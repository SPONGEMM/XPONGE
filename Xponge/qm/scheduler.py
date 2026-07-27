"""QM backend selection and dispatch."""

from __future__ import annotations

import sys

from .backends import psi4_backend, pyscf_backend
from ._esp_memory import normalize_chunk_policy, normalize_safety_factor, parse_memory_limit_bytes
from .capabilities import QMCapabilitySet
from .errors import QMBackendImportError, QMBackendSelectionError, QMCapabilityError
from .models import (
    ESPGridRequest,
    HessianResult,
    OptimizationResult,
    QMMolecule,
    QMRunOptions,
    resolve_scf_reference,
)


_BACKENDS = {
    "pyscf": pyscf_backend,
    "psi4": psi4_backend,
}

_SCF_STRATEGIES = frozenset({
    "direct",
    "density_fit",
    "newton",
    "density_fit_newton",
})


def default_backend_name():
    if sys.platform.startswith("win"):
        return "psi4"
    return "pyscf"


def normalize_backend_name(backend):
    if backend is None:
        return default_backend_name()
    backend_name = str(backend).strip().lower()
    if backend_name not in _BACKENDS:
        supported = ", ".join(sorted(_BACKENDS))
        raise QMBackendSelectionError(f"QM backend should be one of: {supported}")
    return backend_name


def backend_import_or_hint(backend_name, exc):
    message = str(exc)
    if backend_name == "pyscf" and sys.platform.startswith("win"):
        message += " On Windows, install Psi4 via conda-forge or the official Psi4 installer and call calculate_charge('resp', backend='psi4', ...)."
    elif backend_name == "psi4":
        message += " On Windows, Psi4 is not installed through pip by default; install it via conda-forge or the official Psi4 installer and retry."
    raise QMBackendImportError(message) from exc


def get_backend(backend=None):
    return _BACKENDS[normalize_backend_name(backend)]


def normalize_scf_strategy(backend, strategy):
    backend_name = normalize_backend_name(backend)
    strategy_name = str(strategy).strip().lower()
    if strategy_name not in _SCF_STRATEGIES:
        supported = ", ".join(sorted(_SCF_STRATEGIES))
        raise ValueError(f"SCF strategy should be one of: {supported}")
    if backend_name != "pyscf" and strategy_name != "direct":
        raise ValueError(
            f"SCF strategy {strategy_name!r} requires the PySCF backend"
        )
    return strategy_name


def get_capabilities(backend=None) -> QMCapabilitySet:
    return get_backend(backend).capabilities()


def qmmolecule_from_assign(assign, charge, spin):
    atom_names = []
    if hasattr(assign, "names"):
        atom_names = list(assign.names)
    formal_charges = []
    if hasattr(assign, "formal_charge"):
        formal_charges = [int(x) for x in assign.formal_charge]
    bonds = [{int(k): int(v) for k, v in bond_map.items()} for bond_map in assign.bonds.values()]
    return QMMolecule(
        atom_symbols=list(assign.atoms),
        coordinates_angstrom=[tuple(float(x) for x in coord) for coord in assign.coordinate],
        total_charge=int(charge),
        spin=int(spin),
        atom_names=atom_names or None,
        formal_charges=formal_charges or None,
        bonds=bonds,
        metadata={"assign_name": getattr(assign, "name", None)},
    )


def run_scf(
    assign,
    *,
    backend=None,
    basis="6-31g*",
    ecp=None,
    cart=None,
    charge=0,
    spin=0,
    optimize_geometry=False,
    return_timings=False,
    scf_strategy="direct",
    scf_reference="auto",
):
    backend_name = normalize_backend_name(backend)
    scf_strategy = normalize_scf_strategy(backend_name, scf_strategy)
    backend_module = get_backend(backend_name)
    molecule = qmmolecule_from_assign(assign, charge, spin)
    resolved_reference = resolve_scf_reference(scf_reference, molecule.spin)
    options = QMRunOptions(
        backend=backend_name,
        basis=basis,
        ecp=ecp,
        cart=cart,
        method="scf",
        reference=resolved_reference,
        optimize_geometry=optimize_geometry,
        scf_strategy=scf_strategy,
    )
    try:
        if not backend_module.capabilities().supports_scf:
            raise QMCapabilityError(f"{backend_name} does not support SCF")
        return backend_module.run_scf(molecule, options, assign=assign, return_timings=return_timings)
    except ImportError as exc:
        backend_import_or_hint(backend_name, exc)


def compute_esp_on_grid(scf_result, grid_points_bohr, *, memory_limit=None, chunk_policy="auto", safety_factor=0.8):
    backend_module = get_backend(scf_result.backend_name)
    try:
        if not backend_module.capabilities().supports_esp:
            raise QMCapabilityError(f"{scf_result.backend_name} does not support ESP")
        return backend_module.compute_esp(
            scf_result,
            ESPGridRequest(
                grid_points_bohr=grid_points_bohr,
                memory_limit_bytes=parse_memory_limit_bytes(memory_limit),
                chunk_policy=normalize_chunk_policy(chunk_policy),
                safety_factor=normalize_safety_factor(safety_factor),
            ),
        )
    except ImportError as exc:
        backend_import_or_hint(scf_result.backend_name, exc)


def optimize_geometry(
    assign,
    *,
    backend=None,
    basis="6-31g*",
    ecp=None,
    cart=None,
    charge=0,
    spin=0,
    return_timings=False,
    scf_reference="auto",
) -> OptimizationResult:
    backend_name = normalize_backend_name(backend)
    backend_module = get_backend(backend_name)
    molecule = qmmolecule_from_assign(assign, charge, spin)
    resolved_reference = resolve_scf_reference(scf_reference, molecule.spin)
    options = QMRunOptions(
        backend=backend_name,
        basis=basis,
        ecp=ecp,
        cart=cart,
        method="scf",
        reference=resolved_reference,
        optimize_geometry=True,
    )
    try:
        if not backend_module.capabilities().supports_geometry_optimization:
            raise QMCapabilityError(f"{backend_name} does not support geometry optimization")
        return backend_module.optimize_geometry(molecule, options, assign=assign, return_timings=return_timings)
    except ImportError as exc:
        backend_import_or_hint(backend_name, exc)


def compute_hessian(
    assign,
    *,
    backend=None,
    basis="6-31g*",
    ecp=None,
    cart=None,
    charge=0,
    spin=0,
    threads=None,
    memory_limit_bytes=None,
    scf_convergence_tolerance=None,
    scf_max_cycles=None,
    return_timings=False,
    scf_reference="auto",
) -> HessianResult:
    backend_name = normalize_backend_name(backend)
    backend_module = get_backend(backend_name)
    molecule = qmmolecule_from_assign(assign, charge, spin)
    resolved_reference = resolve_scf_reference(scf_reference, molecule.spin)
    options = QMRunOptions(
        backend=backend_name,
        basis=basis,
        ecp=ecp,
        cart=cart,
        method="scf",
        reference=resolved_reference,
        optimize_geometry=False,
        threads=threads,
        memory_limit_bytes=memory_limit_bytes,
        scf_convergence_tolerance=scf_convergence_tolerance,
        scf_max_cycles=scf_max_cycles,
    )
    try:
        if not backend_module.capabilities().supports_hessian:
            raise QMCapabilityError(f"{backend_name} does not support Hessian")
        return backend_module.compute_hessian(molecule, options, assign=assign, return_timings=return_timings)
    except ImportError as exc:
        backend_import_or_hint(backend_name, exc)
