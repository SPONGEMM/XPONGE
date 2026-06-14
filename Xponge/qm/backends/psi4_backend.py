"""Psi4 QM backend adapter."""

from __future__ import annotations

import time

from ..capabilities import QMCapabilitySet
from ..errors import QMBackendImportError, QMCapabilityError
from ..models import ESPGridRequest, ESPResult, HessianResult, OptimizationResult, QMMolecule, QMRunOptions, SCFResult


name = "psi4"
ANGSTROM_PER_BOHR = 0.52918


def require_numpy_psi4():
    try:
        import numpy as np
    except ImportError as exc:
        raise QMBackendImportError("NumPy is required for QM calculations") from exc
    try:
        import psi4
    except ImportError as exc:
        raise QMBackendImportError("Psi4 is required for QM calculations") from exc
    return np, psi4


def capabilities():
    return QMCapabilitySet(
        supports_scf=True,
        supports_esp=True,
        supports_geometry_optimization=True,
        supports_hessian=False,
        supports_open_shell=True,
    )


def _build_geometry_block(molecule: QMMolecule):
    multiplicity = int(molecule.spin) + 1
    geometry_lines = [f"{int(molecule.total_charge)} {multiplicity}"]
    for atom, (x, y, z) in zip(molecule.atom_symbols, molecule.coordinates_angstrom):
        geometry_lines.append(f"{atom} {x:.16f} {y:.16f} {z:.16f}")
    return "\n".join(geometry_lines)


def _update_assign_coordinates(assign, coordinates_angstrom):
    if assign is None:
        return
    for i, coord in enumerate(coordinates_angstrom):
        assign.coordinate[i][0] = float(coord[0])
        assign.coordinate[i][1] = float(coord[1])
        assign.coordinate[i][2] = float(coord[2])


def run_scf(molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> SCFResult:
    np, psi4 = require_numpy_psi4()
    timings = {}
    total_start = time.perf_counter()
    psi4.core.be_quiet()
    start = time.perf_counter()
    mol = psi4.geometry(_build_geometry_block(molecule))
    psi4.core.clean_options()
    psi4.set_options(
        {
            "basis": options.basis,
            "reference": options.reference or ("rhf" if int(molecule.spin) == 0 else "uhf"),
        }
    )
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    if options.optimize_geometry:
        energy, wavefunction = psi4.optimize(options.method, molecule=mol, return_wfn=True)
    else:
        energy, wavefunction = psi4.energy(options.method, molecule=mol, return_wfn=True)
    timings["scf"] = time.perf_counter() - start
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
    out_mol = wavefunction.molecule()
    atom_coordinates_bohr = np.array(out_mol.geometry().np, dtype=float)
    optimized_coordinates = [tuple(float(x) for x in row) for row in atom_coordinates_bohr * ANGSTROM_PER_BOHR]
    if options.optimize_geometry:
        _update_assign_coordinates(assign, atom_coordinates_bohr * ANGSTROM_PER_BOHR)
    return SCFResult(
        backend_name=name,
        total_energy=float(energy),
        converged=True,
        coordinates_bohr=atom_coordinates_bohr,
        nuclear_charges=np.array([out_mol.Z(i) for i in range(out_mol.natom())], dtype=float),
        charge=molecule.total_charge,
        spin=molecule.spin,
        atom_symbols=list(molecule.atom_symbols),
        backend_handle={"wavefunction": wavefunction},
        timings=timings,
        optimized_coordinates_angstrom=optimized_coordinates if options.optimize_geometry else None,
    )


def compute_esp(scf_result: SCFResult, request: ESPGridRequest) -> ESPResult:
    np, psi4 = require_numpy_psi4()
    grid_points_bohr = np.asarray(request.grid_points_bohr, dtype=float)
    wavefunction = scf_result.backend_handle["wavefunction"]
    start = time.perf_counter()
    grid_matrix = psi4.core.Matrix.from_array(grid_points_bohr * ANGSTROM_PER_BOHR)
    prop = psi4.core.ESPPropCalc(wavefunction)
    total_esp = np.array(prop.compute_esp_over_grid_in_memory(grid_matrix), dtype=float).reshape(-1)
    vnuc = np.zeros(len(grid_points_bohr), dtype=float)
    for coord, charge in zip(scf_result.coordinates_bohr, scf_result.nuclear_charges):
        rp = grid_points_bohr - coord
        vnuc += charge / np.linalg.norm(rp, axis=1)
    timings = {"esp": time.perf_counter() - start}
    return ESPResult(
        grid_points_bohr=grid_points_bohr,
        electronic_esp_au=vnuc - total_esp,
        total_esp_au=total_esp,
        nuclear_esp_au=vnuc,
        timings=timings,
    )


def optimize_geometry(molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> OptimizationResult:
    optimized = run_scf(
        molecule,
        QMRunOptions(
            backend=options.backend,
            basis=options.basis,
            method=options.method,
            reference=options.reference,
            optimize_geometry=True,
            threads=options.threads,
            memory=options.memory,
            properties=options.properties,
        ),
        assign=assign,
        return_timings=return_timings,
    )
    return OptimizationResult(
        optimized_coordinates_angstrom=optimized.optimized_coordinates_angstrom
        or [tuple(float(x) for x in row) for row in molecule.coordinates_angstrom],
        converged=optimized.converged,
        final_energy=optimized.total_energy,
        timings=dict(optimized.timings),
    )


def compute_hessian(molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> HessianResult:
    raise QMCapabilityError(f"{name} does not support Hessian")
