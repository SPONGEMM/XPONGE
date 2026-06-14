"""PySCF QM backend adapter."""

from __future__ import annotations

import time

from ..capabilities import QMCapabilitySet
from ..errors import QMBackendImportError
from ..models import ESPGridRequest, ESPResult, HessianResult, OptimizationResult, QMMolecule, QMRunOptions, SCFResult


name = "pyscf"
ANGSTROM_PER_BOHR = 0.52918


def require_numpy_pyscf():
    try:
        import numpy as np
    except ImportError as exc:
        raise QMBackendImportError("NumPy is required for QM calculations") from exc
    try:
        from pyscf import gto, scf
    except ImportError as exc:
        raise QMBackendImportError("PySCF is required for QM calculations") from exc
    return np, gto, scf


def capabilities():
    return QMCapabilitySet(
        supports_scf=True,
        supports_esp=True,
        supports_geometry_optimization=True,
        supports_hessian=True,
        supports_open_shell=True,
    )


def _build_atom_block(molecule: QMMolecule):
    lines = []
    for atom, (x, y, z) in zip(molecule.atom_symbols, molecule.coordinates_angstrom):
        lines.append(f"{atom} {x:f} {y:f} {z:f}")
    return "\n".join(lines)


def _build_molecule(molecule: QMMolecule, options: QMRunOptions, gto):
    return gto.M(
        atom=_build_atom_block(molecule),
        verbose=0,
        basis=options.basis,
        charge=molecule.total_charge,
        spin=molecule.spin,
    )


def _build_wavefunction(mol, molecule: QMMolecule, scf):
    return scf.RHF(mol) if molecule.spin == 0 else scf.UHF(mol)


def _collapse_density_matrix(dm, np):
    dm = np.asarray(dm)
    if dm.ndim == 3:
        return dm.sum(axis=0)
    return dm


def _update_assign_coordinates(assign, coordinates_angstrom):
    if assign is None:
        return
    for i, coord in enumerate(coordinates_angstrom):
        assign.coordinate[i][0] = float(coord[0])
        assign.coordinate[i][1] = float(coord[1])
        assign.coordinate[i][2] = float(coord[2])


def run_scf(molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> SCFResult:
    np, gto, scf = require_numpy_pyscf()
    timings = {}
    total_start = time.perf_counter()
    start = time.perf_counter()
    mol = _build_molecule(molecule, options, gto)
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    wavefunction = _build_wavefunction(mol, molecule, scf)
    if options.optimize_geometry:
        from pyscf.geomopt.geometric_solver import optimize as geometric_opt

        mol = geometric_opt(wavefunction)
        wavefunction = _build_wavefunction(mol, molecule, scf)
        _update_assign_coordinates(assign, mol.atom_coords() * ANGSTROM_PER_BOHR)
    wavefunction.run()
    timings["scf"] = time.perf_counter() - start
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
    total_energy = None
    if hasattr(wavefunction, "e_tot"):
        total_energy = float(wavefunction.e_tot)
    optimized_coordinates = None
    if options.optimize_geometry:
        optimized_coordinates = [tuple(float(x) for x in row) for row in mol.atom_coords() * ANGSTROM_PER_BOHR]
    return SCFResult(
        backend_name=name,
        total_energy=total_energy,
        converged=bool(getattr(wavefunction, "converged", True)),
        coordinates_bohr=np.array([mol.atom_coord(i) for i in range(mol.natm)]),
        nuclear_charges=np.array([mol.atom_charge(i) for i in range(mol.natm)], dtype=float),
        charge=molecule.total_charge,
        spin=molecule.spin,
        atom_symbols=list(molecule.atom_symbols),
        backend_handle={"mol": mol, "wavefunction": wavefunction},
        timings=timings,
        optimized_coordinates_angstrom=optimized_coordinates,
    )


def compute_esp(scf_result: SCFResult, request: ESPGridRequest) -> ESPResult:
    np, gto, _ = require_numpy_pyscf()
    grid_points_bohr = np.asarray(request.grid_points_bohr, dtype=float)
    mol = scf_result.backend_handle["mol"]
    wavefunction = scf_result.backend_handle["wavefunction"]
    start = time.perf_counter()
    try:
        from pyscf import df

        fakemol = gto.fakemol_for_charges(grid_points_bohr)
        dm = _collapse_density_matrix(wavefunction.make_rdm1(), np)
        electronic = np.einsum("ijp,ij->p", df.incore.aux_e2(mol, fakemol), dm)
    except MemoryError:
        dm = _collapse_density_matrix(wavefunction.make_rdm1(), np)
        vele = []
        for point in grid_points_bohr:
            mol.set_rinv_orig_(point)
            vele.append(np.einsum("ij,ij", mol.intor("int1e_rinv"), dm))
        electronic = np.array(vele)
    timings = {"esp": time.perf_counter() - start}
    return ESPResult(
        grid_points_bohr=grid_points_bohr,
        electronic_esp_au=electronic,
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
    np, gto, scf = require_numpy_pyscf()
    timings = {}
    total_start = time.perf_counter()
    start = time.perf_counter()
    mol = _build_molecule(molecule, options, gto)
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    wavefunction = _build_wavefunction(mol, molecule, scf)
    wavefunction.run()
    timings["scf"] = time.perf_counter() - start
    start = time.perf_counter()
    hessian = np.asarray(wavefunction.Hessian().kernel(), dtype=float)
    timings["hessian"] = time.perf_counter() - start
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
    return HessianResult(
        cartesian_hessian_au=hessian,
        coordinates_angstrom=[tuple(float(x) for x in row) for row in mol.atom_coords() * ANGSTROM_PER_BOHR],
        atom_symbols=list(molecule.atom_symbols),
        timings=timings,
    )
