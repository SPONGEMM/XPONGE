"""PySCF QM backend adapter."""

from __future__ import annotations

import math
import time

from .._esp_memory import choose_dual_chunk_layout, estimate_aux_tensor_bytes, iter_chunk_slices
from ..capabilities import QMCapabilitySet
from ..errors import QMBackendImportError
from ..models import ESPGridRequest, ESPResult, HessianResult, OptimizationResult, QMMolecule, QMRunOptions, SCFResult


name = "pyscf"
ANGSTROM_PER_BOHR = 0.529177210903
_ESP_TENSOR_ITEMSIZE = 8


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
    kwargs = {
        "atom": _build_atom_block(molecule),
        "verbose": 0,
        "basis": options.basis,
        "charge": molecule.total_charge,
        "spin": molecule.spin,
    }
    if options.ecp is not None:
        kwargs["ecp"] = options.ecp
    if options.cart is not None:
        kwargs["cart"] = bool(options.cart)
    if options.memory_limit_bytes is not None:
        if options.memory_limit_bytes <= 0:
            raise ValueError("QM memory limit should be positive")
        kwargs["max_memory"] = max(1, int(options.memory_limit_bytes) // (1024 * 1024))
    return gto.M(**kwargs)


def _build_wavefunction(mol, molecule: QMMolecule, scf):
    return scf.RHF(mol) if molecule.spin == 0 else scf.UHF(mol)


def _configure_wavefunction(wavefunction, options: QMRunOptions) -> None:
    if options.scf_convergence_tolerance is not None:
        tolerance = float(options.scf_convergence_tolerance)
        if not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError("SCF convergence tolerance should be positive and finite")
        wavefunction.conv_tol = tolerance
    if options.scf_max_cycles is not None:
        if isinstance(options.scf_max_cycles, bool) or int(options.scf_max_cycles) <= 0:
            raise ValueError("SCF maximum cycle count should be a positive integer")
        wavefunction.max_cycle = int(options.scf_max_cycles)


def _configure_threads(options: QMRunOptions) -> None:
    if options.threads is None:
        return
    if isinstance(options.threads, bool) or int(options.threads) <= 0:
        raise ValueError("QM thread count should be a positive integer")
    from pyscf import lib

    lib.num_threads(int(options.threads))


def _collapse_density_matrix(dm, np):
    dm = np.asarray(dm)
    if dm.ndim == 3:
        return dm.sum(axis=0)
    return dm


def _fakemol_for_charges_like_mol(gto, mol, grid_points_bohr):
    fakemol = gto.fakemol_for_charges(grid_points_bohr)
    fakemol.cart = bool(getattr(mol, "cart", False))
    return fakemol


def _compute_esp_full(np, gto, df_module, mol, dm, grid_points_bohr):
    fakemol = _fakemol_for_charges_like_mol(gto, mol, grid_points_bohr)
    return np.einsum("ijp,ij->p", df_module.incore.aux_e2(mol, fakemol), dm)


def _compute_esp_grid_chunked(np, gto, df_module, mol, dm, grid_points_bohr, chunk_size):
    electronic = np.zeros(len(grid_points_bohr), dtype=float)
    for start, stop in iter_chunk_slices(len(grid_points_bohr), chunk_size):
        fakemol = _fakemol_for_charges_like_mol(gto, mol, grid_points_bohr[start:stop])
        electronic[start:stop] = np.einsum("ijp,ij->p", df_module.incore.aux_e2(mol, fakemol), dm)
    return electronic


def _compute_esp_pointwise(np, mol, dm, grid_points_bohr):
    electronic = np.zeros(len(grid_points_bohr), dtype=float)
    for index, point in enumerate(grid_points_bohr):
        mol.set_rinv_orig_(point)
        electronic[index] = np.einsum("ij,ij", mol.intor("int1e_rinv"), dm)
    return electronic


def _compute_esp_shell_grid_chunked(np, gto, df_module, mol, dm, grid_points_bohr, shell_blocks, grid_chunk_size):
    electronic = np.zeros(len(grid_points_bohr), dtype=float)
    for start, stop in iter_chunk_slices(len(grid_points_bohr), grid_chunk_size):
        grid_chunk = grid_points_bohr[start:stop]
        fakemol = _fakemol_for_charges_like_mol(gto, mol, grid_chunk)
        chunk_esp = np.zeros(len(grid_chunk), dtype=float)
        for ish0, ish1, iao0, iao1 in shell_blocks:
            dm_rows = dm[iao0:iao1]
            for jsh0, jsh1, jao0, jao1 in shell_blocks:
                tensor = df_module.incore.aux_e2(
                    mol,
                    fakemol,
                    shls_slice=(ish0, ish1, jsh0, jsh1, 0, len(grid_chunk)),
                )
                chunk_esp += np.einsum("ijp,ij->p", tensor, dm_rows[:, jao0:jao1])
        electronic[start:stop] = chunk_esp
    return electronic


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
    dm = _collapse_density_matrix(wavefunction.make_rdm1(), np)
    nao = int(mol.nao_nr())
    grid_count = len(grid_points_bohr)
    memory_limit_bytes = int(request.memory_limit_bytes)
    usable_bytes = max(1, int(memory_limit_bytes * float(request.safety_factor)))
    estimated_full_bytes = estimate_aux_tensor_bytes(nao, nao, grid_count, _ESP_TENSOR_ITEMSIZE)
    estimated_per_grid_bytes = estimate_aux_tensor_bytes(nao, nao, 1, _ESP_TENSOR_ITEMSIZE)
    start = time.perf_counter()
    diagnostics = {
        "chunk_policy": request.chunk_policy,
        "estimated_full_bytes": estimated_full_bytes,
        "memory_limit_bytes": memory_limit_bytes,
    }
    if grid_count == 0:
        timings = {"esp": time.perf_counter() - start}
        diagnostics.update({"mode": "full", "grid_chunk_count": 0, "shell_block_count": 0})
        return ESPResult(
            grid_points_bohr=grid_points_bohr,
            electronic_esp_au=np.array([], dtype=float),
            timings=timings,
            diagnostics=diagnostics,
        )

    from pyscf import df

    if request.chunk_policy == "full":
        if estimated_full_bytes > usable_bytes:
            raise ValueError("Requested full ESP evaluation exceeds the configured memory budget")
        electronic = _compute_esp_full(np, gto, df, mol, dm, grid_points_bohr)
        diagnostics.update({"mode": "full", "grid_chunk_count": 1, "shell_block_count": 0})
    elif request.chunk_policy == "grid":
        if estimated_per_grid_bytes > usable_bytes:
            raise ValueError("Grid chunking with full AO blocks exceeds the configured memory budget; use chunk_policy='dual' or 'auto'")
        grid_chunk_size = max(1, min(grid_count, usable_bytes // estimated_per_grid_bytes))
        electronic = _compute_esp_grid_chunked(np, gto, df, mol, dm, grid_points_bohr, grid_chunk_size)
        diagnostics.update(
            {
                "mode": "grid_chunk",
                "grid_chunk_size": grid_chunk_size,
                "grid_chunk_count": math.ceil(grid_count / grid_chunk_size),
                "shell_block_count": 0,
            }
        )
    elif request.chunk_policy == "pointwise":
        electronic = _compute_esp_pointwise(np, mol, dm, grid_points_bohr)
        diagnostics.update({"mode": "pointwise", "grid_chunk_count": grid_count, "shell_block_count": 0})
    elif request.chunk_policy == "dual":
        shell_blocks, grid_chunk_size, largest_block_ao = choose_dual_chunk_layout(
            mol.ao_loc_nr(),
            grid_count,
            usable_bytes,
        )
        electronic = _compute_esp_shell_grid_chunked(
            np,
            gto,
            df,
            mol,
            dm,
            grid_points_bohr,
            shell_blocks,
            grid_chunk_size,
        )
        diagnostics.update(
            {
                "mode": "shell_grid_chunk",
                "grid_chunk_size": grid_chunk_size,
                "grid_chunk_count": math.ceil(grid_count / grid_chunk_size),
                "shell_block_count": len(shell_blocks),
                "largest_shell_block_ao": largest_block_ao,
            }
        )
    else:
        if estimated_full_bytes <= usable_bytes:
            electronic = _compute_esp_full(np, gto, df, mol, dm, grid_points_bohr)
            diagnostics.update({"mode": "full", "grid_chunk_count": 1, "shell_block_count": 0})
        elif estimated_per_grid_bytes <= usable_bytes:
            grid_chunk_size = max(1, min(grid_count, usable_bytes // estimated_per_grid_bytes))
            electronic = _compute_esp_grid_chunked(np, gto, df, mol, dm, grid_points_bohr, grid_chunk_size)
            diagnostics.update(
                {
                    "mode": "grid_chunk",
                    "grid_chunk_size": grid_chunk_size,
                    "grid_chunk_count": math.ceil(grid_count / grid_chunk_size),
                    "shell_block_count": 0,
                }
            )
        else:
            shell_blocks, grid_chunk_size, largest_block_ao = choose_dual_chunk_layout(
                mol.ao_loc_nr(),
                grid_count,
                usable_bytes,
            )
            electronic = _compute_esp_shell_grid_chunked(
                np,
                gto,
                df,
                mol,
                dm,
                grid_points_bohr,
                shell_blocks,
                grid_chunk_size,
            )
            diagnostics.update(
                {
                    "mode": "shell_grid_chunk",
                    "grid_chunk_size": grid_chunk_size,
                    "grid_chunk_count": math.ceil(grid_count / grid_chunk_size),
                    "shell_block_count": len(shell_blocks),
                    "largest_shell_block_ao": largest_block_ao,
                }
            )
    timings = {"esp": time.perf_counter() - start}
    return ESPResult(
        grid_points_bohr=grid_points_bohr,
        electronic_esp_au=electronic,
        timings=timings,
        diagnostics=diagnostics,
    )


def optimize_geometry(molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> OptimizationResult:
    optimized = run_scf(
        molecule,
        QMRunOptions(
            backend=options.backend,
            basis=options.basis,
            ecp=options.ecp,
            cart=options.cart,
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
    _configure_threads(options)
    timings = {}
    total_start = time.perf_counter()
    start = time.perf_counter()
    mol = _build_molecule(molecule, options, gto)
    timings["build"] = time.perf_counter() - start
    start = time.perf_counter()
    wavefunction = _build_wavefunction(mol, molecule, scf)
    _configure_wavefunction(wavefunction, options)
    wavefunction.run()
    timings["scf"] = time.perf_counter() - start
    converged = bool(getattr(wavefunction, "converged", False))
    if not converged:
        raise RuntimeError("PySCF SCF did not converge before Hessian evaluation")
    start = time.perf_counter()
    hessian = np.asarray(wavefunction.Hessian().kernel(), dtype=float)
    timings["hessian"] = time.perf_counter() - start
    atom_count = len(molecule.atom_symbols)
    if hessian.shape != (atom_count, atom_count, 3, 3):
        raise RuntimeError(f"PySCF returned an invalid Hessian shape: {hessian.shape!r}")
    flat_hessian = np.transpose(hessian, (0, 2, 1, 3)).reshape(atom_count * 3, atom_count * 3)
    if not np.isfinite(flat_hessian).all():
        raise RuntimeError("PySCF returned a non-finite Cartesian Hessian")
    if not np.allclose(flat_hessian, flat_hessian.T, atol=1.0e-8, rtol=1.0e-7):
        raise RuntimeError("PySCF returned an asymmetric Cartesian Hessian")
    if return_timings:
        timings["total"] = time.perf_counter() - total_start
    return HessianResult(
        cartesian_hessian_au=hessian,
        coordinates_angstrom=[tuple(float(x) for x in row) for row in mol.atom_coords() * ANGSTROM_PER_BOHR],
        atom_symbols=list(molecule.atom_symbols),
        scf_converged=converged,
        total_energy=float(wavefunction.e_tot) if hasattr(wavefunction, "e_tot") else None,
        timings=timings,
    )
