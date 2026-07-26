"""
This **module** is used to calculate the RESP charge numerically.
"""

from __future__ import annotations

import numpy as np

from ..helper import get_fibonacci_grid
from ..qm.resp_parameters import RESP_PARAMETERS_BY_SYMBOL, get_resp_mk_radius

# Pay Attention To !!!UNIT!!!
default_radius = {symbol: get_resp_mk_radius(symbol) for symbol in RESP_PARAMETERS_BY_SYMBOL}


class _RadiusDict(dict):
    def __missing__(self, key):
        raise KeyError("Radius for element {} not found. Maybe you need to set a dict ({element : radii}) "
                       "to the argument variable radius when calling this function.".format(key))


class _RespMol:
    def __init__(self, atom_numbers):
        self.natm = atom_numbers


def _prepare_linear_constraints(
        atom_count, charge, extra_equivalence=None,
        constraint_matrix=None, constraint_targets=None):
    """Build a deterministic independent ``C q = d`` constraint system."""
    if extra_equivalence is None:
        extra_equivalence = []
    rows = [np.ones(atom_count, dtype=float)]
    targets = [float(charge)]
    if constraint_matrix is not None:
        matrix = np.asarray(constraint_matrix, dtype=float)
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)
        values = np.asarray(constraint_targets, dtype=float).reshape(-1)
        if matrix.ndim != 2 or matrix.shape[1] != atom_count or matrix.shape[0] != values.size:
            raise ValueError("RESP linear constraints have incompatible dimensions")
        if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(values)):
            raise ValueError("RESP linear constraints should be finite")
        rows.extend(matrix)
        targets.extend(values)
    elif constraint_targets is not None:
        raise ValueError("RESP constraint targets require a constraint matrix")
    for group_index, group in enumerate(extra_equivalence):
        indices = [int(index) for index in group]
        if len(indices) < 2 or len(indices) != len(set(indices)):
            raise ValueError(f"RESP equivalence group {group_index} should contain unique atom indices")
        if min(indices) < 0 or max(indices) >= atom_count:
            raise ValueError(f"RESP equivalence group {group_index} contains an out-of-range atom index")
        reference = indices[0]
        for index in indices[1:]:
            row = np.zeros(atom_count, dtype=float)
            row[index] = 1.0
            row[reference] = -1.0
            rows.append(row)
            targets.append(0.0)

    matrix = np.asarray(rows, dtype=float)
    values = np.asarray(targets, dtype=float)
    matrix_rank = np.linalg.matrix_rank(matrix)
    augmented_rank = np.linalg.matrix_rank(np.column_stack((matrix, values)))
    if augmented_rank != matrix_rank:
        raise ValueError("RESP linear constraints are inconsistent")

    independent_rows = []
    independent_targets = []
    current = np.empty((0, atom_count), dtype=float)
    current_rank = 0
    for row, target in zip(matrix, values):
        candidate = np.vstack((current, row))
        candidate_rank = np.linalg.matrix_rank(candidate)
        if candidate_rank == current_rank:
            continue
        independent_rows.append(row)
        independent_targets.append(target)
        current = candidate
        current_rank = candidate_rank
    return (
        np.asarray(independent_rows, dtype=float),
        np.asarray(independent_targets, dtype=float),
        {
            "input_constraint_count": int(matrix.shape[0]),
            "constraint_rank": int(matrix_rank),
            "dropped_dependent_constraint_count": int(matrix.shape[0] - matrix_rank),
        },
    )


def _solve_constrained_quadratic(matrix_a, matrix_b, constraints, targets):
    """Solve a quadratic system with equality constraints without inversion."""
    matrix_a = np.asarray(matrix_a, dtype=float)
    matrix_b = np.asarray(matrix_b, dtype=float).reshape(-1)
    constraints = np.asarray(constraints, dtype=float)
    targets = np.asarray(targets, dtype=float).reshape(-1)
    atom_count = matrix_b.size
    if matrix_a.shape != (atom_count, atom_count):
        raise ValueError("RESP quadratic matrix has incompatible dimensions")
    if constraints.ndim != 2 or constraints.shape[1] != atom_count:
        raise ValueError("RESP constraint matrix has incompatible dimensions")
    if constraints.shape[0] != targets.size:
        raise ValueError("RESP constraint targets have incompatible dimensions")
    if not all(np.all(np.isfinite(value)) for value in (matrix_a, matrix_b, constraints, targets)):
        raise ValueError("RESP quadratic system should be finite")

    constraint_count = constraints.shape[0]
    kkt = np.block([
        [matrix_a, constraints.T],
        [constraints, np.zeros((constraint_count, constraint_count), dtype=float)],
    ])
    rhs = np.concatenate((matrix_b, targets))
    solution, _, rank, singular_values = np.linalg.lstsq(kkt, rhs, rcond=None)
    residual = kkt @ solution - rhs
    scale = max(1.0, float(np.max(np.abs(rhs))))
    if rank < kkt.shape[0] or float(np.max(np.abs(residual))) > 1.0e-10 * scale:
        raise ValueError("RESP constrained quadratic system is singular or inconsistent")
    charges = solution[:atom_count]
    constraint_residual = constraints @ charges - targets
    condition_number = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and singular_values[-1] > 0
        else float("inf")
    )
    if not np.isfinite(condition_number):
        condition_number = None
    return charges, {
        "kkt_rank": int(rank),
        "kkt_size": int(kkt.shape[0]),
        "condition_number": condition_number,
        "max_constraint_residual": (
            float(np.max(np.abs(constraint_residual))) if constraint_residual.size else 0.0
        ),
    }


def _fit_constrained_resp_stage(
        assign, matrix_a0, matrix_b, constraints, targets, initial_charges,
        restraint, restrained_atom_indices, convergence=1.0e-8, max_iterations=1000):
    charges = np.asarray(initial_charges, dtype=float).reshape(-1)
    restrained_atom_indices = frozenset(restrained_atom_indices)
    for iteration in range(1, max_iterations + 1):
        restrained_matrix = np.array(matrix_a0, dtype=float, copy=True)
        for index in restrained_atom_indices:
            if assign.atoms[index] != "H":
                restrained_matrix[index, index] += (
                    restraint / np.sqrt(charges[index] * charges[index] + 0.01)
                )
        updated, diagnostics = _solve_constrained_quadratic(
            restrained_matrix, matrix_b, constraints, targets
        )
        if float(np.max(np.abs(updated - charges))) <= convergence:
            diagnostics = dict(diagnostics)
            diagnostics["iterations"] = iteration
            return updated, diagnostics
        charges = updated
    raise RuntimeError("RESP constrained iteration did not converge")


def _fit_resp_with_linear_constraints(
        assign, matrix_a0, matrix_b, charge, extra_equivalence,
        constraint_matrix, constraint_targets, a1, a2, two_stage, only_esp):
    atom_count = len(assign.atoms)
    constraints, targets, constraint_diagnostics = _prepare_linear_constraints(
        atom_count,
        charge,
        extra_equivalence=extra_equivalence,
        constraint_matrix=constraint_matrix,
        constraint_targets=constraint_targets,
    )
    charges, initial_diagnostics = _solve_constrained_quadratic(
        matrix_a0, matrix_b, constraints, targets
    )
    diagnostics = {
        **constraint_diagnostics,
        "initial": initial_diagnostics,
        "stage1": None,
        "stage2": None,
    }
    if only_esp:
        diagnostics["max_constraint_residual"] = initial_diagnostics["max_constraint_residual"]
        return charges, diagnostics

    charges, stage1_diagnostics = _fit_constrained_resp_stage(
        assign, matrix_a0, matrix_b, constraints, targets, charges,
        a1, range(atom_count),
    )
    diagnostics["stage1"] = stage1_diagnostics
    if not two_stage:
        diagnostics["max_constraint_residual"] = stage1_diagnostics["max_constraint_residual"]
        return charges, diagnostics

    stage2_groups, _, _ = _find_tofit_second(_RespMol(atom_count), assign)
    if not stage2_groups:
        diagnostics["max_constraint_residual"] = stage1_diagnostics["max_constraint_residual"]
        return charges, diagnostics
    stage2_rows = list(constraints)
    stage2_targets = list(targets)
    active_atoms = set()
    for group in stage2_groups:
        active_atoms.update(group)
        reference = group[0]
        for index in group[1:]:
            row = np.zeros(atom_count, dtype=float)
            row[index] = 1.0
            row[reference] = -1.0
            stage2_rows.append(row)
            stage2_targets.append(0.0)
    for index in range(atom_count):
        if index in active_atoms:
            continue
        row = np.zeros(atom_count, dtype=float)
        row[index] = 1.0
        stage2_rows.append(row)
        stage2_targets.append(charges[index])
    stage2_constraints, stage2_values, stage2_constraint_diagnostics = _prepare_linear_constraints(
        atom_count,
        charge,
        constraint_matrix=np.asarray(stage2_rows, dtype=float),
        constraint_targets=np.asarray(stage2_targets, dtype=float),
    )
    charges, stage2_diagnostics = _fit_constrained_resp_stage(
        assign, matrix_a0, matrix_b, stage2_constraints, stage2_values, charges,
        a2, active_atoms,
    )
    stage2_diagnostics = dict(stage2_diagnostics)
    stage2_diagnostics.update(stage2_constraint_diagnostics)
    diagnostics["stage2"] = stage2_diagnostics
    diagnostics["max_constraint_residual"] = stage2_diagnostics["max_constraint_residual"]
    return charges, diagnostics


# Pay Attention To !!!UNIT!!!
def _get_mk_grid(assign, crd, area_density=1.0, layer=4, radius=None):
    """

    :param assign:
    :param crd:
    :param area_density:
    :param layer:
    :param radius:
    :return:
    """
    grids = []
    factor = area_density * 0.52918 * 0.52918 * 4 * np.pi
    real_radius = _RadiusDict(default_radius)
    if radius:
        real_radius.update(radius)

    lists0 = np.array([1.4 + 0.2 * i for i in range(layer)])
    for i, atom in enumerate(assign.atoms):
        r0 = real_radius[atom] / 0.52918
        for r in r0 * lists0:
            grids.extend([*get_fibonacci_grid(int(factor * r * r), crd[i], r)])
    grids = np.array(grids).reshape(-1, 3)
    for i, atom in enumerate(assign.atoms):
        r0 = 1.39 * real_radius[atom] / 0.52918
        t = np.linalg.norm(grids - crd[i], axis=1)
        grids = grids[t >= r0, :]
    return grids


def get_mk_grid(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
    return _get_mk_grid(assign, atom_coordinates_bohr, area_density=area_density, layer=layer, radius=radius)


def _force_equivalence_q(q, extra_equivalence):
    """

    :param q:
    :param extra_equivalence:
    :return:
    """
    q = np.asarray(q).reshape(-1)
    for eq_group in extra_equivalence:
        q_mean = np.mean(q[eq_group])
        q[eq_group] = q_mean
    return q


def _resp_scf_kernel(mol, assign, a, b, matrix_a, matrix_a0, matrix_b, q):
    """

    :param mol:
    :param assign:
    :param a:
    :param b:
    :param matrix_a:
    :param matrix_a0:
    :param matrix_b:
    :param q:
    :return:
    """
    q = np.asarray(q).reshape(-1)
    matrix_b = np.asarray(matrix_b).reshape(-1)
    step = 0
    q_last_step = q
    while step == 0 or np.max(np.abs(q - q_last_step)) > 1e-4:
        step += 1
        q_last_step = q
        for i in range(mol.natm):
            if assign.atoms[i] != "H":
                matrix_a[i][i] = matrix_a0[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)

        ainv = np.linalg.inv(matrix_a)
        q = np.dot(ainv, matrix_b)
        q = q[:-1]

    return q


def _find_tofit_second(mol, assign):
    """

    :param mol:
    :param assign:
    :return:
    """
    tofit_second = []
    fit_group = {i: -1 for i in range(mol.natm)}
    sublength = 0
    for i in range(mol.natm):
        if assign.Atom_Judge(i, "C4"):
            fit_group[i] = len(tofit_second)
            tofit_second.append([i])
            temp = []
            for j in assign.bonds[i].keys():
                if assign.atoms[j] == "H":
                    temp.append(j)
            if temp:
                for j in temp:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(temp)
                sublength += len(temp) - 1

        if assign.Atom_Judge(i, "C3"):
            temp = []
            for j in assign.bonds[i].keys():
                if assign.atoms[j] == "H":
                    temp.append(j)
            if len(temp) == 2:
                fit_group[i] = len(tofit_second)
                tofit_second.append([i])
                for j in temp:
                    fit_group[j] = len(tofit_second)
                tofit_second.append(temp)
                sublength += 1
    return tofit_second, fit_group, sublength


def _correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence, atom_numbers):
    """

    :param tofit_second:
    :param fit_group:
    :param sublength:
    :param extra_equivalence:
    :param atom_numbers:
    :return:
    """
    if extra_equivalence:
        equi_group = [set() for _ in extra_equivalence]
        for i, eq in enumerate(extra_equivalence):
            for eq_atom in eq:
                if fit_group[eq_atom] != -1:
                    equi_group[i].add(fit_group[eq_atom])
            equi_group[i] = list(equi_group[i])
            equi_group[i].sort()

        all_groups = set()
        for atom in range(atom_numbers):
            all_groups.add(fit_group[atom])
        all_groups_list = list(all_groups)
        all_groups_list.sort()

        group_map = {i: i for i in all_groups_list}
        for eq in equi_group:
            for group in eq:
                group_map[group] = eq[0]

        temp_max = 0
        for group in all_groups_list:
            if group == -1:
                continue
            if group_map[group] == group:
                group_map[group] = temp_max
                temp_max += 1
            else:
                group_map[group] = group_map[group_map[group]]

        temp = tofit_second
        tofit_second = [[] for _ in range(temp_max)]
        for i, group in enumerate(temp):
            tofit_second[group_map[i]].extend(group)
            sublength -= len(group) - 1

        for group in tofit_second:
            sublength += len(group) - 1

        for atom in range(atom_numbers):
            fit_group[atom] = group_map[fit_group[atom]]

    return tofit_second, fit_group, sublength


def _find_restrained_second_stage_groups(assign, tofit_second):
    """

    :param assign:
    :param tofit_second:
    :return:
    """
    restrained = []
    for i, group in enumerate(tofit_second):
        if any(assign.atoms[j] != "H" for j in group):
            restrained.append(i)
    return restrained


def _get_a20_and_b20(total_length, tofit_second, fit_group, sublength, mol, matrix_a0, matrix_b, charge, q):
    """

    :param total_length:
    :param tofit_second:
    :param fit_group:
    :param sublength:
    :param mol:
    :param matrix_a0:
    :param matrix_b:
    :param charge:
    :param q:
    :return:
    """
    q = np.asarray(q).reshape(-1)
    matrix_b = np.asarray(matrix_b).reshape(-1)
    a20 = np.zeros((total_length, total_length))
    count = len(tofit_second)
    for i in range(mol.natm):
        if fit_group[i] == -1:
            fit_group[i] = count
            count += 1
        a20[mol.natm - sublength][fit_group[i]] += 1
        a20[fit_group[i]][mol.natm - sublength] += 1

    b20 = np.zeros(total_length)
    for i in range(mol.natm):
        b20[fit_group[i]] += matrix_b[i]
        for j in range(mol.natm):
            a20[fit_group[i]][fit_group[j]] += matrix_a0[i][j]

    b20[mol.natm - sublength] = charge
    count = 0
    for i in range(mol.natm):
        if fit_group[i] >= len(tofit_second):
            b20[mol.natm - sublength + count + 1] = q[i]
            a20[mol.natm - sublength + count + 1][len(tofit_second) + count] = 1
            a20[len(tofit_second) + count][mol.natm - sublength + count + 1] = 1
            count += 1
    return a20, b20


def fit_resp_from_esp(assign, atom_coordinates_bohr, nuclear_charges, grid_points_bohr, esp_values_au, charge,
                      extra_equivalence=None, a1=0.0005, a2=0.001, two_stage=True, only_esp=False,
                      constraint_matrix=None, constraint_targets=None, return_diagnostics=False):
    """
    Fit RESP charges from backend-neutral ESP inputs.

    :param assign:
    :param atom_coordinates_bohr:
    :param nuclear_charges:
    :param grid_points_bohr:
    :param esp_values_au:
    :param charge:
    :param extra_equivalence:
    :param a1:
    :param a2:
    :param two_stage:
    :param only_esp:
    :param constraint_matrix: optional coefficient matrix for ``C q = d``
    :param constraint_targets: target vector ``d`` for ``constraint_matrix``
    :param return_diagnostics: return charges and numerical diagnostics
    :return:
    """
    if extra_equivalence is None:
        extra_equivalence = []

    atom_coordinates_bohr = np.asarray(atom_coordinates_bohr, dtype=float)
    nuclear_charges = np.asarray(nuclear_charges, dtype=float)
    grid_points_bohr = np.asarray(grid_points_bohr, dtype=float)
    electronic_esp = np.asarray(esp_values_au, dtype=float).reshape(-1)
    mol = _RespMol(len(assign.atoms))

    vnuc = np.zeros(len(grid_points_bohr), dtype=float)
    matrix_a0 = np.zeros((mol.natm, mol.natm))
    for i in range(mol.natm):
        r = atom_coordinates_bohr[i]
        z = nuclear_charges[i]
        rp = r - grid_points_bohr
        for j in range(mol.natm):
            rpj = atom_coordinates_bohr[j] - grid_points_bohr
            matrix_a0[i][j] = np.sum(1.0 / np.linalg.norm(rp, axis=1) / np.linalg.norm(rpj, axis=1))

        vnuc += z / np.einsum("xi,xi->x", rp, rp) ** .5

    mep = vnuc - electronic_esp
    matrix_b_atoms = np.zeros(mol.natm)
    for i in range(mol.natm):
        r = atom_coordinates_bohr[i]
        rp = np.linalg.norm(r - grid_points_bohr, axis=1)
        matrix_b_atoms[i] = np.sum(mep / rp)

    if constraint_matrix is not None or return_diagnostics:
        charges, diagnostics = _fit_resp_with_linear_constraints(
            assign,
            matrix_a0,
            matrix_b_atoms,
            charge,
            extra_equivalence,
            constraint_matrix,
            constraint_targets,
            a1,
            a2,
            two_stage,
            only_esp,
        )
        if return_diagnostics:
            fitted_esp = np.zeros_like(mep, dtype=float)
            for atom_index in range(mol.natm):
                distance = np.linalg.norm(
                    atom_coordinates_bohr[atom_index] - grid_points_bohr,
                    axis=1,
                )
                fitted_esp += float(charges[atom_index]) / distance
            esp_residual = fitted_esp - mep
            esp_rmse = float(np.sqrt(np.mean(esp_residual ** 2)))
            reference_rms = float(np.sqrt(np.mean(mep ** 2)))
            diagnostics.update({
                "esp_point_count": int(mep.size),
                "esp_rmse_au": esp_rmse,
                "esp_relative_rmse": (
                    esp_rmse / reference_rms
                    if reference_rms > np.finfo(float).eps
                    else 0.0
                ),
                "esp_mae_au": float(np.mean(np.abs(esp_residual))),
                "esp_max_abs_error_au": float(np.max(np.abs(esp_residual))),
            })
            return {"charges": charges, "diagnostics": diagnostics}
        return charges

    matrix_a0 = np.hstack((matrix_a0, np.ones(mol.natm).reshape(-1, 1)))
    temp = np.ones(mol.natm + 1)
    temp[-1] = 0
    matrix_a0 = np.vstack((matrix_a0, temp.reshape(1, -1)))
    matrix_a = np.zeros_like(matrix_a0)
    matrix_a[:] = matrix_a0
    ainv = np.linalg.inv(matrix_a)

    matrix_b = np.zeros(mol.natm + 1)
    matrix_b[:mol.natm] = matrix_b_atoms
    matrix_b[-1] = charge

    q = np.dot(ainv, matrix_b.reshape(-1, 1))[:-1]

    if only_esp:
        return _force_equivalence_q(q, extra_equivalence)

    q = _resp_scf_kernel(mol, assign, a1, 0.1, matrix_a, matrix_a0, matrix_b, q)

    if not two_stage:
        return _force_equivalence_q(q, extra_equivalence)

    tofit_second, fit_group, sublength = _find_tofit_second(mol, assign)
    tofit_second, fit_group, sublength = _correct_extra_equivalence(
        tofit_second, fit_group, sublength, extra_equivalence, mol.natm
    )

    if tofit_second:
        total_length = mol.natm - sublength + 1 + mol.natm - sublength - len(tofit_second)

        a20, b20 = _get_a20_and_b20(
            total_length, tofit_second, fit_group, sublength, mol, matrix_a0, matrix_b, charge, q
        )
        restrained_groups = _find_restrained_second_stage_groups(assign, tofit_second)

        matrix_a = np.zeros_like(a20)
        matrix_a[:] = a20[:]
        matrix_b = b20.reshape(-1, 1)
        ainv = np.linalg.inv(matrix_a)
        q_temp = np.dot(ainv, matrix_b)[:-1].reshape(-1)

        a = a2
        b = 0.1
        step = 0
        q_last_step = q_temp
        while step == 0 or np.max(np.abs(q_temp - q_last_step)) > 1e-4:
            step += 1
            q_last_step = q_temp
            for i in restrained_groups:
                matrix_a[i][i] = a20[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
            ainv = np.linalg.inv(matrix_a)
            q_temp = np.dot(ainv, matrix_b)[:-1].reshape(-1)

        q = np.asarray(q).reshape(-1)
        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]

    return _force_equivalence_q(q, extra_equivalence)
