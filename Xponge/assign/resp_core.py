"""
This **module** is used to calculate the RESP charge numerically.
"""

from __future__ import annotations

import numpy as np

from ..helper import get_fibonacci_grid

# Pay Attention To !!!UNIT!!!
default_radius = {"H": 1.2, "C": 1.5, "N": 1.5,
                  "O": 1.4, "P": 1.8, "S": 1.75,
                  "F": 1.35, "Cl": 1.7, "Br": 2.3}


class _RadiusDict(dict):
    def __missing__(self, key):
        raise KeyError("Radius for element {} not found. Maybe you need to set a dict ({element : radii}) "
                       "to the argument variable radius when calling this function.".format(key))


class _RespMol:
    def __init__(self, atom_numbers):
        self.natm = atom_numbers


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
                      extra_equivalence=None, a1=0.0005, a2=0.001, two_stage=True, only_esp=False):
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

    matrix_a0 = np.hstack((matrix_a0, np.ones(mol.natm).reshape(-1, 1)))
    temp = np.ones(mol.natm + 1)
    temp[-1] = 0
    matrix_a0 = np.vstack((matrix_a0, temp.reshape(1, -1)))
    matrix_a = np.zeros_like(matrix_a0)
    matrix_a[:] = matrix_a0
    ainv = np.linalg.inv(matrix_a)

    mep = vnuc - electronic_esp

    matrix_b = np.zeros(mol.natm + 1)
    for i in range(mol.natm):
        r = atom_coordinates_bohr[i]
        rp = np.linalg.norm(r - grid_points_bohr, axis=1)
        matrix_b[i] = np.sum(mep / rp)
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
            for i in range(mol.natm - sublength):
                if assign.atoms[i] != "H":
                    matrix_a[i][i] = a20[i][i] + a / np.sqrt(q_last_step[i] * q_last_step[i] + b * b)
            ainv = np.linalg.inv(matrix_a)
            q_temp = np.dot(ainv, matrix_b)[:-1].reshape(-1)

        q = np.asarray(q).reshape(-1)
        for i, group in enumerate(tofit_second):
            for j in group:
                q[j] = q_temp[i]

    return _force_equivalence_q(q, extra_equivalence)
