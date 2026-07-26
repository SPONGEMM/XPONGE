"""Numerical tests for general RESP linear charge constraints."""

from __future__ import annotations

import numpy as np
import pytest

from Xponge.assign.resp_core import (
    fit_resp_from_esp,
    _prepare_linear_constraints,
    _solve_constrained_quadratic,
)


class _ThreeAtomAssign:
    atoms = ["C", "H", "H"]
    bonds = [{1: 1.0, 2: 1.0}, {0: 1.0}, {0: 1.0}]

    @staticmethod
    def Atom_Judge(index, atom_type):
        return index == 0 and atom_type == "C4"


def test_constrained_quadratic_recovers_known_group_targets():
    constraints, targets, constraint_report = _prepare_linear_constraints(
        3,
        1,
        constraint_matrix=[
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        constraint_targets=[0.0, 1.0],
    )
    charges, solve_report = _solve_constrained_quadratic(
        np.eye(3),
        np.asarray([1.0, 2.0, 3.0]),
        constraints,
        targets,
    )

    assert charges == pytest.approx([-0.5, 0.5, 1.0], abs=1.0e-12)
    assert constraint_report["constraint_rank"] == 2
    assert constraint_report["dropped_dependent_constraint_count"] == 1
    assert solve_report["max_constraint_residual"] <= 1.0e-12


def test_equivalence_is_part_of_the_constraint_system():
    constraints, targets, report = _prepare_linear_constraints(
        3,
        0,
        extra_equivalence=[[0, 1]],
    )
    charges, solve_report = _solve_constrained_quadratic(
        np.eye(3),
        np.asarray([1.0, 3.0, -2.0]),
        constraints,
        targets,
    )

    assert charges[0] == pytest.approx(charges[1], abs=1.0e-12)
    assert sum(charges) == pytest.approx(0.0, abs=1.0e-12)
    assert report["constraint_rank"] == 2
    assert solve_report["max_constraint_residual"] <= 1.0e-12


def test_conflicting_constraints_are_rejected():
    with pytest.raises(ValueError, match="inconsistent"):
        _prepare_linear_constraints(
            2,
            0,
            constraint_matrix=[
                [1.0, 0.0],
                [1.0, 0.0],
            ],
            constraint_targets=[0.0, 1.0],
        )


def test_large_partition_and_reference_constraints_are_reduced_stably():
    atom_count = 90
    core_count = 54
    matrix = []
    targets = []
    core = np.zeros(atom_count)
    core[:core_count] = 1.0
    matrix.append(core)
    targets.append(0.0)
    for fragment_index in range(6):
        cap = np.zeros(atom_count)
        start = core_count + fragment_index * 6
        cap[start:start + 6] = 1.0
        matrix.append(cap)
        targets.append(0.0)
    for atom_index, target in enumerate((
        -0.4157, 0.2719, 0.0337, 0.0823, 0.5973, -0.5679,
        -0.3479, 0.2747, -0.2637, 0.1560, 0.7341, -0.5894,
        -0.4157, 0.2719, 0.0337, 0.0823, 0.5973, -0.5679,
    )):
        fixed = np.zeros(atom_count)
        fixed[atom_index] = 1.0
        matrix.append(fixed)
        targets.append(target)

    constraints, reduced_targets, report = _prepare_linear_constraints(
        atom_count,
        0,
        constraint_matrix=matrix,
        constraint_targets=targets,
    )

    assert report["input_constraint_count"] == 26
    assert report["constraint_rank"] == 25
    assert report["dropped_dependent_constraint_count"] == 1
    assert constraints.shape == (25, atom_count)
    candidate, _, _, _ = np.linalg.lstsq(constraints, reduced_targets, rcond=None)
    assert constraints @ candidate == pytest.approx(reduced_targets, abs=1.0e-10)


def test_non_finite_constraints_are_rejected():
    with pytest.raises(ValueError, match="finite"):
        _prepare_linear_constraints(
            2,
            0,
            constraint_matrix=[[1.0, np.nan]],
            constraint_targets=[0.0],
        )


def test_two_stage_fit_preserves_fixed_group_and_equivalence_constraints():
    coordinates = np.asarray([
        [0.0, 0.0, 0.0],
        [1.8, 0.0, 0.0],
        [-1.8, 0.0, 0.0],
    ])
    grid = np.asarray([
        [0.0, 3.0, 0.0],
        [0.0, -3.0, 0.0],
        [0.0, 0.0, 3.0],
        [2.5, 2.5, 0.0],
        [-2.5, 2.5, 0.0],
        [0.0, -2.5, 2.5],
    ])
    target = np.asarray([-0.2, 0.1, 0.1])
    nuclear_charges = np.asarray([6.0, 1.0, 1.0])
    inverse_distance = np.asarray([
        1.0 / np.linalg.norm(grid - coordinate, axis=1)
        for coordinate in coordinates
    ])
    nuclear_esp = nuclear_charges @ inverse_distance
    target_mep = target @ inverse_distance
    electronic_esp = nuclear_esp - target_mep

    result = fit_resp_from_esp(
        _ThreeAtomAssign(),
        coordinates,
        nuclear_charges,
        grid,
        electronic_esp,
        charge=0,
        extra_equivalence=[[1, 2]],
        constraint_matrix=[
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
        ],
        constraint_targets=[-0.2, 0.2],
        two_stage=True,
        return_diagnostics=True,
    )

    assert result["charges"] == pytest.approx(target, abs=1.0e-8)
    assert result["diagnostics"]["stage1"]["max_constraint_residual"] <= 1.0e-10
    assert result["diagnostics"]["stage2"]["max_constraint_residual"] <= 1.0e-10
    assert result["diagnostics"]["esp_point_count"] == len(grid)
    assert result["diagnostics"]["esp_rmse_au"] <= 1.0e-8
    assert result["diagnostics"]["esp_relative_rmse"] <= 1.0e-8
    assert result["diagnostics"]["esp_mae_au"] <= 1.0e-8
    assert result["diagnostics"]["esp_max_abs_error_au"] <= 1.0e-8


def test_constrained_fit_is_repeatable_and_atom_order_invariant():
    coordinates = np.asarray([
        [0.0, 0.0, 0.0],
        [1.8, 0.0, 0.0],
        [-1.8, 0.0, 0.0],
    ])
    grid = np.asarray([
        [0.0, 3.0, 0.0],
        [0.0, -3.0, 0.0],
        [0.0, 0.0, 3.0],
        [2.5, 2.5, 0.0],
        [-2.5, 2.5, 0.0],
        [0.0, -2.5, 2.0],
    ])
    target = np.asarray([-0.2, 0.1, 0.1])
    nuclear_charges = np.asarray([6.0, 1.0, 1.0])
    inverse_distance = np.asarray([
        1.0 / np.linalg.norm(grid - coordinate, axis=1)
        for coordinate in coordinates
    ])
    electronic_esp = nuclear_charges @ inverse_distance - target @ inverse_distance

    def fit(assign, atom_coordinates, atom_nuclear_charges, carbon_index, hydrogen_indices):
        fixed_carbon = np.zeros(3)
        fixed_carbon[carbon_index] = 1.0
        hydrogen_group = np.zeros(3)
        hydrogen_group[list(hydrogen_indices)] = 1.0
        return np.asarray(fit_resp_from_esp(
            assign,
            atom_coordinates,
            atom_nuclear_charges,
            grid,
            electronic_esp,
            charge=0,
            extra_equivalence=[list(hydrogen_indices)],
            constraint_matrix=[fixed_carbon, hydrogen_group],
            constraint_targets=[-0.2, 0.2],
            two_stage=True,
        ))

    first = fit(_ThreeAtomAssign(), coordinates, nuclear_charges, 0, (1, 2))
    repeated = fit(_ThreeAtomAssign(), coordinates, nuclear_charges, 0, (1, 2))
    assert repeated == pytest.approx(first, abs=1.0e-12)

    permutation = np.asarray([2, 0, 1])
    inverse_permutation = np.argsort(permutation)

    class _PermutedAssign:
        atoms = [_ThreeAtomAssign.atoms[index] for index in permutation]
        bonds = []
        for old_index in permutation:
            bonds.append({
                int(inverse_permutation[old_neighbor]): order
                for old_neighbor, order in _ThreeAtomAssign.bonds[old_index].items()
            })

        @staticmethod
        def Atom_Judge(index, atom_type):
            return _PermutedAssign.atoms[index] == "C" and atom_type == "C4"

    permuted = fit(
        _PermutedAssign(),
        coordinates[permutation],
        nuclear_charges[permutation],
        int(inverse_permutation[0]),
        (int(inverse_permutation[1]), int(inverse_permutation[2])),
    )
    mapped_back = permuted[inverse_permutation]
    assert mapped_back == pytest.approx(first, abs=1.0e-12)
