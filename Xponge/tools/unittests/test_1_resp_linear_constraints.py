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
