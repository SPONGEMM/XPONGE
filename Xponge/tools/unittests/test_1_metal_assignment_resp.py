"""RESP provider contracts without a Chemcore or Mokda runtime dependency."""

from __future__ import annotations

from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
from unittest import mock
import unittest

from Xponge.metal_assignment import (
    RespEquivalenceGroup,
    RespFitInput,
    ValidationError,
    resp_fit_input_dumps,
    resp_fit_input_loads,
    validate_resp_fit_input,
)
from Xponge.metal_assignment._resp_worker import _execute
from Xponge.metal_assignment.resp_provider import _model_basis_context
from Xponge.metal_assignment._worker_runtime import worker_command


def _fit_input(**updates):
    values = {
        "schema_version": 2,
        "backend": "pyscf",
        "basis_family": "sto-3g",
        "metal_basis_policy": "require_ecp",
        "basis_source": "unit-test-explicit-basis",
        "optimize_geometry": False,
        "grid_density": 1.0,
        "grid_cell_layer": 1,
        "radius_overrides": {},
        "restraint_a1": 0.0005,
        "restraint_a2": 0.001,
        "two_stage": True,
        "only_esp": False,
        "esp_memory_limit_bytes": 64 * 1024 * 1024,
        "esp_chunk_policy": "pointwise",
        "esp_safety_factor": 0.8,
        "equivalence_groups": (),
        "source": "unit-test-explicit-resp-protocol",
    }
    values.update(updates)
    return RespFitInput(**values).with_computed_hash()


def _worker_request(fit_input=None):
    fit_input = fit_input or _fit_input()
    model_id = "large:water"
    return {
        "protocol_version": 1,
        "model": {
            "external_id": model_id,
            "model_hash": "1" * 64,
            "purpose": "large",
            "coordinate_unit": "angstrom",
            "atomic_charge_role": "initial",
            "electronic_state": {
                "selection_id": "water-site",
                "net_charge": 0,
                "spin_multiplicity": 1,
                "source": "unit-test-explicit-state",
            },
            "atoms": [
                {"model_atom_id": "O1", "element": "O", "coordinates": [0.0, 0.0, 0.0], "initial_charge": None, "is_metal": False},
                {"model_atom_id": "H1", "element": "H", "coordinates": [0.7586, 0.0, 0.5043], "initial_charge": None, "is_metal": False},
                {"model_atom_id": "H2", "element": "H", "coordinates": [-0.7586, 0.0, 0.5043], "initial_charge": None, "is_metal": False},
            ],
            "bonds": [
                {"model_atom_ids": ["O1", "H1"], "order": 1.0},
                {"model_atom_ids": ["O1", "H2"], "order": 1.0},
            ],
            "links": [],
            "linear_constraints": [{
                "constraint_id": "core:water",
                "role": "core",
                "atom_ids": ["O1", "H1", "H2"],
                "coefficients": [1.0, 1.0, 1.0],
                "target_charge": 0.0,
                "source": "unit-test-charge-ledger",
            }],
        },
        "fit_protocol": {
            "fit_input_hash": fit_input.fit_input_hash,
            "backend": fit_input.backend,
            "basis_family": fit_input.basis_family,
            "metal_basis_policy": fit_input.metal_basis_policy,
            "basis_source": fit_input.basis_source,
            "optimize_geometry": fit_input.optimize_geometry,
            "grid_density": fit_input.grid_density,
            "grid_cell_layer": fit_input.grid_cell_layer,
            "radius_overrides": dict(fit_input.radius_overrides),
            "restraint_a1": fit_input.restraint_a1,
            "restraint_a2": fit_input.restraint_a2,
            "two_stage": fit_input.two_stage,
            "only_esp": fit_input.only_esp,
            "esp_memory_limit_bytes": fit_input.esp_memory_limit_bytes,
            "esp_chunk_policy": fit_input.esp_chunk_policy,
            "esp_safety_factor": fit_input.esp_safety_factor,
            "equivalence_groups": [
                {"model_id": group.model_id, "atom_ids": list(group.atom_ids), "source": group.source}
                for group in fit_input.equivalence_groups
            ],
            "source": fit_input.source,
        },
    }


def _fake_resp_result(assign, **kwargs):
    return {
        "charges": [-0.8, 0.4, 0.4],
        "metadata": {
            "method": "RESP",
            "backend": kwargs["backend"],
            "basis_family": kwargs["basis"],
            "scf_converged": True,
            "total_energy_hartree": -75.0,
        },
        "diagnostics": {
            "constraint_rank": 1,
            "max_constraint_residual": 0.0,
        },
    }


def _metal_worker_request(*, basis_family="SDD", backend="pyscf"):
    fit_input = _fit_input(
        backend=backend,
        basis_family=basis_family,
        basis_source=f"unit-test-{basis_family}-basis",
    )
    request = _worker_request(fit_input)
    request["model"]["external_id"] = "large:zinc"
    request["model"]["electronic_state"] = {
        "selection_id": "zinc-site",
        "net_charge": 2,
        "spin_multiplicity": 1,
        "source": "unit-test-explicit-zinc-state",
    }
    request["model"]["atoms"] = [
        {
            "model_atom_id": "Zn1",
            "element": "Zn",
            "coordinates": [0.0, 0.0, 0.0],
            "initial_charge": None,
            "is_metal": True,
        }
    ]
    request["model"]["bonds"] = []
    request["model"]["links"] = []
    request["model"]["linear_constraints"] = [{
        "constraint_id": "core:zinc",
        "role": "core",
        "atom_ids": ["Zn1"],
        "coefficients": [1.0],
        "target_charge": 2.0,
        "source": "unit-test-charge-ledger",
    }]
    return request


def _zinc_water_worker_request():
    request = _metal_worker_request()
    request["model"]["external_id"] = "large:zinc-water"
    request["model"]["atoms"] = [
        {
            "model_atom_id": "Zn1",
            "element": "Zn",
            "coordinates": [0.0, 0.0, 0.0],
            "initial_charge": None,
            "is_metal": True,
        },
        {
            "model_atom_id": "O1",
            "element": "O",
            "coordinates": [1.9, 0.0, 0.0],
            "initial_charge": None,
            "is_metal": False,
        },
        {
            "model_atom_id": "H1",
            "element": "H",
            "coordinates": [2.49, 0.76, 0.0],
            "initial_charge": None,
            "is_metal": False,
        },
        {
            "model_atom_id": "H2",
            "element": "H",
            "coordinates": [2.49, -0.76, 0.0],
            "initial_charge": None,
            "is_metal": False,
        },
    ]
    request["model"]["bonds"] = [
        {"model_atom_ids": ["O1", "H1"], "order": 1.0},
        {"model_atom_ids": ["O1", "H2"], "order": 1.0},
    ]
    request["model"]["links"] = [{"model_atom_ids": ["Zn1", "O1"]}]
    request["model"]["linear_constraints"] = [{
        "constraint_id": "core:zinc-water",
        "role": "core",
        "atom_ids": ["Zn1", "O1", "H1", "H2"],
        "coefficients": [1.0, 1.0, 1.0, 1.0],
        "target_charge": 2.0,
        "source": "unit-test-charge-ledger",
    }]
    return request


class RespFitInputTests(unittest.TestCase):
    def test_hash_closed_round_trip_and_explicit_equivalence(self):
        value = _fit_input(
            equivalence_groups=(RespEquivalenceGroup("large:water", ("H1", "H2"), "symmetry"),)
        )
        validate_resp_fit_input(value)
        self.assertEqual(resp_fit_input_loads(resp_fit_input_dumps(value)), value)
        self.assertEqual(resp_fit_input_dumps(resp_fit_input_loads(resp_fit_input_dumps(value))), resp_fit_input_dumps(value))

    def test_geometry_optimization_and_stale_hash_are_rejected(self):
        with self.assertRaisesRegex(ValidationError, "resp_geometry_must_remain_locked"):
            validate_resp_fit_input(_fit_input(optimize_geometry=True))
        value = _fit_input()
        with self.assertRaisesRegex(ValidationError, "stale_resp_fit_input_hash"):
            validate_resp_fit_input(replace(value, basis_family="6-31G*"))

    def test_basis_source_and_metal_policy_are_explicit(self):
        with self.assertRaisesRegex(ValidationError, "missing_resp_basis_source"):
            validate_resp_fit_input(_fit_input(basis_source=""))
        with self.assertRaisesRegex(ValidationError, "unsupported_resp_metal_basis_policy"):
            validate_resp_fit_input(_fit_input(metal_basis_policy="implicit_default"))

    def test_overlapping_equivalence_groups_are_rejected(self):
        with self.assertRaisesRegex(ValidationError, "overlapping_resp_equivalence_groups"):
            validate_resp_fit_input(
                _fit_input(
                    equivalence_groups=(
                        RespEquivalenceGroup("large:water", ("O1", "H1"), "first"),
                        RespEquivalenceGroup("large:water", ("H1", "H2"), "second"),
                    )
                )
            )

    def test_parent_provider_requires_ecp_for_topology_marked_metal(self):
        model = SimpleNamespace(
            external_id="large:zinc",
            atoms=(SimpleNamespace(element="Zn", external_id="zinc:1"),),
        )
        with self.assertRaisesRegex(ValidationError, "resp_metal_ecp_required"):
            _model_basis_context(
                model,
                _fit_input(basis_family="6-31G*"),
                frozenset({"zinc:1"}),
            )
        resolved, metal_elements = _model_basis_context(
            model,
            _fit_input(basis_family="SDD"),
            frozenset({"zinc:1"}),
        )
        self.assertEqual(metal_elements, ("Zn",))
        self.assertEqual(resolved.ecp, {"Zn": "stuttgart_rsc"})


class RespWorkerTests(unittest.TestCase):
    def test_worker_consumes_locked_graph_and_emits_hash_closed_artifact(self):
        with mock.patch("Xponge.assign.resp.resp_fit", side_effect=_fake_resp_result) as fit:
            response = _execute(_worker_request())
        self.assertTrue(response["ok"])
        self.assertEqual(response["artifact"]["atom_order"], ["O1", "H1", "H2"])
        self.assertAlmostEqual(sum(response["artifact"]["charges"]), 0.0, places=12)
        self.assertEqual(response["fit_report"]["backend_spin_2s"], 0)
        self.assertTrue(response["fit_report"]["geometry_locked"])
        self.assertTrue(response["artifact"]["artifact_hash"])
        self.assertTrue(response["response_hash"])
        kwargs = fit.call_args.kwargs
        self.assertFalse(kwargs["opt"])
        self.assertEqual(kwargs["charge"], 0)
        self.assertEqual(kwargs["spin"], 0)
        self.assertEqual(kwargs["constraint_matrix"], [[1.0, 1.0, 1.0]])
        self.assertEqual(kwargs["constraint_targets"], [0.0])

    def test_worker_converts_multiplicity_to_twice_spin(self):
        request = _worker_request()
        request["model"]["atoms"] = [
            {"model_atom_id": "H1", "element": "H", "coordinates": [0.0, 0.0, 0.0], "initial_charge": None, "is_metal": False}
        ]
        request["model"]["bonds"] = []
        request["model"]["linear_constraints"] = [{
            "constraint_id": "core:hydrogen",
            "role": "core",
            "atom_ids": ["H1"],
            "coefficients": [1.0],
            "target_charge": 0.0,
            "source": "unit-test-charge-ledger",
        }]
        request["model"]["electronic_state"]["spin_multiplicity"] = 2
        with mock.patch(
            "Xponge.assign.resp.resp_fit",
            return_value={
                "charges": [0.0],
                "metadata": {"method": "RESP", "scf_converged": True},
                "diagnostics": {
                    "constraint_rank": 1,
                    "max_constraint_residual": 0.0,
                },
            },
        ) as fit:
            response = _execute(request)
        self.assertEqual(response["fit_report"]["backend_spin_2s"], 1)
        self.assertEqual(fit.call_args.kwargs["spin"], 1)

    def test_worker_rejects_duplicate_or_inferred_edges(self):
        request = _worker_request()
        request["model"]["links"] = [{"model_atom_ids": ["H1", "O1"]}]
        with self.assertRaisesRegex(ValueError, "duplicate edge endpoints"):
            _execute(request)

    def test_worker_rejects_metal_basis_without_resolved_ecp(self):
        with self.assertRaisesRegex(ValueError, "resolved ECP does not cover Zn"):
            _execute(_metal_worker_request(basis_family="6-31G*"))

    def test_worker_accepts_explicit_sdd_ecp_for_topology_marked_metal(self):
        with mock.patch(
            "Xponge.assign.resp.resp_fit",
            return_value={
                "charges": [2.0],
                "metadata": {
                    "method": "RESP",
                    "basis_family": "SDD",
                    "basis": {"Zn": "stuttgart_rsc"},
                    "ecp": {"Zn": "stuttgart_rsc"},
                    "cart": True,
                    "scf_converged": True,
                },
                "diagnostics": {
                    "constraint_rank": 1,
                    "max_constraint_residual": 0.0,
                },
            },
        ) as fit:
            response = _execute(_metal_worker_request())
        self.assertEqual(response["fit_report"]["metal_elements"], ["Zn"])
        self.assertEqual(response["fit_report"]["metal_basis_policy"], "require_ecp")
        self.assertEqual(fit.call_args.kwargs["basis"], "SDD")

    @unittest.skipUnless(os.environ.get("MOKDA_PYSCF_PYTHON"), "set MOKDA_PYSCF_PYTHON for real RESP")
    def test_real_pyscf_worker_smoke(self):
        executable = str(Path(os.environ["MOKDA_PYSCF_PYTHON"]).absolute())
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._resp_worker", python_executable=executable),
            input=json.dumps(_worker_request(), sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=180,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)
        response = json.loads(completed.stdout)
        self.assertTrue(response["ok"])
        self.assertTrue(response["fit_report"]["setup_metadata"]["scf_converged"])
        self.assertLess(abs(response["fit_report"]["charge_residual"]), 1.0e-6)

    @unittest.skipUnless(os.environ.get("MOKDA_PYSCF_PYTHON"), "set MOKDA_PYSCF_PYTHON for real RESP")
    def test_real_pyscf_sdd_metal_ecp_smoke(self):
        executable = str(Path(os.environ["MOKDA_PYSCF_PYTHON"]).absolute())
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._resp_worker", python_executable=executable),
            input=json.dumps(_zinc_water_worker_request(), sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=600,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)
        response = json.loads(completed.stdout)
        self.assertTrue(response["ok"])
        self.assertEqual(response["fit_report"]["metal_elements"], ["Zn"])
        metadata = response["fit_report"]["setup_metadata"]
        self.assertEqual(metadata["ecp"], {"Zn": "stuttgart_rsc"})
        self.assertEqual(
            metadata["basis"],
            {"H": "6-31g*", "O": "6-31g*", "Zn": "stuttgart_rsc"},
        )
        self.assertEqual(response["fit_report"]["atom_count"], 4)
        self.assertLess(abs(response["fit_report"]["charge_residual"]), 1.0e-6)


if __name__ == "__main__":
    unittest.main()
