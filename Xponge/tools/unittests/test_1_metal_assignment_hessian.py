"""Analytical-Hessian provider contracts without Mokda or Chemcore imports."""

from __future__ import annotations

from dataclasses import replace
import os
import signal
import subprocess
from types import SimpleNamespace
from unittest import mock
import unittest

import numpy as np

from Xponge.metal_assignment import (
    DerivedModel,
    HESSIAN_FIT_INPUT_SCHEMA_VERSION,
    HessianFitInput,
    ModelAtom,
    ModelBond,
    ModelElectronicState,
    ValidationError,
    fit_model_hessians,
    hessian_fit_input_dumps,
    hessian_fit_input_loads,
    validate_hessian_fit_input,
)
from Xponge.metal_assignment._hessian_worker import _execute
from Xponge.metal_assignment._ion_worker import _canonical_hash
from Xponge.metal_assignment._worker_runtime import run_worker_subprocess
from Xponge.metal_assignment.hessian_provider import _invoke_hessian_worker
from Xponge.qm import scheduler as qm_scheduler
from Xponge.qm.models import HessianResult, QMMolecule, QMRunOptions


def _fit_input(**updates):
    values = {
        "schema_version": HESSIAN_FIT_INPUT_SCHEMA_VERSION,
        "backend": "pyscf",
        "method": "scf",
        "basis_family": "sto-3g",
        "metal_basis_policy": "require_ecp",
        "basis_source": "unit-test-explicit-basis",
        "optimize_geometry": False,
        "scf_reference": "auto",
        "scf_convergence_tolerance": 1.0e-9,
        "scf_max_cycles": 80,
        "threads": 1,
        "memory_limit_bytes": 256 * 1024 * 1024,
        "source": "unit-test-explicit-hessian-protocol",
    }
    values.update(updates)
    return HessianFitInput(**values).with_computed_hash()


def _raw_identity_hessian(atom_count):
    flat = np.eye(atom_count * 3, dtype=float)
    return np.transpose(
        flat.reshape(atom_count, 3, atom_count, 3), (0, 2, 1, 3)
    )


def _worker_request(fit_input=None):
    fit_input = fit_input or _fit_input()
    return {
        "protocol_version": 1,
        "model": {
            "external_id": "small:water",
            "model_hash": "1" * 64,
            "purpose": "small",
            "coordinate_unit": "angstrom",
            "atomic_charge_role": "absent",
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
            "linear_constraints": [],
        },
        "fit_protocol": {
            "fit_input_hash": fit_input.fit_input_hash,
            "backend": fit_input.backend,
            "method": fit_input.method,
            "basis_family": fit_input.basis_family,
            "metal_basis_policy": fit_input.metal_basis_policy,
            "basis_source": fit_input.basis_source,
            "optimize_geometry": fit_input.optimize_geometry,
            "scf_reference": fit_input.scf_reference,
            "scf_convergence_tolerance": fit_input.scf_convergence_tolerance,
            "scf_max_cycles": fit_input.scf_max_cycles,
            "threads": fit_input.threads,
            "memory_limit_bytes": fit_input.memory_limit_bytes,
            "source": fit_input.source,
        },
    }


def _fake_hessian_result(assign, **kwargs):
    return HessianResult(
        cartesian_hessian_au=_raw_identity_hessian(3),
        coordinates_angstrom=[tuple(float(item) for item in row) for row in assign.coordinate],
        atom_symbols=list(assign.atoms),
        scf_converged=True,
        total_energy=-75.0,
        reference="rhf",
        timings={"scf": 0.1, "hessian": 0.2},
    )


def _model_atom(model_atom_id, external_id, serial, element, coordinates):
    return ModelAtom(
        model_atom_id=model_atom_id,
        external_id=external_id,
        canonical_atom_id=serial,
        mol2_serial=serial,
        element=element,
        coordinates=coordinates,
        role="core",
        canonical_residue_id="water:1",
        chemical_component_id="water",
        partial_charge=None,
        cap_parent_external_id="",
        cut_edge_id="",
        geometry_provenance="canonical",
        charge_projection_group="",
    )


def _small_model():
    charge_accounting_payload = {
        "schema_version": 1,
        "graph_revision": 1,
        "selection_id": "water-site",
        "net_charge": 0,
        "complete": True,
        "contributions": [],
    }
    model = DerivedModel(
        external_id="small:water",
        site_id="water-site",
        purpose="small",
        coordinate_unit="angstrom",
        atomic_charge_role="absent",
        electronic_state=ModelElectronicState("water-site", 0, 1, "unit-test-explicit-state"),
        atoms=(
            _model_atom("O1", "water:O", 1, "O", (0.0, 0.0, 0.0)),
            _model_atom("H1", "water:H1", 2, "H", (0.7586, 0.0, 0.5043)),
            _model_atom("H2", "water:H2", 3, "H", (-0.7586, 0.0, 0.5043)),
        ),
        bonds=(
            ModelBond("bond:OH1", ("O1", "H1"), 1.0, "covalent", "unit-test"),
            ModelBond("bond:OH2", ("O1", "H2"), 1.0, "covalent", "unit-test"),
        ),
        links=(),
        cut_edges=(),
        charge_accounting={
            **charge_accounting_payload,
            "proof_hash": _canonical_hash(charge_accounting_payload),
        },
        mol2_text="@<TRIPOS>MOLECULE\nsmall:water\n",
        model_hash="",
    )
    return replace(model, model_hash=model.computed_hash())


def _package(*models):
    topology_atoms = tuple(
        SimpleNamespace(external_id=atom.external_id, is_metal=False)
        for model in models
        for atom in model.atoms
    )
    return SimpleNamespace(
        request=SimpleNamespace(
            request_id="request:hessian",
            input_hash="2" * 64,
            graph_revision=1,
            topology=SimpleNamespace(atoms=topology_atoms, topology_hash="3" * 64),
            projection_hash="4" * 64,
        ),
        package_hash="5" * 64,
        prepared_artifacts=SimpleNamespace(
            derived_models=SimpleNamespace(models=models)
        ),
    )


class HessianFitInputTests(unittest.TestCase):
    def test_hash_closed_round_trip_and_operation_registration(self):
        value = _fit_input()
        validate_hessian_fit_input(value)
        self.assertEqual(hessian_fit_input_loads(hessian_fit_input_dumps(value)), value)

    def test_geometry_resource_and_stale_hash_fail_closed(self):
        with self.assertRaisesRegex(ValidationError, "hessian_geometry_must_remain_locked"):
            validate_hessian_fit_input(_fit_input(optimize_geometry=True))
        with self.assertRaisesRegex(ValidationError, "invalid_hessian_threads"):
            validate_hessian_fit_input(_fit_input(threads=0))
        with self.assertRaisesRegex(ValidationError, "invalid_hessian_memory_limit_bytes"):
            validate_hessian_fit_input(_fit_input(memory_limit_bytes=0))
        value = _fit_input()
        with self.assertRaisesRegex(ValidationError, "stale_hessian_fit_input_hash"):
            validate_hessian_fit_input(replace(value, scf_max_cycles=81))

    def test_scf_reference_is_hash_closed_and_explicitly_validated(self):
        for reference in ("auto", "rhf", "rohf", "uhf"):
            validate_hessian_fit_input(_fit_input(scf_reference=reference))
        with self.assertRaisesRegex(ValidationError, "unsupported_hessian_scf_reference"):
            validate_hessian_fit_input(_fit_input(scf_reference="implicit"))


class HessianWorkerTests(unittest.TestCase):
    def test_worker_consumes_locked_small_model_and_emits_hash_closed_response(self):
        with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result) as compute:
            response = _execute(_worker_request())
        self.assertTrue(response["ok"])
        self.assertEqual(response["artifact"]["atom_order"], ["O1", "H1", "H2"])
        self.assertEqual(response["response_hash"], _canonical_hash({
            key: value for key, value in response.items() if key != "response_hash"
        }))
        kwargs = compute.call_args.kwargs
        self.assertEqual(kwargs["spin"], 0)
        self.assertEqual(kwargs["threads"], 1)
        self.assertEqual(kwargs["memory_limit_bytes"], 256 * 1024 * 1024)
        self.assertEqual(kwargs["scf_convergence_tolerance"], 1.0e-9)
        self.assertEqual(kwargs["scf_reference"], "auto")
        self.assertTrue(response["fit_report"]["setup_metadata"]["scf_converged"])

    def test_unconverged_result_is_rejected(self):
        result = _fake_hessian_result(SimpleNamespace(
            coordinate=[(0.0, 0.0, 0.0), (0.7586, 0.0, 0.5043), (-0.7586, 0.0, 0.5043)],
            atoms=["O", "H", "H"],
        ))
        result.scf_converged = False
        with mock.patch("Xponge.qm.scheduler.compute_hessian", return_value=result):
            with self.assertRaisesRegex(RuntimeError, "unconverged"):
                _execute(_worker_request())


class HessianParentProviderTests(unittest.TestCase):
    def test_parent_validates_worker_identity_coordinates_and_output_hash(self):
        model = _small_model()

        def invoke(payload, **_kwargs):
            with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result):
                return _execute(payload)

        with mock.patch("Xponge.metal_assignment.hessian_provider.validate_package"), mock.patch(
            "Xponge.metal_assignment.hessian_provider._invoke_hessian_worker", side_effect=invoke
        ):
            output = fit_model_hessians(_package(model), _fit_input())
        self.assertEqual(len(output.hessian_artifacts), 1)
        self.assertEqual(output.hessian_artifacts[0].model_hash, model.model_hash)
        self.assertEqual(output.output_hash, output.computed_hash())
        self.assertEqual(len(output.worker_response_hashes[0]), 64)

    def test_parent_covers_multiple_small_models_in_stable_site_order(self):
        first = _small_model()
        second = replace(
            first,
            external_id="small:water-2",
            site_id="water-site-2",
            model_hash="",
        )
        second = replace(second, model_hash=second.computed_hash())

        def invoke(payload, **_kwargs):
            with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result):
                return _execute(payload)

        with mock.patch("Xponge.metal_assignment.hessian_provider.validate_package"), mock.patch(
            "Xponge.metal_assignment.hessian_provider._invoke_hessian_worker", side_effect=invoke
        ):
            output = fit_model_hessians(_package(second, first), _fit_input())
        self.assertEqual(
            [artifact.model_id for artifact in output.hessian_artifacts],
            ["small:water", "small:water-2"],
        )
        self.assertEqual(len(output.fit_reports), 2)
        self.assertEqual(len(output.worker_response_hashes), 2)

    def test_parent_rejects_forged_artifact_identity_and_malformed_response(self):
        model = _small_model()

        def forged(payload, **_kwargs):
            with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result):
                response = _execute(payload)
            response["artifact"]["model_hash"] = "0" * 64
            response["response_hash"] = _canonical_hash({
                key: value for key, value in response.items() if key != "response_hash"
            })
            return response

        with mock.patch("Xponge.metal_assignment.hessian_provider.validate_package"), mock.patch(
            "Xponge.metal_assignment.hessian_provider._invoke_hessian_worker", side_effect=forged
        ):
            with self.assertRaisesRegex(ValidationError, "hessian_model_identity_mismatch"):
                fit_model_hessians(_package(model), _fit_input())

        def missing_field(payload, **_kwargs):
            with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result):
                response = _execute(payload)
            response.pop("fit_report")
            response["response_hash"] = _canonical_hash({
                key: value for key, value in response.items() if key != "response_hash"
            })
            return response

        with mock.patch("Xponge.metal_assignment.hessian_provider.validate_package"), mock.patch(
            "Xponge.metal_assignment.hessian_provider._invoke_hessian_worker", side_effect=missing_field
        ):
            with self.assertRaisesRegex(ValidationError, "invalid_hessian_worker_response"):
                fit_model_hessians(_package(model), _fit_input())

    def test_parent_rejects_worker_coordinate_drift_even_with_rehashed_response(self):
        model = _small_model()

        def invoke(payload, **_kwargs):
            with mock.patch("Xponge.qm.scheduler.compute_hessian", side_effect=_fake_hessian_result):
                response = _execute(payload)
            response["artifact"]["coordinates_angstrom"][0][0] = 0.01
            response["response_hash"] = _canonical_hash({
                key: value for key, value in response.items() if key != "response_hash"
            })
            return response

        with mock.patch("Xponge.metal_assignment.hessian_provider.validate_package"), mock.patch(
            "Xponge.metal_assignment.hessian_provider._invoke_hessian_worker", side_effect=invoke
        ):
            with self.assertRaisesRegex(ValidationError, "hessian_coordinate_identity_mismatch"):
                fit_model_hessians(_package(model), _fit_input())

    def test_parent_translates_worker_timeout_to_stable_error(self):
        with mock.patch(
            "Xponge.metal_assignment.hessian_provider.run_worker_subprocess",
            side_effect=subprocess.TimeoutExpired(["python"], 1.0),
        ):
            with self.assertRaisesRegex(ValidationError, "hessian_worker_timeout"):
                _invoke_hessian_worker(
                    _worker_request(), timeout_seconds=1.0, python_executable=None
                )

class WorkerRuntimeTests(unittest.TestCase):
    @unittest.skipUnless(os.name == "posix", "POSIX process-group cleanup contract")
    def test_parent_interrupt_kills_and_reaps_worker_process_group(self):
        process = mock.Mock()
        process.pid = 4242
        process.communicate.side_effect = [
            KeyboardInterrupt(),
            ("", "phase-progress"),
        ]
        with mock.patch("subprocess.Popen", return_value=process), mock.patch(
            "os.killpg"
        ) as killpg:
            with self.assertRaises(KeyboardInterrupt):
                run_worker_subprocess(
                    ["python", "-c", "pass"], input_text="{}", timeout_seconds=1.0
                )
        killpg.assert_called_once_with(4242, signal.SIGKILL)
        self.assertEqual(process.communicate.call_count, 2)

    @unittest.skipUnless(os.name == "posix", "POSIX process-group cleanup contract")
    def test_timeout_kills_and_reaps_worker_process_group(self):
        process = mock.Mock()
        process.pid = 4242
        process.communicate.side_effect = [
            subprocess.TimeoutExpired(["python"], 1.0),
            ("", "phase-progress"),
        ]
        with mock.patch("subprocess.Popen", return_value=process) as popen, mock.patch(
            "os.killpg"
        ) as killpg:
            with self.assertRaises(subprocess.TimeoutExpired) as captured:
                run_worker_subprocess(
                    ["python", "-c", "pass"], input_text="{}", timeout_seconds=1.0
                )
        self.assertTrue(popen.call_args.kwargs["start_new_session"])
        killpg.assert_called_once_with(4242, signal.SIGKILL)
        self.assertEqual(process.communicate.call_count, 2)
        self.assertEqual(captured.exception.stderr, "phase-progress")


class PySCFBackendGateTests(unittest.TestCase):
    def test_scheduler_rejects_backend_incompatible_scf_strategy(self):
        with self.assertRaisesRegex(ValueError, "requires the PySCF backend"):
            qm_scheduler.run_scf(
                SimpleNamespace(),
                backend="psi4",
                scf_strategy="density_fit",
            )
        with self.assertRaisesRegex(ValueError, "SCF strategy should be one of"):
            qm_scheduler.run_scf(
                SimpleNamespace(),
                backend="pyscf",
                scf_strategy="implicit",
            )

    def test_wavefunction_strategy_dispatch_is_explicit_and_ordered(self):
        from Xponge.qm.backends import pyscf_backend

        calls = []

        class FakeWavefunction:
            def density_fit(self):
                calls.append("density_fit")
                return self

            def newton(self):
                calls.append("newton")
                return self

        wavefunction = FakeWavefunction()

        class FakeSCF:
            @staticmethod
            def RHF(_mol):
                calls.append("rhf")
                return wavefunction

            @staticmethod
            def ROHF(_mol):
                calls.append("rohf")
                return wavefunction

            @staticmethod
            def UHF(_mol):
                calls.append("uhf")
                return wavefunction

        molecule = QMMolecule(["He"], [(0.0, 0.0, 0.0)], 0, 0)
        for strategy, expected in (
            ("direct", ["rhf"]),
            ("density_fit", ["rhf", "density_fit"]),
            ("newton", ["rhf", "newton"]),
            ("density_fit_newton", ["rhf", "density_fit", "newton"]),
        ):
            calls.clear()
            result = pyscf_backend._build_wavefunction(
                object(),
                molecule,
                FakeSCF,
                QMRunOptions(backend="pyscf", scf_strategy=strategy),
            )
            self.assertIs(result, wavefunction)
            self.assertEqual(calls, expected)

        open_shell = QMMolecule(["H"], [(0.0, 0.0, 0.0)], 0, 1)
        calls.clear()
        result = pyscf_backend._build_wavefunction(
            object(),
            open_shell,
            FakeSCF,
            QMRunOptions(
                backend="pyscf",
                reference="auto",
                scf_strategy="density_fit",
            ),
        )
        self.assertIs(result, wavefunction)
        self.assertEqual(calls, ["rohf", "density_fit"])

        calls.clear()
        pyscf_backend._build_wavefunction(
            object(),
            open_shell,
            FakeSCF,
            QMRunOptions(backend="pyscf", reference="uhf"),
        )
        self.assertEqual(calls, ["uhf"])

    def test_backend_rejects_unconverged_scf_before_hessian(self):
        from Xponge.qm.backends import pyscf_backend

        class FakeMol:
            pass

        class FakeGTO:
            @staticmethod
            def M(**_kwargs):
                return FakeMol()

        class FakeWavefunction:
            converged = False

            def run(self):
                return self

        class FakeSCF:
            @staticmethod
            def RHF(_mol):
                return FakeWavefunction()

        molecule = QMMolecule(["He"], [(0.0, 0.0, 0.0)], 0, 0)
        options = QMRunOptions(backend="pyscf", basis="sto-3g")
        with mock.patch.object(
            pyscf_backend, "require_numpy_pyscf", return_value=(np, FakeGTO, FakeSCF)
        ):
            with self.assertRaisesRegex(RuntimeError, "did not converge"):
                pyscf_backend.compute_hessian(molecule, options)


class RealPySCFHessianSmokeTests(unittest.TestCase):
    def test_real_water_sto3g_worker(self):
        python_executable = os.environ.get("MOKDA_PYSCF_PYTHON")
        if python_executable:
            response = _invoke_hessian_worker(
                _worker_request(),
                timeout_seconds=180.0,
                python_executable=python_executable,
            )
        else:
            try:
                import pyscf  # noqa: F401
            except ImportError:
                self.skipTest("PySCF is not installed and MOKDA_PYSCF_PYTHON is not configured")
            response = _execute(_worker_request())
        self.assertTrue(response["ok"])
        self.assertTrue(response["fit_report"]["setup_metadata"]["scf_converged"])
        self.assertEqual(len(response["artifact"]["cartesian_hessian_au"]), 3)


if __name__ == "__main__":
    unittest.main()
