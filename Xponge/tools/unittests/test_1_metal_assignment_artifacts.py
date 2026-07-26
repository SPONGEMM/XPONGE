"""Cross-repository contract tests for Chemcore-produced artifacts.

The production Xponge package remains independent of Chemcore.  These tests
only import a configured Chemcore build to prove that the neutral wire contract
can consume its real output without topology reconstruction.
"""

from __future__ import annotations

from dataclasses import replace
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from unittest import mock
import unittest

from Xponge.metal_assignment import (
    BaseChargeInput,
    ChargeAssignmentContract,
    ChargePolicy,
    BondedMetalParameterSpec,
    HessianArtifact,
    MetalAtomParameterSpec,
    ModelChargeArtifact,
    RespFitInput,
    ValidationError,
    assign_base_force_field,
    assign_nonbonded_metal_ions,
    assign_standard_biomolecular,
    compose_base_force_field,
    compose_bonded_fit,
    fit_resp_charges,
    metal_assignment_package_dumps,
    metal_assignment_package_loads,
    parameterize,
    project_model_charges,
    seminario_bonded_terms,
    validate_package,
)
from Xponge.metal_assignment.assigned_system import (
    ForceRealizationProtocol,
    apply_parameterization,
)


class ChemcoreArtifactContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        configured_root = os.environ.get("MOKDA_CHEMCORE_ROOT")
        if not configured_root:
            raise unittest.SkipTest("set MOKDA_CHEMCORE_ROOT to run the live Chemcore contract test")
        cls.mokda_root = Path(configured_root).resolve()
        backend_src = cls.mokda_root / "backend" / "src"
        workflow_src = cls.mokda_root / "backend" / "external" / "metal-assignment" / "src"
        test_src = cls.mokda_root / "backend" / "native" / "chemcore" / "tests"
        if not backend_src.is_dir() or not workflow_src.is_dir() or not test_src.is_dir():
            raise unittest.SkipTest(f"invalid MOKDA_CHEMCORE_ROOT: {cls.mokda_root}")
        sys.path[:0] = [str(backend_src), str(workflow_src), str(test_src)]
        cls.chemcore = importlib.import_module("spongui._chemcore")
        cls.workflow_api = importlib.import_module("metal_assignment")
        cls.workflow_xponge = cls.workflow_api.InProcessXpongeProvider()
        fixtures = importlib.import_module("test_v2_topology_projection")
        cls.structure = fixtures.embedded_fe_mmcif()

    def prepare(self, interaction_model: str):
        engine = self.chemcore.ChemcoreEngine()
        session_id = f"xponge-contract-{interaction_model}"
        state = engine.v2_create_session({
            "session_id": session_id,
            "source_format": "mmcif",
            "source_filename": f"{session_id}.cif",
            "structure_text": self.structure,
            "session_mode": "full",
        })
        try:
            detected = engine.v2_detect_metal_sites({
                "session_id": session_id,
                "max_distance": 3.0,
                "include_distance_candidates": False,
                "donor_elements": ["N", "O", "S", "P"],
            })
            candidates = list(detected.get("candidates") or [])
            metal_atom_ids = [
                int((candidate.get("metalAtom") or {}).get("atomId") or 0)
                for candidate in candidates
            ]
            coordination_edge_ids = [
                str(coordination.get("edgeId") or "")
                for candidate in candidates
                for coordination in candidate.get("confirmedCoordination") or []
                if str(coordination.get("edgeId") or "")
            ]
            model_requests = [
                {
                    "selection_id": str(candidate.get("siteId") or ""),
                    "purpose": purpose,
                    "seed_atom_ids": [
                        int((candidate.get("metalAtom") or {}).get("atomId") or 0),
                    ],
                    "selection_strategy": "seed_connected_complete_residues",
                    "cap_strategy": "none" if purpose == "site" else "hydrogen_link",
                    "hydrogen_cap_bond_length_angstrom": 1.09,
                    "spin_multiplicity": 1,
                    "spin_source": "contract-test-explicit-state",
                }
                for candidate in candidates
                for purpose in ("site", "small", "large")
            ] if interaction_model == "bonded" else []
            prepared = engine.v2_prepare_structure_models({
                "session_id": session_id,
                "graph_revision": state["revision"],
                "input_hash": hashlib.sha256(self.structure.encode()).hexdigest(),
                "base_force_field_provider": "gaff2",
                "atomic_center_provider": "metal",
                "component_class_providers": {},
                "focus_atom_ids": metal_atom_ids,
                "assignment_atomic_center_atom_ids": metal_atom_ids,
                "included_interaction_edge_ids": (
                    coordination_edge_ids if interaction_model == "bonded" else []
                ),
                "removed_interaction_edge_ids": (
                    coordination_edge_ids if interaction_model == "nonbonded_12_6" else []
                ),
                "scopes_by_atom_id": {},
                "default_scopes": ["baseForceField", "simulation"],
                "focus_atom_scopes": ["qm", "metalOverlay", "simulation"],
                "connected_residue_scopes": [
                    "baseForceField", "qm", "metalOverlay", "simulation",
                ],
                "capability_snapshot": {
                    "schema_version": 1,
                    "provider_name": "xponge",
                    "provider_version": "1.7b3-dev",
                    "provider_revision": "chemcore-contract-test",
                    "projection_schema_version": 1,
                    "embedded_metal": True,
                    "standalone_metal": True,
                    "multi_metal": False,
                    "bonded": True,
                    "nonbonded_12_6": True,
                    "base_force_field_providers": ["gaff", "gaff2", "standard_biomolecular"],
                },
                "model_requests": model_requests,
            })
        finally:
            engine.v2_close_session({"session_id": session_id})
        request = self.workflow_api.WorkflowRequest(
            schema_version=1,
            request_id=f"contract-{interaction_model}",
            workdir=str(self.mokda_root),
            session_id=session_id,
            graph_revision=state["revision"],
            input_hash=hashlib.sha256(self.structure.encode()).hexdigest(),
            interaction_model=(
                self.workflow_api.InteractionModel.BONDED
                if interaction_model == "bonded"
                else self.workflow_api.InteractionModel.NONBONDED_12_6
            ),
            base_force_field_provider="gaff2",
            water_model="tip3p",
            electronic_state=self.workflow_api.ElectronicState(2, 1),
            capability_snapshot=self.workflow_xponge.capability_snapshot(),
            target_metal_atom_ids=tuple(metal_atom_ids),
            coordination_edge_ids=tuple(coordination_edge_ids),
            parameter_source="qm_fit" if interaction_model == "bonded" else "ion_12_6",
            site_electronic_states=(
                tuple(
                    self.workflow_api.SiteElectronicState(
                        str(candidate.get("siteId") or ""),
                        (int((candidate.get("metalAtom") or {}).get("atomId") or 0),),
                        1,
                        "contract-test-explicit-state",
                    )
                    for candidate in candidates
                )
                if interaction_model == "bonded"
                else ()
            ),
        ).with_computed_hash()
        package_payload = self.workflow_xponge.prepare_package(
            prepared,
            request,
            self.mokda_root,
        )
        return metal_assignment_package_loads(
            json.dumps(package_payload, sort_keys=True, separators=(",", ":"))
        )

    def adapt(self, interaction_model: str):
        return self.prepare(interaction_model)

    @staticmethod
    def rehash_package(package, *, artifacts=None, models=None):
        prepared = package.prepared_artifacts
        if artifacts is not None:
            prepared = replace(prepared, structural_artifacts=artifacts)
        if models is not None:
            prepared = replace(prepared, derived_models=models)
        updated = replace(package, prepared_artifacts=prepared, package_hash="")
        return replace(updated, package_hash=updated.computed_hash())

    @staticmethod
    def rehash_structural_bundle(bundle, *, prepared_system):
        artifact = replace(prepared_system, artifact_hash="")
        artifact = replace(artifact, artifact_hash=artifact.computed_hash())
        updated = replace(bundle, prepared_system=artifact, bundle_hash="")
        return replace(updated, bundle_hash=updated.computed_hash())

    @staticmethod
    def rehash_model_bundle(bundle, *, models=None, sites=None, active_edges=None):
        updated = replace(
            bundle,
            models=bundle.models if models is None else models,
            sites=bundle.sites if sites is None else sites,
            active_metal_interaction_edge_ids=(
                bundle.active_metal_interaction_edge_ids if active_edges is None else active_edges
            ),
            bundle_hash="",
        )
        return replace(updated, bundle_hash=updated.computed_hash())

    def test_bonded_output_is_consumed_and_round_trips(self):
        package = self.adapt("bonded")
        validate_package(package)
        self.assertTrue(package.prepared_artifacts.derived_models.models_generated)
        self.assertEqual(
            tuple(model.purpose for model in package.prepared_artifacts.derived_models.models),
            ("site", "small", "large"),
        )
        site, small, large = package.prepared_artifacts.derived_models.models
        self.assertIsNone(site.electronic_state)
        self.assertEqual((small.electronic_state.net_charge, large.electronic_state.net_charge), (2, 2))
        self.assertEqual(
            (small.electronic_state.source, large.electronic_state.source),
            ("contract-test-explicit-state", "contract-test-explicit-state"),
        )
        metal_canonical_id = next(
            atom.canonical_atom_id for atom in package.request.topology.atoms if atom.is_metal
        )
        expected_selection_id = f"metal:Fe:{metal_canonical_id}"
        self.assertEqual(
            (small.electronic_state.selection_id, large.electronic_state.selection_id),
            (expected_selection_id, expected_selection_id),
        )
        payload = metal_assignment_package_dumps(package)
        self.assertEqual(metal_assignment_package_loads(payload), package)
        self.assertEqual(metal_assignment_package_dumps(metal_assignment_package_loads(payload)), payload)

    def test_nonbonded_parameterize_operation_is_closed(self):
        nonbonded = self.adapt("nonbonded_12_6")
        result = parameterize(
            nonbonded,
            operation="assign-nonbonded-complete",
            provider="gaff2",
            water_model="tip3p",
        )
        self.assertEqual(result.status, "overlay_validated")
        self.assertEqual(
            result.topology_hash,
            nonbonded.request.topology.topology_hash,
        )

    def test_nonbonded_output_proves_model_generation_was_skipped(self):
        package = self.adapt("nonbonded_12_6")
        models = package.prepared_artifacts.derived_models
        self.assertFalse(models.models_generated)
        self.assertEqual(models.skipped_reason, "nonbonded_12_6_fast_path")
        self.assertFalse(models.models)
        self.assertFalse(models.mapping.records)
        ion_output = assign_nonbonded_metal_ions(package, water_model="tip3p")
        metal_ids = tuple(atom.external_id for atom in package.request.topology.atoms if atom.is_metal)
        self.assertEqual(ion_output.overlay.covered_atom_ids, metal_ids)
        self.assertEqual(set(ion_output.overlay.atom_types), set(metal_ids))
        self.assertEqual(set(ion_output.overlay.charges), set(metal_ids))
        self.assertEqual(set(ion_output.overlay.masses), set(metal_ids))
        self.assertEqual(set(ion_output.overlay.lj_parameters), set(metal_ids))
        self.assertFalse(ion_output.overlay.bonded_parameters)
        self.assertEqual(ion_output.report.report_hash, ion_output.report.computed_hash())

    def test_native_gaff2_assignment_consumes_locked_chemcore_component(self):
        package = self.adapt("bonded")
        output = assign_base_force_field(package, "gaff2", timeout_seconds=60.0)
        expected_components = tuple(
            component
            for component in package.request.assignment_components
            if component.base_force_field and component.provider == "gaff2"
        )
        expected_atom_ids = tuple(
            atom_id for component in expected_components for atom_id in component.atom_ids
        )
        self.assertEqual(output.overlay.covered_atom_ids, expected_atom_ids)
        self.assertEqual(set(output.overlay.atom_types), set(expected_atom_ids))
        self.assertEqual(set(output.overlay.masses), set(expected_atom_ids))
        self.assertEqual(set(output.overlay.lj_parameters), set(expected_atom_ids))
        self.assertFalse(output.overlay.charges)
        atom_by_id = {atom.external_id: atom for atom in package.request.topology.atoms}
        self.assertFalse(any(atom_by_id[atom_id].is_metal for atom_id in output.overlay.covered_atom_ids))
        self.assertEqual(output.report.component_ids, tuple(component.external_id for component in expected_components))
        self.assertEqual(output.report.report_hash, output.report.computed_hash())

    def test_native_gaff2_assignment_can_explicitly_assign_component_charges(self):
        package = self.adapt("bonded")
        charge_input = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="xponge-component-charge-test",
        ).with_computed_hash()
        output = assign_base_force_field(
            package,
            "gaff2",
            timeout_seconds=60.0,
            base_charge_input=charge_input,
        )
        expected_components = tuple(
            component
            for component in package.request.assignment_components
            if component.base_force_field and component.provider == "gaff2"
        )
        expected_atom_ids = {
            atom_id for component in expected_components for atom_id in component.atom_ids
        }
        self.assertEqual(set(output.overlay.charges), expected_atom_ids)
        for component in expected_components:
            self.assertAlmostEqual(
                sum(output.overlay.charges[atom_id] for atom_id in component.atom_ids),
                component.net_formal_charge,
                places=6,
            )
        self.assertEqual(output.report.charge_method, "tpacm4")
        self.assertEqual(output.report.charge_source, charge_input.source)
        self.assertEqual(output.report.charge_input_hash, charge_input.charge_input_hash)
    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_3GOU") and os.environ.get("MOKDA_METAL_FIXTURE_4EWL"),
        "set repaired real-fixture paths to run native GAFF2 integration tests",
    )
    def test_native_gaff2_assignment_real_3gou_and_repaired_4ewl(self):
        charge_input = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="real-fixture-empirical-component-charge",
        ).with_computed_hash()
        cases = (
            (Path(os.environ["MOKDA_METAL_FIXTURE_3GOU"]), 4, 288, {-4}),
            (Path(os.environ["MOKDA_METAL_FIXTURE_4EWL"]), 13, 256, {-1, 0}),
        )
        for fixture, expected_component_count, expected_atom_count, expected_charges in cases:
            with self.subTest(fixture=fixture.name):
                self.structure = fixture.read_text()
                package = self.adapt("bonded")
                output = assign_base_force_field(
                    package,
                    "gaff2",
                    timeout_seconds=120.0,
                    base_charge_input=charge_input,
                )
                components = tuple(
                    component
                    for component in package.request.assignment_components
                    if component.base_force_field and component.provider == "gaff2"
                )
                self.assertEqual(len(components), expected_component_count)
                self.assertEqual(output.report.atom_count, expected_atom_count)
                self.assertEqual({component.net_formal_charge for component in components}, expected_charges)
                self.assertTrue(output.report.bonded_term_count)
                self.assertEqual(len(output.report.frcmod_hashes), expected_component_count)
                self.assertTrue(all(output.report.frcmod_hashes.values()))
                self.assertEqual(set(output.overlay.charges), set(output.overlay.covered_atom_ids))
                for component in components:
                    self.assertAlmostEqual(
                        sum(output.overlay.charges[atom_id] for atom_id in component.atom_ids),
                        component.net_formal_charge,
                        places=6,
                    )
                atom_by_id = {atom.external_id: atom for atom in package.request.topology.atoms}
                self.assertFalse(any(atom_by_id[atom_id].is_metal for atom_id in output.overlay.covered_atom_ids))

    def _assert_real_nonbonded_component_charges_apply_and_export(
        self,
        fixture: Path,
        expected_atom_count: int,
        formats: tuple[str, ...],
    ) -> None:
        charge_input = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="real-fixture-empirical-component-charge",
        ).with_computed_hash()
        protocol = ForceRealizationProtocol(
            family="amber",
            exclusion_bond_depth=3,
            one_four_lj_scale=0.5,
            one_four_electrostatic_scale=0.8333333333333334,
            source="amber-default",
        ).with_computed_hash()
        with TemporaryDirectory(prefix="xponge-real-nonbonded-export-") as tempdir:
            self.structure = fixture.read_text()
            package = self.adapt("nonbonded_12_6")
            result = parameterize(
                package,
                operation="assign-nonbonded-complete",
                provider="gaff2",
                water_model="tip3p",
                base_charge_input=charge_input,
                timeout_seconds=180.0,
            )
            applied = apply_parameterization(
                package.request,
                result,
                protocol,
                timeout_seconds=180.0,
            )
            self.assertTrue(applied.assigned_system.complete)
            self.assertEqual(len(applied.assigned_system.atom_parameters), expected_atom_count)
            self.assertEqual(
                applied.assigned_system.topology.topology_hash,
                package.request.topology.topology_hash,
            )
            self.assertLess(
                applied.assigned_system.application_audit["materialized_lj_type_count"],
                expected_atom_count,
            )
            from Xponge import save_sponge_input
            from Xponge.metal_assignment._apply_worker import materialize_assigned_system

            molecule, atoms_by_external_id, _ = materialize_assigned_system(
                applied.assigned_system
            )
            molecule.box_length = [500.0, 500.0, 500.0]
            molecule.box_angle = [90.0, 90.0, 90.0]
            source_ids = {
                atom: external_id
                for external_id, atom in atoms_by_external_id.items()
            }
            mappings = []
            for format_name in formats:
                output_root = Path(tempdir) / format_name
                output_root.mkdir(parents=True, exist_ok=True)
                _, mapping = save_sponge_input(
                    molecule,
                    "system",
                    str(output_root),
                    format=format_name,
                    source_atom_ids=source_ids,
                    return_mapping=True,
                )
                mappings.append(mapping)
                self.assertTrue(any(output_root.iterdir()))
            self.assertTrue(mappings)
            self.assertTrue(all(mapping == mappings[0] for mapping in mappings))

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_3GOU"),
        "set repaired 3GOU fixture path to run nonbonded apply/export integration",
    )
    def test_real_3gou_nonbonded_component_charges_apply_and_export_raw(self):
        self._assert_real_nonbonded_component_charges_apply_and_export(
            Path(os.environ["MOKDA_METAL_FIXTURE_3GOU"]),
            9206,
            ("raw",),
        )

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_3GOU"),
        "set repaired 3GOU fixture path to run nonbonded apply/export integration",
    )
    def test_real_3gou_nonbonded_component_charges_apply_and_export_bundle(self):
        self._assert_real_nonbonded_component_charges_apply_and_export(
            Path(os.environ["MOKDA_METAL_FIXTURE_3GOU"]),
            9206,
            ("bundle",),
        )

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_4EWL"),
        "set repaired 4EWL fixture path to run nonbonded apply/export integration",
    )
    def test_real_4ewl_nonbonded_component_charges_apply_and_export_raw(self):
        self._assert_real_nonbonded_component_charges_apply_and_export(
            Path(os.environ["MOKDA_METAL_FIXTURE_4EWL"]),
            10652,
            ("raw",),
        )

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_4EWL"),
        "set repaired 4EWL fixture path to run nonbonded apply/export integration",
    )
    def test_real_4ewl_nonbonded_component_charges_apply_and_export_bundle(self):
        self._assert_real_nonbonded_component_charges_apply_and_export(
            Path(os.environ["MOKDA_METAL_FIXTURE_4EWL"]),
            10652,
            ("bundle",),
        )

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_3GOU") and os.environ.get("MOKDA_METAL_FIXTURE_4EWL"),
        "set repaired real-fixture paths to run complete base-assignment integration tests",
    )
    def test_real_bonded_empirical_component_charges_apply_for_3gou_and_repaired_4ewl(self):
        charge_input = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="real-fixture-empirical-component-charge",
        ).with_computed_hash()
        cases = (
            (Path(os.environ["MOKDA_METAL_FIXTURE_3GOU"]), 604, 8914),
            (Path(os.environ["MOKDA_METAL_FIXTURE_4EWL"]), 1154, 10394),
        )
        for fixture, expected_standard_components, expected_standard_atoms in cases:
            with self.subTest(fixture=fixture.name):
                self.structure = fixture.read_text()
                package = self.adapt("bonded")
                gaff = assign_base_force_field(
                    package,
                    "gaff2",
                    timeout_seconds=120.0,
                    base_charge_input=charge_input,
                )
                standard = assign_standard_biomolecular(package, timeout_seconds=120.0)
                complete = compose_base_force_field(package, (standard, gaff))
                expected_base_ids = {
                    atom_id
                    for component in package.request.assignment_components
                    if component.base_force_field
                    for atom_id in component.atom_ids
                }
                self.assertEqual(len(standard.report.component_ids), expected_standard_components)
                self.assertEqual(standard.report.atom_count, expected_standard_atoms)
                self.assertEqual(set(complete.overlay.covered_atom_ids), expected_base_ids)
                self.assertEqual(set(complete.overlay.atom_types), expected_base_ids)
                self.assertEqual(set(complete.overlay.charges), expected_base_ids)
                self.assertEqual(set(complete.overlay.masses), expected_base_ids)
                self.assertEqual(set(complete.overlay.lj_parameters), expected_base_ids)
                self.assertEqual(
                    complete.report.providers,
                    ("gaff2", "standard_biomolecular"),
                )
                self.assertEqual(complete.report.report_hash, complete.report.computed_hash())
                metal_atoms = tuple(atom for atom in package.request.topology.atoms if atom.is_metal)
                specs = []
                for atom in metal_atoms:
                    if atom.element == "Fe":
                        atom_type, mass, epsilon, rmin = "Fe2+", 55.85, 0.00941798, 1.353
                    else:
                        atom_type, mass, epsilon, rmin = "Zn2+", 65.4, 0.00330286, 1.271
                    specs.append(MetalAtomParameterSpec(
                        atom.external_id,
                        atom.element,
                        atom.formal_charge,
                        atom_type,
                        float(atom.formal_charge),
                        mass,
                        epsilon,
                        rmin,
                    ))
                metal_spec = BondedMetalParameterSpec(
                    schema_version=package.request.schema_version,
                    topology_hash=package.request.topology.topology_hash,
                    metal_atoms=tuple(specs),
                    donor_atom_types={},
                    parameter_source="unit-test:explicit-metal-nonbonded-v1",
                    provenance={"fixture": fixture.name, "purpose": "deterministic-fit-composition"},
                ).with_computed_hash()
                if fixture.name.startswith("3GOU"):
                    import numpy as np
                    hessians = tuple(
                        HessianArtifact(
                            model.external_id,
                            model.model_hash,
                            tuple(atom.model_atom_id for atom in model.atoms),
                            tuple(tuple(atom.coordinates) for atom in model.atoms),
                            np.eye(len(model.atoms) * 3, dtype=float) * 0.01,
                            "deterministic-test",
                            "1",
                        )
                        for model in package.prepared_artifacts.derived_models.models
                        if model.purpose == "small"
                    )
                    fitted = compose_bonded_fit(
                        package,
                        complete,
                        metal_spec,
                        hessian_artifacts=hessians,
                        force_method="seminario",
                    )
                else:
                    fitted = compose_bonded_fit(
                        package,
                        complete,
                        metal_spec,
                        force_method="empirical_zn_nos",
                    )
                self.assertEqual(fitted.status, "overlay_validated")
                self.assertFalse(fitted.complete)
                self.assertTrue(fitted.bonded_overlay.terms)
                self.assertTrue(all(
                    any(atom_id in {atom.external_id for atom in metal_atoms} for atom_id in term["atom_ids"])
                    for term in fitted.bonded_overlay.terms.values()
                ))
                protocol = ForceRealizationProtocol(
                    "amber", 3, 0.5, 0.833333, "xponge:amber-explicit-v1",
                ).with_computed_hash()
                applied = apply_parameterization(
                    package.request,
                    fitted,
                    protocol,
                    timeout_seconds=180.0,
                )
                self.assertTrue(applied.assigned_system.complete)
                self.assertEqual(
                    applied.assigned_system.topology.topology_hash,
                    package.request.topology.topology_hash,
                )
                self.assertEqual(len(applied.assigned_system.atom_parameters), len(package.request.topology.atoms))
                self.assertLess(
                    applied.assigned_system.application_audit["materialized_lj_type_count"],
                    len(package.request.topology.atoms),
                )

    @unittest.skipUnless(
        os.environ.get("MOKDA_METAL_FIXTURE_3GOU") and os.environ.get("MOKDA_METAL_FIXTURE_4EWL"),
        "set repaired real-fixture paths to run native ion integration tests",
    )
    def test_native_nonbonded_ion_assignment_real_3gou_and_repaired_4ewl(self):
        for fixture, element, expected_count in (
            (Path(os.environ["MOKDA_METAL_FIXTURE_3GOU"]), "Fe", 4),
            (Path(os.environ["MOKDA_METAL_FIXTURE_4EWL"]), "Zn", 2),
        ):
            with self.subTest(fixture=fixture.name):
                self.structure = fixture.read_text()
                package = self.adapt("nonbonded_12_6")
                output = assign_nonbonded_metal_ions(package, water_model="tip3p")
                metal_atoms = tuple(atom for atom in package.request.topology.atoms if atom.is_metal)
                self.assertEqual(len(metal_atoms), expected_count)
                self.assertEqual({atom.element for atom in metal_atoms}, {element})
                self.assertEqual(set(output.overlay.covered_atom_ids), {atom.external_id for atom in metal_atoms})
                self.assertTrue(all(output.overlay.masses[atom.external_id] > 0 for atom in metal_atoms))

    def test_tampered_structural_artifact_is_rejected(self):
        package = self.adapt("bonded")
        artifacts = package.prepared_artifacts.structural_artifacts
        corrupted_prepared = replace(
            artifacts.prepared_system,
            mol2_text=artifacts.prepared_system.mol2_text.replace("NO_CHARGES", "USER_CHARGES", 1),
        )
        corrupted_artifacts = replace(artifacts, prepared_system=corrupted_prepared)
        corrupted_package = replace(
            package,
            prepared_artifacts=replace(package.prepared_artifacts, structural_artifacts=corrupted_artifacts),
        )
        with self.assertRaisesRegex(ValidationError, "stale_structural_bundle_hash|stale_artifact_hash"):
            validate_package(corrupted_package)

    def test_hash_closed_model_edge_with_unknown_endpoint_is_rejected(self):
        package = self.adapt("bonded")
        bundle = package.prepared_artifacts.derived_models
        model = bundle.models[0]
        bond = model.bonds[0]
        corrupted_bond = replace(bond, model_atom_ids=("missing-model-atom", bond.model_atom_ids[1]))
        corrupted_model = replace(model, bonds=(corrupted_bond, *model.bonds[1:]), model_hash="")
        corrupted_model = replace(corrupted_model, model_hash=corrupted_model.computed_hash())
        corrupted_bundle = self.rehash_model_bundle(
            bundle,
            models=(corrupted_model, *bundle.models[1:]),
        )
        corrupted_package = self.rehash_package(package, models=corrupted_bundle)
        with self.assertRaisesRegex(ValidationError, "model_edge_unknown_atom"):
            validate_package(corrupted_package)

    def test_hash_closed_site_with_unknown_metal_is_rejected(self):
        package = self.adapt("bonded")
        bundle = package.prepared_artifacts.derived_models
        corrupted_site = replace(bundle.sites[0], metal_atom_ids=("missing-metal",))
        corrupted_bundle = self.rehash_model_bundle(
            bundle,
            sites=(corrupted_site, *bundle.sites[1:]),
        )
        corrupted_package = self.rehash_package(package, models=corrupted_bundle)
        with self.assertRaisesRegex(ValidationError, "invalid_site_metal_membership"):
            validate_package(corrupted_package)

    def test_nonbonded_hash_closed_active_edge_mismatch_is_rejected(self):
        package = self.adapt("nonbonded_12_6")
        bundle = package.prepared_artifacts.derived_models
        self.assertTrue(bundle.removed_interaction_edge_ids)
        corrupted_bundle = self.rehash_model_bundle(
            bundle,
            active_edges=(bundle.removed_interaction_edge_ids[0],),
        )
        corrupted_package = self.rehash_package(package, models=corrupted_bundle)
        with self.assertRaisesRegex(ValidationError, "active_metal_interaction_edge_mismatch"):
            validate_package(corrupted_package)

    def test_hash_closed_mol2_without_atom_section_is_rejected(self):
        package = self.adapt("bonded")
        bundle = package.prepared_artifacts.structural_artifacts
        corrupted_artifact = replace(
            bundle.prepared_system,
            mol2_text=bundle.prepared_system.mol2_text.replace("@<TRIPOS>ATOM", "@<TRIPOS>ATXM", 1),
        )
        corrupted_bundle = self.rehash_structural_bundle(
            bundle,
            prepared_system=corrupted_artifact,
        )
        corrupted_package = self.rehash_package(package, artifacts=corrupted_bundle)
        with self.assertRaisesRegex(ValidationError, "invalid_mol2"):
            validate_package(corrupted_package)

    def test_seminario_provider_consumes_explicit_small_model_hessian(self):
        package = self.adapt("bonded")
        small = next(
            model for model in package.prepared_artifacts.derived_models.models if model.purpose == "small"
        )
        atom_count = len(small.atoms)
        hessian = tuple(
            tuple(1.0 if row == column else 0.0 for column in range(atom_count * 3))
            for row in range(atom_count * 3)
        )
        artifact = HessianArtifact(
            model_id=small.external_id,
            model_hash=small.model_hash,
            atom_order=tuple(atom.model_atom_id for atom in small.atoms),
            coordinates_angstrom=tuple(atom.coordinates for atom in small.atoms),
            cartesian_hessian_au=hessian,
            provider="deterministic-test",
            provider_version="1",
        )
        terms, report = seminario_bonded_terms(small, artifact)
        self.assertTrue(terms)
        self.assertTrue(all("cap:" not in atom_id for term in terms.values() for atom_id in term["atom_ids"]))
        self.assertEqual(report["model_hash"], small.model_hash)

    def test_charge_provider_projects_to_parent_and_enforces_site_target(self):
        package = self.adapt("bonded")
        topology = package.request.topology
        large = next(
            model for model in package.prepared_artifacts.derived_models.models if model.purpose == "large"
        )
        fit_atom_ids = tuple(str(atom.external_id) for atom in large.atoms if atom.role != "cap")
        component_targets = {
            component.external_id: component.net_formal_charge
            for component in topology.components
            if set(component.atom_ids) & set(fit_atom_ids)
        }
        contract = ChargeAssignmentContract(
            policy=ChargePolicy.SITE_JOINT_FIT,
            fit_atom_ids=fit_atom_ids,
            component_formal_charges=component_targets,
            site_target=2.0,
        )
        model_charges = {
            atom.model_atom_id: 0.0
            for atom in large.atoms
        }
        for group in large.charge_accounting["constraint_groups"]:
            model_charges[group["model_atom_ids"][0]] = float(group["target_charge"])
        artifact = ModelChargeArtifact(
            model_id=large.external_id,
            model_hash=large.model_hash,
            atom_order=tuple(atom.model_atom_id for atom in large.atoms),
            charges=tuple(model_charges[atom.model_atom_id] for atom in large.atoms),
            atomic_charge_role="fitted",
            provider="deterministic-test",
            provider_version="1",
        ).with_computed_hash()
        overlay, report = project_model_charges(topology, large, artifact, contract)
        self.assertAlmostEqual(sum(overlay.charges.values()), 2.0, places=8)
        self.assertEqual(set(overlay.charges), set(fit_atom_ids))
        self.assertEqual(report["cap_atom_count"], 0)
        self.assertTrue(report["constrained_model_charges"])
        self.assertAlmostEqual(report["cap_charge_projected"], 0.0, places=8)
        self.assertAlmostEqual(report["cap_charge_discarded"], 0.0, places=8)
        self.assertLessEqual(report["max_model_constraint_residual"], 1.0e-8)
        self.assertEqual(overlay.overlay_hash, overlay.computed_hash())
        self.assertEqual(tuple(overlay.artifact_hashes), (report["projected_artifact_hash"],))

    def test_isolated_resp_provider_consumes_each_large_model_and_executor_closes_output(self):
        package = self.adapt("bonded")
        topology = package.request.topology
        large_models = tuple(
            model for model in package.prepared_artifacts.derived_models.models if model.purpose == "large"
        )
        fit_atom_ids = tuple(
            str(atom.external_id)
            for model in large_models
            for atom in model.atoms
            if atom.role != "cap"
        )
        component_targets = {
            component.external_id: component.net_formal_charge
            for component in topology.components
            if set(component.atom_ids) & set(fit_atom_ids)
        }
        contract = ChargeAssignmentContract(
            policy=ChargePolicy.SITE_JOINT_FIT,
            fit_atom_ids=fit_atom_ids,
            component_formal_charges=component_targets,
            site_target=2.0,
            source="live-chemcore-resp-contract",
        )
        request = replace(package.request, charge_contract=contract, projection_hash="").with_computed_hash()
        package = replace(package, request=request, package_hash="")
        package = replace(package, package_hash=package.computed_hash())
        fit_input = RespFitInput(
            schema_version=2,
            backend="pyscf",
            basis_family="SDD",
            metal_basis_policy="require_ecp",
            basis_source="xponge:sdd-stuttgart-dz-v1",
            optimize_geometry=False,
            grid_density=1.0,
            grid_cell_layer=1,
            radius_overrides={},
            restraint_a1=0.0005,
            restraint_a2=0.001,
            two_stage=True,
            only_esp=False,
            esp_memory_limit_bytes=64 * 1024 * 1024,
            esp_chunk_policy="pointwise",
            esp_safety_factor=0.8,
            equivalence_groups=(),
            source="live-chemcore-deterministic-resp",
        ).with_computed_hash()

        def fake_resp(assign, **kwargs):
            from Xponge.assign.resp import get_resp_setup_metadata
            from Xponge.assign.resp_core import (
                _prepare_linear_constraints,
                _solve_constrained_quadratic,
            )
            import numpy as np

            charge = float(kwargs["charge"])
            metadata = get_resp_setup_metadata(assign, kwargs["basis"])
            metadata.update({
                "backend": kwargs["backend"],
                "scf_converged": True,
                "total_energy_hartree": -1.0,
            })
            constraints, targets, constraint_report = _prepare_linear_constraints(
                assign.atom_numbers,
                charge,
                extra_equivalence=kwargs["extra_equivalence"],
                constraint_matrix=kwargs["constraint_matrix"],
                constraint_targets=kwargs["constraint_targets"],
            )
            charges, solve_report = _solve_constrained_quadratic(
                np.eye(assign.atom_numbers),
                np.zeros(assign.atom_numbers),
                constraints,
                targets,
            )
            return {
                "charges": charges,
                "metadata": metadata,
                "diagnostics": {
                    **constraint_report,
                    **solve_report,
                    "esp_point_count": 8,
                    "esp_rmse_au": 0.0,
                    "esp_relative_rmse": 0.0,
                    "esp_mae_au": 0.0,
                    "esp_max_abs_error_au": 0.0,
                },
            }

        def invoke(payload, **_kwargs):
            from Xponge.metal_assignment._resp_worker import _execute

            with mock.patch("Xponge.assign.resp.resp_fit", side_effect=fake_resp):
                return _execute(payload)

        with mock.patch("Xponge.metal_assignment.resp_provider._invoke_resp_worker", side_effect=invoke):
            output = fit_resp_charges(package, fit_input)
        self.assertEqual(len(output.model_artifacts), len(large_models))
        self.assertTrue(output.output_hash)


if __name__ == "__main__":
    unittest.main()
