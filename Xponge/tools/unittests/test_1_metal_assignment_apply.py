"""Transactional assigned-system materialization tests."""

from dataclasses import replace
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from Xponge import save_sponge_input
from Xponge.helper import AtomType, Molecule, ResidueType
from Xponge.metal_assignment import (
    AssignmentComponent,
    BaseForceFieldOverlay,
    BondedParameterOverlay,
    ChemicalTopologyProof,
    MetalParameterOverlay,
    ParameterizationResult,
    ProviderProjection,
    ValidationError,
    compose_nonbonded_assignment,
    apply as apply_to_molecule,
)
from Xponge.metal_assignment.apply_input import (
    ParameterApplyInput,
    parameter_apply_input_dumps,
    parameter_apply_input_loads,
)
from Xponge.metal_assignment.assigned_system import (
    ForceRealizationProtocol,
    apply_parameterization,
    assigned_system_dumps,
    assigned_system_loads,
)
from Xponge.metal_assignment._apply_worker import materialize_assigned_system
from Xponge.metal_assignment.base_composition import BaseCompositionOutput, BaseCompositionReport
from Xponge.metal_assignment.metal_overlay import MetalOverlayBuildOutput, MetalOverlayBuildReport
from Xponge.tools.unittests.test_1_metal_assignment_contracts import (
    _bonded_request,
    _nonbonded_request,
    _overlay_result,
)


def _complete_request():
    request = _bonded_request()
    topology = request.topology
    atoms = tuple(
        replace(atom, scopes=("baseForceField", "qm", "simulation"))
        if atom.external_id == "his-ne2" else atom
        for atom in topology.atoms
    )
    topology = replace(topology, atoms=atoms, topology_hash="").with_computed_hash()
    proof = ChemicalTopologyProof(
        1, request.graph_revision, "fixture-perception-v1", ("his-ne2",),
        {"his-ne2": 0}, {"his-ne2": (0,)}, {"his-ne2": 3.0},
        {"his-ne2": "valid_unique"}, (), "complete",
        "component_charge_constrained_valence", True,
    ).with_computed_hash()
    component = AssignmentComponent(
        "protein-component", ("his-ne2",), "polymer_donor_fragment", "standard_biomolecular",
        0, "atom_formal_charge_sum", True, proof,
    )
    projection = ProviderProjection(
        "protein-component", "protein-component", "HIS:A:2", "protein-component",
        "protein-component", "HIS:A:2", ("his-ne2",),
        {"his-ne2": ("baseForceField", "qm", "simulation")},
        parent_component_id="protein-component",
    )
    return replace(
        request,
        topology=topology,
        projections=(*request.projections, projection),
        assignment_components=(*request.assignment_components, component),
        capability_snapshot=replace(
            request.capability_snapshot,
            base_force_field_providers=("gaff", "gaff2", "standard_biomolecular"),
        ),
        projection_hash="",
    ).with_computed_hash()


def _lj(source, epsilon, rmin):
    return {
        "epsilon": epsilon,
        "rmin": rmin,
        "energy_unit": "kcal/mol",
        "length_unit": "angstrom",
        "source": source,
    }


def _complete_result(request):
    topology_hash = request.topology.topology_hash
    base_ids = ("heme-na", "heme-c", "his-ne2")
    base = BaseForceFieldOverlay(
        topology_hash=topology_hash,
        covered_atom_ids=base_ids,
        atom_types={"heme-na": "n", "heme-c": "c", "his-ne2": "NA"},
        charges={"heme-na": -0.4, "heme-c": 0.4, "his-ne2": 0.0},
        masses={"heme-na": 14.01, "heme-c": 12.01, "his-ne2": 14.01},
        lj_parameters={
            "heme-na": _lj("fixture-base", 0.17, 1.824),
            "heme-c": _lj("fixture-base", 0.086, 1.908),
            "his-ne2": _lj("fixture-base", 0.17, 1.824),
        },
        bonded_parameters={
            "base:na-c": {
                "kind": "bond", "atom_ids": ("heme-na", "heme-c"),
                "parameters": {"k": 300.0, "equilibrium": 1.3}, "source": "fixture-base",
            },
        },
        parameter_source="fixture-base",
    )
    metal = MetalParameterOverlay(
        topology_hash=topology_hash,
        covered_atom_ids=("heme-fe",),
        atom_types={"heme-fe": "FE-site"},
        charges={"heme-fe": 2.0},
        masses={"heme-fe": 55.845},
        lj_parameters={"heme-fe": _lj("fixture-metal", 0.01, 1.2)},
        parameter_source="fixture-metal",
    )
    bonded = BondedParameterOverlay(
        topology_hash,
        {
            "metal:fe-na": {
                "kind": "bond", "atom_ids": ("heme-fe", "heme-na"),
                "parameters": {"k": 100.0, "equilibrium": 1.9}, "source": "fixture-fit",
            },
            "metal:fe-his": {
                "kind": "bond", "atom_ids": ("heme-fe", "his-ne2"),
                "parameters": {"k": 80.0, "equilibrium": 2.1}, "source": "fixture-fit",
            },
        },
        "fixture-fit",
    )
    return ParameterizationResult(
        schema_version=1,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology_hash=topology_hash,
        projection_hash=request.projection_hash,
        base_overlay=base,
        metal_overlay=metal,
        bonded_overlay=bonded,
        fit_reports={"base": {"provider": "fixture"}, "fit": {"provider": "fixture"}},
        provenance={"provider_revision": "fixture", "protocol": "fixture"},
        status="overlay_validated",
    ).with_computed_hash()


def _force_protocol():
    return ForceRealizationProtocol(
        family="amber",
        exclusion_bond_depth=3,
        one_four_lj_scale=0.5,
        one_four_electrostatic_scale=0.833333,
        source="xponge:amber-explicit-v1",
    ).with_computed_hash()


class MetalAssignmentApplyTests(unittest.TestCase):
    def test_molecule_api_applies_only_local_overlay_and_preserves_input(self):
        request = _complete_request()
        result = _complete_result(request)
        artifact = apply_parameterization(request, result, _force_protocol()).assigned_system
        molecule, atoms_by_external_id, _ = materialize_assigned_system(artifact)
        molecule.get_atoms()
        mapping = {
            external_id: molecule.atom_index[atom]
            for external_id, atom in atoms_by_external_id.items()
        }
        molecule.box_length = [50.0, 50.0, 50.0]
        molecule.box_angle = [90.0, 90.0, 90.0]
        for index, atom in enumerate(molecule.atoms):
            if atom.x is None:
                atom.x = float(index)
            if atom.y is None:
                atom.y = 0.0
            if atom.z is None:
                atom.z = 0.0
        with TemporaryDirectory(prefix="metal-assignment-checkpoint-save-") as tempdir:
            _, checkpoint_mapping = save_sponge_input(
                molecule,
                "system",
                tempdir,
                source_atom_ids={
                    atom: external_id
                    for external_id, atom in atoms_by_external_id.items()
                },
                return_mapping=True,
            )
        self.assertEqual(
            {row["source_atom_id"] for row in checkpoint_mapping},
            set(mapping),
        )
        atoms_by_external_id["heme-c"].charge = 0.123

        output = apply_to_molecule(molecule, request, result, mapping)

        self.assertIsNot(output.molecule, molecule)
        self.assertAlmostEqual(atoms_by_external_id["heme-c"].charge, 0.123)
        assigned_heme_c = output.molecule.atoms[mapping["heme-c"]]
        assigned_fe = output.molecule.atoms[mapping["heme-fe"]]
        self.assertAlmostEqual(assigned_heme_c.charge, 0.123)
        self.assertAlmostEqual(assigned_fe.charge, 2.0)
        self.assertTrue(output.application_report["topology_preserved"])
        self.assertTrue(output.application_report["ordinary_force_field_preserved"])
        self.assertEqual(
            output.application_report["force_counts"],
            {"bond": 2},
        )
        output.molecule.box_length = [50.0, 50.0, 50.0]
        output.molecule.box_angle = [90.0, 90.0, 90.0]
        for index, atom in enumerate(output.molecule.atoms):
            if atom.x is None:
                atom.x = float(index)
            if atom.y is None:
                atom.y = 0.0
            if atom.z is None:
                atom.z = 0.0
        source_ids = {
            output.molecule.atoms[index]: external_id
            for external_id, index in output.atom_mapping
        }
        with TemporaryDirectory(prefix="metal-assignment-native-save-") as tempdir:
            _, saved_mapping = save_sponge_input(
                output.molecule,
                "system",
                tempdir,
                source_atom_ids=source_ids,
                return_mapping=True,
            )
        self.assertEqual(
            {row["source_atom_id"] for row in saved_mapping},
            set(mapping),
        )

    def test_molecule_api_accepts_site_mapping_inside_larger_molecule(self):
        request = _complete_request()
        result = _complete_result(request)
        result = replace(
            result,
            provenance={**result.provenance, "test_case": "site-inside-larger-molecule"},
            result_hash="",
        ).with_computed_hash()
        artifact = apply_parameterization(request, result, _force_protocol()).assigned_system
        molecule, atoms_by_external_id, _ = materialize_assigned_system(artifact)
        molecule.get_atoms()
        mapping = {
            external_id: molecule.atom_index[atom]
            for external_id, atom in atoms_by_external_id.items()
        }
        original_count = len(molecule.atoms)
        extra_type = ResidueType(name=f"EXTRA_{result.result_hash[:12]}")
        extra_type.add_atom(
            "X",
            molecule.atoms[0].type,
            x=100.0,
            y=100.0,
            z=100.0,
        )
        molecule.add_residue(extra_type)
        molecule.get_atoms()
        molecule.built = True

        output = apply_to_molecule(molecule, request, result, mapping)

        self.assertGreater(len(output.molecule.atoms), len(request.topology.atoms))
        self.assertEqual(len(output.molecule.atoms), len(molecule.atoms))
        self.assertGreater(len(output.molecule.atoms), original_count)
        self.assertEqual(
            dict(output.atom_mapping),
            mapping,
        )

    def test_molecule_api_inplace_commit_is_transactional(self):
        request = _complete_request()
        result = _complete_result(request)
        result = replace(
            result,
            provenance={**result.provenance, "test_case": "transactional-inplace"},
            result_hash="",
        ).with_computed_hash()
        artifact = apply_parameterization(request, result, _force_protocol()).assigned_system
        molecule, atoms_by_external_id, _ = materialize_assigned_system(artifact)
        molecule.get_atoms()
        mapping = {
            external_id: molecule.atom_index[atom]
            for external_id, atom in atoms_by_external_id.items()
        }
        atom_types_before = tuple(atom.type for atom in molecule.atoms)
        charges_before = tuple(atom.charge for atom in molecule.atoms)
        registry_before = tuple(AtomType.get_all_types())

        from Xponge.metal_assignment import molecule_api

        real_get_or_create = molecule_api._get_or_create_type
        calls = 0

        def fail_after_first_registration(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise RuntimeError("synthetic apply failure")
            return real_get_or_create(*args, **kwargs)

        with patch.object(
            molecule_api,
            "_get_or_create_type",
            side_effect=fail_after_first_registration,
        ):
            with self.assertRaisesRegex(RuntimeError, "synthetic apply failure"):
                apply_to_molecule(
                    molecule,
                    request,
                    result,
                    mapping,
                    inplace=True,
                )

        self.assertEqual(tuple(atom.type for atom in molecule.atoms), atom_types_before)
        self.assertEqual(tuple(atom.charge for atom in molecule.atoms), charges_before)
        self.assertEqual(tuple(AtomType.get_all_types()), registry_before)

    def test_native_save_returns_source_mapping_for_raw_and_bundle(self):
        request = _complete_request()
        result = _complete_result(request)
        result = replace(
            result,
            provenance={**result.provenance, "test_case": "native-save-mapping"},
            result_hash="",
        ).with_computed_hash()
        artifact = apply_parameterization(request, result, _force_protocol()).assigned_system
        molecule, atoms_by_external_id, _ = materialize_assigned_system(artifact)
        molecule.get_atoms()
        source_ids = {
            atom: external_id
            for external_id, atom in atoms_by_external_id.items()
        }
        with TemporaryDirectory(prefix="xponge-native-save-") as tempdir:
            returned_raw, raw_mapping = save_sponge_input(
                molecule,
                "fixture",
                tempdir,
                format="raw",
                source_atom_ids=source_ids,
                return_mapping=True,
            )
            returned_bundle, bundle_mapping = save_sponge_input(
                molecule,
                "fixture",
                tempdir,
                format="bundle",
                source_atom_ids=source_ids,
                return_mapping=True,
            )
        self.assertIs(returned_raw, molecule)
        self.assertIs(returned_bundle, molecule)
        self.assertEqual(raw_mapping, bundle_mapping)
        self.assertEqual(
            tuple(record["simulation_index"] for record in raw_mapping),
            tuple(range(len(molecule.atoms))),
        )
        self.assertEqual(
            {record["source_atom_id"] for record in raw_mapping},
            set(atoms_by_external_id),
        )

    def test_nonbonded_provider_outputs_compose_into_apply_ready_result(self):
        request = _nonbonded_request()
        source = _overlay_result(request)
        base_report = BaseCompositionReport(
            schema_version=1,
            topology_hash=request.topology.topology_hash,
            projection_hash=request.projection_hash,
            providers=("fixture-base",),
            provider_report_hashes={"fixture-base": "a" * 64},
            covered_atom_ids=source.base_overlay.covered_atom_ids,
            bonded_term_count=len(source.base_overlay.bonded_parameters),
            parameter_source=source.base_overlay.parameter_source,
        ).with_computed_hash()
        metal_report = MetalOverlayBuildReport(
            schema_version=1,
            topology_hash=request.topology.topology_hash,
            projection_hash=request.projection_hash,
            interaction_model="nonbonded_12_6",
            covered_atom_ids=source.metal_overlay.covered_atom_ids,
            parameter_source=source.metal_overlay.parameter_source,
            precedence=source.metal_overlay.precedence,
            atom_type_count=1,
            charge_count=1,
            mass_count=1,
            lj_count=1,
            bonded_term_count=0,
            provenance={"provider_kind": "fixture-ion"},
        ).with_computed_hash()
        result = compose_nonbonded_assignment(
            request,
            BaseCompositionOutput(source.base_overlay, base_report),
            MetalOverlayBuildOutput(source.metal_overlay, metal_report),
            package_hash="b" * 64,
        )
        self.assertEqual(result.status, "overlay_validated")
        self.assertFalse(result.complete)
        self.assertIsNone(result.charge_overlay)
        self.assertIsNone(result.bonded_overlay)
        self.assertEqual(result.fit_reports["metal_assignment"]["report_hash"], metal_report.report_hash)

    def test_apply_materializes_xponge_molecule_without_parent_registry_mutation(self):
        request = _complete_request()
        result = _complete_result(request)
        atom_types_before = tuple(AtomType.get_all_types())
        molecules_before = tuple(Molecule._all)
        output = apply_parameterization(request, result, _force_protocol())
        self.assertTrue(output.result.complete)
        self.assertEqual(output.result.status, "assigned_system_validated")
        self.assertTrue(output.assigned_system.complete)
        self.assertEqual(output.assigned_system.application_audit["atom_count"], 4)
        self.assertEqual(output.assigned_system.application_audit["residue_count"], 2)
        self.assertEqual(output.assigned_system.application_audit["xponge_residue_link_count"], 1)
        self.assertEqual(output.assigned_system.application_audit["effective_type_count"], 4)
        self.assertEqual(output.assigned_system.application_audit["materialized_lj_type_count"], 3)
        self.assertEqual(output.assigned_system.application_audit["topology_hash_after"], request.topology.topology_hash)
        self.assertEqual(tuple(AtomType.get_all_types()), atom_types_before)
        self.assertEqual(tuple(Molecule._all), molecules_before)
        heme = next(residue for residue in output.assigned_system.topology.residues if residue.external_id == "HEM:A:1")
        self.assertIn("heme-fe", heme.atom_ids)

    def test_apply_accepts_simulation_residue_order_independent_of_stable_atom_order(self):
        request = _complete_request()
        topology = replace(
            request.topology,
            residues=tuple(reversed(request.topology.residues)),
            topology_hash="",
        ).with_computed_hash()
        request = replace(request, topology=topology, projection_hash="").with_computed_hash()

        output = apply_parameterization(request, _complete_result(request), _force_protocol())
        molecule, atoms_by_external_id, _ = materialize_assigned_system(
            output.assigned_system
        )
        source_ids = {
            atom: external_id
            for external_id, atom in atoms_by_external_id.items()
        }
        with TemporaryDirectory(prefix="xponge-native-order-") as tempdir:
            _, mapping = save_sponge_input(
                molecule,
                "fixture",
                tempdir,
                format="raw",
                source_atom_ids=source_ids,
                return_mapping=True,
            )

        self.assertTrue(output.assigned_system.complete)
        self.assertEqual(
            {item["source_atom_id"] for item in mapping},
            {atom.external_id for atom in topology.atoms},
        )
        self.assertEqual(len(mapping), len({item["source_atom_id"] for item in mapping}))
        self.assertEqual(len(output.assigned_system.application_audit["materialized_atom_order_hash"]), 64)
        self.assertEqual(output.assigned_system.application_audit["materialized_lj_type_count"], 3)

    def test_assigned_system_and_apply_input_are_hash_closed(self):
        request = _complete_request()
        result = _complete_result(request)
        output = apply_parameterization(request, result, _force_protocol())
        payload = assigned_system_dumps(output.assigned_system)
        self.assertEqual(assigned_system_loads(payload), output.assigned_system)
        apply_input = ParameterApplyInput(1, result, _force_protocol()).with_computed_hash()
        input_payload = parameter_apply_input_dumps(request, apply_input)
        self.assertEqual(parameter_apply_input_loads(request, input_payload), apply_input)

    def test_apply_rejects_missing_atom_charge_before_worker(self):
        request = _complete_request()
        result = _complete_result(request)
        charges = dict(result.base_overlay.charges)
        charges.pop("heme-c")
        broken = replace(
            result,
            base_overlay=replace(result.base_overlay, charges=charges),
            result_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "incomplete_atom_parameters"):
            apply_parameterization(request, broken, _force_protocol())

    def test_apply_rejects_missing_prepared_link_parameter(self):
        request = _complete_request()
        result = _complete_result(request)
        terms = dict(result.bonded_overlay.terms)
        terms.pop("metal:fe-his")
        broken = replace(
            result,
            bonded_overlay=replace(result.bonded_overlay, terms=terms),
            result_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "missing_bonded_edge_parameter"):
            apply_parameterization(request, broken, _force_protocol())

if __name__ == "__main__":
    unittest.main()
