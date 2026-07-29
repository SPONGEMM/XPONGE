"""Transactional local metal-parameter patch tests."""

from dataclasses import replace
import unittest
from uuid import uuid4

import Xponge.forcefield.amber  # noqa: F401 - register native force types
from Xponge.forcefield.base.lj_base import LJType
from Xponge.helper import AtomType, Molecule, ResidueType
from Xponge.metal_assignment import (
    BaseForceFieldOverlay,
    ValidationError,
    apply,
    build_metal_parameter_patch,
    metal_parameter_patch_dumps,
    metal_parameter_patch_loads,
    prepare_residue_templates,
)
from Xponge.tools.unittests.test_1_metal_assignment_contracts import (
    _bonded_request,
    _overlay_result,
)


def _local_patch():
    request = _bonded_request()
    result = _overlay_result(request)
    result = replace(
        result,
        base_overlay=BaseForceFieldOverlay(
            topology_hash=request.topology.topology_hash,
            parameter_source="xponge:preassigned-molecule",
        ),
        provenance={
            **result.provenance,
            "base_assignment": "preassigned_molecule",
        },
        result_hash="",
    ).with_computed_hash()
    return build_metal_parameter_patch(
        request,
        result,
        site_ids=("site:heme-fe",),
    )


def _ordinary_molecule(patch):
    namespace = f"{patch.patch_hash[:8]}_{uuid4().hex[:8]}"
    lj_name = f"PATCH_TEST_LJ_{namespace}"
    LJType(
        name=f"{lj_name}-{lj_name}",
        epsilon=0.2,
        rmin=1.8,
    )
    residue_type = ResidueType(name=f"PATCH_TEST_RES_{namespace}")
    for index, identity in enumerate(patch.atoms):
        atom_type = AtomType(
            name=f"PATCH_TEST_ATOM_{namespace}_{index}",
            charge=-0.25 if not identity.external_id.endswith("fe") else 0.5,
            mass=12.0,
            LJtype=lj_name,
        )
        residue_type.add_atom(
            identity.atom_name,
            atom_type,
            float(index),
            0.0,
            0.0,
        )
    molecule = Molecule(name=f"patch_test_{namespace}")
    molecule.add_residue(residue_type)
    molecule.get_atoms()
    mapping = {}
    for identity, atom in zip(patch.atoms, molecule.atoms):
        atom.element = identity.element
        mapping[identity.external_id] = molecule.atom_index[atom]
    molecule.built = True
    return molecule, mapping


class MetalAssignmentApplyTests(unittest.TestCase):
    def test_prepare_residue_templates_adds_only_missing_embedded_metal(self):
        patch = _local_patch()
        residue_name = f"PATCH_EMBEDDED_{uuid4().hex[:8]}"
        residue_type = ResidueType(name=residue_name)
        metal_id = patch.target_metal_atom_ids[0]
        metal_identity = next(
            atom for atom in patch.atoms
            if atom.external_id == metal_id
        )

        report = prepare_residue_templates(patch, [{
            "external_id": metal_id,
            "residue_name": residue_name,
            "atom_name": metal_identity.atom_name,
        }])

        self.assertEqual(report["prepared_atom_ids"], (metal_id,))
        self.assertEqual(
            residue_type.name2atom(metal_identity.atom_name).charge,
            patch.parameterization_result.metal_overlay.charges[metal_id],
        )
        repeated = prepare_residue_templates(patch, [{
            "external_id": metal_id,
            "residue_name": residue_name,
            "atom_name": metal_identity.atom_name,
        }])
        self.assertEqual(repeated["prepared_atom_ids"], ())

    def test_hash_closed_patch_applies_to_an_existing_molecule(self):
        patch = _local_patch()
        molecule, mapping = _ordinary_molecule(patch)
        nonmetal_id = next(
            atom.external_id
            for atom in patch.atoms
            if atom.external_id not in patch.target_metal_atom_ids
        )
        molecule.atoms[mapping[nonmetal_id]].charge = -0.123
        original_nonmetal_charge = molecule.atoms[mapping[nonmetal_id]].charge
        original_coordinates = [
            (atom.x, atom.y, atom.z)
            for atom in molecule.atoms
        ]

        output = apply(molecule, patch, mapping)

        self.assertIsNot(output.molecule, molecule)
        self.assertEqual(output.application_report["patch_hash"], patch.patch_hash)
        self.assertTrue(output.application_report["topology_preserved"])
        self.assertTrue(output.application_report["ordinary_force_field_preserved"])
        self.assertEqual(output.application_report["force_counts"], {"bond": 1})
        self.assertAlmostEqual(
            output.molecule.atoms[mapping[nonmetal_id]].charge,
            original_nonmetal_charge,
        )
        self.assertEqual(
            {
                identity.external_id: output.molecule.atoms[
                    mapping[identity.external_id]
                ].name
                for identity in patch.atoms
            },
            {
                identity.external_id: identity.atom_name
                for identity in patch.atoms
            },
        )
        self.assertEqual(
            [(atom.x, atom.y, atom.z) for atom in output.molecule.atoms],
            original_coordinates,
        )
        metal_id = patch.target_metal_atom_ids[0]
        self.assertEqual(
            output.molecule.atoms[mapping[metal_id]].mass,
            patch.parameterization_result.metal_overlay.masses[metal_id],
        )

    def test_identity_mismatch_is_rejected_without_mutating_the_input(self):
        patch = _local_patch()
        molecule, mapping = _ordinary_molecule(patch)
        target_id = patch.atoms[0].external_id
        target_atom = molecule.atoms[mapping[target_id]]
        original_charge = target_atom.charge
        target_atom.name = "WRONG"

        with self.assertRaisesRegex(
            ValidationError,
            "molecule_atom_identity_mismatch",
        ):
            apply(molecule, patch, mapping)

        self.assertEqual(target_atom.charge, original_charge)
        self.assertEqual(target_atom.name, "WRONG")

    def test_patch_wire_round_trip_is_hash_closed(self):
        patch = _local_patch()
        payload = metal_parameter_patch_dumps(patch)
        self.assertEqual(metal_parameter_patch_loads(payload), patch)
        self.assertIn(patch.patch_hash, payload)

    def test_legacy_request_result_apply_signature_is_removed(self):
        patch = _local_patch()
        molecule, mapping = _ordinary_molecule(patch)
        request = _bonded_request()
        result = _overlay_result(request)
        with self.assertRaises(TypeError):
            apply(molecule, request, result, mapping)


if __name__ == "__main__":
    unittest.main()
