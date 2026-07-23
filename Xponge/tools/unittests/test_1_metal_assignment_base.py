"""Native base-force-field assignment tests."""

from __future__ import annotations

from copy import deepcopy
import math
import unittest

from Xponge.assign import AssignRule
from Xponge.metal_assignment import ValidationError
from Xponge.metal_assignment.base_assignment import (
    WORKER_PROTOCOL_VERSION,
    _invoke_native_worker,
    _validate_worker_response,
)
from Xponge.metal_assignment.metal_overlay import _invoke_ion_worker
from Xponge.metal_assignment.standard_assignment import _invoke_worker as _invoke_standard_worker
from Xponge.metal_assignment.standard_assignment import _validate_response as _validate_standard_response


def _benzene_payload(provider: str = "gaff2") -> dict:
    carbon_coordinates = [
        [math.cos(index * math.pi / 3), math.sin(index * math.pi / 3), 0.0]
        for index in range(6)
    ]
    atoms = [
        {
            "external_id": f"atom:C{index + 1}",
            "element": "C",
            "coordinates": coordinates,
            "name": f"C{index + 1}",
            "formal_charge": None,
        }
        for index, coordinates in enumerate(carbon_coordinates)
    ]
    atoms.extend(
        {
            "external_id": f"atom:H{index + 1}",
            "element": "H",
            "coordinates": [1.8 * coordinates[0], 1.8 * coordinates[1], 0.0],
            "name": f"H{index + 1}",
            "formal_charge": None,
        }
        for index, coordinates in enumerate(carbon_coordinates)
    )
    bonds = [
        {
            "external_id": f"bond:CC{index + 1}",
            "atom_ids": [f"atom:C{index + 1}", f"atom:C{(index + 1) % 6 + 1}"],
            "order": 1.5,
            "semantic": "covalent",
            "source": "unit-test-locked-graph",
        }
        for index in range(6)
    ]
    bonds.extend(
        {
            "external_id": f"bond:CH{index + 1}",
            "atom_ids": [f"atom:C{index + 1}", f"atom:H{index + 1}"],
            "order": 1.0,
            "semantic": "covalent",
            "source": "unit-test-locked-graph",
        }
        for index in range(6)
    )
    return {
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": provider,
        "topology_hash": "1" * 64,
        "complete_missing_parameters": True,
        "base_charge_input": None,
        "components": [{
            "external_id": "component:benzene",
            "net_formal_charge": 0,
            "proof_hash": "2" * 64,
            "artifact_hash": "3" * 64,
            "atoms": atoms,
            "bonds": bonds,
        }],
    }


def _nala_payload() -> dict:
    names_and_elements = (
        ("N", "N"), ("H1", "H"), ("H2", "H"), ("H3", "H"),
        ("CA", "C"), ("HA", "H"), ("CB", "C"), ("HB1", "H"),
        ("HB2", "H"), ("HB3", "H"), ("C", "C"), ("O", "O"),
    )
    atoms = [
        {"external_id": f"atom:{name}", "atom_name": name, "element": element}
        for name, element in names_and_elements
    ]
    pairs = (
        ("C", "O"), ("CB", "HB1"), ("CB", "HB2"), ("CB", "HB3"),
        ("CA", "HA"), ("CA", "CB"), ("CA", "C"), ("N", "H1"),
        ("N", "H2"), ("N", "H3"), ("N", "CA"),
    )
    return {
        "protocol_version": 1,
        "provider": "standard_biomolecular",
        "force_field": "ff14sb",
        "water_model": "tip3p",
        "topology_hash": "4" * 64,
        "links": [],
        "components": [{
            "external_id": "assignment:ALA:1",
            "proof_hash": "5" * 64,
            "artifact_hash": "6" * 64,
            "residue_name": "ALA",
            "polymer_position": "start",
            "template_id": "ALA",
            "atoms": atoms,
            "bonds": [
                {
                    "external_id": f"bond:{index + 1}",
                    "atom_ids": [f"atom:{atom1}", f"atom:{atom2}"],
                }
                for index, (atom1, atom2) in enumerate(pairs)
            ],
        }],
    }


class NativeBaseAssignmentTests(unittest.TestCase):
    def test_isolated_gaff2_worker_preserves_parent_registries_and_locked_aromatic_order(self):
        payload = _benzene_payload()
        before = tuple(AssignRule.all.keys())
        response = _invoke_native_worker(payload, timeout_seconds=30.0)
        _validate_worker_response(payload, response)
        self.assertEqual(tuple(AssignRule.all.keys()), before)

        component = response["components"][0]
        self.assertEqual(
            [component["atom_types"][f"atom:C{index + 1}"] for index in range(6)],
            ["ca"] * 6,
        )
        self.assertEqual(
            [component["atom_types"][f"atom:H{index + 1}"] for index in range(6)],
            ["ha"] * 6,
        )
        self.assertEqual(set(component["lj_parameters"]), {atom["external_id"] for atom in payload["components"][0]["atoms"]})
        self.assertTrue(component["frcmod_sha256"])
        self.assertTrue(component["bonded_parameters"])
        self.assertEqual(
            {term["kind"] for term in component["bonded_parameters"].values()},
            {"bond", "angle", "proper_dihedral", "improper_dihedral"},
        )

    def test_worker_is_deterministic(self):
        payload = _benzene_payload("gaff")
        first = _invoke_native_worker(payload, timeout_seconds=30.0)
        second = _invoke_native_worker(payload, timeout_seconds=30.0)
        self.assertEqual(first, second)

    def test_worker_rejects_invalid_locked_bond_order(self):
        payload = deepcopy(_benzene_payload())
        payload["components"][0]["bonds"][0]["order"] = -1
        with self.assertRaisesRegex(ValidationError, "base_assignment_worker_failed"):
            _invoke_native_worker(payload, timeout_seconds=30.0)

    def test_isolated_standard_worker_matches_exact_terminal_template(self):
        payload = _nala_payload()
        before = tuple(AssignRule.all.keys())
        response = _invoke_standard_worker(payload, timeout_seconds=30.0)
        _validate_standard_response(payload, response)
        self.assertEqual(tuple(AssignRule.all.keys()), before)
        component = response["components"][0]
        self.assertEqual(component["template_name"], "NALA")
        expected_atom_ids = {atom["external_id"] for atom in payload["components"][0]["atoms"]}
        self.assertEqual(set(component["atom_types"]), expected_atom_ids)
        self.assertEqual(set(component["charges"]), expected_atom_ids)
        self.assertEqual(set(component["masses"]), expected_atom_ids)
        self.assertEqual(set(component["lj_parameters"]), expected_atom_ids)
        self.assertTrue(component["bonded_parameters"])

    def test_standard_worker_rejects_missing_template_atom(self):
        payload = _nala_payload()
        payload["components"][0]["atoms"] = [
            atom for atom in payload["components"][0]["atoms"] if atom["atom_name"] != "H3"
        ]
        payload["components"][0]["bonds"] = [
            bond for bond in payload["components"][0]["bonds"] if "atom:H3" not in bond["atom_ids"]
        ]
        with self.assertRaisesRegex(ValidationError, "standard_assignment_worker_failed"):
            _invoke_standard_worker(payload, timeout_seconds=30.0)

    def test_isolated_tip3p_ion_worker_resolves_zn_and_fe_without_parent_registry_changes(self):
        payload = {
            "protocol_version": 1,
            "topology_hash": "4" * 64,
            "water_model": "tip3p",
            "metals": [
                {"external_id": "metal:zn", "element": "Zn", "formal_charge": 2},
                {"external_id": "metal:fe2", "element": "Fe", "formal_charge": 2},
                {"external_id": "metal:fe3", "element": "Fe", "formal_charge": 3},
            ],
        }
        before = tuple(AssignRule.all.keys())
        response = _invoke_ion_worker(payload, timeout_seconds=30.0)
        self.assertEqual(tuple(AssignRule.all.keys()), before)
        parameters = response["parameters"]
        self.assertEqual(parameters["metal:zn"]["atom_type"], "Zn2+")
        self.assertEqual(parameters["metal:fe2"]["atom_type"], "Fe2+")
        self.assertEqual(parameters["metal:fe3"]["atom_type"], "Fe3+")
        self.assertAlmostEqual(parameters["metal:zn"]["mass"], 65.4)
        self.assertAlmostEqual(parameters["metal:zn"]["lj"]["rmin"], 1.271)
        self.assertAlmostEqual(parameters["metal:zn"]["lj"]["epsilon"], 0.00330286)

    def test_ion_worker_rejects_unavailable_charge_state(self):
        payload = {
            "protocol_version": 1,
            "topology_hash": "4" * 64,
            "water_model": "tip3p",
            "metals": [{"external_id": "metal:zn", "element": "Zn", "formal_charge": 7}],
        }
        with self.assertRaisesRegex(ValidationError, "ion_assignment_worker_failed"):
            _invoke_ion_worker(payload, timeout_seconds=30.0)


if __name__ == "__main__":
    unittest.main()
