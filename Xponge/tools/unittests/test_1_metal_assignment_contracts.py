"""Contract tests for topology-preserving metal assignment."""

from dataclasses import replace
import hashlib
import json
import unittest

from Xponge.metal_assignment import (
    AssignmentComponent,
    BondedFitInput,
    BondedMetalParameterSpec,
    BaseForceFieldOverlay,
    BondedParameterOverlay,
    ChargeAssignmentContract,
    ChargePolicy,
    ChemicalTopologyProof,
    ElectronicState,
    EMPIRICAL_REGISTRY,
    EmpiricalApplicability,
    EmpiricalCitation,
    EmpiricalRegistryEntry,
    MetalAssignmentInput,
    MetalParameterOverlay,
    MetalAtomParameterSpec,
    ParameterizationResult,
    PartialChargeArtifact,
    PartialChargeArtifactBundle,
    PreparedAtom,
    PreparedBond,
    PreparedChemicalTopology,
    PreparedComponent,
    PreparedLink,
    PreparedResidue,
    ProviderCapabilitySnapshot,
    ProviderProjection,
    ResiduePartitionProof,
    SimulationSplitReason,
    ValidationError,
    metal_assignment_input_dumps,
    metal_assignment_input_from_dict,
    metal_assignment_input_loads,
    metal_assignment_input_to_dict,
    parameterization_result_dumps,
    parameterization_result_loads,
    bonded_fit_input_dumps,
    bonded_fit_input_loads,
    compose_partial_charge_artifacts,
    empirical_registry_descriptor,
    manual_bonded_terms,
    resolve_empirical_registry,
    validate_empirical_registry_entry,
    validate_input,
    validate_result,
)


INPUT_HASH = hashlib.sha256(b"bonded-heme-fixture").hexdigest()


CAPABILITIES = ProviderCapabilitySnapshot(
    schema_version=1,
    provider_name="xponge",
    provider_version="1.7b3-dev",
    provider_revision="fixture-revision",
    projection_schema_version=1,
    embedded_metal=True,
    standalone_metal=True,
    multi_metal=False,
    bonded=True,
    nonbonded_12_6=True,
    base_force_field_providers=("gaff", "gaff2"),
)


def _atom(
    external_id,
    stable_order,
    element,
    coordinates,
    chemical_component_id,
    canonical_residue_id,
    simulation_component_id,
    simulation_residue_id,
    *,
    formal_charge=0,
    is_metal=False,
    scopes=(),
):
    canonical_atom_id = stable_order + 1
    return PreparedAtom(
        external_id,
        canonical_atom_id,
        stable_order,
        element,
        coordinates,
        chemical_component_id,
        canonical_residue_id,
        simulation_component_id,
        simulation_residue_id,
        formal_charge,
        formal_charge_known=formal_charge is not None,
        is_metal=is_metal,
        scopes=scopes,
        atom_name=external_id.upper(),
    )


def _bonded_heme_topology():
    atoms = (
        _atom("heme-fe", 0, "Fe", (0.0, 0.0, 0.0), "HEM-component", "HEM:A:1", "HEM-component", "HEM:A:1", formal_charge=2, is_metal=True, scopes=("qm", "metalOverlay", "simulation")),
        _atom("heme-na", 1, "N", (1.9, 0.0, 0.0), "HEM-component", "HEM:A:1", "HEM-component", "HEM:A:1", scopes=("baseForceField", "qm", "simulation")),
        _atom("heme-c", 2, "C", (3.2, 0.0, 0.0), "HEM-component", "HEM:A:1", "HEM-component", "HEM:A:1", scopes=("baseForceField", "qm", "simulation")),
        _atom("his-ne2", 3, "N", (0.0, 2.1, 0.0), "protein-component", "HIS:A:2", "protein-component", "HIS:A:2", scopes=("qm", "simulation")),
    )
    topology = PreparedChemicalTopology(
        schema_version=1,
        graph_revision=7,
        input_hash=INPUT_HASH,
        coordinate_unit="angstrom",
        atoms=atoms,
        residues=(
            PreparedResidue("HEM:A:1", "HEM:A:1", ("heme-fe", "heme-na", "heme-c"), "HEM-component", "HEM-component", 2, "fixture", True, "HEM", "A", 1),
            PreparedResidue("HIS:A:2", "HIS:A:2", ("his-ne2",), "protein-component", "protein-component", 0, "fixture", True, "HIS", "A", 2, "", "single", "HIS"),
        ),
        components=(
            PreparedComponent("HEM-component", ("heme-fe", "heme-na", "heme-c"), "HEM-component", "HEM-component", 2, "fixture", True),
            PreparedComponent("protein-component", ("his-ne2",), "protein-component", "protein-component", 0, "fixture", True),
        ),
        bonds=(
            PreparedBond("fe-na", ("heme-fe", "heme-na"), semantic="coordination", source="confirmed_input"),
            PreparedBond("na-c", ("heme-na", "heme-c"), source="input"),
        ),
        links=(PreparedLink("fe-his", ("heme-fe", "his-ne2"), source="struct_conn"),),
    )
    return topology.with_computed_hash()


def _bonded_request():
    topology = _bonded_heme_topology()
    request = MetalAssignmentInput(
        schema_version=1,
        request_id="request-heme-1",
        interaction_model="bonded",
        electronic_state=ElectronicState(net_charge=2, spin_multiplicity=1),
        graph_revision=topology.graph_revision,
        input_hash=topology.input_hash,
        topology=topology,
        projections=(
            ProviderProjection(
                "heme-metal", "HEM-component", "HEM:A:1", "heme-metal", "HEM-component", "HEM:A:1",
                ("heme-fe",), {"heme-fe": ("qm", "metalOverlay", "simulation")}, assignment_only_split=True,
                parent_component_id="HEM-component",
            ),
            ProviderProjection(
                "heme-ligand", "HEM-component", "HEM:A:1", "heme-ligand", "HEM-component", "HEM:A:1",
                ("heme-na", "heme-c"), assignment_only_split=True, parent_component_id="HEM-component",
            ),
        ),
        assignment_components=(
            AssignmentComponent(
                "heme-metal", ("heme-fe",), "atomic_center", "metal", 2,
                "atom_formal_charge_sum", True,
                ChemicalTopologyProof(
                    1, 7, "fixture-perception-v1", ("heme-fe",), {"heme-fe": 0},
                    {"heme-fe": ()}, {"heme-fe": 0.0}, {"heme-fe": "not_applicable_metal"},
                    (), "not_applicable", "atomic_metal_not_applicable", True,
                ).with_computed_hash(),
                base_force_field=False,
            ),
            AssignmentComponent(
                "heme-ligand", ("heme-na", "heme-c"), "conjugated_macrocycle", "gaff2", 0,
                "atom_formal_charge_sum", True,
                ChemicalTopologyProof(
                    1, 7, "fixture-perception-v1", ("heme-na", "heme-c"),
                    {"heme-na": 0, "heme-c": 0}, {"heme-na": (0,), "heme-c": (0,)},
                    {"heme-na": 3.0, "heme-c": 4.0},
                    {"heme-na": "valid_unique", "heme-c": "valid_unique"},
                    (), "complete", "component_charge_constrained_valence", True,
                ).with_computed_hash(),
            ),
        ),
        partition_proofs=(
            ResiduePartitionProof(
                "HEM:A:1",
                ("heme-fe", "heme-na", "heme-c"),
                ("fe-na", "na-c"),
                (),
                {"HEM:A:1": ("heme-fe", "heme-na", "heme-c")},
            ).with_computed_hash(),
            ResiduePartitionProof(
                "HIS:A:2",
                ("his-ne2",),
                (),
                (),
                {"HIS:A:2": ("his-ne2",)},
            ).with_computed_hash(),
        ),
        capability_snapshot=CAPABILITIES,
        charge_contract=ChargeAssignmentContract(
            ChargePolicy.SITE_JOINT_FIT,
            ("heme-fe", "heme-na", "heme-c"),
            {"HEM-component": 2},
            site_target=2.0,
        ),
    )
    return request.with_computed_hash()


def _nonbonded_request():
    request = _bonded_request()
    topology = request.topology
    atoms = tuple(
        replace(atom, simulation_component_id="metal-component", simulation_residue_id="FE:A:1")
        if atom.external_id == "heme-fe"
        else replace(atom, simulation_component_id="HEM-ligand-component", simulation_residue_id="HEM-LIG:A:1")
        if atom.external_id.startswith("heme-")
        else atom
        for atom in topology.atoms
    )
    prepared = replace(
        topology,
        atoms=atoms,
        residues=(
            PreparedResidue("FE:A:1", "HEM:A:1", ("heme-fe",), "HEM-component", "metal-component", 2, "fixture", True, "HEM", "A", 1),
            PreparedResidue("HEM-LIG:A:1", "HEM:A:1", ("heme-na", "heme-c"), "HEM-component", "HEM-ligand-component", 0, "fixture", True, "HEM", "A", 1),
            topology.residues[1],
        ),
        components=(
            PreparedComponent("metal-component", ("heme-fe",), "HEM-component", "metal-component", 2, "fixture", True),
            PreparedComponent("HEM-ligand-component", ("heme-na", "heme-c"), "HEM-component", "HEM-ligand-component", 0, "fixture", True),
            topology.components[1],
        ),
        bonds=(PreparedBond("na-c", ("heme-na", "heme-c"), source="input"),),
        links=(),
        topology_hash="",
    ).with_computed_hash()
    projections = (
        replace(
            request.projections[0],
            simulation_component_id="metal-component",
            simulation_residue_id="FE:A:1",
            simulation_split=True,
            simulation_split_reason=SimulationSplitReason.DISCONNECTED_AFTER_INTERACTION_FILTER,
        ),
        replace(
            request.projections[1],
            simulation_component_id="HEM-ligand-component",
            simulation_residue_id="HEM-LIG:A:1",
            simulation_split=True,
            simulation_split_reason=SimulationSplitReason.DISCONNECTED_AFTER_INTERACTION_FILTER,
        ),
    )
    partition_proofs = (
        ResiduePartitionProof(
            "HEM:A:1",
            ("heme-fe", "heme-na", "heme-c"),
            ("na-c",),
            ("fe-na",),
            {"FE:A:1": ("heme-fe",), "HEM-LIG:A:1": ("heme-na", "heme-c")},
        ).with_computed_hash(),
        request.partition_proofs[1],
    )
    return replace(
        request,
        interaction_model="nonbonded_12_6",
        topology=prepared,
        projections=projections,
        partition_proofs=partition_proofs,
        charge_contract=replace(
            request.charge_contract,
            component_formal_charges={"metal-component": 2, "HEM-ligand-component": 0},
        ),
        projection_hash="",
    ).with_computed_hash()


def _overlay_result(request):
    result = ParameterizationResult(
        schema_version=1,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology_hash=request.topology.topology_hash,
        projection_hash=request.projection_hash,
        base_overlay=BaseForceFieldOverlay(
            topology_hash=request.topology.topology_hash,
            covered_atom_ids=("heme-na", "heme-c"),
            atom_types={"heme-na": "n", "heme-c": "c"},
            lj_parameters={
                atom_id: {
                    "epsilon": epsilon,
                    "rmin": rmin,
                    "energy_unit": "kcal/mol",
                    "length_unit": "angstrom",
                    "source": "gaff2.dat",
                }
                for atom_id, epsilon, rmin in (
                    ("heme-na", 0.17, 1.824),
                    ("heme-c", 0.086, 1.908),
                )
            },
            bonded_parameters={
                "base-bond-na-c": {
                    "kind": "bond",
                    "atom_ids": ("heme-na", "heme-c"),
                    "parameters": {
                        "k": 300.0,
                        "equilibrium": 1.3,
                        "k_unit": "kcal/mol/angstrom^2",
                        "equilibrium_unit": "angstrom",
                    },
                    "source": "gaff2.dat",
                },
            },
            parameter_source="gaff2:metal-assignment-base-ff-v1",
        ),
        metal_overlay=MetalParameterOverlay(
            topology_hash=request.topology.topology_hash,
            covered_atom_ids=("heme-fe",),
            atom_types={"heme-fe": "FE-site-1"},
            masses={"heme-fe": 55.845},
            lj_parameters={
                "heme-fe": {
                    "epsilon": 0.01,
                    "rmin": 1.2,
                    "energy_unit": "kcal/mol",
                    "length_unit": "angstrom",
                    "source": "explicit-test",
                },
            },
            parameter_source="explicit-test",
        ),
        bonded_overlay=BondedParameterOverlay(
            request.topology.topology_hash,
            {
                "metal-bond-fe-na": {
                    "kind": "bond",
                    "atom_ids": ("heme-fe", "heme-na"),
                    "parameters": {
                        "k": 100.0,
                        "equilibrium": 1.9,
                        "k_unit": "kcal/mol/angstrom^2",
                        "equilibrium_unit": "angstrom",
                    },
                    "source": "seminario-test",
                },
            },
            "seminario-test",
        ),
        fit_reports={"base_assignment": {"provider": "gaff2"}},
        provenance={"protocol_id": "metal-assignment-base-ff-v1", "provider_revision": "fixture"},
        status="overlay_validated",
    )
    return result.with_computed_hash()


class MetalAssignmentContractTests(unittest.TestCase):
    def test_bonded_embedded_metal_is_assignment_only_split(self):
        request = _bonded_request()
        validate_input(request)
        self.assertTrue(request.projections[0].assignment_only_split)
        self.assertFalse(request.projections[0].simulation_split)
        self.assertEqual(request.topology.atoms[0].simulation_residue_id, "HEM:A:1")

    def test_nonbonded_embedded_metal_uses_separate_simulation_residue(self):
        request = _nonbonded_request()
        validate_input(request)
        self.assertTrue(request.projections[0].simulation_split)
        self.assertNotEqual(request.projections[0].simulation_component_id, request.projections[1].simulation_component_id)

    def test_wire_round_trip_is_canonical_and_strict(self):
        request = _bonded_request()
        payload = metal_assignment_input_dumps(request)
        restored = metal_assignment_input_loads(payload)
        self.assertEqual(restored, request)
        self.assertEqual(metal_assignment_input_dumps(restored), payload)
        self.assertIn('"stable_order":0', payload)
        self.assertIn('"projection_hash"', payload)

    def test_unknown_schema_version_and_field_are_rejected(self):
        payload = metal_assignment_input_to_dict(_bonded_request())
        with self.assertRaisesRegex(ValidationError, "unsupported_schema_version"):
            metal_assignment_input_from_dict({**payload, "schema_version": 2})
        with self.assertRaisesRegex(ValidationError, "unknown_wire_field"):
            metal_assignment_input_from_dict({**payload, "legacy_unknown": True})

    def test_partial_charge_artifacts_are_hash_closed_and_precedence_is_explicit(self):
        request = _nonbonded_request()
        component = next(
            item for item in request.topology.components if item.external_id == "HEM-ligand-component"
        )

        def artifact(artifact_id: str, precedence: int, charges=(-0.4, 0.4)):
            return PartialChargeArtifact(
                schema_version=1,
                artifact_id=artifact_id,
                topology_hash=request.topology.topology_hash,
                graph_revision=request.graph_revision,
                input_hash=request.input_hash,
                atom_ids=component.atom_ids,
                charges=charges,
                scope_kind="component",
                scope_id=component.external_id,
                atomic_charge_role="reference",
                precedence=precedence,
                provider="fixture-charge-provider",
                provider_version="1",
                method="fixture-explicit-charge",
                charge_unit="elementary_charge",
                source="unit-test",
            ).with_computed_hash()

        low = artifact("ligand-reference", 100)
        high = artifact("ligand-override", 150, (-0.5, 0.5))
        bundle = PartialChargeArtifactBundle(
            1,
            request.topology.topology_hash,
            request.graph_revision,
            request.input_hash,
            (low, high),
        ).with_computed_hash()
        request = replace(
            request,
            partial_charge_artifacts=bundle,
            projection_hash="",
        ).with_computed_hash()
        validate_input(request)
        restored = metal_assignment_input_loads(metal_assignment_input_dumps(request))
        self.assertEqual(restored, request)
        overlay, report = compose_partial_charge_artifacts(request)
        self.assertEqual(dict(overlay.charges), {"heme-na": -0.5, "heme-c": 0.5})
        self.assertEqual(overlay.overlay_hash, overlay.computed_hash())
        self.assertEqual(
            {source["artifact_id"] for source in overlay.atom_sources.values()},
            {"ligand-override"},
        )
        self.assertEqual(len(report["overridden"]), 2)

        ambiguous = replace(high, artifact_id="ligand-tie", precedence=100, artifact_hash="").with_computed_hash()
        bad_bundle = replace(bundle, artifacts=(low, ambiguous), bundle_hash="").with_computed_hash()
        bad_request = replace(
            request,
            partial_charge_artifacts=bad_bundle,
            projection_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "ambiguous_partial_charge_precedence"):
            validate_input(bad_request)

    def test_provider_capability_is_a_hard_gate(self):
        request = _bonded_request()
        unsupported = replace(request.capability_snapshot, bonded=False)
        bad = replace(request, capability_snapshot=unsupported, projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "unsupported_provider_capability"):
            validate_input(bad)

    def test_partition_proof_is_complete_and_hash_guarded(self):
        request = _bonded_request()
        incomplete = replace(request, partition_proofs=request.partition_proofs[:1], projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "incomplete_partition_proof"):
            validate_input(incomplete)
        stale_proof = replace(request.partition_proofs[0], proof_hash="stale")
        stale = replace(request, partition_proofs=(stale_proof, request.partition_proofs[1]), projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "stale_partition_proof_hash"):
            validate_input(stale)

    def test_nonfinite_coordinates_are_structured_failures(self):
        request = _bonded_request()
        atom = replace(request.topology.atoms[0], coordinates=(float("nan"), 0.0, 0.0))
        topology = replace(request.topology, atoms=(atom, *request.topology.atoms[1:]), topology_hash="x")
        with self.assertRaisesRegex(ValidationError, "invalid_coordinates"):
            validate_input(replace(request, topology=topology))

    def test_unknown_nonmetal_atom_formal_charge_is_preserved(self):
        request = _bonded_request()
        atom = replace(request.topology.atoms[1], formal_charge=None, formal_charge_known=False)
        residue = replace(request.topology.residues[0], atom_formal_charges_complete=False)
        component = replace(request.topology.components[0], atom_formal_charges_complete=False)
        topology = replace(
            request.topology,
            atoms=(request.topology.atoms[0], atom, *request.topology.atoms[2:]),
            residues=(residue, request.topology.residues[1]),
            components=(component, request.topology.components[1]),
            topology_hash="",
        ).with_computed_hash()
        ligand_component = replace(
            request.assignment_components[1],
            atom_formal_charges_complete=False,
            charge_resolution_method="derived_from_residue_total:fixture",
        )
        validate_input(replace(
            request,
            topology=topology,
            assignment_components=(request.assignment_components[0], ligand_component),
            projection_hash="",
        ).with_computed_hash())

    def test_unknown_metal_formal_charge_is_rejected(self):
        request = _bonded_request()
        atom = replace(request.topology.atoms[0], formal_charge=None, formal_charge_known=False)
        residue = replace(request.topology.residues[0], atom_formal_charges_complete=False)
        component = replace(request.topology.components[0], atom_formal_charges_complete=False)
        topology = replace(
            request.topology,
            atoms=(atom, *request.topology.atoms[1:]),
            residues=(residue, request.topology.residues[1]),
            components=(component, request.topology.components[1]),
            topology_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "missing_metal_formal_charge"):
            validate_input(replace(request, topology=topology, projection_hash="").with_computed_hash())

    def test_inactive_canonical_cross_residue_edge_can_be_retained(self):
        request = _bonded_request()
        inactive = PreparedBond("canonical-fe-his", ("heme-fe", "his-ne2"), active=False, source="canonical")
        topology = replace(request.topology, bonds=(*request.topology.bonds, inactive), topology_hash="").with_computed_hash()
        updated = replace(request, topology=topology, projection_hash="").with_computed_hash()
        validate_input(updated)

    def test_contract_mappings_are_deeply_immutable(self):
        request = _bonded_request()
        with self.assertRaises(TypeError):
            request.partition_proofs[0].simulation_partitions["HEM:A:1"] = ()
        with self.assertRaises(TypeError):
            request.charge_contract.component_formal_charges["HEM-component"] = 1

    def test_base_force_field_assignment_rejects_metal(self):
        request = _bonded_request()
        bad_component = replace(request.assignment_components[0], base_force_field=True, provider="gaff2")
        bad = replace(request, assignment_components=(bad_component, request.assignment_components[1]), projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "metal_in_base_force_field"):
            validate_input(bad)

    def test_gaff_assignment_requires_complete_explicit_hydrogen_proof(self):
        request = _bonded_request()
        ligand = request.assignment_components[1]
        proof = replace(
            ligand.chemical_topology_proof,
            implicit_hydrogen_candidates_by_atom={"heme-na": (1,), "heme-c": (0,)},
            explicit_hydrogen_status="incomplete",
            proof_hash="",
        ).with_computed_hash()
        bad_ligand = replace(ligand, chemical_topology_proof=proof)
        bad = replace(
            request,
            assignment_components=(request.assignment_components[0], bad_ligand),
            projection_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "incomplete_explicit_hydrogen_topology"):
            validate_input(bad)

    def test_stale_topology_projection_revision_and_input_hash_are_rejected(self):
        request = _bonded_request()
        with self.assertRaisesRegex(ValidationError, "stale_topology_hash"):
            validate_input(replace(request, topology=replace(request.topology, topology_hash="stale")))
        with self.assertRaisesRegex(ValidationError, "stale_projection_hash"):
            validate_input(replace(request, projection_hash="stale"))
        with self.assertRaisesRegex(ValidationError, "stale_graph_revision"):
            validate_input(replace(request, graph_revision=8, projection_hash="").with_computed_hash())
        other_hash = hashlib.sha256(b"other").hexdigest()
        with self.assertRaisesRegex(ValidationError, "stale_input_hash"):
            validate_input(replace(request, input_hash=other_hash, projection_hash="").with_computed_hash())

    def test_disconnected_simulation_residue_is_rejected(self):
        request = _bonded_request()
        disconnected = replace(request.topology, bonds=(), topology_hash="").with_computed_hash()
        bad = replace(request, topology=disconnected, projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "disconnected_simulation_residue"):
            validate_input(bad)

    def test_cross_residue_edges_must_be_links(self):
        request = _bonded_request()
        cross_bond = PreparedBond("fe-his-bond", ("heme-fe", "his-ne2"), source="input")
        topology = replace(request.topology, bonds=(*request.topology.bonds, cross_bond), links=(), topology_hash="").with_computed_hash()
        bad = replace(request, topology=topology, projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "cross_residue_bond_misclassified"):
            validate_input(bad)

    def test_projection_and_assignment_atoms_must_match_exactly(self):
        request = _bonded_request()
        projection = replace(request.projections[1], atom_ids=("heme-na",))
        bad = replace(request, projections=(request.projections[0], projection), projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "projection_assignment_mismatch"):
            validate_input(bad)

    def test_overlapping_charge_scopes_are_rejected(self):
        request = _bonded_request()
        bad_charge = replace(request.charge_contract, fixed_atom_ids=("heme-fe",))
        bad = replace(request, charge_contract=bad_charge, projection_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "overlapping_charge_scopes"):
            validate_input(bad)

    def test_bonded_fit_input_round_trip_is_hash_closed(self):
        request = _bonded_request()
        spec = BondedMetalParameterSpec(
            schema_version=1,
            topology_hash=request.topology.topology_hash,
            metal_atoms=(MetalAtomParameterSpec(
                "heme-fe", "Fe", 2, "Fe-site", 1.2, 55.845, 0.01, 1.35,
            ),),
            donor_atom_types={"heme-na": "n-site"},
            parameter_source="unit-test:metal-spec",
            provenance={"protocol": "fixture"},
        ).with_computed_hash()
        fit_input = BondedFitInput(
            schema_version=1,
            metal_parameter_spec=spec,
            charge_artifacts=(),
            hessian_artifacts=(),
            force_method="empirical_zn_nos",
        ).with_computed_hash()
        payload = bonded_fit_input_dumps(fit_input)
        self.assertEqual(bonded_fit_input_loads(payload), fit_input)
        corrupted = json.loads(payload)
        corrupted["scale_factor"] = 2.0
        with self.assertRaisesRegex(ValidationError, "stale_bonded_fit_input_hash"):
            bonded_fit_input_loads(json.dumps(corrupted))

    def test_manual_bonded_uses_confirmed_current_geometry(self):
        topology = _bonded_request().topology
        reference_geometry = {
            "schema_version": 1,
            "graph_revision": topology.graph_revision,
            "input_hash": topology.input_hash,
            "coordinate_unit": "angstrom",
            "angle_unit": "rad",
            "geometry_source": "frozen_current_geometry",
            "selections": [{
                "selection_id": "site:heme-fe",
                "center_external_id": "heme-fe",
                "bonds": [
                    {
                        "edge_id": "fe-na",
                        "center_external_id": "heme-fe",
                        "neighbor_external_id": "heme-na",
                        "equilibrium": 1.9,
                        "unit": "angstrom",
                    },
                    {
                        "edge_id": "fe-his",
                        "center_external_id": "heme-fe",
                        "neighbor_external_id": "his-ne2",
                        "equilibrium": 2.1,
                        "unit": "angstrom",
                    },
                ],
                "angles": [{
                    "edge_id1": "fe-na",
                    "edge_id2": "fe-his",
                    "neighbor1_external_id": "heme-na",
                    "center_external_id": "heme-fe",
                    "neighbor2_external_id": "his-ne2",
                    "equilibrium": 1.5707963267948966,
                    "unit": "rad",
                }],
            }],
        }
        reference_geometry["artifact_hash"] = hashlib.sha256(
            json.dumps(
                reference_geometry,
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        terms, report = manual_bonded_terms(
            topology,
            bond_force_constant=100.0,
            angle_force_constant=20.0,
            reference_geometry_artifact=reference_geometry,
            site_force_constants={
                "site:heme-fe": {
                    "bond_force_constant": 111.0,
                    "angle_force_constant": 22.0,
                },
            },
        )
        bonds = [term for term in terms.values() if term["kind"] == "bond"]
        angles = [term for term in terms.values() if term["kind"] == "angle"]
        self.assertEqual(len(bonds), 2)
        self.assertEqual(len(angles), 1)
        self.assertEqual(
            sorted(term["parameters"]["equilibrium"] for term in bonds),
            [1.9, 2.1],
        )
        self.assertAlmostEqual(angles[0]["parameters"]["equilibrium"], 1.57079633)
        self.assertTrue(all(term["parameters"]["k"] == 111.0 for term in bonds))
        self.assertEqual(angles[0]["parameters"]["k"], 22.0)
        self.assertEqual(report["potential_convention"], "E=k*delta^2")

        request = _bonded_request()
        spec = BondedMetalParameterSpec(
            schema_version=1,
            topology_hash=request.topology.topology_hash,
            metal_atoms=(MetalAtomParameterSpec(
                "heme-fe", "Fe", 2, "Fe-site", 1.2, 55.845, 0.01, 1.35,
            ),),
            donor_atom_types={"heme-na": "n-site", "his-ne2": "n-site"},
            parameter_source="manual-bonded-test",
            provenance={"protocol": "fixture"},
        ).with_computed_hash()
        fit_input = BondedFitInput(
            schema_version=1,
            metal_parameter_spec=spec,
            charge_artifacts=(),
            hessian_artifacts=(),
            force_method="manual_bonded",
            manual_bond_force_constant=100.0,
            manual_angle_force_constant=20.0,
            manual_site_force_constants={
                "site:heme-fe": {
                    "bond_force_constant": 111.0,
                    "angle_force_constant": 22.0,
                },
            },
            equilibrium_geometry_source="frozen_current_geometry",
            reference_geometry_artifact=reference_geometry,
        ).with_computed_hash()
        self.assertEqual(
            bonded_fit_input_loads(bonded_fit_input_dumps(fit_input)),
            fit_input,
        )

    def test_empirical_registry_is_hash_closed_and_citation_bearing(self):
        registry_id = "xponge-zn-nos-experimental-v1"
        entry = EMPIRICAL_REGISTRY[registry_id]
        validate_empirical_registry_entry(entry)
        descriptor = empirical_registry_descriptor(registry_id)
        self.assertEqual(descriptor["registry_id"], registry_id)
        self.assertEqual(descriptor["citations"][0]["doi"], "10.1021/ct1002626")
        self.assertEqual(descriptor["redistribution_status"], "internal-test-only")
        with self.assertRaisesRegex(ValidationError, "stale_empirical_parameter_hash"):
            validate_empirical_registry_entry(replace(entry, parameter_hash="forged"))

    def test_empirical_registry_exactly_matches_topology_identity(self):
        topology = _bonded_request().topology
        atoms = tuple(
            replace(atom, element="Zn", formal_charge=2, formal_charge_known=True)
            if atom.external_id == "heme-fe"
            else atom
            for atom in topology.atoms
        )
        zinc_topology = replace(topology, atoms=atoms, topology_hash="").with_computed_hash()
        match = resolve_empirical_registry(
            zinc_topology,
            registry_id="xponge-zn-nos-experimental-v1",
        )
        self.assertEqual(match.centers[0].element, "Zn")
        self.assertEqual(match.centers[0].oxidation_state, 2)
        self.assertGreater(match.centers[0].coordination_number, 0)
        with self.assertRaisesRegex(ValidationError, "empirical_element_mismatch"):
            resolve_empirical_registry(
                topology,
                registry_id="xponge-zn-nos-experimental-v1",
            )
        wrong_charge = replace(
            zinc_topology,
            atoms=tuple(
                replace(atom, formal_charge=1)
                if atom.external_id == "heme-fe"
                else atom
                for atom in zinc_topology.atoms
            ),
            topology_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "empirical_oxidation_state_mismatch"):
            resolve_empirical_registry(
                wrong_charge,
                registry_id="xponge-zn-nos-experimental-v1",
            )
        with self.assertRaisesRegex(ValidationError, "empirical_geometry_mismatch"):
            resolve_empirical_registry(
                zinc_topology,
                registry_id="xponge-zn-nos-experimental-v1",
                geometry="tetrahedral",
            )

        original = EMPIRICAL_REGISTRY["xponge-zn-nos-experimental-v1"]
        restricted = replace(
            original,
            registry_id="unit-test-zn-oo-v1",
            applicability=replace(
                original.applicability,
                donor_patterns=(("O", "O"),),
                compatible_base_force_fields=("gaff2",),
                compatible_water_models=("tip3p",),
            ),
            source_hash="",
            parameter_hash="",
        ).with_computed_hashes()
        EMPIRICAL_REGISTRY[restricted.registry_id] = restricted
        try:
            with self.assertRaisesRegex(ValidationError, "empirical_base_force_field_mismatch"):
                resolve_empirical_registry(
                    zinc_topology,
                    registry_id=restricted.registry_id,
                    base_force_field="gaff",
                    water_model="tip3p",
                )
            with self.assertRaisesRegex(ValidationError, "empirical_donor_pattern_mismatch"):
                resolve_empirical_registry(
                    zinc_topology,
                    registry_id=restricted.registry_id,
                    base_force_field="gaff2",
                    water_model="tip3p",
                )
        finally:
            del EMPIRICAL_REGISTRY[restricted.registry_id]

    def test_empirical_registry_schema_is_not_zinc_specific(self):
        entry = EmpiricalRegistryEntry(
            registry_id="unit-test-cu-v1",
            version="1",
            status="experimental",
            applicability=EmpiricalApplicability(
                interaction_model="bonded",
                element="Cu",
                oxidation_state=2,
                coordination_numbers=(4,),
                donor_elements=("N",),
                donor_patterns=(("N", "N", "N", "N"),),
                geometries=("square_planar",),
                compatible_base_force_fields=("gaff2",),
                compatible_water_models=("tip3p",),
            ),
            parameters={
                "bond_tables": {"N": ((2.0, 100.0),)},
                "angle_force_constant_default": 50.0,
                "angle_force_constant_with_sulfur": 70.0,
            },
            units={
                "bond_distance": "angstrom",
                "bond_force_constant": "kcal/mol/angstrom^2",
                "angle_force_constant": "kcal/mol/rad^2",
            },
            source="unit test",
            license="test-only",
            redistribution_status="test-only",
            citations=(EmpiricalCitation(
                "Unit-test citation", ("Test Author",), "Test Journal", 2026, "10.0000/test",
            ),),
            validation_cases=("synthetic CuN4",),
        ).with_computed_hashes()
        validate_empirical_registry_entry(entry)
        self.assertEqual(entry.applicability.element, "Cu")

    def test_empirical_registry_fit_input_requires_explicit_registry_id(self):
        request = _bonded_request()
        spec = BondedMetalParameterSpec(
            schema_version=1,
            topology_hash=request.topology.topology_hash,
            metal_atoms=(MetalAtomParameterSpec(
                "heme-fe", "Fe", 2, "Fe-site", 1.2, 55.845, 0.01, 1.35,
            ),),
            donor_atom_types={"heme-na": "n-site"},
            parameter_source="unit-test:metal-spec",
            provenance={"protocol": "fixture"},
        ).with_computed_hash()
        fit_input = BondedFitInput(
            schema_version=1,
            metal_parameter_spec=spec,
            charge_artifacts=(),
            hessian_artifacts=(),
            force_method="empirical_registry",
            empirical_registry_id="xponge-zn-nos-experimental-v1",
        ).with_computed_hash()
        self.assertEqual(bonded_fit_input_loads(bonded_fit_input_dumps(fit_input)), fit_input)
        missing = replace(fit_input, empirical_registry_id="", fit_input_hash="").with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "missing_empirical_registry_id"):
            bonded_fit_input_dumps(missing)

    def test_parameterization_result_round_trip_is_strict_and_hash_closed(self):
        request = _bonded_request()
        result = _overlay_result(request)
        validate_result(request, result)
        payload = parameterization_result_dumps(request, result)
        restored = parameterization_result_loads(payload, request)
        self.assertEqual(restored, result)
        self.assertEqual(parameterization_result_dumps(request, restored), payload)

    def test_parameterization_result_rejects_hash_closed_term_outside_locked_topology(self):
        request = _bonded_request()
        result = _overlay_result(request)
        terms = dict(result.bonded_overlay.terms)
        terms["forged"] = {
            "kind": "bond",
            "atom_ids": ("heme-fe", "heme-c"),
            "parameters": {"k": 1.0, "equilibrium": 2.0},
            "source": "forged",
        }
        corrupted = replace(
            result,
            bonded_overlay=replace(result.bonded_overlay, terms=terms),
            result_hash="",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "bonded_parameter_topology_mismatch"):
            validate_result(request, corrupted)


if __name__ == "__main__":
    unittest.main()
