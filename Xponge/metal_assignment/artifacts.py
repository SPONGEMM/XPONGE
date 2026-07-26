"""Immutable, source-neutral structural and derived-model artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .contracts import (
    ElectronicState,
    MetalAssignmentInput,
    PreparedChemicalTopology,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _sha256,
    _strict_object,
    _strings,
    validate_electronic_state,
)


ATOMIC_CHARGE_ROLES = frozenset({"absent", "initial", "fixed", "reference", "fitted"})


@dataclass(frozen=True, slots=True)
class StructuralAtomMapping:
    artifact_id: str
    artifact_atom_id: str
    mol2_serial: int
    external_id: str
    canonical_atom_id: int
    element: str
    chemical_component_id: str
    canonical_residue_id: str
    parameterization_component_id: str
    simulation_component_id: str
    simulation_residue_id: str


@dataclass(frozen=True, slots=True)
class StructuralArtifact:
    external_id: str
    purpose: str
    classification: str
    provider: str
    base_force_field: bool
    net_formal_charge: int | None
    charge_resolution_method: str
    atom_formal_charges_complete: bool
    coordinate_unit: str
    atomic_charge_role: str
    atom_ids: tuple[str, ...]
    active_edge_ids: tuple[str, ...]
    mapping: tuple[StructuralAtomMapping, ...]
    mol2_text: str
    artifact_hash: str

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "external_id": self.external_id,
            "purpose": self.purpose,
            "classification": self.classification,
            "provider": self.provider,
            "base_force_field": self.base_force_field,
            "net_formal_charge": self.net_formal_charge,
            "charge_resolution_method": self.charge_resolution_method,
            "atom_formal_charges_complete": self.atom_formal_charges_complete,
            "coordinate_unit": self.coordinate_unit,
            "atomic_charge_role": self.atomic_charge_role,
            "atom_ids": _canonicalize(self.atom_ids),
            "active_edge_ids": _canonicalize(self.active_edge_ids),
            "mapping": _canonicalize(self.mapping),
            "mol2_text": self.mol2_text,
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())


@dataclass(frozen=True, slots=True)
class StructuralArtifactBundle:
    schema_version: int
    graph_revision: int
    input_hash: str
    coordinate_unit: str
    prepared_system: StructuralArtifact
    assignment_components: tuple[StructuralArtifact, ...]
    bundle_hash: str

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "graph_revision": self.graph_revision,
            "input_hash": self.input_hash,
            "coordinate_unit": self.coordinate_unit,
            "prepared_system": _canonicalize(self.prepared_system),
            "assignment_components": _canonicalize(self.assignment_components),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())


@dataclass(frozen=True, slots=True)
class ModelAtom:
    model_atom_id: str
    external_id: str | None
    canonical_atom_id: int
    mol2_serial: int
    element: str
    coordinates: tuple[float, float, float]
    role: str
    canonical_residue_id: str
    chemical_component_id: str
    partial_charge: float | None
    cap_parent_external_id: str
    cut_edge_id: str
    geometry_provenance: str
    charge_projection_group: str


@dataclass(frozen=True, slots=True)
class ModelBond:
    external_id: str
    model_atom_ids: tuple[str, str]
    order: float
    semantic: str
    source: str


@dataclass(frozen=True, slots=True)
class ModelLink:
    external_id: str
    model_atom_ids: tuple[str, str]
    kind: str
    source: str


@dataclass(frozen=True, slots=True)
class ModelCutEdge:
    external_id: str
    retained_external_id: str
    excluded_external_id: str
    semantic: str
    source: str


@dataclass(frozen=True, slots=True)
class ModelElectronicState:
    selection_id: str
    net_charge: int
    spin_multiplicity: int
    source: str


@dataclass(frozen=True, slots=True)
class DerivedModel:
    external_id: str
    site_id: str
    purpose: str
    coordinate_unit: str
    atomic_charge_role: str
    electronic_state: ModelElectronicState | None
    atoms: tuple[ModelAtom, ...]
    bonds: tuple[ModelBond, ...]
    links: tuple[ModelLink, ...]
    cut_edges: tuple[ModelCutEdge, ...]
    charge_accounting: Mapping[str, Any]
    mol2_text: str
    model_hash: str
    capped_model_manifest: Mapping[str, Any] | None = None

    def canonical_payload(self) -> dict[str, Any]:
        payload = {
            "external_id": self.external_id,
            "site_id": self.site_id,
            "purpose": self.purpose,
            "coordinate_unit": self.coordinate_unit,
            "atomic_charge_role": self.atomic_charge_role,
            "electronic_state": _canonicalize(self.electronic_state),
            "atoms": _canonicalize(self.atoms),
            "bonds": _canonicalize(self.bonds),
            "links": _canonicalize(self.links),
            "cut_edges": _canonicalize(self.cut_edges),
            "charge_accounting": _canonicalize(self.charge_accounting),
            "mol2_text": self.mol2_text,
        }
        if self.capped_model_manifest is not None:
            payload["capped_model_manifest"] = _canonicalize(self.capped_model_manifest)
        return payload

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())


@dataclass(frozen=True, slots=True)
class ModelMappingRecord:
    site_id: str
    model_id: str
    model_atom_id: str
    mol2_serial: int
    external_id: str | None
    canonical_atom_id: int
    element: str
    role: str
    cap_parent_external_id: str
    cut_edge_id: str


@dataclass(frozen=True, slots=True)
class DerivedModelMapping:
    records: tuple[ModelMappingRecord, ...]
    mapping_hash: str

    def computed_hash(self) -> str:
        return _sha256({"records": _canonicalize(self.records)})


@dataclass(frozen=True, slots=True)
class DerivedSite:
    external_id: str
    metal_atom_ids: tuple[str, ...]
    canonical_residue_ids: tuple[str, ...]
    atom_ids: tuple[str, ...]
    interaction_edge_ids: tuple[str, ...]
    model_ids: tuple[str, ...]
    multi_metal: bool


@dataclass(frozen=True, slots=True)
class DerivedModelBundle:
    schema_version: int
    graph_revision: int
    input_hash: str
    interaction_model: str
    coordinate_unit: str
    models_generated: bool
    skipped_reason: str
    protocol: Mapping[str, Any]
    source_atom_ids: tuple[str, ...]
    active_metal_interaction_edge_ids: tuple[str, ...]
    removed_interaction_edge_ids: tuple[str, ...]
    sites: tuple[DerivedSite, ...]
    models: tuple[DerivedModel, ...]
    mapping: DerivedModelMapping
    bundle_hash: str

    def __post_init__(self) -> None:
        _freeze_field(self, "protocol")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "graph_revision": self.graph_revision,
            "input_hash": self.input_hash,
            "interaction_model": self.interaction_model,
            "coordinate_unit": self.coordinate_unit,
            "models_generated": self.models_generated,
            "skipped_reason": self.skipped_reason,
            "protocol": _canonicalize(self.protocol),
            "source_atom_ids": _canonicalize(self.source_atom_ids),
            "active_metal_interaction_edge_ids": _canonicalize(self.active_metal_interaction_edge_ids),
            "removed_interaction_edge_ids": _canonicalize(self.removed_interaction_edge_ids),
            "sites": _canonicalize(self.sites),
            "models": _canonicalize(self.models),
            "mapping": _canonicalize(self.mapping),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())


@dataclass(frozen=True, slots=True)
class ClosedAtomMapping:
    schema_version: int
    graph_revision: int
    input_hash: str
    coordinate_unit: str
    prepared_topology_hash: str
    parent_records: tuple[StructuralAtomMapping, ...]
    structural_records: tuple[StructuralAtomMapping, ...]
    model_records: tuple[ModelMappingRecord, ...]
    mapping_hash: str

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "graph_revision": self.graph_revision,
            "input_hash": self.input_hash,
            "coordinate_unit": self.coordinate_unit,
            "prepared_topology_hash": self.prepared_topology_hash,
            "parent_records": _canonicalize(self.parent_records),
            "structural_records": _canonicalize(self.structural_records),
            "model_records": _canonicalize(self.model_records),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())


@dataclass(frozen=True, slots=True)
class PreparedArtifactPackage:
    structural_artifacts: StructuralArtifactBundle
    atom_mapping: ClosedAtomMapping
    derived_models: DerivedModelBundle


def _validate_hash(actual: str, expected: str, code: str, path: str) -> None:
    if not actual:
        raise ValidationError(f"missing_{code}", "computed hash is required", path)
    if actual != expected:
        raise ValidationError(f"stale_{code}", "hash does not match canonical payload", path)


def _validate_mol2(text: str, atom_count: int, edge_count: int, charge_role: str, path: str) -> None:
    lines = text.splitlines()
    if len(lines) < 6 or lines[0] != "@<TRIPOS>MOLECULE":
        raise ValidationError("invalid_mol2", "missing MOLECULE header", path)
    try:
        counts = [int(value) for value in lines[2].split()[:3]]
    except (ValueError, IndexError) as exc:
        raise ValidationError("invalid_mol2", "invalid counts line", path) from exc
    if len(counts) != 3 or counts[0] != atom_count or counts[1] != edge_count:
        raise ValidationError("mol2_count_mismatch", "MOL2 counts do not match artifact", path)
    charge_type = lines[4].strip()
    if (charge_role == "absent") != (charge_type == "NO_CHARGES"):
        raise ValidationError(
            "mol2_charge_role_mismatch",
            "NO_CHARGES must be used exactly when atomic charges are absent",
            path,
        )
    section_positions = {
        line: index for index, line in enumerate(lines) if line.startswith("@<TRIPOS>")
    }
    if "@<TRIPOS>ATOM" not in section_positions or "@<TRIPOS>BOND" not in section_positions:
        raise ValidationError("invalid_mol2", "ATOM and BOND sections are required", path)

    def section_rows(name: str) -> list[str]:
        start = section_positions[name] + 1
        end = next(
            (index for index in range(start, len(lines)) if lines[index].startswith("@<TRIPOS>")),
            len(lines),
        )
        return [line for line in lines[start:end] if line.strip() and not line.lstrip().startswith("#")]

    atom_rows = section_rows("@<TRIPOS>ATOM")
    bond_rows = section_rows("@<TRIPOS>BOND")
    if len(atom_rows) != atom_count or len(bond_rows) != edge_count:
        raise ValidationError("mol2_record_count_mismatch", "MOL2 section records do not match counts", path)
    atom_serials: list[int] = []
    for row in atom_rows:
        fields = row.split()
        if len(fields) < 6:
            raise ValidationError("invalid_mol2_atom", row, path)
        try:
            atom_serials.append(int(fields[0]))
            float(fields[2])
            float(fields[3])
            float(fields[4])
        except ValueError as exc:
            raise ValidationError("invalid_mol2_atom", row, path) from exc
    if atom_serials != list(range(1, atom_count + 1)):
        raise ValidationError("noncontiguous_mol2_atom_serial", path)
    bond_serials: list[int] = []
    atom_serial_set = set(atom_serials)
    for row in bond_rows:
        fields = row.split()
        if len(fields) < 4:
            raise ValidationError("invalid_mol2_bond", row, path)
        try:
            serial, atom1, atom2 = map(int, fields[:3])
        except ValueError as exc:
            raise ValidationError("invalid_mol2_bond", row, path) from exc
        if atom1 == atom2 or atom1 not in atom_serial_set or atom2 not in atom_serial_set:
            raise ValidationError("invalid_mol2_bond_endpoint", row, path)
        bond_serials.append(serial)
    if bond_serials != list(range(1, edge_count + 1)):
        raise ValidationError("noncontiguous_mol2_bond_serial", path)


def _validate_structural_mapping(
    artifact: StructuralArtifact,
    topology: PreparedChemicalTopology,
    parameterization_by_atom: Mapping[str, str],
) -> None:
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    records = artifact.mapping
    if len(records) != len(artifact.atom_ids):
        raise ValidationError("incomplete_artifact_mapping", artifact.external_id)
    if tuple(record.external_id for record in records) != artifact.atom_ids:
        raise ValidationError("artifact_mapping_order_mismatch", artifact.external_id)
    if tuple(record.mol2_serial for record in records) != tuple(range(1, len(records) + 1)):
        raise ValidationError("noncontiguous_mol2_serial", artifact.external_id)
    if len({record.artifact_atom_id for record in records}) != len(records):
        raise ValidationError("duplicate_artifact_atom_id", artifact.external_id)
    for record in records:
        atom = atom_by_id.get(record.external_id)
        if atom is None:
            raise ValidationError("mapping_unknown_atom", record.external_id, artifact.external_id)
        expected = (
            artifact.external_id,
            atom.canonical_atom_id,
            atom.element,
            atom.chemical_component_id,
            atom.canonical_residue_id,
            parameterization_by_atom[atom.external_id],
            atom.simulation_component_id,
            atom.simulation_residue_id,
        )
        actual = (
            record.artifact_id,
            record.canonical_atom_id,
            record.element,
            record.chemical_component_id,
            record.canonical_residue_id,
            record.parameterization_component_id,
            record.simulation_component_id,
            record.simulation_residue_id,
        )
        if actual != expected:
            raise ValidationError("artifact_mapping_identity_mismatch", record.external_id, artifact.external_id)


def validate_structural_artifacts(request: MetalAssignmentInput, bundle: StructuralArtifactBundle) -> None:
    topology = request.topology
    if (
        bundle.schema_version != request.schema_version
        or bundle.graph_revision != request.graph_revision
        or bundle.input_hash != request.input_hash
        or bundle.coordinate_unit != topology.coordinate_unit
    ):
        raise ValidationError("artifact_identity_mismatch", "structural bundle belongs to another input")
    _validate_hash(bundle.bundle_hash, bundle.computed_hash(), "structural_bundle_hash", "structural_artifacts")
    parameterization_by_atom: dict[str, str] = {}
    for projection in request.projections:
        for atom_id in projection.atom_ids:
            parameterization_by_atom[atom_id] = projection.parameterization_component_id
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    edge_by_id = {edge.external_id: edge for edge in (*topology.bonds, *topology.links)}
    expected_all_ids = tuple(atom.external_id for atom in topology.atoms)
    prepared = bundle.prepared_system
    if prepared.purpose != "prepared_system" or prepared.atom_ids != expected_all_ids:
        raise ValidationError("prepared_system_atom_mismatch", "prepared system must cover topology in order")
    expected_prepared_edges = {
        edge.external_id for edge in (*topology.bonds, *topology.links) if edge.active
    }
    if len(prepared.active_edge_ids) != len(set(prepared.active_edge_ids)):
        raise ValidationError("duplicate_artifact_edge_id", prepared.external_id)
    if set(prepared.active_edge_ids) != expected_prepared_edges:
        raise ValidationError("prepared_system_edge_mismatch", "active structural edges differ from topology")

    component_by_id = {component.external_id: component for component in request.assignment_components}
    artifact_by_id = {artifact.external_id: artifact for artifact in bundle.assignment_components}
    if len(artifact_by_id) != len(bundle.assignment_components) or set(artifact_by_id) != set(component_by_id):
        raise ValidationError("assignment_artifact_mismatch", "one artifact is required per assignment component")
    assignment_atoms: list[str] = []
    for component_id, component in component_by_id.items():
        artifact = artifact_by_id[component_id]
        if artifact.purpose != "assignment_component" or artifact.atom_ids != component.atom_ids:
            raise ValidationError("assignment_artifact_atom_mismatch", component_id)
        if (
            artifact.classification != component.classification
            or artifact.provider != component.provider
            or artifact.base_force_field != component.base_force_field
            or artifact.net_formal_charge != component.net_formal_charge
            or artifact.charge_resolution_method != component.charge_resolution_method
            or artifact.atom_formal_charges_complete != component.atom_formal_charges_complete
        ):
            raise ValidationError("assignment_artifact_metadata_mismatch", component_id)
        if artifact.base_force_field and any(atom_by_id[atom_id].is_metal for atom_id in artifact.atom_ids):
            raise ValidationError("metal_in_base_force_field", component_id)
        expected_internal_edges = {
            edge_id for edge_id, edge in edge_by_id.items()
            if edge.active and set(edge.atom_ids) <= set(artifact.atom_ids)
        }
        if len(artifact.active_edge_ids) != len(set(artifact.active_edge_ids)):
            raise ValidationError("duplicate_artifact_edge_id", artifact.external_id)
        if set(artifact.active_edge_ids) != expected_internal_edges:
            raise ValidationError("assignment_artifact_edge_mismatch", component_id)
        assignment_atoms.extend(artifact.atom_ids)
    if len(assignment_atoms) != len(set(assignment_atoms)) or set(assignment_atoms) != set(expected_all_ids):
        raise ValidationError("assignment_artifact_partition_mismatch", "assignment artifacts must partition atoms")

    for artifact in (prepared, *bundle.assignment_components):
        if artifact.coordinate_unit != topology.coordinate_unit or artifact.atomic_charge_role not in ATOMIC_CHARGE_ROLES:
            raise ValidationError("invalid_artifact_metadata", artifact.external_id)
        _validate_hash(artifact.artifact_hash, artifact.computed_hash(), "artifact_hash", artifact.external_id)
        _validate_structural_mapping(artifact, topology, parameterization_by_atom)
        _validate_mol2(
            artifact.mol2_text,
            len(artifact.atom_ids),
            len(artifact.active_edge_ids),
            artifact.atomic_charge_role,
            artifact.external_id,
        )


def _validate_derived_sites(request: MetalAssignmentInput, bundle: DerivedModelBundle) -> None:
    topology = request.topology
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    edge_by_id = {edge.external_id: edge for edge in (*topology.bonds, *topology.links)}
    metal_ids = {atom.external_id for atom in topology.atoms if atom.is_metal}
    all_metal_edge_ids = {
        edge.external_id
        for edge in (*topology.bonds, *topology.links)
        if any(atom_by_id[atom_id].is_metal for atom_id in edge.atom_ids)
    }
    active_metal_edge_ids = {edge_id for edge_id in all_metal_edge_ids if edge_by_id[edge_id].active}
    removed_edge_ids = {
        edge_id for proof in request.partition_proofs for edge_id in proof.removed_interaction_edge_ids
    }
    for values, code in (
        (bundle.active_metal_interaction_edge_ids, "duplicate_active_metal_interaction_edge"),
        (bundle.removed_interaction_edge_ids, "duplicate_removed_interaction_edge"),
    ):
        if len(values) != len(set(values)):
            raise ValidationError(code, "interaction edge IDs must be unique")
    if set(bundle.active_metal_interaction_edge_ids) != active_metal_edge_ids:
        raise ValidationError("active_metal_interaction_edge_mismatch", "active metal edges differ from topology")
    if set(bundle.removed_interaction_edge_ids) != removed_edge_ids:
        raise ValidationError("removed_interaction_edge_mismatch", "removed edges differ from partition proofs")
    if request.interaction_model == "nonbonded_12_6" and not bundle.sites:
        return

    site_ids: set[str] = set()
    owned_metals: set[str] = set()
    owned_residues: set[str] = set()
    owned_site_edges: set[str] = set()
    for site in bundle.sites:
        if not site.external_id or site.external_id in site_ids:
            raise ValidationError("duplicate_site_id", site.external_id)
        site_ids.add(site.external_id)
        for values, code in (
            (site.metal_atom_ids, "duplicate_site_metal"),
            (site.canonical_residue_ids, "duplicate_site_residue"),
            (site.atom_ids, "duplicate_site_atom"),
            (site.interaction_edge_ids, "duplicate_site_interaction_edge"),
            (site.model_ids, "duplicate_site_model"),
        ):
            if len(values) != len(set(values)):
                raise ValidationError(code, site.external_id)
        site_metals = set(site.metal_atom_ids)
        site_residues = set(site.canonical_residue_ids)
        site_atoms = set(site.atom_ids)
        site_edges = set(site.interaction_edge_ids)
        if not site_metals or not site_metals <= metal_ids or site.multi_metal != (len(site_metals) > 1):
            raise ValidationError("invalid_site_metal_membership", site.external_id)
        expected_atoms = {
            atom.external_id for atom in topology.atoms if atom.canonical_residue_id in site_residues
        }
        if site_atoms != expected_atoms or not site_metals <= site_atoms:
            raise ValidationError("site_atom_membership_mismatch", site.external_id)
        if {atom_by_id[atom_id].canonical_residue_id for atom_id in site_atoms} != site_residues:
            raise ValidationError("site_residue_membership_mismatch", site.external_id)
        if not site_edges <= all_metal_edge_ids:
            raise ValidationError("site_interaction_edge_unknown", site.external_id)
        for edge_id in site_edges:
            edge = edge_by_id[edge_id]
            if not set(edge.atom_ids) <= site_atoms or not (set(edge.atom_ids) & site_metals):
                raise ValidationError("site_interaction_edge_mismatch", edge_id, site.external_id)
        if owned_metals & site_metals or owned_residues & site_residues:
            raise ValidationError("overlapping_metal_site", site.external_id)
        owned_metals.update(site_metals)
        owned_residues.update(site_residues)
        owned_site_edges.update(site_edges)
    if owned_metals != metal_ids or owned_site_edges != all_metal_edge_ids:
        raise ValidationError("incomplete_metal_site_coverage", "sites must cover every metal and metal interaction")


def _validate_model_topology(
    model: DerivedModel,
    site: DerivedSite,
    topology: PreparedChemicalTopology,
) -> None:
    topology_atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    topology_bond_by_id = {edge.external_id: edge for edge in topology.bonds}
    topology_link_by_id = {edge.external_id: edge for edge in topology.links}
    model_atom_by_id = {atom.model_atom_id: atom for atom in model.atoms}
    core_by_external = {
        str(atom.external_id): atom for atom in model.atoms if atom.role != "cap"
    }
    if set(core_by_external) != set(site.atom_ids) or len(core_by_external) != len(site.atom_ids):
        raise ValidationError("model_core_atom_mismatch", model.external_id)
    edge_ids = [edge.external_id for edge in (*model.bonds, *model.links, *model.cut_edges)]
    if len(edge_ids) != len(set(edge_ids)):
        raise ValidationError("duplicate_model_edge_id", model.external_id)

    expected_internal_bonds: set[str] = set()
    expected_internal_links: set[str] = set()
    expected_cuts: set[str] = set()
    core_ids = set(core_by_external)
    for edge in (*topology.bonds, *topology.links):
        if not edge.active:
            continue
        selected = [atom_id in core_ids for atom_id in edge.atom_ids]
        if all(selected):
            (expected_internal_bonds if edge.external_id in topology_bond_by_id else expected_internal_links).add(
                edge.external_id
            )
        elif any(selected):
            expected_cuts.add(edge.external_id)
    if {edge.external_id for edge in model.bonds} != expected_internal_bonds:
        raise ValidationError("model_bond_set_mismatch", model.external_id)
    if {edge.external_id for edge in model.links} != expected_internal_links:
        raise ValidationError("model_link_set_mismatch", model.external_id)
    if {edge.external_id for edge in model.cut_edges} != expected_cuts:
        raise ValidationError("model_cut_edge_set_mismatch", model.external_id)

    for edge, source_edge, semantic in (
        *(
            (edge, topology_bond_by_id.get(edge.external_id), edge.semantic)
            for edge in model.bonds
        ),
        *(
            (edge, topology_link_by_id.get(edge.external_id), edge.kind)
            for edge in model.links
        ),
    ):
        if source_edge is None or len(edge.model_atom_ids) != 2 or len(set(edge.model_atom_ids)) != 2:
            raise ValidationError("invalid_model_edge", edge.external_id, model.external_id)
        if any(atom_id not in model_atom_by_id or model_atom_by_id[atom_id].role == "cap" for atom_id in edge.model_atom_ids):
            raise ValidationError("model_edge_unknown_atom", edge.external_id, model.external_id)
        external_endpoints = {str(model_atom_by_id[atom_id].external_id) for atom_id in edge.model_atom_ids}
        if external_endpoints != set(source_edge.atom_ids) or semantic != (
            source_edge.semantic if edge.external_id in topology_bond_by_id else source_edge.kind
        ) or edge.source != source_edge.source:
            raise ValidationError("model_edge_identity_mismatch", edge.external_id, model.external_id)
        if edge.external_id in topology_bond_by_id and edge.order != source_edge.order:
            raise ValidationError("model_bond_order_mismatch", edge.external_id, model.external_id)

    cut_by_id = {edge.external_id: edge for edge in model.cut_edges}
    for cut in model.cut_edges:
        source_edge = topology_bond_by_id.get(cut.external_id) or topology_link_by_id.get(cut.external_id)
        if source_edge is None or {cut.retained_external_id, cut.excluded_external_id} != set(source_edge.atom_ids):
            raise ValidationError("model_cut_edge_identity_mismatch", cut.external_id, model.external_id)
        if cut.retained_external_id not in core_ids or cut.excluded_external_id in core_ids:
            raise ValidationError("model_cut_edge_orientation_mismatch", cut.external_id, model.external_id)
        expected_semantic = getattr(source_edge, "semantic", getattr(source_edge, "kind", ""))
        if cut.semantic != expected_semantic or cut.source != source_edge.source:
            raise ValidationError("model_cut_edge_metadata_mismatch", cut.external_id, model.external_id)

    caps = [atom for atom in model.atoms if atom.role == "cap"]
    if model.purpose == "site" and caps:
        raise ValidationError("unexpected_site_model_cap", model.external_id)
    if model.purpose in {"small", "large"} and {atom.cut_edge_id for atom in caps} != set(cut_by_id):
        raise ValidationError("incomplete_model_caps", model.external_id)
    for cap in caps:
        cut = cut_by_id.get(cap.cut_edge_id)
        if cut is None or cap.cap_parent_external_id != cut.retained_external_id:
            raise ValidationError("cap_cut_edge_mismatch", cap.model_atom_id, model.external_id)


def validate_derived_models(request: MetalAssignmentInput, bundle: DerivedModelBundle) -> None:
    topology = request.topology
    if (
        bundle.schema_version != request.schema_version
        or bundle.graph_revision != request.graph_revision
        or bundle.input_hash != request.input_hash
        or bundle.interaction_model != request.interaction_model
        or bundle.coordinate_unit != topology.coordinate_unit
    ):
        raise ValidationError("model_bundle_identity_mismatch", "derived models belong to another input")
    _validate_hash(bundle.bundle_hash, bundle.computed_hash(), "model_bundle_hash", "derived_models")
    _validate_hash(bundle.mapping.mapping_hash, bundle.mapping.computed_hash(), "model_mapping_hash", "derived_models.mapping")
    topology_atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    if tuple(bundle.source_atom_ids) != tuple(atom.external_id for atom in topology.atoms):
        raise ValidationError("model_source_atom_mismatch", "model source atoms must match prepared topology")
    _validate_derived_sites(request, bundle)
    if request.interaction_model == "nonbonded_12_6":
        if (
            bundle.models_generated
            or bundle.models
            or bundle.mapping.records
            or bundle.skipped_reason != "nonbonded_12_6_fast_path"
            or any(site.model_ids for site in bundle.sites)
        ):
            raise ValidationError("nonbonded_model_generation", "nonbonded fast path cannot contain derived models")
        return
    if not bundle.models_generated or bundle.skipped_reason:
        raise ValidationError("missing_bonded_models", "bonded inputs require generated models")
    model_by_id = {model.external_id: model for model in bundle.models}
    if len(model_by_id) != len(bundle.models):
        raise ValidationError("duplicate_model_id", "derived model IDs must be unique")
    site_by_id = {site.external_id: site for site in bundle.sites}
    if set(model.site_id for model in bundle.models) != set(site_by_id):
        raise ValidationError("orphan_model_site", "every bonded site must own models")
    for site in bundle.sites:
        if not site.model_ids or set(site.model_ids) != {
            model.external_id for model in bundle.models if model.site_id == site.external_id
        }:
            raise ValidationError("site_model_mismatch", site.external_id)
    expected_mapping: dict[tuple[str, str], tuple[Any, ...]] = {}
    for model in bundle.models:
        if model.purpose not in {"site", "small", "large"} or model.atomic_charge_role not in ATOMIC_CHARGE_ROLES:
            raise ValidationError("invalid_model_metadata", model.external_id)
        if model.purpose == "site":
            if model.electronic_state is not None:
                raise ValidationError("unexpected_model_electronic_state", model.external_id)
        else:
            state = model.electronic_state
            if state is None:
                raise ValidationError("missing_model_electronic_state", model.external_id)
            if not state.selection_id or not state.source:
                raise ValidationError("missing_model_electronic_state_source", model.external_id)
            validate_electronic_state(
                ElectronicState(state.net_charge, state.spin_multiplicity),
                [atom.element for atom in model.atoms],
            )
        proof = model.charge_accounting
        if not isinstance(proof, Mapping):
            raise ValidationError("invalid_model_charge_accounting", model.external_id)
        legacy_proof_fields = {
            "schema_version", "graph_revision", "selection_id", "net_charge",
            "complete", "contributions", "proof_hash",
        }
        constrained_proof_fields = legacy_proof_fields | {
            "core_target_charge", "cap_target_charge", "model_target_charge",
            "constraint_groups",
        }
        if frozenset(proof) not in {
            frozenset(legacy_proof_fields),
            frozenset(constrained_proof_fields),
        }:
            raise ValidationError("invalid_model_charge_accounting_fields", model.external_id)
        proof_payload = dict(proof)
        proof_hash = proof_payload.pop("proof_hash")
        _validate_hash(
            proof_hash,
            _sha256(proof_payload),
            "model_charge_accounting_hash",
            model.external_id,
        )
        if (
            proof["selection_id"] != model.site_id
            or proof["complete"] is not True
            or not isinstance(proof["net_charge"], int)
            or isinstance(proof["net_charge"], bool)
            or not isinstance(proof["contributions"], list)
        ):
            raise ValidationError("invalid_model_charge_accounting", model.external_id)
        if model.electronic_state is not None and proof["net_charge"] != model.electronic_state.net_charge:
            raise ValidationError("model_charge_accounting_state_mismatch", model.external_id)
        if set(proof) == constrained_proof_fields:
            groups = proof["constraint_groups"]
            if (
                any(
                    isinstance(proof[name], bool) or not isinstance(proof[name], int)
                    for name in ("core_target_charge", "cap_target_charge", "model_target_charge")
                )
                or proof["core_target_charge"] + proof["cap_target_charge"]
                != proof["model_target_charge"]
                or proof["model_target_charge"] != proof["net_charge"]
                or not isinstance(groups, list)
                or not groups
            ):
                raise ValidationError("invalid_constrained_charge_accounting", model.external_id)
            constrained_atom_ids: list[str] = []
            constrained_targets = {"core": 0, "cap": 0}
            for index, group in enumerate(groups):
                if not isinstance(group, Mapping):
                    raise ValidationError("invalid_charge_constraint_group", str(index), model.external_id)
                required_group_fields = {
                    "group_id", "role", "model_atom_ids", "target_charge", "source", "complete",
                }
                atom_ids = group.get("model_atom_ids")
                if (
                    set(group) != required_group_fields
                    or not isinstance(group.get("group_id"), str)
                    or not group["group_id"]
                    or group.get("role") not in {"core", "cap"}
                    or not isinstance(atom_ids, list)
                    or not atom_ids
                    or any(not isinstance(atom_id, str) or not atom_id for atom_id in atom_ids)
                    or len(atom_ids) != len(set(atom_ids))
                    or isinstance(group.get("target_charge"), bool)
                    or not isinstance(group.get("target_charge"), int)
                    or not isinstance(group.get("source"), str)
                    or not group["source"]
                    or group.get("complete") is not True
                ):
                    raise ValidationError("invalid_charge_constraint_group", str(index), model.external_id)
                constrained_atom_ids.extend(atom_ids)
                constrained_targets[group["role"]] += group["target_charge"]
            if len(constrained_atom_ids) != len(set(constrained_atom_ids)):
                raise ValidationError("overlapping_charge_constraint_groups", model.external_id)
            if (
                constrained_targets["core"] != proof["core_target_charge"]
                or constrained_targets["cap"] != proof["cap_target_charge"]
            ):
                raise ValidationError("charge_constraint_target_mismatch", model.external_id)
        manifest = model.capped_model_manifest
        if manifest is not None:
            required_manifest_fields = {
                "schema_version", "model_id", "core_model_atom_ids",
                "environment_model_atom_ids", "cap_model_atom_ids",
                "cut_edge_ids", "complete", "manifest_hash",
            }
            if set(manifest) != required_manifest_fields:
                raise ValidationError("invalid_capped_model_manifest_fields", model.external_id)
            manifest_payload = dict(manifest)
            manifest_hash = manifest_payload.pop("manifest_hash")
            _validate_hash(
                manifest_hash,
                _sha256(manifest_payload),
                "capped_model_manifest_hash",
                model.external_id,
            )
            manifest_atom_fields = (
                "core_model_atom_ids",
                "environment_model_atom_ids",
                "cap_model_atom_ids",
            )
            if any(
                not isinstance(manifest[name], list)
                or any(not isinstance(atom_id, str) or not atom_id for atom_id in manifest[name])
                or len(manifest[name]) != len(set(manifest[name]))
                for name in manifest_atom_fields
            ):
                raise ValidationError("invalid_capped_model_manifest", model.external_id)
            manifest_atom_ids = [
                atom_id
                for name in manifest_atom_fields
                for atom_id in manifest[name]
            ]
            if (
                manifest["model_id"] != model.external_id
                or manifest["complete"] is not True
                or len(manifest_atom_ids) != len(set(manifest_atom_ids))
                or set(manifest_atom_ids) != {atom.model_atom_id for atom in model.atoms}
                or manifest["core_model_atom_ids"] != [
                    atom.model_atom_id for atom in model.atoms if atom.role == "core"
                ]
                or manifest["environment_model_atom_ids"]
                or manifest["cap_model_atom_ids"] != [
                    atom.model_atom_id for atom in model.atoms if atom.role == "cap"
                ]
                or manifest["cut_edge_ids"] != [
                    cut.external_id for cut in model.cut_edges
                ]
            ):
                raise ValidationError("invalid_capped_model_manifest", model.external_id)
            if set(proof) == constrained_proof_fields:
                group_atom_ids_by_role = {
                    role: {
                        atom_id
                        for group in proof["constraint_groups"]
                        if group["role"] == role
                        for atom_id in group["model_atom_ids"]
                    }
                    for role in ("core", "cap")
                }
                if (
                    group_atom_ids_by_role["core"] != set(manifest["core_model_atom_ids"])
                    or group_atom_ids_by_role["cap"] != set(manifest["cap_model_atom_ids"])
                ):
                    raise ValidationError("charge_constraint_atom_coverage_mismatch", model.external_id)
        _validate_hash(model.model_hash, model.computed_hash(), "model_hash", model.external_id)
        _validate_mol2(
            model.mol2_text,
            len(model.atoms),
            len(model.bonds) + len(model.links) + sum(atom.role == "cap" for atom in model.atoms),
            model.atomic_charge_role,
            model.external_id,
        )
        atom_ids = [atom.model_atom_id for atom in model.atoms]
        if len(atom_ids) != len(set(atom_ids)) or [atom.mol2_serial for atom in model.atoms] != list(range(1, len(model.atoms) + 1)):
            raise ValidationError("invalid_model_atom_mapping", model.external_id)
        for atom in model.atoms:
            if atom.role == "cap":
                if (
                    atom.external_id is not None
                    or atom.canonical_atom_id != 0
                    or not atom.cap_parent_external_id
                    or not atom.cut_edge_id
                    or not atom.geometry_provenance
                    or not atom.charge_projection_group
                ):
                    raise ValidationError("invalid_cap_mapping", atom.model_atom_id)
            else:
                if atom.role != "core":
                    raise ValidationError("invalid_model_atom_role", atom.model_atom_id)
                parent = topology_atom_by_id.get(str(atom.external_id))
                if parent is None or atom.canonical_atom_id != parent.canonical_atom_id or atom.element != parent.element:
                    raise ValidationError("model_parent_identity_mismatch", atom.model_atom_id)
            expected_mapping[(model.external_id, atom.model_atom_id)] = (
                model.site_id, atom.mol2_serial, atom.external_id, atom.canonical_atom_id, atom.element,
                atom.role, atom.cap_parent_external_id, atom.cut_edge_id,
            )
        _validate_model_topology(model, site_by_id[model.site_id], topology)
    if len(bundle.mapping.records) != len(expected_mapping):
        raise ValidationError("incomplete_model_mapping", "mapping must cover every model atom")
    for record in bundle.mapping.records:
        expected = expected_mapping.get((record.model_id, record.model_atom_id))
        actual = (
            record.site_id, record.mol2_serial, record.external_id, record.canonical_atom_id,
            record.element, record.role, record.cap_parent_external_id, record.cut_edge_id,
        )
        if expected != actual:
            raise ValidationError("model_mapping_identity_mismatch", record.model_atom_id)


def validate_closed_mapping(
    request: MetalAssignmentInput,
    artifacts: StructuralArtifactBundle,
    models: DerivedModelBundle,
    mapping: ClosedAtomMapping,
) -> None:
    if (
        mapping.schema_version != request.schema_version
        or mapping.graph_revision != request.graph_revision
        or mapping.input_hash != request.input_hash
        or mapping.coordinate_unit != request.topology.coordinate_unit
        or mapping.prepared_topology_hash != request.topology.topology_hash
    ):
        raise ValidationError("closed_mapping_identity_mismatch", "mapping belongs to another topology")
    _validate_hash(mapping.mapping_hash, mapping.computed_hash(), "closed_mapping_hash", "atom_mapping")
    prepared_records = artifacts.prepared_system.mapping
    expected_structural = (
        *prepared_records,
        *(record for artifact in artifacts.assignment_components for record in artifact.mapping),
    )
    if mapping.parent_records != prepared_records:
        raise ValidationError("parent_mapping_mismatch", "parent mapping must equal prepared-system mapping")
    if mapping.structural_records != expected_structural:
        raise ValidationError("structural_mapping_mismatch", "structural mappings differ from artifacts")
    if mapping.model_records != models.mapping.records:
        raise ValidationError("model_mapping_mismatch", "closed mapping differs from derived model mapping")


def validate_artifact_package(request: MetalAssignmentInput, package: PreparedArtifactPackage) -> None:
    validate_structural_artifacts(request, package.structural_artifacts)
    validate_derived_models(request, package.derived_models)
    validate_closed_mapping(request, package.structural_artifacts, package.derived_models, package.atom_mapping)


def _parse_structural_mapping(value: Any, path: str) -> StructuralAtomMapping:
    data = _strict_object(
        value,
        required={
            "artifact_id", "artifact_atom_id", "mol2_serial", "external_id", "canonical_atom_id", "element",
            "chemical_component_id", "canonical_residue_id", "parameterization_component_id",
            "simulation_component_id", "simulation_residue_id",
        },
        path=path,
    )
    return StructuralAtomMapping(**data)


def _parse_structural_artifact(value: Any, path: str) -> StructuralArtifact:
    data = _strict_object(
        value,
        required={
            "external_id", "purpose", "classification", "provider", "base_force_field", "coordinate_unit",
            "net_formal_charge", "charge_resolution_method", "atom_formal_charges_complete",
            "atomic_charge_role", "atom_ids", "active_edge_ids", "mapping", "mol2_text", "artifact_hash",
        },
        path=path,
    )
    if not isinstance(data["mapping"], list):
        raise ValidationError("invalid_wire_type", "expected array", f"{path}.mapping")
    return StructuralArtifact(
        data["external_id"], data["purpose"], data["classification"], data["provider"],
        data["base_force_field"], data["net_formal_charge"], data["charge_resolution_method"],
        data["atom_formal_charges_complete"], data["coordinate_unit"], data["atomic_charge_role"],
        _strings(data["atom_ids"], f"{path}.atom_ids"),
        _strings(data["active_edge_ids"], f"{path}.active_edge_ids"),
        tuple(_parse_structural_mapping(item, f"{path}.mapping[{index}]") for index, item in enumerate(data["mapping"])),
        data["mol2_text"], data["artifact_hash"],
    )


def structural_artifact_bundle_from_dict(value: Any) -> StructuralArtifactBundle:
    data = _strict_object(
        value,
        required={
            "schema_version", "graph_revision", "input_hash", "coordinate_unit", "prepared_system",
            "assignment_components", "bundle_hash",
        },
        path="structural_artifact_bundle",
    )
    if not isinstance(data["assignment_components"], list):
        raise ValidationError("invalid_wire_type", "expected array", "structural_artifact_bundle.assignment_components")
    return StructuralArtifactBundle(
        data["schema_version"], data["graph_revision"], data["input_hash"], data["coordinate_unit"],
        _parse_structural_artifact(data["prepared_system"], "structural_artifact_bundle.prepared_system"),
        tuple(
            _parse_structural_artifact(item, f"structural_artifact_bundle.assignment_components[{index}]")
            for index, item in enumerate(data["assignment_components"])
        ),
        data["bundle_hash"],
    )


def _parse_model_atom(value: Any, path: str) -> ModelAtom:
    data = _strict_object(
        value,
        required={
            "model_atom_id", "external_id", "canonical_atom_id", "mol2_serial", "element", "coordinates",
            "role", "canonical_residue_id", "chemical_component_id", "partial_charge",
            "cap_parent_external_id", "cut_edge_id", "geometry_provenance", "charge_projection_group",
        },
        path=path,
    )
    coordinates = data["coordinates"]
    if not isinstance(coordinates, list) or len(coordinates) != 3:
        raise ValidationError("invalid_wire_type", "coordinates must have three values", f"{path}.coordinates")
    return ModelAtom(
        data["model_atom_id"], data["external_id"], data["canonical_atom_id"], data["mol2_serial"],
        data["element"], tuple(coordinates), data["role"], data["canonical_residue_id"],
        data["chemical_component_id"], data["partial_charge"], data["cap_parent_external_id"],
        data["cut_edge_id"], data["geometry_provenance"], data["charge_projection_group"],
    )


def _parse_model_bond(value: Any, path: str) -> ModelBond:
    data = _strict_object(value, required={"external_id", "model_atom_ids", "order", "semantic", "source"}, path=path)
    return ModelBond(data["external_id"], _strings(data["model_atom_ids"], f"{path}.model_atom_ids"), data["order"], data["semantic"], data["source"])


def _parse_model_link(value: Any, path: str) -> ModelLink:
    data = _strict_object(value, required={"external_id", "model_atom_ids", "kind", "source"}, path=path)
    return ModelLink(data["external_id"], _strings(data["model_atom_ids"], f"{path}.model_atom_ids"), data["kind"], data["source"])


def _parse_model_cut(value: Any, path: str) -> ModelCutEdge:
    data = _strict_object(value, required={"external_id", "retained_external_id", "excluded_external_id", "semantic", "source"}, path=path)
    return ModelCutEdge(**data)


def _parse_derived_model(value: Any, path: str) -> DerivedModel:
    data = _strict_object(
        value,
        required={
            "external_id", "site_id", "purpose", "coordinate_unit", "atomic_charge_role", "atoms", "bonds",
            "links", "cut_edges", "electronic_state", "charge_accounting", "mol2_text", "model_hash",
        },
        optional={"capped_model_manifest"},
        path=path,
    )
    for name in ("atoms", "bonds", "links", "cut_edges"):
        if not isinstance(data[name], list):
            raise ValidationError("invalid_wire_type", "expected array", f"{path}.{name}")
    raw_state = data["electronic_state"]
    state = None
    if raw_state is not None:
        state_data = _strict_object(
            raw_state,
            required={"selection_id", "net_charge", "spin_multiplicity", "source"},
            path=f"{path}.electronic_state",
        )
        state = ModelElectronicState(**state_data)
    if not isinstance(data["charge_accounting"], Mapping):
        raise ValidationError(
            "invalid_wire_type",
            "expected object",
            f"{path}.charge_accounting",
        )
    return DerivedModel(
        data["external_id"], data["site_id"], data["purpose"], data["coordinate_unit"],
        data["atomic_charge_role"], state,
        tuple(_parse_model_atom(item, f"{path}.atoms[{index}]") for index, item in enumerate(data["atoms"])),
        tuple(_parse_model_bond(item, f"{path}.bonds[{index}]") for index, item in enumerate(data["bonds"])),
        tuple(_parse_model_link(item, f"{path}.links[{index}]") for index, item in enumerate(data["links"])),
        tuple(_parse_model_cut(item, f"{path}.cut_edges[{index}]") for index, item in enumerate(data["cut_edges"])),
        data["charge_accounting"], data["mol2_text"], data["model_hash"],
        data.get("capped_model_manifest"),
    )


def _parse_model_mapping(value: Any, path: str) -> ModelMappingRecord:
    data = _strict_object(
        value,
        required={
            "site_id", "model_id", "model_atom_id", "mol2_serial", "external_id", "canonical_atom_id",
            "element", "role", "cap_parent_external_id", "cut_edge_id",
        },
        path=path,
    )
    return ModelMappingRecord(**data)


def _parse_site(value: Any, path: str) -> DerivedSite:
    data = _strict_object(
        value,
        required={
            "external_id", "metal_atom_ids", "canonical_residue_ids", "atom_ids", "interaction_edge_ids",
            "model_ids", "multi_metal",
        },
        path=path,
    )
    return DerivedSite(
        data["external_id"], _strings(data["metal_atom_ids"], f"{path}.metal_atom_ids"),
        _strings(data["canonical_residue_ids"], f"{path}.canonical_residue_ids"),
        _strings(data["atom_ids"], f"{path}.atom_ids"),
        _strings(data["interaction_edge_ids"], f"{path}.interaction_edge_ids"),
        _strings(data["model_ids"], f"{path}.model_ids"), data["multi_metal"],
    )


def derived_model_bundle_from_dict(value: Any) -> DerivedModelBundle:
    data = _strict_object(
        value,
        required={
            "schema_version", "graph_revision", "input_hash", "interaction_model", "coordinate_unit",
            "models_generated", "skipped_reason", "protocol", "source_atom_ids",
            "active_metal_interaction_edge_ids", "removed_interaction_edge_ids", "sites", "models", "mapping",
            "bundle_hash",
        },
        path="derived_model_bundle",
    )
    for name in ("sites", "models"):
        if not isinstance(data[name], list):
            raise ValidationError("invalid_wire_type", "expected array", f"derived_model_bundle.{name}")
    if not isinstance(data["protocol"], Mapping):
        raise ValidationError("invalid_wire_type", "expected object", "derived_model_bundle.protocol")
    mapping = _strict_object(data["mapping"], required={"records", "mapping_hash"}, path="derived_model_bundle.mapping")
    if not isinstance(mapping["records"], list):
        raise ValidationError("invalid_wire_type", "expected array", "derived_model_bundle.mapping.records")
    return DerivedModelBundle(
        data["schema_version"], data["graph_revision"], data["input_hash"], data["interaction_model"],
        data["coordinate_unit"], data["models_generated"], data["skipped_reason"], data["protocol"],
        _strings(data["source_atom_ids"], "derived_model_bundle.source_atom_ids"),
        _strings(data["active_metal_interaction_edge_ids"], "derived_model_bundle.active_metal_interaction_edge_ids"),
        _strings(data["removed_interaction_edge_ids"], "derived_model_bundle.removed_interaction_edge_ids"),
        tuple(_parse_site(item, f"derived_model_bundle.sites[{index}]") for index, item in enumerate(data["sites"])),
        tuple(_parse_derived_model(item, f"derived_model_bundle.models[{index}]") for index, item in enumerate(data["models"])),
        DerivedModelMapping(
            tuple(
                _parse_model_mapping(item, f"derived_model_bundle.mapping.records[{index}]")
                for index, item in enumerate(mapping["records"])
            ),
            mapping["mapping_hash"],
        ),
        data["bundle_hash"],
    )


def closed_atom_mapping_from_dict(value: Any) -> ClosedAtomMapping:
    data = _strict_object(
        value,
        required={
            "schema_version", "graph_revision", "input_hash", "coordinate_unit", "prepared_topology_hash",
            "parent_records", "structural_records", "model_records", "mapping_hash",
        },
        path="atom_mapping",
    )
    for name in ("parent_records", "structural_records", "model_records"):
        if not isinstance(data[name], list):
            raise ValidationError("invalid_wire_type", "expected array", f"atom_mapping.{name}")
    return ClosedAtomMapping(
        data["schema_version"], data["graph_revision"], data["input_hash"], data["coordinate_unit"],
        data["prepared_topology_hash"],
        tuple(_parse_structural_mapping(item, f"atom_mapping.parent_records[{index}]") for index, item in enumerate(data["parent_records"])),
        tuple(_parse_structural_mapping(item, f"atom_mapping.structural_records[{index}]") for index, item in enumerate(data["structural_records"])),
        tuple(_parse_model_mapping(item, f"atom_mapping.model_records[{index}]") for index, item in enumerate(data["model_records"])),
        data["mapping_hash"],
    )


def prepared_artifact_package_from_dict(value: Any) -> PreparedArtifactPackage:
    data = _strict_object(
        value,
        required={"structural_artifacts", "atom_mapping", "derived_models"},
        path="prepared_artifacts",
    )
    return PreparedArtifactPackage(
        structural_artifact_bundle_from_dict(data["structural_artifacts"]),
        closed_atom_mapping_from_dict(data["atom_mapping"]),
        derived_model_bundle_from_dict(data["derived_models"]),
    )


__all__ = [
    "ATOMIC_CHARGE_ROLES", "ClosedAtomMapping", "DerivedModel", "DerivedModelBundle", "DerivedModelMapping",
    "ModelElectronicState",
    "DerivedSite", "ModelAtom", "ModelBond", "ModelCutEdge", "ModelLink", "ModelMappingRecord",
    "PreparedArtifactPackage", "StructuralArtifact", "StructuralArtifactBundle", "StructuralAtomMapping",
    "closed_atom_mapping_from_dict", "derived_model_bundle_from_dict", "prepared_artifact_package_from_dict",
    "structural_artifact_bundle_from_dict", "validate_artifact_package", "validate_closed_mapping",
    "validate_derived_models", "validate_structural_artifacts",
]
