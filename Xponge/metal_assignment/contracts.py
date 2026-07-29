"""Versioned, topology-preserving contracts for metal force-field assignment."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, fields, field, replace
from enum import Enum
import hashlib
import json
import math
import re
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = 1
BASE_FORCE_FIELD_SCOPE = "baseForceField"
COORDINATE_UNIT = "angstrom"


class ValidationError(ValueError):
    """Structured contract validation failure."""

    def __init__(self, code: str, message: str, path: str = "") -> None:
        self.code = code
        self.path = path
        suffix = f" ({path})" if path else ""
        super().__init__(f"{code}: {message}{suffix}")


class SimulationSplitReason(str, Enum):
    NONE_CONNECTED = "none_connected"
    DISCONNECTED_AFTER_INTERACTION_FILTER = "disconnected_after_interaction_filter"


class ChargePolicy(str, Enum):
    PRESERVE_COMPONENT_CHARGES = "preserve_component_charges"
    SITE_JOINT_FIT = "site_joint_fit"
    MIXED_GROUP_CONSTRAINTS = "mixed_group_constraints"


@dataclass(frozen=True, slots=True)
class ProviderCapabilitySnapshot:
    schema_version: int
    provider_name: str
    provider_version: str
    provider_revision: str
    projection_schema_version: int
    embedded_metal: bool
    standalone_metal: bool
    multi_metal: bool
    bonded: bool
    nonbonded_12_6: bool
    base_force_field_providers: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ResiduePartitionProof:
    canonical_residue_id: str
    source_atom_ids: tuple[str, ...]
    active_internal_edge_ids: tuple[str, ...]
    removed_interaction_edge_ids: tuple[str, ...]
    simulation_partitions: Mapping[str, tuple[str, ...]]
    proof_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "simulation_partitions")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "canonical_residue_id": self.canonical_residue_id,
            "source_atom_ids": _canonicalize(self.source_atom_ids),
            "active_internal_edge_ids": _canonicalize(self.active_internal_edge_ids),
            "removed_interaction_edge_ids": _canonicalize(self.removed_interaction_edge_ids),
            "simulation_partitions": _canonicalize(self.simulation_partitions),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ResiduePartitionProof":
        return replace(self, proof_hash=self.computed_hash())


class FrozenMapping(Mapping[str, Any]):
    """Small recursively immutable mapping used by frozen wire contracts."""

    __slots__ = ("_items", "_lookup")

    def __init__(self, value: Mapping[str, Any] | None = None) -> None:
        items = tuple(sorted(((str(key), _freeze(item)) for key, item in (value or {}).items())))
        self._items = items
        self._lookup = dict(items)

    def __getitem__(self, key: str) -> Any:
        return self._lookup[key]

    def __iter__(self) -> Iterator[str]:
        return (key for key, _ in self._items)

    def __len__(self) -> int:
        return len(self._items)

    def __repr__(self) -> str:
        return f"FrozenMapping({dict(self._items)!r})"


def _freeze(value: Any) -> Any:
    if isinstance(value, FrozenMapping):
        return value
    if isinstance(value, Mapping):
        return FrozenMapping(value)
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted(_freeze(item) for item in value))
    return value


def _canonicalize(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if hasattr(value, "__dataclass_fields__"):
        return {item.name: _canonicalize(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, Mapping):
        return {
            str(key): _canonicalize(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_canonicalize(item) for item in value]
    return value


def _sha256(payload: Mapping[str, Any]) -> str:
    try:
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValidationError("noncanonical_payload", str(exc)) from exc
    return hashlib.sha256(encoded).hexdigest()


def _freeze_field(instance: Any, name: str) -> None:
    object.__setattr__(instance, name, _freeze(getattr(instance, name)))


@dataclass(frozen=True, slots=True)
class ElectronicState:
    net_charge: int
    spin_multiplicity: int


@dataclass(frozen=True, slots=True)
class PreparedAtom:
    external_id: str
    canonical_atom_id: int
    stable_order: int
    element: str
    coordinates: tuple[float, float, float]
    chemical_component_id: str
    canonical_residue_id: str
    simulation_component_id: str
    simulation_residue_id: str
    formal_charge: int | None = None
    formal_charge_known: bool = False
    partial_charge: float | None = None
    is_metal: bool = False
    scopes: tuple[str, ...] = ()
    atom_name: str = ""


@dataclass(frozen=True, slots=True)
class PreparedResidue:
    external_id: str
    canonical_residue_id: str
    atom_ids: tuple[str, ...]
    chemical_component_id: str
    simulation_component_id: str
    net_formal_charge: int
    charge_resolution_method: str
    atom_formal_charges_complete: bool
    residue_name: str = ""
    chain_id: str = ""
    res_seq: int = 0
    icode: str = ""
    polymer_position: str = "nonpolymer"
    template_id: str = ""


@dataclass(frozen=True, slots=True)
class PreparedComponent:
    external_id: str
    atom_ids: tuple[str, ...]
    chemical_component_id: str
    simulation_component_id: str
    net_formal_charge: int
    charge_resolution_method: str
    atom_formal_charges_complete: bool


@dataclass(frozen=True, slots=True)
class PreparedBond:
    external_id: str
    atom_ids: tuple[str, str]
    order: float = 1.0
    semantic: str = "covalent"
    active: bool = True
    source: str = "input"


@dataclass(frozen=True, slots=True)
class PreparedLink:
    external_id: str
    atom_ids: tuple[str, str]
    kind: str = "coordination"
    active: bool = True
    confirmed: bool = True
    source: str = "input"


@dataclass(frozen=True, slots=True)
class PreparedChemicalTopology:
    schema_version: int
    graph_revision: int
    input_hash: str
    coordinate_unit: str
    atoms: tuple[PreparedAtom, ...]
    residues: tuple[PreparedResidue, ...]
    components: tuple[PreparedComponent, ...]
    bonds: tuple[PreparedBond, ...] = ()
    links: tuple[PreparedLink, ...] = ()
    topology_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "graph_revision": self.graph_revision,
            "input_hash": self.input_hash,
            "coordinate_unit": self.coordinate_unit,
            "atoms": _canonicalize(self.atoms),
            "residues": _canonicalize(self.residues),
            "components": _canonicalize(self.components),
            "bonds": _canonicalize(self.bonds),
            "links": _canonicalize(self.links),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "PreparedChemicalTopology":
        return replace(self, topology_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class ProviderProjection:
    external_id: str
    chemical_component_id: str
    canonical_residue_id: str
    parameterization_component_id: str
    simulation_component_id: str
    simulation_residue_id: str
    atom_ids: tuple[str, ...]
    included_scopes_by_atom_or_component: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    assignment_only_split: bool = False
    simulation_split: bool = False
    simulation_split_reason: SimulationSplitReason = SimulationSplitReason.NONE_CONNECTED
    parent_component_id: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "included_scopes_by_atom_or_component")


@dataclass(frozen=True, slots=True)
class ChemicalTopologyProof:
    schema_version: int
    graph_revision: int
    perception_rule_version: str
    atom_ids: tuple[str, ...]
    explicit_hydrogen_count_by_atom: Mapping[str, int]
    implicit_hydrogen_candidates_by_atom: Mapping[str, tuple[int, ...]]
    covalent_bond_order_by_atom: Mapping[str, float]
    valence_status_by_atom: Mapping[str, str]
    unresolved_valence_atom_ids: tuple[str, ...]
    explicit_hydrogen_status: str
    hydrogen_resolution_method: str
    valence_complete: bool
    proof_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "explicit_hydrogen_count_by_atom")
        _freeze_field(self, "implicit_hydrogen_candidates_by_atom")
        _freeze_field(self, "covalent_bond_order_by_atom")
        _freeze_field(self, "valence_status_by_atom")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "proof_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ChemicalTopologyProof":
        return replace(self, proof_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class AssignmentComponent:
    external_id: str
    atom_ids: tuple[str, ...]
    classification: str
    provider: str
    net_formal_charge: int
    charge_resolution_method: str
    atom_formal_charges_complete: bool
    chemical_topology_proof: ChemicalTopologyProof
    base_force_field: bool = True


def atomic_center_atom_ids(request: Any) -> frozenset[str]:
    """Return the structural metal centers selected for the metal overlay."""

    return frozenset(
        atom_id
        for component in request.assignment_components
        if component.classification == "atomic_center"
        for atom_id in component.atom_ids
    )


def base_metal_component_supported(component: AssignmentComponent) -> bool:
    """Whether a metal-containing component belongs to ordinary ion assignment."""

    return component.base_force_field and component.classification == "solvent_ion"


@dataclass(frozen=True, slots=True)
class ChargeConstraintGroup:
    external_id: str
    atom_ids: tuple[str, ...]
    target_charge: float
    source: str = "workflow"


@dataclass(frozen=True, slots=True)
class ChargeAssignmentContract:
    policy: ChargePolicy
    fit_atom_ids: tuple[str, ...]
    component_formal_charges: Mapping[str, int]
    component_targets: Mapping[str, float] = field(default_factory=dict)
    site_target: float | None = None
    fixed_atom_ids: tuple[str, ...] = ()
    constraint_groups: tuple[ChargeConstraintGroup, ...] = ()
    charge_unit: str = "elementary_charge"
    source: str = "workflow"
    fit_precedence: int = 200

    def __post_init__(self) -> None:
        _freeze_field(self, "component_formal_charges")
        _freeze_field(self, "component_targets")


@dataclass(frozen=True, slots=True)
class PartialChargeArtifact:
    schema_version: int
    artifact_id: str
    topology_hash: str
    graph_revision: int
    input_hash: str
    atom_ids: tuple[str, ...]
    charges: tuple[float, ...]
    scope_kind: str
    scope_id: str
    atomic_charge_role: str
    precedence: int
    provider: str
    provider_version: str
    method: str
    charge_unit: str
    source: str
    artifact_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "artifact_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "PartialChargeArtifact":
        return replace(self, artifact_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class PartialChargeArtifactBundle:
    schema_version: int
    topology_hash: str
    graph_revision: int
    input_hash: str
    artifacts: tuple[PartialChargeArtifact, ...]
    bundle_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "bundle_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "PartialChargeArtifactBundle":
        return replace(self, bundle_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class BaseForceFieldOverlay:
    topology_hash: str
    covered_atom_ids: tuple[str, ...] = ()
    atom_types: Mapping[str, str] = field(default_factory=dict)
    charges: Mapping[str, float] = field(default_factory=dict)
    masses: Mapping[str, float] = field(default_factory=dict)
    lj_parameters: Mapping[str, Any] = field(default_factory=dict)
    bonded_parameters: Mapping[str, Any] = field(default_factory=dict)
    parameter_source: str = ""

    def __post_init__(self) -> None:
        for name in ("atom_types", "charges", "masses", "lj_parameters", "bonded_parameters"):
            _freeze_field(self, name)


@dataclass(frozen=True, slots=True)
class ChargeOverlay:
    topology_hash: str
    charges: Mapping[str, float]
    source: str
    atom_sources: Mapping[str, Any] = field(default_factory=dict)
    artifact_hashes: tuple[str, ...] = ()
    overlay_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "charges")
        _freeze_field(self, "atom_sources")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "overlay_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ChargeOverlay":
        return replace(self, overlay_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class BondedParameterOverlay:
    topology_hash: str
    terms: Mapping[str, Any]
    source: str

    def __post_init__(self) -> None:
        _freeze_field(self, "terms")


@dataclass(frozen=True, slots=True)
class MetalParameterOverlay:
    topology_hash: str
    covered_atom_ids: tuple[str, ...] = ()
    atom_types: Mapping[str, str] = field(default_factory=dict)
    charges: Mapping[str, float] = field(default_factory=dict)
    masses: Mapping[str, float] = field(default_factory=dict)
    lj_parameters: Mapping[str, Any] = field(default_factory=dict)
    bonded_parameters: Mapping[str, Any] = field(default_factory=dict)
    parameter_source: str = ""
    precedence: int = 100

    def __post_init__(self) -> None:
        for name in ("atom_types", "charges", "masses", "lj_parameters", "bonded_parameters"):
            _freeze_field(self, name)


@dataclass(frozen=True, slots=True)
class ParameterizationResult:
    schema_version: int
    request_id: str
    input_hash: str
    graph_revision: int
    topology_hash: str
    projection_hash: str
    base_overlay: BaseForceFieldOverlay
    metal_overlay: MetalParameterOverlay
    charge_overlay: ChargeOverlay | None = None
    bonded_overlay: BondedParameterOverlay | None = None
    fit_reports: Mapping[str, Any] = field(default_factory=dict)
    provenance: Mapping[str, Any] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    status: str = "overlay_validated"
    application_audit: Mapping[str, Any] = field(default_factory=dict)
    complete: bool = False
    result_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "fit_reports")
        _freeze_field(self, "provenance")
        _freeze_field(self, "application_audit")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "result_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ParameterizationResult":
        return replace(self, result_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class MetalAssignmentInput:
    schema_version: int
    request_id: str
    interaction_model: str
    electronic_state: ElectronicState
    graph_revision: int
    input_hash: str
    topology: PreparedChemicalTopology
    projections: tuple[ProviderProjection, ...]
    assignment_components: tuple[AssignmentComponent, ...]
    partition_proofs: tuple[ResiduePartitionProof, ...]
    capability_snapshot: ProviderCapabilitySnapshot
    charge_contract: ChargeAssignmentContract | None = None
    partial_charge_artifacts: PartialChargeArtifactBundle | None = None
    projection_hash: str = ""

    def canonical_projection_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "request_id": self.request_id,
            "interaction_model": self.interaction_model,
            "electronic_state": _canonicalize(self.electronic_state),
            "graph_revision": self.graph_revision,
            "input_hash": self.input_hash,
            "topology_hash": self.topology.topology_hash,
            "projections": _canonicalize(self.projections),
            "assignment_components": _canonicalize(self.assignment_components),
            "partition_proofs": _canonicalize(self.partition_proofs),
            "charge_contract": _canonicalize(self.charge_contract),
            "partial_charge_artifacts": _canonicalize(self.partial_charge_artifacts),
            "capability_snapshot": _canonicalize(self.capability_snapshot),
        }

    def computed_projection_hash(self) -> str:
        return _sha256(self.canonical_projection_payload())

    def with_computed_hash(self) -> "MetalAssignmentInput":
        return replace(self, projection_hash=self.computed_projection_hash())


def _unique_ids(items: Sequence[Any], *, path: str) -> set[str]:
    ids = [str(item.external_id) for item in items]
    if any(not item_id for item_id in ids):
        raise ValidationError("invalid_external_id", "external IDs must be non-empty strings", path)
    if len(ids) != len(set(ids)):
        raise ValidationError("duplicate_reference", "external IDs must be unique", path)
    return set(ids)


def _validate_partition(groups: Sequence[Any], atom_ids: set[str], *, path: str) -> None:
    members: list[str] = []
    for group in groups:
        if not group.atom_ids:
            raise ValidationError("empty_partition", "partition member has no atoms", f"{path}.{group.external_id}")
        if len(group.atom_ids) != len(set(group.atom_ids)):
            raise ValidationError("duplicate_atom_reference", group.external_id, path)
        if not set(group.atom_ids) <= atom_ids:
            raise ValidationError("missing_atom_reference", group.external_id, path)
        members.extend(group.atom_ids)
    if len(members) != len(set(members)):
        raise ValidationError("overlapping_atom_partition", "partition atom sets overlap", path)
    if set(members) != atom_ids:
        raise ValidationError("incomplete_atom_partition", "partition must cover every atom exactly once", path)


def _validate_residue_connectivity(topology: PreparedChemicalTopology) -> None:
    adjacency = {atom.external_id: set() for atom in topology.atoms}
    for edge in (*topology.bonds, *topology.links):
        if not edge.active:
            continue
        atom1, atom2 = edge.atom_ids
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)
    for residue in topology.residues:
        expected = set(residue.atom_ids)
        visited = {next(iter(expected))}
        pending = list(visited)
        while pending:
            current = pending.pop()
            for neighbor in adjacency[current] & expected:
                if neighbor not in visited:
                    visited.add(neighbor)
                    pending.append(neighbor)
        if visited != expected:
            raise ValidationError(
                "disconnected_simulation_residue",
                "each prepared simulation residue must be internally connected",
                f"residues.{residue.external_id}",
            )


def _valid_sha256(value: str) -> bool:
    return len(value) == 64 and all(character in "0123456789abcdef" for character in value)


def validate_topology(topology: PreparedChemicalTopology) -> None:
    if topology.schema_version != SCHEMA_VERSION:
        raise ValidationError("unsupported_schema_version", str(topology.schema_version), "topology.schema_version")
    if topology.graph_revision < 0:
        raise ValidationError("invalid_graph_revision", str(topology.graph_revision), "topology.graph_revision")
    if not _valid_sha256(topology.input_hash):
        raise ValidationError("invalid_input_hash", topology.input_hash, "topology.input_hash")
    if topology.coordinate_unit != COORDINATE_UNIT:
        raise ValidationError("unsupported_coordinate_unit", topology.coordinate_unit, "topology.coordinate_unit")

    atom_ids = _unique_ids(topology.atoms, path="atoms")
    residue_ids = _unique_ids(topology.residues, path="residues")
    _unique_ids(topology.components, path="components")
    _unique_ids(topology.bonds, path="bonds")
    _unique_ids(topology.links, path="links")
    _validate_partition(topology.residues, atom_ids, path="residues")
    _validate_partition(topology.components, atom_ids, path="components")

    canonical_atom_ids = [atom.canonical_atom_id for atom in topology.atoms]
    if any(not isinstance(value, int) or isinstance(value, bool) or value <= 0 for value in canonical_atom_ids):
        raise ValidationError("invalid_canonical_atom_id", "canonical atom IDs must be positive integers", "atoms")
    if len(canonical_atom_ids) != len(set(canonical_atom_ids)):
        raise ValidationError("duplicate_canonical_atom_id", "canonical atom IDs must be unique", "atoms")

    stable_orders = [atom.stable_order for atom in topology.atoms]
    if any(not isinstance(value, int) or isinstance(value, bool) or value < 0 for value in stable_orders):
        raise ValidationError("invalid_stable_order", "stable atom order must be a non-negative integer", "atoms")
    if len(stable_orders) != len(set(stable_orders)):
        raise ValidationError("duplicate_stable_order", "stable atom order must be unique", "atoms")

    residue_by_id = {item.external_id: item for item in topology.residues}
    component_by_simulation_id = {item.simulation_component_id: item for item in topology.components}
    if len(component_by_simulation_id) != len(topology.components):
        raise ValidationError("duplicate_simulation_component_id", "simulation component IDs must be unique", "components")
    atom_by_id = {item.external_id: item for item in topology.atoms}

    for atom in topology.atoms:
        if (
            not isinstance(atom.external_id, str)
            or not atom.external_id
            or not isinstance(atom.element, str)
            or not atom.element
            or len(atom.coordinates) != 3
            or not all(
                isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
                for value in atom.coordinates
            )
        ):
            raise ValidationError("invalid_coordinates", atom.external_id, f"atoms.{atom.external_id}.coordinates")
        if not isinstance(atom.formal_charge_known, bool):
            raise ValidationError(
                "invalid_formal_charge_known", atom.external_id, f"atoms.{atom.external_id}.formal_charge_known"
            )
        if atom.formal_charge_known != (atom.formal_charge is not None):
            raise ValidationError(
                "formal_charge_presence_mismatch",
                atom.external_id,
                f"atoms.{atom.external_id}.formal_charge_known",
            )
        if atom.formal_charge is not None and (
            not isinstance(atom.formal_charge, int) or isinstance(atom.formal_charge, bool)
        ):
            raise ValidationError("invalid_formal_charge", atom.external_id, f"atoms.{atom.external_id}.formal_charge")
        if atom.is_metal and atom.formal_charge is None:
            raise ValidationError("missing_metal_formal_charge", atom.external_id, f"atoms.{atom.external_id}.formal_charge")
        if atom.partial_charge is not None and (
            not isinstance(atom.partial_charge, (int, float))
            or isinstance(atom.partial_charge, bool)
            or not math.isfinite(atom.partial_charge)
        ):
            raise ValidationError("invalid_partial_charge", atom.external_id, f"atoms.{atom.external_id}.partial_charge")
        if not isinstance(atom.is_metal, bool) or not all(isinstance(scope, str) and scope for scope in atom.scopes):
            raise ValidationError("invalid_atom_metadata", atom.external_id, f"atoms.{atom.external_id}")
        if not isinstance(atom.atom_name, str) or not atom.atom_name.strip():
            raise ValidationError("missing_atom_name", atom.external_id, f"atoms.{atom.external_id}.atom_name")
        if atom.simulation_residue_id not in residue_ids:
            raise ValidationError("missing_residue_reference", atom.simulation_residue_id, f"atoms.{atom.external_id}")
        residue = residue_by_id[atom.simulation_residue_id]
        if atom.external_id not in residue.atom_ids:
            raise ValidationError("residue_membership_mismatch", atom.external_id, f"atoms.{atom.external_id}")
        if atom.canonical_residue_id != residue.canonical_residue_id:
            raise ValidationError("canonical_residue_mismatch", atom.external_id, f"atoms.{atom.external_id}")
        if atom.chemical_component_id != residue.chemical_component_id:
            raise ValidationError("chemical_component_mismatch", atom.external_id, f"atoms.{atom.external_id}")
        if atom.simulation_component_id != residue.simulation_component_id:
            raise ValidationError("simulation_component_mismatch", atom.external_id, f"atoms.{atom.external_id}")
        component = component_by_simulation_id.get(atom.simulation_component_id)
        if component is None or atom.external_id not in component.atom_ids:
            raise ValidationError("component_membership_mismatch", atom.external_id, f"atoms.{atom.external_id}")
        if component.chemical_component_id != atom.chemical_component_id:
            raise ValidationError("chemical_component_mismatch", atom.external_id, f"components.{component.external_id}")
    for group_name, groups in (("residues", topology.residues), ("components", topology.components)):
        for group in groups:
            if group_name == "residues":
                if not isinstance(group.residue_name, str) or not group.residue_name.strip():
                    raise ValidationError(
                        "missing_residue_name", group.external_id, f"residues.{group.external_id}.residue_name"
                    )
                if not isinstance(group.chain_id, str) or not isinstance(group.icode, str):
                    raise ValidationError(
                        "invalid_residue_identity", group.external_id, f"residues.{group.external_id}"
                    )
                if not isinstance(group.res_seq, int) or isinstance(group.res_seq, bool):
                    raise ValidationError(
                        "invalid_residue_sequence", group.external_id, f"residues.{group.external_id}.res_seq"
                    )
                if group.polymer_position not in {
                    "nonpolymer", "unknown", "internal", "start", "end", "single",
                }:
                    raise ValidationError(
                        "invalid_polymer_position",
                        group.polymer_position,
                        f"residues.{group.external_id}.polymer_position",
                    )
            if not isinstance(group.net_formal_charge, int) or isinstance(group.net_formal_charge, bool):
                raise ValidationError("unresolved_component_charge", group.external_id, group_name)
            if not group.charge_resolution_method or not isinstance(group.atom_formal_charges_complete, bool):
                raise ValidationError("invalid_charge_resolution", group.external_id, group_name)
            members = [atom_by_id[atom_id] for atom_id in group.atom_ids]
            complete = all(atom.formal_charge is not None for atom in members)
            if complete != group.atom_formal_charges_complete:
                raise ValidationError("atom_charge_completeness_mismatch", group.external_id, group_name)
            if complete and sum(atom.formal_charge for atom in members) != group.net_formal_charge:
                raise ValidationError("component_charge_mismatch", group.external_id, group_name)

    for bond in topology.bonds:
        if len(bond.atom_ids) != 2 or bond.atom_ids[0] == bond.atom_ids[1] or not set(bond.atom_ids) <= atom_ids:
            raise ValidationError("invalid_edge_endpoint", bond.external_id, "bonds")
        if (
            not isinstance(bond.order, (int, float))
            or isinstance(bond.order, bool)
            or not math.isfinite(bond.order)
            or bond.order <= 0
        ):
            raise ValidationError("invalid_bond_order", bond.external_id, "bonds")
        if not isinstance(bond.active, bool) or not bond.source or not bond.semantic:
            raise ValidationError("missing_edge_source", bond.external_id, "bonds")
        if bond.active and atom_by_id[bond.atom_ids[0]].simulation_residue_id != atom_by_id[bond.atom_ids[1]].simulation_residue_id:
            raise ValidationError("cross_residue_bond_misclassified", bond.external_id, "bonds")
    for link in topology.links:
        if len(link.atom_ids) != 2 or link.atom_ids[0] == link.atom_ids[1] or not set(link.atom_ids) <= atom_ids:
            raise ValidationError("invalid_edge_endpoint", link.external_id, "links")
        if not isinstance(link.active, bool) or not isinstance(link.confirmed, bool):
            raise ValidationError("invalid_link_metadata", link.external_id, "links")
        if link.active and not link.confirmed:
            raise ValidationError("unconfirmed_active_link", link.external_id, "links")
        if not link.source:
            raise ValidationError("missing_edge_source", link.external_id, "links")
        if link.active and atom_by_id[link.atom_ids[0]].simulation_residue_id == atom_by_id[link.atom_ids[1]].simulation_residue_id:
            raise ValidationError("internal_link_misclassified", link.external_id, "links")

    if not topology.topology_hash:
        raise ValidationError("missing_topology_hash", "a computed topology hash is required", "topology_hash")
    if topology.topology_hash != topology.computed_hash():
        raise ValidationError("stale_topology_hash", "topology hash does not match the canonical payload")
    _validate_residue_connectivity(topology)


_ATOMIC_NUMBERS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36, "I": 53,
}


def validate_electronic_state(state: ElectronicState, elements: Sequence[str] | None = None) -> None:
    if (
        not isinstance(state.net_charge, int)
        or isinstance(state.net_charge, bool)
        or not isinstance(state.spin_multiplicity, int)
        or isinstance(state.spin_multiplicity, bool)
    ):
        raise ValidationError("invalid_electronic_state", "net charge and spin multiplicity must be integers")
    if state.spin_multiplicity < 1:
        raise ValidationError("invalid_electronic_state", "spin multiplicity must be positive")
    if elements is None:
        return
    unsupported = sorted(set(elements) - set(_ATOMIC_NUMBERS))
    if unsupported:
        raise ValidationError("unsupported_element", ",".join(unsupported), "topology.atoms")
    electrons = sum(_ATOMIC_NUMBERS[element] for element in elements) - state.net_charge
    unpaired = state.spin_multiplicity - 1
    if electrons < 0 or unpaired > electrons or (electrons - unpaired) % 2:
        raise ValidationError("spin_charge_parity", "spin multiplicity is incompatible with electron parity")


def bookkeeping_electronic_state(elements: Sequence[str], net_charge: int) -> ElectronicState:
    """Return a parity-compatible full-system state for charge bookkeeping."""

    unsupported = sorted(set(elements) - set(_ATOMIC_NUMBERS))
    if unsupported:
        raise ValidationError("unsupported_element", ",".join(unsupported), "topology.atoms")
    electrons = sum(_ATOMIC_NUMBERS[element] for element in elements) - net_charge
    state = ElectronicState(net_charge, 1 if electrons % 2 == 0 else 2)
    validate_electronic_state(state, elements)
    return state


def validate_charge_contract(contract: ChargeAssignmentContract, topology: PreparedChemicalTopology) -> None:
    if not isinstance(contract.policy, ChargePolicy):
        raise ValidationError("invalid_charge_policy", str(contract.policy), "charge_contract.policy")
    if (
        contract.charge_unit != "elementary_charge"
        or not contract.source
        or isinstance(contract.fit_precedence, bool)
        or not isinstance(contract.fit_precedence, int)
        or contract.fit_precedence <= 0
    ):
        raise ValidationError("invalid_charge_contract_metadata", "charge unit and source must be explicit")
    atom_ids = {atom.external_id for atom in topology.atoms}
    component_by_atom = {
        atom_id: component.external_id for component in topology.components for atom_id in component.atom_ids
    }
    component_ids = {component.external_id for component in topology.components}
    fit_ids = set(contract.fit_atom_ids)
    fixed_ids = set(contract.fixed_atom_ids)
    if len(contract.fit_atom_ids) != len(fit_ids) or len(contract.fixed_atom_ids) != len(fixed_ids):
        raise ValidationError("duplicate_charge_target", "charge scopes cannot contain duplicates")
    if not fit_ids or not fit_ids <= atom_ids or not fixed_ids <= atom_ids:
        raise ValidationError("missing_charge_target", "charge scopes must reference known atoms and fit scope cannot be empty")
    if fit_ids & fixed_ids:
        raise ValidationError("overlapping_charge_scopes", "fit and fixed charge scopes overlap")
    referenced_components = set(contract.component_formal_charges) | set(contract.component_targets)
    if not referenced_components <= component_ids:
        raise ValidationError("missing_component_charge_target", "charge contract references an unknown component")
    fit_components = {component_by_atom[atom_id] for atom_id in fit_ids}
    if not fit_components <= set(contract.component_formal_charges):
        raise ValidationError("missing_component_formal_charge", "every fitted component needs an explicit formal charge")

    group_ids = _unique_ids(contract.constraint_groups, path="charge_contract.constraint_groups") if contract.constraint_groups else set()
    del group_ids
    group_atoms: set[str] = set()
    group_targets = 0.0
    for group in contract.constraint_groups:
        atom_set = set(group.atom_ids)
        if (
            not group.atom_ids
            or len(atom_set) != len(group.atom_ids)
            or not atom_set <= fit_ids
            or group_atoms & atom_set
            or not group.source
            or not isinstance(group.target_charge, (int, float))
            or isinstance(group.target_charge, bool)
            or not math.isfinite(group.target_charge)
        ):
            raise ValidationError("invalid_charge_constraint_group", group.external_id)
        group_atoms.update(atom_set)
        group_targets += group.target_charge

    if contract.policy is ChargePolicy.PRESERVE_COMPONENT_CHARGES:
        if set(contract.component_targets) != fit_components:
            raise ValidationError("missing_component_charge_target", "each fitted component needs exactly one target")
        for component_id in fit_components:
            if abs(contract.component_targets[component_id] - contract.component_formal_charges[component_id]) > 1e-8:
                raise ValidationError("component_charge_mismatch", component_id)
        if contract.site_target is not None or contract.constraint_groups:
            raise ValidationError("unexpected_charge_constraint", "preserve-component policy only accepts component targets")
    elif contract.policy is ChargePolicy.SITE_JOINT_FIT:
        if contract.site_target is None:
            raise ValidationError("missing_site_charge_target", "site target is required by this charge policy")
        expected = sum(contract.component_formal_charges[component_id] for component_id in fit_components)
        if abs(contract.site_target - expected) > 1e-8:
            raise ValidationError("site_charge_mismatch", "site target must match fitted component formal charges")
        if contract.constraint_groups:
            raise ValidationError("unexpected_charge_constraint", "site-joint policy does not accept mixed groups")
    elif contract.policy is ChargePolicy.MIXED_GROUP_CONSTRAINTS:
        if contract.site_target is None or not contract.constraint_groups or group_atoms != fit_ids:
            raise ValidationError("missing_group_charge_target", "mixed policy groups must cover every fitted atom")
        if abs(group_targets - contract.site_target) > 1e-8:
            raise ValidationError("group_charge_mismatch", "group targets must sum to the site target")


def validate_partial_charge_artifacts(
    request: MetalAssignmentInput,
    bundle: PartialChargeArtifactBundle,
) -> None:
    if (
        bundle.schema_version != SCHEMA_VERSION
        or bundle.topology_hash != request.topology.topology_hash
        or bundle.graph_revision != request.graph_revision
        or bundle.input_hash != request.input_hash
    ):
        raise ValidationError("partial_charge_bundle_identity_mismatch", request.request_id)
    if not bundle.bundle_hash or bundle.bundle_hash != bundle.computed_hash():
        raise ValidationError("stale_partial_charge_bundle_hash", request.request_id)
    artifact_ids = [artifact.artifact_id for artifact in bundle.artifacts]
    if any(not artifact_id for artifact_id in artifact_ids) or len(artifact_ids) != len(set(artifact_ids)):
        raise ValidationError("duplicate_partial_charge_artifact", "artifact IDs must be non-empty and unique")

    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    component_by_id = {component.external_id: component for component in request.topology.components}
    group_by_id = {
        group.external_id: group
        for group in (request.charge_contract.constraint_groups if request.charge_contract is not None else ())
    }
    precedence_by_atom: dict[str, set[int]] = {}
    tolerance = 1.0e-6
    for artifact in bundle.artifacts:
        path = f"partial_charge_artifacts.{artifact.artifact_id or '<missing>'}"
        if (
            artifact.schema_version != SCHEMA_VERSION
            or artifact.topology_hash != bundle.topology_hash
            or artifact.graph_revision != bundle.graph_revision
            or artifact.input_hash != bundle.input_hash
        ):
            raise ValidationError("partial_charge_artifact_identity_mismatch", artifact.artifact_id, path)
        if not artifact.artifact_hash or artifact.artifact_hash != artifact.computed_hash():
            raise ValidationError("stale_partial_charge_artifact_hash", artifact.artifact_id, path)
        if (
            artifact.scope_kind not in {"component", "fit_scope", "constraint_group", "system"}
            or not artifact.scope_id
            or artifact.atomic_charge_role not in {"fixed", "reference", "fitted", "empirical"}
            or isinstance(artifact.precedence, bool)
            or not isinstance(artifact.precedence, int)
            or artifact.precedence <= 0
            or not artifact.provider
            or not artifact.provider_version
            or not artifact.method
            or artifact.charge_unit != "elementary_charge"
            or not artifact.source
        ):
            raise ValidationError("invalid_partial_charge_artifact_metadata", artifact.artifact_id, path)
        if (
            not artifact.atom_ids
            or len(artifact.atom_ids) != len(set(artifact.atom_ids))
            or len(artifact.atom_ids) != len(artifact.charges)
            or not set(artifact.atom_ids) <= set(atom_by_id)
        ):
            raise ValidationError("invalid_partial_charge_artifact_scope", artifact.artifact_id, path)
        expected_order = tuple(sorted(artifact.atom_ids, key=lambda atom_id: (atom_by_id[atom_id].stable_order, atom_id)))
        if artifact.atom_ids != expected_order:
            raise ValidationError("partial_charge_atom_order_mismatch", artifact.artifact_id, path)
        if any(
            isinstance(charge, bool)
            or not isinstance(charge, (int, float))
            or not math.isfinite(charge)
            for charge in artifact.charges
        ):
            raise ValidationError("invalid_partial_charge_value", artifact.artifact_id, path)

        atom_set = set(artifact.atom_ids)
        total_charge = float(sum(artifact.charges))
        if artifact.scope_kind == "component":
            component = component_by_id.get(artifact.scope_id)
            if component is None or atom_set != set(component.atom_ids):
                raise ValidationError("partial_charge_component_scope_mismatch", artifact.artifact_id, path)
            if abs(total_charge - component.net_formal_charge) > tolerance:
                raise ValidationError("partial_charge_component_conservation", artifact.artifact_id, path)
        elif artifact.scope_kind == "fit_scope":
            contract = request.charge_contract
            if contract is None or atom_set != set(contract.fit_atom_ids):
                raise ValidationError("partial_charge_fit_scope_mismatch", artifact.artifact_id, path)
            scoped = dict(zip(artifact.atom_ids, artifact.charges))
            if contract.policy is ChargePolicy.SITE_JOINT_FIT:
                if abs(total_charge - float(contract.site_target)) > tolerance:
                    raise ValidationError("partial_charge_fit_scope_conservation", artifact.artifact_id, path)
            elif contract.policy is ChargePolicy.PRESERVE_COMPONENT_CHARGES:
                for component_id, target in contract.component_targets.items():
                    component = component_by_id[component_id]
                    if abs(sum(scoped[atom_id] for atom_id in component.atom_ids) - target) > tolerance:
                        raise ValidationError("partial_charge_fit_scope_conservation", component_id, path)
            else:
                for group in contract.constraint_groups:
                    if abs(sum(scoped[atom_id] for atom_id in group.atom_ids) - group.target_charge) > tolerance:
                        raise ValidationError("partial_charge_fit_scope_conservation", group.external_id, path)
        elif artifact.scope_kind == "constraint_group":
            group = group_by_id.get(artifact.scope_id)
            if group is None or atom_set != set(group.atom_ids) or abs(total_charge - group.target_charge) > tolerance:
                raise ValidationError("partial_charge_constraint_group_mismatch", artifact.artifact_id, path)
        elif atom_set != set(atom_by_id) or abs(total_charge - request.electronic_state.net_charge) > tolerance:
            raise ValidationError("partial_charge_system_scope_mismatch", artifact.artifact_id, path)

        for atom_id in artifact.atom_ids:
            priorities = precedence_by_atom.setdefault(atom_id, set())
            if artifact.precedence in priorities:
                raise ValidationError("ambiguous_partial_charge_precedence", atom_id, path)
            priorities.add(artifact.precedence)


def validate_capability_snapshot(snapshot: ProviderCapabilitySnapshot, interaction_model: str) -> None:
    if snapshot.schema_version != SCHEMA_VERSION or snapshot.projection_schema_version != SCHEMA_VERSION:
        raise ValidationError("unsupported_capability_schema", str(snapshot.schema_version), "capability_snapshot")
    if not snapshot.provider_name or not snapshot.provider_version or not snapshot.provider_revision:
        raise ValidationError("invalid_capability_identity", "provider identity and revision are required")
    flags = (
        snapshot.embedded_metal,
        snapshot.standalone_metal,
        snapshot.multi_metal,
        snapshot.bonded,
        snapshot.nonbonded_12_6,
    )
    if not all(isinstance(value, bool) for value in flags):
        raise ValidationError("invalid_capability_flag", "capability flags must be booleans")
    if not snapshot.base_force_field_providers or not all(
        isinstance(provider, str) and provider for provider in snapshot.base_force_field_providers
    ):
        raise ValidationError("missing_base_force_field_provider", "at least one provider is required")
    if interaction_model == "bonded" and not snapshot.bonded:
        raise ValidationError("unsupported_provider_capability", "provider does not support bonded assignment")
    if interaction_model == "nonbonded_12_6" and not snapshot.nonbonded_12_6:
        raise ValidationError("unsupported_provider_capability", "provider does not support nonbonded 12-6 assignment")


def _validate_partition_proofs(request: MetalAssignmentInput) -> None:
    proof_ids = [proof.canonical_residue_id for proof in request.partition_proofs]
    if not proof_ids or len(proof_ids) != len(set(proof_ids)):
        raise ValidationError("duplicate_or_missing_partition_proof", "one proof is required per canonical residue")
    canonical_atom_sets: dict[str, set[str]] = {}
    for atom in request.topology.atoms:
        canonical_atom_sets.setdefault(atom.canonical_residue_id, set()).add(atom.external_id)
    if set(proof_ids) != set(canonical_atom_sets):
        raise ValidationError("incomplete_partition_proof", "proofs must cover every canonical residue")

    active_internal_edges: dict[str, set[str]] = {key: set() for key in canonical_atom_sets}
    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    active_edge_ids = {edge.external_id for edge in (*request.topology.bonds, *request.topology.links) if edge.active}
    for edge in (*request.topology.bonds, *request.topology.links):
        if not edge.active:
            continue
        atom1, atom2 = (atom_by_id[atom_id] for atom_id in edge.atom_ids)
        if atom1.canonical_residue_id == atom2.canonical_residue_id:
            active_internal_edges[atom1.canonical_residue_id].add(edge.external_id)

    proof_by_id = {proof.canonical_residue_id: proof for proof in request.partition_proofs}
    for canonical_residue_id, source_atoms in canonical_atom_sets.items():
        proof = proof_by_id[canonical_residue_id]
        if set(proof.source_atom_ids) != source_atoms or len(proof.source_atom_ids) != len(source_atoms):
            raise ValidationError("partition_source_mismatch", canonical_residue_id)
        partitions = proof.simulation_partitions
        partition_atoms = [atom_id for atom_ids in partitions.values() for atom_id in atom_ids]
        if (
            not partitions
            or len(partition_atoms) != len(set(partition_atoms))
            or set(partition_atoms) != source_atoms
        ):
            raise ValidationError("invalid_simulation_partition", canonical_residue_id)
        expected_residue_ids = {
            atom_by_id[atom_id].simulation_residue_id for atom_id in source_atoms
        }
        if set(partitions) != expected_residue_ids:
            raise ValidationError("simulation_partition_identity_mismatch", canonical_residue_id)
        for residue_id, atom_ids in partitions.items():
            if any(atom_by_id[atom_id].simulation_residue_id != residue_id for atom_id in atom_ids):
                raise ValidationError("simulation_partition_identity_mismatch", canonical_residue_id)
        if set(proof.active_internal_edge_ids) != active_internal_edges[canonical_residue_id]:
            raise ValidationError("active_edge_proof_mismatch", canonical_residue_id)
        if set(proof.removed_interaction_edge_ids) & active_edge_ids:
            raise ValidationError("removed_edge_still_active", canonical_residue_id)
        if not proof.proof_hash or proof.proof_hash != proof.computed_hash():
            raise ValidationError("stale_partition_proof_hash", canonical_residue_id)

    projections_by_residue: dict[str, list[ProviderProjection]] = {}
    for projection in request.projections:
        projections_by_residue.setdefault(projection.canonical_residue_id, []).append(projection)
    for canonical_residue_id, projections in projections_by_residue.items():
        proof = proof_by_id[canonical_residue_id]
        assignment_split = len(projections) > 1
        simulation_split = len(proof.simulation_partitions) > 1
        for projection in projections:
            if projection.assignment_only_split != assignment_split or projection.simulation_split != simulation_split:
                raise ValidationError("projection_treatment_mismatch", projection.external_id)
            if projection.simulation_residue_id not in proof.simulation_partitions:
                raise ValidationError("projection_partition_mismatch", projection.external_id)
        source_atoms = canonical_atom_sets[canonical_residue_id]
        metal_atoms = {atom_id for atom_id in source_atoms if atom_by_id[atom_id].is_metal}
        nonmetal_atoms = source_atoms - metal_atoms
        if request.interaction_model == "nonbonded_12_6" and metal_atoms and nonmetal_atoms:
            for partition in proof.simulation_partitions.values():
                if set(partition) & metal_atoms and set(partition) & nonmetal_atoms:
                    raise ValidationError("nonbonded_embedded_metal_not_split", canonical_residue_id)
            for edge_id in proof.active_internal_edge_ids:
                edge = next(
                    edge for edge in (*request.topology.bonds, *request.topology.links) if edge.external_id == edge_id
                )
                endpoints = set(edge.atom_ids)
                if endpoints & metal_atoms and endpoints & nonmetal_atoms:
                    raise ValidationError("nonbonded_metal_ligand_edge_active", edge_id)


def _validate_projection_mapping(request: MetalAssignmentInput) -> None:
    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    component_ids = _unique_ids(request.assignment_components, path="assignment_components")
    _unique_ids(request.projections, path="projections")
    if not component_ids or not request.projections:
        raise ValidationError("missing_assignment_projection", "assignment components and projections are required")

    assignment_atoms: list[str] = []
    assignment_by_id = {}
    for component in request.assignment_components:
        if not component.atom_ids or len(component.atom_ids) != len(set(component.atom_ids)):
            raise ValidationError("empty_or_duplicate_assignment_component", component.external_id)
        if not set(component.atom_ids) <= set(atom_by_id):
            raise ValidationError("missing_atom_reference", component.external_id, "assignment_components")
        component_contains_metal = any(
            atom_by_id[atom_id].is_metal for atom_id in component.atom_ids
        )
        if component.classification == "atomic_center":
            if component.base_force_field:
                raise ValidationError(
                    "metal_in_base_force_field",
                    component.external_id,
                    "assignment_components",
                )
            if not component_contains_metal:
                raise ValidationError(
                    "invalid_atomic_center_assignment",
                    component.external_id,
                    "assignment_components",
                )
        elif (
            component.base_force_field
            and component_contains_metal
            and not base_metal_component_supported(component)
        ):
            raise ValidationError(
                "metal_in_base_force_field",
                component.external_id,
                "assignment_components",
            )
        if isinstance(component.net_formal_charge, bool) or not isinstance(component.net_formal_charge, int):
            raise ValidationError("invalid_component_net_charge", component.external_id, "assignment_components")
        if not component.charge_resolution_method:
            raise ValidationError("missing_charge_resolution_method", component.external_id, "assignment_components")
        expected_charge_completeness = all(atom_by_id[atom_id].formal_charge_known for atom_id in component.atom_ids)
        if component.classification == "atomic_center":
            if not expected_charge_completeness:
                raise ValidationError("missing_metal_formal_charge", component.external_id)
            expected_metal_charge = sum(
                int(atom_by_id[atom_id].formal_charge) for atom_id in component.atom_ids
            )
            if component.net_formal_charge != expected_metal_charge:
                raise ValidationError("metal_component_charge_mismatch", component.external_id)
        if component.atom_formal_charges_complete != expected_charge_completeness:
            raise ValidationError("atom_formal_charge_completeness_mismatch", component.external_id)
        _validate_chemical_topology_proof(component, atom_by_id, request.graph_revision)
        assignment_atoms.extend(component.atom_ids)
        assignment_by_id[component.external_id] = component
    if len(assignment_atoms) != len(set(assignment_atoms)):
        raise ValidationError("overlapping_atom_partition", "assignment components overlap", "assignment_components")

    projection_atoms: list[str] = []
    seen_parameterization_ids: set[str] = set()
    for projection in request.projections:
        atom_ids = set(projection.atom_ids)
        if not atom_ids or len(projection.atom_ids) != len(atom_ids) or not atom_ids <= set(atom_by_id):
            raise ValidationError("invalid_projection_atoms", projection.external_id, "projections")
        if projection.parameterization_component_id in seen_parameterization_ids:
            raise ValidationError("duplicate_parameterization_component_id", projection.parameterization_component_id)
        seen_parameterization_ids.add(projection.parameterization_component_id)
        component = assignment_by_id.get(projection.parameterization_component_id)
        if component is None or atom_ids != set(component.atom_ids):
            raise ValidationError("projection_assignment_mismatch", projection.external_id, "projections")
        projection_atoms.extend(projection.atom_ids)
        if projection.simulation_split and projection.simulation_split_reason is SimulationSplitReason.NONE_CONNECTED:
            raise ValidationError("missing_simulation_split_reason", projection.external_id, "projections")
        if not projection.simulation_split and projection.simulation_split_reason is not SimulationSplitReason.NONE_CONNECTED:
            raise ValidationError("unexpected_simulation_split_reason", projection.external_id, "projections")
        for atom_id in projection.atom_ids:
            atom = atom_by_id[atom_id]
            if (
                atom.chemical_component_id != projection.chemical_component_id
                or atom.canonical_residue_id != projection.canonical_residue_id
                or atom.simulation_component_id != projection.simulation_component_id
                or atom.simulation_residue_id != projection.simulation_residue_id
            ):
                raise ValidationError("projection_identity_mismatch", atom_id, f"projections.{projection.external_id}")
        allowed_scope_keys = atom_ids | {projection.external_id, projection.parameterization_component_id}
        for key, scopes in projection.included_scopes_by_atom_or_component.items():
            if key not in allowed_scope_keys or not scopes or not all(isinstance(scope, str) and scope for scope in scopes):
                raise ValidationError("invalid_included_scope", key, f"projections.{projection.external_id}")
            if key in atom_ids and not set(scopes) <= set(atom_by_id[key].scopes):
                raise ValidationError("scope_identity_mismatch", key, f"projections.{projection.external_id}")

    if len(projection_atoms) != len(set(projection_atoms)) or set(projection_atoms) != set(assignment_atoms):
        raise ValidationError("projection_partition_mismatch", "projections must partition assignment atoms exactly")
    if seen_parameterization_ids != component_ids:
        raise ValidationError("projection_assignment_mismatch", "each assignment component needs one projection")

    assignment_charge_by_residue: dict[str, int] = {}
    for component in request.assignment_components:
        residue_ids = {atom_by_id[atom_id].canonical_residue_id for atom_id in component.atom_ids}
        if len(residue_ids) != 1:
            raise ValidationError("cross_residue_assignment_component", component.external_id)
        residue_id = next(iter(residue_ids))
        assignment_charge_by_residue[residue_id] = (
            assignment_charge_by_residue.get(residue_id, 0) + component.net_formal_charge
        )
    prepared_charge_by_residue: dict[str, int] = {}
    for residue in request.topology.residues:
        if residue.canonical_residue_id not in assignment_charge_by_residue:
            continue
        prepared_charge_by_residue[residue.canonical_residue_id] = (
            prepared_charge_by_residue.get(residue.canonical_residue_id, 0) + residue.net_formal_charge
        )
    if assignment_charge_by_residue != prepared_charge_by_residue:
        raise ValidationError(
            "assignment_component_charge_partition_mismatch",
            "assignment component charges must conserve each canonical residue charge",
        )


def _validate_chemical_topology_proof(
    component: AssignmentComponent,
    atom_by_id: Mapping[str, PreparedAtom],
    graph_revision: int,
) -> None:
    proof = component.chemical_topology_proof
    path = f"assignment_components.{component.external_id}.chemical_topology_proof"
    if proof.schema_version != SCHEMA_VERSION or proof.graph_revision != graph_revision:
        raise ValidationError("stale_chemical_topology_proof", component.external_id, path)
    if not proof.perception_rule_version:
        raise ValidationError("missing_perception_rule_version", component.external_id, path)
    if not proof.hydrogen_resolution_method:
        raise ValidationError("missing_hydrogen_resolution_method", component.external_id, path)
    if proof.atom_ids != component.atom_ids or len(proof.atom_ids) != len(set(proof.atom_ids)):
        raise ValidationError("chemical_topology_proof_atom_mismatch", component.external_id, path)
    atom_ids = set(component.atom_ids)
    if (
        set(proof.explicit_hydrogen_count_by_atom) != atom_ids
        or set(proof.implicit_hydrogen_candidates_by_atom) != atom_ids
        or set(proof.valence_status_by_atom) != atom_ids
        or not set(proof.covalent_bond_order_by_atom) <= atom_ids
    ):
        raise ValidationError("incomplete_chemical_topology_proof", component.external_id, path)
    if len(proof.unresolved_valence_atom_ids) != len(set(proof.unresolved_valence_atom_ids)) or not set(
        proof.unresolved_valence_atom_ids
    ) <= atom_ids:
        raise ValidationError("invalid_unresolved_valence_atoms", component.external_id, path)
    for atom_id, count in proof.explicit_hydrogen_count_by_atom.items():
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValidationError("invalid_explicit_hydrogen_count", atom_id, path)
    for atom_id, candidates in proof.implicit_hydrogen_candidates_by_atom.items():
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in candidates):
            raise ValidationError("invalid_implicit_hydrogen_candidates", atom_id, path)
    for atom_id, order in proof.covalent_bond_order_by_atom.items():
        if isinstance(order, bool) or not isinstance(order, (int, float)) or not math.isfinite(order) or order < 0:
            raise ValidationError("invalid_covalent_bond_order", atom_id, path)
    if any(not isinstance(status, str) or not status for status in proof.valence_status_by_atom.values()):
        raise ValidationError("invalid_valence_status", component.external_id, path)
    status = proof.explicit_hydrogen_status
    unresolved = set(proof.unresolved_valence_atom_ids)
    metal_only = all(atom_by_id[atom_id].is_metal for atom_id in component.atom_ids)
    if status == "not_applicable":
        valid_status = metal_only and proof.valence_complete and not unresolved
    elif status == "complete":
        valid_status = (
            not metal_only
            and proof.valence_complete
            and not unresolved
            and all(0 in candidates for candidates in proof.implicit_hydrogen_candidates_by_atom.values())
        )
    elif status == "incomplete":
        valid_status = (
            not metal_only
            and proof.valence_complete
            and not unresolved
            and any(
                candidates and 0 not in candidates
                for candidates in proof.implicit_hydrogen_candidates_by_atom.values()
            )
        )
    elif status == "unresolved":
        valid_status = not metal_only and not proof.valence_complete and bool(unresolved)
    else:
        valid_status = False
    if not valid_status:
        raise ValidationError("inconsistent_explicit_hydrogen_status", status, path)
    if not proof.proof_hash or proof.proof_hash != proof.computed_hash():
        raise ValidationError("stale_chemical_topology_proof_hash", component.external_id, path)


def validate_input(request: MetalAssignmentInput) -> None:
    if request.schema_version != SCHEMA_VERSION:
        raise ValidationError("unsupported_schema_version", str(request.schema_version), "schema_version")
    if not request.request_id:
        raise ValidationError("missing_request_id", "request ID is required", "request_id")
    if request.interaction_model not in {"nonbonded_12_6", "bonded"}:
        raise ValidationError("unsupported_interaction_model", request.interaction_model, "interaction_model")
    validate_capability_snapshot(request.capability_snapshot, request.interaction_model)
    validate_topology(request.topology)
    if request.graph_revision != request.topology.graph_revision:
        raise ValidationError("stale_graph_revision", "request and topology graph revisions differ")
    if request.input_hash != request.topology.input_hash:
        raise ValidationError("stale_input_hash", "request and topology input hashes differ")
    validate_electronic_state(request.electronic_state, [atom.element for atom in request.topology.atoms])
    _validate_projection_mapping(request)
    _validate_partition_proofs(request)
    atoms_by_canonical_residue: dict[str, list[PreparedAtom]] = {}
    for atom in request.topology.atoms:
        atoms_by_canonical_residue.setdefault(atom.canonical_residue_id, []).append(atom)
    for atoms in atoms_by_canonical_residue.values():
        metal_count = sum(atom.is_metal for atom in atoms)
        if metal_count and len(atoms) > metal_count and not request.capability_snapshot.embedded_metal:
            raise ValidationError("unsupported_provider_capability", "provider does not support embedded metal")
        if metal_count and len(atoms) == metal_count and not request.capability_snapshot.standalone_metal:
            raise ValidationError("unsupported_provider_capability", "provider does not support standalone metal")
        if metal_count > 1 and not request.capability_snapshot.multi_metal:
            raise ValidationError("unsupported_provider_capability", "provider does not support multi-metal components")
    supported_base_providers = set(request.capability_snapshot.base_force_field_providers)
    for component in request.assignment_components:
        if component.base_force_field and component.provider not in supported_base_providers:
            raise ValidationError("unsupported_base_force_field_provider", component.provider)
        if (
            component.base_force_field
            and component.provider in {"gaff", "gaff2"}
            and component.chemical_topology_proof.explicit_hydrogen_status != "complete"
        ):
            raise ValidationError(
                "incomplete_explicit_hydrogen_topology",
                component.external_id,
                "assignment_components",
            )
    if request.charge_contract is not None:
        validate_charge_contract(request.charge_contract, request.topology)
    if request.partial_charge_artifacts is not None:
        validate_partial_charge_artifacts(request, request.partial_charge_artifacts)
    if not request.projection_hash:
        raise ValidationError("missing_projection_hash", "a computed projection hash is required")
    if request.projection_hash != request.computed_projection_hash():
        raise ValidationError("stale_projection_hash", "projection hash does not match the canonical payload")


def _validate_numeric_tree(value: Any, path: str) -> None:
    if isinstance(value, bool) or value is None:
        raise ValidationError("invalid_parameter_value", str(value), path)
    if isinstance(value, (int, float)):
        if not math.isfinite(value):
            raise ValidationError("nonfinite_parameter", str(value), path)
        return
    if isinstance(value, str):
        if not value:
            raise ValidationError("invalid_parameter_value", "empty string", path)
        return
    if isinstance(value, Mapping):
        if not value:
            raise ValidationError("invalid_parameter_value", "empty parameter object", path)
        for key, item in value.items():
            if not isinstance(key, str) or not key:
                raise ValidationError("invalid_parameter_key", str(key), path)
            _validate_numeric_tree(item, f"{path}.{key}")
        return
    if isinstance(value, (tuple, list)):
        if not value:
            raise ValidationError("invalid_parameter_value", "empty parameter sequence", path)
        for index, item in enumerate(value):
            _validate_numeric_tree(item, f"{path}[{index}]")
        return
    raise ValidationError("invalid_parameter_value", type(value).__name__, path)


def _active_adjacency(topology: PreparedChemicalTopology) -> dict[str, set[str]]:
    adjacency = {atom.external_id: set() for atom in topology.atoms}
    for edge in (*topology.bonds, *topology.links):
        if edge.active:
            atom1, atom2 = edge.atom_ids
            adjacency[atom1].add(atom2)
            adjacency[atom2].add(atom1)
    return adjacency


def _term_follows_topology(kind: str, atom_ids: tuple[str, ...], adjacency: Mapping[str, set[str]]) -> bool:
    expected_size = {
        "bond": 2,
        "distance_constraint": 2,
        "angle": 3,
        "proper_dihedral": 4,
        "improper_dihedral": 4,
    }.get(kind)
    if expected_size is None or len(atom_ids) != expected_size or len(set(atom_ids)) != len(atom_ids):
        return False
    if kind == "bond":
        return atom_ids[1] in adjacency[atom_ids[0]]
    if kind == "distance_constraint":
        return bool(adjacency[atom_ids[0]] & adjacency[atom_ids[1]])
    if kind == "angle":
        return atom_ids[0] in adjacency[atom_ids[1]] and atom_ids[2] in adjacency[atom_ids[1]]
    if kind == "proper_dihedral":
        return all(atom_ids[index + 1] in adjacency[atom_ids[index]] for index in range(3))
    center = atom_ids[1]
    return all(atom_id in adjacency[center] for atom_id in (atom_ids[0], atom_ids[2], atom_ids[3]))


def _validate_bonded_parameter_mapping(
    terms: Mapping[str, Any],
    *,
    topology: PreparedChemicalTopology,
    allowed_atom_ids: set[str],
    require_metal: bool,
    path: str,
) -> None:
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    adjacency = _active_adjacency(topology)
    for term_id, value in terms.items():
        if not isinstance(term_id, str) or not term_id or not isinstance(value, Mapping):
            raise ValidationError("invalid_bonded_parameter", str(term_id), path)
        required = {"kind", "atom_ids", "parameters", "source"}
        if set(value) != required:
            raise ValidationError("invalid_bonded_parameter_fields", term_id, path)
        raw_atom_ids = value["atom_ids"]
        if not isinstance(raw_atom_ids, (tuple, list)) or not all(
            isinstance(atom_id, str) and atom_id for atom_id in raw_atom_ids
        ):
            raise ValidationError("invalid_bonded_parameter_atoms", term_id, path)
        atom_ids = tuple(raw_atom_ids)
        if not set(atom_ids) <= allowed_atom_ids or not set(atom_ids) <= set(atom_by_id):
            raise ValidationError("bonded_parameter_coverage_mismatch", term_id, path)
        if require_metal and not any(atom_by_id[atom_id].is_metal for atom_id in atom_ids):
            raise ValidationError("nonmetal_term_in_metal_overlay", term_id, path)
        if not _term_follows_topology(str(value["kind"]), atom_ids, adjacency):
            raise ValidationError("bonded_parameter_topology_mismatch", term_id, path)
        if not isinstance(value["source"], str) or not value["source"]:
            raise ValidationError("missing_parameter_source", term_id, path)
        if not isinstance(value["parameters"], Mapping):
            raise ValidationError("invalid_bonded_parameter_value", term_id, path)
        _validate_numeric_tree(value["parameters"], f"{path}.{term_id}.parameters")


def _validate_atom_parameter_maps(
    *,
    atom_types: Mapping[str, str],
    charges: Mapping[str, float],
    lj_parameters: Mapping[str, Any],
    covered_atom_ids: set[str],
    path: str,
) -> None:
    if any(not isinstance(atom_id, str) or not atom_id for atom_id in covered_atom_ids):
        raise ValidationError("invalid_overlay_atom", path)
    for atom_id, atom_type in atom_types.items():
        if atom_id not in covered_atom_ids or not isinstance(atom_type, str) or not atom_type:
            raise ValidationError("overlay_atom_type_mismatch", atom_id, path)
    for atom_id, charge in charges.items():
        if atom_id not in covered_atom_ids or isinstance(charge, bool) or not isinstance(charge, (int, float)) or not math.isfinite(charge):
            raise ValidationError("invalid_overlay_charge", atom_id, path)
    for atom_id, parameters in lj_parameters.items():
        if atom_id not in covered_atom_ids or not isinstance(parameters, Mapping):
            raise ValidationError("invalid_lj_parameter", atom_id, path)
        if set(parameters) != {"epsilon", "rmin", "energy_unit", "length_unit", "source"}:
            raise ValidationError("invalid_lj_parameter_fields", atom_id, path)
        epsilon, rmin = parameters["epsilon"], parameters["rmin"]
        if (
            isinstance(epsilon, bool)
            or not isinstance(epsilon, (int, float))
            or not math.isfinite(epsilon)
            or epsilon < 0
            or isinstance(rmin, bool)
            or not isinstance(rmin, (int, float))
            or not math.isfinite(rmin)
            or rmin < 0
            or parameters["energy_unit"] != "kcal/mol"
            or parameters["length_unit"] != "angstrom"
            or not isinstance(parameters["source"], str)
            or not parameters["source"]
        ):
                raise ValidationError("invalid_lj_parameter_value", atom_id, path)


def _validate_mass_parameters(
    masses: Mapping[str, float],
    *,
    covered_atom_ids: set[str],
    path: str,
) -> None:
    for atom_id, mass in masses.items():
        if (
            atom_id not in covered_atom_ids
            or isinstance(mass, bool)
            or not isinstance(mass, (int, float))
            or not math.isfinite(mass)
            or mass <= 0
        ):
            raise ValidationError("invalid_overlay_mass", atom_id, path)


def _validate_charge_overlay_against_contract(
    request: MetalAssignmentInput,
    overlay: ChargeOverlay,
) -> None:
    contract = request.charge_contract
    if contract is None:
        if request.partial_charge_artifacts is None:
            raise ValidationError("unexpected_charge_overlay", "request has no charge artifacts or fit contract")
        return
    charges = overlay.charges
    if not set(contract.fit_atom_ids) <= set(charges):
        raise ValidationError("charge_overlay_scope_mismatch", "charge overlay must cover every fitted atom")
    component_by_id = {component.external_id: component for component in request.topology.components}
    tolerance = 1.0e-6
    if contract.policy is ChargePolicy.PRESERVE_COMPONENT_CHARGES:
        for component_id, target in contract.component_targets.items():
            fitted_ids = set(component_by_id[component_id].atom_ids) & set(contract.fit_atom_ids)
            if abs(sum(charges[atom_id] for atom_id in fitted_ids) - target) > tolerance:
                raise ValidationError("component_charge_conservation", component_id)
    elif contract.policy is ChargePolicy.SITE_JOINT_FIT:
        if abs(sum(charges[atom_id] for atom_id in contract.fit_atom_ids) - float(contract.site_target)) > tolerance:
            raise ValidationError("site_charge_conservation", "fitted charges do not match the site target")
    else:
        for group in contract.constraint_groups:
            if abs(sum(charges[atom_id] for atom_id in group.atom_ids) - group.target_charge) > tolerance:
                raise ValidationError("group_charge_conservation", group.external_id)


def validate_result(request: MetalAssignmentInput, result: ParameterizationResult) -> None:
    validate_input(request)
    if result.schema_version != SCHEMA_VERSION:
        raise ValidationError("unsupported_schema_version", str(result.schema_version), "result.schema_version")
    identity = (result.request_id, result.input_hash, result.graph_revision, result.topology_hash, result.projection_hash)
    expected = (
        request.request_id,
        request.input_hash,
        request.graph_revision,
        request.topology.topology_hash,
        request.projection_hash,
    )
    if identity != expected:
        raise ValidationError("result_identity_mismatch", "result does not belong to the validated request")
    overlays = [result.base_overlay, result.metal_overlay, result.charge_overlay, result.bonded_overlay]
    if any(overlay is not None and overlay.topology_hash != request.topology.topology_hash for overlay in overlays):
        raise ValidationError("overlay_topology_mismatch", "overlay targets a different topology")
    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    base_ids = set(result.base_overlay.covered_atom_ids)
    metal_ids = set(result.metal_overlay.covered_atom_ids)
    if len(base_ids) != len(result.base_overlay.covered_atom_ids) or len(metal_ids) != len(result.metal_overlay.covered_atom_ids):
        raise ValidationError("duplicate_overlay_atom", "overlay coverage must be unique")
    if not base_ids <= set(atom_by_id) or not metal_ids <= set(atom_by_id):
        raise ValidationError("overlay_unknown_atom", "overlay references an unknown atom")
    expected_base_ids = {
        atom_id
        for component in request.assignment_components
        if component.base_force_field
        for atom_id in component.atom_ids
    }
    target_metal_ids = set(atomic_center_atom_ids(request))
    preassigned_molecule = (
        result.provenance.get("base_assignment") == "preassigned_molecule"
    )
    if preassigned_molecule:
        if (
            base_ids
            or result.base_overlay.atom_types
            or result.base_overlay.charges
            or result.base_overlay.masses
            or result.base_overlay.lj_parameters
            or result.base_overlay.bonded_parameters
            or result.base_overlay.parameter_source != "xponge:preassigned-molecule"
        ):
            raise ValidationError(
                "invalid_local_patch_base_overlay",
                "local metal patches must not carry base-system parameters",
            )
    elif base_ids != expected_base_ids:
        raise ValidationError("base_overlay_coverage_mismatch", "base overlay must cover every base assignment atom")
    if not target_metal_ids <= metal_ids:
        raise ValidationError("metal_overlay_coverage_mismatch", "metal overlay must cover every metal atom")
    allowed_metal_overlay_ids = {
        atom_id for atom_id, atom in atom_by_id.items()
        if atom_id in target_metal_ids or "metalOverlay" in atom.scopes
    }
    if not metal_ids <= allowed_metal_overlay_ids:
        raise ValidationError(
            "metal_overlay_scope_mismatch",
            "nonmetal coverage requires explicit metalOverlay scope",
        )
    if not result.base_overlay.parameter_source or not result.metal_overlay.parameter_source:
        raise ValidationError("missing_parameter_source", "overlay sources are required")
    if not isinstance(result.metal_overlay.precedence, int) or isinstance(result.metal_overlay.precedence, bool) or result.metal_overlay.precedence <= 0:
        raise ValidationError("invalid_overlay_precedence", str(result.metal_overlay.precedence))
    _validate_atom_parameter_maps(
        atom_types=result.base_overlay.atom_types,
        charges=result.base_overlay.charges,
        lj_parameters=result.base_overlay.lj_parameters,
        covered_atom_ids=base_ids,
        path="base_overlay",
    )
    _validate_mass_parameters(
        result.base_overlay.masses,
        covered_atom_ids=base_ids,
        path="base_overlay.masses",
    )
    _validate_atom_parameter_maps(
        atom_types=result.metal_overlay.atom_types,
        charges=result.metal_overlay.charges,
        lj_parameters=result.metal_overlay.lj_parameters,
        covered_atom_ids=metal_ids,
        path="metal_overlay",
    )
    _validate_mass_parameters(
        result.metal_overlay.masses,
        covered_atom_ids=metal_ids,
        path="metal_overlay.masses",
    )
    _validate_bonded_parameter_mapping(
        result.base_overlay.bonded_parameters,
        topology=request.topology,
        allowed_atom_ids=base_ids,
        require_metal=False,
        path="base_overlay.bonded_parameters",
    )
    _validate_bonded_parameter_mapping(
        result.metal_overlay.bonded_parameters,
        topology=request.topology,
        allowed_atom_ids=metal_ids,
        require_metal=True,
        path="metal_overlay.bonded_parameters",
    )
    if result.charge_overlay is not None:
        if not result.charge_overlay.source:
            raise ValidationError("missing_parameter_source", "charge overlay source is required")
        if (
            not result.charge_overlay.overlay_hash
            or result.charge_overlay.overlay_hash != result.charge_overlay.computed_hash()
        ):
            raise ValidationError("stale_charge_overlay_hash", request.request_id)
        if set(result.charge_overlay.charges) != set(result.charge_overlay.atom_sources):
            raise ValidationError("charge_overlay_source_coverage_mismatch", request.request_id)
        if (
            len(result.charge_overlay.artifact_hashes) != len(set(result.charge_overlay.artifact_hashes))
            or any(not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None for value in result.charge_overlay.artifact_hashes)
        ):
            raise ValidationError("invalid_charge_overlay_artifact_hash", request.request_id)
        for atom_id, charge in result.charge_overlay.charges.items():
            if atom_id not in atom_by_id or isinstance(charge, bool) or not isinstance(charge, (int, float)) or not math.isfinite(charge):
                raise ValidationError("invalid_overlay_charge", atom_id, "charge_overlay")
            source = result.charge_overlay.atom_sources[atom_id]
            if not isinstance(source, Mapping) or set(source) != {
                "artifact_id", "artifact_hash", "precedence", "scope_kind", "scope_id",
                "atomic_charge_role", "provider", "provider_version", "method", "source",
            }:
                raise ValidationError("invalid_charge_overlay_source", atom_id)
            if (
                not source["artifact_id"]
                or source["artifact_hash"] not in result.charge_overlay.artifact_hashes
                or isinstance(source["precedence"], bool)
                or not isinstance(source["precedence"], int)
                or source["precedence"] <= 0
                or not source["provider"]
                or not source["provider_version"]
                or not source["method"]
                or not source["source"]
            ):
                raise ValidationError("invalid_charge_overlay_source", atom_id)
            if atom_id in result.metal_overlay.charges and source["precedence"] <= result.metal_overlay.precedence:
                raise ValidationError("charge_overlay_precedence_conflict", atom_id)
        if request.partial_charge_artifacts is not None:
            required_hashes = {artifact.artifact_hash for artifact in request.partial_charge_artifacts.artifacts}
            if not required_hashes <= set(result.charge_overlay.artifact_hashes):
                raise ValidationError("missing_partial_charge_artifact_provenance", request.request_id)
            selected_request_artifact: dict[str, PartialChargeArtifact] = {}
            for artifact in request.partial_charge_artifacts.artifacts:
                for atom_id in artifact.atom_ids:
                    previous = selected_request_artifact.get(atom_id)
                    if previous is None or artifact.precedence > previous.precedence:
                        selected_request_artifact[atom_id] = artifact
            for atom_id, artifact in selected_request_artifact.items():
                source = result.charge_overlay.atom_sources.get(atom_id)
                if source is None:
                    raise ValidationError("missing_partial_charge_artifact_atom", atom_id)
                artifact_charge = artifact.charges[artifact.atom_ids.index(atom_id)]
                if source["artifact_hash"] == artifact.artifact_hash:
                    if abs(result.charge_overlay.charges[atom_id] - artifact_charge) > 1.0e-12:
                        raise ValidationError("partial_charge_artifact_value_mismatch", atom_id)
                elif source["precedence"] <= artifact.precedence:
                    raise ValidationError("partial_charge_artifact_override_conflict", atom_id)
        _validate_charge_overlay_against_contract(request, result.charge_overlay)
    if result.bonded_overlay is not None:
        if not result.bonded_overlay.source:
            raise ValidationError("missing_parameter_source", "bonded overlay source is required")
        _validate_bonded_parameter_mapping(
            result.bonded_overlay.terms,
            topology=request.topology,
            allowed_atom_ids=set(atom_by_id),
            require_metal=request.interaction_model == "bonded",
            path="bonded_overlay.terms",
        )
    valid_statuses = {
        "base_assignment_completed", "fit_completed", "overlay_validated",
    }
    if result.status not in valid_statuses:
        raise ValidationError("invalid_result_status", result.status)
    if not result.provenance or not result.fit_reports:
        raise ValidationError("missing_result_provenance", "fit reports and provenance are required")
    if not all(isinstance(warning, str) and warning for warning in result.warnings):
        raise ValidationError("invalid_result_warning", "warnings must be non-empty strings")
    if result.complete:
        raise ValidationError(
            "invalid_result_completion",
            "parameterization results are local overlays, not full-system checkpoints",
        )
    if result.application_audit:
        raise ValidationError(
            "invalid_parameterization_application_audit",
            "Molecule application reports are returned by metal_assignment.apply",
        )
    if not result.result_hash:
        raise ValidationError("missing_result_hash", "a computed result hash is required")
    if result.result_hash != result.computed_hash():
        raise ValidationError("stale_result_hash", "result hash does not match canonical payload")


def _strict_object(
    value: Any,
    *,
    required: set[str],
    optional: set[str] = frozenset(),
    path: str,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValidationError("invalid_wire_type", "expected object", path)
    keys = set(value)
    missing = required - keys
    unknown = keys - required - optional
    if missing:
        raise ValidationError("missing_wire_field", ",".join(sorted(missing)), path)
    if unknown:
        raise ValidationError("unknown_wire_field", ",".join(sorted(unknown)), path)
    return value


def _strings(value: Any, path: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValidationError("invalid_wire_type", "expected string array", path)
    return tuple(value)


def _parse_atom(value: Any, path: str) -> PreparedAtom:
    data = _strict_object(
        value,
        required={
            "external_id", "canonical_atom_id", "stable_order", "element", "atom_name", "coordinates", "chemical_component_id",
            "canonical_residue_id", "simulation_component_id", "simulation_residue_id", "formal_charge",
            "formal_charge_known", "partial_charge", "is_metal", "scopes",
        },
        path=path,
    )
    coordinates = data["coordinates"]
    if not isinstance(coordinates, list) or len(coordinates) != 3:
        raise ValidationError("invalid_wire_type", "coordinates must have three values", f"{path}.coordinates")
    return PreparedAtom(
        external_id=data["external_id"], canonical_atom_id=data["canonical_atom_id"],
        stable_order=data["stable_order"], element=data["element"],
        coordinates=tuple(coordinates), chemical_component_id=data["chemical_component_id"],
        canonical_residue_id=data["canonical_residue_id"], simulation_component_id=data["simulation_component_id"],
        simulation_residue_id=data["simulation_residue_id"], formal_charge=data["formal_charge"],
        formal_charge_known=data["formal_charge_known"],
        partial_charge=data["partial_charge"], is_metal=data["is_metal"], scopes=_strings(data["scopes"], f"{path}.scopes"),
        atom_name=data["atom_name"],
    )


def _parse_residue(value: Any, path: str) -> PreparedResidue:
    data = _strict_object(
        value,
        required={
            "external_id", "canonical_residue_id", "atom_ids", "chemical_component_id",
            "simulation_component_id", "net_formal_charge", "charge_resolution_method",
            "atom_formal_charges_complete", "residue_name", "chain_id", "res_seq", "icode",
            "polymer_position", "template_id",
        },
        path=path,
    )
    return PreparedResidue(
        data["external_id"], data["canonical_residue_id"], _strings(data["atom_ids"], f"{path}.atom_ids"),
        data["chemical_component_id"], data["simulation_component_id"], data["net_formal_charge"],
        data["charge_resolution_method"], data["atom_formal_charges_complete"],
        data["residue_name"], data["chain_id"], data["res_seq"], data["icode"],
        data["polymer_position"], data["template_id"],
    )


def _parse_component(value: Any, path: str) -> PreparedComponent:
    data = _strict_object(
        value,
        required={
            "external_id", "atom_ids", "chemical_component_id", "simulation_component_id",
            "net_formal_charge", "charge_resolution_method", "atom_formal_charges_complete",
        },
        path=path,
    )
    return PreparedComponent(
        data["external_id"], _strings(data["atom_ids"], f"{path}.atom_ids"),
        data["chemical_component_id"], data["simulation_component_id"], data["net_formal_charge"],
        data["charge_resolution_method"], data["atom_formal_charges_complete"],
    )


def _parse_bond(value: Any, path: str) -> PreparedBond:
    data = _strict_object(value, required={"external_id", "atom_ids", "order", "semantic", "active", "source"}, path=path)
    atom_ids = _strings(data["atom_ids"], f"{path}.atom_ids")
    return PreparedBond(data["external_id"], atom_ids, data["order"], data["semantic"], data["active"], data["source"])


def _parse_link(value: Any, path: str) -> PreparedLink:
    data = _strict_object(value, required={"external_id", "atom_ids", "kind", "active", "confirmed", "source"}, path=path)
    atom_ids = _strings(data["atom_ids"], f"{path}.atom_ids")
    return PreparedLink(data["external_id"], atom_ids, data["kind"], data["active"], data["confirmed"], data["source"])


def _parse_topology(value: Any) -> PreparedChemicalTopology:
    data = _strict_object(
        value,
        required={"schema_version", "graph_revision", "input_hash", "coordinate_unit", "atoms", "residues", "components", "bonds", "links", "topology_hash"},
        path="topology",
    )
    if data["schema_version"] != SCHEMA_VERSION:
        raise ValidationError("unsupported_schema_version", str(data["schema_version"]), "topology.schema_version")
    for name in ("atoms", "residues", "components", "bonds", "links"):
        if not isinstance(data[name], list):
            raise ValidationError("invalid_wire_type", "expected array", f"topology.{name}")
    return PreparedChemicalTopology(
        data["schema_version"], data["graph_revision"], data["input_hash"], data["coordinate_unit"],
        tuple(_parse_atom(item, f"topology.atoms[{index}]") for index, item in enumerate(data["atoms"])),
        tuple(_parse_residue(item, f"topology.residues[{index}]") for index, item in enumerate(data["residues"])),
        tuple(_parse_component(item, f"topology.components[{index}]") for index, item in enumerate(data["components"])),
        tuple(_parse_bond(item, f"topology.bonds[{index}]") for index, item in enumerate(data["bonds"])),
        tuple(_parse_link(item, f"topology.links[{index}]") for index, item in enumerate(data["links"])),
        data["topology_hash"],
    )


def _parse_projection(value: Any, path: str) -> ProviderProjection:
    data = _strict_object(
        value,
        required={
            "external_id", "chemical_component_id", "canonical_residue_id", "parameterization_component_id",
            "simulation_component_id", "simulation_residue_id", "atom_ids", "included_scopes_by_atom_or_component",
            "assignment_only_split", "simulation_split", "simulation_split_reason", "parent_component_id",
        },
        path=path,
    )
    scopes = data["included_scopes_by_atom_or_component"]
    if not isinstance(scopes, Mapping):
        raise ValidationError("invalid_wire_type", "expected object", f"{path}.included_scopes_by_atom_or_component")
    try:
        split_reason = SimulationSplitReason(data["simulation_split_reason"])
    except ValueError as exc:
        raise ValidationError("invalid_simulation_split_reason", str(data["simulation_split_reason"]), path) from exc
    return ProviderProjection(
        data["external_id"], data["chemical_component_id"], data["canonical_residue_id"],
        data["parameterization_component_id"], data["simulation_component_id"], data["simulation_residue_id"],
        _strings(data["atom_ids"], f"{path}.atom_ids"),
        {key: _strings(item, f"{path}.included_scopes_by_atom_or_component.{key}") for key, item in scopes.items()},
        data["assignment_only_split"], data["simulation_split"], split_reason, data["parent_component_id"],
    )


def _parse_assignment_component(value: Any, path: str) -> AssignmentComponent:
    data = _strict_object(
        value,
        required={
            "external_id", "atom_ids", "classification", "provider", "net_formal_charge",
            "charge_resolution_method", "atom_formal_charges_complete", "base_force_field",
            "chemical_topology_proof",
        },
        path=path,
    )
    proof_path = f"{path}.chemical_topology_proof"
    proof_data = _strict_object(
        data["chemical_topology_proof"],
        required={
            "schema_version", "graph_revision", "perception_rule_version", "atom_ids",
            "explicit_hydrogen_count_by_atom", "implicit_hydrogen_candidates_by_atom",
            "covalent_bond_order_by_atom", "valence_status_by_atom",
            "unresolved_valence_atom_ids", "explicit_hydrogen_status", "hydrogen_resolution_method",
            "valence_complete", "proof_hash",
        },
        path=proof_path,
    )
    explicit_counts = proof_data["explicit_hydrogen_count_by_atom"]
    implicit_candidates = proof_data["implicit_hydrogen_candidates_by_atom"]
    covalent_orders = proof_data["covalent_bond_order_by_atom"]
    valence_statuses = proof_data["valence_status_by_atom"]
    if not all(isinstance(value, Mapping) for value in (
        explicit_counts, implicit_candidates, covalent_orders, valence_statuses,
    )):
        raise ValidationError("invalid_wire_type", "hydrogen proof maps must be objects", proof_path)
    proof = ChemicalTopologyProof(
        proof_data["schema_version"],
        proof_data["graph_revision"],
        proof_data["perception_rule_version"],
        _strings(proof_data["atom_ids"], f"{proof_path}.atom_ids"),
        dict(explicit_counts),
        {
            key: tuple(values) if isinstance(values, list) else values
            for key, values in implicit_candidates.items()
        },
        dict(covalent_orders),
        dict(valence_statuses),
        _strings(proof_data["unresolved_valence_atom_ids"], f"{proof_path}.unresolved_valence_atom_ids"),
        proof_data["explicit_hydrogen_status"],
        proof_data["hydrogen_resolution_method"],
        proof_data["valence_complete"],
        proof_data["proof_hash"],
    )
    return AssignmentComponent(
        data["external_id"],
        _strings(data["atom_ids"], f"{path}.atom_ids"),
        data["classification"],
        data["provider"],
        data["net_formal_charge"],
        data["charge_resolution_method"],
        data["atom_formal_charges_complete"],
        proof,
        data["base_force_field"],
    )


def _parse_charge_contract(value: Any) -> ChargeAssignmentContract:
    data = _strict_object(
        value,
        required={
            "policy", "fit_atom_ids", "component_formal_charges", "component_targets", "site_target",
            "fixed_atom_ids", "constraint_groups", "charge_unit", "source", "fit_precedence",
        },
        path="charge_contract",
    )
    groups = data["constraint_groups"]
    if not isinstance(groups, list):
        raise ValidationError("invalid_wire_type", "expected array", "charge_contract.constraint_groups")
    parsed_groups = []
    for index, value in enumerate(groups):
        path = f"charge_contract.constraint_groups[{index}]"
        group = _strict_object(value, required={"external_id", "atom_ids", "target_charge", "source"}, path=path)
        parsed_groups.append(ChargeConstraintGroup(group["external_id"], _strings(group["atom_ids"], f"{path}.atom_ids"), group["target_charge"], group["source"]))
    if not isinstance(data["component_formal_charges"], Mapping) or not isinstance(data["component_targets"], Mapping):
        raise ValidationError("invalid_wire_type", "component charge fields must be objects", "charge_contract")
    try:
        policy = ChargePolicy(data["policy"])
    except ValueError as exc:
        raise ValidationError("invalid_charge_policy", str(data["policy"]), "charge_contract.policy") from exc
    return ChargeAssignmentContract(
        policy, _strings(data["fit_atom_ids"], "charge_contract.fit_atom_ids"), data["component_formal_charges"],
        data["component_targets"], data["site_target"], _strings(data["fixed_atom_ids"], "charge_contract.fixed_atom_ids"),
        tuple(parsed_groups), data["charge_unit"], data["source"], data["fit_precedence"],
    )


def _parse_partial_charge_artifact(value: Any, path: str) -> PartialChargeArtifact:
    data = _strict_object(
        value,
        required={
            "schema_version", "artifact_id", "topology_hash", "graph_revision", "input_hash",
            "atom_ids", "charges", "scope_kind", "scope_id", "atomic_charge_role", "precedence",
            "provider", "provider_version", "method", "charge_unit", "source", "artifact_hash",
        },
        path=path,
    )
    if not isinstance(data["charges"], list):
        raise ValidationError("invalid_wire_type", "expected array", f"{path}.charges")
    return PartialChargeArtifact(
        data["schema_version"], data["artifact_id"], data["topology_hash"], data["graph_revision"],
        data["input_hash"], _strings(data["atom_ids"], f"{path}.atom_ids"), tuple(data["charges"]),
        data["scope_kind"], data["scope_id"], data["atomic_charge_role"], data["precedence"],
        data["provider"], data["provider_version"], data["method"], data["charge_unit"],
        data["source"], data["artifact_hash"],
    )


def _parse_partial_charge_bundle(value: Any) -> PartialChargeArtifactBundle:
    data = _strict_object(
        value,
        required={"schema_version", "topology_hash", "graph_revision", "input_hash", "artifacts", "bundle_hash"},
        path="partial_charge_artifacts",
    )
    if not isinstance(data["artifacts"], list):
        raise ValidationError("invalid_wire_type", "expected array", "partial_charge_artifacts.artifacts")
    return PartialChargeArtifactBundle(
        data["schema_version"], data["topology_hash"], data["graph_revision"], data["input_hash"],
        tuple(
            _parse_partial_charge_artifact(item, f"partial_charge_artifacts.artifacts[{index}]")
            for index, item in enumerate(data["artifacts"])
        ),
        data["bundle_hash"],
    )


def partial_charge_artifact_bundle_from_dict(value: Any) -> PartialChargeArtifactBundle:
    """Parse the strict source-neutral partial-charge bundle wire format."""

    return _parse_partial_charge_bundle(value)


def _parse_partition_proof(value: Any, path: str) -> ResiduePartitionProof:
    data = _strict_object(
        value,
        required={
            "canonical_residue_id", "source_atom_ids", "active_internal_edge_ids",
            "removed_interaction_edge_ids", "simulation_partitions", "proof_hash",
        },
        path=path,
    )
    partitions = data["simulation_partitions"]
    if not isinstance(partitions, Mapping):
        raise ValidationError("invalid_wire_type", "simulation partitions must be an object", f"{path}.simulation_partitions")
    return ResiduePartitionProof(
        data["canonical_residue_id"],
        _strings(data["source_atom_ids"], f"{path}.source_atom_ids"),
        _strings(data["active_internal_edge_ids"], f"{path}.active_internal_edge_ids"),
        _strings(data["removed_interaction_edge_ids"], f"{path}.removed_interaction_edge_ids"),
        {key: _strings(item, f"{path}.simulation_partitions.{key}") for key, item in partitions.items()},
        data["proof_hash"],
    )


def _parse_capability_snapshot(value: Any) -> ProviderCapabilitySnapshot:
    data = _strict_object(
        value,
        required={
            "schema_version", "provider_name", "provider_version", "provider_revision", "projection_schema_version",
            "embedded_metal", "standalone_metal", "multi_metal", "bonded", "nonbonded_12_6",
            "base_force_field_providers",
        },
        path="capability_snapshot",
    )
    return ProviderCapabilitySnapshot(
        data["schema_version"], data["provider_name"], data["provider_version"], data["provider_revision"],
        data["projection_schema_version"], data["embedded_metal"], data["standalone_metal"], data["multi_metal"],
        data["bonded"], data["nonbonded_12_6"],
        _strings(data["base_force_field_providers"], "capability_snapshot.base_force_field_providers"),
    )


def metal_assignment_input_to_dict(request: MetalAssignmentInput) -> dict[str, Any]:
    """Return the canonical schema-v1 snake_case wire representation."""

    return _canonicalize(request)


def metal_assignment_input_from_dict(value: Any) -> MetalAssignmentInput:
    """Strictly parse, version-check, and validate a schema-v1 request."""

    data = _strict_object(
        value,
        required={
            "schema_version", "request_id", "interaction_model", "electronic_state", "graph_revision", "input_hash",
            "topology", "projections", "assignment_components", "partition_proofs", "charge_contract",
            "partial_charge_artifacts", "capability_snapshot", "projection_hash",
        },
        path="request",
    )
    if data["schema_version"] != SCHEMA_VERSION:
        raise ValidationError("unsupported_schema_version", str(data["schema_version"]), "schema_version")
    electronic = _strict_object(data["electronic_state"], required={"net_charge", "spin_multiplicity"}, path="electronic_state")
    if (
        not isinstance(data["projections"], list)
        or not isinstance(data["assignment_components"], list)
        or not isinstance(data["partition_proofs"], list)
    ):
        raise ValidationError("invalid_wire_type", "projection fields must be arrays", "request")
    request = MetalAssignmentInput(
        data["schema_version"], data["request_id"], data["interaction_model"],
        ElectronicState(electronic["net_charge"], electronic["spin_multiplicity"]),
        data["graph_revision"], data["input_hash"], _parse_topology(data["topology"]),
        tuple(_parse_projection(item, f"projections[{index}]") for index, item in enumerate(data["projections"])),
        tuple(_parse_assignment_component(item, f"assignment_components[{index}]") for index, item in enumerate(data["assignment_components"])),
        tuple(_parse_partition_proof(item, f"partition_proofs[{index}]") for index, item in enumerate(data["partition_proofs"])),
        _parse_capability_snapshot(data["capability_snapshot"]),
        None if data["charge_contract"] is None else _parse_charge_contract(data["charge_contract"]),
        None if data["partial_charge_artifacts"] is None else _parse_partial_charge_bundle(data["partial_charge_artifacts"]),
        data["projection_hash"],
    )
    validate_input(request)
    return request


def metal_assignment_input_dumps(request: MetalAssignmentInput) -> str:
    validate_input(request)
    return json.dumps(metal_assignment_input_to_dict(request), sort_keys=True, separators=(",", ":"), allow_nan=False)


def metal_assignment_input_loads(payload: str) -> MetalAssignmentInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return metal_assignment_input_from_dict(value)


def _wire_mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or not all(isinstance(key, str) and key for key in value):
        raise ValidationError("invalid_wire_type", "expected string-keyed object", path)
    return value


def _parse_base_overlay(value: Any) -> BaseForceFieldOverlay:
    data = _strict_object(
        value,
        required={
            "topology_hash", "covered_atom_ids", "atom_types", "charges", "masses", "lj_parameters",
            "bonded_parameters", "parameter_source",
        },
        path="result.base_overlay",
    )
    return BaseForceFieldOverlay(
        data["topology_hash"],
        _strings(data["covered_atom_ids"], "result.base_overlay.covered_atom_ids"),
        _wire_mapping(data["atom_types"], "result.base_overlay.atom_types"),
        _wire_mapping(data["charges"], "result.base_overlay.charges"),
        _wire_mapping(data["masses"], "result.base_overlay.masses"),
        _wire_mapping(data["lj_parameters"], "result.base_overlay.lj_parameters"),
        _wire_mapping(data["bonded_parameters"], "result.base_overlay.bonded_parameters"),
        data["parameter_source"],
    )


def _parse_metal_overlay(value: Any) -> MetalParameterOverlay:
    data = _strict_object(
        value,
        required={
            "topology_hash", "covered_atom_ids", "atom_types", "charges", "masses", "lj_parameters",
            "bonded_parameters", "parameter_source", "precedence",
        },
        path="result.metal_overlay",
    )
    return MetalParameterOverlay(
        data["topology_hash"],
        _strings(data["covered_atom_ids"], "result.metal_overlay.covered_atom_ids"),
        _wire_mapping(data["atom_types"], "result.metal_overlay.atom_types"),
        _wire_mapping(data["charges"], "result.metal_overlay.charges"),
        _wire_mapping(data["masses"], "result.metal_overlay.masses"),
        _wire_mapping(data["lj_parameters"], "result.metal_overlay.lj_parameters"),
        _wire_mapping(data["bonded_parameters"], "result.metal_overlay.bonded_parameters"),
        data["parameter_source"],
        data["precedence"],
    )


def _parse_charge_overlay(value: Any) -> ChargeOverlay:
    data = _strict_object(
        value,
        required={"topology_hash", "charges", "source", "atom_sources", "artifact_hashes", "overlay_hash"},
        path="result.charge_overlay",
    )
    return ChargeOverlay(
        data["topology_hash"],
        _wire_mapping(data["charges"], "result.charge_overlay.charges"),
        data["source"],
        _wire_mapping(data["atom_sources"], "result.charge_overlay.atom_sources"),
        _strings(data["artifact_hashes"], "result.charge_overlay.artifact_hashes"),
        data["overlay_hash"],
    )


def _parse_bonded_overlay(value: Any) -> BondedParameterOverlay:
    data = _strict_object(value, required={"topology_hash", "terms", "source"}, path="result.bonded_overlay")
    return BondedParameterOverlay(
        data["topology_hash"],
        _wire_mapping(data["terms"], "result.bonded_overlay.terms"),
        data["source"],
    )


def parameterization_result_to_dict(
    request: MetalAssignmentInput,
    result: ParameterizationResult,
) -> dict[str, Any]:
    validate_result(request, result)
    return _canonicalize(result)


def parameterization_result_from_dict(
    value: Any,
    request: MetalAssignmentInput,
) -> ParameterizationResult:
    result = _parameterization_result_from_dict_unchecked(value)
    validate_result(request, result)
    return result


def _parameterization_result_from_dict_unchecked(
    value: Any,
) -> ParameterizationResult:
    """Parse the result wire shape without requiring a full-system request.

    This private helper exists for the local ``MetalParameterPatch`` parser.
    The patch validator supplies the identity, hash and strictly-local overlay
    checks that ``validate_result`` normally derives from a
    ``MetalAssignmentInput``.
    """

    data = _strict_object(
        value,
        required={
            "schema_version", "request_id", "input_hash", "graph_revision", "topology_hash",
            "projection_hash", "base_overlay", "metal_overlay", "charge_overlay", "bonded_overlay",
            "fit_reports", "provenance", "warnings", "status", "application_audit", "complete",
            "result_hash",
        },
        path="result",
    )
    result = ParameterizationResult(
        schema_version=data["schema_version"],
        request_id=data["request_id"],
        input_hash=data["input_hash"],
        graph_revision=data["graph_revision"],
        topology_hash=data["topology_hash"],
        projection_hash=data["projection_hash"],
        base_overlay=_parse_base_overlay(data["base_overlay"]),
        metal_overlay=_parse_metal_overlay(data["metal_overlay"]),
        charge_overlay=None if data["charge_overlay"] is None else _parse_charge_overlay(data["charge_overlay"]),
        bonded_overlay=None if data["bonded_overlay"] is None else _parse_bonded_overlay(data["bonded_overlay"]),
        fit_reports=_wire_mapping(data["fit_reports"], "result.fit_reports"),
        provenance=_wire_mapping(data["provenance"], "result.provenance"),
        warnings=_strings(data["warnings"], "result.warnings"),
        status=data["status"],
        application_audit=_wire_mapping(data["application_audit"], "result.application_audit"),
        complete=data["complete"],
        result_hash=data["result_hash"],
    )
    return result


def parameterization_result_dumps(request: MetalAssignmentInput, result: ParameterizationResult) -> str:
    return json.dumps(
        parameterization_result_to_dict(request, result),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def parameterization_result_loads(payload: str, request: MetalAssignmentInput) -> ParameterizationResult:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return parameterization_result_from_dict(value, request)


__all__ = [
    "AssignmentComponent", "BASE_FORCE_FIELD_SCOPE", "BaseForceFieldOverlay", "BondedParameterOverlay",
    "COORDINATE_UNIT", "ChargeAssignmentContract", "ChargeConstraintGroup", "ChargeOverlay", "ChargePolicy",
    "ElectronicState", "FrozenMapping", "MetalAssignmentInput", "MetalParameterOverlay", "ParameterizationResult",
    "PartialChargeArtifact", "PartialChargeArtifactBundle",
    "PreparedAtom", "PreparedBond", "PreparedChemicalTopology", "PreparedComponent", "PreparedLink", "PreparedResidue",
    "ProviderCapabilitySnapshot", "ProviderProjection", "ResiduePartitionProof", "SCHEMA_VERSION",
    "SimulationSplitReason", "ValidationError",
    "metal_assignment_input_dumps", "metal_assignment_input_from_dict", "metal_assignment_input_loads",
    "metal_assignment_input_to_dict", "parameterization_result_dumps", "parameterization_result_from_dict",
    "partial_charge_artifact_bundle_from_dict",
    "parameterization_result_loads", "parameterization_result_to_dict", "validate_charge_contract",
    "validate_partial_charge_artifacts",
    "atomic_center_atom_ids", "base_metal_component_supported",
    "validate_electronic_state", "validate_input",
    "validate_capability_snapshot", "validate_result", "validate_topology",
]
