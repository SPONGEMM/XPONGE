"""Topology-preserving application of validated parameter overlays.

The parent process only handles immutable, hash-closed data.  A short-lived
worker imports Xponge force-field modules, materializes a real ``Molecule`` and
verifies that every prepared atom/residue/edge survives unchanged.  This keeps
Xponge's process-global type registries out of callers and makes retries
transactional by construction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
import subprocess
from typing import Any, Mapping

from .contracts import (
    MetalAssignmentInput,
    ParameterizationResult,
    PreparedChemicalTopology,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _parse_topology,
    _sha256,
    _strict_object,
    atomic_center_atom_ids,
    validate_input,
    validate_result,
    validate_topology,
)
from ._worker_runtime import worker_command


ASSIGNED_SYSTEM_SCHEMA_VERSION = 1
APPLY_WORKER_PROTOCOL_VERSION = 1


@dataclass(frozen=True, slots=True)
class ForceRealizationProtocol:
    """Explicit nonbonded/exclusion rules used to realize an assigned graph."""

    family: str
    exclusion_bond_depth: int
    one_four_lj_scale: float
    one_four_electrostatic_scale: float
    source: str
    protocol_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "family": self.family,
            "exclusion_bond_depth": self.exclusion_bond_depth,
            "one_four_lj_scale": self.one_four_lj_scale,
            "one_four_electrostatic_scale": self.one_four_electrostatic_scale,
            "source": self.source,
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ForceRealizationProtocol":
        return replace(self, protocol_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class AssignedAtomParameters:
    external_id: str
    effective_type: str
    charge: float
    mass: float
    epsilon: float
    rmin: float
    sources: Mapping[str, str]

    def __post_init__(self) -> None:
        _freeze_field(self, "sources")


@dataclass(frozen=True, slots=True)
class AssignedForceTerm:
    external_id: str
    kind: str
    atom_ids: tuple[str, ...]
    parameters: Mapping[str, Any]
    source: str

    def __post_init__(self) -> None:
        _freeze_field(self, "parameters")


@dataclass(frozen=True, slots=True)
class AssignedSystemArtifact:
    """Portable recipe for an assigned Xponge system.

    It embeds the prepared topology so later SPONGE export can materialize the
    same system in another isolated process without Chemcore or Mokda imports.
    """

    schema_version: int
    request_id: str
    input_hash: str
    graph_revision: int
    topology: PreparedChemicalTopology
    projection_hash: str
    source_result_hash: str
    atom_parameters: tuple[AssignedAtomParameters, ...]
    force_terms: tuple[AssignedForceTerm, ...]
    force_protocol: ForceRealizationProtocol
    application_audit: Mapping[str, Any]
    complete: bool
    assigned_system_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "application_audit")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in self.__dataclass_fields__.values()
            if item.name != "assigned_system_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "AssignedSystemArtifact":
        return replace(self, assigned_system_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class ParameterApplicationOutput:
    result: ParameterizationResult
    assigned_system: AssignedSystemArtifact


def validate_force_realization_protocol(protocol: ForceRealizationProtocol) -> None:
    if protocol.family != "amber":
        raise ValidationError("unsupported_force_realization_family", protocol.family)
    if protocol.exclusion_bond_depth != 3:
        raise ValidationError("unsupported_exclusion_bond_depth", str(protocol.exclusion_bond_depth))
    for name in ("one_four_lj_scale", "one_four_electrostatic_scale"):
        value = getattr(protocol, name)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValidationError("invalid_force_realization_scale", name)
    if not protocol.source:
        raise ValidationError("missing_force_realization_source", "source")
    if not protocol.protocol_hash or protocol.protocol_hash != protocol.computed_hash():
        raise ValidationError("stale_force_realization_protocol_hash", protocol.family)


def _canonical_term_key(kind: str, atom_ids: tuple[str, ...]) -> tuple[str, tuple[str, ...]]:
    if kind in {"bond", "distance_constraint"}:
        return kind, tuple(sorted(atom_ids))
    if kind == "angle":
        return kind, min(atom_ids, atom_ids[::-1])
    if kind == "proper_dihedral":
        return kind, min(atom_ids, atom_ids[::-1])
    return kind, atom_ids


def _merge_atom_parameters(
    request: MetalAssignmentInput,
    result: ParameterizationResult,
) -> tuple[AssignedAtomParameters, ...]:
    base = result.base_overlay
    metal = result.metal_overlay
    charge_overlay = result.charge_overlay
    atoms: list[AssignedAtomParameters] = []
    missing: dict[str, list[str]] = {}
    namespace = result.result_hash[:12]
    for atom in sorted(request.topology.atoms, key=lambda item: item.stable_order):
        atom_id = atom.external_id
        properties: dict[str, Any] = {}
        sources: dict[str, str] = {}
        for name, values in (
            ("type", base.atom_types), ("charge", base.charges), ("mass", base.masses), ("lj", base.lj_parameters),
        ):
            if atom_id in values:
                properties[name] = values[atom_id]
                sources[name] = base.parameter_source
        for name, values in (
            ("type", metal.atom_types), ("charge", metal.charges), ("mass", metal.masses), ("lj", metal.lj_parameters),
        ):
            if atom_id in values:
                properties[name] = values[atom_id]
                sources[name] = metal.parameter_source
        if charge_overlay is not None and atom_id in charge_overlay.charges:
            properties["charge"] = charge_overlay.charges[atom_id]
            sources["charge"] = charge_overlay.source
        absent = [name for name in ("type", "charge", "mass", "lj") if name not in properties]
        if absent:
            missing[atom_id] = absent
            continue
        lj = properties["lj"]
        atoms.append(AssignedAtomParameters(
            external_id=atom_id,
            effective_type=f"MA_{namespace}_{atom.stable_order}",
            charge=float(properties["charge"]),
            mass=float(properties["mass"]),
            epsilon=float(lj["epsilon"]),
            rmin=float(lj["rmin"]),
            sources=sources,
        ))
    if missing:
        detail = ";".join(f"{atom_id}:{','.join(names)}" for atom_id, names in sorted(missing.items()))
        raise ValidationError("incomplete_atom_parameters", detail, "parameterization_result")
    return tuple(atoms)


def _merge_force_terms(result: ParameterizationResult) -> tuple[tuple[AssignedForceTerm, ...], tuple[str, ...]]:
    merged: dict[tuple[str, tuple[str, ...]], AssignedForceTerm] = {}
    superseded: list[str] = []
    layers = (
        (0, result.base_overlay.bonded_parameters),
        (1, result.metal_overlay.bonded_parameters),
        (2, result.bonded_overlay.terms if result.bonded_overlay is not None else {}),
    )
    for _, values in layers:
        for term_id, value in values.items():
            term = AssignedForceTerm(
                str(term_id), str(value["kind"]), tuple(value["atom_ids"]),
                value["parameters"], str(value["source"]),
            )
            key = _canonical_term_key(term.kind, term.atom_ids)
            previous = merged.get(key)
            if previous is not None:
                superseded.append(previous.external_id)
            merged[key] = term
    ordered = tuple(sorted(merged.values(), key=lambda item: (_canonical_term_key(item.kind, item.atom_ids), item.external_id)))
    return ordered, tuple(sorted(superseded))


def _validate_edge_coverage(topology: PreparedChemicalTopology, terms: tuple[AssignedForceTerm, ...]) -> None:
    covered_pairs = {
        frozenset(term.atom_ids)
        for term in terms
        if term.kind in {"bond", "distance_constraint"} and len(term.atom_ids) == 2
    }
    missing = [
        edge.external_id
        for edge in (*topology.bonds, *topology.links)
        if edge.active and frozenset(edge.atom_ids) not in covered_pairs
    ]
    if missing:
        raise ValidationError("missing_bonded_edge_parameter", ",".join(sorted(missing)), "force_terms")


def validate_assigned_system_artifact(artifact: AssignedSystemArtifact) -> None:
    if artifact.schema_version != ASSIGNED_SYSTEM_SCHEMA_VERSION:
        raise ValidationError("unsupported_assigned_system_schema", str(artifact.schema_version))
    validate_topology(artifact.topology)
    validate_force_realization_protocol(artifact.force_protocol)
    if not artifact.request_id or not artifact.source_result_hash or not artifact.projection_hash:
        raise ValidationError("missing_assigned_system_identity", artifact.request_id)
    if artifact.input_hash != artifact.topology.input_hash or artifact.graph_revision != artifact.topology.graph_revision:
        raise ValidationError("assigned_system_identity_mismatch", artifact.request_id)
    expected_atom_ids = tuple(atom.external_id for atom in sorted(artifact.topology.atoms, key=lambda item: item.stable_order))
    actual_atom_ids = tuple(atom.external_id for atom in artifact.atom_parameters)
    if actual_atom_ids != expected_atom_ids or len(actual_atom_ids) != len(set(actual_atom_ids)):
        raise ValidationError("assigned_atom_order_mismatch", artifact.request_id)
    for atom in artifact.atom_parameters:
        if not atom.effective_type or set(atom.sources) != {"type", "charge", "mass", "lj"}:
            raise ValidationError("incomplete_assigned_atom_provenance", atom.external_id)
        for value in (atom.charge, atom.mass, atom.epsilon, atom.rmin):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValidationError("invalid_assigned_atom_parameter", atom.external_id)
        if atom.mass <= 0 or atom.epsilon < 0 or atom.rmin < 0:
            raise ValidationError("invalid_assigned_atom_parameter", atom.external_id)
    term_ids = [term.external_id for term in artifact.force_terms]
    if len(term_ids) != len(set(term_ids)):
        raise ValidationError("duplicate_assigned_force_term", artifact.request_id)
    _validate_edge_coverage(artifact.topology, artifact.force_terms)
    if artifact.complete != bool(artifact.application_audit):
        raise ValidationError("assigned_system_completion_mismatch", artifact.request_id)
    if not artifact.assigned_system_hash or artifact.assigned_system_hash != artifact.computed_hash():
        raise ValidationError("stale_assigned_system_hash", artifact.request_id)


def _worker_payload(artifact: AssignedSystemArtifact) -> dict[str, Any]:
    return {
        "protocol_version": APPLY_WORKER_PROTOCOL_VERSION,
        "assigned_system": assigned_system_to_dict(artifact),
    }


def _invoke_apply_worker(payload: Mapping[str, Any], timeout_seconds: float) -> Mapping[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._apply_worker"),
            input=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("parameter_apply_worker_timeout", str(timeout_seconds)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError("invalid_parameter_apply_worker_output", completed.stderr[-1000:] or str(exc)) from exc
    if completed.returncode != 0 or not isinstance(response, Mapping) or response.get("ok") is not True:
        error = response.get("error", {}) if isinstance(response, Mapping) else {}
        message = error.get("message", "worker failed") if isinstance(error, Mapping) else "worker failed"
        if isinstance(error, Mapping) and error.get("traceback"):
            message = f"{message}; traceback={str(error['traceback'])[-3000:]}"
        if completed.stderr.strip():
            message = f"{message}; stderr={completed.stderr.strip()[-1000:]}"
        raise ValidationError("parameter_apply_worker_failed", message)
    expected_fields = {"ok", "protocol_version", "assigned_system_hash", "application_audit", "response_hash"}
    if set(response) != expected_fields or response["protocol_version"] != APPLY_WORKER_PROTOCOL_VERSION:
        raise ValidationError("invalid_parameter_apply_worker_response", "unexpected fields or version")
    unhashed = dict(response)
    response_hash = unhashed.pop("response_hash")
    if response_hash != _sha256(unhashed):
        raise ValidationError("stale_parameter_apply_worker_hash", "response hash mismatch")
    if response["assigned_system_hash"] != payload["assigned_system"]["assigned_system_hash"]:
        raise ValidationError("parameter_apply_worker_identity_mismatch", "assigned system hash mismatch")
    if not isinstance(response["application_audit"], Mapping) or not response["application_audit"]:
        raise ValidationError("missing_application_audit", "worker returned no audit")
    return response["application_audit"]


def apply_parameterization(
    request: MetalAssignmentInput,
    result: ParameterizationResult,
    force_protocol: ForceRealizationProtocol,
    *,
    timeout_seconds: float = 120.0,
) -> ParameterApplicationOutput:
    """Apply a validated result in an isolated Xponge process."""

    validate_result(request, result)
    if result.complete or result.status not in {"fit_completed", "overlay_validated"}:
        raise ValidationError("invalid_apply_source_status", result.status)
    validate_force_realization_protocol(force_protocol)
    atoms = _merge_atom_parameters(request, result)
    terms, superseded = _merge_force_terms(result)
    _validate_edge_coverage(request.topology, terms)
    candidate = AssignedSystemArtifact(
        schema_version=ASSIGNED_SYSTEM_SCHEMA_VERSION,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology=request.topology,
        projection_hash=request.projection_hash,
        source_result_hash=result.result_hash,
        atom_parameters=atoms,
        force_terms=terms,
        force_protocol=force_protocol,
        application_audit={},
        complete=False,
    ).with_computed_hash()
    audit = dict(_invoke_apply_worker(_worker_payload(candidate), timeout_seconds))
    audit["superseded_term_ids"] = list(superseded)
    final_artifact = replace(candidate, application_audit=audit, complete=True, assigned_system_hash="").with_computed_hash()
    validate_assigned_system_artifact(final_artifact)
    complete_result = replace(
        result,
        status="assigned_system_validated",
        application_audit={
            "assigned_system_hash": final_artifact.assigned_system_hash,
            **audit,
        },
        complete=True,
        result_hash="",
    ).with_computed_hash()
    validate_result(request, complete_result)
    return ParameterApplicationOutput(complete_result, final_artifact)


def apply_base_force_field(
    request: MetalAssignmentInput,
    base_overlay: Any,
    force_protocol: ForceRealizationProtocol,
    *,
    timeout_seconds: float = 120.0,
) -> AssignedSystemArtifact:
    """Materialize a metal-free base assignment without a metal overlay."""

    validate_input(request)
    validate_force_realization_protocol(force_protocol)
    if atomic_center_atom_ids(request):
        raise ValidationError(
            "metal_in_base_force_field_application",
            "use metal-assignment for systems containing metals",
        )
    topology = request.topology
    ordered_topology_atoms = tuple(sorted(topology.atoms, key=lambda item: item.stable_order))
    atom_ids = tuple(atom.external_id for atom in ordered_topology_atoms)
    expected = set(atom_ids)
    if base_overlay.topology_hash != topology.topology_hash:
        raise ValidationError("base_overlay_topology_mismatch", request.request_id)
    if set(base_overlay.covered_atom_ids) != expected or len(base_overlay.covered_atom_ids) != len(expected):
        raise ValidationError("base_overlay_coverage_mismatch", request.request_id)
    for name in ("atom_types", "charges", "masses", "lj_parameters"):
        if set(getattr(base_overlay, name)) != expected:
            raise ValidationError("incomplete_base_overlay_parameters", name)

    source_hash = _sha256({
        "request_id": request.request_id,
        "projection_hash": request.projection_hash,
        "base_overlay": _canonicalize(base_overlay),
        "force_protocol_hash": force_protocol.protocol_hash,
    })
    namespace = source_hash[:12]
    atoms = []
    for atom in ordered_topology_atoms:
        atom_id = atom.external_id
        lj = base_overlay.lj_parameters[atom_id]
        atoms.append(AssignedAtomParameters(
            external_id=atom_id,
            effective_type=f"FA_{namespace}_{atom.stable_order}",
            charge=float(base_overlay.charges[atom_id]),
            mass=float(base_overlay.masses[atom_id]),
            epsilon=float(lj["epsilon"]),
            rmin=float(lj["rmin"]),
            sources={
                "type": base_overlay.parameter_source,
                "charge": base_overlay.parameter_source,
                "mass": base_overlay.parameter_source,
                "lj": base_overlay.parameter_source,
            },
        ))
    terms = tuple(sorted((
        AssignedForceTerm(
            str(term_id),
            str(value["kind"]),
            tuple(value["atom_ids"]),
            value["parameters"],
            str(value["source"]),
        )
        for term_id, value in base_overlay.bonded_parameters.items()
    ), key=lambda item: (_canonical_term_key(item.kind, item.atom_ids), item.external_id)))
    _validate_edge_coverage(topology, terms)
    candidate = AssignedSystemArtifact(
        schema_version=ASSIGNED_SYSTEM_SCHEMA_VERSION,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology=topology,
        projection_hash=request.projection_hash,
        source_result_hash=source_hash,
        atom_parameters=tuple(atoms),
        force_terms=terms,
        force_protocol=force_protocol,
        application_audit={},
        complete=False,
    ).with_computed_hash()
    audit = dict(_invoke_apply_worker(_worker_payload(candidate), timeout_seconds))
    audit["assignment_mode"] = "base_force_field"
    final_artifact = replace(
        candidate,
        application_audit=audit,
        complete=True,
        assigned_system_hash="",
    ).with_computed_hash()
    validate_assigned_system_artifact(final_artifact)
    return final_artifact


def assigned_system_to_dict(artifact: AssignedSystemArtifact) -> dict[str, Any]:
    if artifact.complete:
        validate_assigned_system_artifact(artifact)
    elif not artifact.assigned_system_hash or artifact.assigned_system_hash != artifact.computed_hash():
        raise ValidationError("stale_assigned_system_hash", artifact.request_id)
    return _canonicalize(artifact)


def _parse_force_protocol(value: Any) -> ForceRealizationProtocol:
    data = _strict_object(
        value,
        required={
            "family", "exclusion_bond_depth", "one_four_lj_scale",
            "one_four_electrostatic_scale", "source", "protocol_hash",
        },
        path="assigned_system.force_protocol",
    )
    protocol = ForceRealizationProtocol(**data)
    validate_force_realization_protocol(protocol)
    return protocol


def assigned_system_from_dict(value: Any) -> AssignedSystemArtifact:
    data = _strict_object(
        value,
        required={
            "schema_version", "request_id", "input_hash", "graph_revision", "topology", "projection_hash",
            "source_result_hash", "atom_parameters", "force_terms", "force_protocol", "application_audit",
            "complete", "assigned_system_hash",
        },
        path="assigned_system",
    )
    if not isinstance(data["atom_parameters"], list) or not isinstance(data["force_terms"], list):
        raise ValidationError("invalid_wire_type", "assigned atom/term fields must be arrays", "assigned_system")
    atom_parameters = []
    for index, value in enumerate(data["atom_parameters"]):
        item = _strict_object(
            value,
            required={"external_id", "effective_type", "charge", "mass", "epsilon", "rmin", "sources"},
            path=f"assigned_system.atom_parameters[{index}]",
        )
        if not isinstance(item["sources"], Mapping):
            raise ValidationError("invalid_wire_type", "sources must be an object")
        atom_parameters.append(AssignedAtomParameters(**item))
    force_terms = []
    for index, value in enumerate(data["force_terms"]):
        item = _strict_object(
            value,
            required={"external_id", "kind", "atom_ids", "parameters", "source"},
            path=f"assigned_system.force_terms[{index}]",
        )
        if not isinstance(item["atom_ids"], list) or not isinstance(item["parameters"], Mapping):
            raise ValidationError("invalid_wire_type", "invalid force term")
        force_terms.append(AssignedForceTerm(
            item["external_id"], item["kind"], tuple(item["atom_ids"]), item["parameters"], item["source"],
        ))
    if not isinstance(data["application_audit"], Mapping):
        raise ValidationError("invalid_wire_type", "application audit must be an object")
    artifact = AssignedSystemArtifact(
        schema_version=data["schema_version"], request_id=data["request_id"], input_hash=data["input_hash"],
        graph_revision=data["graph_revision"], topology=_parse_topology(data["topology"]),
        projection_hash=data["projection_hash"], source_result_hash=data["source_result_hash"],
        atom_parameters=tuple(atom_parameters), force_terms=tuple(force_terms),
        force_protocol=_parse_force_protocol(data["force_protocol"]),
        application_audit=data["application_audit"], complete=data["complete"],
        assigned_system_hash=data["assigned_system_hash"],
    )
    if artifact.complete:
        validate_assigned_system_artifact(artifact)
    elif not artifact.assigned_system_hash or artifact.assigned_system_hash != artifact.computed_hash():
        raise ValidationError("stale_assigned_system_hash", artifact.request_id)
    return artifact


def assigned_system_dumps(artifact: AssignedSystemArtifact) -> str:
    return json.dumps(assigned_system_to_dict(artifact), sort_keys=True, separators=(",", ":"), allow_nan=False)


def assigned_system_loads(payload: str) -> AssignedSystemArtifact:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return assigned_system_from_dict(value)


__all__ = [
    "ASSIGNED_SYSTEM_SCHEMA_VERSION", "AssignedAtomParameters", "AssignedForceTerm",
    "AssignedSystemArtifact", "ForceRealizationProtocol", "ParameterApplicationOutput",
    "apply_parameterization", "assigned_system_dumps", "assigned_system_from_dict",
    "assigned_system_loads", "assigned_system_to_dict", "validate_assigned_system_artifact",
    "validate_force_realization_protocol",
]
