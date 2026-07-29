"""Small, hash-closed metal parameter patches for an existing Molecule."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
import math
from typing import Any, Mapping, Sequence

from .contracts import (
    MetalAssignmentInput,
    ParameterizationResult,
    SCHEMA_VERSION,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _parameterization_result_from_dict_unchecked,
    _sha256,
    validate_result,
)


PATCH_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class PatchAtomIdentity:
    external_id: str
    canonical_atom_id: int
    stable_order: int
    element: str
    atom_name: str
    canonical_residue_id: str
    simulation_residue_id: str


@dataclass(frozen=True, slots=True)
class PatchLinkIdentity:
    external_id: str
    atom_ids: tuple[str, str]
    kind: str


@dataclass(frozen=True, slots=True)
class MetalParameterPatch:
    schema_version: int
    request_id: str
    input_hash: str
    graph_revision: int
    source_topology_hash: str
    source_local_topology_hash: str
    source_local_geometry_hash: str
    target_atom_identity_hash: str
    site_ids: tuple[str, ...]
    target_metal_atom_ids: tuple[str, ...]
    atoms: tuple[PatchAtomIdentity, ...]
    required_links: tuple[PatchLinkIdentity, ...]
    parameterization_result: ParameterizationResult
    provenance: Mapping[str, Any]
    patch_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "provenance")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "patch_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "MetalParameterPatch":
        return replace(self, patch_hash=self.computed_hash())


def _local_force_terms(result: ParameterizationResult) -> tuple[Mapping[str, Any], ...]:
    terms: list[Mapping[str, Any]] = []
    terms.extend(result.metal_overlay.bonded_parameters.values())
    if result.bonded_overlay is not None:
        terms.extend(result.bonded_overlay.terms.values())
    return tuple(terms)


def _target_atom_ids(result: ParameterizationResult) -> set[str]:
    atom_ids = set(result.metal_overlay.covered_atom_ids)
    atom_ids.update(result.metal_overlay.atom_types)
    atom_ids.update(result.metal_overlay.charges)
    atom_ids.update(result.metal_overlay.masses)
    atom_ids.update(result.metal_overlay.lj_parameters)
    if result.charge_overlay is not None:
        atom_ids.update(result.charge_overlay.charges)
    for term in _local_force_terms(result):
        atom_ids.update(str(atom_id) for atom_id in term.get("atom_ids", ()))
    return atom_ids


def _local_topology_hash(
    atoms: Sequence[PatchAtomIdentity],
    links: Sequence[PatchLinkIdentity],
    result: ParameterizationResult,
) -> str:
    terms = [
        {
            "kind": str(term.get("kind") or ""),
            "atom_ids": [str(atom_id) for atom_id in term.get("atom_ids", ())],
        }
        for term in _local_force_terms(result)
    ]
    return _sha256({
        "atoms": [
            {
                "external_id": atom.external_id,
                "canonical_atom_id": atom.canonical_atom_id,
                "element": atom.element,
                "canonical_residue_id": atom.canonical_residue_id,
                "simulation_residue_id": atom.simulation_residue_id,
            }
            for atom in atoms
        ],
        "required_links": _canonicalize(tuple(links)),
        "terms": terms,
    })


def _local_geometry_hash(
    request: MetalAssignmentInput,
    atom_ids: Sequence[str],
) -> str:
    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    ordered = [atom_by_id[atom_id] for atom_id in atom_ids]
    distances: list[float] = []
    for left_index, left in enumerate(ordered):
        for right in ordered[left_index + 1:]:
            squared = sum(
                (float(left.coordinates[axis]) - float(right.coordinates[axis])) ** 2
                for axis in range(3)
            )
            distances.append(round(math.sqrt(squared), 10))
    return _sha256({
        "atom_ids": list(atom_ids),
        "pairwise_distances_angstrom": distances,
    })


def build_metal_parameter_patch(
    request: MetalAssignmentInput,
    result: ParameterizationResult,
    *,
    site_ids: Sequence[str] = (),
    provenance: Mapping[str, Any] | None = None,
) -> MetalParameterPatch:
    """Project a validated parameterization result into a local-only patch."""

    validate_result(request, result)
    if request.interaction_model != "bonded":
        raise ValidationError(
            "invalid_local_patch_interaction_model",
            request.interaction_model,
        )
    if result.provenance.get("base_assignment") != "preassigned_molecule":
        raise ValidationError(
            "local_patch_contains_base_assignment",
            str(result.provenance.get("base_assignment") or ""),
        )
    target_ids = _target_atom_ids(result)
    atom_by_id = {atom.external_id: atom for atom in request.topology.atoms}
    if not target_ids or not target_ids <= set(atom_by_id):
        raise ValidationError(
            "invalid_local_patch_atom_coverage",
            ",".join(sorted(target_ids - set(atom_by_id))),
        )
    ordered_ids = tuple(sorted(
        target_ids,
        key=lambda atom_id: (atom_by_id[atom_id].stable_order, atom_id),
    ))
    atoms = tuple(
        PatchAtomIdentity(
            external_id=atom.external_id,
            canonical_atom_id=atom.canonical_atom_id,
            stable_order=atom.stable_order,
            element=atom.element,
            atom_name=atom.atom_name,
            canonical_residue_id=atom.canonical_residue_id,
            simulation_residue_id=atom.simulation_residue_id,
        )
        for atom in (atom_by_id[atom_id] for atom_id in ordered_ids)
    )
    links = tuple(
        PatchLinkIdentity(link.external_id, link.atom_ids, link.kind)
        for link in request.topology.links
        if link.active and set(link.atom_ids) <= target_ids
    )
    metal_ids = tuple(
        atom.external_id
        for atom in atoms
        if atom_by_id[atom.external_id].is_metal
    )
    identity_hash = _sha256({
        "atoms": _canonicalize(atoms),
        "target_metal_atom_ids": list(metal_ids),
    })
    patch = MetalParameterPatch(
        schema_version=PATCH_SCHEMA_VERSION,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        source_topology_hash=request.topology.topology_hash,
        source_local_topology_hash=_local_topology_hash(atoms, links, result),
        source_local_geometry_hash=_local_geometry_hash(request, ordered_ids),
        target_atom_identity_hash=identity_hash,
        site_ids=tuple(sorted({str(site_id) for site_id in site_ids if str(site_id)})),
        target_metal_atom_ids=metal_ids,
        atoms=atoms,
        required_links=links,
        parameterization_result=result,
        provenance={
            "source_result_hash": result.result_hash,
            "projection_hash": request.projection_hash,
            **dict(provenance or {}),
        },
    ).with_computed_hash()
    validate_metal_parameter_patch(patch)
    return patch


def validate_metal_parameter_patch(patch: MetalParameterPatch) -> None:
    if patch.schema_version != PATCH_SCHEMA_VERSION:
        raise ValidationError(
            "unsupported_metal_parameter_patch_schema",
            str(patch.schema_version),
        )
    result = patch.parameterization_result
    if (
        result.schema_version != SCHEMA_VERSION
        or result.request_id != patch.request_id
        or result.input_hash != patch.input_hash
        or result.graph_revision != patch.graph_revision
        or result.topology_hash != patch.source_topology_hash
        or result.status != "overlay_validated"
        or result.complete
        or not result.result_hash
        or result.result_hash != result.computed_hash()
    ):
        raise ValidationError(
            "metal_parameter_patch_result_identity_mismatch",
            patch.request_id,
        )
    if (
        result.provenance.get("base_assignment") != "preassigned_molecule"
        or result.base_overlay.covered_atom_ids
        or result.base_overlay.atom_types
        or result.base_overlay.charges
        or result.base_overlay.masses
        or result.base_overlay.lj_parameters
        or result.base_overlay.bonded_parameters
    ):
        raise ValidationError(
            "metal_parameter_patch_contains_base_system",
            patch.request_id,
        )
    atom_ids = tuple(atom.external_id for atom in patch.atoms)
    if (
        not atom_ids
        or len(atom_ids) != len(set(atom_ids))
        or tuple(sorted(patch.target_metal_atom_ids)) != tuple(
            sorted(set(patch.target_metal_atom_ids))
        )
        or not set(patch.target_metal_atom_ids) <= set(atom_ids)
        or not _target_atom_ids(result) <= set(atom_ids)
    ):
        raise ValidationError(
            "invalid_metal_parameter_patch_atom_identity",
            patch.request_id,
        )
    if patch.target_atom_identity_hash != _sha256({
        "atoms": _canonicalize(patch.atoms),
        "target_metal_atom_ids": list(patch.target_metal_atom_ids),
    }):
        raise ValidationError(
            "stale_metal_parameter_patch_atom_identity",
            patch.request_id,
        )
    if patch.source_local_topology_hash != _local_topology_hash(
        patch.atoms,
        patch.required_links,
        result,
    ):
        raise ValidationError(
            "stale_metal_parameter_patch_local_topology",
            patch.request_id,
        )
    if not patch.patch_hash or patch.patch_hash != patch.computed_hash():
        raise ValidationError("stale_metal_parameter_patch_hash", patch.request_id)


def metal_parameter_patch_to_dict(
    patch: MetalParameterPatch,
) -> dict[str, Any]:
    validate_metal_parameter_patch(patch)
    return _canonicalize(patch)


def metal_parameter_patch_from_dict(value: Any) -> MetalParameterPatch:
    if not isinstance(value, Mapping):
        raise ValidationError("invalid_metal_parameter_patch", "expected object")
    required = {
        "schema_version", "request_id", "input_hash", "graph_revision",
        "source_topology_hash", "source_local_topology_hash",
        "source_local_geometry_hash", "target_atom_identity_hash", "site_ids",
        "target_metal_atom_ids", "atoms", "required_links",
        "parameterization_result", "provenance", "patch_hash",
    }
    missing = required - set(value)
    unknown = set(value) - required
    if missing or unknown:
        raise ValidationError(
            "invalid_metal_parameter_patch_fields",
            f"missing={sorted(missing)},unknown={sorted(unknown)}",
        )
    atoms = tuple(
        PatchAtomIdentity(
            external_id=item["external_id"],
            canonical_atom_id=item["canonical_atom_id"],
            stable_order=item["stable_order"],
            element=item["element"],
            atom_name=item["atom_name"],
            canonical_residue_id=item["canonical_residue_id"],
            simulation_residue_id=item["simulation_residue_id"],
        )
        for item in value["atoms"]
    )
    links = tuple(
        PatchLinkIdentity(
            external_id=item["external_id"],
            atom_ids=tuple(item["atom_ids"]),
            kind=item["kind"],
        )
        for item in value["required_links"]
    )
    patch = MetalParameterPatch(
        schema_version=value["schema_version"],
        request_id=value["request_id"],
        input_hash=value["input_hash"],
        graph_revision=value["graph_revision"],
        source_topology_hash=value["source_topology_hash"],
        source_local_topology_hash=value["source_local_topology_hash"],
        source_local_geometry_hash=value["source_local_geometry_hash"],
        target_atom_identity_hash=value["target_atom_identity_hash"],
        site_ids=tuple(value["site_ids"]),
        target_metal_atom_ids=tuple(value["target_metal_atom_ids"]),
        atoms=atoms,
        required_links=links,
        parameterization_result=_parameterization_result_from_dict_unchecked(
            value["parameterization_result"]
        ),
        provenance=value["provenance"],
        patch_hash=value["patch_hash"],
    )
    validate_metal_parameter_patch(patch)
    return patch


def metal_parameter_patch_dumps(patch: MetalParameterPatch) -> str:
    return json.dumps(
        metal_parameter_patch_to_dict(patch),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def metal_parameter_patch_loads(payload: str) -> MetalParameterPatch:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return metal_parameter_patch_from_dict(value)


__all__ = [
    "PATCH_SCHEMA_VERSION",
    "MetalParameterPatch",
    "PatchAtomIdentity",
    "PatchLinkIdentity",
    "build_metal_parameter_patch",
    "metal_parameter_patch_dumps",
    "metal_parameter_patch_from_dict",
    "metal_parameter_patch_loads",
    "metal_parameter_patch_to_dict",
    "validate_metal_parameter_patch",
]
