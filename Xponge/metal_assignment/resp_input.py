"""Hash-closed scientific and numerical protocol for RESP charge fitting."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

from .contracts import ValidationError, _canonicalize, _freeze_field, _sha256, _strict_object


RESP_FIT_INPUT_SCHEMA_VERSION = 3
SUPPORTED_RESP_BACKENDS = frozenset({"pyscf", "psi4"})
SUPPORTED_ESP_CHUNK_POLICIES = frozenset({"auto", "full", "grid", "dual", "pointwise"})
SUPPORTED_METAL_BASIS_POLICIES = frozenset({"require_ecp"})


@dataclass(frozen=True, slots=True)
class RespEquivalenceGroup:
    """Model-local atom equivalence supplied by the topology/model producer."""

    model_id: str
    atom_ids: tuple[str, ...]
    source: str


@dataclass(frozen=True, slots=True)
class RespLinearConstraint:
    """Model-local linear constraint expressed only in atom IDs and numbers."""

    model_id: str
    constraint_id: str
    role: str
    atom_ids: tuple[str, ...]
    coefficients: tuple[float, ...]
    target_charge: float
    source: str


@dataclass(frozen=True, slots=True)
class RespFitInput:
    """Versioned RESP settings, independent of a Python environment path."""

    schema_version: int
    backend: str
    basis_family: str
    metal_basis_policy: str
    basis_source: str
    optimize_geometry: bool
    grid_density: float
    grid_cell_layer: int
    radius_overrides: Mapping[str, float]
    restraint_a1: float
    restraint_a2: float
    two_stage: bool
    only_esp: bool
    esp_memory_limit_bytes: int
    esp_chunk_policy: str
    esp_safety_factor: float
    equivalence_groups: tuple[RespEquivalenceGroup, ...]
    linear_constraints: tuple[RespLinearConstraint, ...]
    source: str
    fit_input_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "radius_overrides")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "fit_input_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "RespFitInput":
        return replace(self, fit_input_hash=self.computed_hash())


def validate_resp_fit_input(value: RespFitInput) -> None:
    if value.schema_version != RESP_FIT_INPUT_SCHEMA_VERSION:
        raise ValidationError("unsupported_resp_fit_input_version", str(value.schema_version), "resp_fit_input")
    if value.backend not in SUPPORTED_RESP_BACKENDS:
        raise ValidationError("unsupported_resp_backend", value.backend, "resp_fit_input.backend")
    if not isinstance(value.basis_family, str) or not value.basis_family.strip():
        raise ValidationError(
            "invalid_resp_basis_family", str(value.basis_family), "resp_fit_input.basis_family"
        )
    if value.metal_basis_policy not in SUPPORTED_METAL_BASIS_POLICIES:
        raise ValidationError(
            "unsupported_resp_metal_basis_policy",
            str(value.metal_basis_policy),
            "resp_fit_input.metal_basis_policy",
        )
    if not isinstance(value.basis_source, str) or not value.basis_source.strip():
        raise ValidationError(
            "missing_resp_basis_source", str(value.basis_source), "resp_fit_input.basis_source"
        )
    if value.optimize_geometry is not False:
        raise ValidationError(
            "resp_geometry_must_remain_locked",
            "derived-model coordinates cannot be optimized inside the charge provider",
            "resp_fit_input.optimize_geometry",
        )
    if value.only_esp is not False or value.two_stage is not True:
        raise ValidationError(
            "incomplete_resp_protocol",
            "the metal-assignment RESP provider requires the restrained two-stage fit",
            "resp_fit_input",
        )
    if (
        isinstance(value.grid_density, bool)
        or not isinstance(value.grid_density, (int, float))
        or not math.isfinite(value.grid_density)
        or value.grid_density <= 0
    ):
        raise ValidationError("invalid_resp_grid_density", str(value.grid_density), "resp_fit_input.grid_density")
    if (
        isinstance(value.grid_cell_layer, bool)
        or not isinstance(value.grid_cell_layer, int)
        or value.grid_cell_layer <= 0
    ):
        raise ValidationError("invalid_resp_grid_layer", str(value.grid_cell_layer), "resp_fit_input.grid_cell_layer")
    for name in ("restraint_a1", "restraint_a2"):
        number = getattr(value, name)
        if isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(number) or number <= 0:
            raise ValidationError("invalid_resp_restraint", str(number), f"resp_fit_input.{name}")
    if (
        isinstance(value.esp_memory_limit_bytes, bool)
        or not isinstance(value.esp_memory_limit_bytes, int)
        or value.esp_memory_limit_bytes <= 0
    ):
        raise ValidationError(
            "invalid_resp_memory_limit",
            str(value.esp_memory_limit_bytes),
            "resp_fit_input.esp_memory_limit_bytes",
        )
    if value.esp_chunk_policy not in SUPPORTED_ESP_CHUNK_POLICIES:
        raise ValidationError(
            "invalid_resp_chunk_policy", value.esp_chunk_policy, "resp_fit_input.esp_chunk_policy"
        )
    if (
        isinstance(value.esp_safety_factor, bool)
        or not isinstance(value.esp_safety_factor, (int, float))
        or not math.isfinite(value.esp_safety_factor)
        or not 0 < value.esp_safety_factor <= 1
    ):
        raise ValidationError(
            "invalid_resp_safety_factor", str(value.esp_safety_factor), "resp_fit_input.esp_safety_factor"
        )
    if not isinstance(value.radius_overrides, Mapping):
        raise ValidationError("invalid_resp_radius_overrides", "expected mapping", "resp_fit_input.radius_overrides")
    for element, radius in value.radius_overrides.items():
        if not isinstance(element, str) or not re.fullmatch(r"[A-Z][a-z]?", element):
            raise ValidationError("invalid_resp_radius_element", str(element), "resp_fit_input.radius_overrides")
        if isinstance(radius, bool) or not isinstance(radius, (int, float)) or not math.isfinite(radius) or radius <= 0:
            raise ValidationError("invalid_resp_radius", str(radius), f"resp_fit_input.radius_overrides.{element}")
    group_keys: set[tuple[str, tuple[str, ...]]] = set()
    atoms_by_model: dict[str, set[str]] = {}
    for index, group in enumerate(value.equivalence_groups):
        path = f"resp_fit_input.equivalence_groups[{index}]"
        if not isinstance(group, RespEquivalenceGroup):
            raise ValidationError("invalid_resp_equivalence_group", type(group).__name__, path)
        if not group.model_id or not group.source:
            raise ValidationError("invalid_resp_equivalence_metadata", group.model_id, path)
        if len(group.atom_ids) < 2 or len(group.atom_ids) != len(set(group.atom_ids)) or any(not atom_id for atom_id in group.atom_ids):
            raise ValidationError("invalid_resp_equivalence_group", group.model_id, path)
        key = (group.model_id, group.atom_ids)
        if key in group_keys:
            raise ValidationError("duplicate_resp_equivalence_group", group.model_id, path)
        group_keys.add(key)
        seen = atoms_by_model.setdefault(group.model_id, set())
        overlap = seen & set(group.atom_ids)
        if overlap:
            raise ValidationError("overlapping_resp_equivalence_groups", ",".join(sorted(overlap)), path)
        seen.update(group.atom_ids)
    constraint_keys: set[tuple[str, str]] = set()
    for index, constraint in enumerate(value.linear_constraints):
        path = f"resp_fit_input.linear_constraints[{index}]"
        if not isinstance(constraint, RespLinearConstraint):
            raise ValidationError(
                "invalid_resp_linear_constraint", type(constraint).__name__, path
            )
        if (
            not constraint.model_id
            or not constraint.constraint_id
            or not constraint.role
            or not constraint.source
            or not constraint.atom_ids
            or len(constraint.atom_ids) != len(constraint.coefficients)
            or len(constraint.atom_ids) != len(set(constraint.atom_ids))
            or any(not atom_id for atom_id in constraint.atom_ids)
            or any(
                isinstance(coefficient, bool)
                or not isinstance(coefficient, (int, float))
                or not math.isfinite(coefficient)
                for coefficient in constraint.coefficients
            )
            or not any(abs(float(coefficient)) > 0.0 for coefficient in constraint.coefficients)
            or isinstance(constraint.target_charge, bool)
            or not isinstance(constraint.target_charge, (int, float))
            or not math.isfinite(constraint.target_charge)
        ):
            raise ValidationError(
                "invalid_resp_linear_constraint", constraint.constraint_id, path
            )
        key = (constraint.model_id, constraint.constraint_id)
        if key in constraint_keys:
            raise ValidationError(
                "duplicate_resp_linear_constraint", constraint.constraint_id, path
            )
        constraint_keys.add(key)
    if not isinstance(value.source, str) or not value.source:
        raise ValidationError("missing_resp_protocol_source", "source is required", "resp_fit_input.source")
    if not value.fit_input_hash or value.fit_input_hash != value.computed_hash():
        raise ValidationError("stale_resp_fit_input_hash", "RESP fit input hash mismatch", "resp_fit_input")


def resp_fit_input_from_dict(value: Any) -> RespFitInput:
    data = _strict_object(
        value,
        required={
            "schema_version", "backend", "basis_family", "metal_basis_policy", "basis_source",
            "optimize_geometry", "grid_density",
            "grid_cell_layer", "radius_overrides", "restraint_a1", "restraint_a2", "two_stage",
            "only_esp", "esp_memory_limit_bytes", "esp_chunk_policy", "esp_safety_factor",
            "equivalence_groups", "linear_constraints", "source", "fit_input_hash",
        },
        path="resp_fit_input",
    )
    if not isinstance(data["radius_overrides"], Mapping):
        raise ValidationError("invalid_wire_type", "expected object", "resp_fit_input.radius_overrides")
    if not isinstance(data["equivalence_groups"], list):
        raise ValidationError("invalid_wire_type", "expected array", "resp_fit_input.equivalence_groups")
    if not isinstance(data["linear_constraints"], list):
        raise ValidationError("invalid_wire_type", "expected array", "resp_fit_input.linear_constraints")
    groups = []
    for index, item in enumerate(data["equivalence_groups"]):
        path = f"resp_fit_input.equivalence_groups[{index}]"
        group = _strict_object(item, required={"model_id", "atom_ids", "source"}, path=path)
        if not isinstance(group["atom_ids"], list) or not all(isinstance(atom_id, str) for atom_id in group["atom_ids"]):
            raise ValidationError("invalid_wire_type", "expected string array", f"{path}.atom_ids")
        groups.append(RespEquivalenceGroup(group["model_id"], tuple(group["atom_ids"]), group["source"]))
    constraints = []
    for index, item in enumerate(data["linear_constraints"]):
        path = f"resp_fit_input.linear_constraints[{index}]"
        constraint = _strict_object(
            item,
            required={
                "model_id", "constraint_id", "role", "atom_ids",
                "coefficients", "target_charge", "source",
            },
            path=path,
        )
        if (
            not isinstance(constraint["atom_ids"], list)
            or not all(isinstance(atom_id, str) for atom_id in constraint["atom_ids"])
            or not isinstance(constraint["coefficients"], list)
        ):
            raise ValidationError("invalid_wire_type", "invalid linear constraint arrays", path)
        constraints.append(RespLinearConstraint(
            constraint["model_id"],
            constraint["constraint_id"],
            constraint["role"],
            tuple(constraint["atom_ids"]),
            tuple(constraint["coefficients"]),
            constraint["target_charge"],
            constraint["source"],
        ))
    result = RespFitInput(
        schema_version=data["schema_version"],
        backend=data["backend"],
        basis_family=data["basis_family"],
        metal_basis_policy=data["metal_basis_policy"],
        basis_source=data["basis_source"],
        optimize_geometry=data["optimize_geometry"],
        grid_density=data["grid_density"],
        grid_cell_layer=data["grid_cell_layer"],
        radius_overrides=data["radius_overrides"],
        restraint_a1=data["restraint_a1"],
        restraint_a2=data["restraint_a2"],
        two_stage=data["two_stage"],
        only_esp=data["only_esp"],
        esp_memory_limit_bytes=data["esp_memory_limit_bytes"],
        esp_chunk_policy=data["esp_chunk_policy"],
        esp_safety_factor=data["esp_safety_factor"],
        equivalence_groups=tuple(groups),
        linear_constraints=tuple(constraints),
        source=data["source"],
        fit_input_hash=data["fit_input_hash"],
    )
    validate_resp_fit_input(result)
    return result


def resp_fit_input_to_dict(value: RespFitInput) -> dict[str, Any]:
    validate_resp_fit_input(value)
    return _canonicalize(value)


def resp_fit_input_dumps(value: RespFitInput) -> str:
    return json.dumps(resp_fit_input_to_dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def resp_fit_input_loads(payload: str) -> RespFitInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc), "resp_fit_input") from exc
    return resp_fit_input_from_dict(value)


def load_resp_fit_input(path: str | Path) -> RespFitInput:
    source_path = Path(path).resolve()
    try:
        return resp_fit_input_loads(source_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValidationError("resp_fit_input_read_failed", str(exc), str(source_path)) from exc


__all__ = [
    "RESP_FIT_INPUT_SCHEMA_VERSION", "SUPPORTED_ESP_CHUNK_POLICIES", "SUPPORTED_METAL_BASIS_POLICIES",
    "SUPPORTED_RESP_BACKENDS",
    "RespEquivalenceGroup", "RespLinearConstraint", "RespFitInput", "load_resp_fit_input", "resp_fit_input_dumps",
    "resp_fit_input_from_dict", "resp_fit_input_loads", "resp_fit_input_to_dict", "validate_resp_fit_input",
]
