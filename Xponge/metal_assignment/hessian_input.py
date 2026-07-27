"""Hash-closed protocol for locked small-model Hessian calculations."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
import math
from pathlib import Path
from typing import Any

from .contracts import ValidationError, _canonicalize, _sha256, _strict_object


HESSIAN_FIT_INPUT_SCHEMA_VERSION = 2
SUPPORTED_HESSIAN_BACKENDS = frozenset({"pyscf"})
SUPPORTED_HESSIAN_METHODS = frozenset({"scf"})
SUPPORTED_HESSIAN_METAL_BASIS_POLICIES = frozenset({"require_ecp"})
SUPPORTED_HESSIAN_SCF_REFERENCES = frozenset({"auto", "rhf", "rohf", "uhf"})


@dataclass(frozen=True, slots=True)
class HessianFitInput:
    """Versioned scientific and resource protocol for analytical Hessians."""

    schema_version: int
    backend: str
    method: str
    basis_family: str
    metal_basis_policy: str
    basis_source: str
    optimize_geometry: bool
    scf_reference: str
    scf_convergence_tolerance: float
    scf_max_cycles: int
    threads: int
    memory_limit_bytes: int
    source: str
    fit_input_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "fit_input_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "HessianFitInput":
        return replace(self, fit_input_hash=self.computed_hash())


def validate_hessian_fit_input(value: HessianFitInput) -> None:
    path = "hessian_fit_input"
    if value.schema_version != HESSIAN_FIT_INPUT_SCHEMA_VERSION:
        raise ValidationError(
            "unsupported_hessian_fit_input_version", str(value.schema_version), path
        )
    if value.backend not in SUPPORTED_HESSIAN_BACKENDS:
        raise ValidationError("unsupported_hessian_backend", value.backend, f"{path}.backend")
    if value.method not in SUPPORTED_HESSIAN_METHODS:
        raise ValidationError("unsupported_hessian_method", value.method, f"{path}.method")
    if not isinstance(value.basis_family, str) or not value.basis_family.strip():
        raise ValidationError("invalid_hessian_basis_family", str(value.basis_family), f"{path}.basis_family")
    if value.metal_basis_policy not in SUPPORTED_HESSIAN_METAL_BASIS_POLICIES:
        raise ValidationError(
            "unsupported_hessian_metal_basis_policy",
            str(value.metal_basis_policy),
            f"{path}.metal_basis_policy",
        )
    if not isinstance(value.basis_source, str) or not value.basis_source.strip():
        raise ValidationError("missing_hessian_basis_source", str(value.basis_source), f"{path}.basis_source")
    if value.optimize_geometry is not False:
        raise ValidationError(
            "hessian_geometry_must_remain_locked",
            "derived small-model coordinates cannot be optimized inside the Hessian provider",
            f"{path}.optimize_geometry",
        )
    if value.scf_reference not in SUPPORTED_HESSIAN_SCF_REFERENCES:
        raise ValidationError(
            "unsupported_hessian_scf_reference",
            str(value.scf_reference),
            f"{path}.scf_reference",
        )
    tolerance = value.scf_convergence_tolerance
    if (
        isinstance(tolerance, bool)
        or not isinstance(tolerance, (int, float))
        or not math.isfinite(tolerance)
        or tolerance <= 0
    ):
        raise ValidationError("invalid_hessian_scf_tolerance", str(tolerance), f"{path}.scf_convergence_tolerance")
    for name in ("scf_max_cycles", "threads", "memory_limit_bytes"):
        number = getattr(value, name)
        if isinstance(number, bool) or not isinstance(number, int) or number <= 0:
            raise ValidationError(f"invalid_hessian_{name}", str(number), f"{path}.{name}")
    if not isinstance(value.source, str) or not value.source.strip():
        raise ValidationError("missing_hessian_protocol_source", str(value.source), f"{path}.source")
    if not value.fit_input_hash or value.fit_input_hash != value.computed_hash():
        raise ValidationError("stale_hessian_fit_input_hash", "Hessian fit input hash mismatch", path)


def hessian_fit_input_from_dict(value: Any) -> HessianFitInput:
    data = _strict_object(
        value,
        required={
            "schema_version", "backend", "method", "basis_family", "metal_basis_policy",
            "basis_source", "optimize_geometry", "scf_reference", "scf_convergence_tolerance",
            "scf_max_cycles", "threads", "memory_limit_bytes", "source", "fit_input_hash",
        },
        path="hessian_fit_input",
    )
    result = HessianFitInput(**data)
    validate_hessian_fit_input(result)
    return result


def hessian_fit_input_to_dict(value: HessianFitInput) -> dict[str, Any]:
    validate_hessian_fit_input(value)
    return _canonicalize(value)


def hessian_fit_input_dumps(value: HessianFitInput) -> str:
    return json.dumps(
        hessian_fit_input_to_dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def hessian_fit_input_loads(payload: str) -> HessianFitInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc), "hessian_fit_input") from exc
    return hessian_fit_input_from_dict(value)


def load_hessian_fit_input(path: str | Path) -> HessianFitInput:
    source_path = Path(path).resolve()
    try:
        return hessian_fit_input_loads(source_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValidationError("hessian_fit_input_read_failed", str(exc), str(source_path)) from exc


__all__ = [
    "HESSIAN_FIT_INPUT_SCHEMA_VERSION", "SUPPORTED_HESSIAN_BACKENDS",
    "SUPPORTED_HESSIAN_METHODS", "SUPPORTED_HESSIAN_METAL_BASIS_POLICIES",
    "SUPPORTED_HESSIAN_SCF_REFERENCES",
    "HessianFitInput", "hessian_fit_input_dumps", "hessian_fit_input_from_dict",
    "hessian_fit_input_loads", "hessian_fit_input_to_dict", "load_hessian_fit_input",
    "validate_hessian_fit_input",
]
