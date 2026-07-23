"""Hash-closed file input for transactional parameter application."""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
from pathlib import Path
from typing import Any

from .assigned_system import ForceRealizationProtocol, _parse_force_protocol, validate_force_realization_protocol
from .contracts import (
    MetalAssignmentInput,
    ParameterizationResult,
    ValidationError,
    _canonicalize,
    _sha256,
    _strict_object,
    parameterization_result_from_dict,
    validate_result,
)


APPLY_INPUT_SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class ParameterApplyInput:
    schema_version: int
    parameterization_result: ParameterizationResult
    force_protocol: ForceRealizationProtocol
    apply_input_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "parameterization_result": _canonicalize(self.parameterization_result),
            "force_protocol": _canonicalize(self.force_protocol),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ParameterApplyInput":
        return replace(self, apply_input_hash=self.computed_hash())


def validate_parameter_apply_input(request: MetalAssignmentInput, value: ParameterApplyInput) -> None:
    if value.schema_version != APPLY_INPUT_SCHEMA_VERSION:
        raise ValidationError("unsupported_parameter_apply_input_version", str(value.schema_version))
    validate_result(request, value.parameterization_result)
    validate_force_realization_protocol(value.force_protocol)
    if not value.apply_input_hash or value.apply_input_hash != value.computed_hash():
        raise ValidationError("stale_parameter_apply_input_hash", request.request_id)


def parameter_apply_input_to_dict(
    request: MetalAssignmentInput,
    value: ParameterApplyInput,
) -> dict[str, Any]:
    validate_parameter_apply_input(request, value)
    return _canonicalize(value)


def parameter_apply_input_from_dict(
    request: MetalAssignmentInput,
    value: Any,
) -> ParameterApplyInput:
    data = _strict_object(
        value,
        required={"schema_version", "parameterization_result", "force_protocol", "apply_input_hash"},
        path="parameter_apply_input",
    )
    parsed = ParameterApplyInput(
        schema_version=data["schema_version"],
        parameterization_result=parameterization_result_from_dict(data["parameterization_result"], request),
        force_protocol=_parse_force_protocol(data["force_protocol"]),
        apply_input_hash=data["apply_input_hash"],
    )
    validate_parameter_apply_input(request, parsed)
    return parsed


def parameter_apply_input_dumps(request: MetalAssignmentInput, value: ParameterApplyInput) -> str:
    return json.dumps(
        parameter_apply_input_to_dict(request, value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def parameter_apply_input_loads(request: MetalAssignmentInput, payload: str) -> ParameterApplyInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return parameter_apply_input_from_dict(request, value)


def load_parameter_apply_input(request: MetalAssignmentInput, path: str | Path) -> ParameterApplyInput:
    source = Path(path).resolve()
    try:
        payload = source.read_text(encoding="utf-8")
    except OSError as exc:
        raise ValidationError("parameter_apply_input_read_failed", str(exc), str(source)) from exc
    return parameter_apply_input_loads(request, payload)


__all__ = [
    "APPLY_INPUT_SCHEMA_VERSION", "ParameterApplyInput", "load_parameter_apply_input",
    "parameter_apply_input_dumps", "parameter_apply_input_from_dict", "parameter_apply_input_loads",
    "parameter_apply_input_to_dict", "validate_parameter_apply_input",
]
