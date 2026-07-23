"""Hash-closed charge protocol for non-metal base-force-field components."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
from pathlib import Path
from typing import Any

from .contracts import ValidationError, _canonicalize, _sha256, _strict_object


BASE_CHARGE_INPUT_SCHEMA_VERSION = 1
SUPPORTED_BASE_CHARGE_METHODS = frozenset({"gasteiger", "tpacm4"})


@dataclass(frozen=True, slots=True)
class BaseChargeInput:
    """Explicit empirical charge method for complete projected components.

    Absence of this object means that GAFF/GAFF2 contributes no partial
    charges.  There is deliberately no implicit method or residue-name
    fallback; site RESP or external charge artifacts may override these
    component charges later through the normal overlay precedence rules.
    """

    schema_version: int
    method: str
    source: str
    charge_input_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "charge_input_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "BaseChargeInput":
        return replace(self, charge_input_hash=self.computed_hash())


def validate_base_charge_input(value: BaseChargeInput) -> None:
    if value.schema_version != BASE_CHARGE_INPUT_SCHEMA_VERSION:
        raise ValidationError(
            "unsupported_base_charge_input_version",
            str(value.schema_version),
            "base_charge_input.schema_version",
        )
    if value.method not in SUPPORTED_BASE_CHARGE_METHODS:
        raise ValidationError(
            "unsupported_base_charge_method",
            str(value.method),
            "base_charge_input.method",
        )
    if not isinstance(value.source, str) or not value.source.strip():
        raise ValidationError(
            "missing_base_charge_source",
            "source is required",
            "base_charge_input.source",
        )
    if not value.charge_input_hash or value.charge_input_hash != value.computed_hash():
        raise ValidationError(
            "stale_base_charge_input_hash",
            "base charge input hash mismatch",
            "base_charge_input.charge_input_hash",
        )


def base_charge_input_from_dict(value: Any) -> BaseChargeInput:
    data = _strict_object(
        value,
        required={"schema_version", "method", "source", "charge_input_hash"},
        path="base_charge_input",
    )
    result = BaseChargeInput(
        schema_version=data["schema_version"],
        method=str(data["method"]).strip().lower(),
        source=data["source"],
        charge_input_hash=data["charge_input_hash"],
    )
    validate_base_charge_input(result)
    return result


def base_charge_input_to_dict(value: BaseChargeInput) -> dict[str, Any]:
    validate_base_charge_input(value)
    return _canonicalize(value)


def base_charge_input_dumps(value: BaseChargeInput) -> str:
    return json.dumps(
        base_charge_input_to_dict(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def base_charge_input_loads(payload: str) -> BaseChargeInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc), "base_charge_input") from exc
    return base_charge_input_from_dict(value)


def load_base_charge_input(path: str | Path) -> BaseChargeInput:
    source_path = Path(path).resolve()
    try:
        return base_charge_input_loads(source_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValidationError("base_charge_input_read_failed", str(exc), str(source_path)) from exc


__all__ = [
    "BASE_CHARGE_INPUT_SCHEMA_VERSION", "SUPPORTED_BASE_CHARGE_METHODS", "BaseChargeInput",
    "base_charge_input_dumps", "base_charge_input_from_dict", "base_charge_input_loads",
    "base_charge_input_to_dict", "load_base_charge_input", "validate_base_charge_input",
]
