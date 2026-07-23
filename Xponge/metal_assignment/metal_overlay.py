"""Topology-bound metal overlay builders and native 12-6 ion provider."""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
import subprocess
from typing import Any, Mapping

from .contracts import (
    MetalAssignmentInput,
    MetalParameterOverlay,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _sha256,
    _validate_atom_parameter_maps,
    _validate_bonded_parameter_mapping,
    _validate_mass_parameters,
)
from .input import MetalAssignmentPackage, validate_package
from ._worker_runtime import worker_command


ION_WORKER_PROTOCOL_VERSION = 1
SUPPORTED_ION_WATER_MODELS = frozenset({"tip3p", "tip4pew", "spce", "opc"})


@dataclass(frozen=True, slots=True)
class MetalOverlayBuildReport:
    schema_version: int
    topology_hash: str
    projection_hash: str
    interaction_model: str
    covered_atom_ids: tuple[str, ...]
    parameter_source: str
    precedence: int
    atom_type_count: int
    charge_count: int
    mass_count: int
    lj_count: int
    bonded_term_count: int
    provenance: Mapping[str, Any]
    report_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "provenance")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "topology_hash": self.topology_hash,
            "projection_hash": self.projection_hash,
            "interaction_model": self.interaction_model,
            "covered_atom_ids": _canonicalize(self.covered_atom_ids),
            "parameter_source": self.parameter_source,
            "precedence": self.precedence,
            "atom_type_count": self.atom_type_count,
            "charge_count": self.charge_count,
            "mass_count": self.mass_count,
            "lj_count": self.lj_count,
            "bonded_term_count": self.bonded_term_count,
            "provenance": _canonicalize(self.provenance),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "MetalOverlayBuildReport":
        return replace(self, report_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class MetalOverlayBuildOutput:
    overlay: MetalParameterOverlay
    report: MetalOverlayBuildReport


def build_metal_parameter_overlay(
    request: MetalAssignmentInput,
    *,
    covered_atom_ids: tuple[str, ...],
    atom_types: Mapping[str, str],
    charges: Mapping[str, float],
    masses: Mapping[str, float],
    lj_parameters: Mapping[str, Any],
    bonded_parameters: Mapping[str, Any],
    parameter_source: str,
    precedence: int = 100,
    require_metal_charges: bool = False,
    provenance: Mapping[str, Any] | None = None,
) -> MetalOverlayBuildOutput:
    """Build a validated metal layer without changing prepared topology."""

    topology = request.topology
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    coverage = set(covered_atom_ids)
    if len(coverage) != len(covered_atom_ids):
        raise ValidationError("duplicate_overlay_atom", "metal overlay coverage must be unique")
    metal_ids = {atom_id for atom_id, atom in atom_by_id.items() if atom.is_metal}
    allowed_ids = {
        atom_id for atom_id, atom in atom_by_id.items()
        if atom.is_metal or "metalOverlay" in atom.scopes
    }
    if not metal_ids or not metal_ids <= coverage:
        raise ValidationError("metal_overlay_coverage_mismatch", "every metal atom must be covered")
    if not coverage <= allowed_ids:
        raise ValidationError("metal_overlay_scope_mismatch", "coverage exceeds explicit metalOverlay scope")
    if not parameter_source:
        raise ValidationError("missing_parameter_source", "metal overlay source is required")
    if isinstance(precedence, bool) or not isinstance(precedence, int) or precedence <= 0:
        raise ValidationError("invalid_overlay_precedence", str(precedence))
    _validate_atom_parameter_maps(
        atom_types=atom_types,
        charges=charges,
        lj_parameters=lj_parameters,
        covered_atom_ids=coverage,
        path="metal_overlay",
    )
    _validate_mass_parameters(masses, covered_atom_ids=coverage, path="metal_overlay.masses")
    for name, values in (("atom_types", atom_types), ("masses", masses), ("lj_parameters", lj_parameters)):
        if not metal_ids <= set(values):
            raise ValidationError("incomplete_metal_parameter_map", name, "metal_overlay")
    if require_metal_charges and not metal_ids <= set(charges):
        raise ValidationError("incomplete_metal_parameter_map", "charges", "metal_overlay")
    _validate_bonded_parameter_mapping(
        bonded_parameters,
        topology=topology,
        allowed_atom_ids=coverage,
        require_metal=True,
        path="metal_overlay.bonded_parameters",
    )
    overlay = MetalParameterOverlay(
        topology_hash=topology.topology_hash,
        covered_atom_ids=covered_atom_ids,
        atom_types=atom_types,
        charges=charges,
        masses=masses,
        lj_parameters=lj_parameters,
        bonded_parameters=bonded_parameters,
        parameter_source=parameter_source,
        precedence=precedence,
    )
    report = MetalOverlayBuildReport(
        schema_version=request.schema_version,
        topology_hash=topology.topology_hash,
        projection_hash=request.projection_hash,
        interaction_model=request.interaction_model,
        covered_atom_ids=covered_atom_ids,
        parameter_source=parameter_source,
        precedence=precedence,
        atom_type_count=len(atom_types),
        charge_count=len(charges),
        mass_count=len(masses),
        lj_count=len(lj_parameters),
        bonded_term_count=len(bonded_parameters),
        provenance=provenance or {},
    ).with_computed_hash()
    return MetalOverlayBuildOutput(overlay, report)


def _invoke_ion_worker(payload: Mapping[str, Any], *, timeout_seconds: float) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._ion_worker"),
            input=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("ion_assignment_worker_timeout", str(timeout_seconds)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError(
            "invalid_ion_assignment_worker_output",
            completed.stderr.strip()[-1000:] or str(exc),
        ) from exc
    if not isinstance(response, dict) or response.get("ok") is not True or completed.returncode != 0:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "ion worker failed") if isinstance(error, dict) else "ion worker failed"
        if completed.stderr.strip():
            message = f"{message}; stderr={completed.stderr.strip()[-1000:]}"
        raise ValidationError("ion_assignment_worker_failed", message)
    return response


def assign_nonbonded_metal_ions(
    package: MetalAssignmentPackage,
    *,
    water_model: str,
    timeout_seconds: float = 60.0,
) -> MetalOverlayBuildOutput:
    """Resolve water-model-specific 12-6 ion parameters in an isolated process."""

    validate_package(package)
    request = package.request
    if request.interaction_model != "nonbonded_12_6":
        raise ValidationError("ion_provider_interaction_model_mismatch", request.interaction_model)
    normalized_water_model = water_model.lower()
    if normalized_water_model not in SUPPORTED_ION_WATER_MODELS:
        raise ValidationError("unsupported_ion_water_model", water_model)
    metal_atoms = tuple(atom for atom in request.topology.atoms if atom.is_metal)
    if not metal_atoms:
        raise ValidationError("missing_metal_atom", "ion assignment requires at least one metal")
    for atom in metal_atoms:
        if not atom.formal_charge_known or atom.formal_charge is None:
            raise ValidationError("unresolved_metal_formal_charge", atom.external_id)
    payload = {
        "protocol_version": ION_WORKER_PROTOCOL_VERSION,
        "topology_hash": request.topology.topology_hash,
        "water_model": normalized_water_model,
        "metals": [
            {
                "external_id": atom.external_id,
                "element": atom.element,
                "formal_charge": atom.formal_charge,
            }
            for atom in metal_atoms
        ],
    }
    response = _invoke_ion_worker(payload, timeout_seconds=timeout_seconds)
    required = {
        "ok", "protocol_version", "topology_hash", "water_model", "provider_version", "parameters",
        "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_ion_assignment_worker_response", "unexpected fields")
    response_payload = dict(response)
    response_hash = response_payload.pop("response_hash")
    if response_hash != _sha256(response_payload):
        raise ValidationError("stale_ion_assignment_worker_hash", "response hash mismatch")
    if (
        response["protocol_version"] != ION_WORKER_PROTOCOL_VERSION
        or response["topology_hash"] != request.topology.topology_hash
        or response["water_model"] != normalized_water_model
        or not isinstance(response["provider_version"], str)
        or not response["provider_version"]
        or not isinstance(response["parameters"], Mapping)
    ):
        raise ValidationError("ion_assignment_worker_identity_mismatch", "worker response belongs to another request")
    metal_ids = tuple(atom.external_id for atom in metal_atoms)
    if set(response["parameters"]) != set(metal_ids):
        raise ValidationError("ion_assignment_worker_coverage_mismatch", "metal coverage mismatch")
    source = f"xponge:ion-12-6:{normalized_water_model}:{response['provider_version']}"
    parameters = response["parameters"]
    return build_metal_parameter_overlay(
        request,
        covered_atom_ids=metal_ids,
        atom_types={atom_id: parameters[atom_id]["atom_type"] for atom_id in metal_ids},
        charges={atom_id: parameters[atom_id]["charge"] for atom_id in metal_ids},
        masses={atom_id: parameters[atom_id]["mass"] for atom_id in metal_ids},
        lj_parameters={atom_id: parameters[atom_id]["lj"] for atom_id in metal_ids},
        bonded_parameters={},
        parameter_source=source,
        precedence=100,
        require_metal_charges=True,
        provenance={
            "provider_kind": "ion_12_6",
            "water_model": normalized_water_model,
            "provider_version": response["provider_version"],
            "worker_response_hash": response["response_hash"],
            "template_names": {
                atom_id: parameters[atom_id]["template_name"] for atom_id in metal_ids
            },
        },
    )


__all__ = [
    "MetalOverlayBuildOutput", "MetalOverlayBuildReport", "SUPPORTED_ION_WATER_MODELS",
    "assign_nonbonded_metal_ions", "build_metal_parameter_overlay",
]
