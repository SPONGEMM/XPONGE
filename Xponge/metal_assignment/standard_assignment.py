"""Standard biomolecular base-force-field assignment over locked topology."""

from __future__ import annotations

import json
import math
import subprocess
from typing import Any, Mapping

from .base_assignment import BaseAssignmentOutput, BaseAssignmentReport
from .contracts import (
    BaseForceFieldOverlay,
    ValidationError,
    _sha256,
    base_metal_component_supported,
)
from .input import MetalAssignmentPackage, validate_package
from ._worker_runtime import worker_command


WORKER_PROTOCOL_VERSION = 1
SUPPORTED_STANDARD_FORCE_FIELDS = frozenset({"ff14sb", "ff19sb"})
SUPPORTED_STANDARD_WATER_MODELS = frozenset({"tip3p", "spce", "opc"})
SUPPORTED_TEMPLATE_PROVIDERS = {
    "standard_biomolecular": SUPPORTED_STANDARD_FORCE_FIELDS,
    "glycam": frozenset({"glycam_06j"}),
    "lipid": frozenset({"lipid21"}),
}


def _worker_payload(
    package: MetalAssignmentPackage,
    *,
    provider: str,
    force_field: str,
    water_model: str,
) -> dict[str, Any]:
    request = package.request
    topology = request.topology
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    residue_by_atom: dict[str, Any] = {}
    for residue in topology.residues:
        for atom_id in residue.atom_ids:
            residue_by_atom[atom_id] = residue
    artifact_by_id = {
        artifact.external_id: artifact
        for artifact in package.prepared_artifacts.structural_artifacts.assignment_components
    }
    active_bonds = tuple(bond for bond in topology.bonds if bond.active)
    components = tuple(
        component
        for component in request.assignment_components
        if component.base_force_field and component.provider == provider
    )
    if not components:
        raise ValidationError(
            "missing_provider_scope",
            f"no base assignment component uses provider {provider}",
            "assignment_components",
        )
    standard_atom_ids = {
        atom_id
        for component in components
        for atom_id in component.atom_ids
    }
    cross_component_links = [
        {
            "external_id": link.external_id,
            "atom_ids": list(link.atom_ids),
        }
        for link in topology.links
        if link.active and set(link.atom_ids) <= standard_atom_ids
    ]
    payload_components: list[dict[str, Any]] = []
    for component in components:
        path = f"assignment_components.{component.external_id}"
        residues = {residue_by_atom[atom_id].external_id for atom_id in component.atom_ids}
        if len(residues) != 1:
            raise ValidationError("standard_component_crosses_residues", component.external_id, path)
        residue = residue_by_atom[component.atom_ids[0]]
        if set(residue.atom_ids) != set(component.atom_ids):
            raise ValidationError("standard_component_is_partial_residue", component.external_id, path)
        if (
            any(atom_by_id[atom_id].is_metal for atom_id in component.atom_ids)
            and not base_metal_component_supported(component)
        ):
            raise ValidationError("metal_in_base_force_field", component.external_id, path)
        artifact = artifact_by_id.get(component.external_id)
        if artifact is None:
            raise ValidationError("missing_structural_artifact", component.external_id, path)
        component_atom_ids = set(component.atom_ids)
        bonds = [bond for bond in active_bonds if set(bond.atom_ids) <= component_atom_ids]
        if {bond.external_id for bond in bonds} != set(artifact.active_edge_ids):
            raise ValidationError("base_assignment_edge_mismatch", component.external_id, path)
        payload_components.append({
            "external_id": component.external_id,
            "proof_hash": component.chemical_topology_proof.proof_hash,
            "artifact_hash": artifact.artifact_hash,
            "residue_name": residue.residue_name,
            "polymer_position": residue.polymer_position,
            "template_id": residue.template_id,
            "atoms": [
                {
                    "external_id": atom_id,
                    "atom_name": atom_by_id[atom_id].atom_name,
                    "element": atom_by_id[atom_id].element,
                }
                for atom_id in component.atom_ids
            ],
            "bonds": [
                {"external_id": bond.external_id, "atom_ids": list(bond.atom_ids)}
                for bond in bonds
            ],
        })
    return {
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": provider,
        "force_field": force_field,
        "water_model": water_model,
        "topology_hash": topology.topology_hash,
        "components": payload_components,
        "links": cross_component_links,
    }


def _invoke_worker(payload: Mapping[str, Any], *, timeout_seconds: float) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._biomolecular_worker"),
            input=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("standard_assignment_worker_timeout", str(timeout_seconds)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        stderr = completed.stderr.strip()
        raise ValidationError(
            "invalid_standard_assignment_worker_output",
            stderr[-1000:] if stderr else str(exc),
        ) from exc
    if not isinstance(response, dict) or response.get("ok") is not True:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "worker failed") if isinstance(error, dict) else "worker failed"
        stderr = completed.stderr.strip()
        if stderr:
            message = f"{message}; stderr={stderr[-1000:]}"
        raise ValidationError("standard_assignment_worker_failed", message)
    if completed.returncode != 0:
        raise ValidationError("standard_assignment_worker_failed", f"exit status {completed.returncode}")
    return response


def _validate_response(payload: Mapping[str, Any], response: Mapping[str, Any]) -> None:
    required = {
        "ok", "protocol_version", "provider", "topology_hash", "provider_version", "force_field",
        "water_model", "components", "cross_component_bonded_parameters", "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_standard_assignment_worker_response", "unexpected fields")
    if (
        response["protocol_version"] != WORKER_PROTOCOL_VERSION
        or response["provider"] != payload["provider"]
        or response["topology_hash"] != payload["topology_hash"]
        or response["force_field"] != payload["force_field"]
        or response["water_model"] != payload["water_model"]
        or not isinstance(response["provider_version"], str)
        or not response["provider_version"]
        or not isinstance(response["components"], list)
        or not isinstance(response["cross_component_bonded_parameters"], dict)
    ):
        raise ValidationError("standard_assignment_worker_identity_mismatch", "response belongs to another request")
    response_payload = dict(response)
    response_hash = response_payload.pop("response_hash")
    if response_hash != _sha256(response_payload):
        raise ValidationError("stale_standard_assignment_worker_hash", "response hash mismatch")
    expected_by_id = {component["external_id"]: component for component in payload["components"]}
    component_by_atom = {
        atom["external_id"]: component["external_id"]
        for component in payload["components"]
        for atom in component["atoms"]
    }
    actual_ids = [component.get("external_id") for component in response["components"] if isinstance(component, dict)]
    if len(actual_ids) != len(response["components"]) or set(actual_ids) != set(expected_by_id):
        raise ValidationError("standard_assignment_worker_component_mismatch", "component coverage mismatch")
    for component in response["components"]:
        required_component = {
            "external_id", "proof_hash", "artifact_hash", "template_name", "atom_types", "charges",
            "masses", "lj_parameters", "bonded_parameters", "ignored_force_kinds", "output_hash",
        }
        if set(component) != required_component:
            raise ValidationError("invalid_standard_assignment_worker_component", str(component.get("external_id")))
        expected = expected_by_id[component["external_id"]]
        if component["proof_hash"] != expected["proof_hash"] or component["artifact_hash"] != expected["artifact_hash"]:
            raise ValidationError("standard_assignment_worker_component_identity_mismatch", component["external_id"])
        component_payload = dict(component)
        output_hash = component_payload.pop("output_hash")
        if output_hash != _sha256(component_payload):
            raise ValidationError("stale_standard_assignment_worker_component_hash", component["external_id"])
        expected_atom_ids = {atom["external_id"] for atom in expected["atoms"]}
        for name in ("atom_types", "charges", "masses", "lj_parameters"):
            if set(component[name]) != expected_atom_ids:
                raise ValidationError("standard_assignment_worker_atom_coverage_mismatch", component["external_id"])
    for term_id, raw_term in response["cross_component_bonded_parameters"].items():
        if not isinstance(term_id, str) or not term_id or not isinstance(raw_term, Mapping):
            raise ValidationError("invalid_cross_component_bonded_term", str(term_id))
        if set(raw_term) != {"kind", "atom_ids", "parameters", "source"}:
            raise ValidationError("invalid_cross_component_bonded_term", term_id)
        atom_ids = raw_term["atom_ids"]
        expected_atom_count = {
            "bond": 2,
            "angle": 3,
            "proper_dihedral": 4,
            "improper_dihedral": 4,
        }.get(raw_term["kind"])
        if (
            expected_atom_count is None
            or not isinstance(atom_ids, list)
            or len(atom_ids) != expected_atom_count
            or len(set(atom_ids)) != len(atom_ids)
            or any(atom_id not in component_by_atom for atom_id in atom_ids)
            or len({component_by_atom[atom_id] for atom_id in atom_ids}) < 2
            or not isinstance(raw_term["parameters"], Mapping)
            or not isinstance(raw_term["source"], str)
            or not raw_term["source"]
        ):
            raise ValidationError("invalid_cross_component_bonded_term", term_id)


def assign_template_force_field(
    package: MetalAssignmentPackage,
    *,
    provider: str,
    force_field: str = "ff14sb",
    water_model: str = "tip3p",
    timeout_seconds: float = 120.0,
) -> BaseAssignmentOutput:
    """Assign exact standard templates without mutating the prepared topology."""

    validate_package(package)
    provider = provider.lower()
    force_field = force_field.lower()
    water_model = water_model.lower()
    supported_force_fields = SUPPORTED_TEMPLATE_PROVIDERS.get(provider)
    if supported_force_fields is None:
        raise ValidationError("unsupported_template_provider", provider)
    if force_field not in supported_force_fields:
        raise ValidationError("unsupported_template_force_field", f"{provider}:{force_field}")
    if water_model not in SUPPORTED_STANDARD_WATER_MODELS:
        raise ValidationError("unsupported_standard_water_model", water_model)
    if provider not in package.request.capability_snapshot.base_force_field_providers:
        raise ValidationError("provider_not_in_capability_snapshot", provider)
    payload = _worker_payload(
        package,
        provider=provider,
        force_field=force_field,
        water_model=water_model,
    )
    response = _invoke_worker(payload, timeout_seconds=timeout_seconds)
    _validate_response(payload, response)

    atom_types: dict[str, str] = {}
    charges: dict[str, float] = {}
    masses: dict[str, float] = {}
    lj_parameters: dict[str, Any] = {}
    bonded_parameters: dict[str, Any] = {}
    for component in response["components"]:
        for target, source, code in (
            (atom_types, component["atom_types"], "duplicate_base_atom_type"),
            (charges, component["charges"], "duplicate_base_charge"),
            (masses, component["masses"], "duplicate_base_mass"),
            (lj_parameters, component["lj_parameters"], "duplicate_base_lj_parameter"),
            (bonded_parameters, component["bonded_parameters"], "duplicate_base_bonded_term"),
        ):
            overlap = set(target) & set(source)
            if overlap:
                raise ValidationError(code, ",".join(sorted(overlap)))
            target.update(source)
    overlap = set(bonded_parameters) & set(response["cross_component_bonded_parameters"])
    if overlap:
        raise ValidationError("duplicate_base_bonded_term", ",".join(sorted(overlap)))
    bonded_parameters.update(response["cross_component_bonded_parameters"])
    covered_atom_ids = tuple(
        atom["external_id"]
        for component in payload["components"]
        for atom in component["atoms"]
    )
    parameter_source = (
        f"xponge:{provider}:{force_field}:{water_model}:{response['provider_version']}"
    )
    overlay = BaseForceFieldOverlay(
        topology_hash=package.request.topology.topology_hash,
        covered_atom_ids=covered_atom_ids,
        atom_types=atom_types,
        charges=charges,
        masses=masses,
        lj_parameters=lj_parameters,
        bonded_parameters=bonded_parameters,
        parameter_source=parameter_source,
    )
    components = response["components"]
    report = BaseAssignmentReport(
        schema_version=package.request.schema_version,
        topology_hash=package.request.topology.topology_hash,
        projection_hash=package.request.projection_hash,
        provider=provider,
        provider_version=response["provider_version"],
        component_ids=tuple(component["external_id"] for component in components),
        proof_hashes={component["external_id"]: component["proof_hash"] for component in components},
        artifact_hashes={component["external_id"]: component["artifact_hash"] for component in components},
        component_output_hashes={component["external_id"]: component["output_hash"] for component in components},
        frcmod_hashes={},
        ignored_force_kinds={
            component["external_id"]: tuple(component["ignored_force_kinds"])
            for component in components
        },
        atom_count=len(covered_atom_ids),
        bonded_term_count=len(bonded_parameters),
        parameter_source=parameter_source,
    ).with_computed_hash()
    return BaseAssignmentOutput(overlay, report)


def assign_standard_biomolecular(
    package: MetalAssignmentPackage,
    *,
    force_field: str = "ff14sb",
    water_model: str = "tip3p",
    timeout_seconds: float = 120.0,
) -> BaseAssignmentOutput:
    return assign_template_force_field(
        package,
        provider="standard_biomolecular",
        force_field=force_field,
        water_model=water_model,
        timeout_seconds=timeout_seconds,
    )


__all__ = [
    "SUPPORTED_STANDARD_FORCE_FIELDS", "SUPPORTED_STANDARD_WATER_MODELS",
    "SUPPORTED_TEMPLATE_PROVIDERS", "assign_standard_biomolecular", "assign_template_force_field",
]
