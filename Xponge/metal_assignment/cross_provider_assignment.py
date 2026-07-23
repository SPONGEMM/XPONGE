"""Bonded completion for explicit links crossing template-provider scopes."""

from __future__ import annotations

import json
import math
import subprocess
from typing import Any, Mapping

from .base_assignment import BaseAssignmentOutput, BaseAssignmentReport
from .contracts import BaseForceFieldOverlay, ValidationError, _sha256
from .input import MetalAssignmentPackage, validate_package
from ._worker_runtime import worker_command
from .standard_assignment import WORKER_PROTOCOL_VERSION, _worker_payload


TEMPLATE_PROVIDERS = frozenset({"standard_biomolecular", "glycam", "lipid"})


def _build_payload(
    package: MetalAssignmentPackage,
    *,
    provider_force_fields: Mapping[str, str],
    water_model: str,
) -> dict[str, Any] | None:
    component_by_atom = {
        atom_id: component
        for component in package.request.assignment_components
        if component.base_force_field and component.provider in TEMPLATE_PROVIDERS
        for atom_id in component.atom_ids
    }
    cross_links = [
        link
        for link in package.request.topology.links
        if link.active
        and all(atom_id in component_by_atom for atom_id in link.atom_ids)
        and component_by_atom[link.atom_ids[0]].provider != component_by_atom[link.atom_ids[1]].provider
    ]
    if not cross_links:
        return None
    included_component_ids = {
        component_by_atom[atom_id].external_id
        for link in cross_links
        for atom_id in link.atom_ids
    }
    components = []
    component_providers: dict[str, str] = {}
    force_fields = []
    for provider in sorted({
        component.provider
        for component in component_by_atom.values()
        if component.external_id in included_component_ids
    }):
        force_field = str(provider_force_fields.get(provider) or "").lower()
        if not force_field:
            raise ValidationError("missing_template_provider_force_field", provider)
        provider_payload = _worker_payload(
            package,
            provider=provider,
            force_field=force_field,
            water_model=water_model,
        )
        force_fields.append(force_field)
        for component in provider_payload["components"]:
            if component["external_id"] in included_component_ids:
                components.append(component)
                component_providers[component["external_id"]] = provider
    return {
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": "template_cross_provider",
        "topology_hash": package.request.topology.topology_hash,
        "force_fields": list(dict.fromkeys(force_fields)),
        "components": components,
        "component_providers": component_providers,
        "links": [
            {"external_id": link.external_id, "atom_ids": list(link.atom_ids)}
            for link in cross_links
        ],
    }


def _invoke(payload: Mapping[str, Any], *, timeout_seconds: float) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        completed = subprocess.run(
            worker_command("Xponge.metal_assignment._template_link_worker"),
            input=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("template_link_worker_timeout", str(timeout_seconds)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError(
            "invalid_template_link_worker_output",
            completed.stderr.strip()[-1000:] or str(exc),
        ) from exc
    if completed.returncode != 0 or not isinstance(response, dict) or response.get("ok") is not True:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "worker failed") if isinstance(error, dict) else "worker failed"
        raise ValidationError("template_link_worker_failed", message)
    return response


def _validate_response(payload: Mapping[str, Any], response: Mapping[str, Any]) -> None:
    required = {
        "ok", "protocol_version", "provider", "topology_hash", "component_ids",
        "bonded_parameters", "parameter_source", "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_template_link_worker_response", "unexpected fields")
    if (
        response["protocol_version"] != WORKER_PROTOCOL_VERSION
        or response["provider"] != payload["provider"]
        or response["topology_hash"] != payload["topology_hash"]
        or response["component_ids"] != [item["external_id"] for item in payload["components"]]
        or not isinstance(response["bonded_parameters"], dict)
        or not isinstance(response["parameter_source"], str)
        or not response["parameter_source"]
    ):
        raise ValidationError("template_link_worker_identity_mismatch", "response belongs to another request")
    unhashed = dict(response)
    response_hash = unhashed.pop("response_hash")
    if response_hash != _sha256(unhashed):
        raise ValidationError("stale_template_link_worker_hash", "response hash mismatch")
    known_atoms = {
        atom["external_id"]
        for component in payload["components"]
        for atom in component["atoms"]
    }
    for term_id, term in response["bonded_parameters"].items():
        if not isinstance(term_id, str) or not term_id or not isinstance(term, Mapping):
            raise ValidationError("invalid_template_cross_provider_term", str(term_id))
        if set(term) != {"kind", "atom_ids", "parameters", "source"}:
            raise ValidationError("invalid_template_cross_provider_term", term_id)
        expected_count = {
            "bond": 2,
            "angle": 3,
            "proper_dihedral": 4,
            "improper_dihedral": 4,
        }.get(term["kind"])
        atom_ids = term["atom_ids"]
        if (
            expected_count is None
            or not isinstance(atom_ids, list)
            or len(atom_ids) != expected_count
            or len(set(atom_ids)) != len(atom_ids)
            or any(atom_id not in known_atoms for atom_id in atom_ids)
            or not isinstance(term["parameters"], Mapping)
            or not isinstance(term["source"], str)
            or not term["source"]
        ):
            raise ValidationError("invalid_template_cross_provider_term", term_id)


def assign_template_cross_provider_links(
    package: MetalAssignmentPackage,
    *,
    provider_force_fields: Mapping[str, str],
    water_model: str = "tip3p",
    timeout_seconds: float = 120.0,
) -> BaseAssignmentOutput | None:
    validate_package(package)
    payload = _build_payload(
        package,
        provider_force_fields=provider_force_fields,
        water_model=water_model,
    )
    if payload is None:
        return None
    response = _invoke(payload, timeout_seconds=timeout_seconds)
    _validate_response(payload, response)
    overlay = BaseForceFieldOverlay(
        topology_hash=package.request.topology.topology_hash,
        covered_atom_ids=(),
        atom_types={},
        charges={},
        masses={},
        lj_parameters={},
        bonded_parameters=response["bonded_parameters"],
        parameter_source=response["parameter_source"],
    )
    proof_hashes = {
        component["external_id"]: component["proof_hash"]
        for component in payload["components"]
    }
    artifact_hashes = {
        component["external_id"]: component["artifact_hash"]
        for component in payload["components"]
    }
    report = BaseAssignmentReport(
        schema_version=package.request.schema_version,
        topology_hash=package.request.topology.topology_hash,
        projection_hash=package.request.projection_hash,
        provider="template_cross_provider",
        provider_version="1",
        component_ids=tuple(component["external_id"] for component in payload["components"]),
        proof_hashes=proof_hashes,
        artifact_hashes=artifact_hashes,
        component_output_hashes={},
        frcmod_hashes={},
        ignored_force_kinds={},
        atom_count=0,
        bonded_term_count=len(response["bonded_parameters"]),
        parameter_source=response["parameter_source"],
    ).with_computed_hash()
    return BaseAssignmentOutput(overlay, report)


__all__ = ["TEMPLATE_PROVIDERS", "assign_template_cross_provider_links"]
