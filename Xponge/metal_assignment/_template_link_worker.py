"""Isolated bonded-term completion across locked template-provider scopes."""

from __future__ import annotations

from contextlib import redirect_stdout
import importlib
import json
import sys
from typing import Any, Mapping

from ._biomolecular_worker import (
    WORKER_PROTOCOL_VERSION,
    WorkerInputError,
    _cross_component_terms,
    _validate_component,
)
from ._gaff_worker import _canonical_hash


def _load_force_field(name: str) -> None:
    if name == "glycam_06j":
        importlib.import_module("Xponge.forcefield.amber.ff19sb")
        for submodule in (
            "d_pyranose", "l_pyranose", "d_furanose", "l_furanose", "glycoprotein",
        ):
            importlib.import_module(f"Xponge.forcefield.amber.glycam_06j.{submodule}")
        return
    importlib.import_module(f"Xponge.forcefield.amber.{name}")


def _validate_request(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {
        "protocol_version", "provider", "topology_hash", "force_fields",
        "components", "component_providers", "links",
    }:
        raise WorkerInputError("request: invalid template-link payload")
    if value["protocol_version"] != WORKER_PROTOCOL_VERSION or value["provider"] != "template_cross_provider":
        raise WorkerInputError("request: unsupported protocol or provider")
    if not isinstance(value["topology_hash"], str) or len(value["topology_hash"]) != 64:
        raise WorkerInputError("request.topology_hash: expected sha256")
    if not isinstance(value["force_fields"], list) or not value["force_fields"]:
        raise WorkerInputError("request.force_fields: expected non-empty array")
    force_fields = []
    for force_field in value["force_fields"]:
        if not isinstance(force_field, str) or not force_field:
            raise WorkerInputError("request.force_fields: expected strings")
        force_fields.append(force_field)
    if not isinstance(value["components"], list) or not value["components"]:
        raise WorkerInputError("request.components: expected non-empty array")
    components = [_validate_component(component, index) for index, component in enumerate(value["components"])]
    component_ids = {component["external_id"] for component in components}
    providers = value["component_providers"]
    if (
        not isinstance(providers, Mapping)
        or set(providers) != component_ids
        or not all(isinstance(provider, str) and provider for provider in providers.values())
    ):
        raise WorkerInputError("request.component_providers: coverage mismatch")
    atom_component = {
        atom["external_id"]: component["external_id"]
        for component in components
        for atom in component["atoms"]
    }
    if not isinstance(value["links"], list) or not value["links"]:
        raise WorkerInputError("request.links: expected non-empty array")
    links = []
    seen = set()
    for index, link in enumerate(value["links"]):
        if not isinstance(link, Mapping) or set(link) != {"external_id", "atom_ids"}:
            raise WorkerInputError(f"request.links[{index}]: invalid object")
        atom_ids = link["atom_ids"]
        if (
            not isinstance(link["external_id"], str)
            or not link["external_id"]
            or link["external_id"] in seen
            or not isinstance(atom_ids, list)
            or len(atom_ids) != 2
            or atom_ids[0] == atom_ids[1]
            or any(atom_id not in atom_component for atom_id in atom_ids)
            or providers[atom_component[atom_ids[0]]] == providers[atom_component[atom_ids[1]]]
        ):
            raise WorkerInputError(f"request.links[{index}]: expected a cross-provider link")
        seen.add(link["external_id"])
        links.append({"external_id": link["external_id"], "atom_ids": list(atom_ids)})
    return {
        **value,
        "force_fields": force_fields,
        "components": components,
        "component_providers": dict(providers),
        "links": links,
    }


def _execute(value: Any) -> dict[str, Any]:
    request = _validate_request(value)
    for force_field in request["force_fields"]:
        _load_force_field(force_field)
    source = "xponge:template_cross_provider:" + "+".join(request["force_fields"])
    placeholder_outputs = [
        {"atom_types": {atom["external_id"]: "" for atom in component["atoms"]}}
        for component in request["components"]
    ]
    terms = _cross_component_terms(request, placeholder_outputs, source)
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": request["provider"],
        "topology_hash": request["topology_hash"],
        "component_ids": [component["external_id"] for component in request["components"]],
        "bonded_parameters": terms,
        "parameter_source": source,
    }
    response["response_hash"] = _canonical_hash(response)
    return response


def main() -> int:
    try:
        request = json.load(sys.stdin)
        with redirect_stdout(sys.stderr):
            response = _execute(request)
        json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 0
    except Exception as exc:
        json.dump(
            {"ok": False, "error": {"type": type(exc).__name__, "message": str(exc)}},
            sys.stdout,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
