"""Isolated RESP worker consuming only a locked derived large-model graph."""

from __future__ import annotations

from contextlib import redirect_stdout
import importlib.metadata
import json
import math
import re
import sys
from typing import Any, Mapping

from ._ion_worker import _canonical_hash, _strict_object


WORKER_PROTOCOL_VERSION = 1


class WorkerInputError(ValueError):
    """Invalid private RESP worker request."""


def _package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "source-tree"


def _finite_number(value: Any, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise WorkerInputError(f"{path}: expected finite number")
    return float(value)


def _validate_fit_protocol(
    value: Any,
    model_id: str,
    model_atom_ids: tuple[str, ...],
    elements: tuple[str, ...],
    metal_elements: tuple[str, ...],
) -> Mapping[str, Any]:
    protocol = _strict_object(
        value,
        {
            "fit_input_hash", "backend", "basis_family", "metal_basis_policy", "basis_source",
            "optimize_geometry", "grid_density",
            "grid_cell_layer", "radius_overrides", "restraint_a1", "restraint_a2", "two_stage",
            "only_esp", "esp_memory_limit_bytes", "esp_chunk_policy", "esp_safety_factor",
            "equivalence_groups", "source",
        },
        "fit_protocol",
    )
    if not isinstance(protocol["fit_input_hash"], str) or len(protocol["fit_input_hash"]) != 64:
        raise WorkerInputError("fit_protocol.fit_input_hash: expected sha256")
    if protocol["backend"] not in {"pyscf", "psi4"}:
        raise WorkerInputError("fit_protocol.backend: expected pyscf or psi4")
    if not isinstance(protocol["basis_family"], str) or not protocol["basis_family"].strip():
        raise WorkerInputError("fit_protocol.basis_family: expected non-empty string")
    if protocol["metal_basis_policy"] != "require_ecp":
        raise WorkerInputError("fit_protocol.metal_basis_policy: expected require_ecp")
    if not isinstance(protocol["basis_source"], str) or not protocol["basis_source"].strip():
        raise WorkerInputError("fit_protocol.basis_source: expected non-empty string")
    if protocol["optimize_geometry"] is not False:
        raise WorkerInputError("fit_protocol.optimize_geometry: locked derived-model coordinates are required")
    if protocol["two_stage"] is not True or protocol["only_esp"] is not False:
        raise WorkerInputError("fit_protocol: restrained two-stage RESP is required")
    _finite_number(protocol["grid_density"], "fit_protocol.grid_density")
    if protocol["grid_density"] <= 0:
        raise WorkerInputError("fit_protocol.grid_density: expected positive value")
    if (
        isinstance(protocol["grid_cell_layer"], bool)
        or not isinstance(protocol["grid_cell_layer"], int)
        or protocol["grid_cell_layer"] <= 0
    ):
        raise WorkerInputError("fit_protocol.grid_cell_layer: expected positive integer")
    for name in ("restraint_a1", "restraint_a2"):
        if _finite_number(protocol[name], f"fit_protocol.{name}") <= 0:
            raise WorkerInputError(f"fit_protocol.{name}: expected positive value")
    if (
        isinstance(protocol["esp_memory_limit_bytes"], bool)
        or not isinstance(protocol["esp_memory_limit_bytes"], int)
        or protocol["esp_memory_limit_bytes"] <= 0
    ):
        raise WorkerInputError("fit_protocol.esp_memory_limit_bytes: expected positive integer")
    if protocol["esp_chunk_policy"] not in {"auto", "full", "grid", "dual", "pointwise"}:
        raise WorkerInputError("fit_protocol.esp_chunk_policy: unsupported value")
    safety = _finite_number(protocol["esp_safety_factor"], "fit_protocol.esp_safety_factor")
    if not 0 < safety <= 1:
        raise WorkerInputError("fit_protocol.esp_safety_factor: expected value in (0, 1]")
    if not isinstance(protocol["source"], str) or not protocol["source"]:
        raise WorkerInputError("fit_protocol.source: expected non-empty string")
    if not isinstance(protocol["radius_overrides"], Mapping):
        raise WorkerInputError("fit_protocol.radius_overrides: expected object")
    for element, radius in protocol["radius_overrides"].items():
        if not isinstance(element, str) or not re.fullmatch(r"[A-Z][a-z]?", element):
            raise WorkerInputError(f"fit_protocol.radius_overrides: invalid element {element!r}")
        if _finite_number(radius, f"fit_protocol.radius_overrides.{element}") <= 0:
            raise WorkerInputError(f"fit_protocol.radius_overrides.{element}: expected positive value")
    if not isinstance(protocol["equivalence_groups"], list):
        raise WorkerInputError("fit_protocol.equivalence_groups: expected array")
    known_ids = set(model_atom_ids)
    seen_ids: set[str] = set()
    for index, item in enumerate(protocol["equivalence_groups"]):
        group = _strict_object(item, {"model_id", "atom_ids", "source"}, f"equivalence_groups[{index}]")
        if group["model_id"] != model_id or not isinstance(group["source"], str) or not group["source"]:
            raise WorkerInputError(f"equivalence_groups[{index}]: model identity/source mismatch")
        if (
            not isinstance(group["atom_ids"], list)
            or len(group["atom_ids"]) < 2
            or not all(isinstance(atom_id, str) for atom_id in group["atom_ids"])
            or len(group["atom_ids"]) != len(set(group["atom_ids"]))
            or not set(group["atom_ids"]) <= known_ids
        ):
            raise WorkerInputError(f"equivalence_groups[{index}].atom_ids: invalid model atom scope")
        overlap = seen_ids & set(group["atom_ids"])
        if overlap:
            raise WorkerInputError(f"equivalence_groups[{index}].atom_ids: overlapping scope")
        seen_ids.update(group["atom_ids"])
    from Xponge.qm.resp_basis import resolve_resp_basis

    try:
        resolved_basis = resolve_resp_basis(protocol["basis_family"], elements)
    except ValueError as exc:
        raise WorkerInputError(f"fit_protocol.basis_family: {exc}") from exc
    if metal_elements:
        missing = sorted(set(metal_elements) - set(resolved_basis.ecp or {}))
        if missing:
            raise WorkerInputError(
                "fit_protocol.metal_basis_policy: resolved ECP does not cover " + ",".join(missing)
            )
    if protocol["backend"] == "psi4" and resolved_basis.ecp:
        raise WorkerInputError("Psi4 RESP cannot consume element-mapped ECP basis families; use PySCF")
    return protocol


def _validate_model(
    value: Any,
    *,
    required_purpose: str = "large",
) -> tuple[Mapping[str, Any], tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    model = _strict_object(
        value,
        {
            "external_id", "model_hash", "purpose", "coordinate_unit", "atomic_charge_role",
            "electronic_state", "atoms", "bonds", "links", "linear_constraints",
        },
        "model",
    )
    if not isinstance(model["external_id"], str) or not model["external_id"]:
        raise WorkerInputError("model.external_id: expected non-empty string")
    if not isinstance(model["model_hash"], str) or len(model["model_hash"]) != 64:
        raise WorkerInputError("model.model_hash: expected sha256")
    if model["purpose"] != required_purpose or model["coordinate_unit"] != "angstrom":
        raise WorkerInputError(
            f"model: expected a {required_purpose} model in angstrom"
        )
    if model["atomic_charge_role"] not in {"absent", "initial", "fixed", "reference", "fitted"}:
        raise WorkerInputError("model.atomic_charge_role: unsupported value")
    state = _strict_object(
        model["electronic_state"], {"selection_id", "net_charge", "spin_multiplicity", "source"},
        "model.electronic_state",
    )
    if not isinstance(state["selection_id"], str) or not state["selection_id"] or not isinstance(state["source"], str) or not state["source"]:
        raise WorkerInputError("model.electronic_state: selection/source required")
    if (
        isinstance(state["net_charge"], bool)
        or not isinstance(state["net_charge"], int)
        or isinstance(state["spin_multiplicity"], bool)
        or not isinstance(state["spin_multiplicity"], int)
        or state["spin_multiplicity"] < 1
    ):
        raise WorkerInputError("model.electronic_state: invalid charge or multiplicity")
    if not isinstance(model["atoms"], list) or not model["atoms"]:
        raise WorkerInputError("model.atoms: expected non-empty array")
    atom_ids: list[str] = []
    elements: list[str] = []
    metal_elements: set[str] = set()
    for index, raw_atom in enumerate(model["atoms"]):
        atom = _strict_object(
            raw_atom,
            {"model_atom_id", "element", "coordinates", "initial_charge", "is_metal"},
            f"model.atoms[{index}]",
        )
        if not isinstance(atom["model_atom_id"], str) or not atom["model_atom_id"] or atom["model_atom_id"] in atom_ids:
            raise WorkerInputError(f"model.atoms[{index}].model_atom_id: invalid or duplicate")
        if not isinstance(atom["element"], str) or not re.fullmatch(r"[A-Z][a-z]?", atom["element"]):
            raise WorkerInputError(f"model.atoms[{index}].element: invalid")
        if not isinstance(atom["coordinates"], list) or len(atom["coordinates"]) != 3:
            raise WorkerInputError(f"model.atoms[{index}].coordinates: expected xyz array")
        for axis, coordinate in enumerate(atom["coordinates"]):
            _finite_number(coordinate, f"model.atoms[{index}].coordinates[{axis}]")
        if atom["initial_charge"] is not None:
            _finite_number(atom["initial_charge"], f"model.atoms[{index}].initial_charge")
        if not isinstance(atom["is_metal"], bool):
            raise WorkerInputError(f"model.atoms[{index}].is_metal: expected boolean")
        atom_ids.append(atom["model_atom_id"])
        elements.append(atom["element"])
        if atom["is_metal"]:
            metal_elements.add(atom["element"])
    atom_id_set = set(atom_ids)
    endpoint_keys: set[tuple[str, str]] = set()
    for collection_name, required in (("bonds", {"model_atom_ids", "order"}), ("links", {"model_atom_ids"})):
        collection = model[collection_name]
        if not isinstance(collection, list):
            raise WorkerInputError(f"model.{collection_name}: expected array")
        for index, raw_edge in enumerate(collection):
            edge = _strict_object(raw_edge, required, f"model.{collection_name}[{index}]")
            endpoints = edge["model_atom_ids"]
            if (
                not isinstance(endpoints, list)
                or len(endpoints) != 2
                or endpoints[0] == endpoints[1]
                or not set(endpoints) <= atom_id_set
            ):
                raise WorkerInputError(f"model.{collection_name}[{index}].model_atom_ids: invalid endpoints")
            key = tuple(sorted(endpoints))
            if key in endpoint_keys:
                raise WorkerInputError(f"model.{collection_name}[{index}]: duplicate edge endpoints")
            endpoint_keys.add(key)
            if collection_name == "bonds" and _finite_number(edge["order"], f"model.bonds[{index}].order") <= 0:
                raise WorkerInputError(f"model.bonds[{index}].order: expected positive value")
    if (
        not isinstance(model["linear_constraints"], list)
        or (required_purpose == "large" and not model["linear_constraints"])
    ):
        raise WorkerInputError(
            "model.linear_constraints: expected an array, non-empty for RESP models"
        )
    constraint_ids: set[str] = set()
    for index, raw_constraint in enumerate(model["linear_constraints"]):
        path = f"model.linear_constraints[{index}]"
        constraint = _strict_object(
            raw_constraint,
            {
                "constraint_id", "role", "atom_ids", "coefficients",
                "target_charge", "source",
            },
            path,
        )
        atom_scope = constraint["atom_ids"]
        coefficients = constraint["coefficients"]
        if (
            not isinstance(constraint["constraint_id"], str)
            or not constraint["constraint_id"]
            or constraint["constraint_id"] in constraint_ids
            or not isinstance(constraint["role"], str)
            or not constraint["role"]
            or not isinstance(constraint["source"], str)
            or not constraint["source"]
            or not isinstance(atom_scope, list)
            or not atom_scope
            or len(atom_scope) != len(set(atom_scope))
            or not set(atom_scope) <= atom_id_set
            or not isinstance(coefficients, list)
            or len(coefficients) != len(atom_scope)
        ):
            raise WorkerInputError(f"{path}: invalid metadata or atom scope")
        for coefficient_index, coefficient in enumerate(coefficients):
            _finite_number(coefficient, f"{path}.coefficients[{coefficient_index}]")
        _finite_number(constraint["target_charge"], f"{path}.target_charge")
        constraint_ids.add(constraint["constraint_id"])
    from .contracts import ElectronicState, validate_electronic_state

    validate_electronic_state(ElectronicState(state["net_charge"], state["spin_multiplicity"]), elements)
    return model, tuple(atom_ids), tuple(elements), tuple(sorted(metal_elements))


def _execute(value: Any) -> dict[str, Any]:
    request = _strict_object(value, {"protocol_version", "model", "fit_protocol"}, "request")
    if request["protocol_version"] != WORKER_PROTOCOL_VERSION:
        raise WorkerInputError("unsupported protocol_version")
    model, atom_ids, elements, metal_elements = _validate_model(request["model"])
    protocol = _validate_fit_protocol(
        request["fit_protocol"], model["external_id"], atom_ids, elements, metal_elements
    )
    atom_index = {atom_id: index for index, atom_id in enumerate(atom_ids)}

    from Xponge.assign import Assign
    from Xponge.assign.resp import resp_fit
    from .charge_fit import ModelChargeArtifact
    from .contracts import _canonicalize

    assign = Assign(name=model["external_id"])
    for atom in model["atoms"]:
        x, y, z = atom["coordinates"]
        assign.add_atom(
            atom["element"], x, y, z, name=atom["model_atom_id"],
            charge=0.0 if atom["initial_charge"] is None else atom["initial_charge"],
        )
    for edge in model["bonds"]:
        first, second = edge["model_atom_ids"]
        assign.add_bond(atom_index[first], atom_index[second], edge["order"])
    for edge in model["links"]:
        first, second = edge["model_atom_ids"]
        assign.add_bond(atom_index[first], atom_index[second], 1.0)
    equivalence = [
        [atom_index[atom_id] for atom_id in group["atom_ids"]]
        for group in protocol["equivalence_groups"]
    ]
    constraint_matrix = []
    constraint_targets = []
    for constraint in model["linear_constraints"]:
        row = [0.0] * len(atom_ids)
        for atom_id, coefficient in zip(
            constraint["atom_ids"], constraint["coefficients"]
        ):
            row[atom_index[atom_id]] = coefficient
        constraint_matrix.append(row)
        constraint_targets.append(constraint["target_charge"])
    state = model["electronic_state"]
    result = resp_fit(
        assign,
        basis=protocol["basis_family"],
        opt=False,
        charge=state["net_charge"],
        spin=state["spin_multiplicity"] - 1,
        extra_equivalence=equivalence,
        grid_density=protocol["grid_density"],
        grid_cell_layer=protocol["grid_cell_layer"],
        radius=dict(protocol["radius_overrides"]),
        a1=protocol["restraint_a1"],
        a2=protocol["restraint_a2"],
        two_stage=True,
        only_esp=False,
        backend=protocol["backend"],
        core="python",
        esp_memory_limit=protocol["esp_memory_limit_bytes"],
        esp_chunk_policy=protocol["esp_chunk_policy"],
        esp_safety_factor=protocol["esp_safety_factor"],
        constraint_matrix=constraint_matrix,
        constraint_targets=constraint_targets,
        return_metadata=True,
        return_diagnostics=True,
    )
    if not isinstance(result, Mapping) or set(result) != {"charges", "metadata", "diagnostics"}:
        raise RuntimeError("RESP backend returned an invalid result object")
    raw_charges = list(result["charges"])
    if len(raw_charges) != len(atom_ids):
        raise RuntimeError("RESP backend returned the wrong charge count")
    charges = tuple(_finite_number(charge, f"result.charges[{index}]") for index, charge in enumerate(raw_charges))
    backend_version = _package_version(protocol["backend"])
    xponge_version = _package_version("Xponge")
    provider = f"xponge-resp:{protocol['backend']}"
    provider_version = f"xponge={xponge_version};{protocol['backend']}={backend_version}"
    artifact = ModelChargeArtifact(
        model_id=model["external_id"],
        model_hash=model["model_hash"],
        atom_order=atom_ids,
        charges=charges,
        atomic_charge_role="fitted",
        provider=provider,
        provider_version=provider_version,
    ).with_computed_hash()
    request_hash = _canonical_hash(request)
    fit_report = {
        "model_id": model["external_id"],
        "model_hash": model["model_hash"],
        "fit_input_hash": protocol["fit_input_hash"],
        "request_hash": request_hash,
        "backend": protocol["backend"],
        "backend_version": backend_version,
        "xponge_version": xponge_version,
        "basis_family": protocol["basis_family"],
        "metal_basis_policy": protocol["metal_basis_policy"],
        "basis_source": protocol["basis_source"],
        "metal_elements": list(metal_elements),
        "net_charge": state["net_charge"],
        "spin_multiplicity": state["spin_multiplicity"],
        "backend_spin_2s": state["spin_multiplicity"] - 1,
        "coordinate_unit": "angstrom",
        "geometry_locked": True,
        "atom_count": len(atom_ids),
        "equivalence_group_count": len(equivalence),
        "fitted_charge_sum": float(sum(charges)),
        "target_charge": float(state["net_charge"]),
        "charge_residual": float(sum(charges) - state["net_charge"]),
        "constraint_count": len(constraint_matrix),
        "constraint_diagnostics": result["diagnostics"],
        "setup_metadata": result["metadata"],
        "source": protocol["source"],
    }
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "model_id": model["external_id"],
        "model_hash": model["model_hash"],
        "fit_input_hash": protocol["fit_input_hash"],
        "request_hash": request_hash,
        "artifact": _canonicalize(artifact),
        "fit_report": fit_report,
    }
    response["response_hash"] = _canonical_hash(response)
    return response


def main() -> int:
    try:
        value = json.load(sys.stdin)
        with redirect_stdout(sys.stderr):
            response = _execute(value)
        json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 0
    except Exception as exc:
        json.dump(
            {
                "ok": False,
                "protocol_version": WORKER_PROTOCOL_VERSION,
                "error": {"type": type(exc).__name__, "message": str(exc)},
            },
            sys.stdout,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
