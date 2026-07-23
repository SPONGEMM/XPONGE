"""Isolated GAFF/GAFF2 worker for topology-locked assignment.

This module is intentionally private.  The public process imports neither the
Amber force-field modules nor their process-global registries; it communicates
with this worker through one JSON document on stdin and one JSON document on
stdout.
"""

from __future__ import annotations

from contextlib import redirect_stdout
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from typing import Any, Mapping


WORKER_PROTOCOL_VERSION = 2
SUPPORTED_PROVIDERS = frozenset({"gaff", "gaff2"})
SUPPORTED_CHARGE_METHODS = frozenset({"gasteiger", "tpacm4"})


class WorkerInputError(ValueError):
    """Invalid private worker request."""


def _canonical_hash(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _strict_object(
    value: Any,
    *,
    required: set[str],
    path: str,
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise WorkerInputError(f"{path}: expected object")
    missing = required - set(value)
    unknown = set(value) - required
    if missing:
        raise WorkerInputError(f"{path}: missing {','.join(sorted(missing))}")
    if unknown:
        raise WorkerInputError(f"{path}: unknown {','.join(sorted(unknown))}")
    return value


def _finite_number(value: Any, path: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise WorkerInputError(f"{path}: expected finite number")
    return float(value)


def _validate_component(value: Any, index: int) -> dict[str, Any]:
    path = f"components[{index}]"
    component = _strict_object(
        value,
        required={
            "external_id", "net_formal_charge", "proof_hash", "artifact_hash", "atoms", "bonds",
        },
        path=path,
    )
    if not isinstance(component["external_id"], str) or not component["external_id"]:
        raise WorkerInputError(f"{path}.external_id: expected non-empty string")
    if isinstance(component["net_formal_charge"], bool) or not isinstance(component["net_formal_charge"], int):
        raise WorkerInputError(f"{path}.net_formal_charge: expected integer")
    for name in ("proof_hash", "artifact_hash"):
        if not isinstance(component[name], str) or len(component[name]) != 64:
            raise WorkerInputError(f"{path}.{name}: expected sha256")
    if not isinstance(component["atoms"], list) or not component["atoms"]:
        raise WorkerInputError(f"{path}.atoms: expected non-empty array")
    if not isinstance(component["bonds"], list):
        raise WorkerInputError(f"{path}.bonds: expected array")

    atom_ids: list[str] = []
    atoms: list[dict[str, Any]] = []
    for atom_index, raw_atom in enumerate(component["atoms"]):
        atom_path = f"{path}.atoms[{atom_index}]"
        atom = _strict_object(
            raw_atom,
            required={"external_id", "element", "coordinates", "name", "formal_charge"},
            path=atom_path,
        )
        if not isinstance(atom["external_id"], str) or not atom["external_id"]:
            raise WorkerInputError(f"{atom_path}.external_id: expected non-empty string")
        if not isinstance(atom["element"], str) or not atom["element"]:
            raise WorkerInputError(f"{atom_path}.element: expected non-empty string")
        if not isinstance(atom["name"], str):
            raise WorkerInputError(f"{atom_path}.name: expected string")
        if not isinstance(atom["coordinates"], list) or len(atom["coordinates"]) != 3:
            raise WorkerInputError(f"{atom_path}.coordinates: expected three values")
        coordinates = [
            _finite_number(coordinate, f"{atom_path}.coordinates[{coordinate_index}]")
            for coordinate_index, coordinate in enumerate(atom["coordinates"])
        ]
        formal_charge = atom["formal_charge"]
        if formal_charge is not None and (
            isinstance(formal_charge, bool) or not isinstance(formal_charge, int)
        ):
            raise WorkerInputError(f"{atom_path}.formal_charge: expected integer or null")
        atom_ids.append(atom["external_id"])
        atoms.append({**atom, "coordinates": coordinates})
    if len(atom_ids) != len(set(atom_ids)):
        raise WorkerInputError(f"{path}.atoms: duplicate external_id")

    atom_id_set = set(atom_ids)
    bond_ids: set[str] = set()
    bonds: list[dict[str, Any]] = []
    adjacency = {atom_id: set() for atom_id in atom_ids}
    for bond_index, raw_bond in enumerate(component["bonds"]):
        bond_path = f"{path}.bonds[{bond_index}]"
        bond = _strict_object(
            raw_bond,
            required={"external_id", "atom_ids", "order", "semantic", "source"},
            path=bond_path,
        )
        if not isinstance(bond["external_id"], str) or not bond["external_id"]:
            raise WorkerInputError(f"{bond_path}.external_id: expected non-empty string")
        if bond["external_id"] in bond_ids:
            raise WorkerInputError(f"{bond_path}.external_id: duplicate")
        bond_ids.add(bond["external_id"])
        endpoints = bond["atom_ids"]
        if (
            not isinstance(endpoints, list)
            or len(endpoints) != 2
            or not all(isinstance(atom_id, str) for atom_id in endpoints)
            or endpoints[0] == endpoints[1]
            or not set(endpoints) <= atom_id_set
        ):
            raise WorkerInputError(f"{bond_path}.atom_ids: invalid endpoints")
        order = _finite_number(bond["order"], f"{bond_path}.order")
        if order <= 0 or order > 3:
            raise WorkerInputError(f"{bond_path}.order: expected 0 < order <= 3")
        for name in ("semantic", "source"):
            if not isinstance(bond[name], str) or not bond[name]:
                raise WorkerInputError(f"{bond_path}.{name}: expected non-empty string")
        adjacency[endpoints[0]].add(endpoints[1])
        adjacency[endpoints[1]].add(endpoints[0])
        bonds.append({**bond, "order": order})
    visited = {atom_ids[0]}
    pending = [atom_ids[0]]
    while pending:
        current = pending.pop()
        for neighbor in adjacency[current]:
            if neighbor not in visited:
                visited.add(neighbor)
                pending.append(neighbor)
    if visited != atom_id_set:
        raise WorkerInputError(f"{path}: component graph is disconnected")
    return {**component, "atoms": atoms, "bonds": bonds}


def _validate_base_charge_input(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    charge_input = _strict_object(
        value,
        required={"schema_version", "method", "source", "charge_input_hash"},
        path="request.base_charge_input",
    )
    if charge_input["schema_version"] != 1:
        raise WorkerInputError("request.base_charge_input.schema_version: unsupported")
    if charge_input["method"] not in SUPPORTED_CHARGE_METHODS:
        raise WorkerInputError("request.base_charge_input.method: unsupported")
    if not isinstance(charge_input["source"], str) or not charge_input["source"]:
        raise WorkerInputError("request.base_charge_input.source: expected non-empty string")
    payload = dict(charge_input)
    supplied_hash = payload.pop("charge_input_hash")
    if not isinstance(supplied_hash, str) or supplied_hash != _canonical_hash(payload):
        raise WorkerInputError("request.base_charge_input.charge_input_hash: stale")
    return dict(charge_input)


def _validate_request(value: Any) -> dict[str, Any]:
    request = _strict_object(
        value,
        required={
            "protocol_version", "provider", "topology_hash", "complete_missing_parameters",
            "base_charge_input", "components",
        },
        path="request",
    )
    if request["protocol_version"] != WORKER_PROTOCOL_VERSION:
        raise WorkerInputError("request.protocol_version: unsupported")
    provider = request["provider"]
    if provider not in SUPPORTED_PROVIDERS:
        raise WorkerInputError(f"request.provider: unsupported {provider!r}")
    if not isinstance(request["topology_hash"], str) or len(request["topology_hash"]) != 64:
        raise WorkerInputError("request.topology_hash: expected sha256")
    if not isinstance(request["complete_missing_parameters"], bool):
        raise WorkerInputError("request.complete_missing_parameters: expected boolean")
    if not isinstance(request["components"], list) or not request["components"]:
        raise WorkerInputError("request.components: expected non-empty array")
    components = [_validate_component(component, index) for index, component in enumerate(request["components"])]
    component_ids = [component["external_id"] for component in components]
    if len(component_ids) != len(set(component_ids)):
        raise WorkerInputError("request.components: duplicate external_id")
    return {
        **request,
        "base_charge_input": _validate_base_charge_input(request["base_charge_input"]),
        "components": components,
    }


def _provider_version() -> str:
    try:
        return importlib.metadata.version("Xponge")
    except importlib.metadata.PackageNotFoundError:
        return "source-tree"


def _canonical_term_atoms(kind: str, atom_ids: list[str]) -> tuple[str, ...]:
    values = tuple(atom_ids)
    if kind in {"bond", "angle", "proper_dihedral"}:
        reverse = values[::-1]
        return min(values, reverse)
    # Xponge improper entities place the central atom at index 2.  The neutral
    # contract uses (outer, center, outer, outer).
    center = values[2]
    outers = sorted((values[0], values[1], values[3]))
    return (outers[0], center, outers[1], outers[2])


def _bonded_records(restype: Any, atom_id_by_object: Mapping[Any, str], source: str) -> tuple[list[dict[str, Any]], list[str]]:
    records: list[dict[str, Any]] = []
    ignored: list[str] = []
    for xponge_kind, entities in restype.bonded_forces.items():
        if xponge_kind not in {"bond", "angle", "dihedral", "improper"}:
            if entities:
                ignored.append(xponge_kind)
            continue
        kind = {
            "bond": "bond",
            "angle": "angle",
            "dihedral": "proper_dihedral",
            "improper": "improper_dihedral",
        }[xponge_kind]
        for entity in entities:
            atom_ids = _canonical_term_atoms(kind, [atom_id_by_object[atom] for atom in entity.atoms])
            if kind == "bond":
                parameters = {
                    "force_constant": float(entity.k),
                    "equilibrium_distance": float(entity.b),
                    "force_constant_unit": "kcal/mol/angstrom^2",
                    "distance_unit": "angstrom",
                    "type_name": entity.type.name,
                }
            elif kind == "angle":
                parameters = {
                    "force_constant": float(entity.k),
                    "equilibrium_angle": float(entity.b),
                    "force_constant_unit": "kcal/mol/radian^2",
                    "angle_unit": "radian",
                    "type_name": entity.type.name,
                }
            elif kind == "proper_dihedral":
                parameters = {
                    "force_constants": [float(value) for value in entity.ks],
                    "phases": [float(value) for value in entity.phi0s],
                    "periodicities": [int(value) for value in entity.periodicitys],
                    "energy_unit": "kcal/mol",
                    "angle_unit": "radian",
                    "type_name": entity.type.name,
                }
            else:
                parameters = {
                    "force_constant": float(entity.k),
                    "phase": float(entity.phi0),
                    "periodicity": int(entity.periodicity),
                    "energy_unit": "kcal/mol",
                    "angle_unit": "radian",
                    "type_name": entity.type.name,
                }
            records.append({
                "kind": kind,
                "atom_ids": list(atom_ids),
                "parameters": parameters,
                "source": source,
            })
    records.sort(key=lambda item: (item["kind"], item["atom_ids"], _canonical_hash(item["parameters"])))
    return records, sorted(set(ignored))


def _assign_component(
    component: Mapping[str, Any],
    *,
    provider: str,
    provider_module: Any,
    complete_missing_parameters: bool,
    base_charge_input: Mapping[str, Any] | None,
    component_index: int,
) -> dict[str, Any]:
    from Xponge.assign import Assign
    from Xponge.build import build_bonded_force
    from Xponge.forcefield.base.lj_base import LJType

    assign = Assign(f"metal_assignment_{component_index}")
    atom_index: dict[str, int] = {}
    for index, atom in enumerate(component["atoms"]):
        atom_index[atom["external_id"]] = index
        assign.add_atom(
            atom["element"],
            *atom["coordinates"],
            name=atom["name"] or f"A{index + 1}",
            charge=0.0,
        )
        if atom["formal_charge"] is not None:
            assign.formal_charge[index] = atom["formal_charge"]
    for bond in component["bonds"]:
        atom1, atom2 = bond["atom_ids"]
        # The locked Chemcore order is used verbatim, including delocalized 1.5
        # orders.  No connectivity or bond-order inference API is called.
        assign.add_bond(atom_index[atom1], atom_index[atom2], bond["order"])
    assign.determine_atom_type(provider)
    charges: dict[str, float] = {}
    if base_charge_input is not None:
        assign.calculate_charge(
            base_charge_input["method"],
            charge=component["net_formal_charge"],
        )
        raw_charges = list(assign.charge)
        if len(raw_charges) != len(component["atoms"]) or any(
            isinstance(charge, bool)
            or not isinstance(charge, (int, float))
            or not math.isfinite(charge)
            for charge in raw_charges
        ):
            raise WorkerInputError(
                f"components[{component_index}]: charge provider returned invalid coverage"
            )
        if not math.isclose(
            sum(float(charge) for charge in raw_charges),
            float(component["net_formal_charge"]),
            abs_tol=1e-6,
        ):
            raise WorkerInputError(
                f"components[{component_index}]: charge provider did not conserve component charge"
            )
        charges = {
            atom["external_id"]: float(raw_charges[index])
            for index, atom in enumerate(component["atoms"])
        }
    restype = assign.to_residuetype(f"MA_{component_index}")

    frcmod_sha256 = ""
    if complete_missing_parameters:
        completion = getattr(provider_module, f"parmchk2_{provider}")
        with TemporaryDirectory(prefix="xponge-metal-assignment-") as tempdir:
            frcmod_path = Path(tempdir) / f"{provider}.frcmod"
            completion(restype, str(frcmod_path), direct_load=True, keep=True)
            frcmod_sha256 = hashlib.sha256(frcmod_path.read_bytes()).hexdigest()
    build_bonded_force(restype)

    source = f"xponge:{provider}:{_provider_version()}"
    atom_types: dict[str, str] = {}
    masses: dict[str, float] = {}
    lj_parameters: dict[str, dict[str, Any]] = {}
    atom_id_by_object: dict[Any, str] = {}
    for index, atom in enumerate(restype.atoms):
        external_id = component["atoms"][index]["external_id"]
        atom_id_by_object[atom] = external_id
        atom_types[external_id] = atom.type.name
        masses[external_id] = float(atom.mass)
        lj_type = LJType.get_type(f"{atom.LJtype}-{atom.LJtype}")
        lj_parameters[external_id] = {
            "epsilon": float(lj_type.epsilon),
            "rmin": float(lj_type.rmin),
            "energy_unit": "kcal/mol",
            "length_unit": "angstrom",
            "source": source,
        }

    records, ignored_force_kinds = _bonded_records(restype, atom_id_by_object, source)
    term_counts: dict[str, int] = {}
    bonded_parameters: dict[str, Any] = {}
    for record in records:
        stem = f"base:{provider}:{component['external_id']}:{record['kind']}:{'-'.join(record['atom_ids'])}"
        occurrence = term_counts.get(stem, 0)
        term_counts[stem] = occurrence + 1
        term_id = stem if occurrence == 0 else f"{stem}:{occurrence + 1}"
        bonded_parameters[term_id] = record

    output = {
        "external_id": component["external_id"],
        "proof_hash": component["proof_hash"],
        "artifact_hash": component["artifact_hash"],
        "net_formal_charge": component["net_formal_charge"],
        "atom_types": atom_types,
        "charges": charges,
        "masses": masses,
        "lj_parameters": lj_parameters,
        "bonded_parameters": bonded_parameters,
        "frcmod_sha256": frcmod_sha256,
        "ignored_force_kinds": ignored_force_kinds,
    }
    output["output_hash"] = _canonical_hash(output)
    return output


def _execute(value: Any) -> dict[str, Any]:
    request = _validate_request(value)
    provider = request["provider"]
    provider_module = importlib.import_module(f"Xponge.forcefield.amber.{provider}")
    components = [
        _assign_component(
            component,
            provider=provider,
            provider_module=provider_module,
            complete_missing_parameters=request["complete_missing_parameters"],
            base_charge_input=request["base_charge_input"],
            component_index=index,
        )
        for index, component in enumerate(request["components"])
    ]
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": provider,
        "topology_hash": request["topology_hash"],
        "provider_version": _provider_version(),
        "components": components,
    }
    response["response_hash"] = _canonical_hash(response)
    return response


def main() -> int:
    try:
        request = json.load(sys.stdin)
        # Xponge force-field imports print references.  Keep stdout as a strict
        # one-document machine channel and route all provider chatter to stderr.
        with redirect_stdout(sys.stderr):
            response = _execute(request)
        json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 0
    except Exception as exc:  # structured boundary; parent retains stderr too
        error = {
            "ok": False,
            "protocol_version": WORKER_PROTOCOL_VERSION,
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
            },
        }
        json.dump(error, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
