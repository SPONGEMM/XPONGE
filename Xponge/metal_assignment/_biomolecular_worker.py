"""Isolated standard-biomolecular base-assignment worker.

The parent process supplies a Chemcore-locked residue graph.  This worker only
matches that graph to an existing Xponge force-field template; it never repairs
atoms, guesses connectivity, or changes residue boundaries.
"""

from __future__ import annotations

from contextlib import redirect_stdout
import importlib
import importlib.metadata
import json
import sys
from typing import Any, Mapping

from ._gaff_worker import _bonded_records, _canonical_hash


WORKER_PROTOCOL_VERSION = 1
SUPPORTED_PROVIDER_FORCE_FIELDS = {
    "standard_biomolecular": frozenset({"ff14sb", "ff19sb"}),
    "glycam": frozenset({"glycam_06j"}),
    "lipid": frozenset({"lipid21"}),
}
SUPPORTED_WATER_MODELS = frozenset({"tip3p", "spce", "opc"})
WATER_NAMES = frozenset({"HOH", "H2O", "WAT"})


class WorkerInputError(ValueError):
    """Invalid private worker request."""


def _strict_object(value: Any, *, required: set[str], path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise WorkerInputError(f"{path}: expected object")
    missing = required - set(value)
    unknown = set(value) - required
    if missing:
        raise WorkerInputError(f"{path}: missing {','.join(sorted(missing))}")
    if unknown:
        raise WorkerInputError(f"{path}: unknown {','.join(sorted(unknown))}")
    return value


def _validate_component(value: Any, index: int) -> dict[str, Any]:
    path = f"components[{index}]"
    component = _strict_object(
        value,
        required={
            "external_id", "proof_hash", "artifact_hash", "residue_name", "polymer_position",
            "template_id", "atoms", "bonds",
        },
        path=path,
    )
    for name in ("external_id", "proof_hash", "artifact_hash", "residue_name", "polymer_position", "template_id"):
        if not isinstance(component[name], str):
            raise WorkerInputError(f"{path}.{name}: expected string")
    if not component["external_id"] or not component["residue_name"]:
        raise WorkerInputError(f"{path}: external_id and residue_name must be non-empty")
    if component["polymer_position"] not in {
        "nonpolymer", "unknown", "internal", "start", "end", "single",
    }:
        raise WorkerInputError(f"{path}.polymer_position: invalid value")
    if not isinstance(component["atoms"], list) or not component["atoms"]:
        raise WorkerInputError(f"{path}.atoms: expected non-empty array")
    atoms: list[dict[str, str]] = []
    external_ids: set[str] = set()
    atom_names: set[str] = set()
    for atom_index, raw_atom in enumerate(component["atoms"]):
        atom_path = f"{path}.atoms[{atom_index}]"
        atom = _strict_object(raw_atom, required={"external_id", "atom_name", "element"}, path=atom_path)
        if not all(isinstance(atom[name], str) and atom[name] for name in atom):
            raise WorkerInputError(f"{atom_path}: atom fields must be non-empty strings")
        if atom["external_id"] in external_ids or atom["atom_name"] in atom_names:
            raise WorkerInputError(f"{atom_path}: duplicate atom identity")
        external_ids.add(atom["external_id"])
        atom_names.add(atom["atom_name"])
        atoms.append(dict(atom))
    if not isinstance(component["bonds"], list):
        raise WorkerInputError(f"{path}.bonds: expected array")
    bonds: list[dict[str, Any]] = []
    bond_ids: set[str] = set()
    for bond_index, raw_bond in enumerate(component["bonds"]):
        bond_path = f"{path}.bonds[{bond_index}]"
        bond = _strict_object(raw_bond, required={"external_id", "atom_ids"}, path=bond_path)
        endpoints = bond["atom_ids"]
        if (
            not isinstance(bond["external_id"], str)
            or not bond["external_id"]
            or bond["external_id"] in bond_ids
            or not isinstance(endpoints, list)
            or len(endpoints) != 2
            or endpoints[0] == endpoints[1]
            or not all(isinstance(atom_id, str) and atom_id in external_ids for atom_id in endpoints)
        ):
            raise WorkerInputError(f"{bond_path}: invalid bond")
        bond_ids.add(bond["external_id"])
        bonds.append({"external_id": bond["external_id"], "atom_ids": list(endpoints)})
    return {**component, "atoms": atoms, "bonds": bonds}


def _validate_request(value: Any) -> dict[str, Any]:
    request = _strict_object(
        value,
        required={
            "protocol_version", "provider", "force_field", "water_model", "topology_hash",
            "components", "links",
        },
        path="request",
    )
    if request["protocol_version"] != WORKER_PROTOCOL_VERSION:
        raise WorkerInputError("request: unsupported protocol or provider")
    supported_force_fields = SUPPORTED_PROVIDER_FORCE_FIELDS.get(request["provider"])
    if supported_force_fields is None or request["force_field"] not in supported_force_fields:
        raise WorkerInputError(f"request.force_field: unsupported {request['force_field']!r}")
    if request["water_model"] not in SUPPORTED_WATER_MODELS:
        raise WorkerInputError(f"request.water_model: unsupported {request['water_model']!r}")
    if not isinstance(request["topology_hash"], str) or len(request["topology_hash"]) != 64:
        raise WorkerInputError("request.topology_hash: expected sha256")
    if not isinstance(request["components"], list) or not request["components"]:
        raise WorkerInputError("request.components: expected non-empty array")
    components = [_validate_component(component, index) for index, component in enumerate(request["components"])]
    ids = [component["external_id"] for component in components]
    if len(ids) != len(set(ids)):
        raise WorkerInputError("request.components: duplicate external_id")
    if not isinstance(request["links"], list):
        raise WorkerInputError("request.links: expected array")
    atom_ids = {
        atom["external_id"]
        for component in components
        for atom in component["atoms"]
    }
    links = []
    link_ids = set()
    for index, raw_link in enumerate(request["links"]):
        path = f"request.links[{index}]"
        link = _strict_object(raw_link, required={"external_id", "atom_ids"}, path=path)
        endpoints = link["atom_ids"]
        if (
            not isinstance(link["external_id"], str)
            or not link["external_id"]
            or link["external_id"] in link_ids
            or not isinstance(endpoints, list)
            or len(endpoints) != 2
            or endpoints[0] == endpoints[1]
            or not all(isinstance(atom_id, str) and atom_id in atom_ids for atom_id in endpoints)
        ):
            raise WorkerInputError(f"{path}: invalid link")
        link_ids.add(link["external_id"])
        links.append({"external_id": link["external_id"], "atom_ids": list(endpoints)})
    return {**request, "components": components, "links": links}


def _provider_version() -> str:
    try:
        return importlib.metadata.version("Xponge")
    except importlib.metadata.PackageNotFoundError:
        return "source-tree"


def _candidate_template_names(component: Mapping[str, Any]) -> tuple[str, ...]:
    residue_name = component["residue_name"].upper()
    if residue_name in WATER_NAMES:
        return ("WAT", residue_name)
    bases: list[str] = []
    template_id = component["template_id"].upper()
    if template_id:
        bases.append(template_id)
    bases.append(residue_name)
    if residue_name == "HIS":
        bases.extend(("HID", "HIE", "HIP"))
    elif residue_name == "CYS":
        bases.append("CYX")
    prefix = {"start": "N", "end": "C"}.get(component["polymer_position"], "")
    candidates: list[str] = []
    for base in bases:
        candidates.append(base)
        if prefix and not base.startswith(prefix):
            candidates.insert(0, prefix + base)
    return tuple(dict.fromkeys(candidates))


def _template_bond_names(template: Any) -> set[tuple[str, str]]:
    bonds: set[tuple[str, str]] = set()
    for atom1, neighbors in template.connectivity.items():
        for atom2 in neighbors:
            bonds.add(tuple(sorted((atom1.name, atom2.name))))
    return bonds


def _distance_constraint_names(
    component: Mapping[str, Any],
    template_bonds: set[tuple[str, str]],
    observed_bonds: set[tuple[str, str]],
) -> set[tuple[str, str]]:
    element_by_name = {atom["atom_name"]: atom["element"].upper() for atom in component["atoms"]}
    constraints: set[tuple[str, str]] = set()
    for atom1, atom2 in template_bonds - observed_bonds:
        if element_by_name.get(atom1) != "H" or element_by_name.get(atom2) != "H":
            continue
        neighbors1 = {
            right if left == atom1 else left
            for left, right in template_bonds
            if atom1 in (left, right)
        }
        neighbors2 = {
            right if left == atom2 else left
            for left, right in template_bonds
            if atom2 in (left, right)
        }
        if neighbors1 & neighbors2:
            constraints.add((atom1, atom2))
    return constraints


def _select_template(component: Mapping[str, Any]) -> Any:
    from Xponge.helper import ResidueType

    observed_names = {atom["atom_name"] for atom in component["atoms"]}
    matches: list[Any] = []
    diagnostics: list[str] = []
    for candidate_name in _candidate_template_names(component):
        try:
            candidate = ResidueType.get_type(candidate_name)
        except (KeyError, AssertionError):
            continue
        expected_names = {atom.name for atom in candidate.atoms}
        if expected_names == observed_names:
            if all(candidate is not item for item in matches):
                matches.append(candidate)
        else:
            missing = sorted(expected_names - observed_names)
            extra = sorted(observed_names - expected_names)
            diagnostics.append(f"{candidate_name}:missing={missing}:extra={extra}")
    if len(matches) != 1:
        raise WorkerInputError(
            f"component {component['external_id']}: expected one exact residue template, "
            f"matched={len(matches)}, candidates={_candidate_template_names(component)}, diagnostics={diagnostics}"
        )
    return matches[0]


def _assign_component(component: Mapping[str, Any], source: str) -> dict[str, Any]:
    from Xponge.build import build_bonded_force
    from Xponge.forcefield.base.lj_base import LJType

    template = _select_template(component)
    external_id_by_name = {atom["atom_name"]: atom["external_id"] for atom in component["atoms"]}
    observed_bonds = {
        tuple(sorted((
            next(atom["atom_name"] for atom in component["atoms"] if atom["external_id"] == bond["atom_ids"][0]),
            next(atom["atom_name"] for atom in component["atoms"] if atom["external_id"] == bond["atom_ids"][1]),
        )))
        for bond in component["bonds"]
    }
    expected_bonds = _template_bond_names(template)
    constraint_names = _distance_constraint_names(component, expected_bonds, observed_bonds)
    if observed_bonds != expected_bonds - constraint_names:
        raise WorkerInputError(
            f"component {component['external_id']}: template connectivity mismatch "
            f"missing={sorted(expected_bonds - observed_bonds)}, extra={sorted(observed_bonds - expected_bonds)}"
        )
    build_bonded_force(template)
    atom_types: dict[str, str] = {}
    charges: dict[str, float] = {}
    masses: dict[str, float] = {}
    lj_parameters: dict[str, dict[str, Any]] = {}
    atom_id_by_object: dict[Any, str] = {}
    for atom in template.atoms:
        external_id = external_id_by_name[atom.name]
        atom_id_by_object[atom] = external_id
        atom_types[external_id] = atom.type.name
        charges[external_id] = float(atom.charge)
        masses[external_id] = float(atom.mass)
        lj_type = LJType.get_type(f"{atom.LJtype}-{atom.LJtype}")
        lj_parameters[external_id] = {
            "epsilon": float(lj_type.epsilon),
            "rmin": float(lj_type.rmin),
            "energy_unit": "kcal/mol",
            "length_unit": "angstrom",
            "source": source,
        }
    records, ignored_force_kinds = _bonded_records(template, atom_id_by_object, source)
    active_records: list[dict[str, Any]] = []
    zero_force_terms = 0
    for record in records:
        parameters = record["parameters"]
        if (
            parameters.get("force_constant") == 0
            or ("force_constants" in parameters and not any(parameters["force_constants"]))
        ):
            zero_force_terms += 1
            continue
        active_records.append(record)
    records = active_records
    if zero_force_terms:
        ignored_force_kinds = sorted(set((*ignored_force_kinds, "zero_force_term")))
    constraint_atom_ids = {
        tuple(sorted((external_id_by_name[atom1], external_id_by_name[atom2])))
        for atom1, atom2 in constraint_names
    }
    for record in records:
        if record["kind"] == "bond" and tuple(sorted(record["atom_ids"])) in constraint_atom_ids:
            record["kind"] = "distance_constraint"
    records.sort(key=lambda item: (item["kind"], item["atom_ids"], _canonical_hash(item["parameters"])))
    term_counts: dict[str, int] = {}
    bonded_parameters: dict[str, Any] = {}
    for record in records:
        stem = (
            f"base:{source}:{component['external_id']}:"
            f"{record['kind']}:{'-'.join(record['atom_ids'])}"
        )
        occurrence = term_counts.get(stem, 0)
        term_counts[stem] = occurrence + 1
        term_id = stem if occurrence == 0 else f"{stem}:{occurrence + 1}"
        bonded_parameters[term_id] = record
    output = {
        "external_id": component["external_id"],
        "proof_hash": component["proof_hash"],
        "artifact_hash": component["artifact_hash"],
        "template_name": template.name,
        "atom_types": atom_types,
        "charges": charges,
        "masses": masses,
        "lj_parameters": lj_parameters,
        "bonded_parameters": bonded_parameters,
        "ignored_force_kinds": ignored_force_kinds,
    }
    output["output_hash"] = _canonical_hash(output)
    return output


def _cross_component_terms(
    request: Mapping[str, Any],
    outputs: list[dict[str, Any]],
    source: str,
) -> dict[str, Any]:
    if not request["links"]:
        return {}
    from Xponge.build import build_bonded_force
    from Xponge.helper import Molecule, Residue

    molecule = Molecule("forcefield_assignment_standard_links")
    atom_by_external_id: dict[str, Any] = {}
    component_by_atom: dict[str, str] = {}
    for component, output in zip(request["components"], outputs):
        template = _select_template(component)
        residue = Residue(template, directly_copy=True)
        molecule.add_residue(residue)
        external_by_name = {
            atom["atom_name"]: atom["external_id"]
            for atom in component["atoms"]
        }
        for atom in residue.atoms:
            external_id = external_by_name[atom.name]
            atom_by_external_id[external_id] = atom
            component_by_atom[external_id] = component["external_id"]
        if set(output["atom_types"]) != set(external_by_name.values()):
            raise WorkerInputError(
                f"component {component['external_id']}: assigned atom coverage changed before link completion"
            )
    for link in request["links"]:
        atom1, atom2 = (atom_by_external_id[atom_id] for atom_id in link["atom_ids"])
        molecule.add_residue_link(atom1, atom2)
    build_bonded_force(molecule)
    external_id_by_object = {atom: external_id for external_id, atom in atom_by_external_id.items()}
    records, _ignored = _bonded_records(molecule, external_id_by_object, source)
    records = [
        record
        for record in records
        if len({component_by_atom[atom_id] for atom_id in record["atom_ids"]}) > 1
    ]
    term_counts: dict[str, int] = {}
    terms: dict[str, Any] = {}
    for record in records:
        stem = (
            f"base:{request['provider']}:cross-component:"
            f"{record['kind']}:{'-'.join(record['atom_ids'])}"
        )
        occurrence = term_counts.get(stem, 0)
        term_counts[stem] = occurrence + 1
        term_id = stem if occurrence == 0 else f"{stem}:{occurrence + 1}"
        terms[term_id] = record
    return terms


def _execute(value: Any) -> dict[str, Any]:
    request = _validate_request(value)
    if request["provider"] == "glycam":
        importlib.import_module("Xponge.forcefield.amber.ff19sb")
        for submodule in (
            "d_pyranose", "l_pyranose", "d_furanose", "l_furanose", "glycoprotein",
        ):
            importlib.import_module(f"Xponge.forcefield.amber.glycam_06j.{submodule}")
    else:
        importlib.import_module(f"Xponge.forcefield.amber.{request['force_field']}")
    importlib.import_module(f"Xponge.forcefield.amber.{request['water_model']}")
    source = (
        f"xponge:{request['provider']}:{request['force_field']}:"
        f"{request['water_model']}:{_provider_version()}"
    )
    components = [_assign_component(component, source) for component in request["components"]]
    cross_component_bonded_parameters = _cross_component_terms(request, components, source)
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": request["provider"],
        "topology_hash": request["topology_hash"],
        "provider_version": _provider_version(),
        "force_field": request["force_field"],
        "water_model": request["water_model"],
        "components": components,
        "cross_component_bonded_parameters": cross_component_bonded_parameters,
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
