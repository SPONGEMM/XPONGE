"""Isolated Xponge materialization worker for assigned-system artifacts."""

from __future__ import annotations

from collections import Counter, deque
import json
import sys
import traceback
from typing import Any, Mapping

from .assigned_system import (
    APPLY_WORKER_PROTOCOL_VERSION,
    AssignedForceTerm,
    AssignedSystemArtifact,
    assigned_system_from_dict,
)
from .contracts import ValidationError, _sha256


def _bond_values(term: AssignedForceTerm) -> tuple[float, float]:
    values = term.parameters
    if "force_constant" in values and "equilibrium_distance" in values:
        return float(values["force_constant"]), float(values["equilibrium_distance"])
    if "k" in values and "equilibrium" in values:
        return float(values["k"]), float(values["equilibrium"])
    raise ValidationError("unsupported_bond_parameter_shape", term.external_id)


def _angle_values(term: AssignedForceTerm) -> tuple[float, float]:
    values = term.parameters
    if "force_constant" in values and "equilibrium_angle" in values:
        return float(values["force_constant"]), float(values["equilibrium_angle"])
    if "k" in values and "equilibrium" in values:
        return float(values["k"]), float(values["equilibrium"])
    raise ValidationError("unsupported_angle_parameter_shape", term.external_id)


def _proper_values(term: AssignedForceTerm) -> tuple[list[float], list[float], list[int]]:
    values = term.parameters
    if {"force_constants", "phases", "periodicities"} <= set(values):
        ks = [float(value) for value in values["force_constants"]]
        phases = [float(value) for value in values["phases"]]
        periodicities = [int(value) for value in values["periodicities"]]
    elif {"k", "phase", "periodicity"} <= set(values):
        ks = [float(values["k"])]
        phases = [float(values["phase"])]
        periodicities = [int(values["periodicity"])]
    else:
        raise ValidationError("unsupported_proper_parameter_shape", term.external_id)
    if not ks or len(ks) != len(phases) or len(ks) != len(periodicities):
        raise ValidationError("invalid_proper_parameter_multiplicity", term.external_id)
    return ks, phases, periodicities


def _improper_values(term: AssignedForceTerm) -> tuple[float, float, int]:
    values = term.parameters
    if {"force_constant", "phase", "periodicity"} <= set(values):
        return float(values["force_constant"]), float(values["phase"]), int(values["periodicity"])
    if {"k", "phase", "periodicity"} <= set(values):
        return float(values["k"]), float(values["phase"]), int(values["periodicity"])
    raise ValidationError("unsupported_improper_parameter_shape", term.external_id)


def _graph_distances(adjacency: Mapping[str, set[str]], max_depth: int) -> dict[tuple[str, str], int]:
    distances: dict[tuple[str, str], int] = {}
    for start in sorted(adjacency):
        queue = deque([(start, 0)])
        seen = {start}
        while queue:
            current, depth = queue.popleft()
            if depth == max_depth:
                continue
            for neighbor in sorted(adjacency[current]):
                if neighbor in seen:
                    continue
                seen.add(neighbor)
                next_depth = depth + 1
                pair = tuple(sorted((start, neighbor)))
                distances[pair] = min(next_depth, distances.get(pair, next_depth))
                queue.append((neighbor, next_depth))
    return distances


def materialize_assigned_system(artifact: AssignedSystemArtifact) -> tuple[Any, dict[str, Any], dict[str, Any]]:
    """Return a real Xponge Molecule, external-ID map, and apply audit.

    This function is private to isolated worker processes.  Importing the
    module is harmless; Xponge force-field registries are populated only when
    this function is called.
    """
    # Importing this package registers the Amber exclusion/1-4 conventions and
    # serializer property types.  It happens only in this short-lived process.
    import Xponge.forcefield.amber  # noqa: F401
    from Xponge.forcefield.base.angle_base import AngleType
    from Xponge.forcefield.base.bond_base import BondType
    from Xponge.forcefield.base.dihedral_base import ImproperType, ProperType
    from Xponge.forcefield.base.lj_base import LJType
    from Xponge.forcefield.base.nb14_base import NB14Type
    from Xponge.helper import AtomType, Molecule, Residue, ResidueType

    topology = artifact.topology
    parameter_by_id = {item.external_id: item for item in artifact.atom_parameters}
    prepared_atom_by_id = {item.external_id: item for item in topology.atoms}
    atom_by_external_id: dict[str, Any] = {}
    lj_type_name_by_parameters: dict[tuple[float, float], str] = {}
    residue_id_by_atom_id = {
        atom_id: residue.external_id
        for residue in topology.residues
        for atom_id in residue.atom_ids
    }
    namespace = artifact.assigned_system_hash[:12]
    molecule = Molecule(f"metal_assignment_{namespace}")

    for residue_index, prepared_residue in enumerate(topology.residues):
        residue_type = ResidueType(name=f"MA_RES_{namespace}_{residue_index}")
        used_names: set[str] = set()
        for local_index, atom_id in enumerate(prepared_residue.atom_ids):
            prepared_atom = prepared_atom_by_id[atom_id]
            assigned = parameter_by_id[atom_id]
            lj_parameters = (float(assigned.epsilon), float(assigned.rmin))
            lj_type_name = lj_type_name_by_parameters.get(lj_parameters)
            if lj_type_name is None:
                lj_type_name = f"MA_LJ_{namespace}_{len(lj_type_name_by_parameters)}"
                LJType(
                    name=f"{lj_type_name}-{lj_type_name}",
                    epsilon=assigned.epsilon,
                    rmin=assigned.rmin,
                )
                lj_type_name_by_parameters[lj_parameters] = lj_type_name
            atom_type = AtomType(
                name=assigned.effective_type,
                charge=assigned.charge,
                mass=assigned.mass,
                LJtype=lj_type_name,
            )
            atom_name = prepared_atom.atom_name or f"A{local_index + 1}"
            if atom_name in used_names:
                atom_name = f"{atom_name}_{local_index + 1}"
            used_names.add(atom_name)
            residue_type.add_atom(atom_name, atom_type, *prepared_atom.coordinates)
        residue = Residue(residue_type, directly_copy=True)
        residue.name = prepared_residue.residue_name or prepared_residue.external_id
        molecule.add_residue(residue)
        for atom_id, atom in zip(prepared_residue.atom_ids, residue.atoms):
            coordinates = prepared_atom_by_id[atom_id].coordinates
            atom.x, atom.y, atom.z = map(float, coordinates)
            atom_by_external_id[atom_id] = atom

    active_edges = [edge for edge in (*topology.bonds, *topology.links) if edge.active]
    internal_bond_ids: set[str] = set()
    cross_edge_ids: set[str] = set()
    adjacency = {atom.external_id: set() for atom in topology.atoms}
    for edge in active_edges:
        atom1_id, atom2_id = edge.atom_ids
        adjacency[atom1_id].add(atom2_id)
        adjacency[atom2_id].add(atom1_id)
        same_residue = residue_id_by_atom_id[atom1_id] == residue_id_by_atom_id[atom2_id]
        if same_residue:
            if edge in topology.links:
                raise ValidationError("internal_prepared_link", edge.external_id)
            atom1, atom2 = atom_by_external_id[atom1_id], atom_by_external_id[atom2_id]
            atom1.residue.add_connectivity(atom1, atom2)
            internal_bond_ids.add(edge.external_id)
        else:
            molecule.add_residue_link(atom_by_external_id[atom1_id], atom_by_external_id[atom2_id])
            cross_edge_ids.add(edge.external_id)

    molecule.get_atoms()
    molecule.atom_index = {atom: index for index, atom in enumerate(molecule.atoms)}

    force_counts: Counter[str] = Counter()
    for term in artifact.force_terms:
        atoms = [atom_by_external_id[atom_id] for atom_id in term.atom_ids]
        if term.kind in {"bond", "distance_constraint"}:
            k, equilibrium = _bond_values(term)
            type_name = BondType.Get_Type_Name(atoms)
            force_type = BondType(name=type_name, k=k, b=equilibrium)
            force = BondType.entity(atoms, force_type, term.external_id)
            molecule.add_bonded_force(force)
            force_counts[term.kind] += 1
        elif term.kind == "angle":
            k, equilibrium = _angle_values(term)
            type_name = AngleType.Get_Type_Name(atoms)
            force_type = AngleType(name=type_name, k=k, b=equilibrium)
            force = AngleType.entity(atoms, force_type, term.external_id)
            molecule.add_bonded_force(force)
            force_counts[term.kind] += 1
        elif term.kind == "proper_dihedral":
            ks, phases, periodicities = _proper_values(term)
            type_name = ProperType.Get_Type_Name(atoms)
            force_type = ProperType(name=type_name, k=ks[0], phi0=phases[0], periodicity=periodicities[0])
            force_type.ks = ks
            force_type.phi0s = phases
            force_type.periodicitys = periodicities
            force_type.multiple_numbers = len(ks)
            force = ProperType.entity(atoms, force_type, term.external_id)
            molecule.add_bonded_force(force)
            force_counts[term.kind] += 1
        elif term.kind == "improper_dihedral":
            k, phase, periodicity = _improper_values(term)
            # ImproperType's symmetry registry parses four hyphen-separated
            # atom-type names to identify the central atom.  Assigned-system
            # effective types are unique, so this name is also collision-free.
            improper_type_name = "-".join(atom.type.name for atom in atoms)
            force_type = ImproperType(
                name=improper_type_name,
                k=k,
                phi0=phase,
                periodicity=periodicity,
            )
            force = ImproperType.entity(atoms, force_type, term.external_id)
            molecule.add_bonded_force(force)
            force_counts[term.kind] += 1
        else:
            raise ValidationError("unsupported_assigned_force_kind", term.kind, term.external_id)

    distances = _graph_distances(adjacency, artifact.force_protocol.exclusion_bond_depth)
    for (atom1_id, atom2_id), distance in distances.items():
        atom1, atom2 = atom_by_external_id[atom1_id], atom_by_external_id[atom2_id]
        atom1.link_atom(distance + 1, atom2)
        atom2.link_atom(distance + 1, atom1)

    one_four_pairs = [pair for pair, distance in distances.items() if distance == 3]
    for index, (atom1_id, atom2_id) in enumerate(sorted(one_four_pairs)):
        force_type = NB14Type(
            name=f"MA_NB14_{namespace}_{index}",
            kLJ=artifact.force_protocol.one_four_lj_scale,
            kee=artifact.force_protocol.one_four_electrostatic_scale,
        )
        molecule.add_bonded_force(NB14Type.entity(
            [atom_by_external_id[atom1_id], atom_by_external_id[atom2_id]],
            force_type,
            f"one_four:{atom1_id}:{atom2_id}",
        ))

    molecule.built = True
    expected_atom_ids = {atom.external_id for atom in topology.atoms}
    prepared_residue_order = [
        atom_id
        for residue in topology.residues
        for atom_id in residue.atom_ids
    ]
    if (
        len(prepared_residue_order) != len(expected_atom_ids)
        or len(prepared_residue_order) != len(set(prepared_residue_order))
        or set(prepared_residue_order) != expected_atom_ids
    ):
        raise ValidationError("prepared_residue_atom_bijection_mismatch", artifact.request_id)
    external_id_by_atom = {atom: atom_id for atom_id, atom in atom_by_external_id.items()}
    materialized_order = [external_id_by_atom[atom] for atom in molecule.atoms]
    if materialized_order != prepared_residue_order:
        raise ValidationError("xponge_atom_insertion_order_mismatch", artifact.request_id)
    if len(molecule.atoms) != len(expected_atom_ids) or len(atom_by_external_id) != len(expected_atom_ids):
        raise ValidationError("xponge_atom_materialization_mismatch", artifact.request_id)
    if len(molecule.residue_links) != len(cross_edge_ids):
        raise ValidationError("xponge_residue_link_materialization_mismatch", artifact.request_id)
    if internal_bond_ids | cross_edge_ids != {edge.external_id for edge in active_edges}:
        raise ValidationError("xponge_topology_materialization_mismatch", artifact.request_id)

    audit = {
        "topology_hash_before": topology.topology_hash,
        "topology_hash_after": topology.computed_hash(),
        "atom_count": len(molecule.atoms),
        "residue_count": len(molecule.residues),
        "component_count": len(topology.components),
        "prepared_bond_count": sum(1 for edge in topology.bonds if edge.active),
        "prepared_link_count": sum(1 for edge in topology.links if edge.active),
        "xponge_residue_link_count": len(molecule.residue_links),
        "materialized_force_counts": dict(sorted(force_counts.items())),
        "one_four_pair_count": len(one_four_pairs),
        "excluded_pair_count": len(distances),
        "effective_type_count": len({item.effective_type for item in artifact.atom_parameters}),
        "materialized_lj_type_count": len(lj_type_name_by_parameters),
        "partial_charge_sum": sum(item.charge for item in artifact.atom_parameters),
        "materialized_atom_order_hash": _sha256({"ordered_atom_ids": materialized_order}),
        "external_id_mapping_complete": True,
        "no_inferred_edges": True,
        "registry_isolation": "subprocess",
    }
    return molecule, atom_by_external_id, audit


def _materialize(artifact: AssignedSystemArtifact) -> dict[str, Any]:
    _, _, audit = materialize_assigned_system(artifact)
    return audit


def _execute(value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping) or set(value) != {"protocol_version", "assigned_system"}:
        raise ValidationError("invalid_parameter_apply_request", "unexpected fields")
    if value["protocol_version"] != APPLY_WORKER_PROTOCOL_VERSION:
        raise ValidationError("unsupported_parameter_apply_protocol", str(value["protocol_version"]))
    artifact = assigned_system_from_dict(value["assigned_system"])
    if artifact.complete:
        raise ValidationError("already_complete_assigned_system", artifact.request_id)
    audit = _materialize(artifact)
    response = {
        "ok": True,
        "protocol_version": APPLY_WORKER_PROTOCOL_VERSION,
        "assigned_system_hash": artifact.assigned_system_hash,
        "application_audit": audit,
    }
    response["response_hash"] = _sha256(response)
    return response


def main() -> None:
    try:
        response = _execute(json.load(sys.stdin))
    except Exception as exc:  # worker boundary must always be machine-readable
        response = {
            "ok": False,
            "error": {
                "code": getattr(exc, "code", "parameter_apply_failed"),
                "type": type(exc).__name__,
                "message": str(exc),
                "path": getattr(exc, "path", ""),
                "traceback": traceback.format_exc(limit=12),
            },
        }
    json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
    sys.stdout.write("\n")
    if response.get("ok") is not True:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
