"""Local charge-refit helpers for MCPB workflows."""

from __future__ import annotations

from .. import Assign
from ..assign import resp as resp_module
from ._compat import atom_element, atom_id, atom_to_residue_map, refresh_molecule_indices
from .selection import infer_element_symbol

_RESP_DEFAULT_RADII = {
    "H": 1.2,
    "C": 1.5,
    "N": 1.5,
    "O": 1.4,
    "P": 1.8,
    "S": 1.75,
    "F": 1.35,
    "Cl": 1.7,
    "Br": 2.3,
}

_MCPB_METAL_RESP_RADII = {
    "Li": 1.82,
    "Na": 2.27,
    "K": 2.75,
    "Mg": 1.73,
    "Ca": 2.31,
    "Mn": 1.97,
    "Fe": 1.94,
    "Co": 1.92,
    "Ni": 1.63,
    "Cu": 1.40,
    "Zn": 1.39,
}


def _local_total_charge(request, local_model) -> int:
    refresh_molecule_indices(request.molecule)
    base = sum(float(request.molecule.atoms[source_atom_id].charge) for source_atom_id in local_model.source_atom_ids)
    charge = base
    for info in request.ion_info:
        if info.formal_charge is None or info.atom_id not in local_model.atom_id_map:
            continue
        charge += float(info.formal_charge) - float(request.molecule.atoms[info.atom_id].charge)
    return int(round(charge))


def _local_spin(request, local_model) -> int:
    spin = 0
    for info in request.ion_info:
        if info.spin is None or info.atom_id not in local_model.atom_id_map:
            continue
        spin += int(info.spin)
    return int(spin)


def _resp_radius_overrides(request, local_model) -> dict[str, float]:
    radius = {}
    ion_info_by_atom = {info.atom_id: info for info in request.ion_info}
    source_atom_to_residue = atom_to_residue_map(request.molecule)
    refresh_molecule_indices(request.molecule)
    for source_atom_id in local_model.source_atom_ids:
        atom = request.molecule.atoms[source_atom_id]
        residue = request.molecule.residues[source_atom_to_residue[source_atom_id]]
        info = ion_info_by_atom.get(source_atom_id)
        element = info.element if info is not None else infer_element_symbol(atom_element(atom), atom.name, residue.name)
        if element in _RESP_DEFAULT_RADII:
            continue
        if info is not None and "resp_radius" in info.metadata:
            radius[element] = float(info.metadata["resp_radius"])
            continue
        if element in _MCPB_METAL_RESP_RADII:
            radius[element] = _MCPB_METAL_RESP_RADII[element]
            continue
        raise ValueError(
            f"MCPB local RESP requires an explicit resp_radius for element {element}. "
            "Pass it in ion_info[i]['metadata']['resp_radius']."
        )
    return radius


def build_local_resp_assignment(request, local_model):
    source_atom_to_residue = atom_to_residue_map(request.molecule)
    ion_info_by_atom = {info.atom_id: info for info in request.ion_info}
    refresh_molecule_indices(request.molecule)
    refresh_molecule_indices(local_model.molecule)
    source_by_local = {local_id: source_id for source_id, local_id in local_model.atom_id_map.items()}

    assignment = Assign("MCPB")
    for local_atom_id, atom in enumerate(local_model.molecule.atoms):
        source_atom_id = source_by_local[local_atom_id]
        source_atom = request.molecule.atoms[source_atom_id]
        residue = request.molecule.residues[source_atom_to_residue[source_atom_id]]
        info = ion_info_by_atom.get(source_atom_id)
        element = info.element if info is not None else infer_element_symbol(
            atom_element(source_atom),
            source_atom.name,
            residue.name,
        )
        assignment.add_atom(
            element,
            float(atom.x),
            float(atom.y),
            float(atom.z),
            atom.name,
            float(atom.charge),
        )
        if info is not None and info.formal_charge is not None:
            assignment.formal_charge[local_atom_id] = int(info.formal_charge)

    for residue in local_model.molecule.residues:
        for template_atom, bonded_atoms in residue.type.connectivity.items():
            atom1 = residue.atoms[residue.type.atom2index(template_atom)]
            atom1_id = int(local_model.molecule.atom_index[atom1])
            for template_neighbor in bonded_atoms:
                atom2 = residue.atoms[residue.type.atom2index(template_neighbor)]
                atom2_id = int(local_model.molecule.atom_index[atom2])
                if atom1_id < atom2_id:
                    assignment.add_bond(atom1_id, atom2_id, 1)
    for link in local_model.molecule.residue_links:
        atom1_id = atom_id(link.atom1, local_model.molecule)
        atom2_id = atom_id(link.atom2, local_model.molecule)
        assignment.add_bond(atom1_id, atom2_id, 1)
    return assignment


def run_local_charge_refit(request, local_model):
    assignment = build_local_resp_assignment(request, local_model)
    radius = _resp_radius_overrides(request, local_model)
    total_charge = _local_total_charge(request, local_model)
    spin = _local_spin(request, local_model)
    charges = resp_module.resp_fit(
        assignment,
        basis=request.basis or "6-31g*",
        charge=total_charge,
        spin=spin,
        grid_density=1,
        grid_cell_layer=1,
        two_stage=False,
        radius=radius,
        backend=request.qm_backend,
        core="python",
    )
    refresh_molecule_indices(request.molecule)
    updated_atom_ids: list[int] = []
    for source_atom_id, local_atom_id in local_model.atom_id_map.items():
        request.molecule.atoms[source_atom_id].charge = float(charges[local_atom_id])
        updated_atom_ids.append(int(source_atom_id))
    return {
        "charges": list(charges),
        "updated_atom_ids": updated_atom_ids,
        "assignment": assignment,
        "total_charge": total_charge,
        "spin": spin,
        "radius": radius,
    }
