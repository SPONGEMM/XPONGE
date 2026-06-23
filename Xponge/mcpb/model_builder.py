"""Local model extraction helpers for MCPB workflows."""

from __future__ import annotations

import time

from .. import Molecule
from ._compat import atom_id, atom_to_residue_map, get_residue_link_by_atom_ids, refresh_molecule_indices
from .models import MCPBLocalModel


def _collect_parent_connect_pairs(molecule) -> list[tuple[int, int]]:
    refresh_molecule_indices(molecule)
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    for link in molecule.residue_links:
        atom1 = atom_id(link.atom1, molecule)
        atom2 = atom_id(link.atom2, molecule)
        pair = (atom1, atom2) if atom1 < atom2 else (atom2, atom1)
        if pair in seen:
            continue
        seen.add(pair)
        pairs.append(pair)
    return pairs


def build_local_model(
    molecule,
    residue_ids: tuple[int, ...],
    *,
    name: str,
    extra_connect_pairs: tuple[tuple[int, int], ...] = (),
) -> MCPBLocalModel:
    refresh_molecule_indices(molecule)
    model = Molecule(name)
    source_atom_ids: list[int] = []
    atom_id_map: dict[int, int] = {}
    residue_id_map: dict[int, int] = {}
    forcopy = ("mcpb_local_model", name, tuple(residue_ids), time.time_ns())

    for residue_id in residue_ids:
        source_residue = molecule.residues[int(residue_id)]
        copied_residue = source_residue.deepcopy(forcopy)
        model.add_residue(copied_residue)
        residue_id_map[int(residue_id)] = len(model.residues) - 1

    refresh_molecule_indices(model)
    for residue_id in residue_ids:
        source_residue = molecule.residues[int(residue_id)]
        for source_atom in source_residue.atoms:
            copied_atom = source_atom.copied[forcopy]
            source_index = atom_id(source_atom, molecule)
            mapped_atom_id = int(model.atom_index[copied_atom])
            atom_id_map[source_index] = mapped_atom_id
            source_atom_ids.append(source_index)

    connect_pairs = _collect_parent_connect_pairs(molecule)
    connect_pairs.extend(tuple(pair) for pair in extra_connect_pairs)
    seen_connect: set[tuple[int, int]] = set()
    for atom1, atom2 in connect_pairs:
        if atom1 not in atom_id_map or atom2 not in atom_id_map:
            continue
        mapped = (atom_id_map[atom1], atom_id_map[atom2])
        mapped = mapped if mapped[0] < mapped[1] else (mapped[1], mapped[0])
        if mapped in seen_connect:
            continue
        seen_connect.add(mapped)
        if get_residue_link_by_atom_ids(model, mapped[0], mapped[1]) is None:
            model.add_residue_link(model.atoms[mapped[0]], model.atoms[mapped[1]])

    for residue_id in residue_ids:
        for source_atom in molecule.residues[int(residue_id)].atoms:
            source_atom.copied.pop(forcopy, None)

    return MCPBLocalModel(
        molecule=model,
        source_atom_ids=tuple(source_atom_ids),
        atom_id_map=atom_id_map,
        residue_id_map=residue_id_map,
    )


def build_small_and_large_models(request, selection) -> tuple[MCPBLocalModel, MCPBLocalModel]:
    source_atom_to_residue = atom_to_residue_map(request.molecule)
    ion_residue_ids = {source_atom_to_residue[atom_id] for atom_id in selection.ion_atom_ids}
    coordinating_residue_ids = {source_atom_to_residue[atom_id] for atom_id in selection.coordinating_atom_ids}
    small_residue_ids = tuple(sorted(ion_residue_ids | coordinating_residue_ids))
    large_residue_ids = tuple(sorted(set(selection.selected_residue_ids) | set(small_residue_ids)))
    small_model = build_local_model(
        request.molecule,
        small_residue_ids,
        name="MCPB_SMALL",
        extra_connect_pairs=selection.bonded_pairs,
    )
    large_model = build_local_model(
        request.molecule,
        large_residue_ids,
        name="MCPB_LARGE",
        extra_connect_pairs=selection.bonded_pairs,
    )
    return small_model, large_model
