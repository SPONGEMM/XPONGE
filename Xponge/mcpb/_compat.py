"""Compatibility helpers for MCPB on top of Xponge-origin data structures."""

from __future__ import annotations

from pathlib import Path


def refresh_molecule_indices(molecule):
    molecule.atoms = []
    for residue in molecule.residues:
        molecule.atoms.extend(residue.atoms)
    molecule.atom_index.clear()
    molecule.atom_index.update({atom: i for i, atom in enumerate(molecule.atoms)})
    return molecule


def atom_count(molecule) -> int:
    refresh_molecule_indices(molecule)
    return len(molecule.atoms)


def residue_count(molecule) -> int:
    return len(molecule.residues)


def atom_id(atom, molecule) -> int:
    refresh_molecule_indices(molecule)
    return int(molecule.atom_index[atom])


def residue_id(residue, molecule) -> int:
    return int(molecule.residues.index(residue))


def atom_to_residue_map(molecule) -> dict[int, int]:
    refresh_molecule_indices(molecule)
    mapping = {}
    for residue in molecule.residues:
        res_id = residue_id(residue, molecule)
        for atom in residue.atoms:
            mapping[atom_id(atom, molecule)] = res_id
    return mapping


def atom_ids_for_residue(residue, molecule) -> tuple[int, ...]:
    refresh_molecule_indices(molecule)
    return tuple(atom_id(atom, molecule) for atom in residue.atoms)


def ensure_atom(atom_or_id, molecule):
    refresh_molecule_indices(molecule)
    if isinstance(atom_or_id, int):
        return molecule.atoms[int(atom_or_id)]
    return atom_or_id


def atom_type_name(atom) -> str:
    atom_type = getattr(atom, "type", None)
    return getattr(atom_type, "name", str(atom_type))


def atom_element(atom) -> str:
    return str(getattr(atom, "element", ""))


def add_residue_link_by_atom_ids(molecule, atom1_id: int, atom2_id: int):
    refresh_molecule_indices(molecule)
    molecule.add_residue_link(molecule.atoms[int(atom1_id)], molecule.atoms[int(atom2_id)])


def get_residue_link_by_atom_ids(molecule, atom1_id: int, atom2_id: int):
    refresh_molecule_indices(molecule)
    return molecule.get_residue_link(molecule.atoms[int(atom1_id)], molecule.atoms[int(atom2_id)])


def save_pdb(molecule, filename) -> str:
    from .. import save_pdb as _save_pdb

    path = Path(filename)
    _save_pdb(molecule, str(path))
    return str(path)


def save_sponge_input(molecule, prefix: str, dirname) -> None:
    from .. import save_sponge_input as _save_sponge_input

    _save_sponge_input(molecule, prefix, dirname=str(dirname))
