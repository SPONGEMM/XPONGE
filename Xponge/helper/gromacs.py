"""
    This **module** contains functions to interact with gromacs
"""
from .. import Molecule, set_global_alternative_names
from . import Xdict

def _check_new_residue_when_sorting(mol, mol_res_id,
  sort_map, cri, ri, rname):
    """ check whether it is a new residue """
    if ri != cri:
        cri = ri
        if mol_res_id >= 0:
            residue = mol.residues[mol_res_id]
            restype = residue.type
            if len(restype.name) == 4 and restype.name[0] == "C":
                sort_map["O"] = sort_map.pop("OT1")
                sort_map["OXT"] = sort_map.pop("OT2")
            for atom in residue.atoms:
                if atom.name.startswith("H") or atom.name in ("1H", "2H", "3H"):
                    atom_connected = next(iter(restype.connectivity[restype.name2atom(atom.name)]))
                    find_name = atom_connected.name if atom_connected.name in sort_map else atom_connected.name[:-1]
                    sort_map[atom.name] = sort_map[find_name] + 0.5
            residue.atoms.sort(key=lambda atom: sort_map[atom.name] if atom.name in sort_map else sort_map[atom.name[:-1]])
        sort_map.clear()
        sort_map.not_found_message = "The name of the atom '{}' in the %d-th residue in \
Xponge.Molecule can not be found in the gro file (%5d%-5s)"%(mol_res_id + 1, ri, rname)
        mol_res_id += 1
    return mol_res_id, cri

def sort_atoms_by_gro(mol, gro):
    """
       This **function** sorts the atoms in a Xponge.Molecule according to the index in a gro file

       :param mol: a Xponge.Molecule
       :param gro: the gro file
    """
    if not isinstance(mol, Molecule):
        raise TypeError("The input for sorting should be an Xponge.Molecule")
    mol.get_atoms()
    current_residue_index = 0
    mol_res_id = -1
    current_sort_map = Xdict()
    with open(gro) as f:
        for li, line in enumerate(f):
            if li == 1:
                atom_numbers = int(line)
                if atom_numbers != len(mol.atoms):
                    raise ValueError("The number of atoms is not equal \
in the gro file and in the Xponge.Molecule instance")
            elif 1 < li < len(mol.atoms) + 2:
                residue_index = int(line[:5])
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                if residue_name in ("WAT", "H2O", "HOH", "SOL"):
                    break
                mol_res_id, current_residue_index = _check_new_residue_when_sorting(mol, mol_res_id,
                   current_sort_map, current_residue_index, residue_index, residue_name)
                if atom_name.startswith("H") or atom_name in ("1H", "2H", "3H"):
                    continue
                current_sort_map[atom_name] = len(current_sort_map)
    _check_new_residue_when_sorting(mol, mol_res_id,
       current_sort_map, current_residue_index, current_residue_index + 1, "")


set_global_alternative_names()
