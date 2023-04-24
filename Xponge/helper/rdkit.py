"""
This **module** gives the interface to the package RDKit
"""
try:
    from rdkit import Chem
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "'rdkit' package needed. Maybe you need 'conda install -c rdkit rdkit'") from exc

from . import Xdict, Xprint
from .namespace import set_global_alternative_names


def assign_to_rdmol(assign, ignore_bond_type=False):
    """
    This **function** is used to convert an Xponge.Assign to a RDKit.rdMol

    :param assign: the Assign instance
    :param ignore_bond_type: set the bond type always to UNSPECIFIED
    :return: the RDKit.rdMol instance
    """
    mol_a = Chem.RWMol()
    for atom in assign.atoms:
        mol_a.AddAtom(Chem.Atom(atom))
    for atom, bonds in assign.bonds.items():
        for aton, n in bonds.items():
            if aton < atom:
                continue
            if ignore_bond_type or n == -1:
                temp_bond = Chem.BondType.UNSPECIFIED
            elif n == 1:
                temp_bond = Chem.BondType.SINGLE
            elif n == 2:
                temp_bond = Chem.BondType.DOUBLE
            elif n == 3:
                temp_bond = Chem.BondType.TRIPLE
            else:
                raise NotImplementedError
            mol_a.AddBond(atom, aton, temp_bond)
    conf = Chem.Conformer(assign.atom_numbers)
    for i in range(assign.atom_numbers):
        conf.SetAtomPosition(i, assign.coordinate[i])
    mol = mol_a.GetMol()
    mol.AddConformer(conf)
    for i, atom in enumerate(mol.GetAtoms()):
        atom.SetFormalCharge(assign.formal_charge[i])
    flags = Chem.rdmolops.SanitizeFlags
    Chem.SanitizeMol(mol, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
    return mol


def rdmol_to_assign(rdmol):
    """
    This **function** is used to convert a RDKit.rdMol to an Xponge.Assign

    :param rdmol: the RDKit.rdMol instance
    :return: the Xponge.Assign instance
    """
    from .. import Assign
    assign = Assign()

    for atom in rdmol.GetAtoms():
        x, y, z = rdmol.GetConformer().GetAtomPosition(atom.GetIdx())
        assign.add_atom(atom.GetSymbol(), x=x, y=y, z=z)
        assign.formal_charge[-1] = atom.GetFormalCharge()

    has_unknown_bond = False
    for bond in rdmol.GetBonds():
        temp_bond = bond.GetBondType()
        if temp_bond == Chem.BondType.UNSPECIFIED:
            temp_bond = -1
        elif temp_bond == Chem.BondType.SINGLE:
            temp_bond = 1
        elif temp_bond == Chem.BondType.DOUBLE:
            temp_bond = 2
        elif temp_bond == Chem.BondType.TRIPLE:
            temp_bond = 3
        elif temp_bond == Chem.BondType.AROMATIC:
            temp_bond = -1
            has_unknown_bond = True
        else:
            Xprint(f"Unknown bond type {temp_bond}", "ERROR")
            raise NotImplementedError
        assign.add_bond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), temp_bond)
    assign.determine_ring_and_bond_type()
    if has_unknown_bond:
        assign.determine_bond_order()

    return assign

def insert_atom_type_to_rdmol(mol, res, assign, atom_type_dict=None):
    """
    This **function** inserts the atom types in the force field to the RDKit.rdmol instance.
    This is done by setting the isotope of the atom, so it may not be compatible with the other packages.

    :param mol: the RDKit.rdmol instance
    :param res: the Residue instance corresponding to ``mol``
    :param assign: the Assign instance corresponding to ``mol``
    :param atom_type_dict: the dict mapping the atom type to the isotope number
    :return:
    """
    if atom_type_dict is None:
        atom_type_dict = Xdict()
    for i, a in enumerate(mol.GetAtoms()):
        atom_type = res.name2atom(assign.names[i]).type.name
        if atom_type not in atom_type_dict:
            atom_type_dict[atom_type] = len(atom_type_dict)
        a.SetIsotope(atom_type_dict[atom_type])


def find_equal_atoms(assign):
    """
    This **function** is used to find the chemical equivalent atoms in the molecule

    :param assign: the Assign instance
    :return: a list of equalvalent atom index lists
    """
    mols = []
    canon_smiles = []
    for i in range(len(assign.atoms)):
        mols.append(assign_to_rdmol(assign))
        mols[-1].GetAtoms()[i].SetIsotope(1)
        canon_smiles.append(Chem.MolToSmiles(mols[-1], isomericSmiles=True))
    group = {i: i for i in range(len(assign.atoms))}
    for i in range(len(assign.atoms)):
        if group[i] == i:
            for j in range(i + 1, len(assign.atoms)):
                if canon_smiles[j] == canon_smiles[i]:
                    group[j] = i
    ret = []
    realmap = Xdict()
    for i in group:
        if group[i] == i:
            ret.append([i])
            realmap[i] = len(realmap)
        else:
            ret[realmap[group[i]]].append(i)
    return list(filter(lambda x: len(x) > 1, ret))


set_global_alternative_names()
