"""
This **module** is the basic setting for the force field format of harmonic bond
"""
from ... import Generate_New_Bonded_Force_Type
from ...helper import Molecule, set_global_alternative_names, Xdict, GlobalSetting

# pylint: disable=invalid-name
BondType = Generate_New_Bonded_Force_Type("bond", "1-2", {"k": float, "b": float}, True)

BondType.Set_Property_Unit("k", "energy·distance^-2", "kcal/mol·A^-2")
BondType.Set_Property_Unit("b", "distance", "A")


@Molecule.Set_Save_SPONGE_Input("bond")
def write_bond(self):
    """
    This **function** is used to write SPONGE input file

    :param self: the Molecule instance
    :return: the string to write
    """
    bonds = []
    for bond in self.bonded_forces.get("bond", []):
        order = list(range(2))
        if bond.k != 0:
            if self.atom_index[bond.atoms[order[0]]] > self.atom_index[bond.atoms[order[-1]]]:
                temp_order = order[::-1]
            else:
                temp_order = order
            bonds.append("%d %d %f %f" % (self.atom_index[bond.atoms[temp_order[0]]]
                                          , self.atom_index[bond.atoms[temp_order[1]]], bond.k, bond.b))

    if bonds:
        towrite = "%d\n" % len(bonds)
        bonds.sort(key=lambda x: list(map(int, x.split()[:2])))
        towrite += "\n".join(bonds)

        return towrite
    return None


@Molecule.Set_MindSponge_Todo("bond")
def _do(self, sys_kwarg, ene_kwarg, use_pbc):
    """

    :return:
    """
    from mindsponge.potential import BondEnergy
    if "bond" not in sys_kwarg:
        sys_kwarg["bond"] = []
    if "bond" not in ene_kwarg:
        ene_kwarg["bond"] = Xdict()
        ene_kwarg["bond"]["function"] = lambda system, ene_kwarg: BondEnergy(
            index=system.bond, use_pbc=use_pbc,
            force_constant=ene_kwarg["bond"]["force_constant"],
            bond_length=ene_kwarg["bond"]["bond_length"],
            length_unit="A", energy_unit="kcal/mol")
        ene_kwarg["bond"]["force_constant"] = []
        ene_kwarg["bond"]["bond_length"] = []
    bonds = []
    force_constants = []
    bond_lengths = []
    for bond in self.bonded_forces.get("bond", []):
        if bond.k == 0:
            continue
        bonds.append([self.atom_index[bond.atoms[0]], self.atom_index[bond.atoms[1]]])
        force_constants.append(bond.k * 2)
        bond_lengths.append(bond.b)
    sys_kwarg["bond"].append(bonds)
    ene_kwarg["bond"]["force_constant"].append(force_constants)
    ene_kwarg["bond"]["bond_length"].append(bond_lengths)

@GlobalSetting.set_gmx_bonded_type_parser("bond", 1)
def _gmx_parser(words, mol, stat):
    """ parsing gmx """
    def _type_name(type_):
        return type_ if isinstance(type_, str) else type_.name

    atom1 = stat[int(words[0])]
    atom2 = stat[int(words[1])]
    if len(words) == 3:
        type1 = _type_name(atom1.type)
        type2 = _type_name(atom2.type)
        string = f"{type1}-{type2}"
        reversed_string = f"{type2}-{type1}"
        if string in BondType.get_all_types():
            mol.add_bonded_force(BondType.entity([atom1, atom2], BondType.getType(string)))
        elif reversed_string in BondType.get_all_types():
            mol.add_bonded_force(BondType.entity([atom1, atom2], BondType.getType(reversed_string)))
        else:
            raise KeyError(
                "Bonded Force Type "
                f"{string} (or {reversed_string}) not found. "
                "Did you import the proper force field / load bondtypes, "
                "or provide explicit bond parameters (5 columns) in the [ bonds ] section?"
            )
    elif len(words) == 5:
        new_force = BondType.entity([atom1, atom2], BondType.getType("UNKNOWNS"))
        new_force.k = float(words[3])
        new_force.b = float(words[4])
        mol.add_bonded_force(new_force)
    else:
        raise ValueError(f"Only 3 or 5 words should be in the line '{' '.join(words)}'")

set_global_alternative_names()
