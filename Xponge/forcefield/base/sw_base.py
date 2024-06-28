"""
This **module** is the basic setting for the force field property of Stillinger-Weber
"""
import numpy as np
from ...helper import Molecule, AtomType, Generate_New_Pairwise_Force_Type, Xdict
from itertools import product

AtomType.Add_Property({"SWType": str})

# pylint: disable=invalid-name
SWType = Generate_New_Pairwise_Force_Type("SW",
                                          {"A": float, "B": float, "epsilon": float, "p": float, "q": float, "a": float, "gamma": float, "sigma": float,
                                           "l": float, "b": float})

SWType.Set_Property_Unit("sigma", "distance", "A")
SWType.Set_Property_Unit("epsilon", "energy", "kcal/mol")

def _find_all_types(mol):
    sw_types = []
    sw_type_index = Xdict()
    atom_lines = []
    for atom in mol.atoms:
        if atom.SWType not in sw_type_index.keys():
            sw_type_index[atom.SWType] = len(sw_types)
            sw_types.append(atom.SWType)
        atom_lines.append(f"{sw_type_index[atom.SWType]}")
    return sw_types, sw_type_index, atom_lines

def _get_twobody_outputs(sw_types, sw_type_index):
    results = []
    for type1, type2 in product(sw_types, sw_types):
        sw_type = SWType.get_type(type1 + "-" + type2)
        results.append(f"{sw_type_index[type1]} {sw_type_index[type2]} {sw_type.A} {sw_type.B} " +
            f"{sw_type.epsilon} {sw_type.p} {sw_type.q} {sw_type.a} {sw_type.gamma} {sw_type.sigma}")
    return results

def _get_threebody_outputs(sw_types, sw_type_index):
    results = []
    for type1, type2, type3 in product(sw_types, sw_types, sw_types):
        sw_type = SWType.get_type(type1 + "-" + type2 + "-" + type3)
        results.append(f"{sw_type_index[type1]} {sw_type_index[type2]} {sw_type_index[type3]} " +
            f"{sw_type.l} {sw_type.epsilon} {sw_type.b}")
    return results

#pylint: disable=unused-argument
@Molecule.Set_Save_SPONGE_Input("SW")
def _write_sw(self):
    sw_types, sw_type_index, atom_lines = _find_all_types(self)
    sw_twobody_outputs = _get_twobody_outputs(sw_types, sw_type_index)
    sw_threebody_outputs = _get_threebody_outputs(sw_types, sw_type_index)
    output = f"{len(self.atoms)} {len(sw_types)}"
    output += "\n# type1 type2 A B epsilon[kcal/mol] p q a gamma sigma[Angstrom] (This is the first required comment line)\n"
    output += "\n".join(sw_twobody_outputs)
    output += "\n# type1 type2 type3 lambda epsilon[kcal/mol] b (This is the second required comment line)\n"
    output += "\n".join(sw_threebody_outputs)
    output += "\n# atom type from the zeroth atom (This is the third required comment line)\n"
    output += "\n".join(atom_lines)
    return output

