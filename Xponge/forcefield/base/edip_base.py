"""
This **module** is the basic setting for the force field property of EDIP
"""
import numpy as np
from ...helper import Molecule, AtomType, Generate_New_Pairwise_Force_Type, Xdict
from itertools import product

AtomType.Add_Property({"EDIPType": str})

# pylint: disable=invalid-name
EDIPType = Generate_New_Pairwise_Force_Type("EDIP",
                                          {"A": float, "B": float, "a": float, "c": float, "alpha": float, "beta": float, "eta": float, "gamma": float,
                                           "l": float, "mu": float, "rho": float, "sigma": float, "Q0": float, "u1": float, "u2": float, "u3": float, "u4": float})

EDIPType.Set_Property_Unit("B", "distance", "A")
EDIPType.Set_Property_Unit("a", "distance", "A")
EDIPType.Set_Property_Unit("c", "distance", "A")
EDIPType.Set_Property_Unit("gamma", "distance", "A")
EDIPType.Set_Property_Unit("sigma", "distance", "A")
EDIPType.Set_Property_Unit("A", "energy", "kcal/mol")
EDIPType.Set_Property_Unit("l", "energy", "kcal/mol")

def _find_all_types(mol):
    edip_types = []
    edip_type_index = Xdict()
    atom_lines = []
    for atom in mol.atoms:
        if atom.EDIPType not in edip_type_index.keys():
            edip_type_index[atom.EDIPType] = len(edip_types)
            edip_types.append(atom.EDIPType)
        atom_lines.append(f"{edip_type_index[atom.EDIPType]}")
    return edip_types, edip_type_index, atom_lines

def _get_twobody_outputs(edip_types, edip_type_index):
    results = []
    for type1, type2 in product(edip_types, edip_types):
        edip_type = EDIPType.get_type(type1 + "-" + type2)
        results.append(f"{edip_type_index[type1]} {edip_type_index[type2]} {edip_type.alpha} {edip_type.c} " +
            f"{edip_type.a} {edip_type.A} {edip_type.B} {edip_type.rho} {edip_type.beta} {edip_type.sigma}")
    return results

def _get_threebody_outputs(edip_types, edip_type_index):
    results = []
    for type1, type2, type3 in product(edip_types, edip_types, edip_types):
        edip_type = EDIPType.get_type(type1 + "-" + type2 + "-" + type3)
        results.append(f"{edip_type_index[type1]} {edip_type_index[type2]} {edip_type_index[type3]} " +
            f"{edip_type.eta} {edip_type.gamma} {edip_type.l} {edip_type.Q0} {edip_type.mu} " +
            f"{edip_type.u1} {edip_type.u2} {edip_type.u3} {edip_type.u4}")
    return results

#pylint: disable=unused-argument
@Molecule.Set_Save_SPONGE_Input("EDIP")
def _write_edip(self):
    edip_types, edip_type_index, atom_lines = _find_all_types(self)
    edip_twobody_outputs = _get_twobody_outputs(edip_types, edip_type_index)
    edip_threebody_outputs = _get_threebody_outputs(edip_types, edip_type_index)
    output = f"{len(self.atoms)} {len(edip_types)}"
    output += "\n# type1 type2 alpha c[A] a[A] A[kcal/mol] B[A] rho beta sigma[A] (This is the first required comment line)\n"
    output += "\n".join(edip_twobody_outputs)
    output += "\n# type1 type2 type3 eta gamma[A] l[kcal/mol] Q0 mu u1 u2 u3 u4 (This is the second required comment line)\n"
    output += "\n".join(edip_threebody_outputs)
    output += "\n# atom type from the zeroth atom (This is the third required comment line)\n"
    output += "\n".join(atom_lines)
    return output

