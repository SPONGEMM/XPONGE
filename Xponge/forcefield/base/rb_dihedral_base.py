"""
This **module** is the basic setting for the force field format of Ryckaert-Bellemans dihedral
"""
from ... import Generate_New_Bonded_Force_Type
from ...helper import Molecule, set_global_alternative_names
from .listed_force_base import LISTED_FORCE_DEFINITION

# pylint: disable=invalid-name
ProperType = Generate_New_Bonded_Force_Type("Ryckaert_Bellemans", "1-2-3-4", {"c0": float, "c1": float, "c2": float,
                                                                              "c3": float, "c4": float, "c5": float},
                                            True)

ProperType.Set_Property_Unit("c0", "energy", "kcal/mol")
ProperType.Set_Property_Unit("c1", "energy", "kcal/mol")
ProperType.Set_Property_Unit("c2", "energy", "kcal/mol")
ProperType.Set_Property_Unit("c3", "energy", "kcal/mol")
ProperType.Set_Property_Unit("c4", "energy", "kcal/mol")
ProperType.Set_Property_Unit("c5", "energy", "kcal/mol")

LISTED_FORCE_DEFINITION["Ryckaert_Bellemans"] = """[[[ Ryckaert_Bellemans ]]]
[[ parameters ]]
int atom_a, int atom_b, int atom_c, int atom_d, float c0, float c1, float c2, float c3, float c4, float c5
[[ potential ]]
SADfloat<15> cphi = cosf(phi_abcd - CONSTANT_Pi);
SADfloat<15> cphi2 = cphi * cphi;
SADfloat<15> cphi3 = cphi2 * cphi;
SADfloat<15> cphi4 = cphi3 * cphi;
SADfloat<15> cphi5 = cphi4 * cphi;
E = c0 + c1 * cphi + c2 * cphi2 + c3 * cphi3 + c4 * cphi4 + c5 * cphi5;
[[ end ]]
"""

@Molecule.Set_Save_SPONGE_Input("Ryckaert_Bellemans")
def write_dihedral(self):
    """
    This **function** is used to write SPONGE input file

    :param self: the Molecule instance
    :return: the string to write
    """
    dihedrals = []
    for dihedral in self.bonded_forces.get("Ryckaert_Bellemans", []):
        order = list(range(4))
        if self.atom_index[dihedral.atoms[order[0]]] > self.atom_index[dihedral.atoms[order[-1]]]:
            temp_order = order[::-1]
        else:
            temp_order = order
        if dihedral.c0 != 0 or dihedral.c1 != 0 or dihedral.c2 != 0 or \
           dihedral.c3 != 0 or dihedral.c4 != 0 or dihedral.c5 != 0:
            dihedrals.append("%d %d %d %d %f %f %f %f %f %f" % (self.atom_index[dihedral.atoms[temp_order[0]]],
                                                                self.atom_index[dihedral.atoms[temp_order[1]]],
                                                                self.atom_index[dihedral.atoms[temp_order[2]]],
                                                                self.atom_index[dihedral.atoms[temp_order[3]]],
                                                                dihedral.c0, dihedral.c1,
                                                                dihedral.c2, dihedral.c3,
                                                                dihedral.c4, dihedral.c5))

    if dihedrals:
        towrite = "%d\n" % len(dihedrals)
        dihedrals.sort(key=lambda x: list(map(float, x.split())))
        towrite += "\n".join(dihedrals)

        return towrite
    return None


set_global_alternative_names()
