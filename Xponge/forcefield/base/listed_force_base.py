"""
This **module** is the basic setting for the force field property of general listed forces
"""
from ...helper import Molecule, set_global_alternative_names

LISTED_FORCE_DEFINITION = {}

@Molecule.Set_Save_SPONGE_Input("listed_forces")
def write_listed_forces(self):
    """
    This **function** is used to write SPONGE input file

    :param self: the Molecule instance
    :return: the string to write
    """
    towrite = ""
    for frcs in self.bonded_forces.values():
        if not frcs or frcs[0].get_class_name() not in LISTED_FORCE_DEFINITION:
            continue
        towrite += LISTED_FORCE_DEFINITION[frcs[0].get_class_name()]
    if towrite:
        return towrite
    return None

set_global_alternative_names()
