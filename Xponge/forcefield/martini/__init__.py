"""
This **package** sets the basic configuration of OPLS force field
"""
import os
from ... import GlobalSetting, load_ffitp, AtomType, set_global_alternative_names
from ..base import charge_base, mass_base, lj_base, bond_base

lj_base.LJType.combining_method_A = lj_base.Lorentz_Berthelot_For_A
lj_base.LJType.combining_method_B = lj_base.Lorentz_Berthelot_For_B

def load_parameter_from_ffitp(filename, folder):
    """
    This **function** is used to get Martini force field parameters from GROMACS ffitp

    :param filename: the name of the input file
    :param prefix: the folder of the file
    :return: None
    """
    filename = os.path.join(folder, filename)
    output = load_ffitp(filename)

    AtomType.New_From_String(output["atomtypes"])
    lj_base.LJType.New_From_String(output["LJ"])
