"""
This **package** sets the basic configuration of Martini force field
"""
import os
from ... import GlobalSetting, load_ffitp, AtomType, set_global_alternative_names
from ..base import charge_base, mass_base, lj_base, bond_base

lj_base.LJType.combining_method_A = lj_base.Lorentz_Berthelot_For_A
lj_base.LJType.combining_method_B = lj_base.Lorentz_Berthelot_For_B

MARTINI_DATA_DIR = os.path.join(os.path.dirname(__file__), "martini_v300")

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


for _itp in (
    "martini_v3.0.0.itp",
    "martini_v3.0.0_ions_v1.itp",
    "martini_v3.0.0_solvents_v1.itp",
    "martini_v3.0.0_small_molecules_v1.itp",
    "martini_v3.0.0_sugars_v1.itp",
    "martini_v3.0.0_phospholipids_v1.itp",
    "martini_v3.0.0_nucleobases_v1.itp",
):
    load_parameter_from_ffitp(_itp, MARTINI_DATA_DIR)
