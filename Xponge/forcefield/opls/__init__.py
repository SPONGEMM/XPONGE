"""
This **package** sets the basic configuration of OPLS force field
"""
import os
from ... import GlobalSetting, load_ffitp, AtomType, ResidueType, set_global_alternative_names
from ..base import charge_base, mass_base, lj_base, bond_base, angle_base, \
    dihedral_base, nb14_base, rb_dihedral_base, \
    virtual_atom_base, cmap_base, exclude_base

lj_base.LJType.combining_method_A = lj_base.Good_Hope_For_A
lj_base.LJType.combining_method_B = lj_base.Good_Hope_For_B

GlobalSetting.Set_Invisible_Bonded_Forces(["dihedral"])

OPLS_BOND_TYPE_MAP = {}

def opls_type_name(atoms):
    atom_names = [OPLS_BOND_TYPE_MAP.get(atom.type.name, atom.type.name) for atom in atoms]
    return "-".join(atom_names)

bond_base.BondType.Type_Name_Getter(opls_type_name)
angle_base.AngleType.Type_Name_Getter(opls_type_name)
dihedral_base.ImproperType.Type_Name_Getter(opls_type_name)
rb_dihedral_base.ProperType.Type_Name_Getter(opls_type_name)

exclude_base.Exclude(4)

def load_parameter_from_ffitp(filename, folder):
    """
    This **function** is used to get opls force field parameters from GROMACS ffitp

    :param filename: the name of the input file
    :param prefix: the folder of the file
    :return: None
    """
    filename = os.path.join(folder, filename)
    output = load_ffitp(filename)

    AtomType.New_From_String(output["atomtypes"])
    bond_base.BondType.New_From_String(output["bonds"])
    dihedral_base.ProperType.New_From_String(output["dihedrals"])
    lj_base.LJType.New_From_String(output["LJ"])
    angle_base.AngleType.New_From_String(output["angles"])
    dihedral_base.ImproperType.New_From_String(output["periodic_impropers"])
    rb_dihedral_base.ProperType.New_From_String(output["RB_dihedrals"])
    nb14_base.NB14Type.New_From_String(output["nb14"])
    OPLS_BOND_TYPE_MAP.update(output["bond_type_names"])
