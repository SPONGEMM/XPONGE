"""
This **package** sets the DNA configuration of charmm27 force field
"""
from ....helper import source
from . import CHARMM27_DATA_DIR

source(".....")

load_mol2(os.path.join(CHARMM27_DATA_DIR, "dna.mol2"), as_template=True)

for i in "ATCG":
    i = "D" + i
    res = ResidueType.get_type(i)
    res5 = ResidueType.get_type(i + "5")
    res3 = ResidueType.get_type(i + "3")
    res.tail = "O3'"
    res.tail_next = "C3'"
    res5.tail = "O3'"
    res5.tail_next = "C3'"

    res.head = "P"
    res.head_next = "OP2"
    res3.head = "P"
    res3.head_next = "OP2"

    res.tail_link_conditions.append({"atoms": ["C3'", "O3'"], "parameter": 120 / 180 * np.pi})
    res.tail_link_conditions.append({"atoms": ["H3'", "C3'", "O3'"], "parameter": -54 / 180 * np.pi})
    res5.tail_link_conditions.append({"atoms": ["C3'", "O3'"], "parameter": 120 / 180 * np.pi})
    res5.tail_link_conditions.append({"atoms": ["H3'", "C3'", "O3'"], "parameter": -54 / 180 * np.pi})

    res.head_link_conditions.append({"atoms": ["OP2", "P"], "parameter": 108 / 180 * np.pi})
    res.head_link_conditions.append({"atoms": ["O5'", "OP2", "P"], "parameter": 113 / 180 * np.pi})
    res3.head_link_conditions.append({"atoms": ["OP2", "P"], "parameter": 108 / 180 * np.pi})
    res3.head_link_conditions.append({"atoms": ["O5'", "OP2", "P"], "parameter": 113 / 180 * np.pi})

    GlobalSetting.Add_PDB_Residue_Name_Mapping("head", res.name, res5.name)
    GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", res.name, res3.name)

# pylint:disable=undefined-variable
