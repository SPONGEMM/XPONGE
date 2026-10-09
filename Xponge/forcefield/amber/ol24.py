"""
This **module** set the basic configuration for OL24 (DNA)
"""
from ...helper import source, Xprint

source("....")
amber = source("...amber")

amber.load_parameters_from_parmdat("parm10.dat")
amber.load_parameters_from_frcmod("OL24.frcmod")

load_mol2(os.path.join(AMBER_DATA_DIR, "ol24.mol2"), as_template=True)

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

from ._terminal_monophosphate import register_terminal_monophosphate
register_terminal_monophosphate('dna', ('DA', 'DC', 'DG', 'DT'))

Xprint("""Reference for OL24:
  Refinement of the Sugar Puckering Torsion Potential in the AMBER DNA Force Field.
  DOI: 10.1021/acs.jctc.4c01100
""")
