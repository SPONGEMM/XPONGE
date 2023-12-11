"""
This **package** sets the protein configuration of charmm36 force field
"""
from ....helper import source, set_real_global_variable
from . import CHARMM36_DATA_DIR

source(".....")

load_mol2(os.path.join(CHARMM36_DATA_DIR, "protein.mol2"), as_template=True)

ResidueType.set_type("HIS", ResidueType.get_type("HIE"))
ResidueType.set_type("NHIS", ResidueType.get_type("NHIE"))
ResidueType.set_type("CHIS", ResidueType.get_type("CHIE"))

set_real_global_variable("HIS", ResidueType.get_type("HIS"))
set_real_global_variable("NHIS", ResidueType.get_type("NHIS"))
set_real_global_variable("CHIS", ResidueType.get_type("CHIS"))

residues = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()

for resname in residues:
    res = ResidueType.get_type(resname)
    nres = ResidueType.get_type("N" + resname)
    cres = ResidueType.get_type("C" + resname)
    res.head, cres.head = "N", "N"
    res.head_length, cres.head_length = 1.3, 1.3
    res.head_next, cres.head_next = "CA", "CA"
    res.tail, nres.tail = "C", "C"
    res.tail_next, nres.tail_next = "CA", "CA"
    res.tail_length, nres.tail_length = 1.3, 1.3

    res.head_link_conditions.append({"atoms": ["CA", "N"], "parameter": 120 / 180 * np.pi})
    cres.head_link_conditions.append({"atoms": ["CA", "N"], "parameter": 120 / 180 * np.pi})

    if resname != "PRO":
        res.head_link_conditions.append({"atoms": ["H", "CA", "N"], "parameter": -np.pi})
        cres.head_link_conditions.append({"atoms": ["H", "CA", "N"], "parameter": -np.pi})
    else:
        res.head_link_conditions.append({"atoms": ["HA", "CA", "N"], "parameter": 0})
        cres.head_link_conditions.append({"atoms": ["HA", "CA", "N"], "parameter": 0})

    res.tail_link_conditions.append({"atoms": ["CA", "C"], "parameter": 120 / 180 * np.pi})
    res.tail_link_conditions.append({"atoms": ["O", "CA", "C"], "parameter": -np.pi})
    nres.tail_link_conditions.append({"atoms": ["CA", "C"], "parameter": 120 / 180 * np.pi})
    nres.tail_link_conditions.append({"atoms": ["O", "CA", "C"], "parameter": -np.pi})

    GlobalSetting.Add_PDB_Residue_Name_Mapping("head", resname, "N" + resname)
    GlobalSetting.Add_PDB_Residue_Name_Mapping("tail", resname, "C" + resname)

ResidueType.get_type("ACE").tail_next = "CH3"
ResidueType.get_type("ACE").tail_length = 1.3
ResidueType.get_type("ACE").tail_link_conditions.append({"atoms": ["CH3", "C"], "parameter": 120 / 180 * np.pi})
ResidueType.get_type("ACE").tail_link_conditions.append({"atoms": ["O", "CH3", "C"], "parameter": -np.pi})

ResidueType.get_type("NME").head_next = "CH3"
ResidueType.get_type("NME").head_length = 1.3
ResidueType.get_type("NME").head_link_conditions.append({"atoms": ["CH3", "N"], "parameter": 120 / 180 * np.pi})
ResidueType.get_type("NME").head_link_conditions.append({"atoms": ["H", "CH3", "N"], "parameter": -np.pi})

GlobalSetting.HISMap["DeltaH"] = "HD1"
GlobalSetting.HISMap["EpsilonH"] = "HE2"
GlobalSetting.HISMap["HIS"].update({"HIS": {"HID": "HID", "HIE": "HIE", "HIP": "HIP"},
                                    "CHIS": {"HID": "CHID", "HIE": "CHIE", "HIP": "CHIP"},
                                    "NHIS": {"HID": "NHID", "HIE": "NHIE", "HIP": "NHIP"}})

ResidueType.get_type("CYX").connect_atoms["ssbond"] = "SG"

# pylint:disable=undefined-variable
