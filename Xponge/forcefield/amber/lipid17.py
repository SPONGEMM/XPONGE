"""
This **module** set the basic configuration for lipid17
"""
from ._forcefield_family import activate_forcefield_family

activate_forcefield_family("lipid", "lipid17")

import json

from ...helper import source, Xprint

source("....")
amber = source("...amber")
amber.load_parameters_from_parmdat("lipid17.dat")
load_mol2(os.path.join(AMBER_DATA_DIR, "lipid17.mol2"), as_template=True)
amber.load_parameters_from_parmdat(os.path.join("glycam_06j", "GLYCAM_06j.dat"))
amber.load_parameters_from_frcmod("frcmod.lipid_ext")
load_mol2(os.path.join(AMBER_DATA_DIR, "lipid_ext.mol2"), as_template=True)

for res in "LAL PA MY OL ST AR DHA".split():
    ResidueType.get_type(res).head = "C12"
    ResidueType.get_type(res).tail = "C12"
    ResidueType.get_type(res).head_next = "C13"
    ResidueType.get_type(res).tail_next = "C13"
    ResidueType.get_type(res).head_length = 1.5
    ResidueType.get_type(res).tail_length = 1.5
    ResidueType.get_type(res).head_link_conditions.append({"atoms": ["H2R", "C12"], "parameter": 109.5 / 180 * np.pi})
    ResidueType.get_type(res).head_link_conditions.append(
        {"atoms": ["H2S", "H2R", "C12"], "parameter": -120 / 180 * np.pi})
    ResidueType.get_type(res).tail_link_conditions.append({"atoms": ["H2R", "C12"], "parameter": 109.5 / 180 * np.pi})
    ResidueType.get_type(res).tail_link_conditions.append(
        {"atoms": ["H2S", "H2R", "C12"], "parameter": -120 / 180 * np.pi})

for res in "PC PE PS PGR PH-".split():
    ResidueType.get_type(res).head = "C11"
    ResidueType.get_type(res).tail = "C21"
    ResidueType.get_type(res).head_next = "O11"
    ResidueType.get_type(res).tail_next = "O21"
    ResidueType.get_type(res).head_length = 1.5
    ResidueType.get_type(res).tail_length = 1.5
    ResidueType.get_type(res).head_link_conditions.append({"atoms": ["O11", "C11"], "parameter": 120 / 180 * np.pi})
    ResidueType.get_type(res).head_link_conditions.append({"atoms": ["O12", "O11", "C11"], "parameter": np.pi})
    ResidueType.get_type(res).tail_link_conditions.append({"atoms": ["O21", "C21"], "parameter": 120 / 180 * np.pi})
    ResidueType.get_type(res).tail_link_conditions.append({"atoms": ["O22", "O21", "C21"], "parameter": np.pi})

with open(os.path.join(AMBER_DATA_DIR, "lipid_ext_manifest.json"), encoding="utf-8") as manifest_file:
    lipid_ext_manifest = json.load(manifest_file)

for entry in lipid_ext_manifest["templates"]:
    res = ResidueType.get_type(entry["template"])
    if entry["head_atom"] is not None:
        res.head = entry["head_atom"]
        res.head_next = entry["head_next_atom"]
        res.head_length = 1.5
        res.head_link_conditions.append(
            {"atoms": [entry["head_next_atom"], entry["head_atom"]], "parameter": 120 / 180 * np.pi}
        )
        res.head_link_conditions.append(
            {
                "atoms": [entry["head_reference_atom"], entry["head_next_atom"], entry["head_atom"]],
                "parameter": np.pi,
            }
        )
    if entry["tail_atom"] is not None:
        res.tail = entry["tail_atom"]
        res.tail_next = entry["tail_next_atom"]
        res.tail_length = 1.5
        res.tail_link_conditions.append(
            {"atoms": [entry["tail_next_atom"], entry["tail_atom"]], "parameter": 120 / 180 * np.pi}
        )
        res.tail_link_conditions.append(
            {
                "atoms": [entry["tail_reference_atom"], entry["tail_next_atom"], entry["tail_atom"]],
                "parameter": np.pi,
            }
        )

Xprint("""Reference for Lipid17 and Xponge Lipid17 extensions:
1. Lipid14 / Lipid17 base parameters
  Dickson, C.J., Madej, B.D., Skjevik, A.A., Betz, R.M., Teigen, K., Gould, I.R., Walker, R.C.
    Lipid14: The Amber Lipid Force Field.
    Journal of Chemical Theory and Computation 2014 10(2), 865-879,
    DOI: 10.1021/ct4010307

2. PI, phosphoinositide, and LysoPL extensions
  Schott-Verdugo, S.; Gohlke, H.
    PACKMOL-Memgen: A Simple-To-Use, Generalized Workflow for Membrane-Protein-Lipid-Bilayer System Building.
    Journal of Chemical Information and Modeling 2019 59(6), 2522-2528.
    DOI: 10.1021/acs.jcim.9b00269

3. Supporting parameter families
  Kirschner, K.N. et al.
    GLYCAM06: A generalizable biomolecular force field. Carbohydrates.
    Journal of Computational Chemistry 2008 29(4), 622-655.
    DOI: 10.1002/jcc.20820

  Homeyer, N.; Horn, A.H.C.; Lanig, H.; Sticht, H.
    Revised AMBER Parameters for Bioorganic Phosphates.
    Journal of Chemical Theory and Computation 2012 8(11), 4405-4412.
    DOI: 10.1021/ct300613v

  Selected extension terms are derived from GAFF2 as documented in the source frcmod and Amber Reference Manual.
  GAFF lineage: Wang, J. et al., Journal of Computational Chemistry 2004 25, 1157-1174.
    DOI: 10.1002/jcc.20035
""")
# pylint:disable=undefined-variable
