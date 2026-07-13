"""Reusable PACKMOL-Memgen lipid extension registration."""

import os

from ...helper import Xprint
from ...load import load_mol2
from . import AMBER_DATA_DIR, load_parameters_from_frcmod, load_parameters_from_parmdat
from ._forcefield_family import require_forcefield_family
from ._lipid_common import configure_manifest


LIPID_EXTENSION_REFERENCE_TEXT = """Reference for Xponge lipid extensions:
1. PI, phosphoinositide, and LysoPL extensions
  Schott-Verdugo, S.; Gohlke, H.
    PACKMOL-Memgen: A Simple-To-Use, Generalized Workflow for Membrane-Protein-Lipid-Bilayer System Building.
    Journal of Chemical Information and Modeling 2019 59(6), 2522-2528.
    DOI: 10.1021/acs.jcim.9b00269

2. Supporting parameter families
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
"""


_REGISTERED = False


def register_lipid_extension():
    """Register the shared extension once after a supported lipid base."""
    global _REGISTERED
    require_forcefield_family("lipid", {"lipid17", "lipid21"})
    if _REGISTERED:
        return None
    load_parameters_from_parmdat(os.path.join("glycam_06j", "GLYCAM_06j.dat"))
    load_parameters_from_frcmod("frcmod.lipid_ext")
    load_mol2(os.path.join(AMBER_DATA_DIR, "lipid_ext.mol2"), as_template=True)
    manifest = configure_manifest(os.path.join(AMBER_DATA_DIR, "lipid_ext_manifest.json"))
    Xprint(LIPID_EXTENSION_REFERENCE_TEXT)
    _REGISTERED = True
    return manifest

