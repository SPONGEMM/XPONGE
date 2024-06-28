"""
This **module** sets the basic configuration of coarse-grained model of water (mW)
"""
import os
from .. import sw_base
from ....helper import Xprint, AtomType
from ....load import load_mol2

MW_DATA_DIR = os.path.dirname(__file__)

AtomType.New_From_String("""name mass SWType charge
mW 18.02 mW 0
""")

sw_base.SWType.New_From_String(""" name epsilon sigma a l gamma b A B p q
mW-mW 6.189  2.3925  1.80  23.15  1.20  -0.333333333333 7.049556277  0.6022245584  4.0  0.0
mW-mW-mW 6.189  2.3925  1.80  23.15  1.20  -0.333333333333 7.049556277  0.6022245584  4.0  0.0
""")

load_mol2(os.path.join(MW_DATA_DIR, "mW.mol2"), as_template=True)

Xprint("""Reference for mW:
    Valeria Molinero* and Emily B. Moore
    Water Modeled As an Intermediate Element between Carbon and Silicon.
    J. Phys. Chem. B 2009, 113, 13, 4008–4016
""")
