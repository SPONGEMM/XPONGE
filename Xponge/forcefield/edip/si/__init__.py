"""
This **module** sets the basic configuration of coarse-grained model of water (mW)
"""
import os
from .. import edip_base
from ....helper import Xprint, AtomType
from ....load import load_mol2

SI_DATA_DIR = os.path.dirname(__file__)

AtomType.New_From_String("""name mass EDIPType charge
Si 28.09 Si 0
""")

edip_base.EDIPType.New_From_String(""" name A[eV]   B   a   c   alpha   beta   eta gamma   l[eV]    mu   rho   sigma   Q0  u1   u2   u3   u4
Si-Si 7.9821730 1.5075463 3.1213820 2.5609104 3.1083847 0.0070975 0.2523244 1.1247945 1.4533108 0.6966326 1.2085196 0.5774108 312.1341346 -0.165799 32.557 0.286198 0.66
Si-Si-Si 7.9821730 1.5075463 3.1213820 2.5609104 3.1083847 0.0070975 0.2523244 1.1247945 1.4533108 0.6966326 1.2085196 0.5774108 312.1341346 -0.165799 32.557 0.286198 0.66
""")

load_mol2(os.path.join(SI_DATA_DIR, "Si.mol2"), as_template=True)

Xprint("""Reference for EDIP Si:
    Justo, Bazant, Kaxiras, Bulatov and Yip
    Interatomic potential for silicon defects and disordered phases.
    Phys. Rev, B 1998, 58, 5, 2539–2550
""")
