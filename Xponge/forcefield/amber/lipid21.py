"""Register Amber Lipid21 and the shared Xponge lipid extensions."""

import os

from ...helper import Xprint
from ...load import load_mol2
from . import AMBER_DATA_DIR, load_parameters_from_parmdat
from ._forcefield_family import activate_forcefield_family
from ._lipid_common import configure_manifest
from ._lipid_ext import register_lipid_extension


activate_forcefield_family("lipid", "lipid21")

load_parameters_from_parmdat("lipid21.dat")
load_mol2(os.path.join(AMBER_DATA_DIR, "lipid21.mol2"), as_template=True)
lipid21_manifest = configure_manifest(os.path.join(AMBER_DATA_DIR, "lipid21_manifest.json"))

Xprint("""Reference for Lipid21:
  Dickson, C.J.; Walker, R.C.; Gould, I.R.
    Lipid21: Complex Lipid Membrane Simulations with AMBER.
    Journal of Chemical Theory and Computation 2022 18(3), 1726-1736.
    DOI: 10.1021/acs.jctc.1c01217
""")

register_lipid_extension()
