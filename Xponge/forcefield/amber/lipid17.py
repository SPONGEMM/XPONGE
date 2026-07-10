"""Register Amber Lipid17 and the shared Xponge lipid extensions."""

import os

from ...helper import Xprint
from ...load import load_mol2
from . import AMBER_DATA_DIR, load_parameters_from_parmdat
from ._forcefield_family import activate_forcefield_family
from ._lipid_common import configure_standard_chain, configure_standard_headgroup
from ._lipid_ext import register_lipid_extension


activate_forcefield_family("lipid", "lipid17")

load_parameters_from_parmdat("lipid17.dat")
load_mol2(os.path.join(AMBER_DATA_DIR, "lipid17.mol2"), as_template=True)

for residue_name in "LAL PA MY OL ST AR DHA".split():
    configure_standard_chain(residue_name)

for residue_name in "PC PE PS PGR PH-".split():
    configure_standard_headgroup(residue_name)

Xprint("""Reference for Lipid17:
  Dickson, C.J., Madej, B.D., Skjevik, A.A., Betz, R.M., Teigen, K., Gould, I.R., Walker, R.C.
    Lipid14: The Amber Lipid Force Field.
    Journal of Chemical Theory and Computation 2014 10(2), 865-879,
    DOI: 10.1021/ct4010307
""")

register_lipid_extension()
