"""
This **package** sets the basic configuration of charmm36-july2022 force field
"""
import os
from ....helper import source, Xprint
from .. import load_parameter_from_ffitp

CHARMM36_DATA_DIR = os.path.dirname(__file__)

load_parameter_from_ffitp("forcefield.itp", CHARMM36_DATA_DIR)

source(".protein")
source(".dna")
source(".rna")

Xprint("""Reference for charmm36-july2022:
    Huang, J., Rauscher, S., Nawrocki, G., Ran, T., Feig, M., de Groot, B.L.,
Grubmuller, H., and MacKerell, Jr., A.D. CHARMM36m: an improved force field
for folded and intrinsically disordered proteins, Nature Methods, 2016, DOI:
10.1038/nmeth.4067.
""")
