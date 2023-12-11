"""
This **package** sets the basic configuration of charmm27 force field
"""
import os
from ....helper import source
from .. import load_parameter_from_ffitp

CHARMM27_DATA_DIR = os.path.dirname(__file__)

load_parameter_from_ffitp("forcefield.itp", CHARMM27_DATA_DIR)

source(".protein")
source(".dna")
source(".rna")

Xprint("""Reference for charmm27:
    MacKerell, Jr., A. D., Feig, M., Brooks, C.L., III, Extending the
    treatment of backbone energetics in protein force fields: limitations
    of gas-phase quantum mechanics in reproducing protein conformational
    distributions in molecular dynamics simulations, Journal of
    Computational Chemistry, 25: 1400-1415, 2004.

and 

    MacKerell, Jr., A. D.,  et al. All-atom
    empirical potential for molecular modeling and dynamics Studies of
    proteins.  Journal of Physical Chemistry B, 1998, 102, 3586-3616.

""")
