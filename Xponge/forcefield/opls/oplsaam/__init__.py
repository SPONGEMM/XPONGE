"""
This **package** sets the basic configuration of charmm36-july2022 force field
"""
import os
from ....helper import source, Xprint
from .. import load_parameter_from_ffitp

OPLSAAM_DATA_DIR = os.path.dirname(__file__)

load_parameter_from_ffitp("forcefield.itp", OPLSAAM_DATA_DIR)

source(".protein")

Xprint("""Reference for OPLS-AA/m:
        Robertson, M. J.; Tirado-Rives, J.; Jorgensen, W. L.
    Improved Peptide and Protein Torsional Energetics with the OPLS-AA Force Field.
    J. Chem. Theory Comput. 2015, 11 (7), 3499-3509 doi:10.1021/acs.jctc.5b00356.
""")
