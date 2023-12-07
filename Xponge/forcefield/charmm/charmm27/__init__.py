"""
This **package** sets the basic configuration of charmm27 force field
"""
import os
from ....helper import source
from .. import load_parameter_from_ffitp

CHARMM27_DATA_DIR = os.path.dirname(__file__)

load_parameter_from_ffitp("forcefield.itp", CHARMM27_DATA_DIR)

source(".protein")
