import os
from ....helper import source, Xprint
from .. import load_parameter_from_ffitp

MARTINI_DATA_DIR = os.path.join(os.path.dirname(__file__), "martini300")

for _itp in (
    "martini_v3.0.0.itp",
    "martini_v3.0.0_ions_v1.itp",
    "martini_v3.0.0_solvents_v1.itp",
    "martini_v3.0.0_small_molecules_v1.itp",
    "martini_v3.0.0_sugars_v1.itp",
    "martini_v3.0.0_phospholipids_v1.itp",
    "martini_v3.0.0_nucleobases_v1.itp",
):
    load_parameter_from_ffitp(_itp, MARTINI_DATA_DIR)
    
XPrint("""Reference for Martini 3.0.0:
        Souza, P. C. T.; Alessandri, R.; Barnoud, J.; Thallmair, S.; Faustino, I.; Grünewald, F.; Patmanidis, I.; Abdizadeh, H.; Bruininks, B. M. H.; Wassenaar, T. A.; et al.
    Martini 3: a general purpose force field for coarse-grained molecular dynamics.
    Nat. Methods 2021, 18 (4), 382-388 doi:10.1038/s41592-021-01098-3.
""")