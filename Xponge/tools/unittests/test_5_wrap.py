"""
    This **module** gives the unit tests of the speed of SPONGE NVT simulation
"""

__all__ = ["test_unperiodic"]

def test_unperiodic():
    """ test the speed of NVT simulation """
    import os
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.mdrun import run
    from io import StringIO
    from Xponge.analysis import MdoutReader
    import numpy as np
    from Xponge.helper.cv import CVSystem

    ala = Xponge.ResidueType.Get_Type("ALA")
    up = Xponge.ResidueType.Get_Type("ACE") + ala * 10 + Xponge.ResidueType.Get_Type("NME")
    Xponge.add_solvent_box(up, Xponge.ResidueType.Get_Type("WAT"), 20)
    Xponge.save_sponge_input(up, "up")
    cv = CVSystem(up)
    cv.add_center("c", "protein")
    cv.add_cv_position("x", "c", "x")
    cv.print("x")
    cv.steer("x", 1)
    cv.output("cv_up.txt")
    step_limit = 50000
    assert run(f"SPONGE -mode minimization -default_in_file_prefix up -rst min >> test_unperiodic.log") == 0
    assert run(f"SPONGE -mode NPT -thermostat middle_langevin -barostat andersen_barostat \
-default_in_file_prefix up -step_limit {step_limit} -dt 2e-3 -constrain_mode SHAKE -cutoff 8 \
-cv_in_file cv_up.txt -coordinate_in_file min_coordinate.txt >> test_unperiodic.log") == 0

