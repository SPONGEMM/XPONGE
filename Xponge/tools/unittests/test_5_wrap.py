"""
    This **module** gives the unit tests of the wrapping molecules
"""

__all__ = ["test_unperiodic", "test_periodic"]

def test_unperiodic():
    """ test the wrapping of unperiodic molecules """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.mdrun import run
    from Xponge.helper.cv import CVSystem

    ala = Xponge.ResidueType.Get_Type("ALA")
    up = Xponge.ResidueType.Get_Type("ACE") + ala * 10 + Xponge.ResidueType.Get_Type("NME")
    Xponge.add_solvent_box(up, Xponge.ResidueType.Get_Type("WAT"), 10)
    Xponge.save_sponge_input(up, "up")
    cv = CVSystem(up)
    cv.add_center("c", "protein")
    cv.add_cv_position("x", "c", "x", scaled=False)
    cv.print("x")
    cv.add_cv_rmsd("r", "protein")
    cv.print("r")
    cv.steer("x", 2)
    cv.output("cv_up.txt")
    step_limit = 50000
    assert run("SPONGE -mode minimization -default_in_file_prefix up -rst min > min.out") == 0
    assert run(f"SPONGE -mode NPT -thermostat middle_langevin -barostat andersen_barostat \
-default_in_file_prefix up -step_limit {step_limit} -dt 2e-3 -constrain_mode SHAKE -cutoff 8 \
-cv_in_file cv_up.txt -coordinate_in_file min_coordinate.txt -crd up.dat -box up.box > run.out") == 0

def test_periodic():
    """ test the wrapping of periodic molecules """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.mdrun import run
    from Xponge.helper.cv import CVSystem

    ala = Xponge.ResidueType.Get_Type("ALA")
    p = ala * 10
    p.add_residue_link(p.residues[0].N, p.residues[-1].C)
    Xponge.add_solvent_box(p, Xponge.ResidueType.Get_Type("WAT"), 10)
    Xponge.save_sponge_input(p, "p")
    cv = CVSystem(p)
    cv.add_center("c", "protein")
    cv.add_cv_position("x", "c", "x", scaled=False)
    cv.print("x")
    cv.add_cv_rmsd("r", "protein")
    cv.print("r")
    cv.steer("x", 20)
    cv.output("cv_up.txt")
    step_limit = 50000
    assert run("SPONGE -mode minimization -default_in_file_prefix p -step_limit 5000 -rst min > min.out") == 0
    assert run(f"SPONGE -mode NPT -thermostat middle_langevin -barostat andersen_barostat \
-default_in_file_prefix p -step_limit {step_limit} -dt 2e-3 -constrain_mode SHAKE -cutoff 8 \
-cv_in_file cv_up.txt -coordinate_in_file min_coordinate.txt -crd p.dat -box p.box > run.out") == 0
