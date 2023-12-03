"""
    This **module** gives the unit tests of the simulation for nopbc
"""

__all__ = ["test_rerun", "test_cv_run"]

def _check_one_energy(amber_mdout, amber_name, sponge_out):
    """check one energy term"""
    import re
    import matplotlib.pyplot as plt
    import numpy as np
    sponge_out = sponge_out[:-1]
    with open(amber_mdout) as f:
        t = f.read()
        matches = re.findall(rf"{amber_name}\s*=\s*(\d*.\d*)", t)
        matches = np.array([float(match) for match in matches[2:-2]])
    unit = "kcal/mol"
    if abs(np.mean(matches)) > 1000:
        unit = "Mcal/mol"
        matches /= 1000
        sponge_out /= 1000
    k, b = np.polyfit(matches, sponge_out, 1)
    r = np.corrcoef(matches, sponge_out)[0][1]
    plt.plot([np.min(matches), np.max(matches)],
             [k * np.min(matches) + b, k * np.max(matches) + b],
             label=f"y={k:.3f}x{b:+.3f},r={r:.3f}")
    plt.plot(matches, sponge_out, "o")
    plt.xlabel(f"Result from AMBER [{unit}]")
    plt.ylabel(f"Result from SPONGE [{unit}]")
    plt.legend()
    plt.savefig(f"{amber_name}.png")
    plt.clf()

def test_rerun():
    """
        test the single point energy for a general Born system
    """
    import os
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.forcefield.special import gb
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run


    s = "ALA ALA ALA ALA ALA " + \
        "ALA ALA ALA ALA ALA "

    with open("leaprc", "w") as f:
        f.write(f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p
set default PBRadii mbondi
t = sequence {{ACE {s} NME}}
saveamberparm t t.parm7 t.rst7
quit""")
    with open("mdin", "w") as f:
        f.write("""test ff14SB
&cntrl
  nstlim = 1000
  igb = 1
  cut = 999
  ntt = 3
  temp0 = 300
  ntwx = 1
  ntpr = 1
/
""")
    assert os.system("tleap > tleap.out 2> tleap.out") == 0
    assert os.system("pmemd.cuda -i mdin -p t.parm7 -c t.rst7 -x amber.nc -O > pmemd.out 2> pmemd.out") == 0
    assert os.system("Xponge converter -p t.parm7 -c amber.nc -o amber.dat -of sponge_traj") == 0

    with open("amber.box", "w") as f:
        f.write("999 999 999 90 90 90")
    mol = Xponge.ResidueType.get_type("ACE")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("NME")
    gb.set_gb_radius(mol)
    Xponge.save_sponge_input(mol, "gb")
    assert run("SPONGE -mode rerun -default_in_file_prefix gb " + \
               "-cutoff 999 -pbc 0 -crd amber.dat -box amber.box > rerun.out ") == 0

    t = MdoutReader("mdout.txt")
    _check_one_energy("mdout", " EPtot", t.potential)
    _check_one_energy("mdout", " BOND", t.bond)
    _check_one_energy("mdout", " ANGLE", t.angle)
    _check_one_energy("mdout", " DIHED", t.dihedral)
    _check_one_energy("mdout", " VDWAALS", t.LJ)
    _check_one_energy("mdout", " EELEC", t.Coulomb)
    _check_one_energy("mdout", " 1-4 NB", t.nb14_LJ)
    _check_one_energy("mdout", " 1-4 EEL", t.nb14_EE)
    _check_one_energy("mdout", " EGB", t.gb)

def test_cv_run():
    """
        test the steer MD simulation without pbc
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.forcefield.special import gb
    from Xponge.helper.cv import CVSystem
    from Xponge.mdrun import run


    s = "ALA ALA ALA ALA ALA " + \
        "ALA ALA ALA ALA ALA "

    mol = Xponge.ResidueType.get_type("ACE")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("NME")
    gb.set_gb_radius(mol)
    Xponge.save_sponge_input(mol, "cv")
    cv = CVSystem(mol)
    cv.add_center("c", "protein")
    cv.add_cv_position("x", "c", "x", False)
    cv.steer("x", 100)
    cv.output("cv.txt")
    assert run("SPONGE -mode NVT -thermostat andersen_thermostat -default_in_file_prefix cv " + \
               "-cutoff 999 -pbc 0 -cv_in_file cv.txt -step_limit 100000 > cv.out ") == 0
