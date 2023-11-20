"""
    This **module** includes unittests of the Xponge.forcefield.amber.ol3
"""
import os
import re

__all__ = ["test_ol3"]

def _check_one_energy(amber_mdout, amber_name, sponge_out):
    """check one energy term"""
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

def test_ol3():
    """
        test the single point energy for residues
    """
    import Xponge
    import Xponge.forcefield.amber.ol3
    import Xponge.forcefield.amber.tip3p
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run


    s = "A U C G C U A"

    with open("leaprc", "w") as f:
        f.write(f"""source leaprc.RNA.OL3
source leaprc.water.tip3p
t = sequence {{A5 {s} U3}}
solvatebox t WAT 10
saveamberparm t t.parm7 t.rst7
quit""")
    with open("mdin", "w") as f:
        f.write("""test bsc1
&cntrl
  nstlim = 1000
  ntt = 3
  temp0 = 300
  ntwx = 1
  ntpr = 1
/
""")
    assert os.system("tleap > tleap.out 2> tleap.out") == 0
    assert run("SPONGE -mode minimization -amber_parm7 t.parm7 -amber_rst7 t.rst7 -rst min.rst7 \
-step_limit 2000 -cutoff 8 -minimization_dt_factor 5e-3 > min.out") == 0
    assert os.system("pmemd.cuda -i mdin -p t.parm7 -c min.rst7 -x amber.nc -O > pmemd.out 2> pmemd.out") == 0
    assert os.system("Xponge converter -p t.parm7 -c amber.nc -o amber.dat -of sponge_traj") == 0

    with open("t.parm7") as f:
        n_solvent = f.read().count("WAT ")
    mol = Xponge.ResidueType.get_type("A5")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("U3")
    Xponge.add_solvent_box(mol, Xponge.ResidueType.get_type("WAT"), 20, n_solvent=n_solvent)
    Xponge.save_sponge_input(mol, "ol3")
    assert run("SPONGE -mode rerun -default_in_file_prefix ol3 " + \
               "-cutoff 8 -crd amber.dat -box amber.box > rerun.out ") == 0

    t = MdoutReader("mdout.txt")
    _check_one_energy("mdout", " EPtot", t.potential)
    _check_one_energy("mdout", " BOND", t.bond)
    _check_one_energy("mdout", " ANGLE", t.angle)
    _check_one_energy("mdout", " DIHED", t.dihedral)
    _check_one_energy("mdout", " VDWAALS", t.LJ)
    _check_one_energy("mdout", " EELEC", t.PME)
    _check_one_energy("mdout", " 1-4 NB", t.nb14_LJ)
    _check_one_energy("mdout", " 1-4 EEL", t.nb14_EE)
