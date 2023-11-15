"""
    This **module** includes unittests of the Xponge.forcefield.amber.ff14sb
"""

__all__ = ["test_ff14sb"]

def test_ff14sb():
    """
        test the single point energy for residues
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run
    import matplotlib.pyplot as plt
    import numpy as np
    import os
    import re

    s = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU " + \
        "LYS MET PHE PRO SER THR TRP TYR VAL HIS"

    with open("leaprc", "w") as f:
        f.write(f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p
t = sequence {{ACE {s} NME}}
solvatebox t WAT 10
saveamberparm t t.parm7 t.rst7
quit""")
    with open("mdin", "w") as f:
        f.write(f"""test ff14SB
&cntrl
  nstlim = 1000
  ntt = 3
  temp0 = 300
  ntwx = 100
  ntpr = 100
/
""")
    assert os.system("tleap > tleap.out 2> tleap.out") == 0
    assert run(f"SPONGE -mode minimization -amber_parm7 t.parm7 -amber_rst7 t.rst7 -rst min.rst7 > {os.devnull}") == 0
    assert os.system("pmemd.cuda -i mdin -p t.parm7 -c min.rst7 -x amber.nc -O > pmemd.out 2> pmemd.out") == 0
    assert os.system("Xponge converter -p t.parm7 -c amber.nc -o amber.dat -of sponge_traj") == 0

    with open("t.parm7") as f:
        n_solvent = f.read().count("WAT ")
    mol = Xponge.ResidueType.get_type("ACE")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("NME")
    Xponge.add_solvent_box(mol, Xponge.ResidueType.get_type("WAT"), 20, n_solvent=n_solvent)
    Xponge.save_sponge_input(mol, "ff14sb")
    assert run(f"SPONGE -mode rerun -default_in_file_prefix ff14sb " + \
               "-cutoff 8 -crd amber.dat -box amber.box > {os.devnull} ") == 0

    with open("mdout") as f:
        matches = re.findall(r"EPtot\s*=\s*(\d*.\d*)", f.read())
        matches = np.array([float(match) for match in matches[1:-2]])
    t = MdoutReader("mdout.txt")
    plt.plot(matches / 1000, t.potential / 1000, "o--")
    plt.xlabel("Energy from AMBER [Mcal/mol]")
    plt.ylabel("Energy from SPONGE [Mcal/mol]")
    plt.savefig("ff14sb.png")

