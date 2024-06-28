"""
    This **module** includes unittests of the Xponge.forcefield.sw.mw
"""
import os
import re

__all__ = ["test_mw"]

def _check_force():
    import numpy as np
    lines = []
    with open("lammps/dump.atom") as f:
        start = False
        for line in f:
            if start:
                line = line.split()
                lines.append([int(line[0]), float(line[4]), float(line[5]), float(line[6])])
            elif line.startswith("ITEM: ATOMS id x y z fx fy fz"):
                start = True
    lines.sort(key=lambda x: x[0])
    lmp = np.array([line[1:] for line in lines])
    sponge = np.fromfile("sponge/mW_force.dat", dtype=np.float32).reshape(-1,3)
    assert np.mean(sponge - lmp) < 0.01

def test_mw():
    """
        test the single point energy
    """
    import Xponge
    import Xponge.forcefield.sw.mw
    from Xponge.mdrun import run
    import os
    import shutil

    assert os.system("zenodo_get  10.5281/zenodo.12577736 > download.log 2>&1") == 0
    os.makedirs("lammps", exist_ok=True)
    shutil.move("lammps.in", "lammps/lammps.in")
    shutil.move("mW_real.sw", "lammps/mW_real.sw")
    shutil.move("water.lmp", "lammps/water.lmp")
    assert os.system("cd lammps && lmp -i lammps.in > lammps.log 2>&1") == 0
    os.makedirs("sponge", exist_ok=True)
    shutil.move("build.py", "sponge/build.py")
    shutil.move("mdin.txt", "sponge/mdin.txt")
    shutil.move("lammps_coordinate.txt", "sponge/lammps_coordinate.txt")
    os.chdir("sponge")
    os.system("python build.py > xponge.log 2>&1")
    run("SPONGE > sponge.log 2>&1")
    os.chdir("..")
    _check_force()

