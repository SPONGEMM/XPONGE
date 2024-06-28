"""
    This **module** includes unittests of the basic functions of Xponge.process
"""
import os

__all__ = ["test_impose",
           "test_solvent_process"]

def test_impose():
    """
        test imposing bonds, angles and dihedrals
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import numpy as np
    globals().update(Xponge.ResidueType.get_all_types())

    t1 = ACE + ALA + NME
    ala = t1.residues[1]
    pi = 3.141592654
    Xponge.impose_bond(t1, ala.C, ala.CA, 1.2)
    Xponge.impose_angle(t1, ala.C, ala.CA, ala.N, pi / 2)
    Xponge.impose_dihedral(t1, ala.O, ala.C, ala.CA, ala.N, -pi)
    c = np.array([ala.C.x, ala.C.y, ala.C.z])
    ca = np.array([ala.CA.x, ala.CA.y, ala.CA.z])
    o = np.array([ala.O.x, ala.O.y, ala.O.z])
    n = np.array([ala.N.x, ala.N.y, ala.N.z])
    assert abs(np.linalg.norm(c - ca) - 1.2) < 0.01
    v1 = c - ca
    v2 = ca - n
    assert abs(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)) - pi / 2) < 0.01
    v1, v2, v3 = c - o, ca - c, n - ca
    n1, n2 = np.cross(v1, v2), np.cross(v2, v3)
    k = np.sign(np.dot(n1, v3))
    assert abs(k * np.arccos(np.dot(n1, n2) / np.linalg.norm(n1) / np.linalg.norm(n2)) + pi) < 0.01

def test_solvent_process():
    """
        test the function to add solvents box
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    import Xponge.forcefield.amber.tip3p
    from Xponge.mdrun import run
    globals().update(Xponge.ResidueType.get_all_types())
    t = Xponge.get_peptide_from_sequence("AAAAAAAAAA")
    Xponge.main_axis_rotate(t)
    t = Xponge.add_solvent_box(t, WAT, 20)
    Xponge.solvent_replace(t, WAT, {NA: 5, CL: 5})
    Xponge.h_mass_repartition(t)
    Xponge.save_sponge_input(t, "ALA")
    assert run(f"SPONGE -default_in_file_prefix ALA -mode min -cutoff 8 \
-default_out_file_prefix ALA > sol.log 2>&1") == 0
    assert run(f"SPONGE -default_in_file_prefix ALA -mode NPT -thermostat middle_langevin \
-dt 4e-3 -constrain_mode SHAKE -barostat andersen_barostat -cutoff 8 \
-default_out_file_prefix sol > sol.log 2>&1") == 0

