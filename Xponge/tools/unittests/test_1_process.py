"""
    This **module** includes unittests of the basic functions of Xponge.process
"""
import os

__all__ = ["test_impose",
           "test_solvent_process",
           "test_lattice"]

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
    t = Xponge.add_solvent_box(t, WAT, 20)
    Xponge.main_axis_rotate(t)
    Xponge.solvent_replace(t, WAT, {NA: 5, CL: 5})
    Xponge.h_mass_repartition(t)
    Xponge.save_sponge_input(t, "AlA")
    assert run(f"SPONGE -default_in_file_prefix ALA -mode NPT -thermostat middle_langevin \
-dt 1e-3 -constrain_mode SHAKE -barostat andersen_barostat -write_information_interval 100 \
> {os.devnull} 2> {os.devnull}") == 0

def test_lattice():
    """
        test the lattice building system
    """
    import Xponge
    import Xponge.forcefield.amber.tip3p
    box = Xponge.BlockRegion(0, 0, 0, 60, 60, 60)
    region_1 = Xponge.BlockRegion(0, 0, 20, 20, 20, 40)
    region_2 = Xponge.BlockRegion(0, 0, 40, 20, 20, 60)
    region_3 = Xponge.BlockRegion(0, 0, 0, 20, 20, 20)
    region_4 = Xponge.SphereRegion(20, 10, 30, 10)
    region_5 = Xponge.BlockRegion(0, 0, 0, 20, 20, 40, side="out")
    region_2or3 = Xponge.UnionRegion(region_2, region_3)
    region_4and5 = Xponge.IntersectRegion(region_4, region_5)
    region_6 = Xponge.FrustumRegion(10, 40, 0, 15, 10, 40, 60, 1)
    region_7 = Xponge.PrismRegion(30, 30, 0, 20, 0, 0, 0, 20, 0, 10, 10, 20)
    t = Xponge.Lattice("bcc", basis_molecule=CL, scale=4)
    t2 = Xponge.Lattice("fcc", basis_molecule=K, scale=3)
    t3 = Xponge.Lattice("sc", basis_molecule=NA, scale=3)
    t4 = Xponge.Lattice("hcp", basis_molecule=MG2, scale=4)
    t5 = Xponge.Lattice("diamond", basis_molecule=AL3, scale=5)
    mol = t.Create(box, region_1)
    mol = t2.create(box, region_2or3, mol)
    mol = t3.create(box, region_4and5, mol)
    mol = t4.create(box, region_6, mol)
    mol = t5.create(box, region_7, mol)
