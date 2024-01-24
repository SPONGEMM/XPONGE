"""
    This **module** includes unittests of the basic functions of Xponge.process
"""
import os

__all__ = ["test_all"]

def test_all():
    """
        test the building of all functions
    """
    import Xponge
    import Xponge.forcefield.amber.tip3p

    box = Xponge.BlockRegion(0, 0, 0, 60, 60, 60)
    region_1 = Xponge.BlockRegion(0, 0, 20, 20, 20, 40)
    region_2 = Xponge.BlockRegion(0, 0, 40, 20, 20, 60)
    region_3 = Xponge.BlockRegion(0, 0, 0, 20, 20, 20)
    region_4 = Xponge.SphereRegion(20, 10, 30, 10)
    region_5 = Xponge.BlockRegion(20, 0, 20, 60, 60, 60)
    region_2or3 = Xponge.UnionRegion(region_2, region_3)
    region_4and5 = Xponge.IntersectRegion(region_4, region_5)
    t = Xponge.Lattice("bcc", basis_molecule=Xponge.ResidueType.get_type("CL"), scale=4)
    t2 = Xponge.Lattice("fcc", basis_molecule=Xponge.ResidueType.get_type("K"), scale=3)
    t3 = Xponge.Lattice("sc", basis_molecule=Xponge.ResidueType.get_type("NA"), scale=3)
    mol = t.Create(box, region_1)
    mol = t2.create(box, region_2or3, mol)
    mol = t3.create(box, region_4and5, mol)
    Xponge.Save_PDB(mol, "out.pdb")
