"""
    This **module** includes unittests of the basic functions of Xponge.process
"""
import os

__all__ = ["test_fcc"]

def test_fcc():
    """
        test the building of a fcc latice
    """
    import Xponge
    import Xponge.forcefield.amber.tip3p
    globals().update(Xponge.ResidueType.get_all_types())
    mol0 = NA + CL
    mol0.residues[-1].CL.x = 4
    region = Xponge.BlockRegion(0, 0, 0, 32, 32, 32)
    box = Xponge.BlockRegion(0, 0, 0, 40, 40, 40)
    lattice = Xponge.Lattice("fcc", mol0, 8)
    mol = lattice.Create(box, region)
    Save_PDB(mol, "NaCl.pdb")
    Save_SPONGE_Input(mol, "NaCl")
