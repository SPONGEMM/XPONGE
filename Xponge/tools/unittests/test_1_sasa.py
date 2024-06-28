"""
    This **module** includes unittests of the Solvent Accessible Surface Area
"""
import os

__all__ = ["test_surface"]

def test_surface():
    """
        test the calculation of the surface
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    from Xponge.analysis.sasa import SASA
    from Xponge.analysis.md_analysis import XpongeMoleculeReader
    import MDAnalysis as mda

    t = Xponge.ResidueType.get_type("ALA") * 3
    Xponge.save_pdb(t, "protein.pdb")
    u = mda.Universe(t, format=XpongeMoleculeReader)
    sasa = SASA(u, n_points=1000)
    sasa.main()
    sasa.write_surface_xyz("surface.xyz")
    sasa.write_surface_xyz("surface.txt", headlines=False)

