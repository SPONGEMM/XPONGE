"""
    This **module** includes unittests of the Xponge.forcefield.base
"""

__all__ = ["test_lj",
           "test_atomwise"]

def test_lj():
    """
        test the unit convertion
    """
    from Xponge.forcefield.base.lj_base import LJType
    import numpy as np

    LJType.New_From_String(r"""
name      A         B
AG-AG     2021000   6072
AL-AL     1577000   5035
AU-AU     2307000	6987
""")

    LJType.New_From_String(f"""
name    epsilon   rmin
ag-ag     4.56    {2.955 / 2}
al-al     4.02    {2.925 / 2}
au-au     5.29    {2.951 / 2}
""")

    for name in ["Ag-Ag", "Al-Al", "Au-Au"]:
        er = LJType.get_type(name.lower())
        ab = LJType.get_type(name.upper())
        assert abs(er.epsilon - ab.epsilon) < 0.01, f"{name} epsilon does not match"
        assert abs(er.rmin - ab.rmin) < 0.01, f"{name} rmin does not match"

    LJType.New_From_String(r"""
name    epsilon[eV]   sigma[nm]
y-y     0.0017345     0.32
""")
    yy = LJType.get_type("y-y")
    assert abs(yy.epsilon - 0.03999851) < 0.01
    assert abs(yy.rmin - np.power(2, 1/6) * 3.2 / 2) < 0.01

def test_atomwise():
    """
        test the atomwise forcefield base
    """
    import Xponge
    import Xponge.forcefield.base.mass_base
    import Xponge.forcefield.base.charge_base
    import Xponge.forcefield.base.lj_base

    Xponge.AtomType.New_From_String(r"""
name    mass    charge[e]   LJtype
H       1.008   1.000       HW
""")
