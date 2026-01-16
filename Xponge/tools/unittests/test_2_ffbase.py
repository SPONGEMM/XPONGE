"""
    This **module** includes unittests of the Xponge.forcefield.base
"""

__all__ = ["test_lj",
           "test_atomwise",
           "test_bond_gmx_parser"]

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

def test_bond_gmx_parser():
    """
        test the gmx bond parser uses `mol` not an undefined name
    """
    from Xponge.forcefield.base.bond_base import BondType, _gmx_parser

    class DummyAtomType:
        def __init__(self, name):
            self.name = name
        def __repr__(self):
            return f"Type of Atom: {self.name}"

    class DummyAtom:
        def __init__(self, type_):
            self.type = type_

    class DummyMol:
        def __init__(self):
            self.forces = []
        def add_bonded_force(self, force):
            self.forces.append(force)

    BondType.New_From_String("""
name k b
XTBOND_A-XTBOND_B 100.0 1.5
XTBOND_D-XTBOND_C 200.0 1.6
""")

    # Direct match branch
    mol1 = DummyMol()
    stat1 = {1: DummyAtom(DummyAtomType("XTBOND_A")), 2: DummyAtom(DummyAtomType("XTBOND_B"))}
    _gmx_parser(["1", "2", "1"], mol1, stat1)
    assert len(mol1.forces) == 1

    # Reverse match branch
    mol2 = DummyMol()
    stat2 = {1: DummyAtom(DummyAtomType("XTBOND_C")), 2: DummyAtom(DummyAtomType("XTBOND_D"))}
    _gmx_parser(["1", "2", "1"], mol2, stat2)
    assert len(mol2.forces) == 1

    # Missing type should raise a clear error
    mol3 = DummyMol()
    stat3 = {1: DummyAtom(DummyAtomType("XTBOND_MISSING")), 2: DummyAtom(DummyAtomType("XTBOND_ALSO_MISSING"))}
    try:
        _gmx_parser(["1", "2", "1"], mol3, stat3)
    except KeyError:
        pass
    else:
        raise AssertionError("Expected KeyError for missing bond type")
