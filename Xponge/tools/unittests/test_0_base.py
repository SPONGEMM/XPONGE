"""
    This **module** includes unittests of the basic functions
"""
__all__ = ["test_import",
           "test_assign",
           "test_molecule"]

def test_import():
    """test importing modules"""
    import Xponge
    import Xponge.forcefield.amber.tip3p
    water = Xponge.ResidueType.get_type("WAT")
    Xponge.Xprint(water.name)

def test_assign():
    """test creating an assign by hand"""
    import Xponge
    assign = Xponge.Assign()
    assign.add_atom("O", 0, 0, 0)
    assign.addAtom("H1", 0, 1, 0)
    assign.Add_Atom("H2", 1, 0, 0)

def test_molecule():
    """test creating an assign by hand"""
    import Xponge
    import Xponge.forcefield.base.mass_base
    import Xponge.forcefield.base.charge_base
    Xponge.AtomType.New_From_String(r"""
name  mass   charge[e]
H     1.008  +0.25
C     12.00  -1.00
""")
    h = Xponge.AtomType.get_type("H")
    c = Xponge.AtomType.get_type("C")
    xyj = Xponge.ResidueType(name="XYJ")
    xyj.add_atom("H1", h, 1, 0, 0)
    xyj.add_atom("H2", "H", 0, 1, 0)
    xyj.add_atom("H3", h, -1, 0, 0)
    xyj.add_atom("H4", h, 0, -1, 0)
    xyj.add_atom("C", c, 0, 0, 0)
    mol = Xponge.Molecule(name="XYJ2")
    mol.add_residue(xyj)
