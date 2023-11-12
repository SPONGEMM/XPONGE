"""
    This **module** includes unittests of the Xponge.Assign operators
"""
from io import StringIO

__all__ = ["test_get_assign",
           "test_residuetype_to_assign",
           "test_objective_assign_rule",
           "test_string_assign_rule",
           "test_ring_system",
           "test_atom_deletion",
           "test_bond_deletion",
           "test_equal_atom_search",
           "test_atom_type_determination",
           "test_uff_optimization",
           "test_saving"]

def _check(assign):
    """
        Check whether a Xponge.assignment represents a molecule of benzene
    """
    if assign.atom_numbers != 12:
        return "#atom != 12"
    count1 = 0
    count2 = 0
    for atom in assign.atoms:
        if atom == "C":
            count1 += 1
        elif atom == "H":
            count2 += 1
        else:
            return "Elements except H and C exist"
    if count1 != 6 or count2 != 6:
        return "#H != 6"
    for atom_i, bonds in assign.bonds.items():
        if atom_i == "H":
            if len(bonds) != 1:
                return f"H has {len(bonds)} bonds"
            atom_j, order = next(bonds.items())
            if atom_j != "C" or order != 1:
                return "wrong C-H bond"
        if atom_i == "C":
            if len(bonds) != 3:
                return "#bonds of carbon != 3"
            count1 = 0
            count2 = 0
            for atom_j, order in bonds.items():
                count1 += order
                if atom_j == "C":
                    count2 += 1
            if count1 != 4 or count2 != 2:
                return "Value of C is not right"
    return None


def test_get_assign():
    """
        Test the functions to get assignment
    """
    import Xponge
    error = _check(Xponge.get_assignment_from_smiles("c1ccccc1"))
    if error is not None:
        raise ValueError("smiles", error)
    error = _check(Xponge.get_assignment_from_pubchem("benzene", "name"))
    if error is not None:
        raise ValueError("pubchem", error)
    s = StringIO(r"""
ATOM      0    C ASN     0       0.000   0.000   0.000                      C
ATOM      1    C BEN     1      -1.213  -0.688   0.000                      C
ATOM      2   C1 BEN     1      -1.203   0.706   0.000                      C
ATOM      3   C2 BEN     1      -0.010  -1.395   0.000                      C
ATOM      4   C3 BEN     1       0.010   1.395  -0.000                      C
ATOM      5   C4 BEN     1       1.203  -0.706   0.000                      C
ATOM      6   C5 BEN     1       1.213   0.688   0.000                      C
ATOM      7    H BEN     1      -2.158  -1.224   0.000                      H
ATOM      8   H1 BEN     1      -2.139   1.256   0.000                      H
ATOM      9   H2 BEN     1      -0.018  -2.481  -0.000                      H
ATOM     10   H3 BEN     1       0.018   2.481   0.000                      H
ATOM     11   H4 BEN     1       2.139  -1.256   0.000                      H
ATOM     12   H5 BEN     1       2.158   1.224   0.000                      H
""")
    error = _check(Xponge.get_assignment_from_pdb(s, "BEN"))
    if error is not None:
        raise ValueError("pdb", error)
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    error = _check(Xponge.get_assignment_from_xyz(s))
    if error is not None:
        raise ValueError("xyz", error)
    s = StringIO(r"""
@<TRIPOS>MOLECULE
ASN
 12 12 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    C  -1.2131  -0.6884   0.0000   C.ar         1      ASN   0.000000
     2   C1  -1.2028   0.7064   0.0001   C.ar         1      ASN   0.000000
     3   C2  -0.0103  -1.3948   0.0000   C.ar         1      ASN   0.000000
     4   C3   0.0104   1.3948  -0.0001   C.ar         1      ASN   0.000000
     5   C4   1.2028  -0.7063   0.0000   C.ar         1      ASN   0.000000
     6   C5   1.2131   0.6884   0.0000   C.ar         1      ASN   0.000000
     7    H  -2.1577  -1.2244   0.0000   H            1      ASN   0.000000
     8   H1  -2.1393   1.2564   0.0001   H            1      ASN   0.000000
     9   H2  -0.0184  -2.4809  -0.0001   H            1      ASN   0.000000
    10   H3   0.0184   2.4808   0.0000   H            1      ASN   0.000000
    11   H4   2.1394  -1.2563   0.0001   H            1      ASN   0.000000
    12   H5   2.1577   1.2245   0.0000   H            1      ASN   0.000000
@<TRIPOS>BOND
     1      1      2 ar
     2      1      3 ar
     3      1      7 1
     4      2      4 ar
     5      2      8 1
     6      3      5 ar
     7      3      9 1
     8      4      6 ar
     9      4     10 1
    10      5      6 ar
    11      5     11 1
    12      6     12 1
@<TRIPOS>SUBSTRUCTURE
    1      ASN      1 ****               0 ****  **** 
""")
    error = _check(Xponge.get_assignment_from_mol2(s))
    if error is not None:
        raise ValueError("mol2", error)

def test_objective_assign_rule():
    """
        test creating a new rule to assign the Xponge.AtomType of one atom
        and the usage of the rule
    """
    import Xponge
    from Xponge.assign import AssignRule
    Xponge.AtomType.New_From_String(r"""
name
H
C
O
""")
    rule = AssignRule("myrule")
    def _new_rule(element):
        return lambda i, assign: assign.atoms[i] == element
    for element in ["H", "O", "C"]:
        rule.add_rule(element)(_new_rule(element))

    def _pre_action(assign):
        assign.atoms[1] = "O"

    def _post_action(assign):
        assign.atom_types[2] = Xponge.AtomType.get_type("C")

    rule.set_pre_action(_pre_action)
    rule.set_post_action(_post_action)
    assign = Xponge.Assign()
    assign.add_atom("H", 0, 0, 0)
    assign.add_atom("C", 1, 0, 0)
    assign.add_atom("O", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.add_bond(2, 1, 2)
    assign.determine_atom_type("myrule")
    assert assign.atom_types[0] == Xponge.AtomType.get_type("H")
    assert assign.atom_types[1] == Xponge.AtomType.get_type("O")
    assert assign.atom_types[2] == Xponge.AtomType.get_type("C")

def test_residuetype_to_assign():
    """
        test convert an Xponge.ResidueType to Xponge.Assign
    """
    import Xponge
    import Xponge.forcefield.amber.gaff

    s = StringIO("""
@<TRIPOS>MOLECULE
ASN
 12 12 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    C  -1.2131  -0.6884   0.0000   ca            1      ASN   0.000000
     2   C1  -1.2028   0.7064   0.0001   ca            1      ASN   0.000000
     3   C2  -0.0103  -1.3948   0.0000   ca            1      ASN   0.000000
     4   C3   0.0104   1.3948  -0.0001   ca            1      ASN   0.000000
     5   C4   1.2028  -0.7063   0.0000   ca            1      ASN   0.000000
     6   C5   1.2131   0.6884   0.0000   ca            1      ASN   0.000000
     7    H  -2.1577  -1.2244   0.0000   ha            1      ASN   0.000000
     8   H1  -2.1393   1.2564   0.0001   ha            1      ASN   0.000000
     9   H2  -0.0184  -2.4809  -0.0001   ha            1      ASN   0.000000
    10   H3   0.0184   2.4808   0.0000   ha            1      ASN   0.000000
    11   H4   2.1394  -1.2563   0.0001   ha            1      ASN   0.000000
    12   H5   2.1577   1.2245   0.0000   ha            1      ASN   0.000000
@<TRIPOS>BOND
     1      1      2 ar
     2      1      3 ar
     3      1      7 1
     4      2      4 ar
     5      2      8 1
     6      3      5 ar
     7      3      9 1
     8      4      6 ar
     9      4     10 1
    10      5      6 ar
    11      5     11 1
    12      6     12 1
@<TRIPOS>SUBSTRUCTURE
    1      ASN      1 ****               0 ****  **** 
""")
    ben0 = Xponge.load_mol2(s)
    ben = Xponge.get_assignment_from_residuetype(ben0.residues[0].type)
    assert _check(ben) is None

def test_string_assign_rule():
    """
        test creating a new rule to assign the string of one atom
        and the usage of the rule
    """
    import Xponge
    from Xponge.assign import AssignRule
    rule = AssignRule("myrule", pure_string=True)
    @rule.add_rule("A", -1)
    def _(i, a): #pylint: disable=unused-argument
        return True

    @rule.add_rule("B", 1)
    def _(i, a): #pylint: disable=unused-argument
        return True

    assign = Xponge.Assign()
    assign.add_atom("H", 0, 0, 0)
    assign.add_atom("C", 1, 0, 0)
    assign.add_atom("O", 1, 0, 0)
    assign.add_bond(0, 1, 1)
    assign.add_bond(2, 1, 2)
    results = assign.determine_atom_type("myrule")
    assert results[0] == "B"
    assert results[1] == "B"
    assert results[2] == "B"

def test_ring_system():
    """
        test the basic functions the the helper class _RING
    """
    import Xponge
    nal = Xponge.assign.get_assignment_from_smiles("c1ccc2ccccc2c1")
    assert len(nal.rings) == 2
    for ring in nal.rings:
        for atom in ring.atoms:
            assert "RG6" in nal.atom_marker[atom]
            assert "AR1" in nal.atom_marker[atom]

def test_atom_deletion():
    """
        test the function to delete an atom from an Xponge.Assign
    """
    import Xponge
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    ben.add_atom("O", 0, 0, 0)
    ben.add_bond(0, 12, 1)
    ben.delete_atom(12)
    assert _check(ben) is None
    ben.delete_atom(0)
    ben.add_atom("C", -1.213, -0.688, 0.000)
    ben.add_bond(12, 0, 4 - sum(ben.bonds[0].values()))
    ben.add_bond(12, 1, 4 - sum(ben.bonds[1].values()))
    ben.add_bond(12, 7, 1)
    assert _check(ben) is None

def test_bond_deletion():
    """
        test the function to delete a bond from an Xponge.Assign
    """
    import Xponge
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    ben.add_bond(0, 11, 1)
    ben.delete_bond(0, 11)
    assert _check(ben) is None

def test_equal_atom_search():
    """
        test the function to find the equal atoms in a molecule
    """
    import Xponge
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    results = ben.determine_equal_atoms()
    assert len(results) == 2
    for result in results:
        result1 = set(result) - set(range(6))
        result2 = set(result) - set(range(6, 12))
        assert not result1 or not result2

def test_atom_type_determination():
    """
        test the function to find the atom type
    """
    import Xponge
    import Xponge.forcefield.amber.gaff
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    ben.determine_atom_type("gaff")
    for i in range(6):
        assert ben.atom_types[i].name == "ca"
    for i in range(6, 12):
        assert ben.atom_types[i].name == "ha"

def test_uff_optimization():
    """
        test the optimization of the molecule using UFF
    """
    import Xponge
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    ben.uff_optimize()

def test_saving():
    """
        test the function to save a Xponge.assign
    """
    import Xponge
    s = StringIO(r"""12
BEN
C -1.213  -0.688   0.000
C -1.203   0.706   0.000
C -0.010  -1.395   0.000
C  0.010   1.395  -0.000 
C  1.213   0.688   0.000
C  1.203  -0.706   0.000
H  0.018   2.481   0.000
H -2.158  -1.224   0.000
H -2.139   1.256   0.000
H -0.018  -2.481  -0.000
H  2.139  -1.256   0.000
H  2.158   1.224   0.000
""")
    ben = Xponge.get_assignment_from_xyz(s)
    ben.save_as_pdb("ben.pdb")
    assert _check(Xponge.get_assignment_from_pdb("ben.pdb")) is None
    ben.save_as_mol2("ben.mol2")
    assert _check(Xponge.get_assignment_from_mol2("ben.mol2")) is None
