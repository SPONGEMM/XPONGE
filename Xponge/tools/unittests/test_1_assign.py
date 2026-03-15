"""
    This **module** includes unittests of the Xponge.Assign operators
"""
from io import StringIO
from pathlib import Path
import numpy as np

__all__ = ["test_get_assign",
           "test_residuetype_to_assign",
           "test_objective_assign_rule",
           "test_string_assign_rule",
           "test_ring_system",
           "test_atom_deletion",
           "test_bond_deletion",
           "test_equal_atom_search",
           "test_atom_type_determination",
           "test_gaff2_atom_type_determination",
           "test_gaff2_mismatch_samples_alignment",
           "test_gaff2_remaining_n2_alignment",
           "test_gaff2_remaining_sx_alignment",
           "test_gaff_assignment_state_and_residuetype_regression",
           "test_gaff2_special_nitrogen_types",
           "test_gaff2_amide_and_thiocarbonyl_types",
           "test_gaff2_ring_sp3_carbons",
           "test_gaff2_parmchk2",
           "test_uff_optimization",
           "test_saving",
           "test_assign_pdb_hybrid36_save",
           "test_assign_pdb_hybrid36_load"]

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


def test_gaff2_atom_type_determination():
    """
        test the function to find the atom type with gaff2
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2
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
    ben.determine_atom_type("gaff2")
    for i in range(6):
        assert ben.atom_types[i].name == "ca"
    for i in range(6, 12):
        assert ben.atom_types[i].name == "ha"


def test_gaff2_mismatch_samples_alignment():
    """
        regression against known mismatch samples from benchmark set
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2  # noqa: F401

    base = Path(__file__).resolve().parents[3] / "test" / "gaff2_benchmark_1000_v2"
    samples = [
        (7, "CHEMBL404"),      # nd/nc flip
        (15, "CHEMBL407"),     # cd/cc flip
        (451, "CHEMBL32479"),  # cf/ce flip
    ]

    def _parse_mol2_types(path):
        atom_names = []
        atom_types = []
        in_atom = False
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                words = line.split()
                if len(words) >= 6:
                    atom_names.append(words[1])
                    atom_types.append(words[5])
        return atom_names, atom_types

    for idx, chembl_id in samples:
        init_mol2 = base / "mol2_xponge" / f"{idx:04d}_{chembl_id}.mol2"
        amber_mol2 = base / "mol2_amber" / f"{idx:04d}_{chembl_id}.mol2"
        assign = Xponge.get_assignment_from_mol2(str(init_mol2))
        assign.determine_atom_type("gaff2")
        x_types = [t.name for t in assign.atom_types.values()]
        amber_names, amber_types = _parse_mol2_types(amber_mol2)
        assert len(x_types) == len(amber_types)
        for i, (x_type, amber_type) in enumerate(zip(x_types, amber_types)):
            assert x_type == amber_type, (
                f"{chembl_id} atom#{i + 1} {amber_names[i]}: "
                f"amber={amber_type}, xponge={x_type}"
            )


def test_gaff2_remaining_cc_cd_alignment():
    """
        regression for the remaining pure cc/cd orientation mismatches
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2  # noqa: F401

    base = Path(__file__).resolve().parents[3] / "test" / "gaff2_benchmark_1000_v2"
    samples = [
        (194, "CHEMBL534"),
        (638, "CHEMBL64391"),
        (992, "CHEMBL1397"),
    ]

    def _parse_mol2_types(path):
        atom_names = []
        atom_types = []
        in_atom = False
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                words = line.split()
                if len(words) >= 6:
                    atom_names.append(words[1])
                    atom_types.append(words[5])
        return atom_names, atom_types

    for idx, chembl_id in samples:
        init_mol2 = base / "mol2_xponge" / f"{idx:04d}_{chembl_id}.mol2"
        amber_mol2 = base / "mol2_amber" / f"{idx:04d}_{chembl_id}.mol2"
        assign = Xponge.get_assignment_from_mol2(str(init_mol2))
        assign.determine_atom_type("gaff2")
        x_types = [t.name for t in assign.atom_types.values()]
        amber_names, amber_types = _parse_mol2_types(amber_mol2)
        assert len(x_types) == len(amber_types)
        for i, (x_type, amber_type) in enumerate(zip(x_types, amber_types)):
            assert x_type == amber_type, (
                f"{chembl_id} atom#{i + 1} {amber_names[i]}: "
                f"amber={amber_type}, xponge={x_type}"
            )


def test_gaff2_remaining_sx_alignment():
    """
        regression for residual sx mismatch category
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2  # noqa: F401

    base = Path(__file__).resolve().parents[3] / "test" / "gaff2_benchmark_1000_v2"
    samples = [
        (124, "CHEMBL480", 15), # atom#15 S: amber sx, xponge s4
    ]

    def _parse_mol2_types(path):
        atom_types = []
        in_atom = False
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                words = line.split()
                if len(words) >= 6:
                    atom_types.append(words[5])
        return atom_types

    for idx, chembl_id, atom_index in samples:
        init_mol2 = base / "mol2_xponge" / f"{idx:04d}_{chembl_id}.mol2"
        amber_mol2 = base / "mol2_amber" / f"{idx:04d}_{chembl_id}.mol2"
        assign = Xponge.get_assignment_from_mol2(str(init_mol2))
        assign.determine_atom_type("gaff2")
        x_type = assign.atom_types[atom_index - 1].name
        amber_type = _parse_mol2_types(amber_mol2)[atom_index - 1]
        assert x_type == amber_type, f"{chembl_id} atom#{atom_index}: amber={amber_type}, xponge={x_type}"


def test_gaff2_remaining_n2_alignment():
    """
        regression for residual n2 mismatch category
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2  # noqa: F401

    base = Path(__file__).resolve().parents[3] / "test" / "gaff2_benchmark_1000_v2"
    samples = [
        (31, "CHEMBL19", 12),   # amber n2, xponge nd
        (992, "CHEMBL1397", 8), # amber n2, xponge nc
    ]

    def _parse_mol2_types(path):
        atom_types = []
        in_atom = False
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                words = line.split()
                if len(words) >= 6:
                    atom_types.append(words[5])
        return atom_types

    for idx, chembl_id, atom_index in samples:
        init_mol2 = base / "mol2_xponge" / f"{idx:04d}_{chembl_id}.mol2"
        amber_mol2 = base / "mol2_amber" / f"{idx:04d}_{chembl_id}.mol2"
        assign = Xponge.get_assignment_from_mol2(str(init_mol2))
        assign.determine_atom_type("gaff2")
        x_type = assign.atom_types[atom_index - 1].name
        amber_type = _parse_mol2_types(amber_mol2)[atom_index - 1]
        assert x_type == amber_type, f"{chembl_id} atom#{atom_index}: amber={amber_type}, xponge={x_type}"


def test_gaff2_special_nitrogen_types():
    """
        test special gaff2 nitrogen typing
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2
    s = StringIO(r"""
@<TRIPOS>MOLECULE
MOL
   14    13     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          3.5400    1.4200    0.0000 C.3       1 MOL     -0.444095
      2 H1          2.5460    1.2490   -0.3830 H         1 MOL      0.214061
      3 H2          4.2190    0.6730   -0.3830 H         1 MOL      0.214061
      4 H3          3.5320    1.3950    1.0800 H         1 MOL      0.214061
      5 N1          4.0060    2.7710   -0.4390 N.4       1 MOL      0.063003
      6 H4          4.0060    2.7710   -1.4480 H         1 MOL      0.342733
      7 C2          3.0680    3.8500    0.0010 C.3       1 MOL     -0.444095
      8 H5          3.4180    4.7960   -0.3830 H         1 MOL      0.214061
      9 H6          2.0830    3.6360   -0.3830 H         1 MOL      0.214061
     10 H7          3.0510    3.8700    1.0800 H         1 MOL      0.214061
     11 C3          5.4100    3.0440    0.0010 C.3       1 MOL     -0.444095
     12 H8          5.7170    4.0040   -0.3830 H         1 MOL      0.214061
     13 H9          6.0550    2.2680   -0.3830 H         1 MOL      0.214061
     14 H10         5.4360    3.0480    1.0800 H         1 MOL      0.214061
@<TRIPOS>BOND
     1    2    1 1
     2    3    1 1
     3    4    1 1
     4    5    1 1
     5    6    5 1
     6    7    5 1
     7    8    7 1
     8    9    7 1
     9   10    7 1
    10   11    5 1
    11   12   11 1
    12   13   11 1
    13   14   11 1
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
""")
    me3nh = Xponge.get_assignment_from_mol2(s)
    me3nh.determine_atom_type("gaff2")
    assert me3nh.atom_types[4].name == "nx"
    assert me3nh.atom_types[5].name == "hn"


def test_gaff2_amide_and_thiocarbonyl_types():
    """
        test gaff2-specific ns and cs atom typing
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2

    nma = StringIO(r"""
@<TRIPOS>MOLECULE
NMA
    9     8     1     0     0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.3       1 NMA      0.000000
      2 H1         -0.3600   -0.5800    0.8900 H         1 NMA      0.000000
      3 H2         -0.3600   -0.5800   -0.8900 H         1 NMA      0.000000
      4 H3         -0.5400    0.9800    0.0000 H         1 NMA      0.000000
      5 C2          1.5200   -0.0400    0.0000 C.2       1 NMA      0.000000
      6 O1          2.1800   -1.0700    0.0000 O.2       1 NMA      0.000000
      7 N1          2.2500    1.0900    0.0000 N.am      1 NMA      0.000000
      8 H4          1.7600    1.9600    0.0000 H         1 NMA      0.000000
      9 C3          3.6900    1.1800    0.0000 C.3       1 NMA      0.000000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
     5     5     6 2
     6     5     7 1
     7     7     8 1
     8     7     9 1
@<TRIPOS>SUBSTRUCTURE
     1 NMA         1 TEMP              0 ****  ****    0 ROOT
""")
    nma = Xponge.get_assignment_from_mol2(nma)
    nma.determine_atom_type("gaff2")
    assert nma.atom_types[6].name == "ns"

    thiocarbonyl = StringIO(r"""
@<TRIPOS>MOLECULE
THI
    4     3     1     0     0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
      1 C1          0.0000    0.0000    0.0000 C.2       1 THI      0.000000
      2 S1          1.6200    0.0000    0.0000 S.2       1 THI      0.000000
      3 H1         -0.5800    0.9300    0.0000 H         1 THI      0.000000
      4 H2         -0.5800   -0.9300    0.0000 H         1 THI      0.000000
@<TRIPOS>BOND
     1     1     2 2
     2     1     3 1
     3     1     4 1
@<TRIPOS>SUBSTRUCTURE
     1 THI         1 TEMP              0 ****  ****    0 ROOT
""")
    thiocarbonyl = Xponge.get_assignment_from_mol2(thiocarbonyl)
    thiocarbonyl.determine_atom_type("gaff2")
    assert thiocarbonyl.atom_types[0].name == "cs"


def test_gaff2_ring_sp3_carbons():
    """
        test gaff2 c5/c6 typing for saturated five- and six-membered rings
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2

    cyclopentane = StringIO(r"""
@<TRIPOS>MOLECULE
CYP
   15    15     1     0     0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
      1 C1          1.2140    0.0000    0.0000 C.3       1 CYP      0.000000
      2 C2          0.3750    1.1540    0.0000 C.3       1 CYP      0.000000
      3 C3         -0.9820    0.7130    0.0000 C.3       1 CYP      0.000000
      4 C4         -0.9820   -0.7130    0.0000 C.3       1 CYP      0.000000
      5 C5          0.3750   -1.1540    0.0000 C.3       1 CYP      0.000000
      6 H1          2.2940    0.0000    0.0000 H         1 CYP      0.000000
      7 H2          0.7880    1.9970    0.0000 H         1 CYP      0.000000
      8 H3         -1.8510    1.3440    0.0000 H         1 CYP      0.000000
      9 H4         -1.8510   -1.3440    0.0000 H         1 CYP      0.000000
     10 H5          0.7880   -1.9970    0.0000 H         1 CYP      0.000000
     11 H6          1.2140    0.0000    1.0900 H         1 CYP      0.000000
     12 H7          0.3750    1.1540    1.0900 H         1 CYP      0.000000
     13 H8         -0.9820    0.7130    1.0900 H         1 CYP      0.000000
     14 H9         -0.9820   -0.7130    1.0900 H         1 CYP      0.000000
     15 H10         0.3750   -1.1540    1.0900 H         1 CYP      0.000000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 1
     3     3     4 1
     4     4     5 1
     5     5     1 1
     6     1     6 1
     7     2     7 1
     8     3     8 1
     9     4     9 1
    10     5    10 1
    11     1    11 1
    12     2    12 1
    13     3    13 1
    14     4    14 1
    15     5    15 1
@<TRIPOS>SUBSTRUCTURE
     1 CYP         1 TEMP              0 ****  ****    0 ROOT
""")
    cyclopentane = Xponge.get_assignment_from_mol2(cyclopentane)
    cyclopentane.determine_atom_type("gaff2")
    for i in range(5):
        assert cyclopentane.atom_types[i].name == "c5"

    cyclohexane = StringIO(r"""
@<TRIPOS>MOLECULE
CYH
   18    18     1     0     0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
      1 C1          1.2140    0.7010    0.0000 C.3       1 CYH      0.000000
      2 C2          0.0000    1.4020    0.0000 C.3       1 CYH      0.000000
      3 C3         -1.2140    0.7010    0.0000 C.3       1 CYH      0.000000
      4 C4         -1.2140   -0.7010    0.0000 C.3       1 CYH      0.000000
      5 C5          0.0000   -1.4020    0.0000 C.3       1 CYH      0.000000
      6 C6          1.2140   -0.7010    0.0000 C.3       1 CYH      0.000000
      7 H1          2.1570    1.2450    0.0000 H         1 CYH      0.000000
      8 H2          1.2140    0.7010    1.0900 H         1 CYH      0.000000
      9 H3          0.0000    2.4900    0.0000 H         1 CYH      0.000000
     10 H4          0.0000    1.4020    1.0900 H         1 CYH      0.000000
     11 H5         -2.1570    1.2450    0.0000 H         1 CYH      0.000000
     12 H6         -1.2140    0.7010    1.0900 H         1 CYH      0.000000
     13 H7         -2.1570   -1.2450    0.0000 H         1 CYH      0.000000
     14 H8         -1.2140   -0.7010    1.0900 H         1 CYH      0.000000
     15 H9          0.0000   -2.4900    0.0000 H         1 CYH      0.000000
     16 H10         0.0000   -1.4020    1.0900 H         1 CYH      0.000000
     17 H11         2.1570   -1.2450    0.0000 H         1 CYH      0.000000
     18 H12         1.2140   -0.7010    1.0900 H         1 CYH      0.000000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 1
     3     3     4 1
     4     4     5 1
     5     5     6 1
     6     6     1 1
     7     1     7 1
     8     1     8 1
     9     2     9 1
    10     2    10 1
    11     3    11 1
    12     3    12 1
    13     4    13 1
    14     4    14 1
    15     5    15 1
    16     5    16 1
    17     6    17 1
    18     6    18 1
@<TRIPOS>SUBSTRUCTURE
     1 CYH         1 TEMP              0 ****  ****    0 ROOT
""")
    cyclohexane = Xponge.get_assignment_from_mol2(cyclohexane)
    cyclohexane.determine_atom_type("gaff2")
    for i in range(6):
        assert cyclohexane.atom_types[i].name == "c6"


def test_gaff2_parmchk2(tmp_path):
    """
        test the parmchk2 interface with gaff2
    """
    import Xponge
    import Xponge.forcefield.amber.gaff2 as gaff2
    s = StringIO(r"""
@<TRIPOS>MOLECULE
MOL
   14    13     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          3.5400    1.4200    0.0000 C.3       1 MOL     -0.444095
      2 H1          2.5460    1.2490   -0.3830 H         1 MOL      0.214061
      3 H2          4.2190    0.6730   -0.3830 H         1 MOL      0.214061
      4 H3          3.5320    1.3950    1.0800 H         1 MOL      0.214061
      5 N1          4.0060    2.7710   -0.4390 N.4       1 MOL      0.063003
      6 H4          4.0060    2.7710   -1.4480 H         1 MOL      0.342733
      7 C2          3.0680    3.8500    0.0010 C.3       1 MOL     -0.444095
      8 H5          3.4180    4.7960   -0.3830 H         1 MOL      0.214061
      9 H6          2.0830    3.6360   -0.3830 H         1 MOL      0.214061
     10 H7          3.0510    3.8700    1.0800 H         1 MOL      0.214061
     11 C3          5.4100    3.0440    0.0010 C.3       1 MOL     -0.444095
     12 H8          5.7170    4.0040   -0.3830 H         1 MOL      0.214061
     13 H9          6.0550    2.2680   -0.3830 H         1 MOL      0.214061
     14 H10         5.4360    3.0480    1.0800 H         1 MOL      0.214061
@<TRIPOS>BOND
     1    2    1 1
     2    3    1 1
     3    4    1 1
     4    5    1 1
     5    6    5 1
     6    7    5 1
     7    8    7 1
     8    9    7 1
     9   10    7 1
    10   11    5 1
    11   12   11 1
    12   13   11 1
    13   14   11 1
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
    """)
    me3nh = Xponge.get_assignment_from_mol2(s)
    me3nh.determine_atom_type("gaff2")
    me3nh = me3nh.to_residuetype("ME3")
    frcmod = tmp_path / "me3nh_gaff2.frcmod"
    gaff2.parmchk2_gaff2(me3nh, str(frcmod), direct_load=False)
    content = frcmod.read_text()
    assert "c3-nx" in content
    assert "hn-nx" in content


def test_gaff_assignment_state_and_residuetype_regression():
    """
        regression:
        1) bond-order updates must invalidate built-state cache
        2) gaff/gaff2 typing should not degrade cc/cd to c2 on default path
        3) to_residuetype should preserve corrected atom types
    """
    import Xponge
    import Xponge.forcefield.amber.gaff  # noqa: F401
    import Xponge.forcefield.amber.gaff2  # noqa: F401
    s = StringIO(r"""
@<TRIPOS>MOLECULE
*****
 49 52 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C          -5.0279    0.1289   -3.6674 C.3     1  UNL1        0.0790
      2 O          -5.3788   -0.4045   -2.3945 O.3     1  UNL1       -0.4914
      3 C          -4.5173   -0.4542   -1.2868 C.ar    1  UNL1        0.1635
      4 C          -3.2000    0.0321   -1.3447 C.ar    1  UNL1        0.0099
      5 C          -2.3695   -0.0330   -0.2234 C.ar    1  UNL1        0.0816
      6 N          -1.1016    0.4431   -0.3066 N.ar    1  UNL1       -0.2141
      7 C          -0.2561    0.4015    0.7552 C.ar    1  UNL1        0.2201
      8 N           1.0848    0.9190    0.6348 N.pl3   1  UNL1       -0.2966
      9 C           2.0519    0.7649    1.7317 C.3     1  UNL1        0.0307
     10 C           3.4509    0.4867    1.1717 C.3     1  UNL1        0.0309
     11 N           3.8109    1.4396    0.1093 N.am    1  UNL1       -0.2920
     12 C           5.1991    1.6674   -0.1941 C.2     1  UNL1        0.2820
     13 O           5.6311    2.8484   -0.2818 O.2     1  UNL1       -0.2666
     14 C           6.1395    0.5544   -0.4037 C.ar    1  UNL1        0.1892
     15 C           5.8198   -0.7390   -0.7717 C.ar    1  UNL1       -0.0113
     16 C           7.0220   -1.4055   -0.8700 C.ar    1  UNL1       -0.0231
     17 C           7.9992   -0.4810   -0.5636 C.ar    1  UNL1        0.0928
     18 O           7.4427    0.6775   -0.2948 O.2     1  UNL1       -0.4583
     19 C           2.8100    2.3820   -0.4034 C.3     1  UNL1        0.0309
     20 C           1.4526    1.6893   -0.5640 C.3     1  UNL1        0.0307
     21 N          -0.6894   -0.1269    1.9299 N.ar    1  UNL1       -0.1992
     22 C          -1.9460   -0.6220    2.0847 C.ar    1  UNL1        0.1291
     23 N          -2.3033   -1.1577    3.3656 N.pl3   1  UNL1       -0.3425
     24 C          -2.8413   -0.5883    0.9865 C.ar    1  UNL1        0.0437
     25 C          -4.1632   -1.0751    1.0392 C.ar    1  UNL1       -0.0038
     26 C          -5.0002   -1.0089   -0.0932 C.ar    1  UNL1        0.1623
     27 O          -6.3216   -1.4849   -0.0730 O.3     1  UNL1       -0.4914
     28 C          -6.9550   -2.0733    1.0587 C.3     1  UNL1        0.0790
     29 H          -4.7609    1.2025   -3.5712 H       1  UNL1        0.0660
     30 H          -5.8976    0.0393   -4.3500 H       1  UNL1        0.0660
     31 H          -4.1776   -0.4404   -4.0988 H       1  UNL1        0.0660
     32 H          -2.8103    0.4633   -2.2575 H       1  UNL1        0.0677
     33 H           1.7741   -0.0780    2.4004 H       1  UNL1        0.0480
     34 H           2.0682    1.6966    2.3374 H       1  UNL1        0.0480
     35 H           3.4693   -0.5455    0.7672 H       1  UNL1        0.0480
     36 H           4.1850    0.5516    2.0045 H       1  UNL1        0.0480
     37 H           4.8408   -1.1449   -0.9847 H       1  UNL1        0.0657
     38 H           7.1705   -2.4403   -1.1478 H       1  UNL1        0.0649
     39 H           9.0649   -0.6656   -0.5476 H       1  UNL1        0.1029
     40 H           2.7137    3.2346    0.3036 H       1  UNL1        0.0480
     41 H           3.1255    2.7820   -1.3922 H       1  UNL1        0.0480
     42 H           1.4933    1.0028   -1.4374 H       1  UNL1        0.0480
     43 H           0.6833    2.4672   -0.7626 H       1  UNL1        0.0480
     44 H          -3.2383   -1.5553    3.5907 H       1  UNL1        0.1437
     45 H          -1.5965   -1.1556    4.1344 H       1  UNL1        0.1437
     46 H          -4.5564   -1.5055    1.9429 H       1  UNL1        0.0663
     47 H          -6.9926   -1.3457    1.8968 H       1  UNL1        0.0660
     48 H          -7.9929   -2.3547    0.7870 H       1  UNL1        0.0660
     49 H          -6.4089   -2.9890    1.3691 H       1  UNL1        0.0660
@<TRIPOS>BOND
     1     1     2    1
     2     2     3    1
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     6     7   ar
     7     7     8    1
     8     8     9    1
     9     9    10    1
    10    10    11    1
    11    11    12   am
    12    12    13    2
    13    12    14    1
    14    14    15   ar
    15    15    16   ar
    16    16    17   ar
    17    17    18   ar
    18    11    19    1
    19    19    20    1
    20     7    21   ar
    21    21    22   ar
    22    22    23    1
    23    22    24   ar
    24    24    25   ar
    25    25    26   ar
    26    26    27    1
    27    27    28    1
    28    26     3   ar
    29    24     5   ar
    30    20     8    1
    31    18    14   ar
    32     1    29    1
    33     1    30    1
    34     1    31    1
    35     4    32    1
    36     9    33    1
    37     9    34    1
    38    10    35    1
    39    10    36    1
    40    15    37    1
    41    16    38    1
    42    17    39    1
    43    19    40    1
    44    19    41    1
    45    20    42    1
    46    20    43    1
    47    23    44    1
    48    23    45    1
    49    25    46    1
    50    28    47    1
    51    28    48    1
    52    28    49    1
""")
    assign = Xponge.get_assignment_from_mol2(s)
    assert assign.built is False

    assign.determine_atom_type("gaff2")
    assert [assign.atom_types[i].name for i in [13, 14, 15, 16]] == ["cc", "cd", "cd", "cc"]

    assign.determine_atom_type("gaff")
    assert [assign.atom_types[i].name for i in [13, 14, 15, 16]] == ["cc", "cd", "cd", "cc"]

    restype = assign.to_residuetype("PRZ")
    assert [restype.atoms[i].type.name for i in [13, 14, 15, 16]] == ["cc", "cd", "cd", "cc"]

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

def test_assign_pdb_hybrid36_save(tmp_path):
    """
        test saving Assign as PDB with hybrid-36 atom serials and CONECT
    """
    import Xponge
    from Xponge.helper import Xdict
    atom_numbers = 100001
    assign = Xponge.Assign("H36")
    assign.atom_numbers = atom_numbers
    assign.atoms = ["C"] * atom_numbers
    assign.names = ["C"] * atom_numbers
    assign.element_details = [""] * atom_numbers
    assign.coordinate = np.zeros((atom_numbers, 3), dtype=np.float32)
    assign.charge = np.zeros(atom_numbers, dtype=np.float32)
    assign.formal_charge = [0] * atom_numbers
    assign.bonds = Xdict({i: Xdict() for i in range(atom_numbers)})
    assign.bonds[99999][100000] = 1
    assign.bonds[100000][99999] = 1

    outfile = tmp_path / "assign_h36.pdb"
    assign.save_as_pdb(str(outfile))
    lines = outfile.read_text().splitlines()
    atom_lines = [line for line in lines if line.startswith("ATOM")]
    conect_lines = [line for line in lines if line.startswith("CONECT")]

    assert atom_lines[99999][6:11] == "A0000"
    assert atom_lines[100000][6:11] == "A0001"
    assert any(line[6:11] == "A0000" and line[11:16] == "A0001" for line in conect_lines)

def test_assign_pdb_hybrid36_load():
    """
        test loading Assign from PDB with hybrid-36 atom serials and CONECT
    """
    import Xponge
    s = StringIO(
        "ATOM  A0000 C1   BEN     1       0.000   0.000   0.000                      C\n"
        "ATOM  A0001 C2   BEN     1       1.200   0.000   0.000                      C\n"
        "CONECTA0000A0001\n"
    )
    assign = Xponge.get_assignment_from_pdb(s, "BEN")
    assert assign.atom_numbers == 2
    assert 1 in assign.bonds[0]
    assert 0 in assign.bonds[1]
