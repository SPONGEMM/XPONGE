"""
    This **module** includes unittests of the charge calculations
"""

__all__ = ["test_tpacm4"]

from io import StringIO
from types import SimpleNamespace

import pytest


def test_tpacm4():
    """
        Test the functions to calculate the tpacm4 charge
    """
    import Xponge
    assign = Xponge.get_assignment_from_smiles("c1ccccc1")
    assign.calculate_charge("tpacm4")
    Xponge.Xprint(assign.charge)
    assert abs(assign.charge[0] + 0.155) < 0.01
    assert abs(assign.charge[-1] - 0.155) < 0.01 
    assign = Xponge.get_assignment_from_smiles("OC1=C(C(O)=O)C=C(N)C=C1")
    assign.calculate_charge("tpacm4")
    Xponge.Xprint(assign.charge)


def test_resp_second_stage_hydrogen_groups_are_not_restrained():
    """
        Regression:
        second-stage RESP groups must be classified by their grouped atoms instead of raw atom indexes
    """
    pytest.importorskip("pyscf")

    import Xponge
    from Xponge.assign import resp

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

    assign = Xponge.get_assignment_from_mol2(cyclohexane)
    tofit_second, fit_group, sublength = resp._find_tofit_second(SimpleNamespace(natm=assign.atom_numbers), assign)
    restrained = resp._find_restrained_second_stage_groups(assign, tofit_second)

    assert restrained == [0, 2, 4, 6, 8, 10]
    assert all(any(assign.atoms[j] != "H" for j in tofit_second[i]) for i in restrained)
    assert all(all(assign.atoms[j] == "H" for j in tofit_second[i]) for i in range(len(tofit_second)) if i not in restrained)


def test_resp_scf_kernel_accepts_column_vector_charges():
    """
        Regression:
        RESP iterations should accept legacy ``(n, 1)`` charge vectors even on newer numpy
    """
    pytest.importorskip("pyscf")

    import numpy as np
    from Xponge.assign import resp

    mol = SimpleNamespace(natm=2)
    assign = SimpleNamespace(atoms=["C", "H"])
    matrix_a = np.eye(3)
    matrix_a0 = np.eye(3)
    matrix_b = np.array([[0.2], [0.0], [0.0]])
    q = resp._resp_scf_kernel(mol, assign, 0.0005, 0.1, matrix_a, matrix_a0, matrix_b, np.array([[0.1], [0.0]]))

    assert q.shape == (2,)
