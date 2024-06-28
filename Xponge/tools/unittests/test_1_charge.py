"""
    This **module** includes unittests of the charge calculations
"""

__all__ = ["test_tpacm4"]

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
