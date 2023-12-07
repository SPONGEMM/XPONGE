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
    assign = Xponge.get_assignment_from_smiles("OC1=C(C(O)=O)C=C(N)C=C1")
    assign.calculate_charge("tpacm4")
