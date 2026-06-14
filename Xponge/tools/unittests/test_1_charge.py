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


def test_resp_scf_kernel_accepts_column_vector_charges():
    """
        Test the RESP SCF kernel still accepts legacy column-vector charges
    """
    import numpy as np
    from types import SimpleNamespace
    from Xponge.assign import resp

    mol = SimpleNamespace(natm=2)
    assign = SimpleNamespace(atoms=["C", "H"])
    matrix_a0 = np.array([[2.0, 0.0, 1.0],
                          [0.0, 2.0, 1.0],
                          [1.0, 1.0, 0.0]])
    matrix_a = matrix_a0.copy()
    matrix_b = np.array([[0.2], [-0.2], [0.0]])
    q = np.array([[0.0], [0.0]])

    result = resp._resp_scf_kernel(mol, assign, 0.0005, 0.1, matrix_a, matrix_a0, matrix_b, q)

    assert result.shape == (2,)


def test_qm_capabilities_report_psi4_hessian_as_unsupported():
    """
        Test the QM capability declaration for the optional Psi4 backend
    """
    import Xponge

    caps = Xponge.get_qm_capabilities("psi4")
    assert caps.supports_scf
    assert caps.supports_esp
    assert caps.supports_geometry_optimization
    assert not caps.supports_hessian


def test_calculate_charge_resp_passes_backend_and_core():
    """
        Test the RESP entrypoint forwards backend and core options
    """
    from unittest.mock import patch
    import numpy as np
    import Xponge

    assign = Xponge.Assign("RESP_OPTIONS")
    assign.add_atom("O", 0.0, 0.0, 0.0, name="O1")
    assign.charge = np.array([0.0])
    assign.formal_charge[0] = 0

    with patch("Xponge.assign.resp.RESP_Fit", return_value=np.array([0.0])) as mock_resp:
        assign.calculate_charge("RESP", backend="psi4", core="python")

    assert assign.charge.shape == (1,)
    assert mock_resp.call_args.kwargs["backend"] == "psi4"
    assert mock_resp.call_args.kwargs["core"] == "python"
