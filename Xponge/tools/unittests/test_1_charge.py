"""
    This **module** includes unittests of the charge calculations
"""

__all__ = ["test_tpacm4"]

import pytest


def _resp_smoke_assignment():
    import numpy as np
    import Xponge

    assign = Xponge.Assign("RESP_SMOKE")
    assign.add_atom("O", 0.0, 0.0, 0.0, name="O1")
    assign.add_atom("H", 0.0, 0.757, 0.586, name="H1")
    assign.add_atom("H", 0.0, -0.757, 0.586, name="H2")
    assign.add_bond(0, 1, 1)
    assign.add_bond(0, 2, 1)
    assign.charge = np.zeros(3)
    return assign

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
        assign.calculate_charge(
            "RESP",
            backend="psi4",
            core="python",
            esp_memory_limit="256MB",
            esp_chunk_policy="grid",
            esp_safety_factor=0.5,
        )

    assert assign.charge.shape == (1,)
    assert mock_resp.call_args.kwargs["backend"] == "psi4"
    assert mock_resp.call_args.kwargs["core"] == "python"
    assert mock_resp.call_args.kwargs["esp_memory_limit"] == "256MB"
    assert mock_resp.call_args.kwargs["esp_chunk_policy"] == "grid"
    assert mock_resp.call_args.kwargs["esp_safety_factor"] == 0.5


def test_qm_scheduler_normalizes_esp_request_options():
    """
        Test the QM scheduler parses ESP request options before backend dispatch
    """
    import numpy as np
    from Xponge.qm import scheduler as qm_scheduler
    from Xponge.qm.capabilities import QMCapabilitySet
    from Xponge.qm.models import ESPResult, SCFResult

    requests = []

    class FakeESPBackend:
        name = "fakeesp"

        @staticmethod
        def capabilities():
            return QMCapabilitySet(supports_scf=False, supports_esp=True)

        @staticmethod
        def compute_esp(scf_result, request):
            requests.append(request)
            return ESPResult(
                grid_points_bohr=np.asarray(request.grid_points_bohr, dtype=float),
                electronic_esp_au=np.zeros(len(request.grid_points_bohr), dtype=float),
            )

    qm_scheduler._BACKENDS["fakeesp"] = FakeESPBackend
    try:
        scf_result = SCFResult(
            backend_name="fakeesp",
            total_energy=None,
            converged=True,
            coordinates_bohr=np.zeros((0, 3), dtype=float),
            nuclear_charges=np.zeros(0, dtype=float),
            charge=0,
            spin=0,
            atom_symbols=[],
        )
        result = qm_scheduler.compute_esp_on_grid(
            scf_result,
            np.zeros((3, 3), dtype=float),
            memory_limit="512MB",
            chunk_policy="grid",
            safety_factor=0.5,
        )
    finally:
        del qm_scheduler._BACKENDS["fakeesp"]

    assert len(result.electronic_esp_au) == 3
    assert requests[0].memory_limit_bytes == 512 * 1024 * 1024
    assert requests[0].chunk_policy == "grid"
    assert requests[0].safety_factor == pytest.approx(0.5)


def test_qm_scheduler_runs_pyscf_esp_grid_chunk_mode_smoke():
    """
        Test the PySCF backend can execute explicit grid-chunk ESP evaluation
    """
    pytest.importorskip("pyscf")
    from Xponge.qm import compute_esp_on_grid, run_scf

    assign = _resp_smoke_assignment()
    scf_result = run_scf(assign, backend="pyscf", basis="sto-3g", charge=0, spin=0)
    per_grid_bytes = scf_result.backend_handle["mol"].nao_nr() ** 2 * 8
    esp_result = compute_esp_on_grid(
        scf_result,
        scf_result.coordinates_bohr[:2],
        memory_limit=per_grid_bytes + 1,
        chunk_policy="grid",
        safety_factor=1.0,
    )

    assert len(esp_result.electronic_esp_au) == 2
    assert esp_result.diagnostics["mode"] == "grid_chunk"
    assert esp_result.diagnostics["grid_chunk_size"] == 1


def test_qm_scheduler_runs_pyscf_esp_dual_chunk_mode_smoke():
    """
        Test the PySCF backend can execute dual shell/grid chunk ESP evaluation
    """
    pytest.importorskip("pyscf")
    from Xponge.qm import compute_esp_on_grid, run_scf

    assign = _resp_smoke_assignment()
    scf_result = run_scf(assign, backend="pyscf", basis="sto-3g", charge=0, spin=0)
    per_grid_bytes = scf_result.backend_handle["mol"].nao_nr() ** 2 * 8
    esp_result = compute_esp_on_grid(
        scf_result,
        scf_result.coordinates_bohr[:2],
        memory_limit=max(1, per_grid_bytes - 1),
        chunk_policy="auto",
        safety_factor=1.0,
    )

    assert len(esp_result.electronic_esp_au) == 2
    assert esp_result.diagnostics["mode"] == "shell_grid_chunk"
    assert esp_result.diagnostics["shell_block_count"] >= 1
