"""
    This **module** includes unittests of the charge calculations
"""

__all__ = ["test_tpacm4"]

from io import StringIO
from types import SimpleNamespace

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


def test_resp_parameters_select_standard_basis_and_mk_radius():
    from Xponge.qm.resp_parameters import get_resp_mk_radius, select_resp_basis_family

    assert select_resp_basis_family({"C", "H", "O"}) == "6-31G*"
    assert select_resp_basis_family({"C", "H", "I"}) == "CEP-31G"
    assert select_resp_basis_family({"C", "H", "U"}) == "SDD"
    assert get_resp_mk_radius("I") == pytest.approx(1.99)
    assert get_resp_mk_radius("Be") == pytest.approx(1.7)


def test_resp_basis_resolver_maps_iodine_cep31g_to_sbkjc():
    from Xponge.qm.resp_basis import resolve_default_resp_basis

    resolved = resolve_default_resp_basis({"C", "H", "I"})

    assert resolved.label == "CEP-31G"
    assert resolved.basis == {"H": "dz", "C": "sbkjc", "I": "sbkjc"}
    assert resolved.ecp == {"C": "sbkjc", "I": "sbkjc"}
    assert resolved.cart is True


def test_resp_basis_resolver_keeps_sdd_separate_from_cep31g():
    from Xponge.qm.resp_basis import resolve_resp_basis

    cep = resolve_resp_basis("CEP-31G", {"I"})
    sdd = resolve_resp_basis("SDD", {"I"})

    assert cep.basis == {"I": "sbkjc"}
    assert cep.ecp == {"I": "sbkjc"}
    assert sdd.basis == {"I": "stuttgart"}
    assert sdd.ecp == {"I": "stuttgart"}


def test_resp_default_setup_passes_resolved_basis_ecp_cart_and_mk_radii(monkeypatch):
    import numpy as np
    import Xponge
    from Xponge.assign import resp

    assign = Xponge.Assign("METHYL_IODIDE")
    for index, element in enumerate(["C", "H", "H", "H", "I"]):
        assign.add_atom(element, float(index), 0.0, 0.0, name=f"{element}{index + 1}")
    assign.add_bond(0, 1, 1)
    assign.add_bond(0, 2, 1)
    assign.add_bond(0, 3, 1)
    assign.add_bond(0, 4, 1)
    assign.charge = np.zeros(5)
    calls = {}

    class FakeSCF:
        backend_name = "fake"
        coordinates_bohr = np.zeros((5, 3))
        nuclear_charges = np.array([6.0, 1.0, 1.0, 1.0, 53.0])

    class FakeESP:
        electronic_esp_au = np.zeros(2)

    def fake_run_scf(assign, *, backend=None, basis=None, ecp=None, cart=None, charge=0, spin=0, optimize_geometry=False):
        calls["scf"] = (basis, ecp, cart, charge, spin, optimize_geometry)
        return FakeSCF()

    def fake_grid(assign, atom_coordinates_bohr, area_density=1.0, layer=4, radius=None):
        calls["radius"] = dict(radius)
        return np.zeros((2, 3))

    monkeypatch.setattr(resp.qm_scheduler, "run_scf", fake_run_scf)
    monkeypatch.setattr(resp.qm_scheduler, "compute_esp_on_grid", lambda *args, **kwargs: FakeESP())
    monkeypatch.setattr(resp.resp_core, "get_mk_grid", fake_grid)
    monkeypatch.setattr(resp.resp_core, "fit_resp_from_esp", lambda *args, **kwargs: np.zeros(5))

    result = resp.resp_fit(assign, backend="pyscf", charge=0, return_metadata=True)

    assert calls["scf"] == (
        {"C": "sbkjc", "H": "dz", "I": "sbkjc"},
        {"C": "sbkjc", "I": "sbkjc"},
        True,
        0,
        0,
        False,
    )
    assert calls["radius"]["I"] == pytest.approx(1.99)
    assert calls["radius"]["C"] == pytest.approx(1.61)
    assert result["charges"].shape == (5,)
    assert result["metadata"]["basis_family"] == "CEP-31G"


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
        Test the RESP SCF kernel still accepts legacy column-vector charges
    """
    import numpy as np
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


def test_qm_scheduler_default_backend_matches_platform(monkeypatch):
    """
        Test the default QM backend selection matches the current platform policy
    """
    import Xponge
    from Xponge.qm import scheduler as qm_scheduler

    monkeypatch.setattr(qm_scheduler.sys, "platform", "linux")
    assert qm_scheduler.normalize_backend_name(None) == "pyscf"
    assert qm_scheduler.get_backend(None).name == "pyscf"

    monkeypatch.setattr(qm_scheduler.sys, "platform", "win32")
    assert qm_scheduler.normalize_backend_name(None) == "psi4"
    assert qm_scheduler.get_backend(None).name == "psi4"


def test_qm_scheduler_windows_import_hint_mentions_external_psi4_install(monkeypatch):
    """
        Test the Windows Psi4 import hint points users to external installation.
    """
    from Xponge.qm import scheduler as qm_scheduler

    monkeypatch.setattr(qm_scheduler.sys, "platform", "win32")

    with pytest.raises(qm_scheduler.QMBackendImportError, match="conda-forge"):
        qm_scheduler.backend_import_or_hint("psi4", ImportError("Psi4 is required"))


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
