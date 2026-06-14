"""
    This **module** includes lightweight unittests for MCPB scaffolding
"""


def _ensure_mcpb_test_types():
    import Xponge
    import Xponge.forcefield.base.mass_base
    import Xponge.forcefield.base.charge_base

    if "MCPB_ZN" not in Xponge.AtomType.get_all_types():
        Xponge.AtomType.New_From_String(r"""
name     mass   charge[e]
MCPB_ZN  65.38  2.0
MCPB_N   14.01 -0.5
""")

    if "MCPB_ZN_RES" not in Xponge.ResidueType.get_all_types():
        residue_type = Xponge.ResidueType(name="MCPB_ZN_RES")
        residue_type.add_atom("ZN", "MCPB_ZN", 0.0, 0.0, 0.0)

    if "MCPB_LIG_RES" not in Xponge.ResidueType.get_all_types():
        residue_type = Xponge.ResidueType(name="MCPB_LIG_RES")
        residue_type.add_atom("N1", "MCPB_N", 2.0, 0.0, 0.0)


def _build_mcpb_test_molecule():
    import Xponge

    _ensure_mcpb_test_types()

    mol = Xponge.Molecule(name="MCPB_TEST")
    mol.add_residue(Xponge.ResidueType.get_type("MCPB_ZN_RES"))
    mol.add_residue(Xponge.ResidueType.get_type("MCPB_LIG_RES"))
    zn_atom = mol.residues[0].atoms[0]
    lig_atom = mol.residues[1].atoms[0]
    zn_atom.element = "Zn"
    lig_atom.element = "N"
    zn_atom.charge = 2.0
    lig_atom.charge = -0.5
    return mol


def test_mcpb_selection_normalizes_origin_molecule():
    """
        test MCPB request normalization and selection on an origin Molecule
    """
    from Xponge.mcpb.selection import normalize_request, validate_and_select_environment

    mol = _build_mcpb_test_molecule()
    request = normalize_request(
        mol,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "zn", "formal_charge": 2, "spin": 0}],
        bonded_pairs=[(0, 1)],
        method="blank",
    )
    request, selection, warnings = validate_and_select_environment(request)

    assert request.ion_ids == (0,)
    assert selection.ion_atom_ids == (0,)
    assert selection.coordinating_atom_ids == (1,)
    assert selection.selected_residue_ids == (0, 1)
    assert warnings == []


def test_mcpb_blank_frcmod_uses_atom_type_names():
    """
        test blank frcmod generation uses clean atom type names on origin objects
    """
    from Xponge.mcpb.frcmod import build_blank_frcmod_text
    from Xponge.mcpb.selection import normalize_request, validate_and_select_environment

    mol = _build_mcpb_test_molecule()
    request = normalize_request(
        mol,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 0}],
        bonded_pairs=[(0, 1)],
        method="blank",
    )
    request, selection, _ = validate_and_select_environment(request)
    text = build_blank_frcmod_text(request, selection)

    assert "MCPB_N-MCPB_ZN" in text or "MCPB_ZN-MCPB_N" in text
    assert "Type of Atom" not in text


def test_mcpb_artifact_manifest_records_clean_charge_override_types(tmp_path):
    """
        test artifact manifest uses clean type names in charge override records
    """
    from Xponge.mcpb.export import write_mcpb_artifacts
    from Xponge.mcpb.models import MCPBRequest, MCPBResult, MCPBSelection

    mol = _build_mcpb_test_molecule()
    result = MCPBResult(
        molecule=mol,
        request=MCPBRequest(molecule=mol, ion_ids=(0,), ion_info=()),
        selection=MCPBSelection(ion_atom_ids=(0,), coordinating_atom_ids=(1,), selected_residue_ids=(0, 1), bonded_pairs=((0, 1),)),
        updated_charge_atoms=[0, 1],
    )

    manifest = write_mcpb_artifacts(result, tmp_path, write_sponge_input_files=False, write_local_models=False)

    assert manifest["pdb_path"].endswith(".pdb")
    assert manifest["charge_override_path"].endswith(".json")
    charge_overrides = (tmp_path / "mcpb.charge_overrides.json").read_text()
    assert "MCPB_ZN" in charge_overrides
    assert "Type of Atom" not in charge_overrides


def test_mcpb_model_builder_and_charge_refit_patch_parent_charges():
    """
        test local model extraction and patched charge refit on origin objects
    """
    from unittest.mock import patch
    import numpy as np
    from Xponge.mcpb.charge_refit import run_local_charge_refit
    from Xponge.mcpb.model_builder import build_small_and_large_models
    from Xponge.mcpb.selection import normalize_request, validate_and_select_environment

    mol = _build_mcpb_test_molecule()
    request = normalize_request(
        mol,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 0}],
        bonded_pairs=[(0, 1)],
        method="blank",
        basis="6-31g*",
    )
    request, selection, _ = validate_and_select_environment(request)
    small_model, large_model = build_small_and_large_models(request, selection)

    assert small_model.source_atom_ids == (0, 1)
    assert large_model.source_atom_ids == (0, 1)

    with patch("Xponge.mcpb.charge_refit.resp_module.resp_fit", return_value=np.array([1.25, -1.25])):
        summary = run_local_charge_refit(request, large_model)

    assert summary["updated_atom_ids"] == [0, 1]
    assert mol.atoms[0].charge == 1.25
    assert mol.atoms[1].charge == -1.25


def test_mcpb_api_returns_blank_scaffold_result():
    """
        test top-level MCPB API returns a structured scaffold result in blank mode
    """
    import Xponge

    mol = _build_mcpb_test_molecule()
    result = Xponge.MCPB(
        mol,
        ion_ids=[0],
        ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 0}],
        bonded_pairs=[(0, 1)],
        method="blank",
        charge_mode="fixed",
    )

    assert result.selection.ion_atom_ids == (0,)
    assert result.connect_records == [(0, 1)]
    assert result.frcmod_path is not None
    assert result.metadata["frcmod_generated"] is True
    assert any("fixed-charge MCPB mode selected" in item for item in result.pending_requirements)


def test_mcpb_api_runs_local_resp_refit_when_basis_is_provided():
    """
        test top-level MCPB API executes local RESP refit through the new QM plumbing
    """
    from unittest.mock import patch
    import numpy as np
    import Xponge

    mol = _build_mcpb_test_molecule()
    with patch("Xponge.mcpb.charge_refit.resp_module.resp_fit", return_value=np.array([0.75, -0.75])):
        result = Xponge.MCPB(
            mol,
            ion_ids=[0],
            ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 0}],
            bonded_pairs=[(0, 1)],
            method="blank",
            basis="6-31g*",
        )

    assert result.updated_charge_atoms == [0, 1]
    assert mol.atoms[0].charge == 0.75
    assert mol.atoms[1].charge == -0.75
    assert result.metadata["charge_refit"]["total_charge"] == 2


def test_mcpb_api_runs_mocked_seminario_path():
    """
        test top-level MCPB API can drive the Seminario path through a mocked Hessian backend
    """
    from unittest.mock import patch
    import numpy as np
    import Xponge

    mol = _build_mcpb_test_molecule()
    mocked_hessian = Xponge.qm.HessianResult(
        cartesian_hessian_au=np.eye(6),
        coordinates_angstrom=[(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        atom_symbols=["Zn", "N"],
        timings={"hessian": 0.0},
    )

    with patch("Xponge.mcpb.seminario.compute_hessian", return_value=mocked_hessian):
        result = Xponge.MCPB(
            mol,
            ion_ids=[0],
            ion_info=[{"atom_id": 0, "element": "Zn", "formal_charge": 2, "spin": 0}],
            bonded_pairs=[(0, 1)],
            method="seminario",
            basis="6-31g*",
            charge_mode="fixed",
        )

    assert result.frcmod_path is not None
    assert result.metadata["frcmod_generated"] is True
    assert result.pending_requirements[0].startswith("fixed-charge MCPB mode selected")
