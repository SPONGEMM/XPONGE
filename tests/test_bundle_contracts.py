from Xponge.io_bundle.contracts import CONTRACTS, RUN_MDIN_KEYS


def test_run_mdin_inventory_covers_typed_barostat_controls():
    expected = {
        "barostat_seed",
        "barostat_isotropy",
        "barostat_target_surface_tensor",
        "barostat_compressibility",
        "barostat_g11",
        "barostat_g21",
        "barostat_g22",
        "barostat_g31",
        "barostat_g32",
        "barostat_g33",
        "barostat_tau",
        "barostat_update_interval",
        "barostat_x_constant",
        "barostat_y_constant",
        "barostat_z_constant",
        "monte_carlo_barostat_seed",
    }

    assert expected <= set(RUN_MDIN_KEYS)
    contract_ids = {contract.contract_id for contract in CONTRACTS}
    assert {f"run_mdin.{key}" for key in expected} <= contract_ids
