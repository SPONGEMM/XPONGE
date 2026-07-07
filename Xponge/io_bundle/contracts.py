"""
Contract metadata for SPONGE legacy-to-bundle conversion.

The table is intentionally explicit: converter output manifests should make it
clear which legacy keys were consumed, where their bundled representation lives,
and which contracts still need implementation.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class IOContract:
    """A single legacy-to-bundle mapping contract."""

    contract_id: str
    legacy_keys: tuple[str, ...]
    bundle_file: str
    bundle_path: str
    direction: str
    modes: tuple[str, ...]
    component: str
    override_policy: str
    comparison_rule: str
    status: str = "supported"
    payload_kind: str = "file"


def _contract(
    contract_id: str,
    legacy_key: str,
    bundle_file: str,
    bundle_path: str,
    component: str,
    *,
    direction: str = "input",
    modes: tuple[str, ...] = ("normal", "rerun"),
    override_policy: str = "allowed",
    comparison_rule: str = "text",
    status: str = "supported",
    payload_kind: str = "file",
) -> IOContract:
    return IOContract(
        contract_id=contract_id,
        legacy_keys=(legacy_key,),
        bundle_file=bundle_file,
        bundle_path=bundle_path,
        direction=direction,
        modes=modes,
        component=component,
        override_policy=override_policy,
        comparison_rule=comparison_rule,
        status=status,
        payload_kind=payload_kind,
    )


RUN_MDIN_KEYS = (
    "input_h5_topology_path",
    "input_h5_protocol_path",
    "input_h5_restart_path",
    "input_h5_restart_load",
    "input_h5_trajectory_path",
    "input_h5_trajectory_particle_stream",
    "workspace",
    "buffer_frame",
    "device",
    "device_optimized_block",
    "dont_check_input",
    "end_pause",
    "plugin",
    "default_in_file_prefix",
    "default_out_file_prefix",
    "mode",
    "dt",
    "step_limit",
    "frame_limit",
    "rerun_frame_limit",
    "target_temperature",
    "target_pressure",
    "target_temperature_schedule_mode",
    "target_temperature_schedule_steps",
    "target_temperature_schedule_file",
    "target_pressure_schedule_mode",
    "target_pressure_schedule_steps",
    "target_pressure_schedule_file",
    "pbc",
    "skin",
    "cutoff",
    "velocity_max",
    "make_output_whole",
    "force_whole_output",
    "write_information_interval",
    "write_trajectory_interval",
    "write_mdout_interval",
    "write_restart_file_interval",
    "print_pressure",
    "print_zeroth_frame",
    "max_restart_export_count",
    "rerun_start",
    "rerun_strip",
    "rerun_need_box_update",
    "thermostat",
    "thermostat_mode",
    "thermostat_tau",
    "thermostat_seed",
    "barostat",
    "barostat_mode",
    "nose_hoover_chain_length",
    "monte_carlo_barostat_update_interval",
    "monte_carlo_barostat_check_interval",
    "monte_carlo_barostat_initial_ratio",
    "monte_carlo_barostat_accept_rate_low",
    "monte_carlo_barostat_accept_rate_high",
    "monte_carlo_barostat_couple_dimension",
    "monte_carlo_barostat_only_direction",
    "monte_carlo_barostat_surface_number",
    "monte_carlo_barostat_surface_tension",
    "DOM_DEC_update_interval",
    "DOM_DEC_split_nx",
    "DOM_DEC_split_ny",
    "DOM_DEC_split_nz",
    "neighbor_list_refresh_interval",
    "neighbor_list_skin_permit",
    "neighbor_list_throw_error_when_overflow",
    "neighbor_list_max_neighbor_numbers",
    "neighbor_list_check_overflow_interval",
    "neighbor_list_max_atom_in_grid_numbers",
    "neighbor_list_max_ghost_in_grid_numbers",
    "minimization_max_move",
    "minimization_momentum_keep",
    "minimization_dynamic_dt",
    "minimization_beta1",
    "minimization_beta2",
    "minimization_epsilon",
    "constrain_mode",
    "constrain_angle",
    "constrain_mass",
    "SHAKE_iteration_numbers",
    "SHAKE_step_length",
    "settle_disable",
    "lambda_lj",
    "soft_core_alpha",
    "soft_core_powfer",
    "soft_core_sigma",
    "soft_core_sigma_min",
    "PM_fftx",
    "PM_ffty",
    "PM_fftz",
    "PM_grid_spacing",
    "PM_Direct_Tolerance",
    "PM_MPI_size",
    "PM_print_detail",
    "gb_epsilon",
    "gb_radii_cutoff",
    "gb_radii_offset",
    "hard_wall_x_low",
    "hard_wall_y_low",
    "hard_wall_z_low",
    "hard_wall_x_high",
    "hard_wall_y_high",
    "hard_wall_z_high",
)

H5_OUTPUT_KEYS = (
    "output_h5_trajectory_path",
    "output_h5_trajectory_vds",
    "output_h5_trajectory_chunk_size",
    "output_h5_trajectory_repair_policy",
    "output_h5_restart_path",
    "output_h5_observable_path",
)

LEGACY_OUTPUT_KEYS = (
    "mdinfo",
    "mdout",
    "crd",
    "box",
    "vel",
    "frc",
    "rst",
    "qc_scf_output",
    "nose_hoover_chain_crd",
    "nose_hoover_chain_vel",
    "nose_hoover_chain_restart_output",
    "meta_potential_out_file",
    "meta_direct_export",
    "meta_hills_log",
    "meta_history_log",
    "meta_edge_log",
    "SITS_nk_traj_file",
    "SITS_nk_rest_file",
)

TOPOLOGY_FILE_CONTRACTS = {
    "mass_in_file": "/atoms/mass",
    "charge_in_file": "/atoms/charge",
    "residue_in_file": "/residues",
    "exclude_in_file": "/topology/exclusions",
    "bond_in_file": "/forcefield/bond",
    "angle_in_file": "/forcefield/angle",
    "dihedral_in_file": "/forcefield/dihedral",
    "improper_in_file": "/forcefield/improper",
    "improper_dihedral_in_file": "/forcefield/improper",
    "LJ_in_file": "/forcefield/lj",
    "nb14_in_file": "/forcefield/nb14",
    "nb14_extra_in_file": "/forcefield/nb14_extra",
    "urey_bradley_in_file": "/forcefield/urey_bradley",
    "cmap_in_file": "/forcefield/cmap",
    "gb_in_file": "/forcefield/gb",
    "virtual_atom_in_file": "/forcefield/virtual_atom",
    "virtual_atoms_in_file": "/forcefield/virtual_atom",
    "LJ_soft_core_in_file": "/forcefield/lj_soft_core",
    "subsys_division_in_file": "/forcefield/subsys_division",
    "EAM_in_file": "/manybody/eam",
    "EAM_atom_type_in_file": "/manybody/eam/atom_type",
    "SW_in_file": "/manybody/sw",
    "EDIP_in_file": "/manybody/edip",
    "TERSOFF_in_file": "/manybody/tersoff",
    "REAXFF_in_file": "/manybody/reaxff/parameters",
    "REAXFF_type_in_file": "/manybody/reaxff/type",
    "qc_type_in_file": "/qc/type",
    "pairwise_force_in_file": "/forcefield/custom_force/pairwise",
    "listed_forces_in_file": "/forcefield/custom_force/listed",
}

PROTOCOL_FILE_CONTRACTS = {
    "cv_in_file": "/cv",
    "constrain_in_file": "/constraint/default/pairs",
    "restrain_in_file": "/restraint",
    "restrain_atom_id": "/restraint/default/atom_indices",
    "restrain_weight_in_file": "/restraint/default/weight",
    "restrain_cv_in_file": "/restraint/cv",
    "steer_cv_in_file": "/steer",
    "soft_walls_in_file": "/wall/soft",
    "SITS_in_file": "/sits",
    "SITS_atom_in_file": "/sits/atom_indices",
    "meta_edge_in_file": "/meta/default/grid",
}

RESTART_FILE_CONTRACTS = {
    "nose_hoover_chain_restart_input": "/parameters/restart/thermostat/nose_hoover_chain",
    "restrain_coordinate_in_file": "/parameters/restart/references/restraint/default/coordinate",
    "restrain_amber_rst7": "/parameters/restart/references/restraint/default/coordinate",
    "SITS_nk_in_file": "/parameters/restart/bias/sits/SITS/nk",
    "meta_potential_in_file": "/parameters/restart/bias/meta/default/potential_export",
    "meta_scatter_in_file": "/parameters/restart/bias/meta/default/scatter",
    "hill_in_file": "/parameters/restart/bias/meta/default/hills",
    "hills_in_file": "/parameters/restart/bias/meta/default/hills",
    "metad_hills_in_file": "/parameters/restart/bias/meta/default/hills",
}

PROTOCOL_RESTART_SIDECAR_KEYS = (
    "cv_in_file",
    "constrain_in_file",
    "restrain_in_file",
    "pairwise_force_in_file",
    "listed_forces_in_file",
    "soft_walls_in_file",
    "SITS_in_file",
    "SITS_atom_in_file",
    "SITS_nk_in_file",
    "restrain_atom_id",
    "restrain_coordinate_in_file",
    "restrain_weight_in_file",
    "meta_edge_in_file",
    "meta_potential_in_file",
    "meta_scatter_in_file",
    "restrain_cv_in_file",
    "steer_cv_in_file",
)

TRAJECTORY_INPUT_CONTRACTS = {
    "crd": "/particles/all/position/value",
    "box": "/particles/all/box/edges/value",
    "vel": "/particles/all/velocity/value",
}


CONTRACTS: tuple[IOContract, ...] = (
    IOContract(
        contract_id="restart.coordinate",
        legacy_keys=("coordinate_in_file",),
        bundle_file="restart.spgr.h5",
        bundle_path="/particles/all/position/value",
        direction="input",
        modes=("normal", "rerun"),
        component="base",
        override_policy="forbidden",
        comparison_rule="float32",
    ),
    IOContract(
        contract_id="restart.velocity",
        legacy_keys=("velocity_in_file",),
        bundle_file="restart.spgr.h5",
        bundle_path="/particles/all/velocity/value",
        direction="input",
        modes=("normal", "rerun"),
        component="base",
        override_policy="forbidden",
        comparison_rule="float32",
    ),
    IOContract(
        contract_id="restart.box",
        legacy_keys=("coordinate_in_file",),
        bundle_file="restart.spgr.h5",
        bundle_path="/particles/all/box/edges/value",
        direction="input",
        modes=("normal", "rerun"),
        component="base",
        override_policy="forbidden",
        comparison_rule="float32",
    ),
    IOContract(
        contract_id="rerun.trajectory",
        legacy_keys=("crd",),
        bundle_file="trajectory.spg.h5md",
        bundle_path="/particles/all/position/value",
        direction="input",
        modes=("rerun",),
        component="base",
        override_policy="forbidden",
        comparison_rule="float32",
        status="unsupported",
    ),
) + tuple(
    _contract(
        f"restart.protocol_sidecar.{key}",
        key,
        "restart.spgr.h5",
        f"/parameters/restart/protocol_sidecars/{key}",
        "restart",
        override_policy="legacy_protocol_sidecar",
        comparison_rule="text",
        status="protocol_restart_sidecar",
        payload_kind="file",
    )
    for key in PROTOCOL_RESTART_SIDECAR_KEYS
) + tuple(
    _contract(
        f"run_mdin.{key}",
        key,
        "run.mdin",
        key,
        "run_policy",
        direction="input",
        override_policy="explicit",
        comparison_rule="exact",
        status="run_mdin",
        payload_kind="scalar",
    )
    for key in RUN_MDIN_KEYS
) + tuple(
    _contract(
        f"output.h5.{key}",
        key,
        "run.mdin",
        key,
        "output",
        direction="output",
        override_policy="explicit",
        comparison_rule="exact",
        status="output_plan",
        payload_kind="scalar",
    )
    for key in H5_OUTPUT_KEYS
) + tuple(
    _contract(
        f"output.legacy_sidecar.{key}",
        key,
        "*.legacy",
        f"/parameters/sponge/files/{key}",
        "output",
        direction="output",
        override_policy="explicit",
        comparison_rule="path",
        status="legacy_output_sidecar",
        payload_kind="path",
    )
    for key in LEGACY_OUTPUT_KEYS
) + tuple(
    _contract(
        f"topology.{key.removesuffix('_in_file')}",
        key,
        "topology.spgt.h5",
        path,
        "topology",
        override_policy="forbidden",
        comparison_rule="typed",
        status="compatibility_import",
        payload_kind="file",
    )
    for key, path in TOPOLOGY_FILE_CONTRACTS.items()
) + tuple(
    _contract(
        f"protocol.{key.removesuffix('_in_file')}",
        key,
        "protocol.spgp.h5",
        path,
        "protocol",
        override_policy="allowed",
        comparison_rule="typed",
        status="compatibility_import",
        payload_kind="file",
    )
    for key, path in PROTOCOL_FILE_CONTRACTS.items()
) + tuple(
    _contract(
        f"restart.{key.removesuffix('_in_file')}",
        key,
        "restart.spgr.h5",
        path,
        "restart",
        override_policy="allowed",
        comparison_rule="typed",
        status="compatibility_import",
        payload_kind="file",
    )
    for key, path in RESTART_FILE_CONTRACTS.items()
) + tuple(
    _contract(
        f"trajectory.{key}",
        key,
        "trajectory.spg.h5md",
        path,
        "trajectory",
        modes=("rerun",),
        override_policy="forbidden",
        comparison_rule="typed",
        status="compatibility_import",
        payload_kind="file",
    )
    for key, path in TRAJECTORY_INPUT_CONTRACTS.items()
)


def contracts_by_legacy_key() -> dict[str, tuple[IOContract, ...]]:
    """Return all known contracts indexed by legacy mdin key."""

    index: dict[str, list[IOContract]] = {}
    for contract in CONTRACTS:
        for key in contract.legacy_keys:
            index.setdefault(key, []).append(contract)
    return {key: tuple(value) for key, value in index.items()}
