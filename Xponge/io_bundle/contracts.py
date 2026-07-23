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
    """A bidirectional mapping contract for one SPONGE I/O binding."""

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
    parser_id: str | None = None
    exporter_id: str | None = None
    required_bundle_paths: tuple[str, ...] = ()
    legacy_filename_stem: str | None = None
    legacy_section: str | None = None
    reverse_policy: str = "not_reversible"
    aliases: tuple[str, ...] = ()
    materialization_group: str | None = None


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
    parser_id: str | None = None,
    exporter_id: str | None = None,
    required_bundle_paths: tuple[str, ...] | None = None,
    legacy_filename_stem: str | None = None,
    legacy_section: str | None = None,
    reverse_policy: str | None = None,
    aliases: tuple[str, ...] = (),
    materialization_group: str | None = None,
) -> IOContract:
    if payload_kind == "file" and direction == "input":
        if parser_id is None and status == "compatibility_import":
            parser_id = legacy_key
        if exporter_id is None and legacy_key not in LEGACY_KEY_ALIASES and status == "compatibility_import":
            exporter_id = LEGACY_KEY_EXPORTER_IDS.get(legacy_key, legacy_key)
        if required_bundle_paths is None and status == "compatibility_import":
            required_bundle_paths = (bundle_path,)
        if legacy_filename_stem is None:
            legacy_filename_stem = _legacy_filename_stem(legacy_key)
        if legacy_section is None:
            legacy_section = LEGACY_KEY_SECTIONS.get(legacy_key)
        if not aliases:
            aliases = tuple(
                alias for alias, canonical in LEGACY_KEY_ALIASES.items() if canonical == legacy_key
            )
        if materialization_group is None:
            materialization_group = MATERIALIZATION_GROUPS.get(legacy_key)
        if reverse_policy is None:
            reverse_policy = (
                "alias" if legacy_key in LEGACY_KEY_ALIASES else "typed_or_sidecar"
            )
    if required_bundle_paths is None:
        required_bundle_paths = ()
    if reverse_policy is None:
        reverse_policy = "not_reversible"
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
        parser_id=parser_id,
        exporter_id=exporter_id,
        required_bundle_paths=required_bundle_paths,
        legacy_filename_stem=legacy_filename_stem,
        legacy_section=legacy_section,
        reverse_policy=reverse_policy,
        aliases=aliases,
        materialization_group=materialization_group,
    )


REVERSE_POLICIES = frozenset(
    {
        "alias",
        "embedded_text",
        "not_reversible",
        "scalar",
        "sidecar_only",
        "typed_or_sidecar",
        "typed_required",
    }
)


LEGACY_KEY_ALIASES = {
    "improper_in_file": "improper_dihedral_in_file",
    "virtual_atoms_in_file": "virtual_atom_in_file",
    "hill_in_file": "hills_in_file",
    "metad_hills_in_file": "hills_in_file",
}


LEGACY_KEY_SECTIONS = {
    "EAM_in_file": "EAM",
    "EAM_atom_type_in_file": "EAM",
    "TERSOFF_in_file": "TERSOFF",
    "REAXFF_in_file": "REAXFF",
    "REAXFF_type_in_file": "REAXFF",
}


LEGACY_KEY_EXPORTER_IDS = {
    "restrain_coordinate_in_file": "restraint_reference",
    "restrain_amber_rst7": "restraint_reference",
}


MATERIALIZATION_GROUPS = {
    "coordinate_in_file": "restart.coordinate_file",
    "constrain_in_file": "protocol.constraint_file",
    "LJ_in_file": "topology.lj_file",
    "LJ_soft_core_in_file": "topology.lj_soft_core_file",
    "cmap_in_file": "topology.cmap_file",
    "restrain_coordinate_in_file": "restart.restraint_reference_file",
    "restrain_amber_rst7": "restart.restraint_reference_file",
}


SPONGE_SERIALIZER_TO_LEGACY_KEY = {
    "coordinate": "coordinate_in_file",
    "mass": "mass_in_file",
    "charge": "charge_in_file",
    "residue": "residue_in_file",
    "exclude": "exclude_in_file",
    "bond": "bond_in_file",
    "bond_soft": "bond_soft_in_file",
    "angle": "angle_in_file",
    "dihedral": "dihedral_in_file",
    "improper_dihedral": "improper_dihedral_in_file",
    "LJ": "LJ_in_file",
    "nb14": "nb14_in_file",
    "nb14_extra": "nb14_extra_in_file",
    "urey_bradley": "urey_bradley_in_file",
    "cmap": "cmap_in_file",
    "gb": "gb_in_file",
    "virtual_atom": "virtual_atom_in_file",
    "LJ_soft_core": "LJ_soft_core_in_file",
    "subsys_division": "subsys_division_in_file",
    "SW": "SW_in_file",
    "EDIP": "EDIP_in_file",
    "listed_forces": "listed_forces_in_file",
}


SPONGE_SERIALIZER_METADATA_KEYS = frozenset({"resname", "atom_name", "atom_type_name"})


SPONGE_SERIALIZER_LISTED_FORCE_KEYS = frozenset({"Ryckaert_Bellemans"})


SPONGE_SERIALIZER_UNSUPPORTED_KEYS = frozenset(
    {"fake_mass", "fake_LJ", "fake_charge"}
)


def _legacy_filename_stem(key: str) -> str:
    """Return the canonical legacy filename stem for a normalized mdin key."""

    return key.removesuffix("_in_file")


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
    "nose_hoover_chain_length",
    "monte_carlo_barostat_seed",
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
    "output_h5_restart_topology_hash",
    "output_h5_restart_atom_order_hash",
    "output_h5_restart_protocol_hash",
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
    "bond_soft_in_file": "/forcefield/bond_soft",
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
        parser_id="coordinate_in_file",
        exporter_id="coordinate_in_file",
        required_bundle_paths=(
            "/particles/all/position/value",
            "/particles/all/box/edges/value",
        ),
        legacy_filename_stem="coordinate",
        reverse_policy="typed_required",
        materialization_group="restart.coordinate_file",
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
        parser_id="velocity_in_file",
        exporter_id="velocity_in_file",
        required_bundle_paths=("/particles/all/velocity/value",),
        legacy_filename_stem="velocity",
        reverse_policy="typed_or_sidecar",
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
        parser_id="coordinate_in_file",
        exporter_id="coordinate_in_file",
        required_bundle_paths=(
            "/particles/all/position/value",
            "/particles/all/box/edges/value",
        ),
        legacy_filename_stem="coordinate",
        reverse_policy="typed_required",
        materialization_group="restart.coordinate_file",
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
        reverse_policy="not_reversible",
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
        legacy_filename_stem=_legacy_filename_stem(key),
        reverse_policy="sidecar_only",
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
        reverse_policy="scalar",
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
        reverse_policy=(
            "typed_required"
            if key in {"improper_dihedral_in_file", "residue_in_file"}
            else None
        ),
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


def contracts_by_bundle_file() -> dict[str, tuple[IOContract, ...]]:
    """Return all known contracts indexed by bundled artifact name."""

    index: dict[str, list[IOContract]] = {}
    for contract in CONTRACTS:
        index.setdefault(contract.bundle_file, []).append(contract)
    return {key: tuple(value) for key, value in index.items()}


def reversible_contracts(mode: str) -> tuple[IOContract, ...]:
    """Return canonical input contracts that participate in reverse conversion."""

    mode_class = "rerun" if mode.strip().lower() == "rerun" else "normal"
    return tuple(
        contract
        for contract in CONTRACTS
        if contract.direction == "input"
        and mode_class in contract.modes
        and contract.reverse_policy not in {"alias", "not_reversible"}
    )


def classify_sponge_serializer_key(key: str) -> tuple[str, str | None]:
    """Classify a registered ``save_sponge_input`` serializer key."""

    if key in SPONGE_SERIALIZER_TO_LEGACY_KEY:
        return "contract", SPONGE_SERIALIZER_TO_LEGACY_KEY[key]
    if key in SPONGE_SERIALIZER_METADATA_KEYS:
        return "metadata", None
    if key in SPONGE_SERIALIZER_LISTED_FORCE_KEYS:
        return "listed_force", None
    if key in SPONGE_SERIALIZER_UNSUPPORTED_KEYS:
        return "unsupported", None
    return "unknown", None


def validate_contract_registry(contracts: tuple[IOContract, ...] = CONTRACTS) -> None:
    """Validate static invariants required by both conversion directions."""

    contract_ids: set[str] = set()
    legacy_keys = {key for contract in contracts for key in contract.legacy_keys}
    canonical_exporters: dict[tuple[str, str], IOContract] = {}

    for contract in contracts:
        if contract.contract_id in contract_ids:
            raise ValueError(f"duplicate I/O contract id: {contract.contract_id}")
        contract_ids.add(contract.contract_id)

        if not contract.legacy_keys:
            raise ValueError(f"{contract.contract_id} has no legacy keys")
        if not contract.bundle_file or not contract.bundle_path:
            raise ValueError(f"{contract.contract_id} has an empty bundle binding")
        if not contract.modes:
            raise ValueError(f"{contract.contract_id} has no supported modes")
        if contract.reverse_policy not in REVERSE_POLICIES:
            raise ValueError(
                f"{contract.contract_id} has unknown reverse policy {contract.reverse_policy!r}"
            )
        if contract.bundle_file not in {"run.mdin", "*.legacy"} and not contract.bundle_path.startswith("/"):
            raise ValueError(
                f"{contract.contract_id} has non-absolute HDF5 path {contract.bundle_path!r}"
            )
        if contract.reverse_policy in {"typed_or_sidecar", "typed_required"}:
            if contract.exporter_id is None:
                raise ValueError(f"{contract.contract_id} has no reverse exporter id")
            if not contract.required_bundle_paths:
                raise ValueError(f"{contract.contract_id} has no required bundle paths")
            if contract.legacy_filename_stem is None:
                raise ValueError(f"{contract.contract_id} has no legacy filename stem")
        if contract.reverse_policy == "alias":
            for key in contract.legacy_keys:
                canonical = LEGACY_KEY_ALIASES.get(key)
                if canonical is None or canonical not in legacy_keys:
                    raise ValueError(
                        f"{contract.contract_id} alias {key!r} has no canonical contract"
                    )
        for alias in contract.aliases:
            if LEGACY_KEY_ALIASES.get(alias) not in contract.legacy_keys:
                raise ValueError(
                    f"{contract.contract_id} does not own declared alias {alias!r}"
                )

        if contract.exporter_id is not None and contract.bundle_file not in {"run.mdin", "*.legacy"}:
            route = (contract.bundle_file, contract.bundle_path)
            previous = canonical_exporters.get(route)
            if previous is not None and (
                previous.materialization_group != contract.materialization_group
                or previous.exporter_id != contract.exporter_id
            ):
                raise ValueError(
                    f"conflicting exporters for {route}: "
                    f"{previous.contract_id} and {contract.contract_id}"
                )
            canonical_exporters[route] = contract
