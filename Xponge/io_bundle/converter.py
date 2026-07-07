"""
Legacy SPONGE input group to bundled H5 input group converter.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil

import numpy as np

from Xponge.load import load_coordinate

from .contracts import CONTRACTS, PROTOCOL_RESTART_SIDECAR_KEYS, IOContract, contracts_by_legacy_key
from .h5_writer import (
    ensure_dataset,
    ensure_group,
    ensure_hard_link,
    set_attrs,
    write_bytes,
    write_string,
    write_string_array,
)
from .legacy_case import LegacyCase, render_mdin_without_keys, scan_legacy_case
from .manifest import ConversionManifest, ManifestEntry
from .state_parsers import parse_protocol_or_restart_file
from .topology_parsers import (
    listed_force_schemas,
    pairwise_force_schema,
    parse_listed_force_data_file,
    parse_pairwise_force_data_file,
    parse_topology_file,
)
from .trajectory_parsers import box_to_edges, parse_trajectory_file


class ConversionError(RuntimeError):
    """Raised when a legacy case cannot be converted."""


class LegacyToBundleConverter:
    """Convert a legacy SPONGE case directory into bundled input files."""

    def __init__(self, case: LegacyCase, output_dir: str | Path):
        self.case = case
        self.output_dir = Path(output_dir).resolve()
        self.bundle_dir = self.output_dir / "bundle"
        self.manifest = ConversionManifest(case_root=str(self.case.root), mode=self.case.mode)
        self._contract_index = contracts_by_legacy_key()
        self._handled_contract_ids: set[str] = set()
        self._handled_keys: set[str] = set()
        self._bundled_mdin_omit_keys: set[str] = set()
        self._bundled_mdin_extra_lines: list[str] = []
        self._bundle_files_touched: set[str] = set()
        self._legacy_sidecars_by_bundle: dict[str, list[tuple[str, str]]] = {}
        self._source_paths_by_key: dict[str, Path] = {}
        self._atom_count: int | None = None
        self._restart_has_structural_state = False
        self._restart_has_dynamic_state = False
        self._restart_has_protocol_state = False

    def convert(self, *, dry_run: bool = False) -> ConversionManifest:
        """Run the conversion and return the generated manifest."""

        self.bundle_dir.mkdir(parents=True, exist_ok=True)
        self._convert_restart_structural_state(dry_run=dry_run)
        self._convert_text_sidecars(dry_run=dry_run)
        self._import_compatibility_payloads(dry_run=dry_run)
        self._convert_protocol_restart_sidecars(dry_run=dry_run)
        self._convert_custom_force_dynamic_payloads(dry_run=dry_run)
        self._ensure_required_h5_input_bindings()
        self._write_legacy_sidecar_tables(dry_run=dry_run)
        self._write_bundle_metadata(dry_run=dry_run)
        self._record_scalar_contracts()
        self._record_default_legacy_outputs()
        self._record_unsupported_present_contracts()
        self._write_bundled_mdin(dry_run=dry_run)
        if not dry_run:
            self.manifest.write(self.output_dir / "manifest.json")
        return self.manifest

    def _find_contract(self, contract_id: str) -> IOContract:
        for contract in CONTRACTS:
            if contract.contract_id == contract_id:
                return contract
        raise ConversionError(f"unknown converter contract: {contract_id}")

    def _contract_applies(self, contract: IOContract) -> bool:
        mode_class = "rerun" if self.case.mode == "rerun" else "normal"
        return mode_class in contract.modes

    def _add_converted(self, contract: IOContract, source_key: str, source_path: Path) -> None:
        self._handled_contract_ids.add(contract.contract_id)
        self._handled_keys.add(source_key)
        self._bundle_files_touched.add(contract.bundle_file)
        if contract.direction == "input" and contract.payload_kind == "file":
            self._bundled_mdin_omit_keys.add(source_key)
        self.manifest.add(
            ManifestEntry(
                contract_id=contract.contract_id,
                status="converted",
                source_key=source_key,
                source_path=str(source_path),
                bundle_file=contract.bundle_file,
                bundle_path=contract.bundle_path,
                direction=contract.direction,
                component=contract.component,
                payload_kind=contract.payload_kind,
                override_policy=contract.override_policy,
                comparison_rule=contract.comparison_rule,
            )
        )

    def _add_manifest_entry(
        self,
        contract: IOContract,
        status: str,
        *,
        source_key: str,
        source_path: Path | None = None,
        message: str | None = None,
    ) -> None:
        self._handled_contract_ids.add(contract.contract_id)
        self._handled_keys.add(source_key)
        if contract.bundle_file not in ("run.mdin", "*.legacy"):
            self._bundle_files_touched.add(contract.bundle_file)
        if contract.direction == "input" and contract.payload_kind == "file":
            self._bundled_mdin_omit_keys.add(source_key)
        self.manifest.add(
            ManifestEntry(
                contract_id=contract.contract_id,
                status=status,
                source_key=source_key,
                source_path=str(source_path) if source_path is not None else None,
                bundle_file=contract.bundle_file,
                bundle_path=contract.bundle_path,
                direction=contract.direction,
                component=contract.component,
                payload_kind=contract.payload_kind,
                override_policy=contract.override_policy,
                comparison_rule=contract.comparison_rule,
                message=message,
            )
        )

    def _convert_restart_structural_state(self, *, dry_run: bool) -> None:
        coordinate_path = self.case.resolve_legacy_input_path("coordinate_in_file")
        if coordinate_path is None:
            return
        if not coordinate_path.exists():
            raise ConversionError(f"coordinate_in_file does not exist: {coordinate_path}")

        coordinate, box = load_coordinate(str(coordinate_path))
        self._atom_count = int(coordinate.shape[0])
        box_edges = _box_to_edges(box)
        restart_path = self.bundle_dir / "restart.spgr.h5"

        if not dry_run:
            ensure_dataset(restart_path, "/particles/all/position/value", coordinate[np.newaxis, :, :])
            ensure_dataset(restart_path, "/particles/all/box/edges/value", box_edges[np.newaxis, :, :])
            ensure_dataset(restart_path, "/particles/all/step", np.asarray([0], dtype=np.int64))
            ensure_dataset(restart_path, "/particles/all/time", np.asarray([0.0], dtype=np.float64))
        self._restart_has_structural_state = True

        self._add_converted(self._find_contract("restart.coordinate"), "coordinate_in_file", coordinate_path)
        self._add_converted(self._find_contract("restart.box"), "coordinate_in_file", coordinate_path)

        velocity_path = self.case.resolve_legacy_input_path("velocity_in_file")
        if velocity_path is not None:
            if not velocity_path.exists():
                raise ConversionError(f"velocity_in_file does not exist: {velocity_path}")
            velocity = _read_vector_file(velocity_path, expected_count=coordinate.shape[0])
            if not dry_run:
                ensure_dataset(restart_path, "/particles/all/velocity/value", velocity[np.newaxis, :, :])
            self._add_converted(self._find_contract("restart.velocity"), "velocity_in_file", velocity_path)

    def _convert_text_sidecars(self, *, dry_run: bool) -> None:
        for key, contracts in self._contract_index.items():
            source_path = self.case.resolve_legacy_input_path(key)
            if source_path is None:
                continue
            for contract in contracts:
                if not self._contract_applies(contract):
                    continue
                if contract.status != "supported" or not contract.bundle_path.startswith("/parameters"):
                    continue
                if _has_typed_state_parser(key, source_path):
                    continue
                if contract.contract_id in self._handled_contract_ids:
                    continue
                if not source_path.exists():
                    raise ConversionError(f"{key} does not exist: {source_path}")
                text = source_path.read_text(encoding="utf-8")
                bundle_path = self.bundle_dir / contract.bundle_file
                if not dry_run:
                    write_string(bundle_path, contract.bundle_path, text)
                self._stage_legacy_sidecar(contract, key, source_path, dry_run=dry_run)
                self._add_converted(contract, key, source_path)

    def _import_compatibility_payloads(self, *, dry_run: bool) -> None:
        for key, contracts in self._contract_index.items():
            source_path = self.case.resolve_legacy_input_path(key)
            if source_path is None:
                continue
            for contract in contracts:
                if not self._contract_applies(contract):
                    continue
                if contract.contract_id in self._handled_contract_ids:
                    continue
                if contract.status != "compatibility_import":
                    continue
                if not source_path.exists():
                    raise ConversionError(f"{key} does not exist: {source_path}")
                typed_datasets = _parse_typed_payload(
                    contract,
                    key,
                    source_path,
                    atom_count=self._infer_atom_count() if contract.component == "trajectory" else None,
                )
                if typed_datasets is not None:
                    bundle_path = self.bundle_dir / contract.bundle_file
                    if not dry_run:
                        for dataset in typed_datasets:
                            ensure_dataset(bundle_path, dataset.path, dataset.data)
                    if contract.bundle_file == "restart.spgr.h5":
                        self._track_restart_state_component(key)
                    self._stage_legacy_sidecar(contract, key, source_path, dry_run=dry_run)
                    self._add_manifest_entry(
                        contract,
                        "typed_converted",
                        source_key=key,
                        source_path=source_path,
                        message=f"legacy {contract.component} payload materialized as typed HDF5 datasets",
                    )
                    continue
                import_path = _compatibility_payload_path(key)
                bundle_path = self.bundle_dir / contract.bundle_file
                if not dry_run:
                    if _looks_text_file(source_path):
                        write_string(bundle_path, import_path + "/raw_text", source_path.read_text(encoding="utf-8"))
                    else:
                        write_bytes(bundle_path, import_path + "/raw_bytes", source_path.read_bytes())
                    write_string(bundle_path, import_path + "/source_path", str(source_path))
                    write_string(bundle_path, import_path + "/canonical_path", contract.bundle_path)
                self._stage_legacy_sidecar(contract, key, source_path, dry_run=dry_run)
                self._add_manifest_entry(
                    contract,
                    "compatibility_imported",
                    source_key=key,
                    source_path=source_path,
                    message=(
                        "legacy payload preserved under /compatibility/legacy_import; "
                        "native typed materialization for this contract is still required"
                    ),
                    )

    def _convert_protocol_restart_sidecars(self, *, dry_run: bool) -> None:
        for key in _PROTOCOL_RESTART_SIDECAR_KEYS:
            source_path = self.case.resolve_legacy_input_path(key)
            if source_path is None:
                continue
            if not source_path.exists():
                raise ConversionError(f"{key} does not exist: {source_path}")
            restart_path = self.bundle_dir / "restart.spgr.h5"
            if not dry_run:
                write_string(
                    restart_path,
                    f"/parameters/restart/protocol_sidecars/{key}",
                    source_path.read_text(encoding="utf-8"),
                )
            self._stage_legacy_sidecar_for_bundle("restart.spgr.h5", key, source_path, dry_run=dry_run)
            self._handled_contract_ids.add(f"restart.protocol_sidecar.{key}")
            self._handled_keys.add(key)
            self._bundled_mdin_omit_keys.add(key)
            self._bundle_files_touched.add("restart.spgr.h5")
            self._restart_has_protocol_state = True
            self.manifest.add(
                ManifestEntry(
                    contract_id=f"restart.protocol_sidecar.{key}",
                    status="sidecar_embedded",
                    source_key=key,
                    source_path=str(source_path),
                    bundle_file="restart.spgr.h5",
                    bundle_path=f"/parameters/restart/protocol_sidecars/{key}",
                    direction="input",
                    component="restart",
                    payload_kind="file",
                    override_policy="legacy_protocol_sidecar",
                    comparison_rule="text",
                    message="legacy protocol sidecar text embedded for SPONGE H5 restart materialization",
                )
            )

    def _track_restart_state_component(self, key: str) -> None:
        if key == "nose_hoover_chain_restart_input":
            self._restart_has_dynamic_state = True
        elif key in {
            "SITS_nk_in_file",
            "meta_potential_in_file",
            "meta_scatter_in_file",
            "hill_in_file",
            "hills_in_file",
            "metad_hills_in_file",
        }:
            self._restart_has_protocol_state = True

    def _convert_custom_force_dynamic_payloads(self, *, dry_run: bool) -> None:
        self._convert_pairwise_force_dynamic_payload(dry_run=dry_run)
        self._convert_listed_force_dynamic_payloads(dry_run=dry_run)

    def _convert_pairwise_force_dynamic_payload(self, *, dry_run: bool) -> None:
        config_path = self.case.resolve_legacy_input_path("pairwise_force_in_file")
        if config_path is None:
            return
        if not config_path.exists():
            raise ConversionError(f"pairwise_force_in_file does not exist: {config_path}")
        force_name, parameter_types, parameter_names, ij_count = pairwise_force_schema(config_path)
        source_key, source_path = self._resolve_dynamic_in_file(force_name)
        typed_datasets = parse_pairwise_force_data_file(
            force_name,
            source_path,
            parameter_types=parameter_types,
            parameter_names=parameter_names,
            ij_parameter_count=ij_count,
        )
        self._write_dynamic_custom_force_payload(
            contract_id=f"topology.pairwise_force_data.{force_name}",
            source_key=source_key,
            source_path=source_path,
            bundle_root=f"/forcefield/custom_force/pairwise/data/{_safe_h5_name(force_name)}",
            typed_datasets=typed_datasets,
            dry_run=dry_run,
        )

    def _convert_listed_force_dynamic_payloads(self, *, dry_run: bool) -> None:
        config_path = self.case.resolve_legacy_input_path("listed_forces_in_file")
        if config_path is None:
            return
        if not config_path.exists():
            raise ConversionError(f"listed_forces_in_file does not exist: {config_path}")
        for force_name, parameter_types, parameter_names in listed_force_schemas(config_path):
            source_key, source_path = self._resolve_dynamic_in_file(force_name)
            typed_datasets = parse_listed_force_data_file(
                force_name,
                source_path,
                parameter_types=parameter_types,
                parameter_names=parameter_names,
            )
            self._write_dynamic_custom_force_payload(
                contract_id=f"topology.listed_force_data.{force_name}",
                source_key=source_key,
                source_path=source_path,
                bundle_root=f"/forcefield/custom_force/listed/data/{_safe_h5_name(force_name)}",
                typed_datasets=typed_datasets,
                dry_run=dry_run,
            )

    def _resolve_dynamic_in_file(self, force_name: str) -> tuple[str, Path]:
        source_key = f"{force_name}_in_file"
        source_value = self.case.commands.get(source_key)
        if source_value:
            source_path = self.case.resolve_value_path(source_value)
            if not source_path.exists():
                raise ConversionError(f"{source_key} does not exist: {source_path}")
            return source_key, source_path

        default_prefix = self.case.commands.get("default_in_file_prefix")
        if default_prefix:
            candidate = self.case.resolve_value_path(f"{default_prefix}_{force_name}.txt")
            if candidate.exists():
                return source_key, candidate

        raise ConversionError(
            f"{source_key} is required by custom force {force_name!r} but was not found"
        )

    def _write_dynamic_custom_force_payload(
        self,
        *,
        contract_id: str,
        source_key: str,
        source_path: Path,
        bundle_root: str,
        typed_datasets,
        dry_run: bool,
    ) -> None:
        bundle_file = "topology.spgt.h5"
        if not dry_run:
            bundle_path = self.bundle_dir / bundle_file
            for dataset in typed_datasets:
                ensure_dataset(bundle_path, dataset.path, dataset.data)
        sidecar_rel = self._stage_dynamic_mdin_sidecar(source_key, source_path, dry_run=dry_run)
        self._source_paths_by_key[source_key] = source_path
        self._handled_keys.add(source_key)
        self._bundled_mdin_omit_keys.add(source_key)
        self._bundle_files_touched.add(bundle_file)
        self._bundled_mdin_extra_lines.append(f'{source_key} = "{sidecar_rel.as_posix()}"')
        self.manifest.add(
            ManifestEntry(
                contract_id=contract_id,
                status="typed_converted",
                source_key=source_key,
                source_path=str(source_path),
                bundle_file=bundle_file,
                bundle_path=bundle_root,
                direction="input",
                component="topology",
                payload_kind="file",
                override_policy="legacy_mdin_sidecar",
                comparison_rule="typed",
                message=(
                    "dynamic custom force payload materialized as typed HDF5 datasets "
                    "and rebound to a bundle-local legacy mdin sidecar"
                ),
            )
        )

    def _stage_dynamic_mdin_sidecar(self, key: str, source_path: Path, *, dry_run: bool) -> Path:
        sidecar_rel = Path("legacy_sidecars") / key / source_path.name
        sidecar_abs = self.bundle_dir / sidecar_rel
        if not dry_run:
            sidecar_abs.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_path, sidecar_abs)
        return sidecar_rel

    def _stage_legacy_sidecar(
        self, contract: IOContract, key: str, source_path: Path, *, dry_run: bool
    ) -> None:
        if contract.direction != "input" or contract.bundle_file == "trajectory.spg.h5md":
            return
        self._stage_legacy_sidecar_for_bundle(contract.bundle_file, key, source_path, dry_run=dry_run)

    def _stage_legacy_sidecar_for_bundle(
        self, bundle_file: str, key: str, source_path: Path, *, dry_run: bool
    ) -> None:
        sidecar_rel = Path("legacy_sidecars") / key / source_path.name
        sidecar_abs = self.bundle_dir / sidecar_rel
        if not dry_run:
            sidecar_abs.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_path, sidecar_abs)
        self._source_paths_by_key[key] = source_path
        sidecars = self._legacy_sidecars_by_bundle.setdefault(bundle_file, [])
        entry = (key, sidecar_rel.as_posix())
        if entry not in sidecars:
            sidecars.append(entry)

    def _ensure_required_h5_input_bindings(self) -> None:
        h5_input_files = {
            "topology.spgt.h5",
            "protocol.spgp.h5",
            "restart.spgr.h5",
            "trajectory.spg.h5md",
        }
        touched = self._bundle_files_touched & h5_input_files
        if not touched:
            return
        if "topology.spgt.h5" not in touched:
            raise ConversionError(
                "H5 input bindings require topology.spgt.h5; provide explicit or default-prefixed topology files"
            )
        if "protocol.spgp.h5" not in touched:
            self._bundle_files_touched.add("protocol.spgp.h5")
            self.manifest.add(
                ManifestEntry(
                    contract_id="protocol.required_empty",
                    status="generated",
                    bundle_file="protocol.spgp.h5",
                    bundle_path="/protocol/topology_compatibility",
                    direction="input",
                    component="protocol",
                    payload_kind="metadata",
                    override_policy="required_binding",
                    comparison_rule="metadata",
                    message="generated minimal protocol bundle because SPONGE requires protocol binding for H5 input",
                )
            )
        has_restart = "restart.spgr.h5" in touched
        has_rerun_trajectory = self.case.mode == "rerun" and "trajectory.spg.h5md" in touched
        if not has_restart and not has_rerun_trajectory:
            raise ConversionError(
                "H5 input bindings require restart.spgr.h5 unless rerun trajectory input is bundled"
            )

    def _write_legacy_sidecar_tables(self, *, dry_run: bool) -> None:
        for bundle_file, sidecars in self._legacy_sidecars_by_bundle.items():
            if not sidecars:
                continue
            keys = [key for key, _ in sidecars]
            paths = [path for _, path in sidecars]
            bundle_path = self.bundle_dir / bundle_file
            if not dry_run:
                write_string_array(bundle_path, "/parameters/sponge/files/legacy_sidecars/key", keys)
                write_string_array(bundle_path, "/parameters/sponge/files/legacy_sidecars/path", paths)

    def _write_bundle_metadata(self, *, dry_run: bool) -> None:
        if dry_run:
            return
        schema_version = "xponge.legacy_to_bundle.v1"
        topology_hash = self._content_hash("topology.spgt.h5")
        atom_order_hash = self._content_hash("topology.spgt.h5", include_keys=_TOPOLOGY_ATOM_ORDER_KEYS)
        forcefield_hash = self._content_hash("topology.spgt.h5", exclude_keys=_TOPOLOGY_ATOM_ORDER_KEYS)
        protocol_hash = self._content_hash("protocol.spgp.h5")

        if "topology.spgt.h5" in self._bundle_files_touched:
            topology_path = self.bundle_dir / "topology.spgt.h5"
            write_string(topology_path, "/schema/name", "sponge.topology.h5")
            write_string(topology_path, "/schema/version", schema_version)
            write_string(topology_path, "/parameters/sponge/schema/name", "sponge.topology.h5")
            write_string(topology_path, "/parameters/sponge/schema/version", schema_version)
            write_string(topology_path, "/topology/atom_order_hash", atom_order_hash)
            write_string(topology_path, "/topology/topology_hash", topology_hash)
            write_string(topology_path, "/topology/forcefield_hash", forcefield_hash)
            ensure_dataset(topology_path, "/topology/atom_count", np.asarray(self._infer_atom_count(), dtype=np.int64))

        if "protocol.spgp.h5" in self._bundle_files_touched:
            protocol_path = self.bundle_dir / "protocol.spgp.h5"
            write_string(protocol_path, "/schema/name", "sponge.protocol.h5")
            write_string(protocol_path, "/schema/version", schema_version)
            write_string(protocol_path, "/parameters/sponge/schema/name", "sponge.protocol.h5")
            write_string(protocol_path, "/parameters/sponge/schema/version", schema_version)
            write_string(protocol_path, "/protocol/topology_compatibility/topology_hash", topology_hash)
            write_string(protocol_path, "/identity/content_hash", protocol_hash)
            cv_count, restraint_count, enhanced_methods = self._protocol_metadata_summary()
            ensure_dataset(protocol_path, "/protocol/cv_count", np.asarray(cv_count, dtype=np.int64))
            ensure_dataset(protocol_path, "/protocol/restraint_count", np.asarray(restraint_count, dtype=np.int64))
            if enhanced_methods:
                write_string(protocol_path, "/protocol/enhanced_sampling/method", ",".join(enhanced_methods))

        if "restart.spgr.h5" in self._bundle_files_touched:
            restart_path = self.bundle_dir / "restart.spgr.h5"
            write_string(restart_path, "/parameters/sponge/schema/name", "sponge.restart.h5")
            write_string(restart_path, "/parameters/sponge/schema/version", schema_version)
            self._finalize_restart_h5_layout(restart_path)

        if "trajectory.spg.h5md" in self._bundle_files_touched:
            trajectory_path = self.bundle_dir / "trajectory.spg.h5md"
            write_string(trajectory_path, "/parameters/sponge/schema/name", "sponge.output.h5md")
            write_string(trajectory_path, "/parameters/sponge/schema/version", schema_version)
            self._finalize_trajectory_h5md_layout(trajectory_path)

    def _finalize_trajectory_h5md_layout(self, trajectory_path: Path) -> None:
        ensure_group(trajectory_path, "/h5md")
        ensure_group(trajectory_path, "/h5md/creator")
        ensure_group(trajectory_path, "/particles/all")
        set_attrs(trajectory_path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        set_attrs(trajectory_path, "/h5md/creator", {"name": "XPONGE", "version": "legacy-to-bundle"})

        self._link_particle_time_axis(
            trajectory_path,
            value_path="/particles/all/position/value",
            step_path="/particles/all/position/step",
            time_path="/particles/all/position/time",
            value_attrs={"unit": "Angstrom"},
        )
        self._link_particle_time_axis(
            trajectory_path,
            value_path="/particles/all/velocity/value",
            step_path="/particles/all/velocity/step",
            time_path="/particles/all/velocity/time",
            value_attrs={"unit": "Angstrom ps-1"},
        )
        self._link_particle_time_axis(
            trajectory_path,
            value_path="/particles/all/box/edges/value",
            step_path="/particles/all/box/edges/step",
            time_path="/particles/all/box/edges/time",
            value_attrs={"unit": "Angstrom"},
        )
        set_attrs(trajectory_path, "/particles/all/time", {"unit": "ps"})
        frame_count, last_step, last_time = _particle_completion_metadata(trajectory_path)
        write_string(trajectory_path, "/parameters/sponge/output/mode", "single")
        write_string(trajectory_path, "/parameters/sponge/output/status", "finalized")
        ensure_dataset(
            trajectory_path,
            "/parameters/sponge/output/frame_count",
            np.asarray([frame_count], dtype=np.int64),
        )
        ensure_dataset(
            trajectory_path,
            "/parameters/sponge/output/last_complete_step",
            np.asarray([last_step], dtype=np.int64),
        )
        ensure_dataset(
            trajectory_path,
            "/parameters/sponge/output/last_complete_time",
            np.asarray([last_time], dtype=np.float64),
        )
        write_string_array(trajectory_path, "/parameters/sponge/output/particle_streams", ["all"])

    def _protocol_metadata_summary(self) -> tuple[int, int, list[str]]:
        cv_count = 1 if self.case.resolve_legacy_input_path("cv_in_file") is not None else 0
        restraint_count = 0
        for key in (
            "restrain_in_file",
            "restrain_atom_id",
            "restrain_coordinate_in_file",
            "restrain_weight_in_file",
            "restrain_cv_in_file",
        ):
            if self.case.resolve_legacy_input_path(key) is not None:
                restraint_count += 1
        enhanced_methods = []
        if any(
            self.case.resolve_legacy_input_path(key) is not None
            for key in ("SITS_in_file", "SITS_atom_in_file", "SITS_nk_in_file")
        ):
            enhanced_methods.append("SITS")
        if any(
            self.case.resolve_legacy_input_path(key) is not None
            for key in (
                "meta_edge_in_file",
                "meta_potential_in_file",
                "meta_scatter_in_file",
                "hill_in_file",
                "hills_in_file",
                "metad_hills_in_file",
            )
        ):
            enhanced_methods.append("metadynamics")
        if self.case.resolve_legacy_input_path("steer_cv_in_file") is not None:
            enhanced_methods.append("steer")
        if self.case.resolve_legacy_input_path("soft_walls_in_file") is not None:
            enhanced_methods.append("soft_walls")
        return cv_count, restraint_count, enhanced_methods

    def _finalize_restart_h5_layout(self, restart_path: Path) -> None:
        ensure_group(restart_path, "/h5md")
        ensure_group(restart_path, "/h5md/creator")
        ensure_group(restart_path, "/run")
        ensure_group(restart_path, "/particles/all")
        ensure_group(restart_path, "/parameters/restart")
        ensure_group(restart_path, "/parameters/restart/rng_state")
        ensure_group(restart_path, "/parameters/restart/integrator_state")
        ensure_group(restart_path, "/parameters/restart/thermostat")
        ensure_group(restart_path, "/parameters/restart/barostat")
        ensure_group(restart_path, "/parameters/restart/protocol_sidecars")
        ensure_group(restart_path, "/parameters/restart/bias")
        ensure_group(restart_path, "/parameters/restart/bias/sits")
        ensure_group(restart_path, "/parameters/restart/bias/meta")
        set_attrs(restart_path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        set_attrs(restart_path, "/h5md/creator", {"name": "XPONGE", "version": "legacy-to-bundle"})

        self._link_particle_time_axis(
            restart_path,
            value_path="/particles/all/position/value",
            step_path="/particles/all/position/step",
            time_path="/particles/all/position/time",
            value_attrs={"unit": "Angstrom"},
        )
        self._link_particle_time_axis(
            restart_path,
            value_path="/particles/all/velocity/value",
            step_path="/particles/all/velocity/step",
            time_path="/particles/all/velocity/time",
            value_attrs={"unit": "Angstrom ps-1"},
        )
        self._link_particle_time_axis(
            restart_path,
            value_path="/particles/all/box/edges/value",
            step_path="/particles/all/box/edges/step",
            time_path="/particles/all/box/edges/time",
            value_attrs={"unit": "Angstrom"},
        )
        try:
            set_attrs(restart_path, "/particles/all/time", {"unit": "ps"})
        except KeyError:
            pass

        frame_count, last_step, last_time = _particle_completion_metadata(restart_path)
        write_string(restart_path, "/parameters/sponge/output/status", "finalized")
        ensure_dataset(
            restart_path,
            "/parameters/sponge/output/frame_count",
            np.asarray([frame_count], dtype=np.int64),
        )
        ensure_dataset(
            restart_path,
            "/parameters/sponge/output/last_complete_step",
            np.asarray([last_step], dtype=np.int64),
        )
        ensure_dataset(
            restart_path,
            "/parameters/sponge/output/last_complete_time",
            np.asarray([last_time], dtype=np.float64),
        )
        if frame_count:
            ensure_dataset(restart_path, "/run/current_step", np.asarray([last_step], dtype=np.int64))
            ensure_dataset(restart_path, "/run/current_time", np.asarray([last_time], dtype=np.float64))
            write_string(restart_path, "/run/state_type", "restart")
            write_string_array(restart_path, "/parameters/sponge/output/particle_streams", ["all"])

    def _link_particle_time_axis(
        self,
        trajectory_path: Path,
        *,
        value_path: str,
        step_path: str,
        time_path: str,
        value_attrs: dict[str, str],
    ) -> None:
        try:
            set_attrs(trajectory_path, value_path, value_attrs)
        except KeyError:
            return
        ensure_hard_link(trajectory_path, "/particles/all/step", step_path)
        ensure_hard_link(trajectory_path, "/particles/all/time", time_path)

    def _content_hash(
        self,
        bundle_file: str,
        *,
        include_keys: set[str] | None = None,
        exclude_keys: set[str] | None = None,
    ) -> str:
        digest = hashlib.sha256()
        digest.update(bundle_file.encode("utf-8"))
        sidecars = sorted(self._legacy_sidecars_by_bundle.get(bundle_file, ()))
        for key, rel_path in sidecars:
            if include_keys is not None and key not in include_keys:
                continue
            if exclude_keys is not None and key in exclude_keys:
                continue
            digest.update(b"\0")
            digest.update(key.encode("utf-8"))
            digest.update(b"\0")
            digest.update(rel_path.encode("utf-8"))
            source_path = self._source_paths_by_key.get(key)
            if source_path is not None and source_path.exists():
                digest.update(b"\0")
                digest.update(source_path.read_bytes())
        return "sha256:" + digest.hexdigest()

    def _record_scalar_contracts(self) -> None:
        for key, contracts in self._contract_index.items():
            if key not in self.case.commands:
                continue
            for contract in contracts:
                if not self._contract_applies(contract):
                    continue
                if contract.contract_id in self._handled_contract_ids:
                    continue
                if contract.status == "run_mdin":
                    if key == "default_in_file_prefix" and self._has_h5_input_bundle():
                        self._bundled_mdin_omit_keys.add(key)
                        self._add_manifest_entry(
                            contract,
                            "legacy_input_override_replaced",
                            source_key=key,
                            message=(
                                "legacy default input prefix was removed because bundled H5 sidecars "
                                "provide explicit input paths"
                            ),
                        )
                    elif key in _input_binding_keys_for_bundle_files(self._bundle_files_touched):
                        self._bundled_mdin_omit_keys.add(key)
                        self._add_manifest_entry(
                            contract,
                            "legacy_input_override_replaced",
                            source_key=key,
                            message="legacy input file was converted and replaces this existing H5 input binding",
                        )
                    else:
                        self._add_manifest_entry(contract, "preserved_in_mdin", source_key=key)
                elif contract.status == "output_plan":
                    self._add_manifest_entry(contract, "output_plan_preserved", source_key=key)
                elif contract.status == "legacy_output_sidecar":
                    self._add_manifest_entry(contract, "legacy_output_sidecar_preserved", source_key=key)

    def _record_default_legacy_outputs(self) -> None:
        prefix = self.case.commands.get("default_out_file_prefix")
        if not prefix or self._has_h5_output_plan():
            return
        for key in _DEFAULT_LEGACY_OUTPUT_KEYS:
            if key in self.case.commands:
                continue
            contract = self._find_contract(f"output.legacy_sidecar.{key}")
            self.manifest.add(
                ManifestEntry(
                    contract_id=contract.contract_id,
                    status="legacy_output_sidecar_default",
                    source_key="default_out_file_prefix",
                    source_path=_default_legacy_output_path(prefix, key),
                    bundle_file=contract.bundle_file,
                    bundle_path=contract.bundle_path,
                    direction=contract.direction,
                    component=contract.component,
                    payload_kind=contract.payload_kind,
                    override_policy=contract.override_policy,
                    comparison_rule=contract.comparison_rule,
                    message="legacy output sidecar is enabled by default_out_file_prefix because no H5 output path is set",
                )
            )

    def _record_unsupported_present_contracts(self) -> None:
        for contract in CONTRACTS:
            if contract.status == "supported":
                continue
            for key in contract.legacy_keys:
                if (
                    key in self.case.commands
                    and key not in self._handled_keys
                    and contract.contract_id not in self._handled_contract_ids
                    and self._contract_applies(contract)
                ):
                    value = self.case.commands[key]
                    source_path = self.case.resolve_value_path(value) if value else None
                    self._add_manifest_entry(
                        contract,
                        "unsupported",
                        source_key=key,
                        source_path=source_path,
                        message="contract is known but not implemented by this converter yet",
                    )

    def _write_bundled_mdin(self, *, dry_run: bool) -> None:
        output_path = self.bundle_dir / "mdin.bundled.spg.toml"
        omit_keys = self._bundled_mdin_omit_keys | _input_binding_keys_for_bundle_files(self._bundle_files_touched)
        append_lines = []
        if "topology.spgt.h5" in self._bundle_files_touched:
            append_lines.append('input_h5_topology_path = "topology.spgt.h5"')
        if "protocol.spgp.h5" in self._bundle_files_touched:
            append_lines.append('input_h5_protocol_path = "protocol.spgp.h5"')
        if "restart.spgr.h5" in self._bundle_files_touched:
            append_lines.append('input_h5_restart_path = "restart.spgr.h5"')
            append_lines.append(f'input_h5_restart_load = "{self._restart_load_policy()}"')
        if "trajectory.spg.h5md" in self._bundle_files_touched:
            append_lines.append('input_h5_trajectory_path = "trajectory.spg.h5md"')
            append_lines.append('input_h5_trajectory_particle_stream = "all"')
        append_lines.extend(self._bundled_mdin_extra_lines)
        self.manifest.bundled_mdin = str(output_path)
        if not dry_run:
            output_path.write_text(
                render_mdin_without_keys(self.case.mdin_text, omit_keys, append_lines),
                encoding="utf-8",
            )

    def _has_h5_input_bundle(self) -> bool:
        return bool(
            self._bundle_files_touched
            & {"topology.spgt.h5", "protocol.spgp.h5", "restart.spgr.h5", "trajectory.spg.h5md"}
        )

    def _has_h5_output_plan(self) -> bool:
        return any(self.case.commands.get(key) for key in _H5_OUTPUT_PATH_KEYS)

    def _restart_load_policy(self) -> str:
        if self._restart_has_structural_state and (
            self._restart_has_dynamic_state or self._restart_has_protocol_state
        ):
            return "full"
        if self._restart_has_dynamic_state and self._restart_has_protocol_state:
            return "full"
        if self._restart_has_dynamic_state:
            return "dynamic"
        if self._restart_has_protocol_state:
            return "protocol"
        return "structural"

    def _infer_atom_count(self) -> int:
        if self._atom_count is not None:
            return self._atom_count
        for key in ("mass_in_file", "charge_in_file"):
            path = self.case.resolve_legacy_input_path(key)
            if path is None:
                continue
            if path.exists():
                lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
                if lines:
                    self._atom_count = int(lines[0][0])
                    return self._atom_count
        raise ConversionError("cannot infer atom count needed to convert rerun trajectory inputs")


def convert_legacy_to_bundle(
    case_root: str | Path,
    output_dir: str | Path,
    *,
    mdin: str | Path = "mdin.spg.toml",
    dry_run: bool = False,
) -> ConversionManifest:
    """Convert a legacy case directory and return its manifest."""

    case = scan_legacy_case(case_root, mdin=mdin)
    return LegacyToBundleConverter(case, output_dir).convert(dry_run=dry_run)


def _read_vector_file(path: Path, *, expected_count: int) -> np.ndarray:
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if lines and len(lines[0]) == 1:
        declared = int(lines[0][0])
        if declared != expected_count:
            raise ConversionError(f"{path} declares {declared} atoms but coordinate has {expected_count}")
        lines = lines[1:]
    values = np.asarray([[float(item) for item in line[:3]] for line in lines], dtype=np.float32)
    if values.shape != (expected_count, 3):
        raise ConversionError(f"{path} has shape {values.shape}, expected ({expected_count}, 3)")
    return values


def _box_to_edges(box: np.ndarray) -> np.ndarray:
    try:
        return box_to_edges(box)
    except ValueError as exc:
        raise ConversionError(str(exc)) from exc


def _compatibility_payload_path(key: str) -> str:
    safe_key = "".join(char if char.isalnum() or char == "_" else "_" for char in key)
    return f"/compatibility/legacy_import/{safe_key}"


_H5_OUTPUT_PATH_KEYS = {
    "output_h5_trajectory_path",
    "output_h5_restart_path",
    "output_h5_observable_path",
}


_PROTOCOL_RESTART_SIDECAR_KEYS = PROTOCOL_RESTART_SIDECAR_KEYS


_DEFAULT_LEGACY_OUTPUT_KEYS = (
    "mdout",
    "mdinfo",
    "crd",
    "box",
    "rst",
)


_TOPOLOGY_ATOM_ORDER_KEYS = {
    "mass_in_file",
    "charge_in_file",
    "residue_in_file",
}


def _default_legacy_output_path(prefix: str, key: str) -> str:
    suffixes = {
        "mdout": ".out",
        "mdinfo": ".info",
        "crd": ".dat",
        "box": ".box",
        "rst": "",
    }
    return prefix + suffixes[key]


def _safe_h5_name(name: str) -> str:
    safe = "".join(char if char.isalnum() or char == "_" else "_" for char in name)
    return safe or "unnamed"


def _input_binding_keys_for_bundle_files(bundle_files: set[str]) -> set[str]:
    bindings = {
        "topology.spgt.h5": "input_h5_topology_path",
        "protocol.spgp.h5": "input_h5_protocol_path",
        "restart.spgr.h5": "input_h5_restart_path",
        "trajectory.spg.h5md": "input_h5_trajectory_path",
    }
    keys = {key for bundle_file, key in bindings.items() if bundle_file in bundle_files}
    if "restart.spgr.h5" in bundle_files:
        keys.add("input_h5_restart_load")
    if "trajectory.spg.h5md" in bundle_files:
        keys.add("input_h5_trajectory_particle_stream")
    return keys


def _particle_completion_metadata(bundle_path: Path) -> tuple[int, int, float]:
    if not bundle_path.exists():
        return 0, -1, 0.0
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise ConversionError("h5py is required to finalize particle H5MD metadata") from exc

    with h5py.File(bundle_path, "r") as handle:
        if "/particles/all/step" not in handle:
            return 0, -1, 0.0
        steps = np.asarray(handle["/particles/all/step"][...])
        if steps.size == 0:
            return 0, -1, 0.0
        frame_count = int(steps.shape[0])
        last_step = int(steps[-1])
        last_time = 0.0
        if "/particles/all/time" in handle:
            times = np.asarray(handle["/particles/all/time"][...], dtype=np.float64)
            if times.size:
                last_time = float(times[-1])
        return frame_count, last_step, last_time


def _parse_typed_payload(contract: IOContract, key: str, source_path: Path, *, atom_count: int | None = None):
    if contract.component == "topology":
        return parse_topology_file(key, source_path)
    if contract.component in {"protocol", "restart"}:
        return parse_protocol_or_restart_file(key, source_path)
    if contract.component == "trajectory":
        if atom_count is None:
            raise ConversionError("atom count is required for trajectory conversion")
        return parse_trajectory_file(key, source_path, atom_count=atom_count)
    return None


def _has_typed_state_parser(key: str, source_path: Path) -> bool:
    if not source_path.exists():
        return False
    return parse_protocol_or_restart_file(key, source_path) is not None


def _looks_text_file(path: Path) -> bool:
    try:
        data = path.read_bytes()[:4096]
    except OSError:
        return False
    if b"\0" in data:
        return False
    try:
        data.decode("utf-8")
    except UnicodeDecodeError:
        return False
    return True
