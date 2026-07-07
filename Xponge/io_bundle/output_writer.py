"""
Legacy SPONGE output files to H5MD bundle writer.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .h5_writer import _import_h5py, ensure_dataset, ensure_group, ensure_hard_link, set_attrs, write_string_array
from .legacy_case import LegacyCase, scan_legacy_case
from .manifest import ConversionManifest, ManifestEntry
from .output_parsers import parse_legacy_output_file


_SCHEMA_VERSION = "xponge.legacy_to_bundle.v1"


class LegacyOutputConversionError(RuntimeError):
    """Raised when legacy outputs cannot be materialized as an H5MD bundle."""


class LegacyOutputBundleWriter:
    """Convert existing legacy SPONGE output files into one H5MD output bundle."""

    def __init__(self, case: LegacyCase, output_path: str | Path, *, atom_count: int | None = None):
        self.case = case
        self.output_path = Path(output_path).resolve()
        self.atom_count = atom_count
        self.manifest = ConversionManifest(case_root=str(self.case.root), mode=self.case.mode)
        self._particle_values_by_bundle: dict[Path, set[str]] = {}
        self._observable_values_by_bundle: dict[Path, set[str]] = {}
        self._restart_particle_values: set[str] = set()
        self._wrote_any_restart_dataset = False
        self._finalized_h5md_bundles: set[Path] = set()
        self._legacy_sources_by_bundle: dict[Path, list[tuple[str, str]]] = {}
        self._vds_bundles: set[Path] = set()

    def convert(self, *, dry_run: bool = False) -> ConversionManifest:
        typed_datasets_by_bundle: dict[Path, list] = {}
        restart_datasets = []
        for key in _LEGACY_OUTPUT_PARSE_KEYS:
            source_path = self._resolve_output_path(key)
            if source_path is None:
                continue
            if not source_path.exists():
                self._record_missing_output(key, source_path)
                continue
            output_datasets = self._parse_output(key, source_path)
            for dataset in output_datasets:
                if key in _RESTART_OUTPUT_KEYS:
                    self._wrote_any_restart_dataset = True
                    self._track_restart_dataset_path(dataset.path)
                    restart_datasets.append(dataset)
                else:
                    output_path = self._dataset_output_path(dataset.path)
                    self._track_dataset_path(output_path, dataset.path)
                    typed_datasets_by_bundle.setdefault(output_path, []).append(dataset)
            self._record_legacy_source_for_bundle(self._recorded_bundle_path_for_key(key), key, source_path)
            self._record_converted_output(key, source_path)

        if not dry_run:
            for output_path, typed_datasets in typed_datasets_by_bundle.items():
                if self._use_vds_output_for(output_path):
                    self._write_vds_bundle(output_path, typed_datasets)
                else:
                    for dataset in typed_datasets:
                        ensure_dataset(output_path, dataset.path, dataset.data)
        if not dry_run and restart_datasets:
            self._write_restart_bundle(restart_datasets)

        if not dry_run:
            for output_path in typed_datasets_by_bundle:
                self._finalize_h5md_layout(output_path)
        if not dry_run and self._wrote_any_restart_dataset:
            self._finalize_restart_layout()
        self.manifest.bundled_mdin = None
        if not dry_run:
            self.manifest.write(self.output_path.with_suffix(self.output_path.suffix + ".manifest.json"))
        return self.manifest

    def _parse_output(self, key: str, source_path: Path):
        try:
            return parse_legacy_output_file(
                key,
                source_path,
                atom_count=self._infer_atom_count() if key in _ATOM_COUNT_KEYS else None,
                sits_nk_width=self._infer_sits_nk_width() if key in _SITS_NK_WIDTH_KEYS else None,
            )
        except ValueError as exc:
            raise LegacyOutputConversionError(str(exc)) from exc

    def _resolve_output_path(self, key: str) -> Path | None:
        explicit = self.case.commands.get(key)
        if explicit:
            return self.case.resolve_value_path(explicit)
        prefix = self.case.commands.get("default_out_file_prefix")
        if prefix and key in _DEFAULT_OUTPUT_SUFFIXES:
            candidate = self.case.resolve_value_path(prefix + _DEFAULT_OUTPUT_SUFFIXES[key])
            if candidate.exists():
                return candidate
        if prefix and key in _DEFAULT_OUTPUT_PREFIXED_NAMES:
            candidate = self.case.resolve_value_path(_DEFAULT_OUTPUT_PREFIXED_NAMES[key].format(prefix=prefix))
            if candidate.exists():
                return candidate
        if key in _HARDCODED_OUTPUT_FILES:
            candidate = self.case.root / _HARDCODED_OUTPUT_FILES[key]
            return candidate if candidate.exists() else None
        return None

    def _use_vds_output(self) -> bool:
        value = self.case.commands.get("output_h5_trajectory_vds", "")
        return value.strip().lower() in {"1", "true", "yes", "on"}

    def _use_vds_output_for(self, output_path: Path) -> bool:
        return (
            self._use_vds_output()
            and output_path == self.output_path
            and bool(self._particle_values_by_bundle.get(output_path))
        )

    def _trajectory_chunk_size(self) -> int:
        raw_value = self.case.commands.get("output_h5_trajectory_chunk_size", "")
        if not raw_value:
            return 20
        chunk_size = int(raw_value)
        if chunk_size <= 0:
            raise LegacyOutputConversionError("output_h5_trajectory_chunk_size must be positive")
        return chunk_size

    def _restart_output_path(self) -> Path:
        configured = self.case.commands.get("output_h5_restart_path")
        if configured:
            return self.case.resolve_value_path(configured)
        return self.output_path.with_name("restart.spgr.h5")

    def _observable_output_path(self) -> Path | None:
        configured = self.case.commands.get("output_h5_observable_path")
        if configured:
            return self.case.resolve_value_path(configured)
        return None

    def _dataset_output_path(self, dataset_path: str) -> Path:
        observable_path = self._observable_output_path()
        if observable_path is not None and not dataset_path.startswith("/particles/all/"):
            return observable_path
        return self.output_path

    def _infer_atom_count(self) -> int:
        if self.atom_count is not None:
            return self.atom_count
        for key in ("mass_in_file", "charge_in_file"):
            source_path = self.case.resolve_legacy_input_path(key)
            if source_path is None or not source_path.exists():
                continue
            lines = [line.split() for line in source_path.read_text(encoding="utf-8").splitlines() if line.strip()]
            if lines:
                self.atom_count = int(lines[0][0])
                return self.atom_count
        coordinate_path = self.case.resolve_legacy_input_path("coordinate_in_file")
        if coordinate_path is not None and coordinate_path.exists():
            lines = [line.split() for line in coordinate_path.read_text(encoding="utf-8").splitlines() if line.strip()]
            if lines:
                self.atom_count = int(lines[0][0])
                return self.atom_count
        raise LegacyOutputConversionError(
            "cannot infer atom count needed to convert legacy vector trajectory outputs"
        )

    def _infer_sits_nk_width(self) -> int:
        restart_path = self._resolve_output_path("SITS_nk_rest_file")
        if restart_path is not None and restart_path.exists():
            values = [token for line in restart_path.read_text(encoding="utf-8").splitlines() for token in line.split()]
            if values:
                return len(values)
        input_path = self.case.resolve_legacy_input_path("SITS_nk_in_file")
        if input_path is not None and input_path.exists():
            values = [token for line in input_path.read_text(encoding="utf-8").splitlines() for token in line.split()]
            if values:
                return len(values)
        raise LegacyOutputConversionError(
            "cannot infer SITS nk width needed to convert SITS_nk_traj_file"
        )

    def _track_dataset_path(self, output_path: Path, dataset_path: str) -> None:
        if dataset_path.startswith("/particles/all/") and dataset_path.endswith("/value"):
            self._particle_values_by_bundle.setdefault(output_path, set()).add(dataset_path)
        if dataset_path.startswith("/observables/all/") and dataset_path.endswith("/value"):
            self._observable_values_by_bundle.setdefault(output_path, set()).add(dataset_path)

    def _track_restart_dataset_path(self, dataset_path: str) -> None:
        if dataset_path.startswith("/particles/all/") and dataset_path.endswith("/value"):
            self._restart_particle_values.add(dataset_path)

    def _record_converted_output(self, key: str, source_path: Path) -> None:
        bundle_path = self._recorded_bundle_path_for_key(key)
        self.manifest.add(
            ManifestEntry(
                contract_id=f"output.legacy_import.{key}",
                status="typed_converted",
                source_key=key,
                source_path=str(source_path),
                bundle_file=bundle_path.name,
                bundle_path=_OUTPUT_BUNDLE_PATHS[key],
                direction="output",
                component="output",
                payload_kind="file",
                override_policy="post_run_import",
                comparison_rule="typed",
                message="legacy output file materialized as H5MD/SPONGE extension datasets",
            )
        )

    def _record_legacy_source_for_bundle(self, bundle_path: Path, key: str, source_path: Path) -> None:
        rel_path = _relative_h5_path(source_path.resolve(), self.case.root.resolve())
        sources = self._legacy_sources_by_bundle.setdefault(bundle_path, [])
        entry = (key, rel_path)
        if entry not in sources:
            sources.append(entry)

    def _record_missing_output(self, key: str, source_path: Path) -> None:
        bundle_path = self._recorded_bundle_path_for_key(key)
        self.manifest.add(
            ManifestEntry(
                contract_id=f"output.legacy_import.{key}",
                status="legacy_output_missing",
                source_key=key,
                source_path=str(source_path),
                bundle_file=bundle_path.name,
                bundle_path=_OUTPUT_BUNDLE_PATHS[key],
                direction="output",
                component="output",
                payload_kind="file",
                override_policy="post_run_import",
                comparison_rule="typed",
                message="legacy output path is configured but the file is absent",
            )
        )

    def _recorded_bundle_path_for_key(self, key: str) -> Path:
        if key in _RESTART_OUTPUT_KEYS:
            return self._restart_output_path()
        if key in _PARTICLE_OUTPUT_KEYS:
            return self.output_path
        observable_path = self._observable_output_path()
        return observable_path if observable_path is not None else self.output_path

    def _finalize_h5md_layout(self, output_path: Path) -> None:
        ensure_group(output_path, "/h5md")
        ensure_group(output_path, "/h5md/creator")
        ensure_group(output_path, "/parameters/sponge")
        set_attrs(output_path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        set_attrs(output_path, "/h5md/creator", {"name": "XPONGE", "version": "legacy-output-import"})
        ensure_dataset(output_path, "/parameters/sponge/schema/name", "sponge.output.h5md")
        ensure_dataset(output_path, "/parameters/sponge/schema/version", _SCHEMA_VERSION)
        ensure_dataset(output_path, "/parameters/sponge/output/status", "finalized")
        self._write_legacy_source_table(output_path)
        self._write_output_completion(output_path)
        self._write_stream_tables(output_path)
        ensure_dataset(
            output_path,
            "/parameters/sponge/output/mode",
            "vds" if output_path in self._vds_bundles else "legacy_import",
        )

        if self._particle_values_by_bundle.get(output_path):
            self._link_particle_axes(output_path)
        if self._observable_values_by_bundle.get(output_path):
            self._link_observable_axes(output_path)
        self._record_generated_h5md_bundle(output_path)

    def _record_generated_h5md_bundle(self, output_path: Path) -> None:
        if output_path in self._finalized_h5md_bundles:
            return
        self._finalized_h5md_bundles.add(output_path)
        is_trajectory = bool(self._particle_values_by_bundle.get(output_path))
        contract_id = "output.h5.trajectory" if is_trajectory else "output.h5.observable"
        self.manifest.add(
            ManifestEntry(
                contract_id=contract_id,
                status="generated",
                bundle_file=output_path.name,
                bundle_path="/",
                direction="output",
                component="output",
                payload_kind="metadata",
                override_policy="output_h5_trajectory_path" if is_trajectory else "output_h5_observable_path",
                comparison_rule="h5md",
                message="legacy output datasets were materialized as an H5MD bundle",
            )
        )

    def _link_particle_axes(self, output_path: Path) -> None:
        self._link_particle_axes_for_values(output_path, self._particle_values_by_bundle[output_path])

    def _link_particle_axes_for_values(self, output_path: Path, value_paths: set[str]) -> None:
        if _h5_path_exists(output_path, "/particles/all/time"):
            set_attrs(output_path, "/particles/all/time", {"unit": "ps"})
        for value_path in sorted(value_paths):
            if not _h5_path_exists(output_path, value_path):
                continue
            if value_path.endswith("/position/value"):
                set_attrs(output_path, value_path, {"unit": "Angstrom"})
            elif value_path.endswith("/velocity/value"):
                set_attrs(output_path, value_path, {"unit": "Angstrom ps-1"})
            elif value_path.endswith("/force/value"):
                set_attrs(output_path, value_path, {"unit": "kcal mol-1 Angstrom-1"})
            elif value_path.endswith("/box/edges/value"):
                set_attrs(output_path, value_path, {"unit": "Angstrom"})
            parent = str(Path(value_path).parent).replace("\\", "/")
            if _h5_path_exists(output_path, "/particles/all/step"):
                ensure_hard_link(output_path, "/particles/all/step", parent + "/step")
            if _h5_path_exists(output_path, "/particles/all/time"):
                ensure_hard_link(output_path, "/particles/all/time", parent + "/time")

    def _link_observable_axes(self, output_path: Path) -> None:
        self._link_observable_axes_for_values(output_path, self._observable_values_by_bundle[output_path])

    def _link_observable_axes_for_values(self, output_path: Path, value_paths: set[str]) -> None:
        if _h5_path_exists(output_path, "/observables/all/time"):
            set_attrs(output_path, "/observables/all/time", {"unit": "ps"})
        for value_path in sorted(value_paths):
            if not _h5_path_exists(output_path, value_path):
                continue
            parent = str(Path(value_path).parent).replace("\\", "/")
            if not _h5_path_exists(output_path, parent + "/step") and _h5_path_exists(
                output_path, "/observables/all/step"
            ):
                ensure_hard_link(output_path, "/observables/all/step", parent + "/step")
            if not _h5_path_exists(output_path, parent + "/time") and _h5_path_exists(
                output_path, "/observables/all/time"
            ):
                ensure_hard_link(output_path, "/observables/all/time", parent + "/time")

    def _write_restart_bundle(self, restart_datasets) -> None:
        restart_path = self._restart_output_path()
        for dataset in restart_datasets:
            ensure_dataset(restart_path, dataset.path, dataset.data)

    def _finalize_restart_layout(self) -> None:
        restart_path = self._restart_output_path()
        ensure_group(restart_path, "/run")
        ensure_group(restart_path, "/parameters/restart")
        ensure_group(restart_path, "/parameters/sponge")
        ensure_dataset(restart_path, "/parameters/sponge/schema/name", "sponge.restart.h5")
        ensure_dataset(restart_path, "/parameters/sponge/schema/version", _SCHEMA_VERSION)
        ensure_dataset(restart_path, "/parameters/sponge/output/status", "finalized")
        self._write_legacy_source_table(restart_path)
        self._write_output_completion(restart_path)
        if _h5_path_exists(restart_path, "/particles/all/time"):
            set_attrs(restart_path, "/particles/all/time", {"unit": "ps"})
        for value_path in sorted(self._restart_particle_values):
            if value_path.endswith("/position/value"):
                set_attrs(restart_path, value_path, {"unit": "Angstrom"})
            elif value_path.endswith("/velocity/value"):
                set_attrs(restart_path, value_path, {"unit": "Angstrom ps-1"})
            elif value_path.endswith("/box/edges/value"):
                set_attrs(restart_path, value_path, {"unit": "Angstrom"})
            parent = str(Path(value_path).parent).replace("\\", "/")
            if _h5_path_exists(restart_path, "/particles/all/step"):
                ensure_hard_link(restart_path, "/particles/all/step", parent + "/step")
            if _h5_path_exists(restart_path, "/particles/all/time"):
                ensure_hard_link(restart_path, "/particles/all/time", parent + "/time")
        self.manifest.add(
            ManifestEntry(
                contract_id="output.h5.restart",
                status="generated",
                bundle_file=restart_path.name,
                bundle_path="/parameters/restart",
                direction="output",
                component="restart",
                payload_kind="metadata",
                override_policy="output_h5_restart_path",
                comparison_rule="h5_restart",
                message="legacy restart output state was materialized as restart.spgr.h5",
            )
        )

    def _write_legacy_source_table(self, bundle_path: Path) -> None:
        sources = self._legacy_sources_by_bundle.get(bundle_path, [])
        if not sources:
            return
        keys = [key for key, _ in sources]
        paths = [path for _, path in sources]
        write_string_array(bundle_path, "/parameters/sponge/files/legacy_sidecars/key", keys)
        write_string_array(bundle_path, "/parameters/sponge/files/legacy_sidecars/path", paths)

    def _write_output_completion(self, bundle_path: Path) -> None:
        frame_count, last_step, last_time = _completion_metadata(bundle_path)
        ensure_dataset(
            bundle_path,
            "/parameters/sponge/output/frame_count",
            np.asarray([frame_count], dtype=np.int64),
        )
        ensure_dataset(
            bundle_path,
            "/parameters/sponge/output/last_complete_step",
            np.asarray([last_step], dtype=np.int64),
        )
        ensure_dataset(
            bundle_path,
            "/parameters/sponge/output/last_complete_time",
            np.asarray([last_time], dtype=np.float64),
        )

    def _write_stream_tables(self, output_path: Path) -> None:
        particle_streams = ["all"] if _h5_path_exists(output_path, "/particles/all") else []
        observable_streams = ["all"] if _h5_path_exists(output_path, "/observables/all") else []
        if particle_streams:
            write_string_array(output_path, "/parameters/sponge/output/particle_streams", particle_streams)
        if observable_streams:
            write_string_array(output_path, "/parameters/sponge/output/observable_streams", observable_streams)

    def _write_vds_bundle(self, output_path: Path, typed_datasets) -> None:
        h5py = _import_h5py()
        chunk_size = self._trajectory_chunk_size()

        vds_datasets = []
        direct_datasets = []
        for dataset in typed_datasets:
            if _is_vds_candidate(dataset.path, dataset.data):
                vds_datasets.append(dataset)
            else:
                direct_datasets.append(dataset)

        for dataset in direct_datasets:
            ensure_dataset(output_path, dataset.path, dataset.data)
        if not vds_datasets:
            return

        shard_dir = _vds_shard_dir(output_path)
        shard_dir.mkdir(parents=True, exist_ok=True)

        shard_entries: dict[int, dict[str, Any]] = {}
        shard_particle_values: dict[int, set[str]] = {}
        shard_observable_values: dict[int, set[str]] = {}
        for dataset in vds_datasets:
            data = np.asarray(dataset.data)
            frame_count = data.shape[0]
            dtype = data.dtype
            layout = h5py.VirtualLayout(shape=data.shape, dtype=dtype)
            for shard_index, frame_start in enumerate(range(0, frame_count, chunk_size)):
                frame_end = min(frame_start + chunk_size, frame_count)
                shard_path = shard_dir / f"segment_{shard_index:06d}.spg.h5md"
                source_data = data[frame_start:frame_end]
                ensure_dataset(shard_path, dataset.path, source_data)
                source = h5py.VirtualSource(
                    _relative_h5_path(shard_path, output_path.parent),
                    dataset.path,
                    shape=source_data.shape,
                )
                layout[frame_start:frame_end, ...] = source
                entry = shard_entries.setdefault(
                    shard_index,
                    {
                        "path": _relative_h5_path(shard_path, output_path.parent),
                        "frame_start": frame_start,
                        "frame_count": 0,
                    },
                )
                entry["frame_count"] = max(entry["frame_count"], frame_end - frame_start)
                if dataset.path.startswith("/particles/all/") and dataset.path.endswith("/value"):
                    shard_particle_values.setdefault(shard_index, set()).add(dataset.path)
                if dataset.path.startswith("/observables/all/") and dataset.path.endswith("/value"):
                    shard_observable_values.setdefault(shard_index, set()).add(dataset.path)
                ensure_dataset(shard_path, "/parameters/sponge/shard/index", np.asarray(shard_index, dtype=np.int64))
                ensure_dataset(shard_path, "/parameters/sponge/shard/frame_start", np.asarray(frame_start, dtype=np.int64))
                ensure_dataset(
                    shard_path,
                    "/parameters/sponge/shard/frame_count",
                    np.asarray(frame_end - frame_start, dtype=np.int64),
                )
                ensure_dataset(shard_path, "/parameters/sponge/shard/status", "complete")
            _create_virtual_dataset(output_path, dataset.path, layout)

        if shard_entries:
            mdout_metadata = [
                dataset
                for dataset in direct_datasets
                if dataset.path.startswith("/parameters/sponge/mdout/columns/")
            ]
            for shard_index, entry in shard_entries.items():
                shard_path = output_path.parent / entry["path"]
                if shard_observable_values.get(shard_index):
                    for dataset in mdout_metadata:
                        ensure_dataset(shard_path, dataset.path, dataset.data)
                self._finalize_vds_shard_layout(
                    shard_path,
                    shard_particle_values.get(shard_index, set()),
                    shard_observable_values.get(shard_index, set()),
                )
            self._vds_bundles.add(output_path)
            ordered = [shard_entries[index] for index in sorted(shard_entries)]
            shard_metadata = [
                _vds_shard_metadata(output_path.parent / entry["path"])
                for entry in ordered
            ]
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/trajectory_chunk_size",
                np.asarray(chunk_size, dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/index",
                np.asarray(sorted(shard_entries), dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/path",
                np.asarray([entry["path"] for entry in ordered], dtype=object),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/frame_start",
                np.asarray([entry["frame_start"] for entry in ordered], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/frame_count",
                np.asarray([entry["frame_count"] for entry in ordered], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/observable_frame_count",
                np.asarray([metadata["observable_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/nhc_frame_count",
                np.asarray([metadata["nhc_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/sits_nk_frame_count",
                np.asarray([metadata["sits_nk_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/metadynamics_scalar_frame_count",
                np.asarray([metadata["metadynamics_scalar_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/qc_frame_count",
                np.asarray([metadata["qc_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/reaxff_frame_count",
                np.asarray([metadata["reaxff_frame_count"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/step_start",
                np.asarray([metadata["step_start"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/step_end",
                np.asarray([metadata["step_end"] for metadata in shard_metadata], dtype=np.int64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/time_start",
                np.asarray([metadata["time_start"] for metadata in shard_metadata], dtype=np.float64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/time_end",
                np.asarray([metadata["time_end"] for metadata in shard_metadata], dtype=np.float64),
            )
            ensure_dataset(
                output_path,
                "/parameters/sponge/output/shard_manifest/status",
                np.asarray(["complete"] * len(ordered), dtype=object),
            )
            self.manifest.add(
                ManifestEntry(
                    contract_id="output.h5.vds",
                    status="generated",
                    bundle_file=output_path.name,
                    bundle_path="/parameters/sponge/output/shard_manifest",
                    direction="output",
                    component="output",
                    payload_kind="metadata",
                    override_policy="output_h5_trajectory_vds",
                    comparison_rule="h5md_vds",
                    message="legacy output time series were materialized as HDF5 VDS shards",
                )
            )

    def _finalize_vds_shard_layout(
        self,
        shard_path: Path,
        particle_values: set[str],
        observable_values: set[str],
    ) -> None:
        ensure_group(shard_path, "/h5md")
        ensure_group(shard_path, "/h5md/creator")
        ensure_group(shard_path, "/parameters/sponge")
        set_attrs(shard_path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        set_attrs(shard_path, "/h5md/creator", {"name": "XPONGE", "version": "legacy-output-import"})
        ensure_dataset(shard_path, "/parameters/sponge/schema/name", "sponge.output.h5md")
        ensure_dataset(shard_path, "/parameters/sponge/schema/version", _SCHEMA_VERSION)
        ensure_dataset(shard_path, "/parameters/sponge/output/status", "finalized")
        self._write_vds_shard_completion(shard_path)
        self._write_output_completion(shard_path)
        self._write_stream_tables(shard_path)
        if particle_values:
            self._link_particle_axes_for_values(shard_path, particle_values)
        if observable_values:
            self._link_observable_axes_for_values(shard_path, observable_values)

    def _write_vds_shard_completion(self, shard_path: Path) -> None:
        frame_count, last_step, last_time = _completion_metadata(shard_path)
        ensure_dataset(
            shard_path,
            "/parameters/sponge/shard/frame_count",
            np.asarray(frame_count, dtype=np.int64),
        )
        ensure_dataset(
            shard_path,
            "/parameters/sponge/shard/last_complete_step",
            np.asarray(last_step, dtype=np.int64),
        )
        ensure_dataset(
            shard_path,
            "/parameters/sponge/shard/last_complete_time",
            np.asarray(last_time, dtype=np.float64),
        )


def convert_legacy_outputs_to_bundle(
    case_root: str | Path,
    output_path: str | Path,
    *,
    mdin: str | Path = "mdin.spg.toml",
    atom_count: int | None = None,
    dry_run: bool = False,
) -> ConversionManifest:
    """Convert existing legacy SPONGE outputs into an H5MD output bundle."""

    case = scan_legacy_case(case_root, mdin=mdin)
    return LegacyOutputBundleWriter(case, output_path, atom_count=atom_count).convert(dry_run=dry_run)


_LEGACY_OUTPUT_PARSE_KEYS = (
    "mdinfo",
    "mdout",
    "crd",
    "box",
    "vel",
    "frc",
    "nose_hoover_chain_crd",
    "nose_hoover_chain_vel",
    "SITS_nk_traj_file",
    "SITS_nk_rest_file",
    "meta_potential_out_file",
    "meta_direct_export",
    "meta_hills_log",
    "meta_history_log",
    "meta_edge_log",
    "qc_scf_output",
    "reaxff_eeq_charges",
    "rst",
    "nose_hoover_chain_restart_output",
)

_ATOM_COUNT_KEYS = {"crd", "vel", "frc"}
_SITS_NK_WIDTH_KEYS = {"SITS_nk_traj_file"}
_RESTART_OUTPUT_KEYS = {"rst", "nose_hoover_chain_restart_output"}
_PARTICLE_OUTPUT_KEYS = {"crd", "box", "vel", "frc"}

_DEFAULT_OUTPUT_SUFFIXES = {
    "mdinfo": ".info",
    "mdout": ".out",
    "crd": ".dat",
    "box": ".box",
}

_DEFAULT_OUTPUT_PREFIXED_NAMES = {
    "meta_potential_out_file": "{prefix}_Meta_Potential.txt",
}

_HARDCODED_OUTPUT_FILES = {
    "meta_potential_out_file": "Meta_Potential.txt",
    "meta_direct_export": "Meta_directly.txt",
    "meta_hills_log": "myhill.log",
    "meta_history_log": "history.log",
    "meta_edge_log": "sumhill.log",
    "reaxff_eeq_charges": "eeq_charges.txt",
}

_OUTPUT_BUNDLE_PATHS = {
    "mdinfo": "/parameters/sponge/log/mdinfo_text",
    "mdout": "/observables/all",
    "crd": "/particles/all/position/value",
    "box": "/particles/all/box/edges/value",
    "vel": "/particles/all/velocity/value",
    "frc": "/particles/all/force/value",
    "nose_hoover_chain_crd": "/observables/all/thermostat/nose_hoover_chain/coordinate/value",
    "nose_hoover_chain_vel": "/observables/all/thermostat/nose_hoover_chain/velocity/value",
    "SITS_nk_traj_file": "/observables/all/sits/SITS/nk/value",
    "SITS_nk_rest_file": "/parameters/sponge/restart_exports/sits/SITS/nk/value",
    "meta_potential_out_file": "/parameters/sponge/metadynamics/default/potential_export",
    "meta_direct_export": "/parameters/sponge/metadynamics/default/direct_export",
    "meta_hills_log": "/parameters/sponge/metadynamics/default/hills",
    "meta_history_log": "/parameters/sponge/metadynamics/default/history",
    "meta_edge_log": "/parameters/sponge/metadynamics/default/edge",
    "qc_scf_output": "/parameters/sponge/qc/scf_output",
    "reaxff_eeq_charges": "/parameters/sponge/reaxff/eeq_charges/value",
    "rst": "/particles/all/position/value",
    "nose_hoover_chain_restart_output": "/parameters/restart/thermostat/nose_hoover_chain",
}


def _h5_path_exists(file_path: Path, object_path: str) -> bool:
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise LegacyOutputConversionError("h5py is required to finalize output H5MD links") from exc
    with h5py.File(file_path, "r") as handle:
        return object_path in handle


def _completion_metadata(file_path: Path) -> tuple[int, int, float]:
    h5py = _import_h5py()
    with h5py.File(file_path, "r") as handle:
        for root in ("/particles/all", "/observables/all"):
            step_path = f"{root}/step"
            time_path = f"{root}/time"
            if step_path not in handle:
                continue
            steps = np.asarray(handle[step_path][...])
            if steps.size == 0:
                continue
            frame_count = int(steps.shape[0])
            last_step = int(steps[-1])
            last_time = 0.0
            if time_path in handle:
                times = np.asarray(handle[time_path][...], dtype=np.float64)
                if times.size:
                    last_time = float(times[-1])
                return frame_count, last_step, last_time
    return 0, -1, 0.0


def _vds_shard_metadata(file_path: Path) -> dict[str, int | float]:
    h5py = _import_h5py()
    with h5py.File(file_path, "r") as handle:
        step_start, step_end, time_start, time_end = _h5_series_extent(
            handle,
            [
                "/particles/all/step",
                "/observables/all/step",
                "/observables/all/thermostat/nose_hoover_chain/coordinate/step",
                "/observables/all/sits/SITS/nk/step",
            ],
        )
        _, _, fallback_time_start, fallback_time_end = _h5_series_extent(
            handle,
            [
                "/particles/all/time",
                "/observables/all/time",
                "/observables/all/thermostat/nose_hoover_chain/coordinate/time",
                "/observables/all/sits/SITS/nk/time",
            ],
        )
        if time_start == 0.0 and time_end == 0.0:
            time_start = fallback_time_start
            time_end = fallback_time_end
        return {
            "observable_frame_count": _h5_frame_count(handle, ["/observables/all/step"]),
            "nhc_frame_count": _h5_frame_count(
                handle,
                [
                    "/observables/all/thermostat/nose_hoover_chain/coordinate/step",
                    "/observables/all/thermostat/nose_hoover_chain/velocity/step",
                ],
            ),
            "sits_nk_frame_count": _h5_frame_count(handle, ["/observables/all/sits/SITS/nk/step"]),
            "metadynamics_scalar_frame_count": _h5_child_step_frame_count(
                handle,
                "/observables/all/metadynamics",
            ),
            "qc_frame_count": _h5_frame_count(
                handle,
                [
                    "/observables/all/qc/energy/step",
                    "/observables/all/qc/spin_square/step",
                ],
            ),
            "reaxff_frame_count": _h5_child_step_frame_count(handle, "/observables/all/reaxff"),
            "step_start": step_start,
            "step_end": step_end,
            "time_start": time_start,
            "time_end": time_end,
        }


def _h5_series_extent(handle, paths: list[str]) -> tuple[int, int, float, float]:
    for path in paths:
        if path not in handle:
            continue
        data = np.asarray(handle[path][...])
        if data.size == 0:
            continue
        first = data.reshape(-1)[0]
        last = data.reshape(-1)[-1]
        if np.issubdtype(data.dtype, np.integer):
            return int(first), int(last), 0.0, 0.0
        return -1, -1, float(first), float(last)
    return -1, -1, 0.0, 0.0


def _h5_frame_count(handle, paths: list[str]) -> int:
    frame_count = 0
    for path in paths:
        if path not in handle:
            continue
        data = np.asarray(handle[path][...])
        if data.ndim == 0:
            continue
        frame_count = max(frame_count, int(data.shape[0]))
    return frame_count


def _h5_child_step_frame_count(handle, root_path: str) -> int:
    if root_path not in handle:
        return 0
    frame_count = 0
    group = handle[root_path]
    for name in group:
        step_path = f"{root_path}/{name}/step"
        if step_path in handle:
            data = np.asarray(handle[step_path][...])
            if data.ndim > 0:
                frame_count = max(frame_count, int(data.shape[0]))
    return frame_count


def _is_vds_candidate(dataset_path: str, data) -> bool:
    if not (
        dataset_path.startswith("/particles/all/")
        or dataset_path.startswith("/observables/all/")
    ):
        return False
    array = np.asarray(data)
    if array.ndim == 0 or array.dtype.kind in {"O", "U", "S"}:
        return False
    return True


def _create_virtual_dataset(file_path: Path, dataset_path: str, layout) -> None:
    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a", libver="latest") as handle:
        if dataset_path in handle:
            del handle[dataset_path]
        parent = str(Path(dataset_path).parent).replace("\\", "/")
        if parent != ".":
            handle.require_group(parent)
        handle.create_virtual_dataset(dataset_path, layout, fillvalue=0)


def _vds_shard_dir(output_path: Path) -> Path:
    name = output_path.name
    if name.endswith(".spg.h5md"):
        shard_name = name[: -len(".spg.h5md")] + ".spg.shards"
    elif name.endswith(".h5md"):
        shard_name = name[: -len(".h5md")] + ".shards"
    else:
        shard_name = name + ".shards"
    return output_path.with_name(shard_name)


def _relative_h5_path(path: Path, base: Path) -> str:
    try:
        return path.relative_to(base).as_posix()
    except ValueError:
        return path.as_posix()
