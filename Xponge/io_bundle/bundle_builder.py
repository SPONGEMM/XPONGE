"""Shared writer and finalizer for SPONGE bundled input artifacts."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any, Callable
from uuid import uuid4

import numpy as np

from . import h5_writer


_RESTART_PARTICLE_STATE_PATHS = (
    "/particles/all/step",
    "/particles/all/time",
    "/particles/all/position/value",
    "/particles/all/velocity/value",
    "/particles/all/box/edges/value",
    "/particles/all/force/value",
)
_CANONICAL_NUMERIC_DTYPES = {
    "float32",
    "float64",
    "int8",
    "uint8",
    "int32",
    "uint32",
    "int64",
    "uint64",
}


def canonical_dataset_hash(
    bundle_file: str,
    datasets: dict[str, Any],
    *,
    path_prefixes: tuple[str, ...] = (),
) -> str:
    """Hash logical HDF5 datasets using the SPONGE v2 canonical encoding."""

    digest = hashlib.sha256()
    digest.update(bundle_file.encode("utf-8"))
    for dataset_path, value in sorted(datasets.items()):
        if path_prefixes and not dataset_path.startswith(path_prefixes):
            continue
        array = np.asarray(value)
        string_dataset = array.dtype.kind in {"O", "U", "S"}
        digest.update(b"\0")
        digest.update(dataset_path.encode("utf-8"))
        digest.update(b"\0")
        # HDF5 variable-length UTF-8 strings are exposed as object arrays by
        # h5py and as std::string by SPONGE. Hash the logical type rather than
        # NumPy's source-width spelling, which is not preserved on disk.
        dtype_name = "object" if string_dataset else array.dtype.name
        # HDF5 boolean datasets are represented as an enum on disk. SPONGE's
        # canonical hash operates on native arithmetic buffers, so normalize
        # logical booleans to uint8 (0/1) instead of depending on an HDF5 enum
        # encoding or NumPy's private bool byte spelling.
        if dtype_name == "bool":
            array = array.astype(np.uint8, copy=False)
            dtype_name = "uint8"
        if not string_dataset and dtype_name not in _CANONICAL_NUMERIC_DTYPES:
            raise TypeError(
                f"unsupported canonical dataset dtype {array.dtype} at {dataset_path}"
            )
        digest.update(dtype_name.encode("ascii"))
        digest.update(b"\0")
        digest.update(repr(array.shape).encode("ascii"))
        digest.update(b"\0")
        if string_dataset:
            values = (
                item.decode("utf-8")
                if isinstance(item, (bytes, np.bytes_))
                else str(item)
                for item in array.reshape(-1)
            )
            digest.update("\0".join(values).encode("utf-8"))
        else:
            # SPONGE reads into native arithmetic types before hashing. Do the
            # same so file endianness does not change a logical state hash.
            native_dtype = np.dtype(dtype_name)
            digest.update(np.ascontiguousarray(array, dtype=native_dtype).tobytes())
    return "sha256:" + digest.hexdigest()


def restart_state_hash_from_h5(file_path: str | Path) -> str:
    """Recompute the canonical state digest from a materialized restart."""

    try:
        import h5py  # type: ignore
    except ImportError as exc:  # pragma: no cover - project dependency
        raise RuntimeError("h5py is required to hash restart bundles") from exc

    datasets: dict[str, Any] = {}
    with h5py.File(file_path, "r") as handle:
        for dataset_path in _RESTART_PARTICLE_STATE_PATHS:
            if dataset_path in handle:
                datasets[dataset_path] = _read_h5_dataset(handle[dataset_path], h5py)
        restart_root = "/parameters/restart"
        if restart_root in handle:
            _collect_h5_datasets(
                handle[restart_root], restart_root, datasets, h5py
            )
    return canonical_dataset_hash("restart.spgr.h5", datasets)


def _collect_h5_datasets(group, group_path: str, datasets: dict[str, Any], h5py) -> None:
    for name in sorted(group.keys()):
        child = group[name]
        child_path = f"{group_path}/{name}"
        if isinstance(child, h5py.Dataset):
            datasets[child_path] = _read_h5_dataset(child, h5py)
        elif isinstance(child, h5py.Group):
            _collect_h5_datasets(child, child_path, datasets, h5py)


def _read_h5_dataset(dataset, h5py):
    if h5py.check_string_dtype(dataset.dtype) is not None:
        return np.asarray(dataset.asstr()[...])
    return np.asarray(dataset[...])


@dataclass(frozen=True)
class BundlePaths:
    """Physical paths for the canonical bundled input artifact roles."""

    topology: Path
    protocol: Path
    restart: Path
    trajectory: Path | None = None

    @classmethod
    def canonical(cls, root: str | Path) -> "BundlePaths":
        root = Path(root)
        return cls(
            topology=root / "topology.spgt.h5",
            protocol=root / "protocol.spgp.h5",
            restart=root / "restart.spgr.h5",
            trajectory=root / "trajectory.spg.h5md",
        )

    def for_bundle_file(self, bundle_file: str) -> Path:
        path = {
            "topology.spgt.h5": self.topology,
            "protocol.spgp.h5": self.protocol,
            "restart.spgr.h5": self.restart,
            "trajectory.spg.h5md": self.trajectory,
        }.get(bundle_file)
        if path is None:
            raise KeyError(f"unknown or disabled bundled artifact: {bundle_file}")
        return path


@dataclass(frozen=True)
class BundleIO:
    """Injectable low-level HDF5 operations used by ``BundleBuilder``."""

    ensure_dataset: Callable = h5_writer.ensure_dataset
    ensure_group: Callable = h5_writer.ensure_group
    ensure_hard_link: Callable = h5_writer.ensure_hard_link
    set_attrs: Callable = h5_writer.set_attrs
    write_string: Callable = h5_writer.write_string
    write_bytes: Callable = h5_writer.write_bytes
    write_string_array: Callable = h5_writer.write_string_array


@dataclass(frozen=True)
class BundleMetadata:
    """Metadata required to finalize a set of bundled inputs."""

    atom_count: int
    topology_hash: str
    atom_order_hash: str
    forcefield_hash: str
    protocol_hash: str
    cv_count: int = 0
    restraint_count: int = 0
    enhanced_methods: tuple[str, ...] = ()
    creator_version: str = "legacy-to-bundle"


@dataclass
class BundleBuilder:
    """Write typed datasets and finalize canonical SPONGE HDF5 layouts."""

    paths: BundlePaths
    io: BundleIO = field(default_factory=BundleIO)
    track_dataset_hashes: bool = True
    _dataset_hashes: dict[str, Any] = field(default_factory=dict, init=False)
    identity_uuid: str = field(default_factory=lambda: str(uuid4()))

    def add_datasets(self, bundle_file: str, datasets) -> None:
        path = self.paths.for_bundle_file(bundle_file)
        for dataset in datasets:
            self.io.ensure_dataset(path, dataset.path, dataset.data)
            self._track_dataset(bundle_file, dataset.path, dataset.data)

    def add_dataset(self, bundle_file: str, dataset_path: str, data: Any) -> None:
        path = self.paths.for_bundle_file(bundle_file)
        self.io.ensure_dataset(path, dataset_path, data)
        self._track_dataset(bundle_file, dataset_path, data)

    def add_string(self, bundle_file: str, dataset_path: str, value: str) -> None:
        """Write a UTF-8 dataset and include it in canonical content hashes."""

        path = self.paths.for_bundle_file(bundle_file)
        self.io.write_string(path, dataset_path, value)
        self._track_dataset(bundle_file, dataset_path, value)

    def _track_dataset(self, bundle_file: str, dataset_path: str, data: Any) -> None:
        track_restart_state = bundle_file == "restart.spgr.h5" and (
            dataset_path.startswith("/particles/")
            or dataset_path.startswith("/parameters/restart/")
        )
        if self.track_dataset_hashes or track_restart_state:
            self._dataset_hashes[f"{bundle_file}:{dataset_path}"] = np.asarray(data)

    def write_legacy_sidecars(
        self,
        sidecars_by_bundle: dict[str, list[tuple[str, str]]],
    ) -> None:
        for bundle_file, sidecars in sidecars_by_bundle.items():
            if not sidecars:
                continue
            path = self.paths.for_bundle_file(bundle_file)
            self.io.write_string_array(
                path,
                "/parameters/sponge/files/legacy_sidecars/key",
                [key for key, _ in sidecars],
            )
            self.io.write_string_array(
                path,
                "/parameters/sponge/files/legacy_sidecars/path",
                [value for _, value in sidecars],
            )

    def content_hash(self, *, bundle_file: str, path_prefixes: tuple[str, ...] = ()) -> str:
        datasets = {}
        for key, value in self._dataset_hashes.items():
            current_bundle, dataset_path = key.split(":", 1)
            if current_bundle != bundle_file:
                continue
            datasets[dataset_path] = value
        return canonical_dataset_hash(
            bundle_file,
            datasets,
            path_prefixes=path_prefixes,
        )

    def finalize(self, touched: set[str], metadata: BundleMetadata) -> None:
        input_schema_version = "sponge.input.v2"
        output_schema_version = "sponge.output.v2"
        if "topology.spgt.h5" in touched:
            topology = self.paths.topology
            self.io.write_string(topology, "/schema/name", "sponge.topology.h5")
            self.io.write_string(topology, "/schema/version", input_schema_version)
            self.io.write_string(topology, "/parameters/sponge/schema/name", "sponge.topology.h5")
            self.io.write_string(topology, "/parameters/sponge/schema/version", input_schema_version)
            self.io.write_string(topology, "/identity/uuid", self.identity_uuid)
            self.io.write_string(topology, "/topology/atom_order_hash", metadata.atom_order_hash)
            self.io.write_string(topology, "/topology/topology_hash", metadata.topology_hash)
            self.io.write_string(topology, "/topology/forcefield_hash", metadata.forcefield_hash)
            self.io.ensure_dataset(
                topology,
                "/topology/atom_count",
                np.asarray(metadata.atom_count, dtype=np.int64),
            )
            try:
                self.io.set_attrs(topology, "/atoms/charge", {"unit": "Amber"})
            except KeyError:
                pass

        if "protocol.spgp.h5" in touched:
            protocol = self.paths.protocol
            self.io.write_string(protocol, "/schema/name", "sponge.protocol.h5")
            self.io.write_string(protocol, "/schema/version", input_schema_version)
            self.io.write_string(protocol, "/parameters/sponge/schema/name", "sponge.protocol.h5")
            self.io.write_string(protocol, "/parameters/sponge/schema/version", input_schema_version)
            self.io.write_string(protocol, "/identity/uuid", self.identity_uuid)
            self.io.write_string(
                protocol,
                "/protocol/topology_compatibility/topology_hash",
                metadata.topology_hash,
            )
            self.io.write_string(protocol, "/identity/content_hash", metadata.protocol_hash)
            self.io.ensure_dataset(
                protocol,
                "/protocol/cv_count",
                np.asarray(metadata.cv_count, dtype=np.int64),
            )
            self.io.ensure_dataset(
                protocol,
                "/protocol/restraint_count",
                np.asarray(metadata.restraint_count, dtype=np.int64),
            )
            if metadata.enhanced_methods:
                self.io.write_string(
                    protocol,
                    "/protocol/enhanced_sampling/method",
                    ",".join(metadata.enhanced_methods),
                )

        if "restart.spgr.h5" in touched:
            restart = self.paths.restart
            restart_state_hash = self.content_hash(
                bundle_file="restart.spgr.h5",
                path_prefixes=("/particles/", "/parameters/restart/"),
            )
            self.io.write_string(restart, "/parameters/sponge/schema/name", "sponge.restart.h5")
            self.io.write_string(restart, "/parameters/sponge/schema/version", input_schema_version)
            self.io.write_string(restart, "/schema/name", "sponge.restart.h5")
            self.io.write_string(restart, "/schema/version", input_schema_version)
            self.io.write_string(restart, "/identity/uuid", self.identity_uuid)
            self._finalize_restart(
                restart,
                metadata.creator_version,
                topology_hash=metadata.topology_hash,
                atom_order_hash=metadata.atom_order_hash,
                protocol_hash=metadata.protocol_hash,
                state_hash=restart_state_hash,
            )

        if "trajectory.spg.h5md" in touched:
            trajectory = self.paths.trajectory
            if trajectory is None:
                raise KeyError("trajectory path is disabled")
            self.io.write_string(
                trajectory,
                "/parameters/sponge/schema/name",
                "sponge.output.h5md",
            )
            self.io.write_string(
                trajectory,
                "/parameters/sponge/schema/version",
                output_schema_version,
            )
            self.io.write_string(trajectory, "/identity/uuid", self.identity_uuid)
            self._finalize_trajectory(trajectory, metadata)

    def _finalize_trajectory(self, path: Path, metadata: BundleMetadata) -> None:
        self.io.ensure_group(path, "/h5md")
        self.io.ensure_group(path, "/h5md/creator")
        self.io.ensure_group(path, "/particles/all")
        self.io.set_attrs(path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        self.io.set_attrs(
            path,
            "/h5md/creator",
            {"name": "XPONGE", "version": metadata.creator_version},
        )
        self.io.write_string(
            path,
            "/parameters/sponge/topology_compatibility/topology_hash",
            metadata.topology_hash,
        )
        self.io.write_string(
            path,
            "/parameters/sponge/topology_compatibility/atom_order_hash",
            metadata.atom_order_hash,
        )
        self._link_particle_axes(path)
        self.io.set_attrs(path, "/particles/all/time", {"unit": "ps"})
        frame_count, last_step, last_time = _particle_completion_metadata(path)
        self.io.write_string(path, "/parameters/sponge/output/mode", "single")
        self.io.write_string(path, "/parameters/sponge/output/status", "finalized")
        self._write_completion(path, frame_count, last_step, last_time)
        self.io.ensure_dataset(
            path,
            "/parameters/sponge/output/publication_epoch",
            np.asarray([0, 1], dtype=np.int64),
        )
        self.io.ensure_dataset(
            path,
            "/parameters/sponge/output/streams/particles/committed_count",
            np.asarray([frame_count], dtype=np.int64),
        )
        self.io.write_string(
            path,
            "/parameters/sponge/output/streams/particles/logical_kind",
            "trajectory_frames",
        )
        self.io.write_string(
            path,
            "/parameters/sponge/output/streams/particles/step_path",
            "/particles/all/step",
        )
        self.io.write_string(
            path,
            "/parameters/sponge/output/streams/particles/time_path",
            "/particles/all/time",
        )
        particle_values = [
            value_path
            for value_path in (
                "/particles/all/position/value",
                "/particles/all/box/edges/value",
                "/particles/all/velocity/value",
                "/particles/all/force/value",
            )
            if _h5_path_exists(path, value_path)
        ]
        self.io.write_string_array(
            path,
            "/parameters/sponge/output/streams/particles/value_paths",
            particle_values,
        )
        self.io.write_string(
            path,
            "/parameters/sponge/output/streams/particles/experimental",
            "false",
        )
        self.io.write_string_array(path, "/parameters/sponge/output/particle_streams", ["all"])

    def _finalize_restart(
        self,
        path: Path,
        creator_version: str,
        *,
        topology_hash: str,
        atom_order_hash: str,
        protocol_hash: str,
        state_hash: str,
    ) -> None:
        for group in (
            "/h5md",
            "/h5md/creator",
            "/run",
            "/particles/all",
            "/parameters/restart",
            "/parameters/restart/rng_state",
            "/parameters/restart/integrator_state",
            "/parameters/restart/thermostat",
            "/parameters/restart/barostat",
            "/parameters/restart/protocol_sidecars",
            "/parameters/restart/bias",
            "/parameters/restart/bias/sits",
            "/parameters/restart/bias/meta",
        ):
            self.io.ensure_group(path, group)
        self.io.set_attrs(path, "/h5md", {"version": np.asarray([1, 1], dtype=np.int32)})
        self.io.set_attrs(path, "/h5md/creator", {"name": "XPONGE", "version": creator_version})
        self.io.write_string(path, "/run/topology_hash", topology_hash)
        self.io.write_string(path, "/run/atom_order_hash", atom_order_hash)
        self.io.write_string(path, "/run/producer_protocol_hash", protocol_hash)
        self.io.write_string(path, "/run/state_hash", state_hash)
        self._link_particle_axes(path)
        try:
            self.io.set_attrs(path, "/particles/all/time", {"unit": "ps"})
        except KeyError:
            pass
        frame_count, last_step, last_time = _particle_completion_metadata(path)
        self.io.write_string(path, "/parameters/sponge/output/status", "finalized")
        self._write_completion(path, frame_count, last_step, last_time)
        if frame_count:
            self.io.ensure_dataset(path, "/run/current_step", np.asarray([last_step], dtype=np.int64))
            self.io.ensure_dataset(path, "/run/current_time", np.asarray([last_time], dtype=np.float64))
            self.io.write_string(path, "/run/state_type", "restart")
            self.io.write_string_array(path, "/parameters/sponge/output/particle_streams", ["all"])

    def _link_particle_axes(self, path: Path) -> None:
        for value_path, unit in (
            ("/particles/all/position/value", "Angstrom"),
            ("/particles/all/velocity/value", "Angstrom ps-1"),
            ("/particles/all/force/value", "kcal mol-1 Angstrom-1"),
            ("/particles/all/box/edges/value", "Angstrom"),
        ):
            try:
                self.io.set_attrs(path, value_path, {"unit": unit})
            except KeyError:
                continue
            parent = str(Path(value_path).parent).replace("\\", "/")
            self.io.ensure_hard_link(path, "/particles/all/step", parent + "/step")
            self.io.ensure_hard_link(path, "/particles/all/time", parent + "/time")
        try:
            self.io.set_attrs(
                path,
                "/particles/all/box",
                {"dimension": 3, "boundary": np.asarray(["periodic"] * 3, dtype="S8")},
            )
        except KeyError:
            pass

    def _write_completion(
        self,
        path: Path,
        frame_count: int,
        last_step: int,
        last_time: float,
    ) -> None:
        self.io.ensure_dataset(
            path,
            "/parameters/sponge/output/frame_count",
            np.asarray([frame_count], dtype=np.int64),
        )
        self.io.ensure_dataset(
            path,
            "/parameters/sponge/output/last_complete_step",
            np.asarray([last_step], dtype=np.int64),
        )
        self.io.ensure_dataset(
            path,
            "/parameters/sponge/output/last_complete_time",
            np.asarray([last_time], dtype=np.float64),
        )


def _particle_completion_metadata(path: Path) -> tuple[int, int, float]:
    if not path.exists():
        return 0, -1, 0.0
    try:
        import h5py  # type: ignore
    except ImportError as exc:  # pragma: no cover - project dependency
        raise RuntimeError("h5py is required to finalize particle H5MD metadata") from exc
    with h5py.File(path, "r") as handle:
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


def _h5_path_exists(path: Path, dataset_path: str) -> bool:
    if not path.exists():
        return False
    import h5py  # type: ignore

    with h5py.File(path, "r") as handle:
        return dataset_path in handle
