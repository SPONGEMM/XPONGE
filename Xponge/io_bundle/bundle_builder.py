"""Shared writer and finalizer for SPONGE bundled input artifacts."""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path
from typing import Any, Callable

import numpy as np

from . import h5_writer


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

    def add_datasets(self, bundle_file: str, datasets) -> None:
        path = self.paths.for_bundle_file(bundle_file)
        for dataset in datasets:
            self.io.ensure_dataset(path, dataset.path, dataset.data)
            self._track_dataset(bundle_file, dataset.path, dataset.data)

    def add_dataset(self, bundle_file: str, dataset_path: str, data: Any) -> None:
        path = self.paths.for_bundle_file(bundle_file)
        self.io.ensure_dataset(path, dataset_path, data)
        self._track_dataset(bundle_file, dataset_path, data)

    def _track_dataset(self, bundle_file: str, dataset_path: str, data: Any) -> None:
        if self.track_dataset_hashes:
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
        digest = hashlib.sha256()
        digest.update(bundle_file.encode("utf-8"))
        for key, value in sorted(self._dataset_hashes.items()):
            current_bundle, dataset_path = key.split(":", 1)
            if current_bundle != bundle_file:
                continue
            if path_prefixes and not dataset_path.startswith(path_prefixes):
                continue
            array = np.asarray(value)
            digest.update(b"\0")
            digest.update(dataset_path.encode("utf-8"))
            digest.update(b"\0")
            digest.update(str(array.dtype).encode("ascii"))
            digest.update(b"\0")
            digest.update(repr(array.shape).encode("ascii"))
            digest.update(b"\0")
            if array.dtype.kind in {"O", "U", "S"}:
                digest.update("\0".join(map(str, array.reshape(-1))).encode("utf-8"))
            else:
                digest.update(np.ascontiguousarray(array).tobytes())
        return "sha256:" + digest.hexdigest()

    def finalize(self, touched: set[str], metadata: BundleMetadata) -> None:
        schema_version = "xponge.legacy_to_bundle.v1"
        if "topology.spgt.h5" in touched:
            topology = self.paths.topology
            self.io.write_string(topology, "/schema/name", "sponge.topology.h5")
            self.io.write_string(topology, "/schema/version", schema_version)
            self.io.write_string(topology, "/parameters/sponge/schema/name", "sponge.topology.h5")
            self.io.write_string(topology, "/parameters/sponge/schema/version", schema_version)
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
            self.io.write_string(protocol, "/schema/version", schema_version)
            self.io.write_string(protocol, "/parameters/sponge/schema/name", "sponge.protocol.h5")
            self.io.write_string(protocol, "/parameters/sponge/schema/version", schema_version)
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
            self.io.write_string(restart, "/parameters/sponge/schema/name", "sponge.restart.h5")
            self.io.write_string(restart, "/parameters/sponge/schema/version", schema_version)
            self._finalize_restart(restart, metadata.creator_version)

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
                schema_version,
            )
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
        self.io.write_string_array(path, "/parameters/sponge/output/particle_streams", ["all"])

    def _finalize_restart(self, path: Path, creator_version: str) -> None:
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
