"""Validated HDF5 reader for bundled SPONGE input artifacts."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .bundle_case import BundleCase
from .errors import BundlePathError, BundleSchemaError, BundleValidationError


_INPUT_SCHEMA_VERSION = "sponge.input.v2"
_OUTPUT_SCHEMA_VERSION = "sponge.output.v2"
_SCHEMA_BINDINGS = {
    "topology.spgt.h5": ("/schema/name", "/schema/version", "sponge.topology.h5"),
    "protocol.spgp.h5": ("/schema/name", "/schema/version", "sponge.protocol.h5"),
    "restart.spgr.h5": (
        "/parameters/sponge/schema/name",
        "/parameters/sponge/schema/version",
        "sponge.restart.h5",
    ),
    "trajectory.spg.h5md": (
        "/parameters/sponge/schema/name",
        "/parameters/sponge/schema/version",
        "sponge.output.h5md",
    ),
}


class BundleReader:
    """Keep bundle HDF5 handles open while validating and exporting datasets."""

    def __init__(self, case: BundleCase, *, strict: bool = True):
        self.case = case
        self.strict = strict
        self._handles: dict[str, Any] = {}
        self.warnings: list[str] = []

    def __enter__(self) -> "BundleReader":
        try:
            import h5py  # type: ignore
        except ImportError as exc:  # pragma: no cover - project dependency
            raise BundleValidationError("h5py is required to read bundled inputs") from exc

        for bundle_file in _SCHEMA_BINDINGS:
            path = self.case.path_for_bundle_file(bundle_file)
            if path is not None and path.is_file():
                self._handles[bundle_file] = h5py.File(path, "r")
        try:
            self.validate()
        except Exception:
            self.close()
            raise
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def close(self) -> None:
        for handle in self._handles.values():
            handle.close()
        self._handles.clear()

    def has_bundle_file(self, bundle_file: str) -> bool:
        return bundle_file in self._handles

    def contains(self, bundle_file: str, dataset_path: str) -> bool:
        handle = self._handles.get(bundle_file)
        return handle is not None and dataset_path in handle

    def read(self, bundle_file: str, dataset_path: str):
        handle = self._require_handle(bundle_file)
        if dataset_path not in handle:
            raise BundleValidationError(f"{bundle_file} is missing {dataset_path}")
        return np.asarray(handle[dataset_path][...])

    def read_scalar(self, bundle_file: str, dataset_path: str):
        value = self.read(bundle_file, dataset_path)
        return value.item() if value.shape == () else value.reshape(-1)[0].item()

    def read_text(self, bundle_file: str, dataset_path: str) -> str:
        handle = self._require_handle(bundle_file)
        if dataset_path not in handle:
            raise BundleValidationError(f"{bundle_file} is missing {dataset_path}")
        return _decode_h5_text(handle[dataset_path][()])

    def read_legacy_sidecars(self, bundle_file: str) -> dict[str, Path]:
        key_path = "/parameters/sponge/files/legacy_sidecars/key"
        value_path = "/parameters/sponge/files/legacy_sidecars/path"
        if not self.contains(bundle_file, key_path) and not self.contains(bundle_file, value_path):
            return {}
        if not self.contains(bundle_file, key_path) or not self.contains(bundle_file, value_path):
            raise BundleValidationError(f"{bundle_file} has an incomplete legacy sidecar table")

        keys = [_decode_h5_text(value) for value in self.read(bundle_file, key_path).reshape(-1)]
        paths = [_decode_h5_text(value) for value in self.read(bundle_file, value_path).reshape(-1)]
        if len(keys) != len(paths):
            raise BundleValidationError(f"{bundle_file} legacy sidecar key/path lengths differ")

        sidecars: dict[str, Path] = {}
        root = self.case.root.resolve()
        for key, raw_path in zip(keys, paths):
            path = (root / raw_path).resolve()
            try:
                path.relative_to(root)
            except ValueError as exc:
                raise BundlePathError(
                    f"legacy sidecar for {key} escapes bundle root: {raw_path}"
                ) from exc
            if not path.is_file():
                if self.strict:
                    raise BundleValidationError(f"legacy sidecar for {key} does not exist: {path}")
                self.warnings.append(f"legacy sidecar for {key} does not exist: {path}")
                continue
            sidecars[key] = path
        return sidecars

    def validate(self) -> list[str]:
        for bundle_file, (name_path, version_path, expected_name) in _SCHEMA_BINDINGS.items():
            if bundle_file not in self._handles:
                continue
            self._validate_schema(bundle_file, name_path, version_path, expected_name)
            self.read_legacy_sidecars(bundle_file)
        self._validate_topology_protocol_hash()
        self._validate_atom_dimensions()
        self._validate_restart_load_policy()
        return list(self.warnings)

    def _validate_schema(
        self,
        bundle_file: str,
        name_path: str,
        version_path: str,
        expected_name: str,
    ) -> None:
        if not self.contains(bundle_file, name_path) or not self.contains(bundle_file, version_path):
            self._schema_problem(f"{bundle_file} is missing schema metadata")
            return
        name = self.read_text(bundle_file, name_path)
        version = self.read_text(bundle_file, version_path)
        if name != expected_name:
            self._schema_problem(
                f"{bundle_file} schema name is {name!r}, expected {expected_name!r}"
            )
        expected_version = (
            _OUTPUT_SCHEMA_VERSION if bundle_file == "trajectory.spg.h5md" else _INPUT_SCHEMA_VERSION
        )
        if version != expected_version:
            self._schema_problem(
                f"{bundle_file} schema version is {version!r}, expected {expected_version!r}"
            )

    def _schema_problem(self, message: str) -> None:
        if self.strict:
            raise BundleSchemaError(message)
        self.warnings.append(message)

    def _validate_topology_protocol_hash(self) -> None:
        topology_hash_path = "/topology/topology_hash"
        protocol_hash_path = "/protocol/topology_compatibility/topology_hash"
        if not self.contains("topology.spgt.h5", topology_hash_path) or not self.contains(
            "protocol.spgp.h5", protocol_hash_path
        ):
            return
        topology_hash = self.read_text("topology.spgt.h5", topology_hash_path)
        protocol_hash = self.read_text("protocol.spgp.h5", protocol_hash_path)
        if topology_hash != protocol_hash:
            raise BundleValidationError(
                "protocol topology compatibility hash does not match topology bundle"
            )

    def _validate_atom_dimensions(self) -> None:
        if not self.contains("topology.spgt.h5", "/topology/atom_count"):
            return
        atom_count = int(self.read_scalar("topology.spgt.h5", "/topology/atom_count"))
        for bundle_file in ("restart.spgr.h5", "trajectory.spg.h5md"):
            particle = "all" if bundle_file == "restart.spgr.h5" else self.case.particle_stream
            for path in (
                f"/particles/{particle}/position/value",
                f"/particles/{particle}/velocity/value",
            ):
                if not self.contains(bundle_file, path):
                    continue
                shape = self._handles[bundle_file][path].shape
                if len(shape) != 3 or shape[1:] != (atom_count, 3):
                    raise BundleValidationError(
                        f"{bundle_file} {path} has shape {shape}, expected (*, {atom_count}, 3)"
                    )
        self._validate_trajectory_frames()

    def _validate_trajectory_frames(self) -> None:
        if not self.has_bundle_file("trajectory.spg.h5md"):
            return
        root = f"/particles/{self.case.particle_stream}"
        handle = self._handles["trajectory.spg.h5md"]
        if root not in handle:
            self._validation_problem(
                f"trajectory particle stream {self.case.particle_stream!r} does not exist"
            )
            return
        frame_counts = {}
        for component in ("position", "velocity", "box/edges"):
            path = f"{root}/{component}/value"
            if path in handle:
                frame_counts[component] = int(handle[path].shape[0])
        if frame_counts and len(set(frame_counts.values())) != 1:
            raise BundleValidationError(
                f"trajectory particle stream {self.case.particle_stream!r} has mismatched frame counts: "
                f"{frame_counts}"
            )
        expected = next(iter(frame_counts.values()), 0)
        for axis in ("step", "time"):
            path = f"{root}/{axis}"
            if path in handle and handle[path].shape[0] != expected:
                raise BundleValidationError(
                    f"trajectory {path} has {handle[path].shape[0]} frames, expected {expected}"
                )

    def _validate_restart_load_policy(self) -> None:
        if not self.has_bundle_file("restart.spgr.h5"):
            return
        policy = (self.case.restart_load or "").strip().lower()
        if policy not in {"structural", "dynamic", "protocol", "full"}:
            self._validation_problem(
                f"input_h5_restart_load must be structural, dynamic, protocol, or full; got {policy!r}"
            )
            return

        has_structural = self.contains(
            "restart.spgr.h5", "/particles/all/position/value"
        )
        has_dynamic = self.contains(
            "restart.spgr.h5", "/parameters/restart/thermostat/nose_hoover_chain"
        )
        handle = self._handles["restart.spgr.h5"]
        sidecar_root = "/parameters/restart/protocol_sidecars"
        has_protocol = sidecar_root in handle and len(handle[sidecar_root]) > 0
        has_protocol = has_protocol or any(
            self.contains("restart.spgr.h5", path)
            for path in (
                "/parameters/restart/bias/sits/SITS/nk",
                "/parameters/restart/bias/meta/default/potential/value",
                "/parameters/restart/bias/meta/default/scatter/value",
                "/parameters/restart/bias/meta/default/hills_typed/value",
            )
        )
        available = {
            "structural": has_structural,
            "dynamic": has_dynamic,
            "protocol": has_protocol,
        }
        if policy == "full":
            if sum(available.values()) < 2:
                self._validation_problem(
                    "input_h5_restart_load='full' requires at least two restart state components"
                )
        elif not available[policy]:
            self._validation_problem(
                f"input_h5_restart_load={policy!r} has no corresponding restart state"
            )

    def _validation_problem(self, message: str) -> None:
        if self.strict:
            raise BundleValidationError(message)
        self.warnings.append(message)

    def _require_handle(self, bundle_file: str):
        handle = self._handles.get(bundle_file)
        if handle is None:
            raise BundleValidationError(f"bundle artifact is not configured: {bundle_file}")
        return handle


def _decode_h5_text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    if isinstance(value, np.ndarray) and value.shape == ():
        return _decode_h5_text(value.item())
    return str(value)
