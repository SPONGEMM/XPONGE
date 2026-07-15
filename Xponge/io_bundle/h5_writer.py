"""
Small HDF5 writing helpers for SPONGE bundled fixture generation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


class MissingH5PyError(RuntimeError):
    """Raised when HDF5 output is requested without h5py installed."""


def _import_h5py() -> Any:
    try:
        import h5py  # type: ignore
    except ImportError as exc:
        raise MissingH5PyError(
            "legacy-to-bundle requires h5py to write .h5/.h5md files; install h5py and rerun"
        ) from exc
    return h5py


def ensure_dataset(file_path: Path, dataset_path: str, data: Any) -> None:
    """Create or replace a dataset in an HDF5 file."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    data, dtype = _normalize_dataset_payload(h5py, data)
    with h5py.File(file_path, "a") as handle:
        if dataset_path in handle:
            del handle[dataset_path]
        parent = str(Path(dataset_path).parent).replace("\\", "/")
        if parent != ".":
            handle.require_group(parent)
        if dtype is None:
            handle.create_dataset(dataset_path, data=data)
        else:
            handle.create_dataset(dataset_path, data=data, dtype=dtype)


def ensure_group(file_path: Path, group_path: str) -> None:
    """Create an HDF5 group if it does not already exist."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a") as handle:
        handle.require_group(group_path)


def ensure_hard_link(file_path: Path, target_path: str, link_path: str) -> None:
    """Create or replace a hard link to an existing HDF5 object."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a") as handle:
        if target_path not in handle:
            raise KeyError(f"hard-link target does not exist: {target_path}")
        parent = str(Path(link_path).parent).replace("\\", "/")
        if parent != ".":
            handle.require_group(parent)
        if link_path in handle:
            del handle[link_path]
        handle[link_path] = handle[target_path]


def set_attrs(file_path: Path, object_path: str, attrs: dict[str, Any]) -> None:
    """Set HDF5 attributes on an existing group or dataset."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a") as handle:
        if object_path not in handle:
            raise KeyError(f"HDF5 object does not exist: {object_path}")
        for key, value in attrs.items():
            handle[object_path].attrs[key] = value


def write_string(file_path: Path, dataset_path: str, value: str) -> None:
    """Create or replace a UTF-8 string dataset."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a") as handle:
        if dataset_path in handle:
            del handle[dataset_path]
        parent = str(Path(dataset_path).parent).replace("\\", "/")
        if parent != ".":
            handle.require_group(parent)
        dtype = h5py.string_dtype(encoding="utf-8")
        handle.create_dataset(dataset_path, data=value, dtype=dtype)


def write_bytes(file_path: Path, dataset_path: str, value: bytes) -> None:
    """Create or replace a byte-vector dataset."""

    ensure_dataset(file_path, dataset_path, np.frombuffer(value, dtype=np.uint8))


def write_string_array(file_path: Path, dataset_path: str, values: list[str]) -> None:
    """Create or replace a UTF-8 string-array dataset."""

    h5py = _import_h5py()
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(file_path, "a") as handle:
        if dataset_path in handle:
            del handle[dataset_path]
        parent = str(Path(dataset_path).parent).replace("\\", "/")
        if parent != ".":
            handle.require_group(parent)
        dtype = h5py.string_dtype(encoding="utf-8")
        handle.create_dataset(dataset_path, data=np.asarray(values, dtype=object), dtype=dtype)


def _normalize_dataset_payload(h5py, data: Any):
    if isinstance(data, str):
        return data, h5py.string_dtype(encoding="utf-8")
    array = np.asarray(data)
    if array.dtype.kind in {"U", "S"}:
        return array.astype(object), h5py.string_dtype(encoding="utf-8")
    if array.dtype.kind == "O" and all(isinstance(value, str) for value in array.flat):
        return array, h5py.string_dtype(encoding="utf-8")
    return data, None
