"""Typed restart and rerun particle data exporters."""

from __future__ import annotations

import numpy as np

from .errors import BundleExportError
from .legacy_materializer import LegacyPayload
from .topology_exporters import _format_float


def export_coordinate(contract, reader, context) -> list[LegacyPayload]:
    position = _last_frame(
        reader.read(contract.bundle_file, "/particles/all/position/value"),
        component="position",
    )
    edges = _last_frame(
        reader.read(contract.bundle_file, "/particles/all/box/edges/value"),
        component="box edges",
    )
    if position.ndim != 2 or position.shape[1] != 3:
        raise BundleExportError(f"restart position has shape {position.shape}, expected (atoms, 3)")
    if edges.shape != (3, 3):
        raise BundleExportError(f"restart box edges have shape {edges.shape}, expected (3, 3)")
    box = edges_to_box(edges)
    lines = [str(position.shape[0])]
    lines.extend(" ".join(_format_float(value) for value in row) for row in position)
    lines.append(" ".join(_format_float(value) for value in box))
    return [LegacyPayload("coordinate_in_file", "\n".join(lines) + "\n")]


def export_velocity(contract, reader, context) -> list[LegacyPayload]:
    velocity = _last_frame(
        reader.read(contract.bundle_file, "/particles/all/velocity/value"),
        component="velocity",
    )
    if velocity.ndim != 2 or velocity.shape[1] != 3:
        raise BundleExportError(f"restart velocity has shape {velocity.shape}, expected (atoms, 3)")
    lines = [str(velocity.shape[0])]
    lines.extend(" ".join(_format_float(value) for value in row) for row in velocity)
    return [LegacyPayload("velocity_in_file", "\n".join(lines) + "\n")]


def export_rerun_position(contract, reader, context) -> list[LegacyPayload]:
    root = f"/particles/{context.particle_stream}"
    values = _particle_trajectory(
        reader.read(contract.bundle_file, root + "/position/value"),
        component="rerun position",
    )
    return [
        LegacyPayload(
            "crd",
            np.asarray(values, dtype=np.float32).tobytes(order="C"),
            filename=f"{context.prefix}_crd.dat",
            binary=True,
        )
    ]


def export_rerun_velocity(contract, reader, context) -> list[LegacyPayload]:
    root = f"/particles/{context.particle_stream}"
    values = _particle_trajectory(
        reader.read(contract.bundle_file, root + "/velocity/value"),
        component="rerun velocity",
    )
    return [
        LegacyPayload(
            "vel",
            np.asarray(values, dtype=np.float32).tobytes(order="C"),
            filename=f"{context.prefix}_vel.dat",
            binary=True,
        )
    ]


def export_rerun_box(contract, reader, context) -> list[LegacyPayload]:
    root = f"/particles/{context.particle_stream}"
    values = np.asarray(reader.read(contract.bundle_file, root + "/box/edges/value"))
    if values.ndim == 2:
        values = values[np.newaxis, ...]
    if values.ndim != 3 or values.shape[1:] != (3, 3):
        raise BundleExportError(f"rerun box edges have unsupported shape {values.shape}")
    lines = [" ".join(_format_float(value) for value in edges_to_box(frame)) for frame in values]
    return [
        LegacyPayload(
            "box",
            "\n".join(lines) + "\n",
            filename=f"{context.prefix}_box.txt",
        )
    ]


def edges_to_box(edges) -> np.ndarray:
    """Convert a 3x3 H5MD cell matrix to lengths and angles in degrees."""

    vectors = np.asarray(edges, dtype=np.float64)
    if vectors.shape != (3, 3):
        raise BundleExportError(f"box edges have shape {vectors.shape}, expected (3, 3)")
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 0) or not np.all(np.isfinite(lengths)):
        raise BundleExportError("box edges contain a zero or non-finite vector")

    def angle(left: int, right: int) -> float:
        cosine = np.dot(vectors[left], vectors[right]) / (lengths[left] * lengths[right])
        return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))

    alpha = angle(1, 2)
    beta = angle(0, 2)
    gamma = angle(0, 1)
    return np.asarray([*lengths, alpha, beta, gamma], dtype=np.float64)


TRAJECTORY_EXPORTERS = {
    "coordinate_in_file": export_coordinate,
    "velocity_in_file": export_velocity,
    "crd": export_rerun_position,
    "vel": export_rerun_velocity,
    "box": export_rerun_box,
}


def _last_frame(values, *, component: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim == 3:
        if array.shape[0] == 0:
            raise BundleExportError(f"{component} has no frames")
        return array[-1]
    if array.ndim == 2:
        return array
    raise BundleExportError(f"{component} has unsupported shape {array.shape}")


def _particle_trajectory(values, *, component: str) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim != 3 or array.shape[0] == 0 or array.shape[2] != 3:
        raise BundleExportError(f"{component} has unsupported shape {array.shape}")
    return array
