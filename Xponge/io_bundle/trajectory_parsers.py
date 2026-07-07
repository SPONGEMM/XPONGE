"""
Typed parsers for legacy SPONGE rerun trajectory inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class TypedDataset:
    """A typed trajectory dataset ready to be written into an H5MD bundle."""

    path: str
    data: np.ndarray


def parse_trajectory_file(key: str, path: Path, *, atom_count: int) -> list[TypedDataset] | None:
    """Parse legacy rerun trajectory files into H5MD particle datasets."""

    if key == "crd":
        position = _read_binary_vector_trajectory(path, atom_count=atom_count)
        frame_count = position.shape[0]
        return [
            TypedDataset("/particles/all/position/value", position),
            TypedDataset("/particles/all/step", np.arange(frame_count, dtype=np.int64)),
            TypedDataset("/particles/all/time", np.arange(frame_count, dtype=np.float64)),
        ]
    if key == "box":
        return [
            TypedDataset(
                "/particles/all/box/edges/value",
                _read_box_trajectory(path),
            )
        ]
    if key == "vel":
        return [
            TypedDataset(
                "/particles/all/velocity/value",
                _read_binary_vector_trajectory(path, atom_count=atom_count),
            )
        ]
    return None


def _read_binary_vector_trajectory(path: Path, *, atom_count: int) -> np.ndarray:
    if atom_count <= 0:
        raise ValueError("atom_count must be positive")
    values = np.frombuffer(path.read_bytes(), dtype=np.float32)
    frame_width = atom_count * 3
    if values.size == 0 or values.size % frame_width != 0:
        raise ValueError(f"{path} has {values.size} float32 values, not a multiple of {frame_width}")
    return values.reshape(values.size // frame_width, atom_count, 3)


def _read_box_trajectory(path: Path) -> np.ndarray:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    edges = []
    for row_index, row in enumerate(rows):
        if len(row) != 6:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected 6")
        edges.append(box_to_edges(np.asarray([float(value) for value in row], dtype=np.float32)))
    return np.asarray(edges, dtype=np.float32)


def box_to_edges(box: np.ndarray) -> np.ndarray:
    """Convert SPONGE length/angle box rows into H5MD edge matrices."""

    if box.shape[0] < 3:
        raise ValueError("box must contain at least 3 length values")
    lx, ly, lz = [float(value) for value in box[:3]]
    if box.shape[0] < 6:
        return np.asarray([[lx, 0.0, 0.0], [0.0, ly, 0.0], [0.0, 0.0, lz]], dtype=np.float32)

    alpha, beta, gamma = np.deg2rad([float(value) for value in box[3:6]])
    sin_gamma = float(np.sin(gamma))
    if abs(sin_gamma) < 1e-7:
        raise ValueError("box gamma angle produces a singular cell")
    cos_alpha = float(np.cos(alpha))
    cos_beta = float(np.cos(beta))
    cos_gamma = float(np.cos(gamma))
    cy = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cz_squared = lz * lz - (lz * cos_beta) ** 2 - cy * cy
    if cz_squared < -1e-5:
        raise ValueError("box angles produce an invalid cell")
    return np.asarray(
        [
            [lx, 0.0, 0.0],
            [ly * cos_gamma, ly * sin_gamma, 0.0],
            [lz * cos_beta, cy, np.sqrt(max(cz_squared, 0.0))],
        ],
        dtype=np.float32,
    )
