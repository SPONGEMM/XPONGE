"""
Typed parsers for legacy SPONGE output files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .trajectory_parsers import box_to_edges


@dataclass(frozen=True)
class TypedDataset:
    """A typed output dataset ready to be written into an H5MD bundle."""

    path: str
    data: np.ndarray | str


def parse_legacy_output_file(
    key: str,
    path: Path,
    *,
    atom_count: int | None = None,
    sits_nk_width: int | None = None,
) -> list[TypedDataset] | None:
    """Parse a supported legacy output file into H5MD/SPONGE extension datasets."""

    if key == "mdinfo":
        return [TypedDataset("/parameters/sponge/log/mdinfo_text", path.read_text(encoding="utf-8"))]
    if key == "mdout":
        return _parse_mdout(path)
    if key in {"crd", "vel", "frc"}:
        if atom_count is None:
            raise ValueError(f"atom_count is required to parse legacy {key} trajectory output")
        value_path = {
            "crd": "/particles/all/position/value",
            "vel": "/particles/all/velocity/value",
            "frc": "/particles/all/force/value",
        }[key]
        return _parse_vector_trajectory(path, atom_count=atom_count, value_path=value_path)
    if key == "box":
        return _parse_box_trajectory(path)
    if key == "nose_hoover_chain_crd":
        return _parse_observable_matrix(
            path,
            root="/observables/all/thermostat/nose_hoover_chain/coordinate",
        )
    if key == "nose_hoover_chain_vel":
        return _parse_observable_matrix(
            path,
            root="/observables/all/thermostat/nose_hoover_chain/velocity",
        )
    if key == "SITS_nk_traj_file":
        if sits_nk_width is None:
            raise ValueError("sits_nk_width is required to parse SITS_nk_traj_file")
        return _parse_binary_observable_matrix(
            path,
            root="/observables/all/sits/SITS/nk",
            width=sits_nk_width,
            dtype=np.float32,
        )
    if key == "SITS_nk_rest_file":
        return _parse_sits_nk_restart_output(path)
    if key == "meta_potential_out_file":
        return _parse_meta_potential_export(path)
    if key in {
        "meta_direct_export",
        "meta_hills_log",
        "meta_history_log",
        "meta_edge_log",
    }:
        component = {
            "meta_direct_export": "direct_export",
            "meta_hills_log": "hills",
            "meta_history_log": "history",
            "meta_edge_log": "edge",
        }[key]
        return _parse_metadynamics_diagnostic_text(path, component=component)
    if key == "qc_scf_output":
        return _parse_qc_scf_output(path)
    if key == "reaxff_eeq_charges":
        return _parse_reaxff_eeq_charges(path)
    if key == "rst":
        return _parse_restart_rst7_output(path)
    if key == "nose_hoover_chain_restart_output":
        return _parse_nose_hoover_chain_restart_output(path)
    return None


def _parse_mdout(path: Path) -> list[TypedDataset]:
    raw_text = path.read_text(encoding="utf-8")
    rows = [_split_mdout_line(line) for line in raw_text.splitlines()]
    rows = [row for row in rows if row]
    if not rows:
        raise ValueError(f"{path} is empty")

    header = None
    data_rows = rows
    if any(_is_not_number(token) for token in rows[0]):
        header = [_safe_column_name(token) for token in rows[0]]
        data_rows = rows[1:]
    if not data_rows:
        raise ValueError(f"{path} does not contain mdout numeric rows")

    width = len(data_rows[0])
    if any(len(row) != width for row in data_rows):
        raise ValueError(f"{path} has inconsistent mdout row widths")
    values = np.asarray([[float(token) for token in row] for row in data_rows], dtype=np.float64)
    if header is None:
        header = [f"col_{index}" for index in range(width)]
    if len(header) != width:
        raise ValueError(f"{path} mdout header has {len(header)} columns but data has {width}")

    lower_header = [name.lower() for name in header]
    step_index = lower_header.index("step") if "step" in lower_header else None
    time_index = lower_header.index("time") if "time" in lower_header else None
    step = (
        values[:, step_index].astype(np.int64)
        if step_index is not None
        else np.arange(values.shape[0], dtype=np.int64)
    )
    time = (
        values[:, time_index].astype(np.float64)
        if time_index is not None
        else np.arange(values.shape[0], dtype=np.float64)
    )

    datasets: list[TypedDataset] = [
        TypedDataset("/observables/all/step", step),
        TypedDataset("/observables/all/time", time),
        TypedDataset("/parameters/sponge/mdout/raw_text", raw_text),
    ]
    original_names = []
    hdf5_names = []
    categories = []
    used_names: set[str] = set()
    for column_index, original_name in enumerate(header):
        if column_index in {step_index, time_index}:
            continue
        hdf5_name = _dedupe_observable_path(_mdout_h5_relative_path(original_name), used_names)
        original_names.append(original_name)
        hdf5_names.append(hdf5_name)
        categories.append(_mdout_category(original_name))
        datasets.append(TypedDataset(f"/observables/all/{hdf5_name}/value", values[:, column_index]))
    datasets.extend(
        [
            TypedDataset("/parameters/sponge/mdout/columns/original_name", np.asarray(original_names, dtype=object)),
            TypedDataset("/parameters/sponge/mdout/columns/hdf5_name", np.asarray(hdf5_names, dtype=object)),
            TypedDataset("/parameters/sponge/mdout/columns/category", np.asarray(categories, dtype=object)),
        ]
    )
    return datasets


def _parse_vector_trajectory(path: Path, *, atom_count: int, value_path: str) -> list[TypedDataset]:
    if atom_count <= 0:
        raise ValueError("atom_count must be positive")
    values = np.frombuffer(path.read_bytes(), dtype=np.float32)
    frame_width = atom_count * 3
    if values.size == 0 or values.size % frame_width != 0:
        raise ValueError(f"{path} has {values.size} float32 values, not a multiple of {frame_width}")
    trajectory = values.reshape(values.size // frame_width, atom_count, 3)
    frame_count = trajectory.shape[0]
    return [
        TypedDataset(value_path, trajectory),
        TypedDataset("/particles/all/step", np.arange(frame_count, dtype=np.int64)),
        TypedDataset("/particles/all/time", np.arange(frame_count, dtype=np.float64)),
    ]


def _parse_box_trajectory(path: Path) -> list[TypedDataset]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    edges = []
    for row_index, row in enumerate(rows):
        if len(row) != 6:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected 6")
        edges.append(box_to_edges(np.asarray([float(value) for value in row], dtype=np.float32)))
    values = np.asarray(edges, dtype=np.float32)
    frame_count = values.shape[0]
    return [
        TypedDataset("/particles/all/box/edges/value", values),
        TypedDataset("/particles/all/step", np.arange(frame_count, dtype=np.int64)),
        TypedDataset("/particles/all/time", np.arange(frame_count, dtype=np.float64)),
    ]


def _parse_observable_matrix(path: Path, *, root: str) -> list[TypedDataset]:
    rows = [[float(token) for token in line.split()] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError(f"{path} has inconsistent row widths")
    values = np.asarray(rows, dtype=np.float64)
    frame_count = values.shape[0]
    return [
        TypedDataset(f"{root}/value", values),
        TypedDataset(f"{root}/step", np.arange(frame_count, dtype=np.int64)),
        TypedDataset(f"{root}/time", np.arange(frame_count, dtype=np.float64)),
    ]


def _parse_binary_observable_matrix(
    path: Path,
    *,
    root: str,
    width: int,
    dtype,
) -> list[TypedDataset]:
    if width <= 0:
        raise ValueError("observable matrix width must be positive")
    values = np.frombuffer(path.read_bytes(), dtype=dtype)
    if values.size == 0 or values.size % width != 0:
        raise ValueError(f"{path} has {values.size} values, not a multiple of {width}")
    matrix = values.reshape(values.size // width, width)
    frame_count = matrix.shape[0]
    return [
        TypedDataset(f"{root}/value", matrix),
        TypedDataset(f"{root}/step", np.arange(frame_count, dtype=np.int64)),
        TypedDataset(f"{root}/time", np.arange(frame_count, dtype=np.float64)),
    ]


def _parse_sits_nk_restart_output(path: Path) -> list[TypedDataset]:
    values = np.asarray(
        [float(token) for line in path.read_text(encoding="utf-8").splitlines() for token in line.split()],
        dtype=np.float32,
    )
    if values.size == 0:
        raise ValueError(f"{path} is empty")
    return [
        TypedDataset("/parameters/sponge/restart_exports/sits/SITS/nk/value", values),
        TypedDataset("/parameters/sponge/restart_exports/sits/SITS/nk/raw_text", path.read_text(encoding="utf-8")),
    ]


def _parse_meta_potential_export(path: Path) -> list[TypedDataset]:
    raw_text = path.read_text(encoding="utf-8")
    lines = [line.strip() for line in raw_text.splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"{path} is too short for a meta potential export")
    title = lines[0]
    ndim = _infer_meta_export_ndim(lines, path)
    axes = []
    for axis_index in range(ndim):
        row = [float(value) for value in lines[1 + axis_index].split()]
        if len(row) != 3:
            raise ValueError(f"{path} meta axis row {axis_index} must have 3 columns")
        axes.append(row)
    header_values = lines[1 + ndim].split()
    if len(header_values) < ndim + 1:
        raise ValueError(f"{path} meta grid row must contain ndim grid counts and row count")
    grid = np.asarray([int(value) for value in header_values[:ndim]], dtype=np.int32)
    row_count = int(header_values[ndim])
    data_lines = lines[2 + ndim :]
    if len(data_lines) != row_count:
        raise ValueError(f"{path} declares {row_count} meta rows but contains {len(data_lines)}")
    rows = np.asarray([[float(value) for value in line.split()] for line in data_lines], dtype=np.float32)
    if rows.ndim != 2 or rows.shape[1] < ndim + 1:
        raise ValueError(f"{path} meta rows must contain coordinates and at least one value")

    root = "/parameters/sponge/metadynamics/default/potential_export"
    axes_array = np.asarray(axes, dtype=np.float32)
    datasets = [
        TypedDataset(f"{root}/raw_text", raw_text),
        TypedDataset(f"{root}/title", title),
        TypedDataset(f"{root}/grid/min", axes_array[:, 0]),
        TypedDataset(f"{root}/grid/max", axes_array[:, 1]),
        TypedDataset(f"{root}/grid/width", axes_array[:, 2]),
        TypedDataset(f"{root}/grid/count", grid),
        TypedDataset(f"{root}/coordinate", rows[:, :ndim]),
        TypedDataset(f"{root}/value", rows[:, -1]),
    ]
    if rows.shape[1] >= ndim + 2:
        datasets.append(TypedDataset(f"{root}/force", rows[:, ndim:-1]))
    return datasets


def _parse_metadynamics_diagnostic_text(path: Path, *, component: str, name: str = "default") -> list[TypedDataset]:
    return [
        TypedDataset(
            f"/parameters/sponge/metadynamics/{name}/{component}",
            path.read_text(encoding="utf-8"),
        )
    ]


def _infer_meta_export_ndim(lines: list[str], path: Path) -> int:
    for ndim in range(1, min(6, len(lines) - 1)):
        header_index = 1 + ndim
        header_values = lines[header_index].split()
        if len(header_values) < ndim + 1:
            continue
        try:
            row_count = int(header_values[ndim])
            for value in header_values[:ndim]:
                int(value)
        except ValueError:
            continue
        if row_count == len(lines[header_index + 1 :]):
            return ndim
    raise ValueError(f"{path} cannot infer meta potential export dimensionality")


def _parse_qc_scf_output(path: Path) -> list[TypedDataset]:
    return [TypedDataset("/parameters/sponge/qc/scf_output", path.read_text(encoding="utf-8"))]


def _parse_reaxff_eeq_charges(path: Path) -> list[TypedDataset]:
    raw_text = path.read_text(encoding="utf-8")
    rows = []
    for line_index, line in enumerate(raw_text.splitlines()):
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if len(tokens) != 2:
            raise ValueError(f"{path} line {line_index + 1} must contain atom index and charge")
        rows.append((int(tokens[0]), float(tokens[1])))
    if not rows:
        raise ValueError(f"{path} is empty")
    return [
        TypedDataset(
            "/parameters/sponge/reaxff/eeq_charges/atom_index",
            np.asarray([row[0] for row in rows], dtype=np.int32),
        ),
        TypedDataset(
            "/parameters/sponge/reaxff/eeq_charges/value",
            np.asarray([row[1] for row in rows], dtype=np.float32),
        ),
        TypedDataset("/parameters/sponge/reaxff/eeq_charges/raw_text", raw_text),
    ]


def _parse_restart_rst7_output(path: Path) -> list[TypedDataset]:
    title, atom_count, time, coordinates, velocities, box = _read_amber_rst7(path)
    datasets = [
        TypedDataset("/particles/all/step", np.asarray([0], dtype=np.int64)),
        TypedDataset("/particles/all/time", np.asarray([time], dtype=np.float64)),
        TypedDataset("/particles/all/position/value", coordinates[np.newaxis, :, :]),
        TypedDataset("/particles/all/box/edges/value", box_to_edges(box)[np.newaxis, :, :]),
        TypedDataset("/parameters/restart/source/title", title),
        TypedDataset("/parameters/restart/source/atom_count", np.asarray(atom_count, dtype=np.int64)),
    ]
    if velocities is not None:
        datasets.append(TypedDataset("/particles/all/velocity/value", velocities[np.newaxis, :, :]))
    return datasets


def _parse_nose_hoover_chain_restart_output(path: Path) -> list[TypedDataset]:
    matrix = _read_numeric_matrix(path)
    if matrix.shape[1] != 2:
        raise ValueError(f"{path} nose-hoover chain restart rows must have 2 columns")
    return [
        TypedDataset(
            "/parameters/restart/thermostat/nose_hoover_chain",
            matrix.astype(np.float32),
        )
    ]


def _read_amber_rst7(path: Path) -> tuple[str, int, float, np.ndarray, np.ndarray | None, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path} is too short for Amber rst7")
    title = lines[0].strip()
    header = lines[1].split()
    if not header:
        raise ValueError(f"{path} missing Amber rst7 atom count")
    atom_count = int(header[0])
    time = float(header[1]) if len(header) >= 2 else 0.0
    has_velocity = len(header) >= 2
    values = [float(token) for line in lines[2:] for token in line.split()]
    coord_count = atom_count * 3
    offset = 0
    if len(values) < coord_count:
        raise ValueError(f"{path} missing Amber rst7 coordinates")
    coordinates = np.asarray(values[offset : offset + coord_count], dtype=np.float32).reshape(atom_count, 3)
    offset += coord_count
    velocities = None
    if has_velocity:
        if len(values) < offset + coord_count:
            raise ValueError(f"{path} missing Amber rst7 velocities")
        velocities = np.asarray(values[offset : offset + coord_count], dtype=np.float32).reshape(atom_count, 3)
        offset += coord_count
    if len(values) < offset + 6:
        raise ValueError(f"{path} missing Amber rst7 box")
    box = np.asarray(values[offset : offset + 6], dtype=np.float32)
    return title, atom_count, time, coordinates, velocities, box


def _read_numeric_matrix(path: Path) -> np.ndarray:
    rows = [[float(token) for token in line.split()] for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError(f"{path} has inconsistent row widths")
    return np.asarray(rows, dtype=np.float64)


def _split_mdout_line(line: str) -> list[str]:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return []
    return stripped.split()


def _is_not_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return True
    return False


def _safe_column_name(name: str) -> str:
    return name.strip().strip(":,")


def _safe_h5_name(name: str) -> str:
    safe = "".join(char if char.isalnum() or char == "_" else "_" for char in name.strip())
    if not safe or safe[0].isdigit():
        safe = "_" + safe
    if safe in {"step", "time"}:
        safe = f"mdout_{safe}"
    return safe


def _dedupe_name(name: str, used_names: set[str]) -> str:
    candidate = name
    suffix = 0
    while candidate in used_names:
        suffix += 1
        candidate = f"{name}_{suffix}"
    used_names.add(candidate)
    return candidate


def _dedupe_observable_path(path: str, used_paths: set[str]) -> str:
    candidate = path
    suffix = 0
    while candidate in used_paths:
        suffix += 1
        if "/" in path:
            prefix, leaf = path.rsplit("/", 1)
            candidate = f"{prefix}/{leaf}_{suffix}"
        else:
            candidate = f"{path}_{suffix}"
    used_paths.add(candidate)
    return candidate


def _mdout_h5_relative_path(name: str) -> str:
    lower = name.lower()
    safe = _safe_h5_name(name)
    if lower in {"meta", "rbias", "rct"}:
        return f"metadynamics/{lower}"
    if lower == "qc":
        return "qc/energy"
    if lower == "qc_s_sq":
        return "qc/spin_square"
    if _is_reaxff_output_key(name):
        return f"reaxff/{safe}"
    return safe


def _is_reaxff_output_key(name: str) -> bool:
    return name == "REAXFF" or name.startswith("REAXFF_")


def _mdout_category(name: str) -> str:
    lower = name.lower()
    if lower in {"temperature", "potential", "kinetic", "total", "density", "volume"}:
        return "core"
    if lower.startswith("p") and lower in {"pressure", "pxx", "pyy", "pzz", "pxy", "pxz", "pyz"}:
        return "pressure"
    if lower in {"lj", "lj_short", "lj_long", "coulomb", "pm", "gb"}:
        return "nonbonded"
    if lower in {"bond", "angle", "dihedral", "improper_dihedral", "cmap", "nb14_lj", "nb14_ee"}:
        return "bonded"
    if lower in {"restrain", "steer_cv", "restrain_cv"}:
        return "restraint"
    if lower in {"meta", "rbias", "rct"}:
        return "metadynamics"
    if _is_reaxff_output_key(name):
        return "reaxff"
    if lower.startswith("qc"):
        return "qc"
    return "custom"
