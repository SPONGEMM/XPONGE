"""
Typed parsers for legacy SPONGE protocol and restart state files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 fallback path.
    tomllib = None

import numpy as np


@dataclass(frozen=True)
class TypedDataset:
    """A typed dataset ready to be written into a bundle."""

    path: str
    data: np.ndarray | str


def parse_protocol_or_restart_file(key: str, path: Path) -> list[TypedDataset] | None:
    """Parse simple protocol/restart legacy files into typed HDF5 datasets."""

    if key == "cv_in_file":
        return _parse_cv_config(path, root="/cv")
    if key == "restrain_in_file":
        return _parse_cv_config(path, root="/restraint")
    if key == "restrain_cv_in_file":
        return _parse_cv_config(path, root="/restraint/cv")
    if key == "steer_cv_in_file":
        return _parse_cv_config(path, root="/steer")
    if key == "SITS_in_file":
        return _parse_cv_config(path, root="/sits")
    if key == "SITS_nk_in_file":
        return [
            TypedDataset(
                "/parameters/restart/bias/sits/SITS/nk",
                _read_flat_float_file(path),
            )
        ]
    if key == "restrain_atom_id":
        return [
            TypedDataset(
                "/restraint/default/atom_indices",
                _read_flat_int_file(path),
            )
        ]
    if key == "restrain_weight_in_file":
        return [
            TypedDataset(
                "/restraint/default/weight",
                _read_matrix_file(path, dtype=np.float32, columns=3),
            )
        ]
    if key == "restrain_coordinate_in_file":
        return [
            TypedDataset(
                "/parameters/restart/references/restraint/default/coordinate",
                _read_counted_xyz(path),
            )
        ]
    if key == "restrain_amber_rst7":
        return _parse_restrain_amber_rst7(path)
    if key == "nose_hoover_chain_restart_input":
        return [
            TypedDataset(
                "/parameters/restart/thermostat/nose_hoover_chain",
                _read_matrix_file(path, dtype=np.float32, columns=2),
            )
        ]
    if key == "constrain_in_file":
        atoms, r0 = _read_counted_pairs_with_r0(path)
        return [
            TypedDataset("/constraint/default/pairs/atoms", atoms),
            TypedDataset("/constraint/default/pairs/r0", r0),
        ]
    if key == "SITS_atom_in_file":
        return [
            TypedDataset(
                "/sits/atom_indices",
                _read_flat_int_file(path),
            )
        ]
    if key == "soft_walls_in_file":
        names, potentials = _read_soft_wall_config(path)
        return [
            TypedDataset("/wall/soft/name", np.asarray(names, dtype=object)),
            TypedDataset("/wall/soft/potential", np.asarray(potentials, dtype=object)),
            TypedDataset("/wall/soft/count", np.asarray(len(names), dtype=np.int64)),
        ]
    if key == "meta_edge_in_file":
        return _parse_meta_edge(path)
    if key in {"meta_potential_in_file", "meta_scatter_in_file"}:
        root = "/parameters/restart/bias/meta/default/potential"
        if key == "meta_scatter_in_file":
            root = "/parameters/restart/bias/meta/default/scatter"
        datasets = _parse_meta_potential(path, root=root)
        if key == "meta_potential_in_file":
            datasets.insert(
                0,
                TypedDataset(
                    "/parameters/restart/bias/meta/default/potential_export",
                    path.read_text(encoding="utf-8"),
                ),
            )
        return datasets
    if key in {"hill_in_file", "hills_in_file", "metad_hills_in_file"}:
        return _parse_meta_hills(path)
    return None


def _parse_cv_config(path: Path, *, root: str) -> list[TypedDataset]:
    if path.suffix.lower() == ".toml":
        sections = _read_toml_sections(path)
    else:
        sections = _read_braced_sections(path)
    if not sections:
        raise ValueError(f"{path} does not contain any sections")

    section_names = []
    key_values: list[tuple[str, str]] = []
    key_offsets = [0]
    for section_name, values in sections:
        section_names.append(section_name)
        for key, value in values:
            key_values.append((key, value))
        key_offsets.append(len(key_values))

    keys = [key for key, _ in key_values]
    values = [value for _, value in key_values]
    return [
        TypedDataset(f"{root}/config/section/name", np.asarray(section_names, dtype=object)),
        TypedDataset(f"{root}/config/section/key_offset", np.asarray(key_offsets, dtype=np.int64)),
        TypedDataset(f"{root}/config/key", np.asarray(keys, dtype=object)),
        TypedDataset(f"{root}/config/value", np.asarray(values, dtype=object)),
        TypedDataset(f"{root}/config/section/count", np.asarray(len(section_names), dtype=np.int64)),
    ]


def _read_toml_sections(path: Path) -> list[tuple[str, list[tuple[str, str]]]]:
    if tomllib is None:
        raise ValueError(f"{path} is TOML but tomllib is unavailable on this Python")
    data = tomllib.loads(path.read_text(encoding="utf-8"))
    sections: list[tuple[str, list[tuple[str, str]]]] = []
    for section_name, value in data.items():
        if isinstance(value, dict):
            flat_values = _flatten_toml_values(value)
            sections.append((section_name, flat_values))
        else:
            sections.append(("", [(section_name, _format_config_value(value))]))
    return sections


def _flatten_toml_values(data: dict, prefix: str = "") -> list[tuple[str, str]]:
    values = []
    for key, value in data.items():
        name = f"{prefix}_{key}" if prefix else key
        if isinstance(value, dict):
            values.extend(_flatten_toml_values(value, name))
        else:
            values.append((name, _format_config_value(value)))
    return values


def _format_config_value(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return " ".join(_format_config_value(item) for item in value)
    return str(value)


def _read_braced_sections(path: Path) -> list[tuple[str, list[tuple[str, str]]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    sections: list[tuple[str, list[tuple[str, str]]]] = []
    index = 0
    while index < len(lines):
        line = _strip_cv_comment(lines[index]).strip()
        index += 1
        if not line:
            continue
        section_name = line
        if index >= len(lines) or _strip_cv_comment(lines[index]).strip() != "{":
            raise ValueError(f"{path} expected '{{' after CV section {section_name}")
        index += 1
        values: list[tuple[str, str]] = []
        while index < len(lines):
            raw_line = _strip_cv_comment(lines[index]).strip()
            index += 1
            if not raw_line:
                continue
            if raw_line == "}":
                break
            if "=" not in raw_line:
                raise ValueError(f"{path} expected key = value in CV section {section_name}")
            key, value = raw_line.split("=", 1)
            values.append((key.strip(), value.strip()))
        else:
            raise ValueError(f"{path} missing '}}' for CV section {section_name}")
        sections.append((section_name, values))
    return sections


def _strip_cv_comment(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return ""
    return line


def _tokens(path: Path) -> list[str]:
    return [token for line in path.read_text(encoding="utf-8").splitlines() for token in line.split()]


def _parse_restrain_amber_rst7(path: Path) -> list[TypedDataset]:
    title, atom_count, time, coordinates, velocities, box = _read_amber_rst7(path)
    datasets = [
        TypedDataset("/parameters/restart/references/restraint/default/source_title", np.asarray(title, dtype=object)),
        TypedDataset("/parameters/restart/references/restraint/default/atom_count", np.asarray(atom_count, dtype=np.int64)),
        TypedDataset("/parameters/restart/references/restraint/default/time", np.asarray(time, dtype=np.float64)),
        TypedDataset("/parameters/restart/references/restraint/default/coordinate", coordinates),
        TypedDataset("/parameters/restart/references/restraint/default/box", box),
    ]
    if velocities is not None:
        datasets.append(TypedDataset("/parameters/restart/references/restraint/default/velocity", velocities))
    return datasets


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


def _parse_meta_edge(path: Path) -> list[TypedDataset]:
    matrix = _read_numeric_matrix(path)
    if matrix.shape[1] < 3 or matrix.shape[1] % 2 != 1:
        raise ValueError(f"{path} meta edge rows must have 2 * ndim + 1 columns")
    ndim = (matrix.shape[1] - 1) // 2
    return [
        TypedDataset("/meta/default/grid/coordinate", matrix[:, :ndim].astype(np.float32)),
        TypedDataset("/meta/default/grid/normalization", matrix[:, ndim].astype(np.float32)),
        TypedDataset("/meta/default/grid/force", matrix[:, ndim + 1 :].astype(np.float32)),
        TypedDataset("/meta/default/grid/ndim", np.asarray(ndim, dtype=np.int32)),
        TypedDataset("/meta/default/grid/count", np.asarray(matrix.shape[0], dtype=np.int64)),
    ]


def _parse_meta_potential(path: Path, *, root: str) -> list[TypedDataset]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 3:
        raise ValueError(f"{path} is too short for a meta potential file")
    title = lines[0]
    ndim = _infer_meta_potential_ndim(lines, path)
    axes = []
    for axis_index in range(ndim):
        row = [float(value) for value in lines[1 + axis_index].split()]
        if len(row) != 3:
            raise ValueError(f"{path} meta axis row {axis_index} must have 3 columns")
        axes.append(row)
    header_values = lines[1 + ndim].split()
    if len(header_values) < ndim + 1:
        raise ValueError(f"{path} meta grid/scatter row must contain ndim grid counts and scatter size")
    grids = np.asarray([int(value) for value in header_values[:ndim]], dtype=np.int32)
    scatter_size = int(header_values[ndim])
    row_lines = lines[2 + ndim :]
    if len(row_lines) != scatter_size:
        raise ValueError(f"{path} declares {scatter_size} meta rows but contains {len(row_lines)}")
    rows = _numeric_rows(row_lines, path)
    coordinates = rows[:, :ndim].astype(np.float32)
    potential = np.zeros(scatter_size, dtype=np.float32)
    force = np.zeros((scatter_size, ndim), dtype=np.float32)
    has_potential = np.zeros(scatter_size, dtype=np.int32)
    has_force = np.zeros(scatter_size, dtype=np.int32)
    for row_index, row in enumerate(rows):
        nwords = row.shape[0]
        if nwords < ndim:
            raise ValueError(f"{path} meta row {row_index} has fewer coordinate columns than ndim")
        if nwords >= ndim + 2:
            potential[row_index] = row[-1]
            has_potential[row_index] = 1
        if nwords == 2 * ndim + 2:
            force[row_index, :] = row[ndim + 1 : ndim + 1 + ndim]
            has_force[row_index] = 1
    return [
        TypedDataset(f"{root}/title", np.asarray(title, dtype=object)),
        TypedDataset(f"{root}/axis/min", np.asarray(axes, dtype=np.float32)[:, 0]),
        TypedDataset(f"{root}/axis/max", np.asarray(axes, dtype=np.float32)[:, 1]),
        TypedDataset(f"{root}/axis/delta", np.asarray(axes, dtype=np.float32)[:, 2]),
        TypedDataset(f"{root}/grid", grids),
        TypedDataset(f"{root}/scatter_size", np.asarray(scatter_size, dtype=np.int64)),
        TypedDataset(f"{root}/coordinate", coordinates),
        TypedDataset(f"{root}/value", potential),
        TypedDataset(f"{root}/force", force),
        TypedDataset(f"{root}/has_value", has_potential),
        TypedDataset(f"{root}/has_force", has_force),
        TypedDataset(f"{root}/ndim", np.asarray(ndim, dtype=np.int32)),
    ]


def _infer_meta_potential_ndim(lines: list[str], path: Path) -> int:
    for ndim in range(1, min(8, len(lines) - 2) + 1):
        axis_ok = True
        for axis_index in range(ndim):
            try:
                row = [float(value) for value in lines[1 + axis_index].split()]
            except ValueError:
                axis_ok = False
                break
            if len(row) != 3:
                axis_ok = False
                break
        if not axis_ok:
            continue
        header = lines[1 + ndim].split()
        if len(header) < ndim + 1:
            continue
        try:
            scatter_size = int(header[ndim])
            for value in header[:ndim]:
                int(value)
        except ValueError:
            continue
        if len(lines[2 + ndim :]) == scatter_size:
            return ndim
    raise ValueError(f"{path} cannot infer meta potential dimensionality")


def _parse_meta_hills(path: Path) -> list[TypedDataset]:
    matrix = _read_numeric_matrix(path)
    return [
        TypedDataset("/parameters/restart/bias/meta/default/hills", path.read_text(encoding="utf-8")),
        TypedDataset("/parameters/restart/bias/meta/default/hills_typed/value", matrix.astype(np.float32)),
        TypedDataset("/parameters/restart/bias/meta/default/hills_typed/count", np.asarray(matrix.shape[0], dtype=np.int64)),
        TypedDataset("/parameters/restart/bias/meta/default/hills_typed/column_count", np.asarray(matrix.shape[1], dtype=np.int32)),
    ]


def _read_numeric_matrix(path: Path) -> np.ndarray:
    rows = [
        [float(value) for value in line.split()]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not rows:
        raise ValueError(f"{path} is empty")
    widths = {len(row) for row in rows}
    if len(widths) != 1:
        raise ValueError(f"{path} has ragged numeric rows")
    return np.asarray(rows, dtype=np.float32)


def _numeric_rows(lines: list[str], path: Path) -> np.ndarray:
    rows = [[float(value) for value in line.split()] for line in lines]
    if not rows:
        return np.empty((0, 0), dtype=np.float32)
    widths = {len(row) for row in rows}
    if len(widths) != 1:
        raise ValueError(f"{path} has ragged numeric rows")
    return np.asarray(rows, dtype=np.float32)


def _read_flat_float_file(path: Path) -> np.ndarray:
    values = _tokens(path)
    if not values:
        raise ValueError(f"{path} is empty")
    return np.asarray([float(value) for value in values], dtype=np.float32)


def _read_flat_int_file(path: Path) -> np.ndarray:
    values = _tokens(path)
    if not values:
        raise ValueError(f"{path} is empty")
    return np.asarray([int(value) for value in values], dtype=np.int32)


def _read_matrix_file(path: Path, *, dtype, columns: int) -> np.ndarray:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    matrix = []
    for row_index, row in enumerate(rows):
        if len(row) != columns:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected {columns}")
        matrix.append([float(value) for value in row])
    return np.asarray(matrix, dtype=dtype)


def _read_counted_xyz(path: Path) -> np.ndarray:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    count = int(rows[0][0])
    xyz_rows = rows[1:]
    if len(xyz_rows) != count:
        raise ValueError(f"{path} declares {count} coordinates but contains {len(xyz_rows)}")
    matrix = []
    for row_index, row in enumerate(xyz_rows):
        if len(row) < 3:
            raise ValueError(f"{path} coordinate row {row_index} has {len(row)} columns, expected at least 3")
        matrix.append([float(value) for value in row[:3]])
    return np.asarray(matrix, dtype=np.float32)


def _read_counted_pairs_with_r0(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not rows:
        raise ValueError(f"{path} is empty")
    count = int(rows[0][0])
    pair_rows = rows[1:]
    if len(pair_rows) != count:
        raise ValueError(f"{path} declares {count} pairs but contains {len(pair_rows)}")
    atoms = []
    r0 = []
    for row_index, row in enumerate(pair_rows):
        if len(row) < 3:
            raise ValueError(f"{path} pair row {row_index} has {len(row)} columns, expected at least 3")
        atoms.append([int(row[0]), int(row[1])])
        r0.append(float(row[2]))
    return np.asarray(atoms, dtype=np.int32), np.asarray(r0, dtype=np.float32)


def _read_soft_wall_config(path: Path) -> tuple[list[str], list[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    sections: list[tuple[str, dict[str, str]]] = []
    index = 0
    while index < len(lines):
        line = lines[index].strip()
        index += 1
        if not line:
            continue
        if not (line.startswith("[[[") and line.endswith("]]]")):
            raise ValueError(f"{path} expected section header '[[[ name ]]]'")
        section = line[3:-3].strip()
        values: dict[str, str] = {}
        while index < len(lines):
            key_line = lines[index].strip()
            index += 1
            if not key_line:
                continue
            if not (key_line.startswith("[[") and key_line.endswith("]]")):
                raise ValueError(f"{path} expected key header '[[ key ]]' in section {section}")
            key = key_line[2:-2].strip()
            if key == "end":
                break
            value_lines = []
            while index < len(lines):
                next_line = lines[index].strip()
                if next_line.startswith("[[") and next_line.endswith("]]"):
                    break
                if next_line:
                    value_lines.append(next_line)
                index += 1
            if not value_lines:
                raise ValueError(f"{path} missing value for section {section} key {key}")
            values[key] = "".join(value_lines)
        else:
            raise ValueError(f"{path} missing '[[ end ]]' for section {section}")
        if "potential" not in values:
            raise ValueError(f"{path} soft wall section {section} missing potential")
        unknown = set(values) - {"potential"}
        if unknown:
            raise ValueError(f"{path} soft wall section {section} has unsupported keys: {sorted(unknown)}")
        sections.append((section, values))
    if not sections:
        raise ValueError(f"{path} does not contain any soft wall sections")
    return [section for section, _ in sections], [values["potential"] for _, values in sections]
