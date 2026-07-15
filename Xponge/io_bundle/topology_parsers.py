"""
Typed parsers for legacy SPONGE topology input files.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


EV_TO_KCAL_MOL = 23.0605480
HARTREE_BOHR_TO_EV_ANGSTROM = 14.3888000


@dataclass(frozen=True)
class TypedDataset:
    """A typed dataset ready to be written into a bundle."""

    path: str
    data: np.ndarray


def parse_topology_file(key: str, path: Path) -> list[TypedDataset] | None:
    """Parse a supported legacy topology file into typed HDF5 datasets."""

    if key == "mass_in_file":
        return [TypedDataset("/atoms/mass", _read_counted_vector(path, np.float32))]
    if key == "charge_in_file":
        return [TypedDataset("/atoms/charge", _read_counted_vector(path, np.float32))]
    if key == "residue_in_file":
        return _parse_residue(path)
    if key == "exclude_in_file":
        return _parse_exclude(path)
    if key == "bond_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=2,
            datasets=(
                ("/forcefield/bond/atoms", np.int32, slice(0, 2)),
                ("/forcefield/bond/k", np.float32, 2),
                ("/forcefield/bond/r0", np.float32, 3),
            ),
        )
    if key == "angle_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=3,
            datasets=(
                ("/forcefield/angle/atoms", np.int32, slice(0, 3)),
                ("/forcefield/angle/k", np.float32, 3),
                ("/forcefield/angle/theta0", np.float32, 4),
            ),
        )
    if key == "dihedral_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=4,
            datasets=(
                ("/forcefield/dihedral/atoms", np.int32, slice(0, 4)),
                ("/forcefield/dihedral/periodicity", np.int32, 4),
                ("/forcefield/dihedral/k", np.float32, 5),
                ("/forcefield/dihedral/phi0", np.float32, 6),
            ),
        )
    if key in {"improper_in_file", "improper_dihedral_in_file"}:
        return _parse_fixed_table(
            path,
            atom_columns=4,
            datasets=(
                ("/forcefield/improper/atoms", np.int32, slice(0, 4)),
                ("/forcefield/improper/pk", np.float32, 4),
                ("/forcefield/improper/phi0", np.float32, 5),
            ),
        )
    if key == "nb14_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=2,
            datasets=(
                ("/forcefield/nb14/atoms", np.int32, slice(0, 2)),
                ("/forcefield/nb14/params", np.float32, slice(2, 4)),
            ),
        )
    if key == "nb14_extra_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=2,
            datasets=(
                ("/forcefield/nb14_extra/atoms", np.int32, slice(0, 2)),
                ("/forcefield/nb14_extra/params", np.float32, slice(2, 5)),
            ),
        )
    if key == "urey_bradley_in_file":
        return _parse_fixed_table(
            path,
            atom_columns=3,
            datasets=(
                ("/forcefield/urey_bradley/atoms", np.int32, slice(0, 3)),
                ("/forcefield/urey_bradley/angle_k", np.float32, 3),
                ("/forcefield/urey_bradley/angle_theta0", np.float32, 4),
                ("/forcefield/urey_bradley/bond_k", np.float32, 5),
                ("/forcefield/urey_bradley/bond_r0", np.float32, 6),
            ),
        )
    if key == "cmap_in_file":
        return _parse_cmap(path)
    if key == "gb_in_file":
        return [TypedDataset("/forcefield/gb/params", _read_counted_matrix(path, np.float32, columns=2))]
    if key in {"virtual_atom_in_file", "virtual_atoms_in_file"}:
        return _parse_virtual_atom(path)
    if key == "LJ_in_file":
        return _parse_lj(path)
    if key == "LJ_soft_core_in_file":
        return _parse_lj_soft_core(path)
    if key == "subsys_division_in_file":
        return [TypedDataset("/forcefield/subsys_division", _read_counted_vector(path, np.int32))]
    if key == "EAM_in_file":
        return _parse_eam(path)
    if key == "EAM_atom_type_in_file":
        return [
            TypedDataset("/manybody/eam/atom_type", _read_flat_int_file(path)),
        ]
    if key == "SW_in_file":
        return _parse_sw(path)
    if key == "EDIP_in_file":
        return _parse_edip(path)
    if key == "TERSOFF_in_file":
        return _parse_tersoff(path)
    if key == "qc_type_in_file":
        return _parse_qc_type(path)
    if key == "REAXFF_in_file":
        return _parse_reaxff_parameters(path)
    if key == "REAXFF_type_in_file":
        return _parse_reaxff_type(path)
    if key == "pairwise_force_in_file":
        return _parse_pairwise_force_config(path)
    if key == "listed_forces_in_file":
        return _parse_listed_forces_config(path)
    return None


def pairwise_force_schema(path: Path) -> tuple[str, list[str], list[str], int]:
    """Return force name, parameter types, parameter names, and *_ij count."""

    sections = _read_configuration_sections(path)
    if len(sections) != 1:
        raise ValueError(f"{path} pairwise force config must contain exactly one section")
    name, values = sections[0]
    _require_config_keys(path, name, values, {"potential", "parameters"})
    parameter_types, parameter_names = _parse_custom_force_parameters(path, name, values["parameters"])
    return name, parameter_types, parameter_names, _pairwise_ij_parameter_count(path, parameter_names)


def listed_force_schemas(path: Path) -> list[tuple[str, list[str], list[str]]]:
    """Return listed force module names and parameter schemas from a config file."""

    schemas = []
    for name, values in _read_configuration_sections(path):
        _require_config_keys(path, name, values, {"potential", "parameters"})
        parameter_types, parameter_names = _parse_custom_force_parameters(path, name, values["parameters"])
        schemas.append((name, parameter_types, parameter_names))
    return schemas


def parse_pairwise_force_data_file(
    force_name: str,
    path: Path,
    *,
    parameter_types: list[str],
    parameter_names: list[str],
    ij_parameter_count: int,
) -> list[TypedDataset]:
    """Parse the dynamic ``<force_name>_in_file`` payload for pairwise forces."""

    tokens = [token for line in _data_lines(path) for token in line]
    if len(tokens) < 2:
        raise ValueError(f"{path} must start with atom and type counts")
    atom_count = int(tokens[0])
    type_count = int(tokens[1])
    pair_count = type_count * (type_count + 1) // 2
    expected = 2 + ij_parameter_count * pair_count + atom_count
    if len(tokens) != expected:
        raise ValueError(f"{path} has {len(tokens)} values, expected {expected}")

    values = np.zeros((ij_parameter_count, pair_count), dtype=np.float32)
    int_values = np.zeros((ij_parameter_count, pair_count), dtype=np.int32)
    float_values = np.full((ij_parameter_count, pair_count), np.nan, dtype=np.float32)
    cursor = 2
    for parameter_index in range(ij_parameter_count):
        row = tokens[cursor : cursor + pair_count]
        cursor += pair_count
        if parameter_types[parameter_index] == "int":
            parsed_ints = np.asarray([int(value) for value in row], dtype=np.int32)
            int_values[parameter_index, :] = parsed_ints
            values[parameter_index, :] = parsed_ints.astype(np.float32)
        else:
            parsed_floats = np.asarray([float(value) for value in row], dtype=np.float32)
            float_values[parameter_index, :] = parsed_floats
            values[parameter_index, :] = parsed_floats
    atom_type = np.asarray([int(value) for value in tokens[cursor:]], dtype=np.int32)
    root = f"/forcefield/custom_force/pairwise/data/{_safe_h5_name(force_name)}"
    return [
        TypedDataset(f"{root}/name", np.asarray(force_name, dtype=object)),
        TypedDataset(f"{root}/atom_count", np.asarray(atom_count, dtype=np.int64)),
        TypedDataset(f"{root}/type_count", np.asarray(type_count, dtype=np.int32)),
        TypedDataset(f"{root}/pair_count", np.asarray(pair_count, dtype=np.int64)),
        TypedDataset(f"{root}/parameter/name", np.asarray(parameter_names[:ij_parameter_count], dtype=object)),
        TypedDataset(f"{root}/parameter/type", np.asarray(parameter_types[:ij_parameter_count], dtype=object)),
        TypedDataset(f"{root}/parameter/value", values),
        TypedDataset(f"{root}/parameter/int_value", int_values),
        TypedDataset(f"{root}/parameter/float_value", float_values),
        TypedDataset(f"{root}/atom_type", atom_type),
    ]


def parse_listed_force_data_file(
    force_name: str,
    path: Path,
    *,
    parameter_types: list[str],
    parameter_names: list[str],
) -> list[TypedDataset]:
    """Parse the dynamic ``<force_name>_in_file`` payload for listed forces."""

    tokens = [token for line in _data_lines(path) for token in line]
    if not tokens:
        raise ValueError(f"{path} must start with item count")
    item_count = int(tokens[0])
    parameter_count = len(parameter_names)
    expected = 1 + item_count * parameter_count
    if len(tokens) != expected:
        raise ValueError(f"{path} has {len(tokens)} values, expected {expected}")

    rows = []
    int_rows = []
    float_rows = []
    cursor = 1
    for _ in range(item_count):
        row = []
        int_row = []
        float_row = []
        for parameter_type in parameter_types:
            value = tokens[cursor]
            cursor += 1
            if parameter_type == "int":
                parsed_int = int(value)
                row.append(parsed_int)
                int_row.append(parsed_int)
                float_row.append(np.nan)
            else:
                parsed_float = float(value)
                row.append(parsed_float)
                int_row.append(0)
                float_row.append(parsed_float)
        rows.append(row)
        int_rows.append(int_row)
        float_rows.append(float_row)
    value = np.asarray(rows, dtype=np.float32).reshape(item_count, parameter_count)
    int_value = np.asarray(int_rows, dtype=np.int32).reshape(item_count, parameter_count)
    float_value = np.asarray(float_rows, dtype=np.float32).reshape(item_count, parameter_count)
    int_mask = np.asarray([parameter_type == "int" for parameter_type in parameter_types], dtype=np.bool_)
    root = f"/forcefield/custom_force/listed/data/{_safe_h5_name(force_name)}"
    return [
        TypedDataset(f"{root}/name", np.asarray(force_name, dtype=object)),
        TypedDataset(f"{root}/item_count", np.asarray(item_count, dtype=np.int64)),
        TypedDataset(f"{root}/parameter/name", np.asarray(parameter_names, dtype=object)),
        TypedDataset(f"{root}/parameter/type", np.asarray(parameter_types, dtype=object)),
        TypedDataset(f"{root}/parameter/is_int", int_mask),
        TypedDataset(f"{root}/parameter/value", value),
        TypedDataset(f"{root}/parameter/int_value", int_value),
        TypedDataset(f"{root}/parameter/float_value", float_value),
    ]


def _data_lines(path: Path) -> list[list[str]]:
    return [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _non_empty_raw_lines(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _read_counted_vector(path: Path, dtype) -> np.ndarray:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    count = int(lines[0][0])
    values = np.asarray([float(line[0]) for line in lines[1:]], dtype=dtype)
    if values.shape != (count,):
        raise ValueError(f"{path} declares {count} values but contains {values.shape[0]}")
    return values


def _read_counted_matrix(path: Path, dtype, *, columns: int) -> np.ndarray:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    count = int(lines[0][0])
    rows = lines[1:]
    if len(rows) != count:
        raise ValueError(f"{path} declares {count} rows but contains {len(rows)}")
    matrix = []
    for row_index, row in enumerate(rows):
        if len(row) != columns:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected {columns}")
        matrix.append([float(value) for value in row])
    return np.asarray(matrix, dtype=dtype)


def _read_flat_int_file(path: Path) -> np.ndarray:
    values = [int(value) for line in _data_lines(path) for value in line]
    if not values:
        raise ValueError(f"{path} is empty")
    return np.asarray(values, dtype=np.int32)


def _parse_residue(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    atom_count, residue_count = int(lines[0][0]), int(lines[0][1])
    counts = np.asarray([int(line[0]) for line in lines[1:]], dtype=np.int64)
    if counts.shape != (residue_count,):
        raise ValueError(f"{path} declares {residue_count} residues but contains {counts.shape[0]}")
    offsets = np.zeros(residue_count + 1, dtype=np.int64)
    offsets[1:] = np.cumsum(counts)
    if int(offsets[-1]) != atom_count:
        raise ValueError(f"{path} residue atom counts sum to {offsets[-1]}, expected {atom_count}")
    residue_index = np.empty(atom_count, dtype=np.int32)
    start = 0
    for index, count in enumerate(counts):
        residue_index[start : start + int(count)] = index
        start += int(count)
    return [
        TypedDataset("/residues/atom_offset", offsets),
        TypedDataset("/atoms/residue_index", residue_index),
    ]


def _parse_exclude(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    atom_count = int(lines[0][0])
    offsets = np.zeros(atom_count + 1, dtype=np.int64)
    values: list[int] = []
    if len(lines[1:]) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atoms but contains {len(lines) - 1} rows")
    for atom_index, line in enumerate(lines[1:]):
        count = int(line[0])
        row = [int(value) for value in line[1:]]
        if len(row) != count:
            raise ValueError(f"{path} exclude row {atom_index} declares {count} values but contains {len(row)}")
        values.extend(row)
        offsets[atom_index + 1] = len(values)
    return [
        TypedDataset("/topology/exclusions/offset", offsets),
        TypedDataset("/topology/exclusions/list", np.asarray(values, dtype=np.int32)),
    ]


def _parse_fixed_table(path: Path, *, atom_columns: int, datasets) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    count = int(lines[0][0])
    rows = lines[1:]
    if len(rows) != count:
        raise ValueError(f"{path} declares {count} records but contains {len(rows)}")
    min_columns = max(
        selector.stop if isinstance(selector, slice) else selector + 1
        for _, _, selector in datasets
    )
    table = []
    for row_index, row in enumerate(rows):
        if len(row) < min_columns:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected at least {min_columns}")
        table.append([float(value) for value in row])
    array = np.asarray(table, dtype=np.float64)
    parsed = []
    for dataset_path, dtype, selector in datasets:
        parsed.append(TypedDataset(dataset_path, array[:, selector].astype(dtype)))
    parsed.append(TypedDataset(str(Path(datasets[0][0]).parent / "count").replace("\\", "/"), np.asarray(count, dtype=np.int64)))
    return parsed


def _parse_lj(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    atom_count, type_count = int(lines[0][0]), int(lines[0][1])
    values = [float(item) for line in lines[1:] for item in line]
    pair_count = type_count * (type_count + 1) // 2
    expected = pair_count * 2 + atom_count
    if len(values) != expected:
        raise ValueError(f"{path} LJ payload has {len(values)} values, expected {expected}")
    pair_a = np.asarray(values[:pair_count], dtype=np.float32)
    pair_b = np.asarray(values[pair_count : pair_count * 2], dtype=np.float32)
    atom_type = np.asarray(values[pair_count * 2 :], dtype=np.int32)
    return [
        TypedDataset("/forcefield/lj/atom_type_count", np.asarray(type_count, dtype=np.int32)),
        TypedDataset("/forcefield/lj/pair_A_12", pair_a),
        TypedDataset("/forcefield/lj/pair_B_6", pair_b),
        TypedDataset("/forcefield/lj/type", atom_type),
    ]


def _parse_cmap(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    cmap_count, type_count = int(lines[0][0]), int(lines[0][1])
    if len(lines) < 2:
        raise ValueError(f"{path} missing CMAP resolution row")
    resolution = np.asarray([int(value) for value in lines[1]], dtype=np.int32)
    if resolution.shape != (type_count,):
        raise ValueError(f"{path} declares {type_count} CMAP types but has {resolution.shape[0]} resolutions")

    flat = [value for line in lines[2:] for value in line]
    grid_count = int(np.sum(resolution.astype(np.int64) ** 2))
    record_count = cmap_count * 6
    expected = grid_count + record_count
    if len(flat) != expected:
        raise ValueError(f"{path} CMAP payload has {len(flat)} values, expected {expected}")
    grid_value = np.asarray([float(value) for value in flat[:grid_count]], dtype=np.float32)
    records = np.asarray([int(value) for value in flat[grid_count:]], dtype=np.int32).reshape(cmap_count, 6)
    return [
        TypedDataset("/forcefield/cmap/atoms", records[:, :5]),
        TypedDataset("/forcefield/cmap/type", records[:, 5]),
        TypedDataset("/forcefield/cmap/resolution", resolution),
        TypedDataset("/forcefield/cmap/grid_value", grid_value),
        TypedDataset("/forcefield/cmap/count", np.asarray(cmap_count, dtype=np.int64)),
    ]


def _parse_virtual_atom(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    types = []
    atoms = []
    from_values: list[int] = []
    from_offsets = [0]
    parameter_values: list[float] = []
    parameter_offsets = [0]
    for row_index, row in enumerate(lines):
        vtype = int(row[0])
        atom = int(row[1])
        if vtype == 0:
            expected_columns = 4
            from_count = 1
        elif vtype == 1:
            expected_columns = 5
            from_count = 2
        elif vtype in (2, 3):
            expected_columns = 7
            from_count = 3
        else:
            raise ValueError(f"{path} row {row_index} has unsupported virtual atom type {vtype}")
        if len(row) != expected_columns:
            raise ValueError(f"{path} row {row_index} has {len(row)} columns, expected {expected_columns}")
        types.append(vtype)
        atoms.append(atom)
        from_values.extend(int(value) for value in row[2 : 2 + from_count])
        parameter_values.extend(float(value) for value in row[2 + from_count :])
        from_offsets.append(len(from_values))
        parameter_offsets.append(len(parameter_values))
    return [
        TypedDataset("/forcefield/virtual_atom/type", np.asarray(types, dtype=np.int32)),
        TypedDataset("/forcefield/virtual_atom/atom", np.asarray(atoms, dtype=np.int32)),
        TypedDataset("/forcefield/virtual_atom/from_offset", np.asarray(from_offsets, dtype=np.int64)),
        TypedDataset("/forcefield/virtual_atom/from", np.asarray(from_values, dtype=np.int32)),
        TypedDataset("/forcefield/virtual_atom/parameter_offset", np.asarray(parameter_offsets, dtype=np.int64)),
        TypedDataset("/forcefield/virtual_atom/parameter", np.asarray(parameter_values, dtype=np.float32)),
        TypedDataset("/forcefield/virtual_atom/count", np.asarray(len(types), dtype=np.int64)),
    ]


def _parse_lj_soft_core(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    atom_count, type_count_a, type_count_b = (int(value) for value in lines[0][:3])
    pair_count_a = type_count_a * (type_count_a + 1) // 2
    pair_count_b = type_count_b * (type_count_b + 1) // 2
    values = [float(item) for line in lines[1:] for item in line]
    expected = pair_count_a * 2 + pair_count_b * 2 + atom_count * 2
    if len(values) != expected:
        raise ValueError(f"{path} LJ soft-core payload has {len(values)} values, expected {expected}")
    offset = 0
    pair_aa = np.asarray(values[offset : offset + pair_count_a], dtype=np.float32)
    offset += pair_count_a
    pair_ab = np.asarray(values[offset : offset + pair_count_a], dtype=np.float32)
    offset += pair_count_a
    pair_ba = np.asarray(values[offset : offset + pair_count_b], dtype=np.float32)
    offset += pair_count_b
    pair_bb = np.asarray(values[offset : offset + pair_count_b], dtype=np.float32)
    offset += pair_count_b
    atom_types = np.asarray(values[offset:], dtype=np.int32).reshape(atom_count, 2)
    return [
        TypedDataset("/forcefield/lj_soft_core/atom_type_A", atom_types[:, 0]),
        TypedDataset("/forcefield/lj_soft_core/atom_type_B", atom_types[:, 1]),
        TypedDataset("/forcefield/lj_soft_core/atom_type_count_A", np.asarray(type_count_a, dtype=np.int32)),
        TypedDataset("/forcefield/lj_soft_core/atom_type_count_B", np.asarray(type_count_b, dtype=np.int32)),
        TypedDataset("/forcefield/lj_soft_core/pair_AA", pair_aa),
        TypedDataset("/forcefield/lj_soft_core/pair_AB", pair_ab),
        TypedDataset("/forcefield/lj_soft_core/pair_BA", pair_ba),
        TypedDataset("/forcefield/lj_soft_core/pair_BB", pair_bb),
    ]


def _parse_manybody_table(
    path: Path,
    *,
    root: str,
    pair_param_columns: int,
    triple_param_columns: int,
) -> list[TypedDataset]:
    lines = _non_empty_raw_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    header = lines[0].split()
    if len(header) < 2:
        raise ValueError(f"{path} header must declare atom and type counts")
    atom_count, type_count = int(header[0]), int(header[1])
    pair_count = type_count * type_count
    triple_count = type_count * type_count * type_count

    cursor = 1
    cursor = _require_comment(lines, cursor, path, "first")
    pair_rows = _parse_manybody_rows(
        lines[cursor : cursor + pair_count],
        path,
        label="pair",
        index_columns=2,
        param_columns=pair_param_columns,
    )
    cursor += pair_count

    cursor = _require_comment(lines, cursor, path, "second")
    triple_rows = _parse_manybody_rows(
        lines[cursor : cursor + triple_count],
        path,
        label="triple",
        index_columns=3,
        param_columns=triple_param_columns,
    )
    cursor += triple_count

    cursor = _require_comment(lines, cursor, path, "third")
    atom_type_values = [int(value) for line in lines[cursor:] for value in line.split()]
    if len(atom_type_values) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atom types but contains {len(atom_type_values)}")
    atom_types = np.asarray(atom_type_values, dtype=np.int32)
    if np.any(atom_types >= type_count) or np.any(atom_types < 0):
        raise ValueError(f"{path} contains atom type outside [0, {type_count})")

    pair_index = pair_rows[:, :2].astype(np.int32)
    pair_params = pair_rows[:, 2:].astype(np.float32)
    triple_index = triple_rows[:, :3].astype(np.int32)
    triple_params = triple_rows[:, 3:].astype(np.float32)
    _validate_manybody_coverage(pair_index, type_count, path, "pair")
    _validate_manybody_coverage(triple_index, type_count, path, "triple")

    return [
        TypedDataset(f"{root}/atom_type_count", np.asarray(type_count, dtype=np.int32)),
        TypedDataset(f"{root}/pair/type", pair_index),
        TypedDataset(f"{root}/pair/parameters", pair_params),
        TypedDataset(f"{root}/triple/type", triple_index),
        TypedDataset(f"{root}/triple/parameters", triple_params),
        TypedDataset(f"{root}/atom_type", atom_types),
    ]


def _require_comment(lines: list[str], index: int, path: Path, label: str) -> int:
    if index >= len(lines) or not lines[index].startswith("#"):
        raise ValueError(f"{path} missing {label} comment line")
    return index + 1


def _parse_manybody_rows(
    rows: list[str],
    path: Path,
    *,
    label: str,
    index_columns: int,
    param_columns: int,
) -> np.ndarray:
    expected_columns = index_columns + param_columns
    if len(rows) == 0:
        raise ValueError(f"{path} missing {label} parameter rows")
    parsed = []
    for row_index, row in enumerate(rows):
        columns = row.split()
        if len(columns) != expected_columns:
            raise ValueError(
                f"{path} {label} row {row_index} has {len(columns)} columns, expected {expected_columns}"
            )
        parsed.append([float(value) for value in columns])
    return np.asarray(parsed, dtype=np.float64)


def _validate_manybody_coverage(index: np.ndarray, type_count: int, path: Path, label: str) -> None:
    multipliers = np.asarray([type_count ** power for power in range(index.shape[1] - 1, -1, -1)], dtype=np.int64)
    if np.any(index < 0) or np.any(index >= type_count):
        raise ValueError(f"{path} contains {label} type index outside [0, {type_count})")
    flattened = index.astype(np.int64) @ multipliers
    if np.unique(flattened).shape[0] != flattened.shape[0]:
        raise ValueError(f"{path} contains duplicate {label} type rows")


def _parse_sw(path: Path) -> list[TypedDataset]:
    return _parse_manybody_table(
        path,
        root="/manybody/sw",
        pair_param_columns=8,
        triple_param_columns=3,
    )


def _parse_edip(path: Path) -> list[TypedDataset]:
    return _parse_manybody_table(
        path,
        root="/manybody/edip",
        pair_param_columns=8,
        triple_param_columns=9,
    )


def _parse_eam(path: Path) -> list[TypedDataset]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path} is too short for an EAM file")
    second = lines[1].split()
    if len(second) >= 2 and _is_int_token(second[0]) and _is_float_token(second[1]):
        return _parse_eam_funcfl(path, lines)
    return _parse_eam_setfl(path, lines)


def _parse_eam_funcfl(path: Path, lines: list[str]) -> list[TypedDataset]:
    header = lines[1].split()
    if len(header) < 4:
        raise ValueError(f"{path} funcfl header must contain atomic number, mass, lattice constant, and lattice type")
    atomic_number = int(header[0])
    mass = float(header[1])
    lattice_constant = float(header[2])
    lattice_type = header[3]
    nrho, drho, nr, dr, cut, cursor = _parse_eam_table_parameters(path, lines, 2)
    values = _float_tokens(lines[cursor:])
    expected = nrho + nr + nr
    if len(values) < expected:
        raise ValueError(f"{path} funcfl payload has {len(values)} values, expected at least {expected}")
    offset = 0
    embed_raw = np.asarray(values[offset : offset + nrho], dtype=np.float32).reshape(1, nrho)
    offset += nrho
    z_raw = np.asarray(values[offset : offset + nr], dtype=np.float32)
    offset += nr
    electron_density = np.asarray(values[offset : offset + nr], dtype=np.float32).reshape(1, nr)
    pair_potential = _eam_funcfl_pair_potential(z_raw, dr).reshape(1, 1, nr)
    return _eam_common_datasets(
        format_name="funcfl",
        type_count=1,
        nrho=nrho,
        drho=drho,
        nr=nr,
        dr=dr,
        cut=cut,
        embed_raw=embed_raw,
        electron_density=electron_density,
        pair_potential=pair_potential,
    ) + [
        TypedDataset("/manybody/eam/atomic_number", np.asarray([atomic_number], dtype=np.int32)),
        TypedDataset("/manybody/eam/mass", np.asarray([mass], dtype=np.float32)),
        TypedDataset("/manybody/eam/lattice_constant", np.asarray([lattice_constant], dtype=np.float32)),
        TypedDataset("/manybody/eam/lattice_type", np.asarray([lattice_type], dtype=object)),
        TypedDataset("/manybody/eam/funcfl/z", z_raw),
    ]


def _parse_eam_setfl(path: Path, lines: list[str]) -> list[TypedDataset]:
    if len(lines) < 5:
        raise ValueError(f"{path} is too short for an EAM setfl file")
    type_line = lines[3].split()
    if not type_line:
        raise ValueError(f"{path} missing setfl type line")
    type_count = int(type_line[0])
    type_names = type_line[1:]
    if len(type_names) != type_count:
        raise ValueError(f"{path} declares {type_count} EAM types but contains {len(type_names)} names")
    nrho, drho, nr, dr, cut, cursor = _parse_eam_table_parameters(path, lines, 4)
    atomic_numbers = []
    masses = []
    lattice_constants = []
    lattice_types = []
    embed_rows = []
    rho_rows = []
    for type_index in range(type_count):
        if cursor >= len(lines):
            raise ValueError(f"{path} missing setfl element header {type_index}")
        header = lines[cursor].split()
        cursor += 1
        if len(header) < 4:
            raise ValueError(f"{path} setfl element header {type_index} must contain 4 columns")
        atomic_numbers.append(int(header[0]))
        masses.append(float(header[1]))
        lattice_constants.append(float(header[2]))
        lattice_types.append(header[3])
        values = _float_tokens_from_cursor(lines, cursor, nrho + nr)
        embed_rows.append(values[:nrho])
        rho_rows.append(values[nrho : nrho + nr])
        cursor = values.cursor
    pair_values = _float_tokens(lines[cursor:])
    pair_count = type_count * (type_count + 1) // 2
    expected_pair_values = pair_count * nr
    if len(pair_values) < expected_pair_values:
        raise ValueError(f"{path} setfl pair payload has {len(pair_values)} values, expected at least {expected_pair_values}")
    pair_potential = np.zeros((type_count, type_count, nr), dtype=np.float32)
    offset = 0
    for i in range(type_count):
        for j in range(i + 1):
            raw = np.asarray(pair_values[offset : offset + nr], dtype=np.float32)
            offset += nr
            converted = _eam_setfl_pair_potential(raw, dr)
            pair_potential[i, j, :] = converted
            pair_potential[j, i, :] = converted
    return _eam_common_datasets(
        format_name="setfl",
        type_count=type_count,
        nrho=nrho,
        drho=drho,
        nr=nr,
        dr=dr,
        cut=cut,
        embed_raw=np.asarray(embed_rows, dtype=np.float32),
        electron_density=np.asarray(rho_rows, dtype=np.float32),
        pair_potential=pair_potential,
    ) + [
        TypedDataset("/manybody/eam/type_name", np.asarray(type_names, dtype=object)),
        TypedDataset("/manybody/eam/atomic_number", np.asarray(atomic_numbers, dtype=np.int32)),
        TypedDataset("/manybody/eam/mass", np.asarray(masses, dtype=np.float32)),
        TypedDataset("/manybody/eam/lattice_constant", np.asarray(lattice_constants, dtype=np.float32)),
        TypedDataset("/manybody/eam/lattice_type", np.asarray(lattice_types, dtype=object)),
    ]


def _parse_eam_table_parameters(path: Path, lines: list[str], index: int) -> tuple[int, float, int, float, float, int]:
    if index >= len(lines):
        raise ValueError(f"{path} missing EAM table parameter line")
    row = lines[index].split()
    if len(row) < 5:
        raise ValueError(f"{path} EAM table parameter line must contain nrho, drho, nr, dr, and cut")
    return int(row[0]), float(row[1]), int(row[2]), float(row[3]), float(row[4]), index + 1


def _eam_common_datasets(
    *,
    format_name: str,
    type_count: int,
    nrho: int,
    drho: float,
    nr: int,
    dr: float,
    cut: float,
    embed_raw: np.ndarray,
    electron_density: np.ndarray,
    pair_potential: np.ndarray,
) -> list[TypedDataset]:
    return [
        TypedDataset("/manybody/eam/format", np.asarray(format_name, dtype=object)),
        TypedDataset("/manybody/eam/atom_type_count", np.asarray(type_count, dtype=np.int32)),
        TypedDataset("/manybody/eam/nrho", np.asarray(nrho, dtype=np.int32)),
        TypedDataset("/manybody/eam/drho", np.asarray(drho, dtype=np.float32)),
        TypedDataset("/manybody/eam/nr", np.asarray(nr, dtype=np.int32)),
        TypedDataset("/manybody/eam/dr", np.asarray(dr, dtype=np.float32)),
        TypedDataset("/manybody/eam/cut", np.asarray(cut, dtype=np.float32)),
        TypedDataset("/manybody/eam/embed/raw_ev", embed_raw.astype(np.float32)),
        TypedDataset("/manybody/eam/embed/value", (embed_raw * EV_TO_KCAL_MOL).astype(np.float32)),
        TypedDataset("/manybody/eam/electron_density/value", electron_density.astype(np.float32)),
        TypedDataset("/manybody/eam/pair_potential/value", pair_potential.astype(np.float32)),
    ]


def _eam_funcfl_pair_potential(z_values: np.ndarray, dr: float) -> np.ndarray:
    converted = np.zeros(z_values.shape, dtype=np.float32)
    for index, z_value in enumerate(z_values):
        r = index * dr
        if index == 0:
            r = 1.0e-8
        converted[index] = (z_value * z_value / r) * HARTREE_BOHR_TO_EV_ANGSTROM * EV_TO_KCAL_MOL
    return converted


def _eam_setfl_pair_potential(values: np.ndarray, dr: float) -> np.ndarray:
    converted = np.zeros(values.shape, dtype=np.float32)
    for index, value in enumerate(values):
        r = index * dr
        if index == 0:
            r = 1.0e-8
        converted[index] = (value / r) * EV_TO_KCAL_MOL
    return converted


class _FloatTokens(list):
    cursor: int


def _float_tokens_from_cursor(lines: list[str], cursor: int, count: int) -> _FloatTokens:
    values = _FloatTokens()
    while cursor < len(lines) and len(values) < count:
        values.extend(float(value) for value in lines[cursor].split())
        cursor += 1
    if len(values) < count:
        raise ValueError("EAM table ended before expected value count")
    del values[count:]
    values.cursor = cursor
    return values


def _float_tokens(lines: list[str]) -> list[float]:
    return [float(value) for line in lines for value in line.split()]


def _parse_qc_type(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    header = lines[0]
    if len(header) < 3:
        raise ValueError(f"{path} header must declare natm, charge, and multiplicity")
    atom_count, charge, multiplicity = int(header[0]), int(header[1]), int(header[2])
    atom_rows = lines[1:]
    if len(atom_rows) != atom_count:
        raise ValueError(f"{path} declares {atom_count} QM atoms but contains {len(atom_rows)}")
    atom_index = []
    symbol = []
    for row_index, row in enumerate(atom_rows):
        if len(row) < 2:
            raise ValueError(f"{path} atom row {row_index} must contain MD index and element symbol")
        atom_index.append(int(row[0]))
        symbol.append(row[1])
    return [
        TypedDataset("/qc/type/atom_index", np.asarray(atom_index, dtype=np.int32)),
        TypedDataset("/qc/type/symbol", np.asarray(symbol, dtype=object)),
        TypedDataset("/qc/type/charge", np.asarray(charge, dtype=np.int32)),
        TypedDataset("/qc/type/multiplicity", np.asarray(multiplicity, dtype=np.int32)),
        TypedDataset("/qc/type/count", np.asarray(atom_count, dtype=np.int64)),
    ]


def _parse_reaxff_type(path: Path) -> list[TypedDataset]:
    lines = _data_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    atom_count = int(lines[0][0])
    type_rows = lines[1:]
    if len(type_rows) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atom types but contains {len(type_rows)}")
    names = []
    for row_index, row in enumerate(type_rows):
        if len(row) < 1:
            raise ValueError(f"{path} atom type row {row_index} is empty")
        names.append(row[0])
    return [
        TypedDataset("/manybody/reaxff/type/name", np.asarray(names, dtype=object)),
        TypedDataset("/manybody/reaxff/type/count", np.asarray(atom_count, dtype=np.int64)),
    ]


def _parse_reaxff_parameters(path: Path) -> list[TypedDataset]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path} is too short for a REAXFF parameter file")

    cursor = 0
    header = lines[cursor].strip()
    cursor += 1

    general_count, label, cursor = _read_reaxff_count(lines, cursor, path, "general parameter")
    general_values = []
    general_labels = []
    for index in range(general_count):
        columns, line_label, cursor = _read_reaxff_numeric_line(lines, cursor, path, f"general parameter {index}")
        if not columns:
            raise ValueError(f"{path} general parameter {index} is empty")
        general_values.append(columns[0])
        general_labels.append(line_label)

    atom_type_count, atom_type_label, cursor = _read_reaxff_count(lines, cursor, path, "atom type")
    atom_header_lines = _read_reaxff_header_lines(lines, cursor, 3, path, "atom type header")
    cursor += 3
    type_names: list[str] = []
    atom_values: list[float] = []
    atom_value_offsets = [0]
    atom_line_value_offsets = [0]
    atom_line_labels: list[str] = []
    for type_index in range(atom_type_count):
        columns, line_label, cursor = _read_reaxff_data_line(lines, cursor, path, f"atom type {type_index} line 1")
        if not columns or _is_float_token(columns[0]):
            raise ValueError(f"{path} atom type {type_index} line 1 must start with an element name")
        type_names.append(columns[0])
        line_values = [float(value) for value in columns[1:]]
        atom_values.extend(line_values)
        atom_line_value_offsets.append(len(atom_values))
        atom_line_labels.append(line_label)
        for line_index in range(1, 4):
            values, line_label, cursor = _read_reaxff_numeric_line(
                lines,
                cursor,
                path,
                f"atom type {type_index} line {line_index + 1}",
            )
            atom_values.extend(values)
            atom_line_value_offsets.append(len(atom_values))
            atom_line_labels.append(line_label)
        atom_value_offsets.append(len(atom_values))

    bond_count, bond_label, cursor = _read_reaxff_count(lines, cursor, path, "bond")
    bond_header_lines = _read_reaxff_header_lines(lines, cursor, 1, path, "bond header")
    cursor += 1
    bond_types: list[list[int]] = []
    bond_values: list[float] = []
    bond_value_offsets = [0]
    bond_line_value_offsets = [0]
    bond_line_labels: list[str] = []
    for bond_index in range(bond_count):
        columns, line_label, cursor = _read_reaxff_numeric_line(lines, cursor, path, f"bond {bond_index} line 1")
        if len(columns) < 2:
            raise ValueError(f"{path} bond {bond_index} line 1 must contain two atom type indices")
        bond_types.append([int(columns[0]), int(columns[1])])
        bond_values.extend(columns[2:])
        bond_line_value_offsets.append(len(bond_values))
        bond_line_labels.append(line_label)
        values, line_label, cursor = _read_reaxff_numeric_line(lines, cursor, path, f"bond {bond_index} line 2")
        bond_values.extend(values)
        bond_line_value_offsets.append(len(bond_values))
        bond_line_labels.append(line_label)
        bond_value_offsets.append(len(bond_values))

    offdiag_count, offdiag_label, cursor = _read_reaxff_count(lines, cursor, path, "off-diagonal")
    offdiag_types, offdiag_values, cursor = _read_reaxff_indexed_rows(
        lines,
        cursor,
        path,
        row_count=offdiag_count,
        index_columns=2,
        section="off-diagonal",
    )

    angle_count, angle_label, cursor = _read_reaxff_count(lines, cursor, path, "angle")
    angle_types, angle_values, cursor = _read_reaxff_indexed_rows(
        lines,
        cursor,
        path,
        row_count=angle_count,
        index_columns=3,
        section="angle",
    )

    torsion_count, torsion_label, cursor = _read_reaxff_count(lines, cursor, path, "torsion")
    torsion_types, torsion_values, cursor = _read_reaxff_indexed_rows(
        lines,
        cursor,
        path,
        row_count=torsion_count,
        index_columns=4,
        section="torsion",
    )

    hydrogen_bond_count, hydrogen_bond_label, cursor = _read_reaxff_count(lines, cursor, path, "hydrogen bond")
    hydrogen_bond_types, hydrogen_bond_values, cursor = _read_reaxff_indexed_rows(
        lines,
        cursor,
        path,
        row_count=hydrogen_bond_count,
        index_columns=3,
        section="hydrogen_bond",
    )

    datasets = [
        TypedDataset("/manybody/reaxff/parameters/header", np.asarray(header, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/general/count", np.asarray(general_count, dtype=np.int64)),
        TypedDataset("/manybody/reaxff/parameters/general/count_label", np.asarray(label, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/general/value", np.asarray(general_values, dtype=np.float32)),
        TypedDataset("/manybody/reaxff/parameters/general/label", np.asarray(general_labels, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/atom/count", np.asarray(atom_type_count, dtype=np.int64)),
        TypedDataset("/manybody/reaxff/parameters/atom/count_label", np.asarray(atom_type_label, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/atom/header", np.asarray(atom_header_lines, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/atom/type_name", np.asarray(type_names, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/atom/value", np.asarray(atom_values, dtype=np.float32)),
        TypedDataset("/manybody/reaxff/parameters/atom/value_offset", np.asarray(atom_value_offsets, dtype=np.int64)),
        TypedDataset(
            "/manybody/reaxff/parameters/atom/line_value_offset",
            np.asarray(atom_line_value_offsets, dtype=np.int64),
        ),
        TypedDataset("/manybody/reaxff/parameters/atom/line_label", np.asarray(atom_line_labels, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/bond/count", np.asarray(bond_count, dtype=np.int64)),
        TypedDataset("/manybody/reaxff/parameters/bond/count_label", np.asarray(bond_label, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/bond/header", np.asarray(bond_header_lines, dtype=object)),
        TypedDataset("/manybody/reaxff/parameters/bond/type", np.asarray(bond_types, dtype=np.int32).reshape(bond_count, 2)),
        TypedDataset("/manybody/reaxff/parameters/bond/value", np.asarray(bond_values, dtype=np.float32)),
        TypedDataset("/manybody/reaxff/parameters/bond/value_offset", np.asarray(bond_value_offsets, dtype=np.int64)),
        TypedDataset(
            "/manybody/reaxff/parameters/bond/line_value_offset",
            np.asarray(bond_line_value_offsets, dtype=np.int64),
        ),
        TypedDataset("/manybody/reaxff/parameters/bond/line_label", np.asarray(bond_line_labels, dtype=object)),
    ]
    datasets.extend(
        _reaxff_row_section_datasets(
            "off_diagonal",
            offdiag_count,
            offdiag_label,
            offdiag_types,
            offdiag_values,
        )
    )
    datasets.extend(_reaxff_row_section_datasets("angle", angle_count, angle_label, angle_types, angle_values))
    datasets.extend(_reaxff_row_section_datasets("torsion", torsion_count, torsion_label, torsion_types, torsion_values))
    datasets.extend(
        _reaxff_row_section_datasets(
            "hydrogen_bond",
            hydrogen_bond_count,
            hydrogen_bond_label,
            hydrogen_bond_types,
            hydrogen_bond_values,
        )
    )
    return datasets


def _read_reaxff_count(lines: list[str], cursor: int, path: Path, section: str) -> tuple[int, str, int]:
    columns, label, cursor = _read_reaxff_data_line(lines, cursor, path, f"{section} count line")
    if not columns or not _is_int_token(columns[0]):
        raise ValueError(f"{path} {section} count line must start with an integer")
    return int(columns[0]), label, cursor


def _read_reaxff_header_lines(lines: list[str], cursor: int, count: int, path: Path, section: str) -> list[str]:
    if cursor + count > len(lines):
        raise ValueError(f"{path} missing {section}")
    return [lines[cursor + index].strip() for index in range(count)]


def _read_reaxff_data_line(lines: list[str], cursor: int, path: Path, section: str) -> tuple[list[str], str, int]:
    if cursor >= len(lines):
        raise ValueError(f"{path} missing {section}")
    content, label = _split_reaxff_label(lines[cursor])
    columns = content.split()
    if not columns:
        raise ValueError(f"{path} {section} is empty")
    return columns, label, cursor + 1


def _read_reaxff_numeric_line(lines: list[str], cursor: int, path: Path, section: str) -> tuple[list[float], str, int]:
    columns, label, cursor = _read_reaxff_data_line(lines, cursor, path, section)
    values = []
    for column in columns:
        if not _is_float_token(column):
            raise ValueError(f"{path} {section} contains non-numeric token {column!r}")
        values.append(float(column))
    return values, label, cursor


def _split_reaxff_label(line: str) -> tuple[str, str]:
    content, separator, label = line.partition("!")
    if not separator:
        return content.strip(), ""
    return content.strip(), label.strip()


def _read_reaxff_indexed_rows(
    lines: list[str],
    cursor: int,
    path: Path,
    *,
    row_count: int,
    index_columns: int,
    section: str,
) -> tuple[np.ndarray, np.ndarray, int]:
    types = []
    values = []
    for row_index in range(row_count):
        row_values, _, cursor = _read_reaxff_numeric_line(lines, cursor, path, f"{section} row {row_index}")
        if len(row_values) < index_columns:
            raise ValueError(f"{path} {section} row {row_index} must contain {index_columns} type indices")
        types.append([int(value) for value in row_values[:index_columns]])
        values.append(row_values[index_columns:])
    return (
        np.asarray(types, dtype=np.int32).reshape(row_count, index_columns),
        _rectangular_float_matrix(values),
        cursor,
    )


def _rectangular_float_matrix(rows: list[list[float]]) -> np.ndarray:
    if not rows:
        return np.empty((0, 0), dtype=np.float32)
    width = max(len(row) for row in rows)
    matrix = np.full((len(rows), width), np.nan, dtype=np.float32)
    for row_index, row in enumerate(rows):
        matrix[row_index, : len(row)] = row
    return matrix


def _reaxff_row_section_datasets(
    name: str,
    count: int,
    label: str,
    types: np.ndarray,
    values: np.ndarray,
) -> list[TypedDataset]:
    root = f"/manybody/reaxff/parameters/{name}"
    return [
        TypedDataset(f"{root}/count", np.asarray(count, dtype=np.int64)),
        TypedDataset(f"{root}/count_label", np.asarray(label, dtype=object)),
        TypedDataset(f"{root}/type", types),
        TypedDataset(f"{root}/value", values),
    ]


def _parse_pairwise_force_config(path: Path) -> list[TypedDataset]:
    sections = _read_configuration_sections(path)
    if len(sections) != 1:
        raise ValueError(f"{path} pairwise force config must contain exactly one section")
    name, values = sections[0]
    _require_config_keys(path, name, values, {"potential", "parameters"})
    unknown = set(values) - {"potential", "parameters", "with_ele", "electrostatic_potential"}
    if unknown:
        raise ValueError(f"{path} pairwise force section {name} has unsupported keys: {sorted(unknown)}")
    parameter_type, parameter_name = _parse_custom_force_parameters(path, name, values["parameters"])
    ij_count = _pairwise_ij_parameter_count(path, parameter_name)
    with_ele_value = values.get("with_ele", "true").strip()
    with_ele = not with_ele_value.lower() in {"false", "0"}
    return [
        TypedDataset("/forcefield/custom_force/pairwise/name", np.asarray(name, dtype=object)),
        TypedDataset("/forcefield/custom_force/pairwise/potential", np.asarray(values["potential"], dtype=object)),
        TypedDataset("/forcefield/custom_force/pairwise/parameters/text", np.asarray(values["parameters"], dtype=object)),
        TypedDataset("/forcefield/custom_force/pairwise/parameters/type", np.asarray(parameter_type, dtype=object)),
        TypedDataset("/forcefield/custom_force/pairwise/parameters/name", np.asarray(parameter_name, dtype=object)),
        TypedDataset("/forcefield/custom_force/pairwise/parameters/ij_count", np.asarray(abs(ij_count), dtype=np.int64)),
        TypedDataset("/forcefield/custom_force/pairwise/with_ele", np.asarray(with_ele, dtype=np.bool_)),
        TypedDataset(
            "/forcefield/custom_force/pairwise/electrostatic_potential",
            np.asarray(values.get("electrostatic_potential", ""), dtype=object),
        ),
    ]


def _pairwise_ij_parameter_count(path: Path, parameter_names: list[str]) -> int:
    ij_count = 0
    for parameter in parameter_names:
        if parameter.endswith("_ij"):
            if ij_count < 0:
                raise ValueError(f"{path} pairwise force parameters must put *_ij parameters first")
            ij_count += 1
        else:
            ij_count = -ij_count
    return abs(ij_count)


def _parse_listed_forces_config(path: Path) -> list[TypedDataset]:
    sections = _read_configuration_sections(path)
    if not sections:
        raise ValueError(f"{path} does not contain any listed force sections")

    names = []
    potentials = []
    parameter_texts = []
    parameter_types = []
    parameter_names = []
    parameter_offsets = [0]
    parameter_is_atom = []
    parameter_is_int = []
    connected_atoms = []
    constrain_distance = []
    for name, values in sections:
        _require_config_keys(path, name, values, {"potential", "parameters"})
        unknown = set(values) - {"potential", "parameters", "connected_atoms", "constrain_distance"}
        if unknown:
            raise ValueError(f"{path} listed force section {name} has unsupported keys: {sorted(unknown)}")
        types, params = _parse_custom_force_parameters(path, name, values["parameters"])
        atom_flags = []
        int_flags = []
        atom_count = 0
        for parameter_type, parameter_name in zip(types, params):
            is_int = parameter_type == "int"
            is_atom = is_int and parameter_name.startswith("atom_") and len(parameter_name) == 6
            if is_atom:
                atom_count += 1
            atom_flags.append(int(is_atom))
            int_flags.append(int(is_int))
        if atom_count < 1 or atom_count > 6:
            raise ValueError(f"{path} listed force section {name} must declare 1 to 6 atom parameters")
        connected = values.get("connected_atoms", "")
        if connected and (len(connected) != 2 or connected[0] == connected[1]):
            raise ValueError(f"{path} listed force section {name} connected_atoms must contain two different labels")

        names.append(name)
        potentials.append(values["potential"])
        parameter_texts.append(values["parameters"])
        parameter_types.extend(types)
        parameter_names.extend(params)
        parameter_is_atom.extend(atom_flags)
        parameter_is_int.extend(int_flags)
        parameter_offsets.append(len(parameter_names))
        connected_atoms.append(connected)
        constrain_distance.append(values.get("constrain_distance", ""))

    return [
        TypedDataset("/forcefield/custom_force/listed/name", np.asarray(names, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/potential", np.asarray(potentials, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/parameters/text", np.asarray(parameter_texts, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/parameters/type", np.asarray(parameter_types, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/parameters/name", np.asarray(parameter_names, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/parameters/offset", np.asarray(parameter_offsets, dtype=np.int64)),
        TypedDataset("/forcefield/custom_force/listed/parameters/is_atom", np.asarray(parameter_is_atom, dtype=np.int32)),
        TypedDataset("/forcefield/custom_force/listed/parameters/is_int", np.asarray(parameter_is_int, dtype=np.int32)),
        TypedDataset("/forcefield/custom_force/listed/connected_atoms", np.asarray(connected_atoms, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/constrain_distance", np.asarray(constrain_distance, dtype=object)),
        TypedDataset("/forcefield/custom_force/listed/count", np.asarray(len(names), dtype=np.int64)),
    ]


def _read_configuration_sections(path: Path) -> list[tuple[str, dict[str, str]]]:
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
            values[key] = "\n".join(value_lines)
        else:
            raise ValueError(f"{path} missing '[[ end ]]' for section {section}")
        sections.append((section, values))
    return sections


def _require_config_keys(path: Path, section: str, values: dict[str, str], required: set[str]) -> None:
    missing = required - set(values)
    if missing:
        raise ValueError(f"{path} section {section} missing required keys: {sorted(missing)}")


def _parse_custom_force_parameters(path: Path, section: str, text: str) -> tuple[list[str], list[str]]:
    parameter_types = []
    parameter_names = []
    for item in text.split(","):
        columns = item.strip().split()
        if len(columns) != 2:
            raise ValueError(f"{path} section {section} parameter {item!r} must be '<type> <name>'")
        if columns[0] not in {"int", "float"}:
            raise ValueError(f"{path} section {section} parameter type must be 'int' or 'float'")
        parameter_types.append(columns[0])
        parameter_names.append(columns[1])
    if not parameter_names:
        raise ValueError(f"{path} section {section} must declare at least one parameter")
    return parameter_types, parameter_names


def _safe_h5_name(name: str) -> str:
    safe = "".join(char if char.isalnum() or char == "_" else "_" for char in name)
    return safe or "unnamed"


def _parse_tersoff(path: Path) -> list[TypedDataset]:
    lines = _non_empty_raw_lines(path)
    if not lines:
        raise ValueError(f"{path} is empty")
    header = lines[0].split()
    if len(header) < 2:
        raise ValueError(f"{path} header must declare atom and type counts")
    atom_count, type_count = int(header[0]), int(header[1])

    cursor = 1
    type_names: list[str] = []
    while cursor < len(lines) and len(type_names) < type_count:
        if lines[cursor].startswith("#"):
            cursor += 1
            continue
        type_names.extend(lines[cursor].split())
        cursor += 1
    if len(type_names) != type_count:
        raise ValueError(f"{path} declares {type_count} TERSOFF types but contains {len(type_names)} names")

    entry_names: list[list[str]] = []
    raw_params: list[list[float]] = []
    while cursor < len(lines):
        line = lines[cursor]
        columns = line.split()
        if line.startswith("#") or not columns:
            cursor += 1
            continue
        if _is_int_token(columns[0]):
            break
        if len(columns) != 17:
            raise ValueError(f"{path} TERSOFF row {len(entry_names)} has {len(columns)} columns, expected 17")
        entry_names.append(columns[:3])
        raw_params.append([float(value) for value in columns[3:]])
        cursor += 1
    if not entry_names:
        raise ValueError(f"{path} does not contain any TERSOFF parameter rows")

    atom_type_values = []
    while cursor < len(lines):
        line = lines[cursor]
        cursor += 1
        if line.startswith("#"):
            continue
        atom_type_values.extend(int(value) for value in line.split())
    if len(atom_type_values) != atom_count:
        raise ValueError(f"{path} declares {atom_count} atom types but contains {len(atom_type_values)}")
    atom_types = np.asarray(atom_type_values, dtype=np.int32)
    if np.any(atom_types < 0) or np.any(atom_types >= type_count):
        raise ValueError(f"{path} contains atom type outside [0, {type_count})")

    raw = np.asarray(raw_params, dtype=np.float32)
    native = np.zeros((raw.shape[0], 18), dtype=np.float32)
    native[:, :14] = raw
    native[:, 9] *= EV_TO_KCAL_MOL
    native[:, 13] *= EV_TO_KCAL_MOL
    for row_index, n_value in enumerate(raw[:, 6]):
        if n_value > 0:
            native[row_index, 14] = np.power(2.0 * n_value * 1.0e-16, -1.0 / n_value)
            native[row_index, 15] = np.power(2.0 * n_value * 1.0e-8, -1.0 / n_value)
            native[row_index, 16] = 1.0 / native[row_index, 15]
            native[row_index, 17] = 1.0 / native[row_index, 14]

    entry_name_array = np.asarray(entry_names, dtype=object)
    type_index = {name: index for index, name in enumerate(type_names)}
    type_triples = np.asarray(
        [[type_index[name] for name in names] for names in entry_names],
        dtype=np.int32,
    )
    map_values = np.full(type_count * type_count * type_count, -1, dtype=np.int32)
    for entry_index, triple in enumerate(type_triples):
        flat_index = int(triple[0]) * type_count * type_count + int(triple[1]) * type_count + int(triple[2])
        map_values[flat_index] = entry_index

    return [
        TypedDataset("/manybody/tersoff/type_name", np.asarray(type_names, dtype=object)),
        TypedDataset("/manybody/tersoff/entry/type_name", entry_name_array),
        TypedDataset("/manybody/tersoff/entry/type", type_triples),
        TypedDataset("/manybody/tersoff/entry/parameters_raw", raw),
        TypedDataset("/manybody/tersoff/entry/parameters", native),
        TypedDataset("/manybody/tersoff/map", map_values),
        TypedDataset("/manybody/tersoff/atom_type", atom_types),
        TypedDataset("/manybody/tersoff/atom_type_count", np.asarray(type_count, dtype=np.int32)),
        TypedDataset("/manybody/tersoff/entry/count", np.asarray(raw.shape[0], dtype=np.int64)),
    ]


def _is_int_token(value: str) -> bool:
    try:
        int(value)
    except ValueError:
        return False
    return True


def _is_float_token(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True
