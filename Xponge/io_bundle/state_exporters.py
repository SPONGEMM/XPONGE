"""Typed protocol and restart bundle to legacy SPONGE exporters."""

from __future__ import annotations

import numpy as np

from .errors import BundleExportError
from .legacy_materializer import LegacyPayload
from .topology_exporters import _format_float, _format_scalar, _text_vector, _vector


_CONFIG_ROOTS = {
    "cv_in_file": "/cv",
    "restrain_in_file": "/restraint",
    "restrain_cv_in_file": "/restraint/cv",
    "steer_cv_in_file": "/steer",
    "SITS_in_file": "/sits",
}


def export_config(contract, reader, context) -> list[LegacyPayload]:
    root = _CONFIG_ROOTS[contract.exporter_id]
    names = _text_vector(reader.read(contract.bundle_file, root + "/config/section/name"))
    offsets = _vector(
        reader.read(contract.bundle_file, root + "/config/section/key_offset"), np.int64
    )
    keys = _text_vector(reader.read(contract.bundle_file, root + "/config/key"))
    values = _text_vector(reader.read(contract.bundle_file, root + "/config/value"))
    if offsets.size != len(names) + 1 or offsets[0] != 0 or offsets[-1] != len(keys):
        raise BundleExportError(f"{contract.contract_id} config section offsets are invalid")
    if len(keys) != len(values):
        raise BundleExportError(f"{contract.contract_id} config key/value lengths differ")
    lines = []
    for section_index, name in enumerate(names):
        lines.extend([name, "{"])
        for item_index in range(offsets[section_index], offsets[section_index + 1]):
            lines.append(f"    {keys[item_index]} = {values[item_index]}")
        lines.append("}")
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_constraint(contract, reader, context) -> list[LegacyPayload]:
    atoms = np.asarray(reader.read(contract.bundle_file, "/constraint/default/pairs/atoms"), dtype=np.int64)
    r0 = _vector(reader.read(contract.bundle_file, "/constraint/default/pairs/r0"), np.float64)
    if atoms.ndim != 2 or atoms.shape[1] != 2 or atoms.shape[0] != r0.size:
        raise BundleExportError("constraint atom and r0 datasets are inconsistent")
    lines = [str(atoms.shape[0])]
    lines.extend(
        f"{int(row[0])} {int(row[1])} {_format_float(value)}"
        for row, value in zip(atoms, r0)
    )
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_flat_int(contract, reader, context) -> list[LegacyPayload]:
    values = _vector(reader.read(contract.bundle_file, contract.bundle_path), np.int64)
    return [
        LegacyPayload(
            contract.legacy_keys[0],
            "\n".join(str(int(value)) for value in values) + "\n",
        )
    ]


def export_flat_float(contract, reader, context) -> list[LegacyPayload]:
    values = _vector(reader.read(contract.bundle_file, contract.bundle_path), np.float64)
    return [
        LegacyPayload(
            contract.legacy_keys[0],
            " ".join(_format_float(value) for value in values) + "\n",
        )
    ]


def export_matrix(contract, reader, context) -> list[LegacyPayload]:
    values = np.asarray(reader.read(contract.bundle_file, contract.bundle_path))
    if values.ndim != 2:
        raise BundleExportError(f"{contract.contract_id} matrix has shape {values.shape}")
    lines = [" ".join(_format_scalar(value) for value in row) for row in values]
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_restraint_reference(contract, reader, context) -> list[LegacyPayload]:
    root = "/parameters/restart/references/restraint/default"
    coordinate = np.asarray(reader.read(contract.bundle_file, root + "/coordinate"))
    if coordinate.ndim != 2 or coordinate.shape[1] != 3:
        raise BundleExportError(f"restraint reference coordinate has shape {coordinate.shape}")
    if reader.contains(contract.bundle_file, root + "/source_title"):
        title = reader.read_text(contract.bundle_file, root + "/source_title")
        atom_count = int(reader.read_scalar(contract.bundle_file, root + "/atom_count"))
        time = float(reader.read_scalar(contract.bundle_file, root + "/time"))
        box = _vector(reader.read(contract.bundle_file, root + "/box"), np.float64)
        if atom_count != coordinate.shape[0] or box.size != 6:
            raise BundleExportError("Amber restraint rst7 metadata is inconsistent")
        lines = [title, f"{atom_count} {_format_float(time)}"]
        lines.extend(_format_xyz_rows(coordinate))
        if reader.contains(contract.bundle_file, root + "/velocity"):
            velocity = np.asarray(reader.read(contract.bundle_file, root + "/velocity"))
            if velocity.shape != coordinate.shape:
                raise BundleExportError("Amber restraint rst7 velocity shape differs from coordinate")
            lines.extend(_format_xyz_rows(velocity))
        lines.append(" ".join(_format_float(value) for value in box))
        return [LegacyPayload("restrain_amber_rst7", "\n".join(lines) + "\n")]

    lines = [str(coordinate.shape[0]), *_format_xyz_rows(coordinate)]
    return [LegacyPayload("restrain_coordinate_in_file", "\n".join(lines) + "\n")]


def export_soft_walls(contract, reader, context) -> list[LegacyPayload]:
    names = _text_vector(reader.read(contract.bundle_file, "/wall/soft/name"))
    potentials = _text_vector(reader.read(contract.bundle_file, "/wall/soft/potential"))
    if len(names) != len(potentials):
        raise BundleExportError("soft wall name and potential lengths differ")
    lines = []
    for name, potential in zip(names, potentials):
        lines.extend([f"[[[ {name} ]]]", "[[ potential ]]", potential, "[[ end ]]"])
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_meta_edge(contract, reader, context) -> list[LegacyPayload]:
    root = "/meta/default/grid"
    coordinate = np.asarray(reader.read(contract.bundle_file, root + "/coordinate"))
    normalization = _vector(reader.read(contract.bundle_file, root + "/normalization"), np.float64)
    force = np.asarray(reader.read(contract.bundle_file, root + "/force"))
    if coordinate.ndim != 2 or force.shape != coordinate.shape or normalization.size != coordinate.shape[0]:
        raise BundleExportError("meta edge datasets have inconsistent shapes")
    lines = []
    for coordinate_row, normalization_value, force_row in zip(coordinate, normalization, force):
        values = [*coordinate_row, normalization_value, *force_row]
        lines.append(" ".join(_format_float(value) for value in values))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_meta_state(contract, reader, context) -> list[LegacyPayload]:
    key = contract.legacy_keys[0]
    root = "/parameters/restart/bias/meta/default/potential"
    if key == "meta_scatter_in_file":
        root = "/parameters/restart/bias/meta/default/scatter"
    title = reader.read_text(contract.bundle_file, root + "/title")
    axis_min = _vector(reader.read(contract.bundle_file, root + "/axis/min"), np.float64)
    axis_max = _vector(reader.read(contract.bundle_file, root + "/axis/max"), np.float64)
    axis_delta = _vector(reader.read(contract.bundle_file, root + "/axis/delta"), np.float64)
    grid = _vector(reader.read(contract.bundle_file, root + "/grid"), np.int64)
    coordinate = np.asarray(reader.read(contract.bundle_file, root + "/coordinate"))
    value = _vector(reader.read(contract.bundle_file, root + "/value"), np.float64)
    force = np.asarray(reader.read(contract.bundle_file, root + "/force"))
    has_value = _vector(reader.read(contract.bundle_file, root + "/has_value"), np.int64)
    has_force = _vector(reader.read(contract.bundle_file, root + "/has_force"), np.int64)
    ndim = axis_min.size
    row_count = coordinate.shape[0]
    if not (
        axis_max.size == axis_delta.size == grid.size == ndim
        and coordinate.shape == (row_count, ndim)
        and force.shape == (row_count, ndim)
        and value.size == has_value.size == has_force.size == row_count
    ):
        raise BundleExportError(f"{contract.contract_id} meta state datasets are inconsistent")
    lines = [title]
    lines.extend(
        f"{_format_float(left)} {_format_float(right)} {_format_float(delta)}"
        for left, right, delta in zip(axis_min, axis_max, axis_delta)
    )
    lines.append(" ".join([*(str(int(item)) for item in grid), str(row_count)]))
    for index in range(row_count):
        row = [_format_float(item) for item in coordinate[index]]
        if has_force[index]:
            row.append("0")
            row.extend(_format_float(item) for item in force[index])
            row.append(_format_float(value[index]))
        elif has_value[index]:
            row.extend(["0", _format_float(value[index])])
        lines.append(" ".join(row))
    return [LegacyPayload(key, "\n".join(lines) + "\n")]


def export_meta_hills(contract, reader, context) -> list[LegacyPayload]:
    root = "/parameters/restart/bias/meta/default/hills_typed"
    values = np.asarray(reader.read(contract.bundle_file, root + "/value"))
    if values.ndim != 2:
        raise BundleExportError(f"meta hills matrix has shape {values.shape}")
    lines = [" ".join(_format_float(item) for item in row) for row in values]
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


STATE_EXPORTERS = {
    "cv_in_file": export_config,
    "restrain_in_file": export_config,
    "restrain_cv_in_file": export_config,
    "steer_cv_in_file": export_config,
    "SITS_in_file": export_config,
    "constrain_in_file": export_constraint,
    "restrain_atom_id": export_flat_int,
    "restrain_weight_in_file": export_matrix,
    "restraint_reference": export_restraint_reference,
    "nose_hoover_chain_restart_input": export_matrix,
    "SITS_atom_in_file": export_flat_int,
    "SITS_nk_in_file": export_flat_float,
    "soft_walls_in_file": export_soft_walls,
    "meta_edge_in_file": export_meta_edge,
    "meta_potential_in_file": export_meta_state,
    "meta_scatter_in_file": export_meta_state,
    "hills_in_file": export_meta_hills,
}


def _format_xyz_rows(values: np.ndarray) -> list[str]:
    return [" ".join(_format_float(value) for value in row) for row in values]
