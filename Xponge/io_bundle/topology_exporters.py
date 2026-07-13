"""Typed topology bundle to legacy SPONGE payload exporters."""

from __future__ import annotations

import numpy as np

from .errors import BundleExportError
from .legacy_materializer import LegacyPayload
from .topology_parsers import EV_TO_KCAL_MOL


def export_mass(contract, reader, context) -> list[LegacyPayload]:
    return [_counted_vector(contract.legacy_keys[0], reader.read(contract.bundle_file, "/atoms/mass"))]


def export_charge(contract, reader, context) -> list[LegacyPayload]:
    return [_counted_vector(contract.legacy_keys[0], reader.read(contract.bundle_file, "/atoms/charge"))]


def export_residue(contract, reader, context) -> list[LegacyPayload]:
    offsets = _vector(reader.read(contract.bundle_file, "/residues/atom_offset"), np.int64)
    if offsets.size < 2 or offsets[0] != 0 or np.any(np.diff(offsets) < 0):
        raise BundleExportError("residue atom offsets are invalid")
    atom_count = int(offsets[-1])
    counts = np.diff(offsets)
    lines = [f"{atom_count} {counts.size}"]
    lines.extend(str(int(value)) for value in counts)
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_exclude(contract, reader, context) -> list[LegacyPayload]:
    offsets = _vector(reader.read(contract.bundle_file, "/topology/exclusions/offset"), np.int64)
    values = _vector(reader.read(contract.bundle_file, "/topology/exclusions/list"), np.int64)
    if offsets.size < 1 or offsets[0] != 0 or offsets[-1] != values.size or np.any(np.diff(offsets) < 0):
        raise BundleExportError("exclusion offsets and values are inconsistent")
    atom_count = offsets.size - 1
    lines = [f"{atom_count} {values.size}"]
    for index in range(atom_count):
        row = values[offsets[index] : offsets[index + 1]]
        suffix = " ".join(str(int(value)) for value in row)
        lines.append(f"{row.size} {suffix}".rstrip())
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_bond(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        ("/forcefield/bond/atoms", "/forcefield/bond/k", "/forcefield/bond/r0"),
    )


def export_angle(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        ("/forcefield/angle/atoms", "/forcefield/angle/k", "/forcefield/angle/theta0"),
    )


def export_dihedral(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        (
            "/forcefield/dihedral/atoms",
            "/forcefield/dihedral/periodicity",
            "/forcefield/dihedral/k",
            "/forcefield/dihedral/phi0",
        ),
    )


def export_improper(contract, reader, context) -> list[LegacyPayload]:
    parameter_path = "/forcefield/improper/pk"
    if not reader.contains(contract.bundle_file, parameter_path):
        parameter_path = "/forcefield/improper/k"
    return _fixed_table(
        contract,
        reader,
        ("/forcefield/improper/atoms", parameter_path, "/forcefield/improper/phi0"),
    )


def export_nb14(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        ("/forcefield/nb14/atoms", "/forcefield/nb14/params"),
    )


def export_nb14_extra(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        ("/forcefield/nb14_extra/atoms", "/forcefield/nb14_extra/params"),
    )


def export_urey_bradley(contract, reader, context) -> list[LegacyPayload]:
    return _fixed_table(
        contract,
        reader,
        (
            "/forcefield/urey_bradley/atoms",
            "/forcefield/urey_bradley/angle_k",
            "/forcefield/urey_bradley/angle_theta0",
            "/forcefield/urey_bradley/bond_k",
            "/forcefield/urey_bradley/bond_r0",
        ),
    )


def export_gb(contract, reader, context) -> list[LegacyPayload]:
    values = np.asarray(reader.read(contract.bundle_file, "/forcefield/gb/params"))
    if values.ndim != 2 or values.shape[1] != 2:
        raise BundleExportError(f"GB parameters have shape {values.shape}, expected (atoms, 2)")
    lines = [str(values.shape[0])]
    lines.extend(" ".join(_format_float(value) for value in row) for row in values)
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_subsys_division(contract, reader, context) -> list[LegacyPayload]:
    return [
        _counted_vector(
            contract.legacy_keys[0],
            reader.read(contract.bundle_file, "/forcefield/subsys_division"),
        )
    ]


def export_cmap(contract, reader, context) -> list[LegacyPayload]:
    atoms = np.asarray(reader.read(contract.bundle_file, "/forcefield/cmap/atoms"), dtype=np.int64)
    cmap_type = _vector(reader.read(contract.bundle_file, "/forcefield/cmap/type"), np.int64)
    resolution = _vector(reader.read(contract.bundle_file, "/forcefield/cmap/resolution"), np.int64)
    grid = _vector(reader.read(contract.bundle_file, "/forcefield/cmap/grid_value"), np.float64)
    if atoms.ndim != 2 or atoms.shape[1] != 5 or atoms.shape[0] != cmap_type.size:
        raise BundleExportError("CMAP atom and type datasets are inconsistent")
    expected_grid = int(np.sum(resolution * resolution))
    if grid.size != expected_grid:
        raise BundleExportError(f"CMAP grid has {grid.size} values, expected {expected_grid}")
    lines = [f"{atoms.shape[0]} {resolution.size}", " ".join(str(int(value)) for value in resolution)]
    offset = 0
    for width in resolution:
        size = int(width) * int(width)
        values = grid[offset : offset + size]
        offset += size
        for row in values.reshape(int(width), int(width)):
            lines.append(" ".join(_format_float(value) for value in row))
    for atom_row, type_value in zip(atoms, cmap_type):
        lines.append(" ".join(str(int(value)) for value in (*atom_row, type_value)))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_virtual_atom(contract, reader, context) -> list[LegacyPayload]:
    types = _vector(reader.read(contract.bundle_file, "/forcefield/virtual_atom/type"), np.int64)
    atoms = _vector(reader.read(contract.bundle_file, "/forcefield/virtual_atom/atom"), np.int64)
    from_offsets = _vector(
        reader.read(contract.bundle_file, "/forcefield/virtual_atom/from_offset"), np.int64
    )
    from_values = _vector(reader.read(contract.bundle_file, "/forcefield/virtual_atom/from"), np.int64)
    parameter_offsets = _vector(
        reader.read(contract.bundle_file, "/forcefield/virtual_atom/parameter_offset"), np.int64
    )
    parameters = _vector(
        reader.read(contract.bundle_file, "/forcefield/virtual_atom/parameter"), np.float64
    )
    count = types.size
    if atoms.size != count or from_offsets.size != count + 1 or parameter_offsets.size != count + 1:
        raise BundleExportError("virtual atom datasets have inconsistent lengths")
    if from_offsets[0] != 0 or from_offsets[-1] != from_values.size:
        raise BundleExportError("virtual atom source offsets are invalid")
    if parameter_offsets[0] != 0 or parameter_offsets[-1] != parameters.size:
        raise BundleExportError("virtual atom parameter offsets are invalid")
    lines = []
    for index in range(count):
        source = from_values[from_offsets[index] : from_offsets[index + 1]]
        parameter = parameters[parameter_offsets[index] : parameter_offsets[index + 1]]
        row = [str(int(types[index])), str(int(atoms[index]))]
        row.extend(str(int(value)) for value in source)
        row.extend(_format_float(value) for value in parameter)
        lines.append(" ".join(row))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_lj(contract, reader, context) -> list[LegacyPayload]:
    type_count = int(reader.read_scalar(contract.bundle_file, "/forcefield/lj/atom_type_count"))
    pair_a = _vector(reader.read(contract.bundle_file, "/forcefield/lj/pair_A_12"), np.float64)
    pair_b = _vector(reader.read(contract.bundle_file, "/forcefield/lj/pair_B_6"), np.float64)
    atom_type = _vector(reader.read(contract.bundle_file, "/forcefield/lj/type"), np.int64)
    pair_count = type_count * (type_count + 1) // 2
    if pair_a.size != pair_count or pair_b.size != pair_count:
        raise BundleExportError(
            f"LJ type count {type_count} requires {pair_count} pair values, got "
            f"{pair_a.size} and {pair_b.size}"
        )

    lines = [f"{atom_type.size} {type_count}", ""]
    offset = 0
    for width in range(1, type_count + 1):
        lines.append(" ".join(_format_float(value) for value in pair_a[offset : offset + width]))
        offset += width
    lines.append("")
    offset = 0
    for width in range(1, type_count + 1):
        lines.append(" ".join(_format_float(value) for value in pair_b[offset : offset + width]))
        offset += width
    lines.append("")
    lines.extend(str(int(value)) for value in atom_type)
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_lj_soft_core(contract, reader, context) -> list[LegacyPayload]:
    type_count_a = int(
        reader.read_scalar(contract.bundle_file, "/forcefield/lj_soft_core/atom_type_count_A")
    )
    type_count_b = int(
        reader.read_scalar(contract.bundle_file, "/forcefield/lj_soft_core/atom_type_count_B")
    )
    pair_aa = _vector(reader.read(contract.bundle_file, "/forcefield/lj_soft_core/pair_AA"), np.float64)
    pair_ab = _vector(reader.read(contract.bundle_file, "/forcefield/lj_soft_core/pair_AB"), np.float64)
    pair_ba = _vector(reader.read(contract.bundle_file, "/forcefield/lj_soft_core/pair_BA"), np.float64)
    pair_bb = _vector(reader.read(contract.bundle_file, "/forcefield/lj_soft_core/pair_BB"), np.float64)
    atom_type_a = _vector(
        reader.read(contract.bundle_file, "/forcefield/lj_soft_core/atom_type_A"), np.int64
    )
    atom_type_b = _vector(
        reader.read(contract.bundle_file, "/forcefield/lj_soft_core/atom_type_B"), np.int64
    )
    if atom_type_a.size != atom_type_b.size:
        raise BundleExportError("LJ soft-core atom type lengths differ")
    expected_a = type_count_a * (type_count_a + 1) // 2
    expected_b = type_count_b * (type_count_b + 1) // 2
    if pair_aa.size != expected_a or pair_ab.size != expected_a:
        raise BundleExportError("LJ soft-core A pair table size is invalid")
    if pair_ba.size != expected_b or pair_bb.size != expected_b:
        raise BundleExportError("LJ soft-core B pair table size is invalid")
    lines = [f"{atom_type_a.size} {type_count_a} {type_count_b}"]
    lines.extend(_format_float(value) for values in (pair_aa, pair_ab, pair_ba, pair_bb) for value in values)
    lines.extend(f"{int(left)} {int(right)}" for left, right in zip(atom_type_a, atom_type_b))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_manybody(contract, reader, context) -> list[LegacyPayload]:
    root = contract.bundle_path
    type_count = int(reader.read_scalar(contract.bundle_file, root + "/atom_type_count"))
    pair_type = np.asarray(reader.read(contract.bundle_file, root + "/pair/type"), dtype=np.int64)
    pair_parameters = np.asarray(reader.read(contract.bundle_file, root + "/pair/parameters"))
    triple_type = np.asarray(reader.read(contract.bundle_file, root + "/triple/type"), dtype=np.int64)
    triple_parameters = np.asarray(reader.read(contract.bundle_file, root + "/triple/parameters"))
    atom_type = _vector(reader.read(contract.bundle_file, root + "/atom_type"), np.int64)
    if pair_type.shape != (type_count * type_count, 2) or pair_parameters.shape[0] != pair_type.shape[0]:
        raise BundleExportError(f"{contract.contract_id} pair table is inconsistent")
    if triple_type.shape != (type_count ** 3, 3) or triple_parameters.shape[0] != triple_type.shape[0]:
        raise BundleExportError(f"{contract.contract_id} triple table is inconsistent")
    lines = [f"{atom_type.size} {type_count}", "# pair"]
    for indices, parameters in zip(pair_type, pair_parameters):
        lines.append(
            " ".join(
                [*(str(int(value)) for value in indices), *(_format_float(value) for value in parameters)]
            )
        )
    lines.append("# triple")
    for indices, parameters in zip(triple_type, triple_parameters):
        lines.append(
            " ".join(
                [*(str(int(value)) for value in indices), *(_format_float(value) for value in parameters)]
            )
        )
    lines.append("# atom types")
    lines.extend(str(int(value)) for value in atom_type)
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_eam_atom_type(contract, reader, context) -> list[LegacyPayload]:
    values = _vector(reader.read(contract.bundle_file, "/manybody/eam/atom_type"), np.int64)
    return [
        LegacyPayload(
            contract.legacy_keys[0],
            "\n".join(str(int(value)) for value in values) + "\n",
        )
    ]


def export_eam(contract, reader, context) -> list[LegacyPayload]:
    root = "/manybody/eam"
    format_name = reader.read_text(contract.bundle_file, root + "/format")
    type_count = int(reader.read_scalar(contract.bundle_file, root + "/atom_type_count"))
    nrho = int(reader.read_scalar(contract.bundle_file, root + "/nrho"))
    drho = float(reader.read_scalar(contract.bundle_file, root + "/drho"))
    nr = int(reader.read_scalar(contract.bundle_file, root + "/nr"))
    dr = float(reader.read_scalar(contract.bundle_file, root + "/dr"))
    cut = float(reader.read_scalar(contract.bundle_file, root + "/cut"))
    embed = np.asarray(reader.read(contract.bundle_file, root + "/embed/raw_ev"))
    density = np.asarray(reader.read(contract.bundle_file, root + "/electron_density/value"))
    if embed.shape != (type_count, nrho) or density.shape != (type_count, nr):
        raise BundleExportError("EAM embedding or density table shape is inconsistent")

    if format_name == "funcfl":
        if type_count != 1:
            raise BundleExportError("funcfl EAM must contain exactly one atom type")
        atomic_number = int(reader.read(contract.bundle_file, root + "/atomic_number").reshape(-1)[0])
        mass = float(reader.read(contract.bundle_file, root + "/mass").reshape(-1)[0])
        lattice_constant = float(
            reader.read(contract.bundle_file, root + "/lattice_constant").reshape(-1)[0]
        )
        lattice_type = _text_vector(reader.read(contract.bundle_file, root + "/lattice_type"))[0]
        z_values = _vector(reader.read(contract.bundle_file, root + "/funcfl/z"), np.float64)
        if z_values.size != nr:
            raise BundleExportError("funcfl EAM z table has an invalid length")
        lines = [
            "Generated by XPONGE bundled input export",
            f"{atomic_number} {_format_float(mass)} {_format_float(lattice_constant)} {lattice_type}",
            f"{nrho} {_format_float(drho)} {nr} {_format_float(dr)} {_format_float(cut)}",
        ]
        lines.extend(_format_float(value) for value in embed[0])
        lines.extend(_format_float(value) for value in z_values)
        lines.extend(_format_float(value) for value in density[0])
    elif format_name == "setfl":
        type_names = _text_vector(reader.read(contract.bundle_file, root + "/type_name"))
        atomic_numbers = _vector(reader.read(contract.bundle_file, root + "/atomic_number"), np.int64)
        masses = _vector(reader.read(contract.bundle_file, root + "/mass"), np.float64)
        lattice_constants = _vector(
            reader.read(contract.bundle_file, root + "/lattice_constant"), np.float64
        )
        lattice_types = _text_vector(reader.read(contract.bundle_file, root + "/lattice_type"))
        if not all(
            len(values) == type_count
            for values in (type_names, atomic_numbers, masses, lattice_constants, lattice_types)
        ):
            raise BundleExportError("setfl EAM type metadata lengths differ")
        pair = np.asarray(reader.read(contract.bundle_file, root + "/pair_potential/value"))
        if pair.shape != (type_count, type_count, nr):
            raise BundleExportError("setfl EAM pair potential shape is invalid")
        lines = [
            "Generated by XPONGE bundled input export",
            "",
            "",
            f"{type_count} " + " ".join(type_names),
            f"{nrho} {_format_float(drho)} {nr} {_format_float(dr)} {_format_float(cut)}",
        ]
        for index in range(type_count):
            lines.append(
                f"{int(atomic_numbers[index])} {_format_float(masses[index])} "
                f"{_format_float(lattice_constants[index])} {lattice_types[index]}"
            )
            lines.extend(_format_float(value) for value in embed[index])
            lines.extend(_format_float(value) for value in density[index])
        for left in range(type_count):
            for right in range(left + 1):
                for index, value in enumerate(pair[left, right]):
                    radius = index * dr if index else 1.0e-8
                    raw = float(value) * radius / EV_TO_KCAL_MOL
                    lines.append(_format_float(raw))
    else:
        raise BundleExportError(f"unsupported EAM format {format_name!r}")
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_tersoff(contract, reader, context) -> list[LegacyPayload]:
    root = "/manybody/tersoff"
    type_names = _text_vector(reader.read(contract.bundle_file, root + "/type_name"))
    entry_names = np.asarray(reader.read(contract.bundle_file, root + "/entry/type_name"))
    raw = np.asarray(reader.read(contract.bundle_file, root + "/entry/parameters_raw"))
    atom_type = _vector(reader.read(contract.bundle_file, root + "/atom_type"), np.int64)
    if entry_names.ndim != 2 or entry_names.shape[1] != 3 or raw.shape != (entry_names.shape[0], 14):
        raise BundleExportError("TERSOFF entry names and raw parameters are inconsistent")
    lines = [f"{atom_type.size} {len(type_names)}", " ".join(type_names)]
    for names, parameters in zip(entry_names, raw):
        rendered_names = [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in names]
        lines.append(
            " ".join([*rendered_names, *(_format_float(value) for value in parameters)])
        )
    lines.append("# Atom types")
    lines.append(" ".join(str(int(value)) for value in atom_type))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_qc_type(contract, reader, context) -> list[LegacyPayload]:
    atom_index = _vector(reader.read(contract.bundle_file, "/qc/type/atom_index"), np.int64)
    symbols = _text_vector(reader.read(contract.bundle_file, "/qc/type/symbol"))
    charge = int(reader.read_scalar(contract.bundle_file, "/qc/type/charge"))
    multiplicity = int(reader.read_scalar(contract.bundle_file, "/qc/type/multiplicity"))
    if atom_index.size != len(symbols):
        raise BundleExportError("QC type atom index and symbol lengths differ")
    lines = [f"{atom_index.size} {charge} {multiplicity}"]
    lines.extend(f"{int(index)} {symbol}" for index, symbol in zip(atom_index, symbols))
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_reaxff_type(contract, reader, context) -> list[LegacyPayload]:
    names = _text_vector(reader.read(contract.bundle_file, "/manybody/reaxff/type/name"))
    return [LegacyPayload(contract.legacy_keys[0], f"{len(names)}\n" + "\n".join(names) + "\n")]


def export_reaxff_parameters(contract, reader, context) -> list[LegacyPayload]:
    bundle_file = contract.bundle_file
    root = "/manybody/reaxff/parameters"
    lines = [reader.read_text(bundle_file, root + "/header")]

    general_count = int(reader.read_scalar(bundle_file, root + "/general/count"))
    general_count_label = reader.read_text(bundle_file, root + "/general/count_label")
    general_values = _vector(reader.read(bundle_file, root + "/general/value"), np.float64)
    general_labels = _text_vector(reader.read(bundle_file, root + "/general/label"))
    if general_values.size != general_count or len(general_labels) != general_count:
        raise BundleExportError("REAXFF general parameter lengths differ")
    lines.append(_labeled_line([general_count], general_count_label))
    lines.extend(
        _labeled_line([value], label) for value, label in zip(general_values, general_labels)
    )

    atom_root = root + "/atom"
    atom_count = int(reader.read_scalar(bundle_file, atom_root + "/count"))
    lines.append(_labeled_line([atom_count], reader.read_text(bundle_file, atom_root + "/count_label")))
    lines.extend(_text_vector(reader.read(bundle_file, atom_root + "/header")))
    atom_names = _text_vector(reader.read(bundle_file, atom_root + "/type_name"))
    atom_values = _vector(reader.read(bundle_file, atom_root + "/value"), np.float64)
    atom_offsets = _vector(reader.read(bundle_file, atom_root + "/value_offset"), np.int64)
    atom_line_offsets = _vector(reader.read(bundle_file, atom_root + "/line_value_offset"), np.int64)
    atom_line_labels = _text_vector(reader.read(bundle_file, atom_root + "/line_label"))
    if len(atom_names) != atom_count or atom_offsets.size != atom_count + 1:
        raise BundleExportError("REAXFF atom metadata lengths differ")
    if atom_line_offsets.size != atom_count * 4 + 1 or len(atom_line_labels) != atom_count * 4:
        raise BundleExportError("REAXFF atom line offsets are invalid")
    for atom_index, name in enumerate(atom_names):
        for line_index in range(4):
            flat_index = atom_index * 4 + line_index
            values = atom_values[atom_line_offsets[flat_index] : atom_line_offsets[flat_index + 1]]
            prefix = [name] if line_index == 0 else []
            lines.append(_labeled_line([*prefix, *values], atom_line_labels[flat_index]))

    bond_root = root + "/bond"
    bond_count = int(reader.read_scalar(bundle_file, bond_root + "/count"))
    lines.append(_labeled_line([bond_count], reader.read_text(bundle_file, bond_root + "/count_label")))
    lines.extend(_text_vector(reader.read(bundle_file, bond_root + "/header")))
    bond_type = np.asarray(reader.read(bundle_file, bond_root + "/type"), dtype=np.int64)
    bond_values = _vector(reader.read(bundle_file, bond_root + "/value"), np.float64)
    bond_line_offsets = _vector(reader.read(bundle_file, bond_root + "/line_value_offset"), np.int64)
    bond_line_labels = _text_vector(reader.read(bundle_file, bond_root + "/line_label"))
    if bond_type.shape != (bond_count, 2) or bond_line_offsets.size != bond_count * 2 + 1:
        raise BundleExportError("REAXFF bond metadata lengths differ")
    for bond_index, types in enumerate(bond_type):
        for line_index in range(2):
            flat_index = bond_index * 2 + line_index
            values = bond_values[bond_line_offsets[flat_index] : bond_line_offsets[flat_index + 1]]
            prefix = list(types) if line_index == 0 else []
            lines.append(_labeled_line([*prefix, *values], bond_line_labels[flat_index]))

    for name, index_columns in (
        ("off_diagonal", 2),
        ("angle", 3),
        ("torsion", 4),
        ("hydrogen_bond", 3),
    ):
        section_root = root + "/" + name
        count = int(reader.read_scalar(bundle_file, section_root + "/count"))
        label = reader.read_text(bundle_file, section_root + "/count_label")
        types = np.asarray(reader.read(bundle_file, section_root + "/type"), dtype=np.int64)
        values = np.asarray(reader.read(bundle_file, section_root + "/value"), dtype=np.float64)
        if types.shape != (count, index_columns) or values.shape[0] != count:
            raise BundleExportError(f"REAXFF {name} table shape is invalid")
        lines.append(_labeled_line([count], label))
        for type_row, value_row in zip(types, values):
            finite_values = value_row[np.isfinite(value_row)]
            lines.append(_labeled_line([*type_row, *finite_values], ""))

    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def export_pairwise_force(contract, reader, context) -> list[LegacyPayload]:
    root = "/forcefield/custom_force/pairwise"
    name = reader.read_text(contract.bundle_file, root + "/name")
    potential = reader.read_text(contract.bundle_file, root + "/potential")
    parameter_text = reader.read_text(contract.bundle_file, root + "/parameters/text")
    with_ele = bool(reader.read_scalar(contract.bundle_file, root + "/with_ele"))
    electrostatic = reader.read_text(contract.bundle_file, root + "/electrostatic_potential")
    config_values = {
        "potential": potential,
        "parameters": parameter_text,
        "with_ele": "true" if with_ele else "false",
    }
    if electrostatic:
        config_values["electrostatic_potential"] = electrostatic
    payloads = [
        LegacyPayload(
            contract.legacy_keys[0],
            _render_configuration_sections([(name, config_values)]),
        )
    ]

    data_root = root + "/data/" + _safe_h5_name(name)
    if reader.contains(contract.bundle_file, data_root + "/parameter/value"):
        type_count = int(reader.read_scalar(contract.bundle_file, data_root + "/type_count"))
        atom_type = _vector(reader.read(contract.bundle_file, data_root + "/atom_type"), np.int64)
        parameter_types = _text_vector(reader.read(contract.bundle_file, data_root + "/parameter/type"))
        int_values = np.asarray(reader.read(contract.bundle_file, data_root + "/parameter/int_value"))
        float_values = np.asarray(reader.read(contract.bundle_file, data_root + "/parameter/float_value"))
        lines = [f"{atom_type.size} {type_count}"]
        for index, parameter_type in enumerate(parameter_types):
            values = int_values[index] if parameter_type == "int" else float_values[index]
            lines.extend(_format_scalar(value) for value in values)
        lines.extend(str(int(value)) for value in atom_type)
        payloads.append(
            LegacyPayload(
                f"{name}_in_file",
                "\n".join(lines) + "\n",
                filename=f"{context.prefix}_{name}.txt",
            )
        )
    return payloads


def export_listed_forces(contract, reader, context) -> list[LegacyPayload]:
    root = "/forcefield/custom_force/listed"
    names = _text_vector(reader.read(contract.bundle_file, root + "/name"))
    potentials = _text_vector(reader.read(contract.bundle_file, root + "/potential"))
    parameter_texts = _text_vector(reader.read(contract.bundle_file, root + "/parameters/text"))
    connected = _text_vector(reader.read(contract.bundle_file, root + "/connected_atoms"))
    constrain = _text_vector(reader.read(contract.bundle_file, root + "/constrain_distance"))
    if not all(len(values) == len(names) for values in (potentials, parameter_texts, connected, constrain)):
        raise BundleExportError("listed force config dataset lengths differ")
    sections = []
    for name, potential, parameters, connected_atoms, constrain_distance in zip(
        names, potentials, parameter_texts, connected, constrain
    ):
        values = {"potential": potential, "parameters": parameters}
        if connected_atoms:
            values["connected_atoms"] = connected_atoms
        if constrain_distance:
            values["constrain_distance"] = constrain_distance
        sections.append((name, values))
    payloads = [LegacyPayload(contract.legacy_keys[0], _render_configuration_sections(sections))]
    for name in names:
        data_root = root + "/data/" + _safe_h5_name(name)
        if not reader.contains(contract.bundle_file, data_root + "/parameter/value"):
            continue
        parameter_types = _text_vector(reader.read(contract.bundle_file, data_root + "/parameter/type"))
        int_values = np.asarray(reader.read(contract.bundle_file, data_root + "/parameter/int_value"))
        float_values = np.asarray(reader.read(contract.bundle_file, data_root + "/parameter/float_value"))
        if int_values.shape != float_values.shape or int_values.shape[1] != len(parameter_types):
            raise BundleExportError(f"listed force {name} data shape is invalid")
        lines = [str(int_values.shape[0])]
        for row_index in range(int_values.shape[0]):
            row = [
                str(int(int_values[row_index, column_index]))
                if parameter_type == "int"
                else _format_float(float_values[row_index, column_index])
                for column_index, parameter_type in enumerate(parameter_types)
            ]
            lines.append(" ".join(row))
        payloads.append(
            LegacyPayload(
                f"{name}_in_file",
                "\n".join(lines) + "\n",
                filename=f"{context.prefix}_{name}.txt",
            )
        )
    return payloads


TOPOLOGY_EXPORTERS = {
    "mass_in_file": export_mass,
    "charge_in_file": export_charge,
    "residue_in_file": export_residue,
    "exclude_in_file": export_exclude,
    "bond_in_file": export_bond,
    "angle_in_file": export_angle,
    "dihedral_in_file": export_dihedral,
    "improper_dihedral_in_file": export_improper,
    "nb14_in_file": export_nb14,
    "nb14_extra_in_file": export_nb14_extra,
    "urey_bradley_in_file": export_urey_bradley,
    "cmap_in_file": export_cmap,
    "gb_in_file": export_gb,
    "virtual_atom_in_file": export_virtual_atom,
    "LJ_in_file": export_lj,
    "LJ_soft_core_in_file": export_lj_soft_core,
    "subsys_division_in_file": export_subsys_division,
    "EAM_in_file": export_eam,
    "EAM_atom_type_in_file": export_eam_atom_type,
    "SW_in_file": export_manybody,
    "EDIP_in_file": export_manybody,
    "TERSOFF_in_file": export_tersoff,
    "qc_type_in_file": export_qc_type,
    "REAXFF_in_file": export_reaxff_parameters,
    "REAXFF_type_in_file": export_reaxff_type,
    "pairwise_force_in_file": export_pairwise_force,
    "listed_forces_in_file": export_listed_forces,
}


def _counted_vector(key: str, values) -> LegacyPayload:
    values = np.asarray(values).reshape(-1)
    lines = [str(values.size)]
    lines.extend(_format_scalar(value) for value in values)
    return LegacyPayload(key, "\n".join(lines) + "\n")


def _fixed_table(contract, reader, paths: tuple[str, ...]) -> list[LegacyPayload]:
    columns = []
    row_count = None
    for path in paths:
        values = np.asarray(reader.read(contract.bundle_file, path))
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        elif values.ndim != 2:
            raise BundleExportError(f"{contract.contract_id} dataset {path} has shape {values.shape}")
        if row_count is None:
            row_count = values.shape[0]
        elif values.shape[0] != row_count:
            raise BundleExportError(f"{contract.contract_id} table column lengths differ")
        columns.append(values)
    table = np.concatenate(columns, axis=1)
    lines = [str(int(row_count or 0))]
    lines.extend(" ".join(_format_scalar(value) for value in row) for row in table)
    return [LegacyPayload(contract.legacy_keys[0], "\n".join(lines) + "\n")]


def _vector(values, dtype) -> np.ndarray:
    return np.asarray(values, dtype=dtype).reshape(-1)


def _text_vector(values) -> list[str]:
    rendered = []
    for value in np.asarray(values).reshape(-1):
        rendered.append(value.decode("utf-8") if isinstance(value, bytes) else str(value))
    return rendered


def _format_scalar(value) -> str:
    if np.issubdtype(np.asarray(value).dtype, np.integer):
        return str(int(value))
    return _format_float(value)


def _format_float(value) -> str:
    number = float(value)
    if not np.isfinite(number):
        raise BundleExportError(f"cannot export non-finite value {number}")
    return format(number, ".9g")


def _labeled_line(values, label: str) -> str:
    rendered = []
    for value in values:
        if isinstance(value, str):
            rendered.append(value)
        else:
            rendered.append(_format_scalar(value))
    line = " ".join(rendered)
    return f"{line} ! {label}" if label else line


def _render_configuration_sections(sections: list[tuple[str, dict[str, str]]]) -> str:
    lines = []
    for name, values in sections:
        lines.append(f"[[[ {name} ]]]")
        for key, value in values.items():
            lines.append(f"[[ {key} ]]")
            lines.extend(str(value).splitlines() or [""])
        lines.append("[[ end ]]")
    return "\n".join(lines) + "\n"


def _safe_h5_name(name: str) -> str:
    safe = "".join(char if char.isalnum() or char == "_" else "_" for char in name)
    return safe or "unnamed"
