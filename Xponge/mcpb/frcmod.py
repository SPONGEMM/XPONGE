"""Minimal frcmod artifact helpers for MCPB workflows."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

from ._compat import atom_element, atom_type_name

_ZN_EMPIRICAL_BOND_TABLES = {
    ("N", "Zn"): [
        (1.926, 124.58), (1.947, 113.59), (1.978, 93.97), (1.982, 90.76), (1.984, 90.80),
        (2.011, 78.18), (2.027, 70.85), (2.028, 70.86), (2.029, 66.61), (2.040, 66.69),
        (2.041, 66.11), (2.047, 62.61), (2.073, 52.77), (2.089, 50.10), (2.133, 36.30),
        (2.145, 35.80), (2.176, 29.24),
    ],
    ("O", "Zn"): [
        (1.860, 169.29), (1.986, 76.81), (2.011, 71.26), (2.054, 56.37), (2.109, 41.86), (2.112, 41.32),
    ],
    ("S", "Zn"): [
        (2.262, 88.50), (2.305, 67.39), (2.353, 50.99), (2.355, 51.79), (2.426, 32.69),
    ],
}


def _distance(atom1, atom2) -> float:
    dx = float(atom1.x) - float(atom2.x)
    dy = float(atom1.y) - float(atom2.y)
    dz = float(atom1.z) - float(atom2.z)
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def _angle(atom1, atom2, atom3) -> float:
    v1 = (
        float(atom1.x) - float(atom2.x),
        float(atom1.y) - float(atom2.y),
        float(atom1.z) - float(atom2.z),
    )
    v2 = (
        float(atom3.x) - float(atom2.x),
        float(atom3.y) - float(atom2.y),
        float(atom3.z) - float(atom2.z),
    )
    dot = sum(a * b for a, b in zip(v1, v2))
    n1 = math.sqrt(sum(a * a for a in v1)) or 1.0
    n2 = math.sqrt(sum(a * a for a in v2)) or 1.0
    cosine = max(-1.0, min(1.0, dot / n1 / n2))
    return math.degrees(math.acos(cosine))


def _sorted_type_pair(type1: str, type2: str) -> tuple[str, str]:
    return (type1, type2) if type1 <= type2 else (type2, type1)


def _sorted_element_pair(atom1, atom2) -> tuple[str, str]:
    return tuple(sorted((atom_element(atom1), atom_element(atom2))))


def _interpolate_force_constant(points: list[tuple[float, float]], distance: float) -> float:
    ordered = sorted(points)
    if distance <= ordered[0][0]:
        return round(float(ordered[0][1]), 1)
    if distance >= ordered[-1][0]:
        return round(float(ordered[-1][1]), 1)
    for (x1, y1), (x2, y2) in zip(ordered, ordered[1:]):
        if x1 <= distance <= x2:
            if x2 == x1:
                return round(float(max(y1, y2)), 1)
            ratio = (distance - x1) / (x2 - x1)
            return round(float(y1 + ratio * (y2 - y1)), 1)
    return round(float(ordered[-1][1]), 1)


def _blank_bond_terms(molecule, selection):
    bond_terms: dict[tuple[str, str], tuple[float, float]] = {}
    neighbors_by_ion: dict[int, list[int]] = {atom_id: [] for atom_id in selection.ion_atom_ids}
    for atom1, atom2 in selection.bonded_pairs:
        ion_id = atom1 if atom1 in neighbors_by_ion else atom2
        neighbor_id = atom2 if ion_id == atom1 else atom1
        neighbors_by_ion[ion_id].append(neighbor_id)
        type1 = atom_type_name(molecule.atoms[atom1])
        type2 = atom_type_name(molecule.atoms[atom2])
        bond_terms[_sorted_type_pair(type1, type2)] = (0.0, _distance(molecule.atoms[atom1], molecule.atoms[atom2]))
    return bond_terms, neighbors_by_ion


def _blank_angle_terms(molecule, neighbors_by_ion):
    angle_terms: dict[tuple[str, str, str], tuple[float, float]] = {}
    for ion_id, neighbors in neighbors_by_ion.items():
        if len(neighbors) < 2:
            continue
        ion_type = atom_type_name(molecule.atoms[ion_id])
        for idx, atom1 in enumerate(neighbors):
            for atom3 in neighbors[idx + 1:]:
                type1 = atom_type_name(molecule.atoms[atom1])
                type3 = atom_type_name(molecule.atoms[atom3])
                key = (type1, ion_type, type3)
                rev = (type3, ion_type, type1)
                angle_terms[min(key, rev)] = (0.0, _angle(molecule.atoms[atom1], molecule.atoms[ion_id], molecule.atoms[atom3]))
    return angle_terms


def build_blank_frcmod_text(request, selection) -> str:
    molecule = request.molecule
    bond_terms, neighbors_by_ion = _blank_bond_terms(molecule, selection)
    angle_terms = _blank_angle_terms(molecule, neighbors_by_ion)

    lines = [
        "Xponge MCPB blank frcmod",
        "MASS",
        "",
        "BOND",
    ]
    for (type1, type2), (force_constant, length) in sorted(bond_terms.items()):
        lines.append(f"{type1:>2}-{type2:<2} {force_constant:5.1f}    {length:7.4f}")
    lines.extend(["", "ANGL"])
    for (type1, type2, type3), (force_constant, theta) in sorted(angle_terms.items()):
        lines.append(f"{type1:>2}-{type2:<2}-{type3:<2} {force_constant:5.2f}    {theta:7.2f}")
    lines.extend(["", "DIHE", "", "IMPR", "", "NONBON", ""])
    return "\n".join(lines)


def build_empirical_frcmod_text(request, selection) -> str:
    molecule = request.molecule
    bond_terms: dict[tuple[str, str], tuple[float, float]] = {}
    neighbors_by_ion: dict[int, list[int]] = {atom_id: [] for atom_id in selection.ion_atom_ids}
    for atom1, atom2 in selection.bonded_pairs:
        ion_id = atom1 if atom1 in neighbors_by_ion else atom2
        neighbor_id = atom2 if ion_id == atom1 else atom1
        neighbors_by_ion[ion_id].append(neighbor_id)
        pair = _sorted_element_pair(molecule.atoms[atom1], molecule.atoms[atom2])
        table = _ZN_EMPIRICAL_BOND_TABLES.get(pair)
        if table is None:
            raise ValueError(
                f"empirical MCPB currently supports only Zn-N, Zn-O, and Zn-S bonds; got {pair[0]}-{pair[1]}"
            )
        type1 = atom_type_name(molecule.atoms[atom1])
        type2 = atom_type_name(molecule.atoms[atom2])
        distance = _distance(molecule.atoms[atom1], molecule.atoms[atom2])
        bond_terms[_sorted_type_pair(type1, type2)] = (_interpolate_force_constant(table, distance), distance)

    angle_terms: dict[tuple[str, str, str], tuple[float, float]] = {}
    for ion_id, neighbors in neighbors_by_ion.items():
        if len(neighbors) < 2:
            continue
        ion_type = atom_type_name(molecule.atoms[ion_id])
        for idx, atom1 in enumerate(neighbors):
            for atom3 in neighbors[idx + 1:]:
                type1 = atom_type_name(molecule.atoms[atom1])
                type3 = atom_type_name(molecule.atoms[atom3])
                key = (type1, ion_type, type3)
                rev = (type3, ion_type, type1)
                pair1 = _sorted_element_pair(molecule.atoms[atom1], molecule.atoms[ion_id])
                pair2 = _sorted_element_pair(molecule.atoms[atom3], molecule.atoms[ion_id])
                if pair1 == ("S", "Zn") or pair2 == ("S", "Zn"):
                    force_constant = 70.0
                elif pair1 in {("N", "Zn"), ("O", "Zn")} or pair2 in {("N", "Zn"), ("O", "Zn")}:
                    force_constant = 50.0
                else:
                    force_constant = 35.0
                angle_terms[min(key, rev)] = (force_constant, _angle(molecule.atoms[atom1], molecule.atoms[ion_id], molecule.atoms[atom3]))

    lines = [
        "Xponge MCPB empirical frcmod",
        "MASS",
        "",
        "BOND",
    ]
    for (type1, type2), (force_constant, length) in sorted(bond_terms.items()):
        lines.append(f"{type1:>2}-{type2:<2} {force_constant:5.1f}    {length:7.4f}")
    lines.extend(["", "ANGL"])
    for (type1, type2, type3), (force_constant, theta) in sorted(angle_terms.items()):
        lines.append(f"{type1:>2}-{type2:<2}-{type3:<2} {force_constant:5.2f}    {theta:7.2f}")
    lines.extend(["", "DIHE", "", "IMPR", "", "NONBON", ""])
    return "\n".join(lines)


def register_blank_parameters(request, selection) -> None:
    from .. import register_amber_angle_parameter, register_amber_bond_parameter

    bond_terms, neighbors_by_ion = _blank_bond_terms(request.molecule, selection)
    angle_terms = _blank_angle_terms(request.molecule, neighbors_by_ion)
    for (type1, type2), (force_constant, length) in bond_terms.items():
        register_amber_bond_parameter(type1, type2, force_constant, length)
    for atom_types, (force_constant, theta) in angle_terms.items():
        register_amber_angle_parameter(atom_types, force_constant, theta)


def register_empirical_parameters(request, selection) -> None:
    from .. import register_amber_angle_parameter, register_amber_bond_parameter

    molecule = request.molecule
    neighbors_by_ion: dict[int, list[int]] = {atom_id: [] for atom_id in selection.ion_atom_ids}
    for atom1, atom2 in selection.bonded_pairs:
        ion_id = atom1 if atom1 in neighbors_by_ion else atom2
        neighbor_id = atom2 if ion_id == atom1 else atom1
        neighbors_by_ion[ion_id].append(neighbor_id)
        pair = _sorted_element_pair(molecule.atoms[atom1], molecule.atoms[atom2])
        table = _ZN_EMPIRICAL_BOND_TABLES.get(pair)
        if table is None:
            raise ValueError(
                f"empirical MCPB currently supports only Zn-N, Zn-O, and Zn-S bonds; got {pair[0]}-{pair[1]}"
            )
        type1 = atom_type_name(molecule.atoms[atom1])
        type2 = atom_type_name(molecule.atoms[atom2])
        distance = _distance(molecule.atoms[atom1], molecule.atoms[atom2])
        register_amber_bond_parameter(type1, type2, _interpolate_force_constant(table, distance), distance)

    for ion_id, neighbors in neighbors_by_ion.items():
        if len(neighbors) < 2:
            continue
        ion_type = atom_type_name(molecule.atoms[ion_id])
        for idx, atom1 in enumerate(neighbors):
            for atom3 in neighbors[idx + 1:]:
                type1 = atom_type_name(molecule.atoms[atom1])
                type3 = atom_type_name(molecule.atoms[atom3])
                pair1 = _sorted_element_pair(molecule.atoms[atom1], molecule.atoms[ion_id])
                pair2 = _sorted_element_pair(molecule.atoms[atom3], molecule.atoms[ion_id])
                if pair1 == ("S", "Zn") or pair2 == ("S", "Zn"):
                    force_constant = 70.0
                elif pair1 in {("N", "Zn"), ("O", "Zn")} or pair2 in {("N", "Zn"), ("O", "Zn")}:
                    force_constant = 50.0
                else:
                    force_constant = 35.0
                register_amber_angle_parameter(
                    (type1, ion_type, type3),
                    force_constant,
                    _angle(molecule.atoms[atom1], molecule.atoms[ion_id], molecule.atoms[atom3]),
                )


def write_blank_frcmod_artifact(request, selection, directory: str | Path | None = None) -> str:
    if directory is None:
        directory = Path(tempfile.mkdtemp(prefix="xponge_mcpb_frcmod_"))
    else:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
    path = directory / "metal_center_blank.frcmod"
    path.write_text(build_blank_frcmod_text(request, selection), encoding="utf-8")
    return str(path)


def write_empirical_frcmod_artifact(request, selection, directory: str | Path | None = None) -> str:
    if directory is None:
        directory = Path(tempfile.mkdtemp(prefix="xponge_mcpb_frcmod_"))
    else:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
    path = directory / "metal_center_empirical.frcmod"
    path.write_text(build_empirical_frcmod_text(request, selection), encoding="utf-8")
    return str(path)
