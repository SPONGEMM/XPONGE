"""
Generate Xponge GLYCAM template mol2 files from Amber PREP units.
"""
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

from audit import EXTERNAL_COVERAGE, FUNCTIONAL_GROUP_TEMPLATES, _parse_mol2_units

_UNIT_RE = re.compile(r"^([A-Za-z0-9]{3,4})\s+INT\s+0\s*$")
_ATOM_TYPE_NORMALIZATION = {
    "HC": "Hc",
    "HP": "Hp",
}


@dataclass
class PrepAtom:
    index: int
    name: str
    atom_type: str
    parent: int
    angle_ref: int
    dihedral_ref: int
    bond_length: float
    angle_deg: float
    dihedral_deg: float
    charge: float


@dataclass
class PrepUnit:
    name: str
    net_charge: float
    atoms: list[PrepAtom]
    loops: list[tuple[str, str]]


def _normalized_ref(current_index: int, reference_index: int) -> int:
    """
    Amber PREP occasionally contains self/forward references in a few modified
    monosaccharide units (e.g. 4AE). For template export we normalize those
    rare records back to the previous atom index.
    """
    if current_index > 3 and reference_index >= current_index:
        return current_index - 1
    return reference_index


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _default_ambertools_root() -> Path:
    return Path("/mnt/data8t/Software/AmberTools26/ambertools26_src")


def parse_prep_units(prep_file: Path) -> dict[str, PrepUnit]:
    units: dict[str, PrepUnit] = {}
    lines = prep_file.read_text().splitlines()
    i = 0
    while i < len(lines):
        match = _UNIT_RE.match(lines[i].strip())
        if not match:
            i += 1
            continue
        name = match.group(1)
        i += 1
        if i >= len(lines):
            break
        i += 1  # CORRECT OMIT ...
        net_charge = float(lines[i].strip())
        i += 1
        atoms: list[PrepAtom] = []
        loops: list[tuple[str, str]] = []
        while i < len(lines):
            stripped = lines[i].strip()
            if not stripped:
                i += 1
                continue
            if stripped == "LOOP":
                i += 1
                while i < len(lines):
                    stripped = lines[i].strip()
                    if not stripped:
                        i += 1
                        continue
                    if stripped == "DONE":
                        i += 1
                        break
                    loops.append(tuple(stripped.split()[:2]))
                    i += 1
                break
            if stripped == "DONE":
                i += 1
                break
            fields = stripped.split()
            if len(fields) < 11 or not fields[0].isdigit():
                i += 1
                continue
            atoms.append(
                PrepAtom(
                    index=int(fields[0]),
                    name=fields[1],
                    atom_type=fields[2],
                    parent=int(fields[4]),
                    angle_ref=int(fields[5]),
                    dihedral_ref=int(fields[6]),
                    bond_length=float(fields[7]),
                    angle_deg=float(fields[8]),
                    dihedral_deg=float(fields[9]),
                    charge=float(fields[10]),
                )
            )
            i += 1
        units[name] = PrepUnit(name=name, net_charge=net_charge, atoms=atoms, loops=loops)
    return units


def _vec_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_scale(a, scalar):
    return (a[0] * scalar, a[1] * scalar, a[2] * scalar)


def _vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _vec_cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vec_norm(a):
    return math.sqrt(_vec_dot(a, a))


def _vec_unit(a):
    norm = _vec_norm(a)
    if norm < 1e-12:
        raise ValueError("zero-length vector")
    return _vec_scale(a, 1.0 / norm)


def _fallback_normal(axis):
    if abs(axis[0]) < 0.9:
        probe = (1.0, 0.0, 0.0)
    else:
        probe = (0.0, 1.0, 0.0)
    normal = _vec_cross(axis, probe)
    return _vec_unit(normal)


def _zmatrix_to_cartesian(unit: PrepUnit) -> dict[int, tuple[float, float, float]]:
    coords: dict[int, tuple[float, float, float]] = {}
    for atom in unit.atoms:
        if atom.index == 1:
            coords[atom.index] = (0.0, 0.0, 0.0)
            continue
        if atom.index == 2:
            coords[atom.index] = (atom.bond_length, 0.0, 0.0)
            continue
        if atom.index == 3:
            theta = math.radians(atom.angle_deg)
            parent = coords[_normalized_ref(atom.index, atom.parent)]
            coords[atom.index] = (
                parent[0] - atom.bond_length * math.cos(theta),
                parent[1] + atom.bond_length * math.sin(theta),
                0.0,
            )
            continue

        parent_index = _normalized_ref(atom.index, atom.parent)
        angle_ref_index = _normalized_ref(atom.index, atom.angle_ref)
        dihedral_ref_index = _normalized_ref(atom.index, atom.dihedral_ref)

        bond = coords[parent_index]
        angle = coords[angle_ref_index]
        dihedral = coords[dihedral_ref_index]

        e1 = _vec_unit(_vec_sub(bond, angle))
        plane_normal = _vec_cross(_vec_sub(angle, dihedral), e1)
        if _vec_norm(plane_normal) < 1e-8:
            plane_normal = _fallback_normal(e1)
        else:
            plane_normal = _vec_unit(plane_normal)
        e2 = _vec_cross(plane_normal, e1)

        theta = math.radians(atom.angle_deg)
        phi = math.radians(atom.dihedral_deg)
        local = _vec_add(
            _vec_scale(e1, -atom.bond_length * math.cos(theta)),
            _vec_add(
                _vec_scale(e2, atom.bond_length * math.sin(theta) * math.cos(phi)),
                _vec_scale(plane_normal, atom.bond_length * math.sin(theta) * math.sin(phi)),
            ),
        )
        coords[atom.index] = _vec_add(bond, local)
    return coords


def _mol2_atom_line(atom_id: int, atom: PrepAtom, coord, residue_id: int, residue_name: str) -> str:
    return (
        f"{atom_id:>6d} {atom.name:>4s} "
        f"{coord[0]:>9.4f} {coord[1]:>9.4f} {coord[2]:>9.4f} "
        f"{_ATOM_TYPE_NORMALIZATION.get(atom.atom_type, atom.atom_type):>5s} "
        f"{residue_id:>5d} {residue_name:>8s} {atom.charge:>10.6f}"
    )


def render_units_as_mol2(units: list[PrepUnit], title: str = "glycam_extended_templates") -> str:
    atom_lines: list[str] = []
    bond_lines: list[str] = []
    substructure_lines: list[str] = []
    atom_id = 1
    bond_id = 1
    residue_id = 1
    for unit in units:
        coords = _zmatrix_to_cartesian(unit)
        real_atoms = [atom for atom in unit.atoms if atom.index > 3]
        local_to_global: dict[int, int] = {}
        first_atom_id = atom_id
        for atom in real_atoms:
            local_to_global[atom.index] = atom_id
            atom_lines.append(_mol2_atom_line(atom_id, atom, coords[atom.index], residue_id, unit.name))
            atom_id += 1
        for atom in real_atoms:
            parent_index = _normalized_ref(atom.index, atom.parent)
            if parent_index > 3:
                bond_lines.append(f"{bond_id:>6d} {local_to_global[atom.index]:>6d} {local_to_global[parent_index]:>6d} 1")
                bond_id += 1
        name_to_index = {atom.name: atom.index for atom in real_atoms}
        seen_loop_bonds = set()
        for atom_a, atom_b in unit.loops:
            if atom_a not in name_to_index or atom_b not in name_to_index:
                continue
            pair = tuple(sorted((name_to_index[atom_a], name_to_index[atom_b])))
            if pair in seen_loop_bonds:
                continue
            seen_loop_bonds.add(pair)
            bond_lines.append(f"{bond_id:>6d} {local_to_global[pair[0]]:>6d} {local_to_global[pair[1]]:>6d} 1")
            bond_id += 1
        substructure_lines.append(f"{residue_id:>6d} {unit.name:>8s} {first_atom_id:>6d} ****               0 ****  ****")
        residue_id += 1

    header = [
        "@<TRIPOS>MOLECULE",
        title,
        f"{len(atom_lines)} {len(bond_lines)} {len(substructure_lines)} 0 1",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    return "\n".join(
        header
        + atom_lines
        + ["@<TRIPOS>BOND"]
        + bond_lines
        + ["@<TRIPOS>SUBSTRUCTURE"]
        + substructure_lines
        + [""]
    )


def gather_existing_xponge_units(repo_root: Path, target_name: str) -> set[str]:
    glycam_dir = repo_root / "Xponge" / "forcefield" / "amber" / "glycam_06j"
    units = set()
    for mol2_file in sorted(glycam_dir.glob("*.mol2")):
        if mol2_file.name == target_name:
            continue
        units |= _parse_mol2_units(mol2_file)
    return units


def build_default_unit_list(repo_root: Path, amber_root: Path, target_name: str) -> list[str]:
    prep_file = amber_root / "dat" / "leap" / "prep" / "GLYCAM_06j-1.prep"
    all_prep_units = parse_prep_units(prep_file)
    existing = gather_existing_xponge_units(repo_root, target_name)
    covered_elsewhere = set(EXTERNAL_COVERAGE) | {"HYP", "NHYP", "CHYP"}
    missing = [
        name for name in sorted(all_prep_units)
        if name not in existing and name not in covered_elsewhere and name not in FUNCTIONAL_GROUP_TEMPLATES
    ]
    return missing


def main():
    parser = argparse.ArgumentParser(description="Export missing Amber GLYCAM PREP units into a mol2 template file.")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=_default_repo_root(),
        help="Xponge repository root",
    )
    parser.add_argument(
        "--ambertools-root",
        type=Path,
        default=_default_ambertools_root(),
        help="AmberTools source root",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output mol2 path. Defaults to glycam_06j/modified_monosaccharides.mol2",
    )
    parser.add_argument(
        "--unit",
        action="append",
        default=None,
        help="Specific PREP unit name to export. May be repeated.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    amber_root = args.ambertools_root.resolve()
    output = args.output or repo_root / "Xponge" / "forcefield" / "amber" / "glycam_06j" / "modified_monosaccharides.mol2"

    prep_file = amber_root / "dat" / "leap" / "prep" / "GLYCAM_06j-1.prep"
    prep_units = parse_prep_units(prep_file)

    requested_units = args.unit or build_default_unit_list(repo_root, amber_root, output.name)
    units = [prep_units[name] for name in requested_units if name in prep_units]
    output.write_text(render_units_as_mol2(units))
    print(f"Wrote {len(units)} PREP units to {output}")


if __name__ == "__main__":
    main()
