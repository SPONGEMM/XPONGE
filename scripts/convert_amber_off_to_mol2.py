#!/usr/bin/env python3
"""Convert single-residue Amber OFF libraries into Xponge template MOL2 data."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path


_HEADER = re.compile(r"^!entry\.([^.]*)\.unit\.([^ ]+)")


@dataclass
class Atom:
    name: str
    atom_type: str
    atomic_number: int
    charge: float
    position: tuple[float, float, float] | None = None


@dataclass
class Template:
    name: str
    atoms: list[Atom] = field(default_factory=list)
    bonds: list[tuple[int, int, int]] = field(default_factory=list)
    connect: tuple[int, int] = (0, 0)
    unit_name: str | None = None
    residue_names: list[str] = field(default_factory=list)


def _section_rows(lines: list[str], start: int) -> tuple[list[str], int]:
    rows = []
    index = start
    while index < len(lines) and not lines[index].startswith("!"):
        if lines[index].strip():
            rows.append(lines[index])
        index += 1
    return rows, index


def parse_off(path: Path) -> list[Template]:
    """Parse the OFF fields required for deterministic Xponge templates."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or not lines[0].startswith("!!index array str"):
        raise ValueError(f"{path} is not an Amber OFF library")

    names = []
    index = 1
    while index < len(lines) and not lines[index].startswith("!"):
        if lines[index].strip():
            names.append(shlex.split(lines[index])[0])
        index += 1
    templates = {name: Template(name) for name in names}

    while index < len(lines):
        match = _HEADER.match(lines[index])
        if not match:
            index += 1
            continue
        name, section = match.groups()
        if name not in templates:
            raise ValueError(f"section for undeclared OFF entry {name!r}")
        rows, index = _section_rows(lines, index + 1)
        template = templates[name]

        if section == "atoms":
            for row in rows:
                words = shlex.split(row)
                template.atoms.append(
                    Atom(words[0], words[1], int(words[6]), float(words[7]))
                )
        elif section == "positions":
            if len(rows) != len(template.atoms):
                raise ValueError(f"{name}: atom/position count mismatch")
            for atom, row in zip(template.atoms, rows):
                xyz = tuple(float(value) for value in row.split())
                if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
                    raise ValueError(f"{name}: invalid atom position {row!r}")
                atom.position = xyz
        elif section == "connectivity":
            template.bonds = [tuple(int(value) for value in row.split()) for row in rows]
        elif section == "connect":
            if len(rows) != 2:
                raise ValueError(f"{name}: expected two unit.connect entries")
            template.connect = tuple(int(row.strip()) for row in rows)
        elif section == "name":
            template.unit_name = shlex.split(rows[0])[0]
        elif section == "residues":
            template.residue_names = [shlex.split(row)[0] for row in rows]

    result = []
    for name in names:
        template = templates[name]
        if not template.atoms or any(atom.position is None for atom in template.atoms):
            raise ValueError(f"{name}: incomplete atom data")
        if template.residue_names != [name]:
            raise ValueError(f"{name}: only one-residue units with matching residue names are supported")
        if not template.unit_name:
            raise ValueError(f"{name}: missing unit.name")
        for atom1, atom2, order in template.bonds:
            if not (1 <= atom1 <= len(template.atoms) and 1 <= atom2 <= len(template.atoms)):
                raise ValueError(f"{name}: bond atom index out of range")
            if order not in (1, 2, 3):
                raise ValueError(f"{name}: unsupported OFF connectivity flag {order}")
        for atom_index in template.connect:
            if not 0 <= atom_index <= len(template.atoms):
                raise ValueError(f"{name}: unit.connect index out of range")
        result.append(template)
    return result


def render_mol2(source: Path, templates: list[Template]) -> str:
    """Render all templates as substructures in one Xponge-compatible MOL2."""
    atom_count = sum(len(template.atoms) for template in templates)
    bond_count = sum(len(template.bonds) for template in templates)
    lines = [
        "@<TRIPOS>MOLECULE",
        source.stem.upper(),
        f"{atom_count:6d} {bond_count:6d} {len(templates):6d} 0 1",
        "SMALL",
        "USER_CHARGES",
        "@<TRIPOS>ATOM",
    ]
    atom_offset = 0
    for residue_index, template in enumerate(templates, 1):
        for local_index, atom in enumerate(template.atoms, 1):
            x, y, z = atom.position
            lines.append(
                f"{atom_offset + local_index:6d} {atom.name:<6s} "
                f"{x:12.6f} {y:12.6f} {z:12.6f} {atom.atom_type:<4s} "
                f"{residue_index:5d} {template.name:<8s} {atom.charge:12.6f}"
            )
        atom_offset += len(template.atoms)

    lines.append("@<TRIPOS>BOND")
    atom_offset = 0
    bond_index = 1
    for template in templates:
        for atom1, atom2, order in template.bonds:
            lines.append(
                f"{bond_index:6d} {atom_offset + atom1:6d} {atom_offset + atom2:6d} {order}"
            )
            bond_index += 1
        atom_offset += len(template.atoms)

    lines.append("@<TRIPOS>SUBSTRUCTURE")
    atom_offset = 0
    for residue_index, template in enumerate(templates, 1):
        lines.append(
            f"{residue_index:6d} {template.name:<8s} {atom_offset + 1:6d} "
            "****               0 ****  ****"
        )
        atom_offset += len(template.atoms)
    return "\n".join(lines) + "\n"


def _connection_metadata(template: Template, position: str, anchor: str | None) -> dict:
    """Derive Xponge orientation metadata for a validated Amber lipid connection."""
    if anchor is None:
        return {
            f"{position}_next_atom": None,
            f"{position}_reference_atom": None,
            f"{position}_link_conditions": [],
            f"{position}_rule_source": "none",
        }

    atom_names = {atom.name for atom in template.atoms}

    def result(next_atom, reference, angle_atoms, angle, dihedral_atoms, dihedral, source):
        required = {anchor, next_atom, reference, *angle_atoms, *dihedral_atoms}
        missing = required - atom_names
        if missing:
            raise ValueError(
                f"{template.name}: {position} connection rule references missing atoms {sorted(missing)}"
            )
        if len({anchor, next_atom, reference}) != 3:
            raise ValueError(
                f"{template.name}: {position} anchor/next/reference atoms must be distinct"
            )
        return {
            f"{position}_next_atom": next_atom,
            f"{position}_reference_atom": reference,
            f"{position}_link_conditions": [
                {"atoms": angle_atoms, "parameter_degrees": angle},
                {"atoms": dihedral_atoms, "parameter_degrees": dihedral},
            ],
            f"{position}_rule_source": source,
        }

    if template.name == "SPM" and position == "head":
        return result("N11", "O12", ["N11", "C11"], 120.0,
                      ["O12", "N11", "C11"], 180.0, "explicit_spm_amide")
    if template.name == "SPM" and position == "tail":
        return result("C2", "O22", ["C2", "C1"], 109.5,
                      ["O22", "C2", "C1"], -120.0, "explicit_spm_sphingosine")
    if template.name == "SA" and anchor == "C12":
        return result("C13", "H2R", ["H2R", "C12"], 120.0,
                      ["C13", "H2R", "C12"], 180.0, "explicit_sa_alkene")

    if "C11" in anchor:
        next_atom = anchor.replace("C11", "O11", 1)
        reference = anchor.replace("C11", "O12", 1)
        if next_atom in atom_names and reference in atom_names:
            return result(next_atom, reference, [next_atom, anchor], 120.0,
                          [reference, next_atom, anchor], 180.0, "carbonyl_head")
    if "C21" in anchor:
        next_atom = anchor.replace("C21", "O21", 1)
        reference = anchor.replace("C21", "O22", 1)
        if next_atom in atom_names and reference in atom_names:
            return result(next_atom, reference, [next_atom, anchor], 120.0,
                          [reference, next_atom, anchor], 180.0, "carbonyl_tail")
    if anchor == "C12" and {"C13", "H2R", "H2S"}.issubset(atom_names):
        return result("C13", "H2S", ["H2R", "C12"], 109.5,
                      ["H2S", "H2R", "C12"],
                      -120.0 if position == "head" else 120.0, "standard_chain")

    raise ValueError(f"{template.name}: cannot derive {position} connection geometry for {anchor}")


def build_manifest(source: Path, templates: list[Template], source_license: str) -> dict:
    entries = []
    for template in templates:
        charge = sum(atom.charge for atom in template.atoms)
        if abs(charge) < 5e-9:
            charge = 0.0
        head_index, tail_index = template.connect
        head_atom = template.atoms[head_index - 1].name if head_index else None
        tail_atom = template.atoms[tail_index - 1].name if tail_index else None
        entry = {
                "template": template.name,
                "source_unit_name": template.unit_name,
                "atom_count": len(template.atoms),
                "bond_count": len(template.bonds),
                "total_charge": round(charge, 8),
                "expected_integer_charge": round(charge),
                "head_atom": head_atom,
                "tail_atom": tail_atom,
                "source_connect_indices": [head_index, tail_index],
                "source_connectivity_flags": sorted({bond[2] for bond in template.bonds}),
        }
        entry.update(_connection_metadata(template, "head", head_atom))
        entry.update(_connection_metadata(template, "tail", tail_atom))
        entries.append(entry)
    return {
        "format_version": 2,
        "source": source.name,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "source_license": source_license,
        "template_count": len(entries),
        "templates": entries,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Amber OFF/lib input")
    parser.add_argument("mol2", type=Path, help="generated MOL2 output")
    parser.add_argument("manifest", type=Path, help="generated JSON manifest")
    parser.add_argument(
        "--source-license", default="Apache-2.0",
        help="license/provenance label recorded in the generated manifest",
    )
    args = parser.parse_args()

    templates = parse_off(args.source)
    args.mol2.write_text(render_mol2(args.source, templates), encoding="utf-8")
    args.manifest.write_text(
        json.dumps(
            build_manifest(args.source, templates, args.source_license),
            indent=2,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
