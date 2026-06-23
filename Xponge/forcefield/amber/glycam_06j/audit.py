"""
Utilities for auditing GLYCAM template coverage against AmberTools.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

FUNCTIONAL_GROUP_TEMPLATES = {"CA2", "MEX", "OME", "ROH", "SO3", "TBT"}
EXTERNAL_COVERAGE = {"CA2": "covered by Amber ion templates outside glycam_06j"}
NAME_EQUIVALENTS = {}

_PREP_UNIT_RE = re.compile(r"^([A-Za-z0-9]{3,4})\s+INT\s+0\s*$")
_LIB_INDEX_RE = re.compile(r'^\s*"([A-Za-z0-9]{3,4})"\s*$')


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def _default_ambertools_root() -> Path:
    return Path("/mnt/data8t/Software/AmberTools26/ambertools26_src")


def _parse_prep_units(prep_file: Path) -> set[str]:
    names = set()
    with prep_file.open() as f:
        for line in f:
            match = _PREP_UNIT_RE.match(line.strip())
            if match:
                names.add(match.group(1))
    return names


def _parse_lib_units(lib_file: Path) -> set[str]:
    names = set()
    in_index = False
    with lib_file.open() as f:
        for line in f:
            stripped = line.strip()
            if stripped == "!!index array str":
                in_index = True
                continue
            if in_index and stripped.startswith("!entry."):
                break
            if in_index:
                match = _LIB_INDEX_RE.match(stripped)
                if match:
                    names.add(match.group(1))
    return names


def _parse_mol2_units(mol2_file: Path) -> set[str]:
    names = set()
    in_substructure = False
    with mol2_file.open() as f:
        for line in f:
            stripped = line.strip()
            if stripped == "@<TRIPOS>SUBSTRUCTURE":
                in_substructure = True
                continue
            if stripped.startswith("@<TRIPOS>") and stripped != "@<TRIPOS>SUBSTRUCTURE":
                in_substructure = False
            elif in_substructure and stripped:
                parts = stripped.split()
                if len(parts) >= 2:
                    names.add(parts[1])
    return names


def gather_xponge_glycam_units(repo_root: Path | None = None) -> set[str]:
    repo_root = repo_root or _default_repo_root()
    glycam_dir = repo_root / "Xponge" / "forcefield" / "amber" / "glycam_06j"
    units = set()
    for mol2_file in sorted(glycam_dir.glob("*.mol2")):
        units |= _parse_mol2_units(mol2_file)
    return units


def audit_glycam_coverage(
    prep_file: Path | None = None,
    lib_files: list[Path] | None = None,
    repo_root: Path | None = None,
) -> dict:
    repo_root = repo_root or _default_repo_root()
    amber_root = _default_ambertools_root()
    prep_file = prep_file or amber_root / "dat" / "leap" / "prep" / "GLYCAM_06j-1.prep"
    lib_files = lib_files or [
        amber_root / "dat" / "leap" / "lib" / "GLYCAM_amino_06j_12SB.lib",
        amber_root / "dat" / "leap" / "lib" / "GLYCAM_aminont_06j_12SB.lib",
        amber_root / "dat" / "leap" / "lib" / "GLYCAM_aminoct_06j_12SB.lib",
    ]

    amber_units = _parse_prep_units(prep_file)
    for lib_file in lib_files:
        amber_units |= _parse_lib_units(lib_file)

    xponge_units = gather_xponge_glycam_units(repo_root)

    covered = sorted(name for name in amber_units if name in xponge_units)
    equivalent = {
        amber_name: real_name
        for amber_name, real_name in NAME_EQUIVALENTS.items()
        if amber_name in amber_units and real_name in xponge_units
    }
    covered_elsewhere = {
        amber_name: note
        for amber_name, note in EXTERNAL_COVERAGE.items()
        if amber_name in amber_units and amber_name not in xponge_units
    }

    accounted = set(covered) | set(equivalent) | set(covered_elsewhere)
    missing = amber_units - accounted
    missing_functional_groups = sorted(name for name in missing if name in FUNCTIONAL_GROUP_TEMPLATES)
    missing_modified_monosaccharides = sorted(name for name in missing if name not in FUNCTIONAL_GROUP_TEMPLATES)

    return {
        "amber_unit_count": len(amber_units),
        "xponge_unit_count": len(xponge_units),
        "covered": covered,
        "equivalent": equivalent,
        "covered_elsewhere": covered_elsewhere,
        "missing_functional_groups": missing_functional_groups,
        "missing_modified_monosaccharides": missing_modified_monosaccharides,
    }


def _main():
    parser = argparse.ArgumentParser(description="Audit Amber GLYCAM coverage in Xponge.")
    parser.add_argument("--prep-file", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--lib-file", type=Path, action="append", default=None)
    args = parser.parse_args()
    report = audit_glycam_coverage(args.prep_file, args.lib_file, args.repo_root)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
