#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a CHARMM-GUI GROMACS directory into Xponge/SPONGE inputs and a PDB."
    )
    parser.add_argument(
        "gromacs_dir",
        nargs="?",
        default="/home/ylj/文档/CHARMM_DEMO/charmm-gui-7398083341/gromacs",
        help="Directory containing topol.top and step3_input.gro",
    )
    parser.add_argument(
        "--prefix",
        default="step3_input",
        help="Output file prefix, without extension",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Directory for generated files. Defaults to the input GROMACS directory.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    gromacs_dir = Path(args.gromacs_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else gromacs_dir
    topol = gromacs_dir / "topol.top"
    gro = gromacs_dir / "step3_input.gro"

    if not topol.exists():
        raise FileNotFoundError(f"Cannot find topology file: {topol}")
    if not gro.exists():
        raise FileNotFoundError(f"Cannot find coordinate file: {gro}")

    output_dir.mkdir(parents=True, exist_ok=True)

    import Xponge
    import Xponge.forcefield.charmm.charmm36
    import Xponge.forcefield.charmm.tip3p_charmm

    system, molecules = Xponge.load_molitp(str(topol), water_replace=False)
    Xponge.load_gro(str(gro), system)

    prefix = output_dir / args.prefix
    Xponge.save_sponge_input(system, prefix.name, str(output_dir))
    Xponge.save_pdb(system, str(prefix.with_suffix(".pdb")))

    print(f"[ok] loaded system: {system.name}")
    print(f"[ok] molecule types: {', '.join(sorted(molecules.keys()))}")
    print(f"[ok] residues: {len(system.residues)}")
    print(f"[ok] atoms: {len(system.get_atoms())}")
    print(f"[ok] box: {system.box_length} angles={getattr(system, 'box_angle', None)}")
    print(f"[ok] SPONGE prefix: {prefix}")
    print(f"[ok] PDB: {prefix.with_suffix('.pdb')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
