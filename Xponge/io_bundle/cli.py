"""
CLI entry point for SPONGE legacy-to-bundle conversion.
"""

from __future__ import annotations

from .converter import convert_legacy_to_bundle
from .output_writer import convert_legacy_outputs_to_bundle


def add_legacy_to_bundle_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "legacy-to-bundle",
        help="convert a legacy SPONGE input directory to bundled H5 input files",
    )
    parser.add_argument("case_root", help="legacy case directory")
    parser.add_argument("-m", "--mdin", default="mdin.spg.toml", help="legacy mdin path relative to case root")
    parser.add_argument("-o", "--output", required=True, help="output directory for bundle and manifest")
    parser.add_argument("--dry-run", action="store_true", help="scan and build manifest without writing H5 files")
    parser.set_defaults(func=main)

    output_parser = subparsers.add_parser(
        "legacy-outputs-to-bundle",
        help="convert existing legacy SPONGE output files to an H5MD output bundle",
    )
    output_parser.add_argument("case_root", help="legacy case directory")
    output_parser.add_argument("-m", "--mdin", default="mdin.spg.toml", help="legacy mdin path relative to case root")
    output_parser.add_argument("-o", "--output", required=True, help="output H5MD bundle path")
    output_parser.add_argument("--atom-count", type=int, default=None, help="atom count for vector trajectory outputs")
    output_parser.add_argument("--dry-run", action="store_true", help="scan and build manifest without writing H5 files")
    output_parser.set_defaults(func=output_main)


def main(args) -> None:
    manifest = convert_legacy_to_bundle(args.case_root, args.output, mdin=args.mdin, dry_run=args.dry_run)
    converted = sum(1 for entry in manifest.entries if entry.status == "converted")
    typed_converted = sum(1 for entry in manifest.entries if entry.status == "typed_converted")
    compatibility_imported = sum(1 for entry in manifest.entries if entry.status == "compatibility_imported")
    print(f"converted entries: {converted}")
    print(f"typed converted entries: {typed_converted}")
    if compatibility_imported:
        print(f"compatibility imported entries: {compatibility_imported}")
    unsupported = [entry for entry in manifest.entries if entry.status == "unsupported"]
    if unsupported:
        print(f"unsupported entries: {len(unsupported)}")


def output_main(args) -> None:
    manifest = convert_legacy_outputs_to_bundle(
        args.case_root,
        args.output,
        mdin=args.mdin,
        atom_count=args.atom_count,
        dry_run=args.dry_run,
    )
    typed_converted = sum(1 for entry in manifest.entries if entry.status == "typed_converted")
    missing = sum(1 for entry in manifest.entries if entry.status == "legacy_output_missing")
    print(f"typed converted output entries: {typed_converted}")
    if missing:
        print(f"missing legacy output entries: {missing}")
