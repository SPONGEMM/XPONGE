"""
CLI entry point for SPONGE legacy-to-bundle conversion.
"""

from __future__ import annotations

from .converter import convert_legacy_to_bundle
from .output_writer import convert_legacy_outputs_to_bundle
from .reverse_converter import convert_bundle_to_legacy


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
    output_parser.add_argument(
        "--input-topology",
        default=None,
        help="sponge.input.v2 topology bundle used to make imported restart output launchable",
    )
    output_parser.add_argument(
        "--input-protocol",
        default=None,
        help="optional matching sponge.input.v2 protocol bundle for restart lineage",
    )
    output_parser.add_argument("--dry-run", action="store_true", help="scan and build manifest without writing H5 files")
    output_parser.set_defaults(func=output_main)

    reverse_parser = subparsers.add_parser(
        "bundle-to-legacy",
        help="convert bundled H5 input files to direct/legacy SPONGE inputs",
    )
    reverse_parser.add_argument("bundle_root", help="directory containing the bundled mdin")
    reverse_parser.add_argument(
        "-m",
        "--mdin",
        default="mdin.bundled.spg.toml",
        help="bundled mdin path relative to bundle root",
    )
    reverse_parser.add_argument("-o", "--output", required=True, help="legacy output directory")
    reverse_parser.add_argument("--prefix", default=None, help="legacy input filename prefix")
    reverse_parser.add_argument(
        "--no-strict",
        action="store_false",
        dest="strict",
        help="record unsupported optional contracts instead of failing",
    )
    reverse_parser.add_argument(
        "--overwrite",
        action="store_true",
        help="allow replacement of planned legacy output files",
    )
    reverse_parser.add_argument("--dry-run", action="store_true", help="validate and plan without writing files")
    reverse_parser.set_defaults(func=reverse_main, strict=True)


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
        input_topology_path=args.input_topology,
        input_protocol_path=args.input_protocol,
        dry_run=args.dry_run,
    )
    typed_converted = sum(1 for entry in manifest.entries if entry.status == "typed_converted")
    missing = sum(1 for entry in manifest.entries if entry.status == "legacy_output_missing")
    print(f"typed converted output entries: {typed_converted}")
    if missing:
        print(f"missing legacy output entries: {missing}")


def reverse_main(args) -> None:
    manifest = convert_bundle_to_legacy(
        args.bundle_root,
        args.output,
        mdin=args.mdin,
        prefix=args.prefix,
        strict=args.strict,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
    )
    typed = sum(1 for entry in manifest.entries if entry.status == "typed_exported")
    fallback = sum(
        1
        for entry in manifest.entries
        if entry.status in {"sidecar_restored", "embedded_exported", "compatibility_restored"}
    )
    unsupported = sum(1 for entry in manifest.entries if entry.status == "unsupported")
    print(f"typed exported entries: {typed}")
    if fallback:
        print(f"fallback restored entries: {fallback}")
    if unsupported:
        print(f"unsupported entries: {unsupported}")
    if manifest.generated_mdin is not None:
        print(f"generated mdin: {manifest.generated_mdin}")
