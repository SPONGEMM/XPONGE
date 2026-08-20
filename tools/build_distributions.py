"""Build the canonical Xponge distribution and the mokda-xponge alias."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
VERSION_FILE = ROOT / "Xponge" / "__init__.py"
ALIAS_TEMPLATE_DIR = ROOT / "packaging" / "mokda-xponge"


def read_version() -> str:
    """Read ``Xponge.__version__`` without importing the package."""
    module = ast.parse(VERSION_FILE.read_text(encoding="utf-8"), VERSION_FILE.as_posix())
    for statement in module.body:
        if not isinstance(statement, (ast.Assign, ast.AnnAssign)):
            continue
        targets = statement.targets if isinstance(statement, ast.Assign) else [statement.target]
        if not any(isinstance(target, ast.Name) and target.id == "__version__" for target in targets):
            continue
        value = ast.literal_eval(statement.value)
        if not isinstance(value, str) or not value:
            break
        return value
    raise RuntimeError(f"Could not find a string __version__ in {VERSION_FILE}")


def build(source: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, "-m", "build", "--outdir", str(output), str(source)],
        check=True,
    )


def build_alias(version: str, output: Path) -> None:
    template = (ALIAS_TEMPLATE_DIR / "pyproject.toml.template").read_text(encoding="utf-8")
    if template.count("@VERSION@") != 2:
        raise RuntimeError("Alias template must contain exactly two @VERSION@ placeholders")

    with tempfile.TemporaryDirectory(prefix="xponge-alias-") as temporary_directory:
        source = Path(temporary_directory)
        (source / "pyproject.toml").write_text(
            template.replace("@VERSION@", version),
            encoding="utf-8",
        )
        shutil.copy2(ALIAS_TEMPLATE_DIR / "README.md", source / "README.md")
        build(source, output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "dist",
        help="Parent directory for the xponge and mokda-xponge artifact directories",
    )
    parser.add_argument(
        "--expected-version",
        help="Fail unless this version (optionally prefixed with v) matches Xponge.__version__",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    version = read_version()
    if args.expected_version and args.expected_version.removeprefix("v") != version:
        raise RuntimeError(
            f"Release version mismatch: tag {args.expected_version!r}, Xponge.__version__ {version!r}"
        )

    output_root = args.outdir.resolve()
    build(ROOT, output_root / "xponge")
    build_alias(version, output_root / "mokda-xponge")


if __name__ == "__main__":
    main()
