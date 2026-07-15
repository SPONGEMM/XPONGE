"""
Legacy SPONGE case scanning and mdin parsing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import shlex


_KEY_VALUE_RE = re.compile(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?)\s*(?:#.*)?$")
_SECTION_RE = re.compile(r"^\s*\[([A-Za-z_][A-Za-z0-9_.]*)\]\s*(?:#.*)?$")


@dataclass
class LegacyCase:
    """A scanned legacy SPONGE case directory."""

    root: Path
    mdin_path: Path
    mdin_text: str
    commands: dict[str, str]

    @property
    def mode(self) -> str:
        return self.commands.get("mode", "normal").strip().lower()

    def resolve_value_path(self, value: str) -> Path:
        path = Path(value)
        if path.is_absolute():
            return path
        return self.root / path

    def resolve_legacy_input_path(self, key: str) -> Path | None:
        """Resolve explicit or ``default_in_file_prefix`` legacy input paths."""

        value = self.commands.get(key)
        if value:
            return self.resolve_value_path(value)
        suffix = "_in_file"
        if not key.endswith(suffix):
            return None
        prefix = self.commands.get("default_in_file_prefix")
        if not prefix:
            return None
        stem = key[: -len(suffix)]
        if stem.endswith("_"):
            stem = stem[:-1]
        candidate = self.resolve_value_path(f"{prefix}_{stem}.txt")
        if candidate.exists():
            return candidate
        return None


def parse_mdin_text(text: str) -> dict[str, str]:
    """Parse the flat key/value subset of SPONGE mdin TOML used by legacy I/O."""

    commands: dict[str, str] = {}
    section: str | None = None
    for line in text.splitlines():
        section_match = _SECTION_RE.match(line)
        if section_match:
            section = section_match.group(1).replace(".", "_")
            continue
        match = _KEY_VALUE_RE.match(line)
        if not match:
            continue
        key, raw_value = match.groups()
        if section:
            key = f"{section}_{key}"
        raw_value = raw_value.strip()
        if not raw_value:
            commands[key] = ""
            continue
        try:
            parts = shlex.split(raw_value, comments=False, posix=True)
        except ValueError:
            parts = []
        commands[key] = parts[0] if len(parts) == 1 else raw_value.strip('"').strip("'")
    return commands


def render_mdin_without_keys(text: str, omit_keys: set[str], append_lines: list[str]) -> str:
    """Render mdin text after removing normalized command keys.

    Sectioned TOML commands such as ``[EAM] in_file`` are normalized with the
    same ``EAM_in_file`` rule used by ``parse_mdin_text``.
    """

    output: list[str] = []
    pending_section: str | None = None
    pending_section_lines: list[str] = []
    section_has_payload = False
    section: str | None = None

    def flush_section() -> None:
        nonlocal pending_section, pending_section_lines, section_has_payload
        if pending_section is not None and section_has_payload:
            output.extend(pending_section_lines)
        pending_section = None
        pending_section_lines = []
        section_has_payload = False

    for line in text.splitlines():
        section_match = _SECTION_RE.match(line)
        if section_match:
            flush_section()
            section = section_match.group(1).replace(".", "_")
            pending_section = section
            pending_section_lines = [line]
            section_has_payload = False
            continue

        match = _KEY_VALUE_RE.match(line)
        normalized_key = None
        if match:
            key = match.group(1)
            normalized_key = f"{section}_{key}" if section else key
        if normalized_key in omit_keys:
            continue

        target = pending_section_lines if pending_section is not None else output
        target.append(line)
        if pending_section is not None:
            if normalized_key is not None or line.strip():
                section_has_payload = True

    flush_section()
    rendered = output + append_lines
    return "\n".join(rendered).rstrip() + "\n"


def scan_legacy_case(case_root: str | Path, mdin: str | Path = "mdin.spg.toml") -> LegacyCase:
    """Read a legacy case directory and parse its mdin file."""

    root = Path(case_root).resolve()
    mdin_path = Path(mdin)
    if not mdin_path.is_absolute():
        mdin_path = root / mdin_path
    if not mdin_path.exists():
        raise FileNotFoundError(f"legacy mdin file does not exist: {mdin_path}")
    mdin_text = mdin_path.read_text(encoding="utf-8")
    commands = parse_mdin_text(mdin_text)
    return LegacyCase(root=root, mdin_path=mdin_path, mdin_text=mdin_text, commands=commands)
