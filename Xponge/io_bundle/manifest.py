"""
Manifest model for SPONGE legacy-to-bundle conversion.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
from pathlib import Path
from typing import Any


@dataclass
class ManifestEntry:
    """One converted or unsupported I/O contract."""

    contract_id: str
    status: str
    source_key: str | None = None
    source_path: str | None = None
    bundle_file: str | None = None
    bundle_path: str | None = None
    direction: str | None = None
    component: str | None = None
    payload_kind: str | None = None
    override_policy: str | None = None
    comparison_rule: str | None = None
    message: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {key: value for key, value in self.__dict__.items() if value is not None}


@dataclass
class ConversionManifest:
    """Machine-readable conversion manifest."""

    schema: str = "xponge.legacy_to_bundle.manifest"
    schema_version: int = 1
    case_root: str = ""
    mode: str = "normal"
    entries: list[ManifestEntry] = field(default_factory=list)
    bundled_mdin: str | None = None

    def add(self, entry: ManifestEntry) -> None:
        self.entries.append(entry)

    def to_dict(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "schema": self.schema,
            "schema_version": self.schema_version,
            "case_root": self.case_root,
            "mode": self.mode,
            "entries": [entry.to_dict() for entry in self.entries],
        }
        if self.bundled_mdin is not None:
            data["bundled_mdin"] = self.bundled_mdin
        return data

    def write(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
