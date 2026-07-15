"""Safe materialization of legacy SPONGE input payloads."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
import os
import tempfile

from .errors import BundleConflictError, BundlePathError


@dataclass(frozen=True)
class LegacyPayload:
    """One legacy input file generated from a bundled contract."""

    key: str
    data: str | bytes
    filename: str | None = None
    binary: bool = False
    source_kind: str = "typed"
    source_path: str | None = None
    bind_in_mdin: bool = True


class LegacyMaterializer:
    """Plan and atomically write reverse-converted legacy files."""

    def __init__(self, output_dir: str | Path, *, overwrite: bool = False):
        self.output_dir = Path(output_dir).resolve()
        self.overwrite = overwrite
        self._planned: dict[Path, LegacyPayload] = {}

    @property
    def planned(self) -> dict[Path, LegacyPayload]:
        return dict(self._planned)

    def plan_payload(self, payload: LegacyPayload, filename: str) -> Path:
        target = self._safe_target(filename)
        payload = replace(payload, filename=filename)
        previous = self._planned.get(target)
        if previous is not None and previous != payload:
            raise BundleConflictError(f"multiple payloads target the same path: {target}")
        self._planned[target] = payload
        return target

    def validate_targets(self) -> None:
        if self.overwrite:
            return
        conflicts = [path for path in self._planned if path.exists()]
        if conflicts:
            rendered = ", ".join(str(path) for path in conflicts)
            raise BundleConflictError(f"legacy output targets already exist: {rendered}")

    def write_all(self) -> None:
        self.validate_targets()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        for target, payload in self._planned.items():
            target.parent.mkdir(parents=True, exist_ok=True)
            data = payload.data if isinstance(payload.data, bytes) else payload.data.encode("utf-8")
            _atomic_write(target, data)

    def _safe_target(self, filename: str) -> Path:
        raw = Path(filename)
        if raw.is_absolute():
            raise BundlePathError(f"legacy target filename must be relative: {filename}")
        target = (self.output_dir / raw).resolve()
        try:
            target.relative_to(self.output_dir)
        except ValueError as exc:
            raise BundlePathError(f"legacy target escapes output root: {filename}") from exc
        return target


def _atomic_write(target: Path, data: bytes) -> None:
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{target.name}.", dir=target.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
        os.replace(temporary_name, target)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise
