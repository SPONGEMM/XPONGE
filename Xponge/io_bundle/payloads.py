"""Path-like payload adapters for SPONGE input codecs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


class InputPayload(Protocol):
    """Minimal path-like interface consumed by legacy input parsers."""

    @property
    def suffix(self) -> str: ...

    def read_text(self, encoding: str = "utf-8") -> str: ...

    def read_bytes(self) -> bytes: ...


@dataclass(frozen=True)
class PathPayload:
    """Adapter around a filesystem path."""

    path: Path

    @property
    def suffix(self) -> str:
        return self.path.suffix

    def read_text(self, encoding: str = "utf-8") -> str:
        return self.path.read_text(encoding=encoding)

    def read_bytes(self) -> bytes:
        return self.path.read_bytes()

    def __str__(self) -> str:
        return str(self.path)


@dataclass(frozen=True)
class MemoryPayload:
    """In-memory text or bytes exposed through the parser path protocol."""

    name: str
    data: str | bytes

    @property
    def suffix(self) -> str:
        return Path(self.name).suffix

    def read_text(self, encoding: str = "utf-8") -> str:
        if isinstance(self.data, str):
            return self.data
        return self.data.decode(encoding)

    def read_bytes(self) -> bytes:
        if isinstance(self.data, bytes):
            return self.data
        return self.data.encode("utf-8")

    def __str__(self) -> str:
        return self.name
