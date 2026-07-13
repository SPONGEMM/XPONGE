"""Registry for bundled input to legacy payload exporters."""

from __future__ import annotations

from dataclasses import dataclass

from .errors import BundleCapabilityError
from .state_exporters import STATE_EXPORTERS
from .topology_exporters import TOPOLOGY_EXPORTERS
from .trajectory_exporters import TRAJECTORY_EXPORTERS


@dataclass(frozen=True)
class ExportContext:
    """Shared options available to typed exporters."""

    mode: str
    prefix: str
    particle_stream: str = "all"


EXPORTERS = {
    **TOPOLOGY_EXPORTERS,
    **STATE_EXPORTERS,
    **TRAJECTORY_EXPORTERS,
}


def get_exporter(exporter_id: str):
    """Return a registered exporter or raise an explicit capability error."""

    try:
        return EXPORTERS[exporter_id]
    except KeyError as exc:
        raise BundleCapabilityError(f"no bundle exporter is registered for {exporter_id}") from exc


def validate_exporter_registry() -> None:
    """Reject malformed exporter registry entries."""

    for exporter_id, exporter in EXPORTERS.items():
        if not exporter_id or not callable(exporter):
            raise ValueError(f"invalid bundle exporter registry entry: {exporter_id!r}")
