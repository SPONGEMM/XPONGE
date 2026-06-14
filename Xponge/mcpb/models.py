"""Structured request/result models for MCPB workflows."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(slots=True)
class MCPBIonInfo:
    """Explicit user-facing metal-ion description for MCPB."""

    atom_id: int
    element: str
    formal_charge: int | None = None
    spin: int | None = None
    resname: str | None = None
    atom_name: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class MCPBSelection:
    """Selected metal-center environment used by the workflow."""

    ion_atom_ids: tuple[int, ...]
    coordinating_atom_ids: tuple[int, ...]
    selected_residue_ids: tuple[int, ...]
    bonded_pairs: tuple[tuple[int, int], ...]


@dataclass(slots=True)
class MCPBLocalModel:
    """Locally extracted MCPB submodel plus source-index mappings."""

    molecule: Any
    source_atom_ids: tuple[int, ...]
    atom_id_map: dict[int, int] = field(default_factory=dict)
    residue_id_map: dict[int, int] = field(default_factory=dict)


@dataclass(slots=True)
class MCPBRequest:
    """Normalized MCPB request."""

    molecule: Any
    ion_ids: tuple[int, ...]
    ion_info: tuple[MCPBIonInfo, ...]
    method: str = "seminario"
    model: str = "bonded"
    cutoff: float = 2.8
    bonded_pairs: tuple[tuple[int, int], ...] = ()
    additional_residue_ids: tuple[int, ...] = ()
    charge_mode: str = "resp"
    qm_backend: str = "pyscf"
    basis: str | None = None
    scale_factor: float = 1.0
    frcmod_files: tuple[str, ...] = ()
    gaff: int | None = None
    force_field: str | None = None
    water_model: str | None = None


@dataclass(slots=True)
class MCPBResult:
    """Structured MCPB result."""

    molecule: Any
    request: MCPBRequest
    selection: MCPBSelection
    small_model: MCPBLocalModel | None = None
    large_model: MCPBLocalModel | None = None
    frcmod_path: str | None = None
    connect_records: list[tuple[int, int]] = field(default_factory=list)
    updated_charge_atoms: list[int] = field(default_factory=list)
    registered_metal_templates: list[str] = field(default_factory=list)
    sponge_ready: bool = False
    pending_requirements: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
