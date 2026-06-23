"""Capability flags for QM backends."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class QMCapabilitySet:
    supports_scf: bool = True
    supports_esp: bool = True
    supports_geometry_optimization: bool = False
    supports_hessian: bool = False
    supports_open_shell: bool = True
    supports_point_charges: bool = False
    supports_constraints: bool = False
