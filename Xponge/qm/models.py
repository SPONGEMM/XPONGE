"""Backend-neutral QM data models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


SUPPORTED_SCF_REFERENCES = frozenset({"auto", "rhf", "rohf", "uhf"})


def resolve_scf_reference(reference: str | None, spin: int) -> str:
    """Resolve a requested SCF reference for a molecule with twice-spin ``spin``."""

    reference_name = "auto" if reference is None else str(reference).strip().lower()
    if reference_name not in SUPPORTED_SCF_REFERENCES:
        supported = ", ".join(sorted(SUPPORTED_SCF_REFERENCES))
        raise ValueError(f"SCF reference should be one of: {supported}")
    if reference_name == "auto":
        return "rhf" if int(spin) == 0 else "rohf"
    if reference_name == "rhf" and int(spin) != 0:
        raise ValueError("RHF requires a closed-shell molecule with spin=0")
    if reference_name == "rohf" and int(spin) == 0:
        raise ValueError("ROHF requires an open-shell molecule with spin>0")
    return reference_name


@dataclass(slots=True)
class QMMolecule:
    atom_symbols: list[str]
    coordinates_angstrom: list[tuple[float, float, float]]
    total_charge: int
    spin: int = 0
    atom_names: list[str] | None = None
    formal_charges: list[int] | None = None
    bonds: list[dict[int, int]] | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class QMRunOptions:
    backend: str | None = None
    basis: Any = "6-31g*"
    ecp: Any = None
    cart: bool | None = None
    method: str = "scf"
    reference: str = "auto"
    optimize_geometry: bool = False
    threads: int | None = None
    memory: str | None = None
    memory_limit_bytes: int | None = None
    scf_convergence_tolerance: float | None = None
    scf_max_cycles: int | None = None
    scf_strategy: str = "direct"
    properties: tuple[str, ...] = ()


@dataclass(slots=True)
class ESPGridRequest:
    grid_points_bohr: Any
    include_nuclear_term: bool = False
    include_electronic_term: bool = True
    memory_limit_bytes: int | None = None
    chunk_policy: str = "auto"
    safety_factor: float = 0.8


@dataclass(slots=True)
class SCFResult:
    backend_name: str
    total_energy: float | None
    converged: bool
    coordinates_bohr: Any
    nuclear_charges: Any
    charge: int
    spin: int
    atom_symbols: list[str]
    reference: str = ""
    backend_handle: Any = None
    timings: dict[str, float] = field(default_factory=dict)
    optimized_coordinates_angstrom: list[tuple[float, float, float]] | None = None


@dataclass(slots=True)
class ESPResult:
    grid_points_bohr: Any
    electronic_esp_au: Any
    total_esp_au: Any | None = None
    nuclear_esp_au: Any | None = None
    timings: dict[str, float] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class OptimizationResult:
    optimized_coordinates_angstrom: list[tuple[float, float, float]]
    converged: bool
    iterations: int | None = None
    final_energy: float | None = None
    timings: dict[str, float] = field(default_factory=dict)


@dataclass(slots=True)
class HessianResult:
    cartesian_hessian_au: Any
    coordinates_angstrom: list[tuple[float, float, float]]
    atom_symbols: list[str]
    scf_converged: bool = False
    total_energy: float | None = None
    reference: str = ""
    timings: dict[str, float] = field(default_factory=dict)
