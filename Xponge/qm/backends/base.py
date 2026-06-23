"""Base protocol for QM backend adapters."""

from __future__ import annotations

from typing import Protocol

from ..capabilities import QMCapabilitySet
from ..models import ESPGridRequest, ESPResult, HessianResult, OptimizationResult, QMMolecule, QMRunOptions, SCFResult


class QMBackend(Protocol):
    name: str

    def capabilities(self) -> QMCapabilitySet: ...

    def run_scf(self, molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> SCFResult: ...

    def compute_esp(self, scf_result: SCFResult, request: ESPGridRequest) -> ESPResult: ...

    def optimize_geometry(
        self, molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False
    ) -> OptimizationResult: ...

    def compute_hessian(self, molecule: QMMolecule, options: QMRunOptions, assign=None, return_timings: bool = False) -> HessianResult: ...
