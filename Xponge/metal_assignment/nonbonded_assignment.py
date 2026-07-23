"""Composition of nonbonded base and ion overlays into an apply-ready result."""

from __future__ import annotations

from typing import Any, Mapping

from .base_composition import BaseCompositionOutput
from .contracts import (
    MetalAssignmentInput,
    ParameterizationResult,
    ValidationError,
    _canonicalize,
    validate_input,
    validate_result,
)
from .metal_overlay import MetalOverlayBuildOutput
from .partial_charge import compose_partial_charge_artifacts


def compose_nonbonded_assignment(
    request: MetalAssignmentInput,
    base: BaseCompositionOutput,
    metal: MetalOverlayBuildOutput,
    *,
    package_hash: str,
    provenance: Mapping[str, Any] | None = None,
) -> ParameterizationResult:
    """Create a validated result without making topology or charge decisions.

    Base and ion providers remain independently auditable.  This function only
    closes their hashes and identities into the common result consumed by the
    transactional ``apply`` operation.
    """

    validate_input(request)
    if request.interaction_model != "nonbonded_12_6":
        raise ValidationError(
            "nonbonded_composition_interaction_model_mismatch",
            request.interaction_model,
        )
    if (
        not isinstance(package_hash, str)
        or len(package_hash) != 64
        or any(character not in "0123456789abcdef" for character in package_hash)
    ):
        raise ValidationError("missing_package_hash", "nonbonded composition requires package identity")
    topology_hash = request.topology.topology_hash
    if (
        base.overlay.topology_hash != topology_hash
        or base.report.topology_hash != topology_hash
        or metal.overlay.topology_hash != topology_hash
        or metal.report.topology_hash != topology_hash
    ):
        raise ValidationError("overlay_topology_mismatch", "nonbonded provider output targets another topology")
    if (
        base.report.projection_hash != request.projection_hash
        or metal.report.projection_hash != request.projection_hash
    ):
        raise ValidationError("overlay_projection_mismatch", "nonbonded provider output targets another projection")
    if not base.report.report_hash or base.report.report_hash != base.report.computed_hash():
        raise ValidationError("stale_base_assignment_report_hash", request.request_id)
    if not metal.report.report_hash or metal.report.report_hash != metal.report.computed_hash():
        raise ValidationError("stale_metal_assignment_report_hash", request.request_id)
    if metal.report.interaction_model != request.interaction_model:
        raise ValidationError("metal_assignment_interaction_model_mismatch", metal.report.interaction_model)

    charge_overlay, charge_report = compose_partial_charge_artifacts(request)

    result = ParameterizationResult(
        schema_version=request.schema_version,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology_hash=topology_hash,
        projection_hash=request.projection_hash,
        base_overlay=base.overlay,
        metal_overlay=metal.overlay,
        charge_overlay=charge_overlay,
        fit_reports={
            "base_assignment": _canonicalize(base.report),
            "metal_assignment": _canonicalize(metal.report),
            "partial_charge_composition": _canonicalize(charge_report),
        },
        provenance={
            "package_hash": package_hash,
            "interaction_model": request.interaction_model,
            "composition": "xponge:nonbonded-overlay-composition:v1",
            **dict(provenance or {}),
        },
        status="overlay_validated",
        complete=False,
    ).with_computed_hash()
    validate_result(request, result)
    return result


__all__ = ["compose_nonbonded_assignment"]
