"""Deterministic composition of disjoint base-force-field provider fragments."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any, Mapping, Sequence

from .base_assignment import BaseAssignmentOutput
from .contracts import BaseForceFieldOverlay, ValidationError, _canonicalize, _freeze_field, _sha256
from .input import MetalAssignmentPackage, validate_package


@dataclass(frozen=True, slots=True)
class BaseCompositionReport:
    schema_version: int
    topology_hash: str
    projection_hash: str
    providers: tuple[str, ...]
    provider_report_hashes: Mapping[str, str]
    covered_atom_ids: tuple[str, ...]
    bonded_term_count: int
    parameter_source: str
    report_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "provider_report_hashes")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "topology_hash": self.topology_hash,
            "projection_hash": self.projection_hash,
            "providers": _canonicalize(self.providers),
            "provider_report_hashes": _canonicalize(self.provider_report_hashes),
            "covered_atom_ids": _canonicalize(self.covered_atom_ids),
            "bonded_term_count": self.bonded_term_count,
            "parameter_source": self.parameter_source,
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "BaseCompositionReport":
        return replace(self, report_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class BaseCompositionOutput:
    overlay: BaseForceFieldOverlay
    report: BaseCompositionReport


def _merge_map(
    target: dict[str, Any],
    source: Mapping[str, Any],
    *,
    code: str,
) -> None:
    overlap = set(target) & set(source)
    if overlap:
        raise ValidationError(code, ",".join(sorted(overlap)))
    target.update(source)


def compose_base_force_field(
    package: MetalAssignmentPackage,
    fragments: Sequence[BaseAssignmentOutput],
) -> BaseCompositionOutput:
    """Merge provider fragments only when their scopes are complete and disjoint."""

    validate_package(package)
    expected_ids = {
        atom_id
        for component in package.request.assignment_components
        if component.base_force_field
        for atom_id in component.atom_ids
    }
    if not fragments and expected_ids:
        raise ValidationError("missing_base_assignment_fragment", "at least one provider fragment is required")
    ordered = sorted(
        fragments,
        key=lambda item: (item.report.provider, item.report.parameter_source, item.report.report_hash),
    )
    topology_hash = package.request.topology.topology_hash
    projection_hash = package.request.projection_hash
    atom_order = {atom.external_id: atom.stable_order for atom in package.request.topology.atoms}
    covered: set[str] = set()
    atom_types: dict[str, str] = {}
    charges: dict[str, float] = {}
    masses: dict[str, float] = {}
    lj_parameters: dict[str, Any] = {}
    bonded_parameters: dict[str, Any] = {}
    provider_report_hashes: dict[str, str] = {}
    for fragment in ordered:
        overlay = fragment.overlay
        report = fragment.report
        if overlay.topology_hash != topology_hash or report.topology_hash != topology_hash:
            raise ValidationError("base_fragment_topology_mismatch", report.provider)
        if report.projection_hash != projection_hash:
            raise ValidationError("base_fragment_projection_mismatch", report.provider)
        if not report.report_hash or report.report_hash != report.computed_hash():
            raise ValidationError("stale_base_assignment_report_hash", report.provider)
        fragment_ids = set(overlay.covered_atom_ids)
        if len(fragment_ids) != len(overlay.covered_atom_ids):
            raise ValidationError("duplicate_overlay_atom", report.provider)
        overlap = covered & fragment_ids
        if overlap:
            raise ValidationError("overlapping_base_provider_scope", ",".join(sorted(overlap)))
        for name, values in (
            ("atom_types", overlay.atom_types),
            ("masses", overlay.masses),
            ("lj_parameters", overlay.lj_parameters),
        ):
            if set(values) != fragment_ids:
                raise ValidationError("incomplete_base_fragment_parameters", f"{report.provider}:{name}")
        if not set(overlay.charges) <= fragment_ids:
            raise ValidationError("base_fragment_charge_scope_mismatch", report.provider)
        if report.provider in provider_report_hashes:
            raise ValidationError("duplicate_base_provider_fragment", report.provider)
        provider_report_hashes[report.provider] = report.report_hash
        covered.update(fragment_ids)
        _merge_map(atom_types, overlay.atom_types, code="duplicate_base_atom_type")
        _merge_map(charges, overlay.charges, code="duplicate_base_charge")
        _merge_map(masses, overlay.masses, code="duplicate_base_mass")
        _merge_map(lj_parameters, overlay.lj_parameters, code="duplicate_base_lj_parameter")
        _merge_map(bonded_parameters, overlay.bonded_parameters, code="duplicate_base_bonded_term")
    if covered != expected_ids:
        missing = sorted(expected_ids - covered)
        extra = sorted(covered - expected_ids)
        raise ValidationError("base_overlay_coverage_mismatch", f"missing={missing},extra={extra}")
    covered_atom_ids = tuple(sorted(covered, key=lambda atom_id: (atom_order[atom_id], atom_id)))
    providers = tuple(fragment.report.provider for fragment in ordered)
    parameter_source = "xponge:base-composition:" + ("+".join(providers) if providers else "none")
    overlay = BaseForceFieldOverlay(
        topology_hash=topology_hash,
        covered_atom_ids=covered_atom_ids,
        atom_types=atom_types,
        charges=charges,
        masses=masses,
        lj_parameters=lj_parameters,
        bonded_parameters=bonded_parameters,
        parameter_source=parameter_source,
    )
    report = BaseCompositionReport(
        schema_version=package.request.schema_version,
        topology_hash=topology_hash,
        projection_hash=projection_hash,
        providers=providers,
        provider_report_hashes=provider_report_hashes,
        covered_atom_ids=covered_atom_ids,
        bonded_term_count=len(bonded_parameters),
        parameter_source=parameter_source,
    ).with_computed_hash()
    return BaseCompositionOutput(overlay, report)


__all__ = ["BaseCompositionOutput", "BaseCompositionReport", "compose_base_force_field"]
