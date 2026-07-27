"""Deterministic composition of bonded metal force-field fit artifacts."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import math
from typing import Any, Mapping, Sequence

from .base_composition import BaseCompositionOutput
from .charge_fit import ModelChargeArtifact, project_model_charge_artifact
from .contracts import (
    BondedParameterOverlay,
    MetalAssignmentInput,
    ParameterizationResult,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _sha256,
    validate_result,
)
from .force_fit import (
    HessianArtifact,
    empirical_registry_bonded_terms,
    empirical_zn_nos_bonded_terms,
    manual_bonded_terms,
    seminario_bonded_terms,
)
from .input import MetalAssignmentPackage, validate_package
from .metal_overlay import build_metal_parameter_overlay
from .partial_charge import compose_partial_charge_artifacts


@dataclass(frozen=True, slots=True)
class MetalAtomParameterSpec:
    external_id: str
    element: str
    formal_charge: int
    atom_type: str
    partial_charge: float
    mass: float
    epsilon: float
    rmin: float


@dataclass(frozen=True, slots=True)
class BondedMetalParameterSpec:
    schema_version: int
    topology_hash: str
    metal_atoms: tuple[MetalAtomParameterSpec, ...]
    donor_atom_types: Mapping[str, str]
    parameter_source: str
    provenance: Mapping[str, Any]
    precedence: int = 100
    spec_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "donor_atom_types")
        _freeze_field(self, "provenance")

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "spec_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "BondedMetalParameterSpec":
        return replace(self, spec_hash=self.computed_hash())


def validate_bonded_metal_parameter_spec(
    request: MetalAssignmentInput,
    spec: BondedMetalParameterSpec,
) -> None:
    topology = request.topology
    if spec.schema_version != request.schema_version:
        raise ValidationError("unsupported_schema_version", str(spec.schema_version), "metal_parameter_spec")
    if spec.topology_hash != topology.topology_hash:
        raise ValidationError("metal_parameter_spec_topology_mismatch", spec.topology_hash)
    if not spec.spec_hash or spec.spec_hash != spec.computed_hash():
        raise ValidationError("stale_metal_parameter_spec_hash", "spec hash mismatch")
    if not spec.parameter_source or not spec.provenance:
        raise ValidationError("missing_metal_parameter_provenance", spec.parameter_source)
    if isinstance(spec.precedence, bool) or not isinstance(spec.precedence, int) or spec.precedence <= 0:
        raise ValidationError("invalid_overlay_precedence", str(spec.precedence))
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    expected_metal_ids = {atom_id for atom_id, atom in atom_by_id.items() if atom.is_metal}
    actual_ids = [atom.external_id for atom in spec.metal_atoms]
    if len(actual_ids) != len(set(actual_ids)) or set(actual_ids) != expected_metal_ids:
        raise ValidationError("metal_parameter_spec_coverage_mismatch", ",".join(actual_ids))
    for parameter in spec.metal_atoms:
        path = f"metal_parameter_spec.{parameter.external_id}"
        atom = atom_by_id[parameter.external_id]
        if parameter.element != atom.element or parameter.formal_charge != atom.formal_charge:
            raise ValidationError("metal_parameter_spec_identity_mismatch", parameter.external_id, path)
        if not parameter.atom_type:
            raise ValidationError("missing_metal_atom_type", parameter.external_id, path)
        for name in ("partial_charge", "mass", "epsilon", "rmin"):
            value = getattr(parameter, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValidationError("invalid_metal_parameter", name, path)
        if parameter.mass <= 0 or parameter.epsilon < 0 or parameter.rmin < 0:
            raise ValidationError("invalid_metal_parameter", parameter.external_id, path)
    allowed_donors = {
        atom_id for atom_id, atom in atom_by_id.items()
        if not atom.is_metal and "metalOverlay" in atom.scopes
    }
    if not set(spec.donor_atom_types) <= allowed_donors:
        raise ValidationError("metal_parameter_spec_donor_scope_mismatch", "donor type scope")
    if any(not isinstance(value, str) or not value for value in spec.donor_atom_types.values()):
        raise ValidationError("invalid_donor_atom_type", "donor atom types must be non-empty")


def _validate_complete_base(
    package: MetalAssignmentPackage,
    base: BaseCompositionOutput,
) -> None:
    report = base.report
    if (
        report.topology_hash != package.request.topology.topology_hash
        or report.projection_hash != package.request.projection_hash
        or not report.report_hash
        or report.report_hash != report.computed_hash()
    ):
        raise ValidationError("invalid_base_composition_report", report.parameter_source)


def _merge_terms(target: dict[str, Mapping[str, Any]], source: Mapping[str, Mapping[str, Any]]) -> None:
    overlap = set(target) & set(source)
    if overlap:
        raise ValidationError("conflicting_bonded_fit_term", ",".join(sorted(overlap)))
    target.update(source)


def compose_bonded_fit(
    package: MetalAssignmentPackage,
    base: BaseCompositionOutput,
    metal_spec: BondedMetalParameterSpec,
    *,
    charge_artifacts: Sequence[ModelChargeArtifact] = (),
    hessian_artifacts: Sequence[HessianArtifact] = (),
    force_method: str = "seminario",
    scale_factor: float = 1.0,
    empirical_registry_id: str = "",
    empirical_geometry: str = "unclassified",
    empirical_base_force_field: str = "",
    empirical_water_model: str = "",
    manual_bond_force_constant: float | None = None,
    manual_angle_force_constant: float | None = None,
    manual_site_force_constants: Mapping[str, Any] | None = None,
    reference_geometry_artifact: Mapping[str, Any] | None = None,
) -> ParameterizationResult:
    """Compose immutable base, charge, metal and bonded overlays for a bonded request."""

    validate_package(package)
    request = package.request
    if request.interaction_model != "bonded":
        raise ValidationError("invalid_bonded_fit_interaction_model", request.interaction_model)
    _validate_complete_base(package, base)
    validate_bonded_metal_parameter_spec(request, metal_spec)

    model_by_id = {model.external_id: model for model in package.prepared_artifacts.derived_models.models}
    projected_charge_artifacts = ()
    charge_projection_report: Mapping[str, Any] = {"status": "not_requested"}
    if request.charge_contract is None:
        if charge_artifacts:
            raise ValidationError("unexpected_charge_artifact", "request has no charge contract")
    else:
        large_models = {model.external_id: model for model in model_by_id.values() if model.purpose == "large"}
        artifact_by_id = {artifact.model_id: artifact for artifact in charge_artifacts}
        if len(artifact_by_id) != len(charge_artifacts) or set(artifact_by_id) != set(large_models):
            raise ValidationError("charge_artifact_model_coverage_mismatch", "large model coverage must be exact")
        projected_artifact, charge_projection_report = project_model_charge_artifact(
            request.topology,
            tuple((large_models[model_id], artifact_by_id[model_id]) for model_id in sorted(large_models)),
            request.charge_contract,
        )
        projected_charge_artifacts = (projected_artifact,)
    charge_overlay, charge_composition_report = compose_partial_charge_artifacts(
        request,
        projected_charge_artifacts,
    )
    charge_report = {
        "projection": _canonicalize(charge_projection_report),
        "composition": _canonicalize(charge_composition_report),
    }

    bonded_terms: dict[str, Mapping[str, Any]] = {}
    bonded_reports: list[Mapping[str, Any]] = []
    normalized_force_method = force_method.lower()
    if normalized_force_method == "seminario":
        small_models = {model.external_id: model for model in model_by_id.values() if model.purpose == "small"}
        artifact_by_id = {artifact.model_id: artifact for artifact in hessian_artifacts}
        if len(artifact_by_id) != len(hessian_artifacts) or set(artifact_by_id) != set(small_models):
            raise ValidationError("hessian_artifact_model_coverage_mismatch", "small model coverage must be exact")
        for model_id in sorted(small_models):
            terms, report = seminario_bonded_terms(
                small_models[model_id],
                artifact_by_id[model_id],
                scale_factor=scale_factor,
            )
            _merge_terms(bonded_terms, terms)
            bonded_reports.append(report)
    elif normalized_force_method in {"empirical_registry", "empirical_zn_nos"}:
        if hessian_artifacts:
            raise ValidationError("unexpected_hessian_artifact", "empirical force method does not consume Hessians")
        if normalized_force_method == "empirical_registry":
            if not empirical_registry_id:
                raise ValidationError("missing_empirical_registry_id", normalized_force_method)
            terms, report = empirical_registry_bonded_terms(
                request.topology,
                registry_id=empirical_registry_id,
                geometry=empirical_geometry,
                base_force_field=empirical_base_force_field,
                water_model=empirical_water_model,
            )
        else:
            terms, report = empirical_zn_nos_bonded_terms(request.topology)
        _merge_terms(bonded_terms, terms)
        bonded_reports.append(report)
    elif normalized_force_method == "manual_bonded":
        if hessian_artifacts:
            raise ValidationError(
                "unexpected_hessian_artifact",
                "manual_bonded does not consume Hessians",
            )
        terms, report = manual_bonded_terms(
            request.topology,
            bond_force_constant=manual_bond_force_constant,
            angle_force_constant=manual_angle_force_constant,
            site_force_constants=manual_site_force_constants,
            reference_geometry_artifact=reference_geometry_artifact,
        )
        _merge_terms(bonded_terms, terms)
        bonded_reports.append(report)
    else:
        raise ValidationError("unsupported_bonded_force_method", force_method)
    if not bonded_terms:
        raise ValidationError("missing_bonded_fit_terms", normalized_force_method)

    metal_by_id = {parameter.external_id: parameter for parameter in metal_spec.metal_atoms}
    stable_order = {atom.external_id: atom.stable_order for atom in request.topology.atoms}
    coverage = set(metal_by_id) | set(metal_spec.donor_atom_types)
    covered_atom_ids = tuple(sorted(coverage, key=lambda atom_id: (stable_order[atom_id], atom_id)))
    atom_types = {atom_id: parameter.atom_type for atom_id, parameter in metal_by_id.items()}
    atom_types.update(metal_spec.donor_atom_types)
    metal_output = build_metal_parameter_overlay(
        request,
        covered_atom_ids=covered_atom_ids,
        atom_types=atom_types,
        charges={atom_id: parameter.partial_charge for atom_id, parameter in metal_by_id.items()},
        masses={atom_id: parameter.mass for atom_id, parameter in metal_by_id.items()},
        lj_parameters={
            atom_id: {
                "epsilon": parameter.epsilon,
                "rmin": parameter.rmin,
                "energy_unit": "kcal/mol",
                "length_unit": "angstrom",
                "source": metal_spec.parameter_source,
            }
            for atom_id, parameter in metal_by_id.items()
        },
        bonded_parameters={},
        parameter_source=metal_spec.parameter_source,
        precedence=metal_spec.precedence,
        require_metal_charges=True,
        provenance={
            "metal_parameter_spec_hash": metal_spec.spec_hash,
            "spec_provenance": metal_spec.provenance,
        },
    )
    bonded_overlay = BondedParameterOverlay(
        topology_hash=request.topology.topology_hash,
        terms=bonded_terms,
        source="bonded-fit:" + normalized_force_method,
    )
    empirical_provenance = {}
    if bonded_reports and normalized_force_method in {"empirical_registry", "empirical_zn_nos"}:
        empirical_provenance = {
            "empirical_registry": _canonicalize(bonded_reports[0].get("registry", {})),
        }
    result = ParameterizationResult(
        schema_version=request.schema_version,
        request_id=request.request_id,
        input_hash=request.input_hash,
        graph_revision=request.graph_revision,
        topology_hash=request.topology.topology_hash,
        projection_hash=request.projection_hash,
        base_overlay=base.overlay,
        metal_overlay=metal_output.overlay,
        charge_overlay=charge_overlay,
        bonded_overlay=bonded_overlay,
        fit_reports={
            "base_composition": _canonicalize(base.report),
            "metal_parameters": _canonicalize(metal_output.report),
            "charge_fit": charge_report,
            "bonded_fit": tuple(bonded_reports),
        },
        provenance={
            "package_hash": package.package_hash,
            "metal_parameter_spec_hash": metal_spec.spec_hash,
            "force_method": normalized_force_method,
            **empirical_provenance,
        },
        status="overlay_validated",
        complete=False,
    ).with_computed_hash()
    validate_result(request, result)
    return result


__all__ = [
    "BondedMetalParameterSpec", "MetalAtomParameterSpec", "compose_bonded_fit",
    "validate_bonded_metal_parameter_spec",
]
