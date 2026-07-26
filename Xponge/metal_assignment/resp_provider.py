"""Parent-side isolated RESP provider for Chemcore-derived large models."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
import math
import subprocess
from typing import Any, Mapping

from ._worker_runtime import run_worker_subprocess, worker_command
from .artifacts import DerivedModel
from .charge_fit import ModelChargeArtifact, model_charge_artifact_from_dict, validate_model_charge_artifact
from .contracts import ValidationError, _canonicalize, _freeze, _sha256
from .input import MetalAssignmentPackage, validate_package
from .resp_input import RespFitInput, validate_resp_fit_input
from ..qm.resp_basis import ResolvedRespBasis, resolve_resp_basis


RESP_WORKER_PROTOCOL_VERSION = 1


@dataclass(frozen=True, slots=True)
class RespFitOutput:
    schema_version: int
    fit_input_hash: str
    model_artifacts: tuple[ModelChargeArtifact, ...]
    fit_reports: tuple[Mapping[str, Any], ...]
    worker_response_hashes: tuple[str, ...]
    output_hash: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "fit_reports", tuple(_freeze(report) for report in self.fit_reports))

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "output_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "RespFitOutput":
        return replace(self, output_hash=self.computed_hash())


def _model_payload(
    model: DerivedModel,
    fit_input: RespFitInput,
    metal_atom_ids: frozenset[str],
) -> dict[str, Any]:
    state = model.electronic_state
    if state is None:
        raise ValidationError("missing_model_electronic_state", model.external_id)
    atom_order = tuple(atom.model_atom_id for atom in model.atoms)
    known_atom_ids = set(atom_order)
    proof = model.charge_accounting
    constraint_groups = proof.get("constraint_groups")
    manifest = model.capped_model_manifest
    if (
        manifest is None
        or not isinstance(constraint_groups, list)
        or proof.get("complete") is not True
        or proof.get("model_target_charge") != state.net_charge
    ):
        raise ValidationError("missing_constrained_charge_ledger", model.external_id)
    linear_constraints = []
    constraint_ids: set[str] = set()
    for index, group in enumerate(constraint_groups):
        if not isinstance(group, Mapping):
            raise ValidationError("invalid_charge_constraint_group", str(index), model.external_id)
        atom_ids = group.get("model_atom_ids")
        if (
            not isinstance(atom_ids, list)
            or not atom_ids
            or len(atom_ids) != len(set(atom_ids))
            or not set(atom_ids) <= known_atom_ids
            or not isinstance(group.get("target_charge"), int)
            or isinstance(group.get("target_charge"), bool)
            or not isinstance(group.get("group_id"), str)
            or not group["group_id"]
            or group.get("role") not in {"core", "cap"}
            or not isinstance(group.get("source"), str)
            or not group["source"]
            or group.get("complete") is not True
        ):
            raise ValidationError("invalid_charge_constraint_group", str(index), model.external_id)
        constraint_id = group["group_id"]
        constraint_ids.add(constraint_id)
        linear_constraints.append({
            "constraint_id": constraint_id,
            "role": group.get("role"),
            "atom_ids": atom_ids,
            "coefficients": [1.0] * len(atom_ids),
            "target_charge": float(group["target_charge"]),
            "source": group.get("source"),
        })
    for constraint in fit_input.linear_constraints:
        if constraint.model_id != model.external_id:
            continue
        if constraint.constraint_id in constraint_ids:
            raise ValidationError(
                "duplicate_resp_constraint_id",
                constraint.constraint_id,
                model.external_id,
            )
        if not set(constraint.atom_ids) <= known_atom_ids:
            raise ValidationError(
                "resp_constraint_atom_mismatch",
                constraint.constraint_id,
                model.external_id,
            )
        constraint_ids.add(constraint.constraint_id)
        linear_constraints.append({
            "constraint_id": constraint.constraint_id,
            "role": constraint.role,
            "atom_ids": list(constraint.atom_ids),
            "coefficients": [float(value) for value in constraint.coefficients],
            "target_charge": float(constraint.target_charge),
            "source": constraint.source,
        })
    return {
        "external_id": model.external_id,
        "model_hash": model.model_hash,
        "purpose": model.purpose,
        "coordinate_unit": model.coordinate_unit,
        "atomic_charge_role": model.atomic_charge_role,
        "electronic_state": _canonicalize(state),
        "atoms": [
            {
                "model_atom_id": atom.model_atom_id,
                "element": atom.element,
                "coordinates": list(atom.coordinates),
                "initial_charge": atom.partial_charge,
                "is_metal": atom.external_id in metal_atom_ids,
            }
            for atom in model.atoms
        ],
        "bonds": [
            {"model_atom_ids": list(bond.model_atom_ids), "order": bond.order}
            for bond in model.bonds
        ],
        "links": [
            {"model_atom_ids": list(link.model_atom_ids)}
            for link in model.links
        ],
        "linear_constraints": linear_constraints,
    }


def _fit_protocol_payload(value: RespFitInput, model_id: str) -> dict[str, Any]:
    return {
        "fit_input_hash": value.fit_input_hash,
        "backend": value.backend,
        "basis_family": value.basis_family,
        "metal_basis_policy": value.metal_basis_policy,
        "basis_source": value.basis_source,
        "optimize_geometry": value.optimize_geometry,
        "grid_density": value.grid_density,
        "grid_cell_layer": value.grid_cell_layer,
        "radius_overrides": dict(value.radius_overrides),
        "restraint_a1": value.restraint_a1,
        "restraint_a2": value.restraint_a2,
        "two_stage": value.two_stage,
        "only_esp": value.only_esp,
        "esp_memory_limit_bytes": value.esp_memory_limit_bytes,
        "esp_chunk_policy": value.esp_chunk_policy,
        "esp_safety_factor": value.esp_safety_factor,
        "equivalence_groups": [
            _canonicalize(group) for group in value.equivalence_groups if group.model_id == model_id
        ],
        "source": value.source,
    }


def _model_basis_context(
    model: DerivedModel,
    fit_input: RespFitInput,
    metal_atom_ids: frozenset[str],
) -> tuple[ResolvedRespBasis, tuple[str, ...]]:
    elements = tuple(atom.element for atom in model.atoms)
    metal_elements = tuple(sorted({
        atom.element for atom in model.atoms if atom.external_id in metal_atom_ids
    }))
    try:
        resolved = resolve_resp_basis(fit_input.basis_family, elements)
    except ValueError as exc:
        raise ValidationError(
            "invalid_resp_basis_mapping", str(exc), f"large_model.{model.external_id}"
        ) from exc
    if metal_elements and fit_input.metal_basis_policy == "require_ecp":
        ecp_elements = set(resolved.ecp or {})
        missing = sorted(set(metal_elements) - ecp_elements)
        if missing:
            raise ValidationError(
                "resp_metal_ecp_required",
                f"basis family {fit_input.basis_family!r} has no ECP for {','.join(missing)}",
                f"large_model.{model.external_id}",
            )
    if resolved.ecp and fit_input.backend == "psi4":
        raise ValidationError(
            "resp_backend_basis_unsupported",
            "Psi4 cannot consume element-mapped ECP RESP basis families; use PySCF",
            f"large_model.{model.external_id}",
        )
    return resolved, metal_elements


def _invoke_resp_worker(
    payload: Mapping[str, Any],
    *,
    timeout_seconds: float,
    python_executable: str | None,
) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        command = worker_command(
            "Xponge.metal_assignment._resp_worker",
            python_executable=python_executable,
        )
    except (FileNotFoundError, OSError) as exc:
        raise ValidationError("invalid_resp_python_executable", str(exc), "python_executable") from exc
    try:
        completed = run_worker_subprocess(
            command,
            input_text=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("resp_worker_timeout", str(timeout_seconds)) from exc
    except OSError as exc:
        raise ValidationError("resp_worker_launch_failed", str(exc)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError(
            "invalid_resp_worker_output", completed.stderr.strip()[-1000:] or str(exc)
        ) from exc
    if not isinstance(response, dict) or response.get("ok") is not True or completed.returncode != 0:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "RESP worker failed") if isinstance(error, dict) else "RESP worker failed"
        if completed.stderr.strip():
            message = f"{message}; stderr={completed.stderr.strip()[-1000:]}"
        raise ValidationError("resp_worker_failed", message)
    return response


def _validate_worker_response(
    response: Mapping[str, Any],
    *,
    payload: Mapping[str, Any],
    model: DerivedModel,
    fit_input: RespFitInput,
    resolved_basis: ResolvedRespBasis,
    metal_elements: tuple[str, ...],
) -> tuple[ModelChargeArtifact, Mapping[str, Any], str]:
    required = {
        "ok", "protocol_version", "model_id", "model_hash", "fit_input_hash", "request_hash",
        "artifact", "fit_report", "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_resp_worker_response", "unexpected fields")
    response_payload = dict(response)
    response_hash = response_payload.pop("response_hash")
    if not isinstance(response_hash, str) or response_hash != _sha256(response_payload):
        raise ValidationError("stale_resp_worker_hash", model.external_id)
    request_hash = _sha256(payload)
    if (
        response["protocol_version"] != RESP_WORKER_PROTOCOL_VERSION
        or response["model_id"] != model.external_id
        or response["model_hash"] != model.model_hash
        or response["fit_input_hash"] != fit_input.fit_input_hash
        or response["request_hash"] != request_hash
    ):
        raise ValidationError("resp_worker_identity_mismatch", model.external_id)
    artifact = model_charge_artifact_from_dict(
        response["artifact"], path=f"resp_worker.{model.external_id}.artifact"
    )
    validate_model_charge_artifact(model, artifact)
    report = response["fit_report"]
    if not isinstance(report, Mapping):
        raise ValidationError("invalid_resp_fit_report", model.external_id)
    expected_report = {
        "model_id": model.external_id,
        "model_hash": model.model_hash,
        "fit_input_hash": fit_input.fit_input_hash,
        "request_hash": request_hash,
        "backend": fit_input.backend,
        "basis_family": fit_input.basis_family,
        "metal_basis_policy": fit_input.metal_basis_policy,
        "basis_source": fit_input.basis_source,
        "metal_elements": list(metal_elements),
        "coordinate_unit": "angstrom",
        "geometry_locked": True,
    }
    for name, expected in expected_report.items():
        if report.get(name) != expected:
            raise ValidationError("resp_fit_report_identity_mismatch", name, model.external_id)
    residual = report.get("charge_residual")
    if isinstance(residual, bool) or not isinstance(residual, (int, float)) or not math.isfinite(residual):
        raise ValidationError("invalid_resp_charge_residual", model.external_id)
    if abs(residual) > 1.0e-6:
        raise ValidationError("resp_charge_not_conserved", repr(residual), model.external_id)
    constraint_diagnostics = report.get("constraint_diagnostics")
    if not isinstance(constraint_diagnostics, Mapping):
        raise ValidationError("missing_resp_constraint_diagnostics", model.external_id)
    expected_constraint_count = (
        len(model.charge_accounting["constraint_groups"])
        + sum(
            1
            for constraint in fit_input.linear_constraints
            if constraint.model_id == model.external_id
        )
    )
    expected_input_constraint_count = (
        1
        + expected_constraint_count
        + sum(
            len(group.atom_ids) - 1
            for group in fit_input.equivalence_groups
            if group.model_id == model.external_id
        )
    )
    constraint_count = report.get("constraint_count")
    if (
        isinstance(constraint_count, bool)
        or not isinstance(constraint_count, int)
        or constraint_count != expected_constraint_count
        or constraint_diagnostics.get("input_constraint_count") != expected_input_constraint_count
    ):
        raise ValidationError("resp_constraint_count_mismatch", model.external_id)
    constraint_residual = constraint_diagnostics.get("max_constraint_residual")
    if (
        isinstance(constraint_residual, bool)
        or not isinstance(constraint_residual, (int, float))
        or not math.isfinite(constraint_residual)
        or abs(constraint_residual) > 1.0e-8
    ):
        raise ValidationError(
            "resp_constraints_not_satisfied",
            repr(constraint_residual),
            model.external_id,
        )
    for name in (
        "esp_rmse_au",
        "esp_relative_rmse",
        "esp_mae_au",
        "esp_max_abs_error_au",
    ):
        value = constraint_diagnostics.get(name)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value < 0.0
        ):
            raise ValidationError("invalid_resp_esp_fit_diagnostic", name, model.external_id)
    point_count = constraint_diagnostics.get("esp_point_count")
    if isinstance(point_count, bool) or not isinstance(point_count, int) or point_count <= 0:
        raise ValidationError("invalid_resp_esp_fit_diagnostic", "esp_point_count", model.external_id)
    metadata = report.get("setup_metadata")
    if not isinstance(metadata, Mapping) or metadata.get("method") != "RESP" or metadata.get("scf_converged") is not True:
        raise ValidationError("incomplete_resp_provenance", model.external_id)
    expected_setup = {
        "basis_family": resolved_basis.label,
        "basis": _canonicalize(resolved_basis.basis),
        "ecp": _canonicalize(resolved_basis.ecp),
        "cart": resolved_basis.cart,
    }
    for name, expected in expected_setup.items():
        if _canonicalize(metadata.get(name)) != expected:
            raise ValidationError("resp_resolved_basis_mismatch", name, model.external_id)
    return artifact, report, response_hash


def fit_resp_charges(
    package: MetalAssignmentPackage,
    fit_input: RespFitInput,
    *,
    timeout_seconds: float = 600.0,
    python_executable: str | None = None,
) -> RespFitOutput:
    """Fit every locked large model in isolated Xponge worker processes."""

    validate_package(package)
    validate_resp_fit_input(fit_input)
    request = package.request
    if request.interaction_model != "bonded":
        raise ValidationError("resp_interaction_model_mismatch", request.interaction_model)
    if request.charge_contract is None:
        raise ValidationError("missing_charge_contract", "RESP requires an explicit charge contract")
    models = package.prepared_artifacts.derived_models
    large_models = tuple(model for model in models.models if model.purpose == "large")
    if not large_models:
        raise ValidationError("missing_large_model", "RESP requires at least one large model")
    site_ids = {site.external_id for site in models.sites}
    large_site_ids = [model.site_id for model in large_models]
    if set(large_site_ids) != site_ids or len(large_site_ids) != len(set(large_site_ids)):
        raise ValidationError("invalid_large_model_site_coverage", "expected exactly one large model per site")
    large_by_id = {model.external_id: model for model in large_models}
    for group in fit_input.equivalence_groups:
        model = large_by_id.get(group.model_id)
        if model is None:
            raise ValidationError("resp_equivalence_model_mismatch", group.model_id)
        model_atom_ids = {atom.model_atom_id for atom in model.atoms}
        if not set(group.atom_ids) <= model_atom_ids:
            raise ValidationError("resp_equivalence_atom_mismatch", group.model_id)
    for constraint in fit_input.linear_constraints:
        model = large_by_id.get(constraint.model_id)
        if model is None:
            raise ValidationError("resp_constraint_model_mismatch", constraint.model_id)
        model_atom_ids = {atom.model_atom_id for atom in model.atoms}
        if not set(constraint.atom_ids) <= model_atom_ids:
            raise ValidationError(
                "resp_constraint_atom_mismatch",
                constraint.constraint_id,
                constraint.model_id,
            )
    artifacts: list[ModelChargeArtifact] = []
    reports: list[Mapping[str, Any]] = []
    response_hashes: list[str] = []
    metal_atom_ids = frozenset(atom.external_id for atom in request.topology.atoms if atom.is_metal)
    for model in sorted(large_models, key=lambda item: (item.site_id, item.external_id)):
        resolved_basis, metal_elements = _model_basis_context(model, fit_input, metal_atom_ids)
        payload = {
            "protocol_version": RESP_WORKER_PROTOCOL_VERSION,
            "model": _model_payload(model, fit_input, metal_atom_ids),
            "fit_protocol": _fit_protocol_payload(fit_input, model.external_id),
        }
        response = _invoke_resp_worker(
            payload,
            timeout_seconds=timeout_seconds,
            python_executable=python_executable,
        )
        artifact, report, response_hash = _validate_worker_response(
            response,
            payload=payload,
            model=model,
            fit_input=fit_input,
            resolved_basis=resolved_basis,
            metal_elements=metal_elements,
        )
        artifacts.append(artifact)
        reports.append(report)
        response_hashes.append(response_hash)
    return RespFitOutput(
        schema_version=RESP_WORKER_PROTOCOL_VERSION,
        fit_input_hash=fit_input.fit_input_hash,
        model_artifacts=tuple(artifacts),
        fit_reports=tuple(reports),
        worker_response_hashes=tuple(response_hashes),
    ).with_computed_hash()


__all__ = ["RESP_WORKER_PROTOCOL_VERSION", "RespFitOutput", "fit_resp_charges"]
