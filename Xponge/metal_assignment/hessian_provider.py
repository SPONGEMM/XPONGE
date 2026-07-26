"""Parent-side isolated Hessian provider for Chemcore-derived small models."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
import math
import subprocess
from typing import Any, Mapping

from ._worker_runtime import run_worker_subprocess, worker_command
from .artifacts import DerivedModel
from .contracts import ValidationError, _canonicalize, _freeze, _sha256
from .force_fit import HessianArtifact, hessian_artifact_from_dict, validate_hessian_artifact
from .hessian_input import HessianFitInput, validate_hessian_fit_input
from .input import MetalAssignmentPackage, validate_package
from ..qm.resp_basis import ResolvedRespBasis, resolve_resp_basis


HESSIAN_WORKER_PROTOCOL_VERSION = 1


@dataclass(frozen=True, slots=True)
class HessianFitOutput:
    schema_version: int
    fit_input_hash: str
    hessian_artifacts: tuple[HessianArtifact, ...]
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

    def with_computed_hash(self) -> "HessianFitOutput":
        return replace(self, output_hash=self.computed_hash())


def _model_payload(model: DerivedModel, metal_atom_ids: frozenset[str]) -> dict[str, Any]:
    state = model.electronic_state
    if state is None:
        raise ValidationError("missing_model_electronic_state", model.external_id)
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
        "linear_constraints": [],
    }


def _fit_protocol_payload(value: HessianFitInput) -> dict[str, Any]:
    return {
        "fit_input_hash": value.fit_input_hash,
        "backend": value.backend,
        "method": value.method,
        "basis_family": value.basis_family,
        "metal_basis_policy": value.metal_basis_policy,
        "basis_source": value.basis_source,
        "optimize_geometry": value.optimize_geometry,
        "scf_convergence_tolerance": value.scf_convergence_tolerance,
        "scf_max_cycles": value.scf_max_cycles,
        "threads": value.threads,
        "memory_limit_bytes": value.memory_limit_bytes,
        "source": value.source,
    }


def _model_basis_context(
    model: DerivedModel,
    fit_input: HessianFitInput,
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
            "invalid_hessian_basis_mapping", str(exc), f"small_model.{model.external_id}"
        ) from exc
    if metal_elements and fit_input.metal_basis_policy == "require_ecp":
        missing = sorted(set(metal_elements) - set(resolved.ecp or {}))
        if missing:
            raise ValidationError(
                "hessian_metal_ecp_required",
                f"basis family {fit_input.basis_family!r} has no ECP for {','.join(missing)}",
                f"small_model.{model.external_id}",
            )
    return resolved, metal_elements


def _invoke_hessian_worker(
    payload: Mapping[str, Any],
    *,
    timeout_seconds: float,
    python_executable: str | None,
) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        command = worker_command(
            "Xponge.metal_assignment._hessian_worker",
            python_executable=python_executable,
        )
    except (FileNotFoundError, OSError) as exc:
        raise ValidationError("invalid_hessian_python_executable", str(exc), "python_executable") from exc
    try:
        completed = run_worker_subprocess(
            command,
            input_text=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("hessian_worker_timeout", str(timeout_seconds)) from exc
    except OSError as exc:
        raise ValidationError("hessian_worker_launch_failed", str(exc)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise ValidationError(
            "invalid_hessian_worker_output", completed.stderr.strip()[-1000:] or str(exc)
        ) from exc
    if not isinstance(response, dict) or response.get("ok") is not True or completed.returncode != 0:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "Hessian worker failed") if isinstance(error, dict) else "Hessian worker failed"
        if completed.stderr.strip():
            message = f"{message}; stderr={completed.stderr.strip()[-1000:]}"
        raise ValidationError("hessian_worker_failed", message)
    return response


def _validate_worker_response(
    response: Mapping[str, Any],
    *,
    payload: Mapping[str, Any],
    model: DerivedModel,
    fit_input: HessianFitInput,
    resolved_basis: ResolvedRespBasis,
    metal_elements: tuple[str, ...],
) -> tuple[HessianArtifact, Mapping[str, Any], str]:
    required = {
        "ok", "protocol_version", "model_id", "model_hash", "fit_input_hash",
        "request_hash", "artifact", "fit_report", "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_hessian_worker_response", "unexpected fields")
    response_payload = dict(response)
    response_hash = response_payload.pop("response_hash")
    if not isinstance(response_hash, str) or response_hash != _sha256(response_payload):
        raise ValidationError("stale_hessian_worker_hash", model.external_id)
    request_hash = _sha256(payload)
    if (
        response["protocol_version"] != HESSIAN_WORKER_PROTOCOL_VERSION
        or response["model_id"] != model.external_id
        or response["model_hash"] != model.model_hash
        or response["fit_input_hash"] != fit_input.fit_input_hash
        or response["request_hash"] != request_hash
    ):
        raise ValidationError("hessian_worker_identity_mismatch", model.external_id)
    artifact = hessian_artifact_from_dict(
        response["artifact"], path=f"hessian_worker.{model.external_id}.artifact"
    )
    validate_hessian_artifact(model, artifact)
    report = response["fit_report"]
    if not isinstance(report, Mapping):
        raise ValidationError("invalid_hessian_fit_report", model.external_id)
    expected_report = {
        "model_id": model.external_id,
        "model_hash": model.model_hash,
        "fit_input_hash": fit_input.fit_input_hash,
        "request_hash": request_hash,
        "backend": fit_input.backend,
        "method": fit_input.method,
        "basis_family": fit_input.basis_family,
        "metal_basis_policy": fit_input.metal_basis_policy,
        "basis_source": fit_input.basis_source,
        "metal_elements": list(metal_elements),
        "coordinate_unit": "angstrom",
        "geometry_locked": True,
    }
    for name, expected in expected_report.items():
        if report.get(name) != expected:
            raise ValidationError("hessian_fit_report_identity_mismatch", name, model.external_id)
    metadata = report.get("setup_metadata")
    if (
        not isinstance(metadata, Mapping)
        or metadata.get("method") != "SCF_HESSIAN"
        or metadata.get("scf_converged") is not True
    ):
        raise ValidationError("incomplete_hessian_provenance", model.external_id)
    expected_setup = {
        "basis_family": resolved_basis.label,
        "basis": _canonicalize(resolved_basis.basis),
        "ecp": _canonicalize(resolved_basis.ecp),
        "cart": resolved_basis.cart,
        "references": list(resolved_basis.references),
        "scf_convergence_tolerance": float(fit_input.scf_convergence_tolerance),
        "scf_max_cycles": fit_input.scf_max_cycles,
        "threads": fit_input.threads,
        "memory_limit_bytes": fit_input.memory_limit_bytes,
    }
    for name, expected in expected_setup.items():
        if _canonicalize(metadata.get(name)) != expected:
            raise ValidationError("hessian_resolved_setup_mismatch", name, model.external_id)
    return artifact, report, response_hash


def fit_model_hessians(
    package: MetalAssignmentPackage,
    fit_input: HessianFitInput,
    *,
    timeout_seconds: float = 1800.0,
    python_executable: str | None = None,
) -> HessianFitOutput:
    """Compute every locked small-model Hessian in isolated worker processes."""

    validate_package(package)
    validate_hessian_fit_input(fit_input)
    small_models = tuple(
        model for model in package.prepared_artifacts.derived_models.models
        if model.purpose == "small"
    )
    if not small_models:
        raise ValidationError("missing_small_models", "Hessian fitting requires at least one small model")
    request = package.request
    metal_atom_ids = frozenset(atom.external_id for atom in request.topology.atoms if atom.is_metal)
    artifacts: list[HessianArtifact] = []
    reports: list[Mapping[str, Any]] = []
    response_hashes: list[str] = []
    for model in sorted(small_models, key=lambda item: (item.site_id, item.external_id)):
        resolved_basis, metal_elements = _model_basis_context(model, fit_input, metal_atom_ids)
        payload = {
            "protocol_version": HESSIAN_WORKER_PROTOCOL_VERSION,
            "model": _model_payload(model, metal_atom_ids),
            "fit_protocol": _fit_protocol_payload(fit_input),
        }
        response = _invoke_hessian_worker(
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
    return HessianFitOutput(
        schema_version=HESSIAN_WORKER_PROTOCOL_VERSION,
        fit_input_hash=fit_input.fit_input_hash,
        hessian_artifacts=tuple(artifacts),
        fit_reports=tuple(reports),
        worker_response_hashes=tuple(response_hashes),
    ).with_computed_hash()


__all__ = [
    "HESSIAN_WORKER_PROTOCOL_VERSION", "HessianFitOutput", "fit_model_hessians",
]
