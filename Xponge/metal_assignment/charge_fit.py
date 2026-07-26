"""Charge-fit projection and constraint enforcement for prepared models."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import math
from typing import Any, Mapping, Sequence

import numpy as np

from .artifacts import DerivedModel
from .contracts import (
    ChargeAssignmentContract,
    ChargeOverlay,
    ChargePolicy,
    PartialChargeArtifact,
    PreparedChemicalTopology,
    ValidationError,
    _canonicalize,
    _sha256,
    _strict_object,
    validate_charge_contract,
)


@dataclass(frozen=True, slots=True)
class ModelChargeArtifact:
    model_id: str
    model_hash: str
    atom_order: tuple[str, ...]
    charges: tuple[float, ...]
    atomic_charge_role: str
    provider: str
    provider_version: str
    artifact_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "artifact_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "ModelChargeArtifact":
        return replace(self, artifact_hash=self.computed_hash())


def validate_model_charge_artifact(model: DerivedModel, artifact: ModelChargeArtifact) -> None:
    if model.purpose != "large":
        raise ValidationError("invalid_charge_model_purpose", model.purpose)
    if artifact.model_id != model.external_id or artifact.model_hash != model.model_hash:
        raise ValidationError("charge_model_identity_mismatch", artifact.model_id)
    if artifact.atom_order != tuple(atom.model_atom_id for atom in model.atoms):
        raise ValidationError("charge_atom_order_mismatch", artifact.model_id)
    if len(artifact.charges) != len(model.atoms):
        raise ValidationError("charge_count_mismatch", artifact.model_id)
    if artifact.atomic_charge_role not in {"reference", "fitted"}:
        raise ValidationError("invalid_charge_artifact_role", artifact.atomic_charge_role)
    if not artifact.provider or not artifact.provider_version:
        raise ValidationError("missing_charge_provenance", artifact.model_id)
    if not artifact.artifact_hash or artifact.artifact_hash != artifact.computed_hash():
        raise ValidationError("stale_model_charge_artifact_hash", artifact.model_id)
    if any(
        isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)
        for value in artifact.charges
    ):
        raise ValidationError("invalid_fitted_charge", artifact.model_id)


def model_charge_artifact_from_dict(value: Any, *, path: str = "model_charge_artifact") -> ModelChargeArtifact:
    """Parse a strict worker/file representation of a model charge artifact."""

    data = _strict_object(
        value,
        required={
            "model_id", "model_hash", "atom_order", "charges", "atomic_charge_role", "provider",
            "provider_version", "artifact_hash",
        },
        path=path,
    )
    if not isinstance(data["atom_order"], list) or not all(isinstance(item, str) for item in data["atom_order"]):
        raise ValidationError("invalid_wire_type", "expected string array", f"{path}.atom_order")
    if not isinstance(data["charges"], list):
        raise ValidationError("invalid_wire_type", "expected array", f"{path}.charges")
    return ModelChargeArtifact(
        model_id=data["model_id"],
        model_hash=data["model_hash"],
        atom_order=tuple(data["atom_order"]),
        charges=tuple(data["charges"]),
        atomic_charge_role=data["atomic_charge_role"],
        provider=data["provider"],
        provider_version=data["provider_version"],
        artifact_hash=data["artifact_hash"],
    )


def _constraint_rows(
    topology: PreparedChemicalTopology,
    contract: ChargeAssignmentContract,
    fit_atom_ids: tuple[str, ...],
) -> tuple[list[list[float]], list[float], list[str]]:
    fit_index = {atom_id: index for index, atom_id in enumerate(fit_atom_ids)}
    component_by_id = {component.external_id: component for component in topology.components}
    rows: list[list[float]] = []
    targets: list[float] = []
    labels: list[str] = []

    def append(atom_ids, target: float, label: str) -> None:
        row = [0.0] * len(fit_atom_ids)
        for atom_id in atom_ids:
            if atom_id not in fit_index:
                raise ValidationError("charge_constraint_outside_fit_scope", atom_id, label)
            row[fit_index[atom_id]] = 1.0
        if not any(row):
            raise ValidationError("empty_charge_constraint", label)
        rows.append(row)
        targets.append(float(target))
        labels.append(label)

    if contract.policy is ChargePolicy.PRESERVE_COMPONENT_CHARGES:
        for component_id, target in sorted(contract.component_targets.items()):
            component = component_by_id[component_id]
            if not set(component.atom_ids) <= set(fit_atom_ids):
                raise ValidationError(
                    "partial_component_charge_fit",
                    "preserve-component policy requires the complete component in fit scope",
                    component_id,
                )
            append(component.atom_ids, target, f"component:{component_id}")
    elif contract.policy is ChargePolicy.SITE_JOINT_FIT:
        append(fit_atom_ids, float(contract.site_target), "site")
    else:
        for group in contract.constraint_groups:
            append(group.atom_ids, group.target_charge, f"group:{group.external_id}")
    return rows, targets, labels


def _constrained_projection(
    raw_charges: Mapping[str, float],
    rows: list[list[float]],
    targets: list[float],
) -> dict[str, float]:
    atom_ids = tuple(raw_charges)
    initial = np.asarray([raw_charges[atom_id] for atom_id in atom_ids], dtype=float)
    matrix = np.asarray(rows, dtype=float)
    target = np.asarray(targets, dtype=float)
    if matrix.shape[1] != len(atom_ids):
        raise ValidationError("invalid_charge_constraint_matrix", repr(matrix.shape))
    gram = matrix @ matrix.T
    residual = target - matrix @ initial
    try:
        correction = matrix.T @ np.linalg.solve(gram, residual)
    except np.linalg.LinAlgError as exc:
        raise ValidationError("unsatisfiable_charge_constraints", "constraint matrix is singular") from exc
    fitted = initial + correction
    if not np.isfinite(fitted).all() or not np.allclose(matrix @ fitted, target, atol=1.0e-8, rtol=0.0):
        raise ValidationError("unsatisfiable_charge_constraints", "charge targets were not satisfied")
    return {atom_id: float(fitted[index]) for index, atom_id in enumerate(atom_ids)}


def project_model_charges(
    topology: PreparedChemicalTopology,
    model: DerivedModel,
    artifact: ModelChargeArtifact,
    contract: ChargeAssignmentContract,
) -> tuple[ChargeOverlay, Mapping[str, Any]]:
    """Remove caps, project their charge to declared parents, and enforce targets."""

    return project_model_charge_artifacts(topology, ((model, artifact),), contract)


def project_model_charge_artifact(
    topology: PreparedChemicalTopology,
    model_artifacts: Sequence[tuple[DerivedModel, ModelChargeArtifact]],
    contract: ChargeAssignmentContract,
) -> tuple[PartialChargeArtifact, Mapping[str, Any]]:
    """Project one or more disjoint large-model charge artifacts as one fit scope."""

    validate_charge_contract(contract, topology)
    if not model_artifacts:
        raise ValidationError("missing_charge_artifact", "at least one large-model charge artifact is required")
    raw_parent_charges: dict[str, float] = {}
    cap_projection: dict[str, float] = {}
    discarded_cap_charge = 0.0
    discarded_environment_charge = 0.0
    max_model_constraint_residual = 0.0
    constrained_models = all(
        model.capped_model_manifest is not None
        and isinstance(model.charge_accounting.get("constraint_groups"), list)
        for model, _ in model_artifacts
    )
    model_reports: list[dict[str, Any]] = []
    for model, artifact in model_artifacts:
        validate_model_charge_artifact(model, artifact)
        charge_by_model_atom = dict(zip(artifact.atom_order, artifact.charges))
        if constrained_models:
            for group in model.charge_accounting["constraint_groups"]:
                fitted_group_charge = sum(
                    charge_by_model_atom[atom_id]
                    for atom_id in group["model_atom_ids"]
                )
                residual = fitted_group_charge - float(group["target_charge"])
                max_model_constraint_residual = max(
                    max_model_constraint_residual,
                    abs(residual),
                )
                if abs(residual) > 1.0e-8:
                    raise ValidationError(
                        "model_charge_constraint_not_satisfied",
                        f"{group['group_id']} residual={residual!r}",
                        model.external_id,
                    )
        for atom in model.atoms:
            charge = float(charge_by_model_atom[atom.model_atom_id])
            if atom.role == "cap":
                if constrained_models:
                    discarded_cap_charge += charge
                    continue
                if not atom.cap_parent_external_id or not atom.charge_projection_group:
                    raise ValidationError("missing_cap_charge_projection", atom.model_atom_id)
                if atom.cap_parent_external_id not in contract.fit_atom_ids:
                    raise ValidationError("cap_parent_outside_fit_scope", atom.cap_parent_external_id)
                cap_projection[atom.cap_parent_external_id] = cap_projection.get(atom.cap_parent_external_id, 0.0) + charge
                continue
            if atom.role == "environment":
                if not constrained_models:
                    raise ValidationError(
                        "legacy_environment_charge_projection_unsupported",
                        atom.model_atom_id,
                    )
                discarded_environment_charge += charge
                continue
            if atom.external_id is None:
                raise ValidationError("invalid_charge_parent_mapping", atom.model_atom_id)
            if atom.external_id in raw_parent_charges:
                raise ValidationError("overlapping_charge_model_scope", atom.external_id)
            raw_parent_charges[atom.external_id] = charge
        model_reports.append({
            "model_id": model.external_id,
            "model_hash": model.model_hash,
            "provider": artifact.provider,
            "provider_version": artifact.provider_version,
            "atomic_charge_role": artifact.atomic_charge_role,
            "cap_atom_count": sum(atom.role == "cap" for atom in model.atoms),
            "environment_atom_count": sum(atom.role == "environment" for atom in model.atoms),
        })
    missing_fit_atoms = set(contract.fit_atom_ids) - set(raw_parent_charges)
    if missing_fit_atoms:
        raise ValidationError("incomplete_charge_fit_mapping", ",".join(sorted(missing_fit_atoms)))
    scoped_charges = {
        atom_id: raw_parent_charges[atom_id] + cap_projection.get(atom_id, 0.0)
        for atom_id in contract.fit_atom_ids
    }
    rows, targets, labels = _constraint_rows(topology, contract, contract.fit_atom_ids)
    if constrained_models:
        fitted = dict(scoped_charges)
        charge_vector = np.asarray(
            [fitted[atom_id] for atom_id in contract.fit_atom_ids],
            dtype=float,
        )
        if not np.allclose(
            np.asarray(rows, dtype=float) @ charge_vector,
            np.asarray(targets, dtype=float),
            atol=1.0e-8,
            rtol=0.0,
        ):
            raise ValidationError(
                "projected_charge_contract_not_satisfied",
                "constrained RESP core charges do not satisfy the parent contract",
            )
    else:
        fitted = _constrained_projection(scoped_charges, rows, targets)
    stable_order = {atom.external_id: atom.stable_order for atom in topology.atoms}
    atom_ids = tuple(sorted(fitted, key=lambda atom_id: (stable_order[atom_id], atom_id)))
    provider = "+".join(sorted({artifact.provider for _, artifact in model_artifacts}))
    provider_version = "+".join(sorted({artifact.provider_version for _, artifact in model_artifacts}))
    artifact_identity = _sha256({
        "topology_hash": topology.topology_hash,
        "model_artifact_hashes": sorted(artifact.artifact_hash for _, artifact in model_artifacts),
        "contract": _canonicalize(contract),
    })
    projected_artifact = PartialChargeArtifact(
        schema_version=topology.schema_version,
        artifact_id=f"projected-model-charge:{artifact_identity[:16]}",
        topology_hash=topology.topology_hash,
        graph_revision=topology.graph_revision,
        input_hash=topology.input_hash,
        atom_ids=atom_ids,
        charges=tuple(fitted[atom_id] for atom_id in atom_ids),
        scope_kind="fit_scope",
        scope_id=contract.source,
        atomic_charge_role="fitted",
        precedence=contract.fit_precedence,
        provider=provider,
        provider_version=provider_version,
        method=(
            f"charge-fit-constrained:{contract.policy.value}"
            if constrained_models
            else f"charge-fit-projection:{contract.policy.value}"
        ),
        charge_unit=contract.charge_unit,
        source=contract.source,
    ).with_computed_hash()
    report = {
        "models": tuple(model_reports),
        "policy": contract.policy.value,
        "constraint_labels": tuple(labels),
        "raw_fit_scope_charge": sum(scoped_charges.values()),
        "final_fit_scope_charge": sum(fitted.values()),
        "cap_charge_projected": sum(cap_projection.values()),
        "cap_charge_discarded": discarded_cap_charge,
        "environment_charge_discarded": discarded_environment_charge,
        "constrained_model_charges": constrained_models,
        "max_model_constraint_residual": max_model_constraint_residual,
        "cap_atom_count": sum(item["cap_atom_count"] for item in model_reports),
        "environment_atom_count": sum(
            item["environment_atom_count"] for item in model_reports
        ),
        "max_charge_correction": max(
            abs(fitted[atom_id] - scoped_charges[atom_id]) for atom_id in contract.fit_atom_ids
        ),
        "projected_artifact_hash": projected_artifact.artifact_hash,
        "projected_artifact_id": projected_artifact.artifact_id,
        "precedence": projected_artifact.precedence,
    }
    return projected_artifact, report


def project_model_charge_artifacts(
    topology: PreparedChemicalTopology,
    model_artifacts: Sequence[tuple[DerivedModel, ModelChargeArtifact]],
    contract: ChargeAssignmentContract,
) -> tuple[ChargeOverlay, Mapping[str, Any]]:
    artifact, report = project_model_charge_artifact(topology, model_artifacts, contract)
    charges = dict(zip(artifact.atom_ids, artifact.charges))
    atom_sources = {
        atom_id: {
            "artifact_id": artifact.artifact_id,
            "artifact_hash": artifact.artifact_hash,
            "precedence": artifact.precedence,
            "scope_kind": artifact.scope_kind,
            "scope_id": artifact.scope_id,
            "atomic_charge_role": artifact.atomic_charge_role,
            "provider": artifact.provider,
            "provider_version": artifact.provider_version,
            "method": artifact.method,
            "source": artifact.source,
        }
        for atom_id in artifact.atom_ids
    }
    overlay = ChargeOverlay(
        topology_hash=topology.topology_hash,
        charges=charges,
        source="partial-charge-artifact-composition:v1",
        atom_sources=atom_sources,
        artifact_hashes=(artifact.artifact_hash,),
    ).with_computed_hash()
    return overlay, report


__all__ = [
    "ModelChargeArtifact", "model_charge_artifact_from_dict", "project_model_charge_artifact",
    "project_model_charge_artifacts", "project_model_charges", "validate_model_charge_artifact",
]
