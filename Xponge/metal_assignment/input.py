"""Unified Metal-assignment input package and Chemcore-result adapter."""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
from typing import Any

from .artifacts import (
    DerivedModelBundle,
    PreparedArtifactPackage,
    derived_model_bundle_from_dict,
    prepared_artifact_package_from_dict,
    validate_artifact_package,
    validate_derived_models,
)
from .contracts import (
    ChargeAssignmentContract,
    ElectronicState,
    MetalAssignmentInput,
    PartialChargeArtifactBundle,
    ProviderCapabilitySnapshot,
    SCHEMA_VERSION,
    ValidationError,
    _canonicalize,
    _parse_assignment_component,
    _parse_partition_proof,
    _parse_projection,
    _parse_topology,
    _sha256,
    _strict_object,
    metal_assignment_input_from_dict,
    validate_input,
)


@dataclass(frozen=True, slots=True)
class MetalAssignmentPackage:
    request: MetalAssignmentInput
    prepared_artifacts: PreparedArtifactPackage
    package_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "request": _canonicalize(self.request),
            "prepared_artifacts": _canonicalize(self.prepared_artifacts),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "MetalAssignmentPackage":
        return replace(self, package_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class MetalLocalModelPackage:
    """Bonded-only package containing no prepared or assigned full system."""

    request: MetalAssignmentInput
    derived_models: DerivedModelBundle
    package_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "request": _canonicalize(self.request),
            "derived_models": _canonicalize(self.derived_models),
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "MetalLocalModelPackage":
        return replace(self, package_hash=self.computed_hash())


def validate_package(package: MetalAssignmentPackage | MetalLocalModelPackage) -> None:
    validate_input(package.request)
    if isinstance(package, MetalLocalModelPackage):
        validate_derived_models(package.request, package.derived_models)
    else:
        validate_artifact_package(package.request, package.prepared_artifacts)
    if not package.package_hash:
        raise ValidationError("missing_package_hash", "a computed package hash is required")
    if package.package_hash != package.computed_hash():
        raise ValidationError("stale_package_hash", "package hash does not match canonical payload")


def package_derived_models(
    package: MetalAssignmentPackage | MetalLocalModelPackage,
) -> DerivedModelBundle:
    """Return the local derived models without exposing a prepared-system API."""

    return (
        package.derived_models
        if isinstance(package, MetalLocalModelPackage)
        else package.prepared_artifacts.derived_models
    )


def metal_assignment_package_from_chemcore_result(
    value: Any,
    *,
    request_id: str,
    electronic_state: ElectronicState,
    capability_snapshot: ProviderCapabilitySnapshot,
    charge_contract: ChargeAssignmentContract | None = None,
    partial_charge_artifacts: PartialChargeArtifactBundle | None = None,
) -> MetalAssignmentPackage:
    """Adapt a source-neutral Chemcore Phase-1 result without changing topology."""

    data = _strict_object(
        value,
        required={"topology", "structural_artifact_bundle", "derived_model_bundle", "atom_mapping"},
        path="chemcore_preparation_result",
    )
    topology_result = _strict_object(
        data["topology"],
        required={
            "schema_version", "graph_revision", "input_hash", "interaction_model", "prepared_topology",
            "projections", "assignment_components", "partition_proofs", "warnings",
        },
        path="chemcore_preparation_result.topology",
    )
    for name in ("projections", "assignment_components", "partition_proofs", "warnings"):
        if not isinstance(topology_result[name], list):
            raise ValidationError("invalid_wire_type", "expected array", f"chemcore_preparation_result.topology.{name}")
    if topology_result["schema_version"] != SCHEMA_VERSION:
        raise ValidationError(
            "unsupported_schema_version",
            str(topology_result["schema_version"]),
            "chemcore_preparation_result.topology.schema_version",
        )
    topology = _parse_topology(topology_result["prepared_topology"])
    request = MetalAssignmentInput(
        schema_version=topology_result["schema_version"],
        request_id=request_id,
        interaction_model=topology_result["interaction_model"],
        electronic_state=electronic_state,
        graph_revision=topology_result["graph_revision"],
        input_hash=topology_result["input_hash"],
        topology=topology,
        projections=tuple(
            _parse_projection(item, f"projections[{index}]")
            for index, item in enumerate(topology_result["projections"])
        ),
        assignment_components=tuple(
            _parse_assignment_component(item, f"assignment_components[{index}]")
            for index, item in enumerate(topology_result["assignment_components"])
        ),
        partition_proofs=tuple(
            _parse_partition_proof(item, f"partition_proofs[{index}]")
            for index, item in enumerate(topology_result["partition_proofs"])
        ),
        capability_snapshot=capability_snapshot,
        charge_contract=charge_contract,
        partial_charge_artifacts=partial_charge_artifacts,
    ).with_computed_hash()
    prepared_artifacts = prepared_artifact_package_from_dict({
        "structural_artifacts": data["structural_artifact_bundle"],
        "atom_mapping": data["atom_mapping"],
        "derived_models": data["derived_model_bundle"],
    })
    package = MetalAssignmentPackage(request, prepared_artifacts).with_computed_hash()
    validate_package(package)
    return package


def metal_local_model_package_from_chemcore_result(
    value: Any,
    *,
    request_id: str,
    electronic_state: ElectronicState,
    capability_snapshot: ProviderCapabilitySnapshot,
    charge_contract: ChargeAssignmentContract | None = None,
    partial_charge_artifacts: PartialChargeArtifactBundle | None = None,
) -> MetalLocalModelPackage:
    """Adapt a Chemcore focus environment into the bonded local-model contract."""

    if not isinstance(value, dict):
        raise ValidationError(
            "invalid_wire_type",
            "expected object",
            "chemcore_preparation_result",
        )
    required = {"topology", "derived_model_bundle"}
    missing = required - set(value)
    if missing:
        raise ValidationError(
            "missing_wire_fields",
            ",".join(sorted(missing)),
            "chemcore_preparation_result",
        )
    topology_result = _strict_object(
        value["topology"],
        required={
            "schema_version", "graph_revision", "input_hash", "interaction_model", "prepared_topology",
            "projections", "assignment_components", "partition_proofs", "warnings",
        },
        path="chemcore_preparation_result.topology",
    )
    for name in ("projections", "assignment_components", "partition_proofs", "warnings"):
        if not isinstance(topology_result[name], list):
            raise ValidationError(
                "invalid_wire_type",
                "expected array",
                f"chemcore_preparation_result.topology.{name}",
            )
    if topology_result["schema_version"] != SCHEMA_VERSION:
        raise ValidationError(
            "unsupported_schema_version",
            str(topology_result["schema_version"]),
            "chemcore_preparation_result.topology.schema_version",
        )
    topology = _parse_topology(topology_result["prepared_topology"])
    request = MetalAssignmentInput(
        schema_version=topology_result["schema_version"],
        request_id=request_id,
        interaction_model=topology_result["interaction_model"],
        electronic_state=electronic_state,
        graph_revision=topology_result["graph_revision"],
        input_hash=topology_result["input_hash"],
        topology=topology,
        projections=tuple(
            _parse_projection(item, f"projections[{index}]")
            for index, item in enumerate(topology_result["projections"])
        ),
        assignment_components=tuple(
            _parse_assignment_component(item, f"assignment_components[{index}]")
            for index, item in enumerate(topology_result["assignment_components"])
        ),
        partition_proofs=tuple(
            _parse_partition_proof(item, f"partition_proofs[{index}]")
            for index, item in enumerate(topology_result["partition_proofs"])
        ),
        capability_snapshot=capability_snapshot,
        charge_contract=charge_contract,
        partial_charge_artifacts=partial_charge_artifacts,
    ).with_computed_hash()
    package = MetalLocalModelPackage(
        request,
        derived_model_bundle_from_dict(value["derived_model_bundle"]),
    ).with_computed_hash()
    validate_package(package)
    return package


def metal_assignment_package_to_dict(package: MetalAssignmentPackage) -> dict[str, Any]:
    validate_package(package)
    return _canonicalize(package)


def metal_assignment_package_from_dict(value: Any) -> MetalAssignmentPackage:
    data = _strict_object(
        value,
        required={"request", "prepared_artifacts", "package_hash"},
        path="package",
    )
    package = MetalAssignmentPackage(
        metal_assignment_input_from_dict(data["request"]),
        prepared_artifact_package_from_dict(data["prepared_artifacts"]),
        data["package_hash"],
    )
    validate_package(package)
    return package


def metal_assignment_package_dumps(package: MetalAssignmentPackage) -> str:
    return json.dumps(
        metal_assignment_package_to_dict(package),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def metal_assignment_package_loads(payload: str) -> MetalAssignmentPackage:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc)) from exc
    return metal_assignment_package_from_dict(value)


def metal_local_model_package_to_dict(package: MetalLocalModelPackage) -> dict[str, Any]:
    validate_package(package)
    return _canonicalize(package)


def metal_local_model_package_from_dict(value: Any) -> MetalLocalModelPackage:
    data = _strict_object(
        value,
        required={"request", "derived_models", "package_hash"},
        path="local_model_package",
    )
    package = MetalLocalModelPackage(
        metal_assignment_input_from_dict(data["request"]),
        derived_model_bundle_from_dict(data["derived_models"]),
        data["package_hash"],
    )
    validate_package(package)
    return package


__all__ = [
    "MetalAssignmentPackage", "MetalLocalModelPackage",
    "metal_assignment_package_dumps", "metal_assignment_package_from_chemcore_result",
    "metal_assignment_package_from_dict", "metal_assignment_package_loads", "metal_assignment_package_to_dict",
    "metal_local_model_package_from_chemcore_result", "metal_local_model_package_from_dict",
    "metal_local_model_package_to_dict", "package_derived_models", "validate_package",
]
