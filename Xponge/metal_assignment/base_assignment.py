"""Source-neutral base-force-field assignment over Chemcore-locked graphs."""

from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
import subprocess
from typing import Any, Mapping

from .contracts import BaseForceFieldOverlay, ValidationError, _canonicalize, _freeze_field, _sha256
from .input import MetalAssignmentPackage, validate_package
from ._worker_runtime import run_worker_subprocess, worker_command
from .base_charge import BaseChargeInput, base_charge_input_to_dict, validate_base_charge_input


WORKER_PROTOCOL_VERSION = 2
SUPPORTED_NATIVE_BASE_PROVIDERS = frozenset({"gaff", "gaff2"})


@dataclass(frozen=True, slots=True)
class BaseAssignmentReport:
    """Hash-closed audit for one provider-specific base overlay fragment."""

    schema_version: int
    topology_hash: str
    projection_hash: str
    provider: str
    provider_version: str
    component_ids: tuple[str, ...]
    proof_hashes: Mapping[str, str]
    artifact_hashes: Mapping[str, str]
    component_output_hashes: Mapping[str, str]
    frcmod_hashes: Mapping[str, str]
    ignored_force_kinds: Mapping[str, tuple[str, ...]]
    atom_count: int
    bonded_term_count: int
    parameter_source: str
    charge_method: str = ""
    charge_source: str = ""
    charge_input_hash: str = ""
    report_hash: str = ""

    def __post_init__(self) -> None:
        for name in (
            "proof_hashes", "artifact_hashes", "component_output_hashes", "frcmod_hashes",
            "ignored_force_kinds",
        ):
            _freeze_field(self, name)

    def canonical_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "topology_hash": self.topology_hash,
            "projection_hash": self.projection_hash,
            "provider": self.provider,
            "provider_version": self.provider_version,
            "charge_method": self.charge_method,
            "charge_source": self.charge_source,
            "charge_input_hash": self.charge_input_hash,
            "component_ids": _canonicalize(self.component_ids),
            "proof_hashes": _canonicalize(self.proof_hashes),
            "artifact_hashes": _canonicalize(self.artifact_hashes),
            "component_output_hashes": _canonicalize(self.component_output_hashes),
            "frcmod_hashes": _canonicalize(self.frcmod_hashes),
            "ignored_force_kinds": _canonicalize(self.ignored_force_kinds),
            "atom_count": self.atom_count,
            "bonded_term_count": self.bonded_term_count,
            "parameter_source": self.parameter_source,
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "BaseAssignmentReport":
        return replace(self, report_hash=self.computed_hash())


@dataclass(frozen=True, slots=True)
class BaseAssignmentOutput:
    """Provider-specific overlay fragment and its audit report."""

    overlay: BaseForceFieldOverlay
    report: BaseAssignmentReport


def _worker_payload(
    package: MetalAssignmentPackage,
    *,
    provider: str,
    complete_missing_parameters: bool,
    base_charge_input: BaseChargeInput | None,
) -> dict[str, Any]:
    request = package.request
    topology = request.topology
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    artifact_by_id = {
        artifact.external_id: artifact
        for artifact in package.prepared_artifacts.structural_artifacts.assignment_components
    }
    components = tuple(
        component
        for component in request.assignment_components
        if component.base_force_field and component.provider == provider
    )
    if not components:
        raise ValidationError(
            "missing_provider_scope",
            f"no base assignment component uses provider {provider}",
            "assignment_components",
        )

    payload_components: list[dict[str, Any]] = []
    active_links = {link.external_id: link for link in topology.links if link.active}
    active_bonds = {bond.external_id: bond for bond in topology.bonds if bond.active}
    for component in components:
        path = f"assignment_components.{component.external_id}"
        component_atom_ids = set(component.atom_ids)
        if any(atom_by_id[atom_id].is_metal for atom_id in component.atom_ids):
            raise ValidationError("metal_in_base_force_field", component.external_id, path)
        proof = component.chemical_topology_proof
        if not proof.valence_complete or proof.explicit_hydrogen_status != "complete":
            raise ValidationError("incomplete_chemical_topology_proof", component.external_id, path)
        internal_links = [
            link for link in active_links.values() if set(link.atom_ids) <= component_atom_ids
        ]
        if internal_links:
            raise ValidationError(
                "base_assignment_link_requires_bond_order",
                ",".join(link.external_id for link in internal_links),
                path,
            )
        bonds = [
            bond for bond in active_bonds.values() if set(bond.atom_ids) <= component_atom_ids
        ]
        expected_edge_ids = set(artifact_by_id[component.external_id].active_edge_ids)
        if {bond.external_id for bond in bonds} != expected_edge_ids:
            raise ValidationError("base_assignment_edge_mismatch", component.external_id, path)
        for bond in bonds:
            if (
                isinstance(bond.order, bool)
                or not isinstance(bond.order, (int, float))
                or not math.isfinite(bond.order)
                or bond.order <= 0
                or bond.order > 3
            ):
                raise ValidationError("invalid_base_assignment_bond_order", bond.external_id, path)
        artifact = artifact_by_id[component.external_id]
        payload_components.append({
            "external_id": component.external_id,
            "net_formal_charge": component.net_formal_charge,
            "proof_hash": proof.proof_hash,
            "artifact_hash": artifact.artifact_hash,
            "atoms": [
                {
                    "external_id": atom_id,
                    "element": atom_by_id[atom_id].element,
                    "coordinates": list(atom_by_id[atom_id].coordinates),
                    "name": f"A{index + 1}",
                    "formal_charge": (
                        atom_by_id[atom_id].formal_charge
                        if atom_by_id[atom_id].formal_charge_known
                        else None
                    ),
                }
                for index, atom_id in enumerate(component.atom_ids)
            ],
            "bonds": [
                {
                    "external_id": bond.external_id,
                    "atom_ids": list(bond.atom_ids),
                    "order": bond.order,
                    "semantic": bond.semantic,
                    "source": bond.source,
                }
                for bond in bonds
            ],
        })
    return {
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "provider": provider,
        "topology_hash": topology.topology_hash,
        "complete_missing_parameters": complete_missing_parameters,
        "base_charge_input": (
            base_charge_input_to_dict(base_charge_input)
            if base_charge_input is not None
            else None
        ),
        "components": payload_components,
    }


def _invoke_native_worker(payload: Mapping[str, Any], *, timeout_seconds: float) -> dict[str, Any]:
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValidationError("invalid_worker_timeout", str(timeout_seconds))
    try:
        completed = run_worker_subprocess(
            worker_command("Xponge.metal_assignment._gaff_worker"),
            input_text=json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False),
            timeout_seconds=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        raise ValidationError("base_assignment_worker_timeout", str(timeout_seconds)) from exc
    try:
        response = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        stderr = completed.stderr.strip()
        raise ValidationError(
            "invalid_base_assignment_worker_output",
            stderr[-1000:] if stderr else str(exc),
        ) from exc
    if not isinstance(response, dict) or response.get("ok") is not True:
        error = response.get("error", {}) if isinstance(response, dict) else {}
        message = error.get("message", "worker failed") if isinstance(error, dict) else "worker failed"
        stderr = completed.stderr.strip()
        if stderr:
            message = f"{message}; stderr={stderr[-1000:]}"
        raise ValidationError("base_assignment_worker_failed", message)
    if completed.returncode != 0:
        raise ValidationError("base_assignment_worker_failed", f"exit status {completed.returncode}")
    return response


def _validate_worker_response(payload: Mapping[str, Any], response: Mapping[str, Any]) -> None:
    required = {
        "ok", "protocol_version", "provider", "topology_hash", "provider_version", "components",
        "response_hash",
    }
    if set(response) != required:
        raise ValidationError("invalid_base_assignment_worker_response", "unexpected fields")
    if (
        response["protocol_version"] != WORKER_PROTOCOL_VERSION
        or response["provider"] != payload["provider"]
        or response["topology_hash"] != payload["topology_hash"]
        or not isinstance(response["provider_version"], str)
        or not response["provider_version"]
        or not isinstance(response["components"], list)
    ):
        raise ValidationError("base_assignment_worker_identity_mismatch", "worker response belongs to another request")
    response_payload = dict(response)
    response_hash = response_payload.pop("response_hash")
    if response_hash != _sha256(response_payload):
        raise ValidationError("stale_base_assignment_worker_hash", "response hash mismatch")
    expected_by_id = {component["external_id"]: component for component in payload["components"]}
    actual_ids = [component.get("external_id") for component in response["components"] if isinstance(component, dict)]
    if len(actual_ids) != len(response["components"]) or set(actual_ids) != set(expected_by_id) or len(actual_ids) != len(set(actual_ids)):
        raise ValidationError("base_assignment_worker_component_mismatch", "worker component coverage mismatch")
    for component in response["components"]:
        required_component = {
            "external_id", "proof_hash", "artifact_hash", "net_formal_charge", "atom_types",
            "charges", "masses", "lj_parameters", "bonded_parameters", "frcmod_sha256",
            "ignored_force_kinds", "output_hash",
        }
        if set(component) != required_component:
            raise ValidationError("invalid_base_assignment_worker_component", str(component.get("external_id")))
        expected = expected_by_id[component["external_id"]]
        if (
            component["proof_hash"] != expected["proof_hash"]
            or component["artifact_hash"] != expected["artifact_hash"]
            or component["net_formal_charge"] != expected["net_formal_charge"]
        ):
            raise ValidationError("base_assignment_worker_component_identity_mismatch", component["external_id"])
        component_payload = dict(component)
        output_hash = component_payload.pop("output_hash")
        if output_hash != _sha256(component_payload):
            raise ValidationError("stale_base_assignment_worker_component_hash", component["external_id"])
        expected_atom_ids = {atom["external_id"] for atom in expected["atoms"]}
        expected_charge_ids = expected_atom_ids if payload["base_charge_input"] is not None else set()
        if (
            set(component["atom_types"]) != expected_atom_ids
            or set(component["charges"]) != expected_charge_ids
            or set(component["masses"]) != expected_atom_ids
            or set(component["lj_parameters"]) != expected_atom_ids
        ):
            raise ValidationError("base_assignment_worker_atom_coverage_mismatch", component["external_id"])
        if component["charges"]:
            charges = tuple(component["charges"].values())
            if any(
                isinstance(charge, bool)
                or not isinstance(charge, (int, float))
                or not math.isfinite(charge)
                for charge in charges
            ):
                raise ValidationError("invalid_base_assignment_charge", component["external_id"])
            if not math.isclose(
                sum(float(charge) for charge in charges),
                float(expected["net_formal_charge"]),
                abs_tol=1e-6,
            ):
                raise ValidationError("base_assignment_charge_not_conserved", component["external_id"])


def assign_base_force_field(
    package: MetalAssignmentPackage,
    provider: str,
    *,
    complete_missing_parameters: bool = True,
    timeout_seconds: float = 120.0,
    base_charge_input: BaseChargeInput | None = None,
) -> BaseAssignmentOutput:
    """Assign one provider scope without mutating or reconstructing topology.

    The returned overlay is a provider-specific fragment.  A complete
    ``ParameterizationResult`` may merge it with other base providers (for
    example a standard biomolecular provider) before result validation.
    """

    validate_package(package)
    normalized_provider = provider.lower()
    if normalized_provider not in SUPPORTED_NATIVE_BASE_PROVIDERS:
        raise ValidationError("unsupported_base_force_field_provider", provider)
    if normalized_provider not in package.request.capability_snapshot.base_force_field_providers:
        raise ValidationError("provider_not_in_capability_snapshot", normalized_provider)
    if base_charge_input is not None:
        validate_base_charge_input(base_charge_input)
    payload = _worker_payload(
        package,
        provider=normalized_provider,
        complete_missing_parameters=complete_missing_parameters,
        base_charge_input=base_charge_input,
    )
    response = _invoke_native_worker(payload, timeout_seconds=timeout_seconds)
    _validate_worker_response(payload, response)

    atom_types: dict[str, str] = {}
    charges: dict[str, float] = {}
    masses: dict[str, float] = {}
    lj_parameters: dict[str, Any] = {}
    bonded_parameters: dict[str, Any] = {}
    for component in response["components"]:
        for target, source, code in (
            (atom_types, component["atom_types"], "duplicate_base_atom_type"),
            (charges, component["charges"], "duplicate_base_charge"),
            (masses, component["masses"], "duplicate_base_mass"),
            (lj_parameters, component["lj_parameters"], "duplicate_base_lj_parameter"),
            (bonded_parameters, component["bonded_parameters"], "duplicate_base_bonded_term"),
        ):
            overlap = set(target) & set(source)
            if overlap:
                raise ValidationError(code, ",".join(sorted(overlap)))
            target.update(source)
    covered_atom_ids = tuple(
        atom["external_id"]
        for component in payload["components"]
        for atom in component["atoms"]
    )
    parameter_source = f"xponge:{normalized_provider}:{response['provider_version']}"
    overlay = BaseForceFieldOverlay(
        topology_hash=package.request.topology.topology_hash,
        covered_atom_ids=covered_atom_ids,
        atom_types=atom_types,
        charges=charges,
        masses=masses,
        lj_parameters=lj_parameters,
        bonded_parameters=bonded_parameters,
        parameter_source=parameter_source,
    )
    components = response["components"]
    report = BaseAssignmentReport(
        schema_version=package.request.schema_version,
        topology_hash=package.request.topology.topology_hash,
        projection_hash=package.request.projection_hash,
        provider=normalized_provider,
        provider_version=response["provider_version"],
        charge_method=base_charge_input.method if base_charge_input is not None else "",
        charge_source=base_charge_input.source if base_charge_input is not None else "",
        charge_input_hash=(
            base_charge_input.charge_input_hash if base_charge_input is not None else ""
        ),
        component_ids=tuple(component["external_id"] for component in components),
        proof_hashes={component["external_id"]: component["proof_hash"] for component in components},
        artifact_hashes={component["external_id"]: component["artifact_hash"] for component in components},
        component_output_hashes={component["external_id"]: component["output_hash"] for component in components},
        frcmod_hashes={component["external_id"]: component["frcmod_sha256"] for component in components},
        ignored_force_kinds={
            component["external_id"]: tuple(component["ignored_force_kinds"])
            for component in components
        },
        atom_count=len(covered_atom_ids),
        bonded_term_count=len(bonded_parameters),
        parameter_source=parameter_source,
    ).with_computed_hash()
    return BaseAssignmentOutput(overlay, report)


__all__ = [
    "BaseAssignmentOutput", "BaseAssignmentReport", "SUPPORTED_NATIVE_BASE_PROVIDERS",
    "assign_base_force_field",
]
