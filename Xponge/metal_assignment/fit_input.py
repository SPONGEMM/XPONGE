"""Hash-closed file contract for bonded fit artifacts."""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import json
from pathlib import Path
from typing import Any, Mapping

from .bonded_fit import BondedMetalParameterSpec, MetalAtomParameterSpec
from .charge_fit import ModelChargeArtifact, model_charge_artifact_from_dict
from .contracts import ValidationError, _canonicalize, _sha256, _strict_object
from .force_fit import HessianArtifact, hessian_artifact_from_dict


@dataclass(frozen=True, slots=True)
class BondedFitInput:
    schema_version: int
    metal_parameter_spec: BondedMetalParameterSpec
    charge_artifacts: tuple[ModelChargeArtifact, ...]
    hessian_artifacts: tuple[HessianArtifact, ...]
    force_method: str
    scale_factor: float = 1.0
    empirical_registry_id: str = ""
    empirical_geometry: str = "unclassified"
    empirical_base_force_field: str = ""
    empirical_water_model: str = ""
    fit_input_hash: str = ""

    def canonical_payload(self) -> dict[str, Any]:
        return {
            item.name: _canonicalize(getattr(self, item.name))
            for item in fields(self)
            if item.name != "fit_input_hash"
        }

    def computed_hash(self) -> str:
        return _sha256(self.canonical_payload())

    def with_computed_hash(self) -> "BondedFitInput":
        return replace(self, fit_input_hash=self.computed_hash())


def validate_bonded_fit_input(value: BondedFitInput) -> None:
    if value.schema_version != 1:
        raise ValidationError("unsupported_schema_version", str(value.schema_version), "bonded_fit_input")
    if value.force_method not in {"seminario", "empirical_registry", "empirical_zn_nos"}:
        raise ValidationError("unsupported_bonded_force_method", value.force_method)
    if value.force_method == "empirical_registry" and not value.empirical_registry_id:
        raise ValidationError("missing_empirical_registry_id", value.force_method)
    if value.force_method == "seminario" and value.empirical_registry_id:
        raise ValidationError("unexpected_empirical_registry_id", value.empirical_registry_id)
    for name in (
        "empirical_registry_id", "empirical_geometry",
        "empirical_base_force_field", "empirical_water_model",
    ):
        if not isinstance(getattr(value, name), str):
            raise ValidationError("invalid_wire_type", "expected string", f"bonded_fit_input.{name}")
    if not value.fit_input_hash or value.fit_input_hash != value.computed_hash():
        raise ValidationError("stale_bonded_fit_input_hash", "fit input hash mismatch")


def default_bonded_fit_input(
    package: Any,
    *,
    parameter_source: str = "qm_fit",
    water_model: str,
    timeout_seconds: float = 60.0,
) -> BondedFitInput:
    """Build the provider-owned metal baseline for a bonded fit.

    Metal mass/LJ identity comes from the same water-model ion registry used by
    the nonbonded path.  RESP/Hessian artifacts remain explicit later inputs.
    """

    from .input import validate_package
    from .metal_overlay import resolve_metal_ion_parameters
    validate_package(package)
    metals, parameters, ion_source, ion_provenance = resolve_metal_ion_parameters(
        package,
        water_model=water_model,
        timeout_seconds=timeout_seconds,
    )
    if parameter_source == "qm_fit":
        force_method = "seminario"
        registry_id = ""
    elif parameter_source.startswith("empirical:") and parameter_source[10:]:
        force_method = "empirical_registry"
        registry_id = parameter_source[10:]
    else:
        raise ValidationError("unsupported_parameter_source", parameter_source)
    metal_spec = BondedMetalParameterSpec(
        schema_version=package.request.schema_version,
        topology_hash=package.request.topology.topology_hash,
        metal_atoms=tuple(
            MetalAtomParameterSpec(
                external_id=atom.external_id,
                element=atom.element,
                formal_charge=int(atom.formal_charge),
                atom_type=str(parameters[atom.external_id]["atom_type"]),
                partial_charge=float(parameters[atom.external_id]["charge"]),
                mass=float(parameters[atom.external_id]["mass"]),
                epsilon=float(parameters[atom.external_id]["lj"]["epsilon"]),
                rmin=float(parameters[atom.external_id]["lj"]["rmin"]),
            )
            for atom in metals
        ),
        donor_atom_types={},
        parameter_source=f"{parameter_source};baseline={ion_source}",
        provenance={
            "parameter_source": parameter_source,
            "baseline_ion_registry": ion_provenance,
        },
    ).with_computed_hash()
    result = BondedFitInput(
        schema_version=package.request.schema_version,
        metal_parameter_spec=metal_spec,
        charge_artifacts=(),
        hessian_artifacts=(),
        force_method=force_method,
        scale_factor=1.0,
        empirical_registry_id=registry_id,
        empirical_geometry="unclassified",
        empirical_base_force_field="",
        empirical_water_model=water_model,
    ).with_computed_hash()
    validate_bonded_fit_input(result)
    return result


def _strings(value: Any, path: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValidationError("invalid_wire_type", "expected string array", path)
    return tuple(value)


def _parse_metal_spec(value: Any) -> BondedMetalParameterSpec:
    data = _strict_object(
        value,
        required={
            "schema_version", "topology_hash", "metal_atoms", "donor_atom_types", "parameter_source",
            "provenance", "precedence", "spec_hash",
        },
        path="bonded_fit_input.metal_parameter_spec",
    )
    if not isinstance(data["metal_atoms"], list):
        raise ValidationError("invalid_wire_type", "expected array", "bonded_fit_input.metal_parameter_spec.metal_atoms")
    metal_atoms = []
    for index, item in enumerate(data["metal_atoms"]):
        path = f"bonded_fit_input.metal_parameter_spec.metal_atoms[{index}]"
        atom = _strict_object(
            item,
            required={
                "external_id", "element", "formal_charge", "atom_type", "partial_charge", "mass",
                "epsilon", "rmin",
            },
            path=path,
        )
        metal_atoms.append(MetalAtomParameterSpec(**atom))
    if not isinstance(data["donor_atom_types"], Mapping) or not isinstance(data["provenance"], Mapping):
        raise ValidationError("invalid_wire_type", "expected object", "bonded_fit_input.metal_parameter_spec")
    return BondedMetalParameterSpec(
        schema_version=data["schema_version"],
        topology_hash=data["topology_hash"],
        metal_atoms=tuple(metal_atoms),
        donor_atom_types=data["donor_atom_types"],
        parameter_source=data["parameter_source"],
        provenance=data["provenance"],
        precedence=data["precedence"],
        spec_hash=data["spec_hash"],
    )


def _parse_charge_artifact(value: Any, index: int) -> ModelChargeArtifact:
    path = f"bonded_fit_input.charge_artifacts[{index}]"
    return model_charge_artifact_from_dict(value, path=path)


def _parse_hessian_artifact(value: Any, index: int) -> HessianArtifact:
    path = f"bonded_fit_input.hessian_artifacts[{index}]"
    return hessian_artifact_from_dict(value, path=path)


def bonded_fit_input_from_dict(value: Any) -> BondedFitInput:
    data = _strict_object(
        value,
        required={
            "schema_version", "metal_parameter_spec", "charge_artifacts", "hessian_artifacts",
            "force_method", "scale_factor", "fit_input_hash",
        },
        optional={
            "empirical_registry_id", "empirical_geometry",
            "empirical_base_force_field", "empirical_water_model",
        },
        path="bonded_fit_input",
    )
    if not isinstance(data["charge_artifacts"], list) or not isinstance(data["hessian_artifacts"], list):
        raise ValidationError("invalid_wire_type", "artifact fields must be arrays", "bonded_fit_input")
    result = BondedFitInput(
        schema_version=data["schema_version"],
        metal_parameter_spec=_parse_metal_spec(data["metal_parameter_spec"]),
        charge_artifacts=tuple(
            _parse_charge_artifact(item, index) for index, item in enumerate(data["charge_artifacts"])
        ),
        hessian_artifacts=tuple(
            _parse_hessian_artifact(item, index) for index, item in enumerate(data["hessian_artifacts"])
        ),
        force_method=data["force_method"],
        scale_factor=data["scale_factor"],
        empirical_registry_id=data.get("empirical_registry_id", ""),
        empirical_geometry=data.get("empirical_geometry", "unclassified"),
        empirical_base_force_field=data.get("empirical_base_force_field", ""),
        empirical_water_model=data.get("empirical_water_model", ""),
        fit_input_hash=data["fit_input_hash"],
    )
    validate_bonded_fit_input(result)
    return result


def bonded_fit_input_to_dict(value: BondedFitInput) -> dict[str, Any]:
    validate_bonded_fit_input(value)
    return _canonicalize(value)


def bonded_fit_input_dumps(value: BondedFitInput) -> str:
    return json.dumps(bonded_fit_input_to_dict(value), sort_keys=True, separators=(",", ":"), allow_nan=False)


def bonded_fit_input_loads(payload: str) -> BondedFitInput:
    try:
        value = json.loads(payload)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValidationError("invalid_json", str(exc), "bonded_fit_input") from exc
    return bonded_fit_input_from_dict(value)


def load_bonded_fit_input(path: str | Path) -> BondedFitInput:
    source_path = Path(path).resolve()
    try:
        return bonded_fit_input_loads(source_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValidationError("fit_input_read_failed", str(exc), str(source_path)) from exc


__all__ = [
    "BondedFitInput", "bonded_fit_input_dumps", "bonded_fit_input_from_dict",
    "bonded_fit_input_loads", "bonded_fit_input_to_dict", "load_bonded_fit_input",
    "validate_bonded_fit_input", "default_bonded_fit_input",
]
