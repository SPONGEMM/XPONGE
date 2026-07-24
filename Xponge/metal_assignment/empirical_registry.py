"""Versioned empirical metal-parameter registry with strict applicability.

Registry records are owned by Xponge and contain no Mokda or Chemcore types.
They describe parameter applicability and provenance; callers must request a
specific registry ID and no resolver path silently falls back to empirical
parameters.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations_with_replacement
from typing import Any, Mapping

from .contracts import (
    PreparedChemicalTopology,
    ValidationError,
    _canonicalize,
    _freeze_field,
    _sha256,
)


@dataclass(frozen=True, slots=True)
class EmpiricalCitation:
    title: str
    authors: tuple[str, ...]
    journal: str
    year: int
    doi: str
    relation: str = "method_context"


@dataclass(frozen=True, slots=True)
class EmpiricalApplicability:
    interaction_model: str
    element: str
    oxidation_state: int
    coordination_numbers: tuple[int, ...]
    donor_elements: tuple[str, ...]
    donor_patterns: tuple[tuple[str, ...], ...]
    geometries: tuple[str, ...]
    compatible_base_force_fields: tuple[str, ...]
    compatible_water_models: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class EmpiricalRegistryEntry:
    registry_id: str
    version: str
    status: str
    applicability: EmpiricalApplicability
    parameters: Mapping[str, Any]
    units: Mapping[str, str]
    source: str
    license: str
    redistribution_status: str
    citations: tuple[EmpiricalCitation, ...]
    validation_cases: tuple[str, ...]
    source_hash: str = ""
    parameter_hash: str = ""

    def __post_init__(self) -> None:
        _freeze_field(self, "parameters")
        _freeze_field(self, "units")

    def source_payload(self) -> dict[str, Any]:
        return {
            "registry_id": self.registry_id,
            "version": self.version,
            "status": self.status,
            "source": self.source,
            "license": self.license,
            "redistribution_status": self.redistribution_status,
            "citations": _canonicalize(self.citations),
            "validation_cases": _canonicalize(self.validation_cases),
        }

    def parameter_payload(self) -> dict[str, Any]:
        return {
            "registry_id": self.registry_id,
            "version": self.version,
            "applicability": _canonicalize(self.applicability),
            "parameters": _canonicalize(self.parameters),
            "units": _canonicalize(self.units),
        }

    def computed_source_hash(self) -> str:
        return _sha256(self.source_payload())

    def computed_parameter_hash(self) -> str:
        return _sha256(self.parameter_payload())

    def with_computed_hashes(self) -> "EmpiricalRegistryEntry":
        return replace(
            self,
            source_hash=self.computed_source_hash(),
            parameter_hash=self.computed_parameter_hash(),
        )


@dataclass(frozen=True, slots=True)
class EmpiricalCenterMatch:
    metal_atom_id: str
    element: str
    oxidation_state: int
    coordination_number: int
    donor_elements: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class EmpiricalRegistryMatch:
    entry: EmpiricalRegistryEntry
    centers: tuple[EmpiricalCenterMatch, ...]
    geometry: str
    base_force_field: str
    water_model: str

    def descriptor(self) -> dict[str, Any]:
        return {
            "registry_id": self.entry.registry_id,
            "version": self.entry.version,
            "status": self.entry.status,
            "centers": _canonicalize(self.centers),
            "geometry": self.geometry,
            "base_force_field": self.base_force_field,
            "water_model": self.water_model,
            "source": self.entry.source,
            "source_hash": self.entry.source_hash,
            "parameter_hash": self.entry.parameter_hash,
            "license": self.entry.license,
            "redistribution_status": self.entry.redistribution_status,
            "citations": _canonicalize(self.entry.citations),
            "validation_cases": _canonicalize(self.entry.validation_cases),
        }


_ZN_EMPIRICAL_BOND_TABLES = {
    "N": (
        (1.926, 124.58), (1.947, 113.59), (1.978, 93.97), (1.982, 90.76), (1.984, 90.80),
        (2.011, 78.18), (2.027, 70.85), (2.028, 70.86), (2.029, 66.61), (2.040, 66.69),
        (2.041, 66.11), (2.047, 62.61), (2.073, 52.77), (2.089, 50.10), (2.133, 36.30),
        (2.145, 35.80), (2.176, 29.24),
    ),
    "O": (
        (1.860, 169.29), (1.986, 76.81), (2.011, 71.26), (2.054, 56.37), (2.109, 41.86),
        (2.112, 41.32),
    ),
    "S": (
        (2.262, 88.50), (2.305, 67.39), (2.353, 50.99), (2.355, 51.79), (2.426, 32.69),
    ),
}


def _donor_patterns(
    donor_elements: tuple[str, ...],
    coordination_numbers: tuple[int, ...],
) -> tuple[tuple[str, ...], ...]:
    return tuple(
        pattern
        for coordination_number in coordination_numbers
        for pattern in combinations_with_replacement(donor_elements, coordination_number)
    )


def _built_in_entries() -> tuple[EmpiricalRegistryEntry, ...]:
    # These are the pre-existing Xponge experimental tables.  The ZAFF paper
    # is recorded as method context only; parameter identity is not claimed,
    # and redistribution remains restricted until the data source is audited.
    return (EmpiricalRegistryEntry(
        registry_id="xponge-zn-nos-experimental-v1",
        version="1",
        status="experimental",
        applicability=EmpiricalApplicability(
            interaction_model="bonded",
            element="Zn",
            oxidation_state=2,
            coordination_numbers=(1, 2, 3, 4, 5, 6),
            donor_elements=("N", "O", "S"),
            donor_patterns=_donor_patterns(("N", "O", "S"), (1, 2, 3, 4, 5, 6)),
            geometries=("unclassified",),
            compatible_base_force_fields=("*",),
            compatible_water_models=("*",),
        ),
        parameters={
            "bond_tables": _ZN_EMPIRICAL_BOND_TABLES,
            "angle_force_constant_default": 50.0,
            "angle_force_constant_with_sulfur": 70.0,
        },
        units={
            "bond_distance": "angstrom",
            "bond_force_constant": "kcal/mol/angstrom^2",
            "angle_force_constant": "kcal/mol/rad^2",
        },
        source="pre-existing Xponge empirical_zn_nos tables; provenance audit required",
        license="unknown",
        redistribution_status="internal-test-only",
        citations=(EmpiricalCitation(
            title="Structural Survey of Zinc-Containing Proteins and Development of the Zinc AMBER Force Field (ZAFF)",
            authors=(
                "Martin B. Peters", "Yue Yang", "Bing Wang", "Laszlo Fusti-Molnar",
                "Michael N. Weaver", "Kenneth M. Merz Jr.",
            ),
            journal="Journal of Chemical Theory and Computation",
            year=2010,
            doi="10.1021/ct1002626",
            relation="method_context_only_parameter_identity_not_claimed",
        ),),
        validation_cases=("Mokda 4EWL engineering fixture",),
    ).with_computed_hashes(),)


EMPIRICAL_REGISTRY = {
    entry.registry_id: entry
    for entry in _built_in_entries()
}


def validate_empirical_registry_entry(entry: EmpiricalRegistryEntry) -> None:
    if not entry.registry_id or not entry.version or not entry.status:
        raise ValidationError("invalid_empirical_registry_identity", entry.registry_id)
    if entry.status not in {"experimental", "supported"}:
        raise ValidationError("invalid_empirical_registry_status", entry.status, entry.registry_id)
    if entry.applicability.interaction_model != "bonded":
        raise ValidationError(
            "invalid_empirical_interaction_model",
            entry.applicability.interaction_model,
            entry.registry_id,
        )
    if (
        not entry.applicability.element
        or not entry.applicability.coordination_numbers
        or not entry.applicability.donor_elements
        or not entry.applicability.donor_patterns
        or not entry.applicability.geometries
        or not entry.applicability.compatible_base_force_fields
        or not entry.applicability.compatible_water_models
    ):
        raise ValidationError("incomplete_empirical_applicability", entry.registry_id)
    if any(value <= 0 for value in entry.applicability.coordination_numbers):
        raise ValidationError("invalid_empirical_coordination_number", entry.registry_id)
    for pattern in entry.applicability.donor_patterns:
        if (
            len(pattern) not in entry.applicability.coordination_numbers
            or tuple(sorted(pattern)) != pattern
            or set(pattern) - set(entry.applicability.donor_elements)
        ):
            raise ValidationError("invalid_empirical_donor_pattern", entry.registry_id)
    if not entry.parameters or not entry.units:
        raise ValidationError("missing_empirical_parameters", entry.registry_id)
    if not entry.source or not entry.license or not entry.redistribution_status or not entry.citations:
        raise ValidationError("missing_empirical_provenance", entry.registry_id)
    if not entry.source_hash or entry.source_hash != entry.computed_source_hash():
        raise ValidationError("stale_empirical_source_hash", entry.registry_id)
    if not entry.parameter_hash or entry.parameter_hash != entry.computed_parameter_hash():
        raise ValidationError("stale_empirical_parameter_hash", entry.registry_id)


def empirical_registry_descriptor(registry_id: str) -> Mapping[str, Any]:
    entry = EMPIRICAL_REGISTRY.get(registry_id)
    if entry is None:
        raise ValidationError("unknown_empirical_registry_id", registry_id)
    validate_empirical_registry_entry(entry)
    return {
        **entry.source_payload(),
        "applicability": _canonicalize(entry.applicability),
        "source_hash": entry.source_hash,
        "parameter_hash": entry.parameter_hash,
    }


def empirical_registry_descriptors() -> tuple[Mapping[str, Any], ...]:
    """Return deterministic, provenance-complete descriptors for discovery."""

    return tuple(
        empirical_registry_descriptor(registry_id)
        for registry_id in sorted(EMPIRICAL_REGISTRY)
    )


def matching_empirical_registry_descriptors(
    *,
    interaction_model: str,
    element: str,
    oxidation_state: int,
    coordination_number: int,
    donor_elements: tuple[str, ...],
    geometry: str = "unclassified",
    base_force_field: str = "",
    water_model: str = "",
) -> tuple[Mapping[str, Any], ...]:
    """List registry records whose declared applicability exactly matches a site.

    This discovery API uses only source-neutral site facts.  It never selects a
    record or falls back to an element-specific default.
    """

    normalized_donors = tuple(sorted(str(value) for value in donor_elements))
    matches: list[Mapping[str, Any]] = []
    for registry_id in sorted(EMPIRICAL_REGISTRY):
        entry = EMPIRICAL_REGISTRY[registry_id]
        validate_empirical_registry_entry(entry)
        applicability = entry.applicability
        if (
            interaction_model != applicability.interaction_model
            or element != applicability.element
            or oxidation_state != applicability.oxidation_state
            or coordination_number not in applicability.coordination_numbers
            or normalized_donors not in applicability.donor_patterns
            or geometry not in applicability.geometries
            or (
                "*" not in applicability.compatible_base_force_fields
                and base_force_field not in applicability.compatible_base_force_fields
            )
            or (
                "*" not in applicability.compatible_water_models
                and water_model not in applicability.compatible_water_models
            )
        ):
            continue
        matches.append(empirical_registry_descriptor(registry_id))
    return tuple(matches)


def resolve_empirical_registry(
    topology: PreparedChemicalTopology,
    *,
    registry_id: str,
    geometry: str = "unclassified",
    base_force_field: str = "",
    water_model: str = "",
) -> EmpiricalRegistryMatch:
    entry = EMPIRICAL_REGISTRY.get(registry_id)
    if entry is None:
        raise ValidationError("unknown_empirical_registry_id", registry_id)
    validate_empirical_registry_entry(entry)
    applicability = entry.applicability
    if geometry not in applicability.geometries:
        raise ValidationError("empirical_geometry_mismatch", geometry, registry_id)
    if (
        "*" not in applicability.compatible_base_force_fields
        and base_force_field not in applicability.compatible_base_force_fields
    ):
        raise ValidationError("empirical_base_force_field_mismatch", base_force_field, registry_id)
    if (
        "*" not in applicability.compatible_water_models
        and water_model not in applicability.compatible_water_models
    ):
        raise ValidationError("empirical_water_model_mismatch", water_model, registry_id)

    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    centers: list[EmpiricalCenterMatch] = []
    for metal in sorted(
        (atom for atom in topology.atoms if atom.is_metal),
        key=lambda atom: (atom.stable_order, atom.external_id),
    ):
        if metal.element != applicability.element:
            raise ValidationError("empirical_element_mismatch", metal.element, metal.external_id)
        if not metal.formal_charge_known or metal.formal_charge != applicability.oxidation_state:
            raise ValidationError(
                "empirical_oxidation_state_mismatch",
                str(metal.formal_charge),
                metal.external_id,
            )
        donor_ids: list[str] = []
        for edge in (*topology.bonds, *topology.links):
            if not edge.active or metal.external_id not in edge.atom_ids:
                continue
            donor_ids.append(edge.atom_ids[1] if edge.atom_ids[0] == metal.external_id else edge.atom_ids[0])
        donor_ids = sorted(set(donor_ids))
        donor_elements = tuple(sorted(atom_by_id[atom_id].element for atom_id in donor_ids))
        if len(donor_ids) not in applicability.coordination_numbers:
            raise ValidationError(
                "empirical_coordination_number_mismatch",
                str(len(donor_ids)),
                metal.external_id,
            )
        unsupported_donors = sorted(set(donor_elements) - set(applicability.donor_elements))
        if unsupported_donors:
            raise ValidationError(
                "empirical_donor_pattern_mismatch",
                ",".join(unsupported_donors),
                metal.external_id,
            )
        if donor_elements not in applicability.donor_patterns:
            raise ValidationError(
                "empirical_donor_pattern_mismatch",
                ",".join(donor_elements),
                metal.external_id,
            )
        centers.append(EmpiricalCenterMatch(
            metal_atom_id=metal.external_id,
            element=metal.element,
            oxidation_state=int(metal.formal_charge),
            coordination_number=len(donor_ids),
            donor_elements=donor_elements,
        ))
    if not centers:
        raise ValidationError("missing_empirical_metal_center", registry_id)
    return EmpiricalRegistryMatch(
        entry=entry,
        centers=tuple(centers),
        geometry=geometry,
        base_force_field=base_force_field,
        water_model=water_model,
    )


__all__ = [
    "EMPIRICAL_REGISTRY", "EmpiricalApplicability", "EmpiricalCenterMatch", "EmpiricalCitation",
    "EmpiricalRegistryEntry", "EmpiricalRegistryMatch", "empirical_registry_descriptor",
    "empirical_registry_descriptors", "matching_empirical_registry_descriptors",
    "resolve_empirical_registry",
    "validate_empirical_registry_entry",
]
