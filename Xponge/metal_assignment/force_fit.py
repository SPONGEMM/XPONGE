"""Topology-preserving metal-related force-constant providers.

The functions in this module are pure with respect to Xponge's global force-
field registries.  They consume locked topology/model records and return
external-ID keyed overlay terms; they never construct or patch a parent model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

from .artifacts import DerivedModel
from .contracts import PreparedChemicalTopology, ValidationError, _freeze, _strict_object
from .empirical_registry import resolve_empirical_registry


HARTREE_TO_KCAL_MOL = 627.5094740631
BOHR_TO_ANGSTROM = 0.529177210903
BOND_FORCE_CONVERSION = HARTREE_TO_KCAL_MOL / (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM)
ANGLE_FORCE_CONVERSION = HARTREE_TO_KCAL_MOL


@dataclass(frozen=True, slots=True)
class HessianArtifact:
    model_id: str
    model_hash: str
    atom_order: tuple[str, ...]
    coordinates_angstrom: tuple[tuple[float, float, float], ...]
    cartesian_hessian_au: Any
    provider: str
    provider_version: str


def _finite_coordinates(value: Sequence[Sequence[float]], expected: int) -> None:
    if len(value) != expected:
        raise ValidationError("hessian_coordinate_count_mismatch", str(len(value)))
    for index, coordinate in enumerate(value):
        if len(coordinate) != 3 or not all(
            isinstance(item, (int, float)) and not isinstance(item, bool) and math.isfinite(item)
            for item in coordinate
        ):
            raise ValidationError("invalid_hessian_coordinate", str(index))


def _flatten_hessian(value: Any, atom_count: int):
    import numpy as np

    hessian = np.asarray(value, dtype=float)
    if hessian.ndim == 4 and hessian.shape == (atom_count, atom_count, 3, 3):
        hessian = np.transpose(hessian, (0, 2, 1, 3)).reshape(atom_count * 3, atom_count * 3)
    if hessian.ndim != 2 or hessian.shape != (atom_count * 3, atom_count * 3):
        raise ValidationError("invalid_hessian_shape", repr(hessian.shape))
    if not np.isfinite(hessian).all():
        raise ValidationError("nonfinite_hessian", "Hessian contains non-finite values")
    if not np.allclose(hessian, hessian.T, atol=1.0e-8, rtol=1.0e-7):
        raise ValidationError("asymmetric_hessian", "Cartesian Hessian must be symmetric")
    return hessian


def _wire_strings(value: Any, path: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValidationError("invalid_wire_type", "expected string array", path)
    return tuple(value)


def hessian_artifact_from_dict(value: Any, *, path: str = "hessian_artifact") -> HessianArtifact:
    data = _strict_object(
        value,
        required={
            "model_id", "model_hash", "atom_order", "coordinates_angstrom",
            "cartesian_hessian_au", "provider", "provider_version",
        },
        path=path,
    )
    if not isinstance(data["coordinates_angstrom"], list) or not isinstance(data["cartesian_hessian_au"], list):
        raise ValidationError("invalid_wire_type", "coordinates and Hessian must be arrays", path)
    return HessianArtifact(
        model_id=data["model_id"],
        model_hash=data["model_hash"],
        atom_order=_wire_strings(data["atom_order"], f"{path}.atom_order"),
        coordinates_angstrom=tuple(tuple(coordinate) for coordinate in data["coordinates_angstrom"]),
        cartesian_hessian_au=_freeze(data["cartesian_hessian_au"]),
        provider=data["provider"],
        provider_version=data["provider_version"],
    )


def validate_hessian_artifact(model: DerivedModel, artifact: HessianArtifact) -> None:
    if model.purpose != "small":
        raise ValidationError("invalid_hessian_model_purpose", model.purpose)
    if artifact.model_id != model.external_id or artifact.model_hash != model.model_hash:
        raise ValidationError("hessian_model_identity_mismatch", artifact.model_id)
    if artifact.atom_order != tuple(atom.model_atom_id for atom in model.atoms):
        raise ValidationError("hessian_atom_order_mismatch", artifact.model_id)
    if not artifact.provider or not artifact.provider_version:
        raise ValidationError("missing_hessian_provenance", artifact.model_id)
    _finite_coordinates(artifact.coordinates_angstrom, len(model.atoms))

    import numpy as np

    artifact_coordinates = np.asarray(artifact.coordinates_angstrom, dtype=float)
    model_coordinates = np.asarray([atom.coordinates for atom in model.atoms], dtype=float)
    if not np.allclose(artifact_coordinates, model_coordinates, atol=1.0e-8, rtol=0.0):
        raise ValidationError("hessian_coordinate_identity_mismatch", artifact.model_id)
    _flatten_hessian(artifact.cartesian_hessian_au, len(model.atoms))


def _bond_b_vector(coords_bohr, atom1: int, atom2: int):
    import numpy as np

    vector = coords_bohr[atom1] - coords_bohr[atom2]
    distance = float(np.linalg.norm(vector))
    if distance <= 1.0e-12:
        raise ValidationError("invalid_force_fit_geometry", "bond atoms coincide")
    unit = vector / distance
    b_vector = np.zeros((coords_bohr.shape[0], 3), dtype=float)
    b_vector[atom1] = unit
    b_vector[atom2] = -unit
    return b_vector.reshape(-1), distance * BOHR_TO_ANGSTROM


def _angle_b_vector(coords_bohr, atom1: int, atom2: int, atom3: int):
    import numpy as np

    vector21 = coords_bohr[atom1] - coords_bohr[atom2]
    vector23 = coords_bohr[atom3] - coords_bohr[atom2]
    norm21 = float(np.linalg.norm(vector21))
    norm23 = float(np.linalg.norm(vector23))
    if norm21 <= 1.0e-12 or norm23 <= 1.0e-12:
        raise ValidationError("invalid_force_fit_geometry", "angle contains a zero-length bond")
    unit21 = vector21 / norm21
    unit23 = vector23 / norm23
    cosine = float(np.clip(np.dot(unit21, unit23), -1.0, 1.0))
    theta = math.acos(cosine)
    sine = math.sin(theta)
    if abs(sine) <= 1.0e-10:
        raise ValidationError("invalid_force_fit_geometry", "linear angle cannot be projected")
    gradient1 = (unit21 * cosine - unit23) / (norm21 * sine)
    gradient3 = (unit23 * cosine - unit21) / (norm23 * sine)
    gradient2 = -(gradient1 + gradient3)
    b_vector = np.zeros((coords_bohr.shape[0], 3), dtype=float)
    b_vector[atom1] = gradient1
    b_vector[atom2] = gradient2
    b_vector[atom3] = gradient3
    return b_vector.reshape(-1), theta


def _project_force_constant(hessian, b_vector) -> float:
    import numpy as np

    norm_squared = float(np.dot(b_vector, b_vector))
    if norm_squared <= 1.0e-16:
        raise ValidationError("invalid_force_fit_coordinate", "internal-coordinate vector has zero norm")
    projected = float(np.dot(b_vector, hessian @ b_vector) / (norm_squared * norm_squared))
    if not math.isfinite(projected):
        raise ValidationError("nonfinite_force_constant", "projected force constant is non-finite")
    return max(0.0, projected)


def seminario_bonded_terms(
    model: DerivedModel,
    artifact: HessianArtifact,
    *,
    scale_factor: float = 1.0,
) -> tuple[dict[str, Mapping[str, Any]], Mapping[str, Any]]:
    """Project metal bond/angle terms from an explicit model Hessian."""

    validate_hessian_artifact(model, artifact)
    if isinstance(scale_factor, bool) or not isinstance(scale_factor, (int, float)) or not math.isfinite(scale_factor) or scale_factor <= 0:
        raise ValidationError("invalid_force_constant_scale", str(scale_factor))
    hessian = _flatten_hessian(artifact.cartesian_hessian_au, len(model.atoms))

    import numpy as np

    coordinates_bohr = np.asarray(artifact.coordinates_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    index_by_model_id = {atom.model_atom_id: index for index, atom in enumerate(model.atoms)}
    atom_by_model_id = {atom.model_atom_id: atom for atom in model.atoms}
    core_external_by_model_id = {
        atom.model_atom_id: atom.external_id for atom in model.atoms if atom.role != "cap"
    }
    metal_model_ids = {
        atom.model_atom_id for atom in model.atoms
        if atom.role != "cap" and atom.element in {"Li", "Na", "K", "Mg", "Ca", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
    }
    neighbor_ids: dict[str, set[str]] = {model_atom_id: set() for model_atom_id in metal_model_ids}
    terms: dict[str, Mapping[str, Any]] = {}
    for edge in (*model.bonds, *model.links):
        atom1, atom2 = edge.model_atom_ids
        if atom1 not in core_external_by_model_id or atom2 not in core_external_by_model_id:
            continue
        if not ({atom1, atom2} & metal_model_ids):
            continue
        b_vector, equilibrium = _bond_b_vector(
            coordinates_bohr,
            index_by_model_id[atom1],
            index_by_model_id[atom2],
        )
        force_constant = _project_force_constant(hessian, b_vector) * BOND_FORCE_CONVERSION * scale_factor
        external_ids = (str(core_external_by_model_id[atom1]), str(core_external_by_model_id[atom2]))
        term_id = f"seminario:bond:{edge.external_id}"
        terms[term_id] = {
            "kind": "bond",
            "atom_ids": external_ids,
            "parameters": {
                "k": round(force_constant, 8),
                "equilibrium": round(equilibrium, 8),
                "k_unit": "kcal/mol/angstrom^2",
                "equilibrium_unit": "angstrom",
            },
            "source": f"seminario:{artifact.provider}:{artifact.provider_version}",
        }
        if atom1 in metal_model_ids:
            neighbor_ids[atom1].add(atom2)
        if atom2 in metal_model_ids:
            neighbor_ids[atom2].add(atom1)

    for metal_model_id, neighbors in sorted(neighbor_ids.items()):
        ordered_neighbors = sorted(neighbors, key=lambda atom_id: str(core_external_by_model_id[atom_id]))
        for index, atom1 in enumerate(ordered_neighbors):
            for atom3 in ordered_neighbors[index + 1:]:
                b_vector, equilibrium = _angle_b_vector(
                    coordinates_bohr,
                    index_by_model_id[atom1],
                    index_by_model_id[metal_model_id],
                    index_by_model_id[atom3],
                )
                force_constant = _project_force_constant(hessian, b_vector) * ANGLE_FORCE_CONVERSION * scale_factor
                external_ids = (
                    str(core_external_by_model_id[atom1]),
                    str(core_external_by_model_id[metal_model_id]),
                    str(core_external_by_model_id[atom3]),
                )
                canonical_ends = sorted((external_ids[0], external_ids[2]))
                term_id = f"seminario:angle:{canonical_ends[0]}:{external_ids[1]}:{canonical_ends[1]}"
                terms[term_id] = {
                    "kind": "angle",
                    "atom_ids": external_ids,
                    "parameters": {
                        "k": round(force_constant, 8),
                        "equilibrium": round(equilibrium, 8),
                        "k_unit": "kcal/mol/rad^2",
                        "equilibrium_unit": "rad",
                    },
                    "source": f"seminario:{artifact.provider}:{artifact.provider_version}",
                }
    report = {
        "provider": "seminario",
        "provider_revision": "external-hessian-v1",
        "hessian_provider": artifact.provider,
        "hessian_provider_version": artifact.provider_version,
        "model_id": model.external_id,
        "model_hash": model.model_hash,
        "term_count": len(terms),
        "scale_factor": float(scale_factor),
        "cap_terms_discarded": True,
    }
    return terms, report


def _distance(atom1, atom2) -> float:
    return math.sqrt(sum((float(left) - float(right)) ** 2 for left, right in zip(atom1.coordinates, atom2.coordinates)))


def _angle(atom1, atom2, atom3) -> float:
    vector1 = tuple(float(left) - float(center) for left, center in zip(atom1.coordinates, atom2.coordinates))
    vector2 = tuple(float(right) - float(center) for right, center in zip(atom3.coordinates, atom2.coordinates))
    norm1 = math.sqrt(sum(value * value for value in vector1))
    norm2 = math.sqrt(sum(value * value for value in vector2))
    if norm1 <= 1.0e-12 or norm2 <= 1.0e-12:
        raise ValidationError("invalid_force_fit_geometry", "empirical angle contains a zero-length bond")
    cosine = max(-1.0, min(1.0, sum(left * right for left, right in zip(vector1, vector2)) / norm1 / norm2))
    return math.acos(cosine)


def manual_bonded_terms(
    topology: PreparedChemicalTopology,
    *,
    bond_force_constant: float,
    angle_force_constant: float,
    reference_geometry_artifact: Mapping[str, Any],
    site_force_constants: Mapping[str, Any] | None = None,
) -> tuple[dict[str, Mapping[str, Any]], Mapping[str, Any]]:
    """Project uniform constants onto an explicit, topology-locked geometry."""

    for name, value in (
        ("bond_force_constant", bond_force_constant),
        ("angle_force_constant", angle_force_constant),
    ):
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or value <= 0
            or value > 10000
        ):
            raise ValidationError("invalid_manual_force_constant", name)
    artifact = _strict_object(
        reference_geometry_artifact,
        required={
            "schema_version", "graph_revision", "input_hash",
            "coordinate_unit", "angle_unit", "geometry_source", "selections",
            "artifact_hash",
        },
        path="reference_geometry_artifact",
    )
    if (
        artifact["schema_version"] != 1
        or artifact["graph_revision"] != topology.graph_revision
        or artifact["input_hash"] != topology.input_hash
        or artifact["coordinate_unit"] != "angstrom"
        or artifact["angle_unit"] != "rad"
        or artifact["geometry_source"] != "frozen_current_geometry"
        or not isinstance(artifact["selections"], (list, tuple))
        or not artifact["selections"]
    ):
        raise ValidationError(
            "reference_geometry_identity_mismatch",
            "reference geometry does not match the prepared topology",
        )
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    expected_edges: dict[str, tuple[str, str]] = {}
    for edge in topology.bonds:
        if edge.active and edge.semantic == "coordination":
            expected_edges[edge.external_id] = edge.atom_ids
    for edge in topology.links:
        if edge.active and edge.confirmed and edge.kind == "coordination":
            if edge.external_id in expected_edges:
                raise ValidationError(
                    "duplicate_manual_coordination_link",
                    edge.external_id,
                )
            expected_edges[edge.external_id] = edge.atom_ids
    if not expected_edges:
        raise ValidationError(
            "missing_manual_coordination_terms",
            "no confirmed active coordination links",
        )
    expected_selection_ids = {
        selection.get("selection_id")
        for selection in artifact["selections"]
        if isinstance(selection, Mapping)
    }
    if site_force_constants is not None and set(site_force_constants) != expected_selection_ids:
        raise ValidationError(
            "manual_site_force_constant_coverage_mismatch",
            f"expected={sorted(expected_selection_ids)},"
            f"actual={sorted(site_force_constants)}",
        )

    seen_edges: set[str] = set()
    seen_angle_pairs: set[tuple[str, str]] = set()
    terms: dict[str, Mapping[str, Any]] = {}
    expected_angle_count = 0
    for selection_index, raw_selection in enumerate(artifact["selections"]):
        selection = _strict_object(
            raw_selection,
            required={"selection_id", "center_external_id", "bonds", "angles"},
            path=f"reference_geometry_artifact.selections[{selection_index}]",
        )
        metal_id = selection["center_external_id"]
        selection_constants = (
            site_force_constants.get(selection["selection_id"], {})
            if site_force_constants is not None
            else {}
        )
        site_bond_force_constant = float(
            selection_constants.get("bond_force_constant", bond_force_constant)
        )
        site_angle_force_constant = float(
            selection_constants.get("angle_force_constant", angle_force_constant)
        )
        metal = atom_by_id.get(metal_id)
        if metal is None or not metal.is_metal:
            raise ValidationError("invalid_reference_geometry_center", str(metal_id))
        if not isinstance(selection["bonds"], (list, tuple)):
            raise ValidationError("invalid_reference_geometry_bonds", str(metal_id))
        donor_by_edge: dict[str, str] = {}
        for raw_bond in selection["bonds"]:
            bond = _strict_object(
                raw_bond,
                required={
                    "edge_id", "center_external_id", "neighbor_external_id",
                    "equilibrium", "unit",
                },
                path="reference_geometry_artifact.bond",
            )
            edge_id = bond["edge_id"]
            donor_id = bond["neighbor_external_id"]
            equilibrium = bond["equilibrium"]
            if (
                edge_id in seen_edges
                or bond["center_external_id"] != metal_id
                or bond["unit"] != "angstrom"
                or donor_id not in atom_by_id
                or atom_by_id[donor_id].is_metal
                or isinstance(equilibrium, bool)
                or not isinstance(equilibrium, (int, float))
                or not math.isfinite(equilibrium)
                or equilibrium <= 1.0e-12
                or set(expected_edges.get(edge_id, ())) != {metal_id, donor_id}
            ):
                raise ValidationError("invalid_reference_geometry_bond", str(edge_id))
            seen_edges.add(edge_id)
            donor_by_edge[edge_id] = donor_id
            terms[f"manual:bond:{edge_id}"] = {
                "kind": "bond",
                "atom_ids": (metal_id, donor_id),
                "parameters": {
                    "k": site_bond_force_constant,
                    "equilibrium": round(float(equilibrium), 8),
                    "k_unit": "kcal/mol/angstrom^2",
                    "equilibrium_unit": "angstrom",
                },
                "source": "manual_bonded:explicit_reference_geometry",
            }
        expected_angle_count += len(donor_by_edge) * (len(donor_by_edge) - 1) // 2
        if not isinstance(selection["angles"], (list, tuple)):
            raise ValidationError("invalid_reference_geometry_angles", str(metal_id))
        for raw_angle in selection["angles"]:
            angle = _strict_object(
                raw_angle,
                required={
                    "edge_id1", "edge_id2", "neighbor1_external_id",
                    "center_external_id", "neighbor2_external_id",
                    "equilibrium", "unit",
                },
                path="reference_geometry_artifact.angle",
            )
            pair = tuple(sorted((angle["edge_id1"], angle["edge_id2"])))
            equilibrium = angle["equilibrium"]
            if (
                pair in seen_angle_pairs
                or pair[0] == pair[1]
                or angle["center_external_id"] != metal_id
                or angle["unit"] != "rad"
                or donor_by_edge.get(angle["edge_id1"]) != angle["neighbor1_external_id"]
                or donor_by_edge.get(angle["edge_id2"]) != angle["neighbor2_external_id"]
                or isinstance(equilibrium, bool)
                or not isinstance(equilibrium, (int, float))
                or not math.isfinite(equilibrium)
                or equilibrium <= 0
                or equilibrium > math.pi
            ):
                raise ValidationError("invalid_reference_geometry_angle", ":".join(pair))
            seen_angle_pairs.add(pair)
            terms[f"manual:angle:{pair[0]}:{pair[1]}"] = {
                "kind": "angle",
                "atom_ids": (
                    angle["neighbor1_external_id"],
                    metal_id,
                    angle["neighbor2_external_id"],
                ),
                "parameters": {
                    "k": site_angle_force_constant,
                    "equilibrium": round(float(equilibrium), 8),
                    "k_unit": "kcal/mol/rad^2",
                    "equilibrium_unit": "rad",
                },
                "source": "manual_bonded:explicit_reference_geometry",
            }
    if seen_edges != set(expected_edges) or len(seen_angle_pairs) != expected_angle_count:
        raise ValidationError(
            "reference_geometry_term_coverage_mismatch",
            f"bonds={len(seen_edges)}/{len(expected_edges)},"
            f"angles={len(seen_angle_pairs)}/{expected_angle_count}",
        )
    return terms, {
        "provider": "manual_bonded",
        "provider_revision": "explicit-reference-geometry-v1",
        "equilibrium_geometry_source": "frozen_current_geometry",
        "bond_force_constant": float(bond_force_constant),
        "bond_force_constant_unit": "kcal/mol/angstrom^2",
        "angle_force_constant": float(angle_force_constant),
        "angle_force_constant_unit": "kcal/mol/rad^2",
        "potential_convention": "E=k*delta^2",
        "site_force_constants": (
            {
                key: {
                    "bond_force_constant": float(value["bond_force_constant"]),
                    "angle_force_constant": float(value["angle_force_constant"]),
                }
                for key, value in sorted(site_force_constants.items())
            }
            if site_force_constants is not None
            else {}
        ),
        "metal_site_count": len(artifact["selections"]),
        "bond_term_count": len(seen_edges),
        "angle_term_count": len(seen_angle_pairs),
        "term_count": len(terms),
    }


def _interpolate(points: Sequence[tuple[float, float]], distance: float) -> float:
    if distance <= points[0][0]:
        return float(points[0][1])
    if distance >= points[-1][0]:
        return float(points[-1][1])
    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        if x1 <= distance <= x2:
            return float(y1 + (distance - x1) * (y2 - y1) / (x2 - x1))
    raise ValidationError("empirical_interpolation_failure", str(distance))


def empirical_registry_bonded_terms(
    topology: PreparedChemicalTopology,
    *,
    registry_id: str,
    geometry: str = "unclassified",
    base_force_field: str = "",
    water_model: str = "",
) -> tuple[dict[str, Mapping[str, Any]], Mapping[str, Any]]:
    """Return an exact-match empirical overlay from Xponge's registry."""

    match = resolve_empirical_registry(
        topology,
        registry_id=registry_id,
        geometry=geometry,
        base_force_field=base_force_field,
        water_model=water_model,
    )
    parameters = match.entry.parameters
    bond_tables = parameters["bond_tables"]
    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    metal_ids = {center.metal_atom_id for center in match.centers}
    neighbors: dict[str, list[str]] = {atom_id: [] for atom_id in metal_ids}
    terms: dict[str, Mapping[str, Any]] = {}
    for edge in (*topology.bonds, *topology.links):
        if not edge.active:
            continue
        metals = set(edge.atom_ids) & metal_ids
        if not metals:
            continue
        if len(metals) != 1:
            raise ValidationError("unsupported_empirical_metal_bridge", edge.external_id)
        metal_id = next(iter(metals))
        donor_id = edge.atom_ids[1] if edge.atom_ids[0] == metal_id else edge.atom_ids[0]
        donor = atom_by_id[donor_id]
        table = bond_tables.get(donor.element)
        if table is None:
            raise ValidationError(
                "unsupported_empirical_donor",
                f"{atom_by_id[metal_id].element}-{donor.element}",
                edge.external_id,
            )
        distance = _distance(atom_by_id[metal_id], donor)
        terms[f"empirical:bond:{edge.external_id}"] = {
            "kind": "bond",
            "atom_ids": (metal_id, donor_id),
            "parameters": {
                "k": round(_interpolate(table, distance), 8),
                "equilibrium": round(distance, 8),
                "k_unit": "kcal/mol/angstrom^2",
                "equilibrium_unit": "angstrom",
            },
            "source": f"{match.entry.registry_id}:{match.entry.version}",
        }
        neighbors[metal_id].append(donor_id)

    for metal_id, donor_ids in sorted(neighbors.items()):
        ordered = sorted(set(donor_ids))
        for index, atom1_id in enumerate(ordered):
            for atom3_id in ordered[index + 1:]:
                elements = {atom_by_id[atom1_id].element, atom_by_id[atom3_id].element}
                force_constant = (
                    float(parameters["angle_force_constant_with_sulfur"])
                    if "S" in elements
                    else float(parameters["angle_force_constant_default"])
                )
                term_id = f"empirical:angle:{atom1_id}:{metal_id}:{atom3_id}"
                terms[term_id] = {
                    "kind": "angle",
                    "atom_ids": (atom1_id, metal_id, atom3_id),
                    "parameters": {
                        "k": force_constant,
                        "equilibrium": round(_angle(atom_by_id[atom1_id], atom_by_id[metal_id], atom_by_id[atom3_id]), 8),
                        "k_unit": "kcal/mol/rad^2",
                        "equilibrium_unit": "rad",
                    },
                    "source": f"{match.entry.registry_id}:{match.entry.version}",
                }
    if not terms:
        raise ValidationError("missing_empirical_terms", f"no active interaction matched {registry_id}")
    return terms, {
        "provider": "empirical_registry",
        "registry": match.descriptor(),
        "status": match.entry.status,
        "metal_count": len(metal_ids),
        "term_count": len(terms),
    }


def empirical_zn_nos_bonded_terms(
    topology: PreparedChemicalTopology,
) -> tuple[dict[str, Mapping[str, Any]], Mapping[str, Any]]:
    """Compatibility adapter for the pre-registry experimental Zn provider."""

    return empirical_registry_bonded_terms(
        topology,
        registry_id="xponge-zn-nos-experimental-v1",
    )


__all__ = [
    "ANGLE_FORCE_CONVERSION", "BOHR_TO_ANGSTROM", "BOND_FORCE_CONVERSION", "HARTREE_TO_KCAL_MOL",
    "HessianArtifact", "empirical_registry_bonded_terms", "empirical_zn_nos_bonded_terms",
    "manual_bonded_terms",
    "hessian_artifact_from_dict", "seminario_bonded_terms", "validate_hessian_artifact",
]
