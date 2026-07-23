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


_ZN_EMPIRICAL_BOND_TABLES = {
    ("N", "Zn"): (
        (1.926, 124.58), (1.947, 113.59), (1.978, 93.97), (1.982, 90.76), (1.984, 90.80),
        (2.011, 78.18), (2.027, 70.85), (2.028, 70.86), (2.029, 66.61), (2.040, 66.69),
        (2.041, 66.11), (2.047, 62.61), (2.073, 52.77), (2.089, 50.10), (2.133, 36.30),
        (2.145, 35.80), (2.176, 29.24),
    ),
    ("O", "Zn"): (
        (1.860, 169.29), (1.986, 76.81), (2.011, 71.26), (2.054, 56.37), (2.109, 41.86),
        (2.112, 41.32),
    ),
    ("S", "Zn"): (
        (2.262, 88.50), (2.305, 67.39), (2.353, 50.99), (2.355, 51.79), (2.426, 32.69),
    ),
}


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


def _interpolate(points: Sequence[tuple[float, float]], distance: float) -> float:
    if distance <= points[0][0]:
        return float(points[0][1])
    if distance >= points[-1][0]:
        return float(points[-1][1])
    for (x1, y1), (x2, y2) in zip(points, points[1:]):
        if x1 <= distance <= x2:
            return float(y1 + (distance - x1) * (y2 - y1) / (x2 - x1))
    raise ValidationError("empirical_interpolation_failure", str(distance))


def empirical_zn_nos_bonded_terms(
    topology: PreparedChemicalTopology,
) -> tuple[dict[str, Mapping[str, Any]], Mapping[str, Any]]:
    """Return the explicitly limited experimental Zn-N/O/S parameter overlay."""

    atom_by_id = {atom.external_id: atom for atom in topology.atoms}
    zinc_ids = {atom.external_id for atom in topology.atoms if atom.is_metal and atom.element == "Zn"}
    other_metal_ids = {atom.external_id for atom in topology.atoms if atom.is_metal and atom.element != "Zn"}
    if other_metal_ids:
        raise ValidationError("unsupported_empirical_metal", ",".join(sorted(other_metal_ids)))
    neighbors: dict[str, list[str]] = {atom_id: [] for atom_id in zinc_ids}
    terms: dict[str, Mapping[str, Any]] = {}
    for edge in (*topology.bonds, *topology.links):
        if not edge.active:
            continue
        metals = set(edge.atom_ids) & zinc_ids
        if not metals:
            continue
        if len(metals) != 1:
            raise ValidationError("unsupported_empirical_metal_bridge", edge.external_id)
        metal_id = next(iter(metals))
        donor_id = edge.atom_ids[1] if edge.atom_ids[0] == metal_id else edge.atom_ids[0]
        donor = atom_by_id[donor_id]
        table = _ZN_EMPIRICAL_BOND_TABLES.get(tuple(sorted(("Zn", donor.element))))
        if table is None:
            raise ValidationError("unsupported_empirical_donor", f"Zn-{donor.element}", edge.external_id)
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
            "source": "empirical_zn_nos_v1",
        }
        neighbors[metal_id].append(donor_id)

    for metal_id, donor_ids in sorted(neighbors.items()):
        ordered = sorted(set(donor_ids))
        for index, atom1_id in enumerate(ordered):
            for atom3_id in ordered[index + 1:]:
                elements = {atom_by_id[atom1_id].element, atom_by_id[atom3_id].element}
                force_constant = 70.0 if "S" in elements else 50.0
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
                    "source": "empirical_zn_nos_v1",
                }
    if not terms:
        raise ValidationError("missing_empirical_terms", "no active Zn-N/O/S interaction was found")
    return terms, {
        "provider": "empirical_zn_nos_v1",
        "status": "experimental",
        "metal_count": len(zinc_ids),
        "term_count": len(terms),
    }


__all__ = [
    "ANGLE_FORCE_CONVERSION", "BOHR_TO_ANGSTROM", "BOND_FORCE_CONVERSION", "HARTREE_TO_KCAL_MOL",
    "HessianArtifact", "empirical_zn_nos_bonded_terms", "hessian_artifact_from_dict",
    "seminario_bonded_terms", "validate_hessian_artifact",
]
