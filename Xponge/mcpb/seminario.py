"""First-pass Seminario-style frcmod generation helpers."""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

from ..qm import compute_hessian
from ._compat import atom_type_name
from .charge_refit import _local_spin, _local_total_charge, build_local_resp_assignment

HARTREE_TO_KCAL_MOL = 627.5094740631
BOHR_TO_ANGSTROM = 0.529177210903
BOND_FORCE_CONVERSION = HARTREE_TO_KCAL_MOL / (BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM)
ANGLE_FORCE_CONVERSION = HARTREE_TO_KCAL_MOL


def _flatten_cartesian_hessian(cartesian_hessian_au):
    import numpy as np

    hessian = np.asarray(cartesian_hessian_au, dtype=float)
    if hessian.ndim == 4:
        natoms = int(hessian.shape[0])
        return np.transpose(hessian, (0, 2, 1, 3)).reshape(natoms * 3, natoms * 3)
    if hessian.ndim == 2:
        return hessian
    raise ValueError(f"unexpected Hessian shape: {hessian.shape!r}")


def _bond_b_vector(coords_bohr, atom1: int, atom2: int):
    import numpy as np

    vector = coords_bohr[atom1] - coords_bohr[atom2]
    distance = float(np.linalg.norm(vector))
    if distance <= 1.0e-12:
        raise ValueError("cannot build Seminario bond coordinate for coincident atoms")
    unit = vector / distance
    b_vector = np.zeros((coords_bohr.shape[0], 3), dtype=float)
    b_vector[atom1] = unit
    b_vector[atom2] = -unit
    return b_vector.reshape(-1), distance * BOHR_TO_ANGSTROM


def _angle_b_vector(coords_bohr, atom1: int, atom2: int, atom3: int):
    import numpy as np

    vec21 = coords_bohr[atom1] - coords_bohr[atom2]
    vec23 = coords_bohr[atom3] - coords_bohr[atom2]
    norm21 = float(np.linalg.norm(vec21))
    norm23 = float(np.linalg.norm(vec23))
    if norm21 <= 1.0e-12 or norm23 <= 1.0e-12:
        raise ValueError("cannot build Seminario angle coordinate with zero-length bond")
    unit21 = vec21 / norm21
    unit23 = vec23 / norm23
    cosine = float(np.clip(np.dot(unit21, unit23), -1.0, 1.0))
    theta = math.acos(cosine)
    sine = math.sin(theta)
    if abs(sine) <= 1.0e-10:
        raise ValueError("cannot build Seminario angle coordinate for a linear triplet")
    grad1 = (unit21 * cosine - unit23) / (norm21 * sine)
    grad3 = (unit23 * cosine - unit21) / (norm23 * sine)
    grad2 = -(grad1 + grad3)
    b_vector = np.zeros((coords_bohr.shape[0], 3), dtype=float)
    b_vector[atom1] = grad1
    b_vector[atom2] = grad2
    b_vector[atom3] = grad3
    return b_vector.reshape(-1), math.degrees(theta)


def _project_force_constant(cartesian_hessian_au, b_vector):
    import numpy as np

    norm_sq = float(np.dot(b_vector, b_vector))
    if norm_sq <= 1.0e-16:
        raise ValueError("cannot project Seminario force constant with a zero-norm internal coordinate")
    projected = float(np.dot(b_vector, cartesian_hessian_au @ b_vector) / (norm_sq * norm_sq))
    return max(0.0, projected)


def _sorted_type_pair(type1: str, type2: str) -> tuple[str, str]:
    return (type1, type2) if type1 <= type2 else (type2, type1)


def _neighbors_by_ion(selection):
    neighbors = {atom_id: [] for atom_id in selection.ion_atom_ids}
    for atom1, atom2 in selection.bonded_pairs:
        ion_id = atom1 if atom1 in neighbors else atom2
        neighbor_id = atom2 if ion_id == atom1 else atom1
        neighbors[ion_id].append(neighbor_id)
    return neighbors


def _compute_seminario_summary(request, selection, local_model):
    assignment = build_local_resp_assignment(request, local_model)
    total_charge = _local_total_charge(request, local_model)
    spin = _local_spin(request, local_model)
    hessian_result = compute_hessian(
        assignment,
        backend=request.qm_backend,
        basis=request.basis or "6-31g*",
        charge=total_charge,
        spin=spin,
        return_timings=True,
    )
    import numpy as np

    coords_bohr = np.asarray(hessian_result.coordinates_angstrom, dtype=float) / BOHR_TO_ANGSTROM
    cartesian_hessian = _flatten_cartesian_hessian(hessian_result.cartesian_hessian_au)
    neighbors_by_ion = _neighbors_by_ion(selection)
    bond_terms: dict[tuple[str, str], tuple[float, float]] = {}
    angle_terms: dict[tuple[str, str, str], tuple[float, float]] = {}
    for atom1, atom2 in selection.bonded_pairs:
        local1 = local_model.atom_id_map[atom1]
        local2 = local_model.atom_id_map[atom2]
        b_vector, distance_angstrom = _bond_b_vector(coords_bohr, local1, local2)
        force_constant = round(
            _project_force_constant(cartesian_hessian, b_vector) * BOND_FORCE_CONVERSION * float(request.scale_factor),
            1,
        )
        type1 = atom_type_name(request.molecule.atoms[atom1])
        type2 = atom_type_name(request.molecule.atoms[atom2])
        bond_terms[_sorted_type_pair(type1, type2)] = (force_constant, distance_angstrom)
    for ion_atom_id, neighbors in neighbors_by_ion.items():
        if len(neighbors) < 2:
            continue
        local_ion = local_model.atom_id_map[ion_atom_id]
        ion_type = atom_type_name(request.molecule.atoms[ion_atom_id])
        for index, atom1 in enumerate(neighbors):
            for atom3 in neighbors[index + 1:]:
                local1 = local_model.atom_id_map[atom1]
                local3 = local_model.atom_id_map[atom3]
                b_vector, theta = _angle_b_vector(coords_bohr, local1, local_ion, local3)
                force_constant = round(
                    _project_force_constant(cartesian_hessian, b_vector)
                    * ANGLE_FORCE_CONVERSION
                    * float(request.scale_factor),
                    2,
                )
                type1 = atom_type_name(request.molecule.atoms[atom1])
                type3 = atom_type_name(request.molecule.atoms[atom3])
                key = (type1, ion_type, type3)
                reverse = (type3, ion_type, type1)
                angle_terms[min(key, reverse)] = (force_constant, theta)
    return {
        "assignment": assignment,
        "total_charge": total_charge,
        "spin": spin,
        "timings": dict(hessian_result.timings),
        "bond_terms": bond_terms,
        "angle_terms": angle_terms,
    }


def register_seminario_parameters(request, selection, local_model):
    from .. import register_amber_angle_parameter, register_amber_bond_parameter

    summary = _compute_seminario_summary(request, selection, local_model)
    for (type1, type2), (force_constant, distance) in summary["bond_terms"].items():
        register_amber_bond_parameter(type1, type2, force_constant, distance)
    for atom_types, (force_constant, theta) in summary["angle_terms"].items():
        register_amber_angle_parameter(atom_types, force_constant, theta)
    return summary


def build_seminario_frcmod_text(seminario_summary) -> str:
    lines = [
        "Xponge MCPB seminario frcmod",
        "MASS",
        "",
        "BOND",
    ]
    for (type1, type2), (force_constant, length) in sorted(seminario_summary["bond_terms"].items()):
        lines.append(f"{type1:>2}-{type2:<2} {force_constant:5.1f}    {length:7.4f}")
    lines.extend(["", "ANGL"])
    for (type1, type2, type3), (force_constant, theta) in sorted(seminario_summary["angle_terms"].items()):
        lines.append(f"{type1:>2}-{type2:<2}-{type3:<2} {force_constant:5.2f}    {theta:7.2f}")
    lines.extend(["", "DIHE", "", "IMPR", "", "NONBON", ""])
    return "\n".join(lines)


def write_seminario_frcmod_artifact(*, seminario_summary, directory: str | Path | None = None) -> str:
    if directory is None:
        directory = Path(tempfile.mkdtemp(prefix="xponge_mcpb_frcmod_"))
    else:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
    path = directory / "metal_center_seminario.frcmod"
    path.write_text(build_seminario_frcmod_text(seminario_summary), encoding="utf-8")
    return str(path)
