"""Composable Metal-assignment operations for an existing Xponge Molecule.

The public boundary in this module follows Xponge's native modeling style:
ordinary force fields are imported and assigned first, then a metal-local
overlay is applied to a transactional Molecule copy.  Hash-closed package and
worker objects remain implementation details of ``parameterize``.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Any, Mapping, Sequence

from .contracts import (
    ParameterizationResult,
    ValidationError,
    _freeze_field,
)
from .bonded_fit import compose_bonded_fit
from .fit_input import validate_bonded_fit_input
from .input import MetalAssignmentPackage, validate_package


_PARAMETERIZATION_OPERATIONS = frozenset({"fit-bonded-local"})


def _bond_values(term: Any) -> tuple[float, float]:
    values = term.parameters
    if "force_constant" in values and "equilibrium_distance" in values:
        return float(values["force_constant"]), float(values["equilibrium_distance"])
    if "k" in values and "equilibrium" in values:
        return float(values["k"]), float(values["equilibrium"])
    raise ValidationError("unsupported_bond_parameter_shape", term.external_id)


def _angle_values(term: Any) -> tuple[float, float]:
    values = term.parameters
    if "force_constant" in values and "equilibrium_angle" in values:
        return float(values["force_constant"]), float(values["equilibrium_angle"])
    if "k" in values and "equilibrium" in values:
        return float(values["k"]), float(values["equilibrium"])
    raise ValidationError("unsupported_angle_parameter_shape", term.external_id)


def _proper_values(term: Any) -> tuple[list[float], list[float], list[int]]:
    values = term.parameters
    if {"force_constants", "phases", "periodicities"} <= set(values):
        ks = [float(value) for value in values["force_constants"]]
        phases = [float(value) for value in values["phases"]]
        periodicities = [int(value) for value in values["periodicities"]]
    elif {"k", "phase", "periodicity"} <= set(values):
        ks = [float(values["k"])]
        phases = [float(values["phase"])]
        periodicities = [int(values["periodicity"])]
    else:
        raise ValidationError("unsupported_proper_parameter_shape", term.external_id)
    if not ks or len(ks) != len(phases) or len(ks) != len(periodicities):
        raise ValidationError("invalid_proper_parameter_multiplicity", term.external_id)
    return ks, phases, periodicities


def _improper_values(term: Any) -> tuple[float, float, int]:
    values = term.parameters
    if {"force_constant", "phase", "periodicity"} <= set(values):
        return float(values["force_constant"]), float(values["phase"]), int(values["periodicity"])
    if {"k", "phase", "periodicity"} <= set(values):
        return float(values["k"]), float(values["phase"]), int(values["periodicity"])
    raise ValidationError("unsupported_improper_parameter_shape", term.external_id)


@dataclass(frozen=True, slots=True)
class MetalAssignmentResult:
    """Result of applying a metal-local overlay to an Xponge Molecule."""

    molecule: Any
    parameterization_result: ParameterizationResult
    atom_mapping: tuple[tuple[str, int], ...]
    application_report: Mapping[str, Any]

    def __post_init__(self) -> None:
        _freeze_field(self, "application_report")


def parameterize(
    package: MetalAssignmentPackage,
    *,
    operation: str,
    **options: Any,
) -> ParameterizationResult:
    """Return a validated metal parameterization result.

    ``fit-bonded-local`` produces only the metal-local patch for application
    to an already assigned Molecule.  Lower-level RESP/Hessian providers
    remain available from their existing modules.
    """

    validate_package(package)
    if operation not in _PARAMETERIZATION_OPERATIONS:
        raise ValidationError(
            "unsupported_molecule_parameterization_operation",
            operation,
        )
    allowed_options = {"bonded_fit_input"}
    unknown = set(options) - allowed_options
    if unknown:
        raise ValidationError(
            "unexpected_molecule_parameterization_option",
            ",".join(sorted(unknown)),
        )
    request = package.request
    if request.interaction_model != "bonded":
        raise ValidationError(
            "invalid_bonded_fit_interaction_model",
            request.interaction_model,
        )
    bonded_fit_input = options.get("bonded_fit_input")
    if bonded_fit_input is None:
        raise ValidationError(
            "missing_bonded_fit_input",
            "fit-bonded requires fit artifacts",
        )
    validate_bonded_fit_input(bonded_fit_input)
    return compose_bonded_fit(
        package,
        None,
        bonded_fit_input.metal_parameter_spec,
        charge_artifacts=bonded_fit_input.charge_artifacts,
        hessian_artifacts=bonded_fit_input.hessian_artifacts,
        force_method=bonded_fit_input.force_method,
        scale_factor=bonded_fit_input.scale_factor,
        empirical_registry_id=bonded_fit_input.empirical_registry_id,
        empirical_geometry=bonded_fit_input.empirical_geometry,
        empirical_base_force_field=bonded_fit_input.empirical_base_force_field,
        empirical_water_model=bonded_fit_input.empirical_water_model,
        manual_bond_force_constant=bonded_fit_input.manual_bond_force_constant,
        manual_angle_force_constant=bonded_fit_input.manual_angle_force_constant,
        manual_site_force_constants=bonded_fit_input.manual_site_force_constants,
        reference_geometry_artifact=bonded_fit_input.reference_geometry_artifact,
    )


def _normalized_patch_mapping(
    molecule: Any,
    patch: Any,
    atom_mapping: Mapping[str, Any] | None,
) -> tuple[tuple[str, int], ...]:
    molecule.get_atoms()
    expected_atoms = {atom.external_id: atom for atom in patch.atoms}
    if atom_mapping is None:
        atom_mapping = {
            atom.external_id: atom.stable_order
            for atom in patch.atoms
        }
    if not isinstance(atom_mapping, Mapping):
        raise ValidationError("invalid_molecule_atom_mapping", "expected mapping")
    if not set(expected_atoms) <= set(atom_mapping):
        missing = sorted(set(expected_atoms) - set(atom_mapping))
        raise ValidationError(
            "molecule_atom_mapping_coverage_mismatch",
            f"missing={missing}",
        )
    molecule_atom_index = {atom: index for index, atom in enumerate(molecule.atoms)}
    normalized: list[tuple[str, int]] = []
    for external_id, identity in expected_atoms.items():
        value = atom_mapping[external_id]
        if isinstance(value, bool):
            raise ValidationError("invalid_molecule_atom_index", external_id)
        if isinstance(value, int):
            index = value
        else:
            try:
                index = molecule_atom_index[value]
            except (KeyError, TypeError) as exc:
                raise ValidationError("molecule_atom_not_found", external_id) from exc
        if index < 0 or index >= len(molecule.atoms):
            raise ValidationError("invalid_molecule_atom_index", external_id)
        molecule_atom = molecule.atoms[index]
        observed_name = str(getattr(molecule_atom, "name", "") or "").strip()
        if observed_name != identity.atom_name:
            raise ValidationError(
                "molecule_atom_identity_mismatch",
                f"{external_id}:expected_name={identity.atom_name},observed_name={observed_name}",
            )
        observed_element = getattr(molecule_atom, "element", None)
        if observed_element is None:
            observed_element = getattr(
                getattr(molecule_atom, "type", None),
                "element",
                None,
            )
        if (
            observed_element is not None
            and str(observed_element).strip()
            and str(observed_element).strip().upper()
            != identity.element.strip().upper()
        ):
            raise ValidationError(
                "molecule_atom_identity_mismatch",
                (
                    f"{external_id}:expected_element={identity.element},"
                    f"observed_element={observed_element}"
                ),
            )
        normalized.append((external_id, index))
    indices = [index for _, index in normalized]
    if len(indices) != len(set(indices)):
        raise ValidationError("duplicate_molecule_atom_mapping", "atom indices")
    return tuple(sorted(
        normalized,
        key=lambda item: (
            expected_atoms[item[0]].stable_order,
            item[0],
        ),
    ))


def _topology_fingerprint(molecule: Any) -> tuple[Any, ...]:
    molecule.get_atoms()
    atom_index = {atom: index for index, atom in enumerate(molecule.atoms)}
    residue_rows = tuple(
        (residue.name, tuple(atom_index[atom] for atom in residue.atoms))
        for residue in molecule.residues
    )
    links = tuple(sorted(
        tuple(sorted((atom_index[link.atom1], atom_index[link.atom2])))
        for link in molecule.residue_links
    ))
    return (
        len(molecule.atoms),
        len(molecule.residues),
        residue_rows,
        links,
    )


def _validate_patch_links(
    molecule: Any,
    patch: Any,
    index_by_external_id: Mapping[str, int],
) -> None:
    for link in patch.required_links:
        atom1 = molecule.atoms[index_by_external_id[link.atom_ids[0]]]
        atom2 = molecule.atoms[index_by_external_id[link.atom_ids[1]]]
        if molecule.get_residue_link(atom1, atom2) is None:
            raise ValidationError("missing_molecule_residue_link", link.external_id)


def _get_or_create_type(type_class: Any, name: str, **parameters: Any) -> Any:
    existing = type_class.get_all_types().get(name)
    normalized = {key: value for key, value in parameters.items() if value is not None}
    if existing is not None:
        for key, value in normalized.items():
            current = getattr(existing, key)
            if isinstance(value, float):
                if not math.isclose(float(current), value, rel_tol=1e-10, abs_tol=1e-12):
                    raise ValidationError("conflicting_metal_assignment_type", name)
            elif current != value:
                raise ValidationError("conflicting_metal_assignment_type", name)
        return existing
    return type_class(name=name, **normalized)


def _short_hash(value: Any) -> str:
    def plain(item: Any) -> Any:
        if isinstance(item, Mapping):
            return {str(key): plain(child) for key, child in item.items()}
        if isinstance(item, (tuple, list)):
            return [plain(child) for child in item]
        return item

    payload = json.dumps(
        plain(value),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:12]


def _get_or_create_proper_type(
    proper_type: Any,
    name: str,
    ks: list[float],
    phases: list[float],
    periodicities: list[int],
) -> Any:
    existing = proper_type.get_all_types().get(name)
    if existing is not None:
        actual = (
            tuple(float(value) for value in existing.ks),
            tuple(float(value) for value in existing.phi0s),
            tuple(int(value) for value in existing.periodicitys),
        )
        expected = (
            tuple(float(value) for value in ks),
            tuple(float(value) for value in phases),
            tuple(int(value) for value in periodicities),
        )
        if actual != expected or int(existing.multiple_numbers) != len(ks):
            raise ValidationError("conflicting_metal_assignment_type", name)
        return existing
    created = proper_type(
        name=name,
        k=ks[0],
        phi0=phases[0],
        periodicity=periodicities[0],
    )
    created.ks = list(ks)
    created.phi0s = list(phases)
    created.periodicitys = list(periodicities)
    created.multiple_numbers = len(ks)
    return created


def _snapshot_registries(type_classes: tuple[Any, ...]) -> tuple[tuple[Any, dict[Any, Any]], ...]:
    snapshots = []
    seen = set()
    for type_class in type_classes:
        for attribute in ("_types", "_types_different_name"):
            registry = getattr(type_class, attribute)
            if id(registry) in seen:
                continue
            seen.add(id(registry))
            snapshots.append((registry, dict(registry)))
    return tuple(snapshots)


def _restore_registries(snapshots: tuple[tuple[Any, dict[Any, Any]], ...]) -> None:
    for registry, values in snapshots:
        registry.clear()
        registry.update(values)


def _local_force_terms(result: ParameterizationResult) -> tuple[tuple[str, Mapping[str, Any]], ...]:
    merged: dict[tuple[str, tuple[str, ...]], tuple[str, Mapping[str, Any]]] = {}

    def canonical(kind: str, atom_ids: tuple[str, ...]) -> tuple[str, tuple[str, ...]]:
        if kind in {"bond", "distance_constraint"}:
            return "bond", tuple(sorted(atom_ids))
        if kind in {"angle", "proper_dihedral"}:
            return kind, min(atom_ids, atom_ids[::-1])
        return kind, atom_ids

    layers = (
        result.metal_overlay.bonded_parameters,
        result.bonded_overlay.terms if result.bonded_overlay is not None else {},
    )
    for layer in layers:
        for term_id, value in layer.items():
            atom_ids = tuple(value["atom_ids"])
            merged[canonical(str(value["kind"]), atom_ids)] = (str(term_id), value)
    return tuple(merged[key] for key in sorted(merged))


def _replace_force(
    molecule: Any,
    force_class: Any,
    atoms: list[Any],
    force_type: Any,
    external_id: str,
    kind: str,
) -> None:
    class_name = force_class.get_class_name()
    target = tuple(atoms)

    def equivalent(entity_atoms: tuple[Any, ...]) -> bool:
        if kind in {"bond", "distance_constraint", "angle", "proper_dihedral"}:
            return entity_atoms in {target, target[::-1]}
        if kind == "improper_dihedral":
            return (
                len(entity_atoms) == 4
                and entity_atoms[2] is target[2]
                and frozenset((entity_atoms[0], entity_atoms[1], entity_atoms[3]))
                == frozenset((target[0], target[1], target[3]))
            )
        return entity_atoms == target

    existing = molecule.bonded_forces.setdefault(class_name, [])
    molecule.bonded_forces[class_name] = [
        entity for entity in existing
        if not equivalent(tuple(entity.atoms))
    ]
    molecule.add_bonded_force(force_class.entity(atoms, force_type, external_id))


def _register_native_build_force_types(
    molecule: Any,
    bond_type: Any,
    angle_type: Any,
    proper_type: Any,
    improper_type: Any,
) -> None:
    """Register lookup aliases needed when the native saver rebuilds links."""

    for force in molecule.bonded_forces.get(bond_type.get_class_name(), ()):
        _get_or_create_type(
            bond_type,
            bond_type.Get_Type_Name(force.atoms),
            k=float(force.k),
            b=float(force.b),
        )
    for force in molecule.bonded_forces.get(angle_type.get_class_name(), ()):
        _get_or_create_type(
            angle_type,
            angle_type.Get_Type_Name(force.atoms),
            k=float(force.k),
            b=float(force.b),
        )
    for force in molecule.bonded_forces.get(proper_type.get_class_name(), ()):
        _get_or_create_proper_type(
            proper_type,
            proper_type.Get_Type_Name(force.atoms),
            [float(value) for value in force.ks],
            [float(value) for value in force.phi0s],
            [int(value) for value in force.periodicitys],
        )
    for force in molecule.bonded_forces.get(improper_type.get_class_name(), ()):
        _get_or_create_type(
            improper_type,
            improper_type.Get_Type_Name(force.atoms),
            k=float(force.k),
            phi0=float(force.phi0),
            periodicity=int(force.periodicity),
        )


def prepare_residue_templates(
    patch: Any,
    atom_records: Sequence[Mapping[str, Any]],
) -> Mapping[str, Any]:
    """Ensure embedded patch metals can be loaded without GAFF typing them."""

    from Xponge.forcefield.base.lj_base import LJType
    from Xponge.helper import AtomType, ResidueType
    from .patch import MetalParameterPatch, validate_metal_parameter_patch

    if not isinstance(patch, MetalParameterPatch):
        raise TypeError(
            "metal_assignment.prepare_residue_templates expects a "
            "MetalParameterPatch"
        )
    validate_metal_parameter_patch(patch)
    records = {
        str(record.get("external_id") or ""): record
        for record in atom_records
        if isinstance(record, Mapping)
        and str(record.get("external_id") or "")
    }
    overlay = patch.parameterization_result.metal_overlay
    prepared: list[str] = []
    for external_id in patch.target_metal_atom_ids:
        record = records.get(external_id)
        if record is None:
            raise ValidationError(
                "metal_template_atom_record_missing",
                external_id,
            )
        residue_name = str(record.get("residue_name") or "").strip()
        atom_name = str(record.get("atom_name") or "").strip()
        if not residue_name or not atom_name:
            raise ValidationError(
                "invalid_metal_template_atom_record",
                external_id,
            )
        try:
            residue_type = ResidueType.get_type(residue_name)
        except KeyError as exc:
            raise ValidationError(
                "metal_template_residue_type_missing",
                residue_name,
            ) from exc
        try:
            residue_type.name2atom(atom_name)
            continue
        except KeyError:
            pass

        lj = overlay.lj_parameters.get(external_id)
        mass = overlay.masses.get(external_id)
        charge = overlay.charges.get(external_id)
        if lj is None or mass is None or charge is None:
            raise ValidationError(
                "metal_template_baseline_parameters_missing",
                external_id,
            )
        namespace = patch.patch_hash[:12]
        lj_name = f"MA_LOAD_LJ_{namespace}_{len(prepared)}"
        _get_or_create_type(
            LJType,
            f"{lj_name}-{lj_name}",
            epsilon=float(lj["epsilon"]),
            rmin=float(lj["rmin"]),
        )
        atom_type = _get_or_create_type(
            AtomType,
            f"MA_LOAD_{namespace}_{len(prepared)}",
            charge=float(charge),
            mass=float(mass),
            LJtype=lj_name,
        )
        residue_type.add_atom(atom_name, atom_type, 0.0, 0.0, 0.0)
        prepared.append(external_id)
    return {
        "patch_hash": patch.patch_hash,
        "prepared_atom_ids": tuple(prepared),
    }


def apply(
    molecule: Any,
    patch: Any,
    atom_mapping: Mapping[str, Any] | None = None,
    *,
    inplace: bool = False,
) -> MetalAssignmentResult:
    """Apply only the metal-local overlay to an existing Xponge Molecule.

    The molecule must already carry its ordinary force field.  The function
    never creates, removes, or reorders atoms, residues, or residue links.
    """

    from Xponge.build import build_bonded_force
    from Xponge.forcefield.base.angle_base import AngleType
    from Xponge.forcefield.base.bond_base import BondType
    from Xponge.forcefield.base.dihedral_base import ImproperType, ProperType
    from Xponge.forcefield.base.lj_base import LJType
    from Xponge.helper import AtomType, Molecule
    from .patch import MetalParameterPatch, validate_metal_parameter_patch

    if not isinstance(molecule, Molecule):
        raise TypeError("metal_assignment.apply expects an Xponge Molecule")
    if not isinstance(patch, MetalParameterPatch):
        raise TypeError(
            "metal_assignment.apply expects a MetalParameterPatch"
        )
    validate_metal_parameter_patch(patch)
    result = patch.parameterization_result
    normalized_mapping = _normalized_patch_mapping(
        molecule,
        patch,
        atom_mapping,
    )
    request_id = patch.request_id
    request_atom_ids = {atom.external_id for atom in patch.atoms}
    index_by_external_id = dict(normalized_mapping)
    before = _topology_fingerprint(molecule)
    _validate_patch_links(molecule, patch, index_by_external_id)

    local_terms = _local_force_terms(result)
    local_atom_ids = set(result.metal_overlay.covered_atom_ids)
    local_atom_ids.update(result.metal_overlay.charges)
    if result.charge_overlay is not None:
        local_atom_ids.update(result.charge_overlay.charges)
    for _, term in local_terms:
        local_atom_ids.update(term["atom_ids"])
    if not local_atom_ids <= request_atom_ids:
        raise ValidationError(
            "metal_overlay_atom_not_in_request",
            ",".join(sorted(local_atom_ids - request_atom_ids)),
        )

    namespace = result.result_hash[:12]
    registry_snapshots = _snapshot_registries(
        (AtomType, LJType, BondType, AngleType, ProperType, ImproperType)
    )
    molecule_registry = dict(Molecule._all)
    try:
        # Always work on a copy.  ``inplace`` is a commit policy, not permission
        # to expose a half-applied overlay when validation fails.
        assigned = molecule.deepcopy()
        # Molecule.deepcopy() preserves residue topology but reconstructs each
        # Residue from its type name.  Instance-level names carry source
        # identity and must survive this transactional boundary.
        for source_residue, assigned_residue in zip(molecule.residues, assigned.residues):
            assigned_residue.name = source_residue.name
        assigned.get_atoms()
        if not assigned.built:
            suspended_patch_links = []
            for link in patch.required_links:
                atom1 = assigned.atoms[index_by_external_id[link.atom_ids[0]]]
                atom2 = assigned.atoms[index_by_external_id[link.atom_ids[1]]]
                residue_link = assigned.get_residue_link(atom1, atom2)
                if residue_link is not None:
                    assigned.del_residue_link(atom1, atom2)
                    suspended_patch_links.append((atom1, atom2))
            try:
                build_bonded_force(assigned)
            except Exception as exc:
                raise ValidationError(
                    "ordinary_force_field_not_assigned",
                    "build the ordinary Xponge force field before applying metal parameters",
                ) from exc
            for atom1, atom2 in suspended_patch_links:
                assigned.add_residue_link(atom1, atom2)

        atom_by_external_id = {
            external_id: assigned.atoms[index]
            for external_id, index in normalized_mapping
        }
        lj_names: dict[tuple[float, float], str] = {}
        for local_index, external_id in enumerate(sorted(local_atom_ids)):
            atom = atom_by_external_id[external_id]
            old_type = atom.type
            original_atom_name = atom.name
            original_coordinates = (atom.x, atom.y, atom.z)
            contents = dict(old_type.contents)
            if getattr(atom, "charge", None) is not None:
                contents["charge"] = float(atom.charge)
            if external_id in result.metal_overlay.charges:
                contents["charge"] = float(result.metal_overlay.charges[external_id])
            if result.charge_overlay is not None and external_id in result.charge_overlay.charges:
                contents["charge"] = float(result.charge_overlay.charges[external_id])
            if external_id in result.metal_overlay.masses:
                contents["mass"] = float(result.metal_overlay.masses[external_id])
            if external_id in result.metal_overlay.lj_parameters:
                values = result.metal_overlay.lj_parameters[external_id]
                signature = (float(values["epsilon"]), float(values["rmin"]))
                lj_name = lj_names.get(signature)
                if lj_name is None:
                    lj_name = f"MA_LJ_{namespace}_{len(lj_names)}"
                    _get_or_create_type(
                        LJType,
                        f"{lj_name}-{lj_name}",
                        epsilon=signature[0],
                        rmin=signature[1],
                    )
                    lj_names[signature] = lj_name
                contents["LJtype"] = lj_name
            requested_type = result.metal_overlay.atom_types.get(external_id, old_type.name)
            parameters = {
                key: value
                for key, value in contents.items()
                if key in AtomType._parameters and key != "name" and value is not None
            }
            effective_name = (
                f"MA_{namespace}_{local_index}_{requested_type}_{_short_hash(parameters)}"
            )
            effective_type = _get_or_create_type(AtomType, effective_name, **parameters)
            atom.type = effective_type
            atom.contents = dict(effective_type.contents)
            atom.name = original_atom_name
            atom.x, atom.y, atom.z = original_coordinates

        force_counts: dict[str, int] = {}
        for term_index, (term_id, value) in enumerate(local_terms):
            kind = str(value["kind"])
            atoms = [atom_by_external_id[atom_id] for atom_id in value["atom_ids"]]
            type_name = f"MA_TERM_{namespace}_{term_index}_{_short_hash(value)}"
            term = type("_Term", (), {
                "parameters": value["parameters"],
                "external_id": term_id,
            })()
            if kind in {"bond", "distance_constraint"}:
                k, equilibrium = _bond_values(term)
                force_type = _get_or_create_type(BondType, type_name, k=k, b=equilibrium)
                _replace_force(assigned, BondType, atoms, force_type, term_id, kind)
            elif kind == "angle":
                k, equilibrium = _angle_values(term)
                force_type = _get_or_create_type(AngleType, type_name, k=k, b=equilibrium)
                _replace_force(assigned, AngleType, atoms, force_type, term_id, kind)
            elif kind == "proper_dihedral":
                ks, phases, periodicities = _proper_values(term)
                force_type = _get_or_create_proper_type(
                    ProperType,
                    type_name,
                    ks,
                    phases,
                    periodicities,
                )
                _replace_force(assigned, ProperType, atoms, force_type, term_id, kind)
            elif kind == "improper_dihedral":
                k, phase, periodicity = _improper_values(term)
                force_type = _get_or_create_type(
                    ImproperType,
                    f"{type_name}_IMPROPER",
                    k=k,
                    phi0=phase,
                    periodicity=periodicity,
                )
                _replace_force(assigned, ImproperType, atoms, force_type, term_id, kind)
            else:
                raise ValidationError("unsupported_assigned_force_kind", kind, term_id)
            force_counts[kind] = force_counts.get(kind, 0) + 1

        _register_native_build_force_types(
            assigned,
            BondType,
            AngleType,
            ProperType,
            ImproperType,
        )
        after = _topology_fingerprint(assigned)
        if after != before:
            raise ValidationError("molecule_topology_changed_during_apply", request_id)
        assigned.built = True
    except Exception:
        _restore_registries(registry_snapshots)
        Molecule._all.clear()
        Molecule._all.update(molecule_registry)
        raise

    if inplace:
        molecule.__dict__.clear()
        molecule.__dict__.update(assigned.__dict__)
        assigned = molecule
        Molecule._all[assigned.name] = assigned
    report = {
        "request_id": request_id,
        "result_hash": result.result_hash,
        "patch_hash": patch.patch_hash,
        "inplace": bool(inplace),
        "atom_count": len(assigned.atoms),
        "residue_count": len(assigned.residues),
        "residue_link_count": len(assigned.residue_links),
        "local_atom_count": len(local_atom_ids),
        "force_counts": force_counts,
        "topology_preserved": True,
        "ordinary_force_field_preserved": True,
    }
    return MetalAssignmentResult(
        molecule=assigned,
        parameterization_result=result,
        atom_mapping=normalized_mapping,
        application_report=report,
    )


__all__ = [
    "MetalAssignmentResult",
    "apply",
    "parameterize",
    "prepare_residue_templates",
]
