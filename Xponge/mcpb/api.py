"""Initial MCPB workflow entrypoint."""

from __future__ import annotations

from ..helper import ResidueType, set_real_global_variable
from ._compat import add_residue_link_by_atom_ids, atom_to_residue_map, get_residue_link_by_atom_ids, refresh_molecule_indices
from .charge_refit import run_local_charge_refit
from .export import audit_sponge_ready
from .frcmod import (
    register_blank_parameters,
    register_empirical_parameters,
    write_blank_frcmod_artifact,
    write_empirical_frcmod_artifact,
)
from .model_builder import build_small_and_large_models
from .models import MCPBResult
from .seminario import register_seminario_parameters, write_seminario_frcmod_artifact
from .selection import normalize_request, validate_and_select_environment

_ION_FORCEFIELD_LOADED: set[str] = set()


def _canonical_single_atom_ion_name(element: str, formal_charge: int | None) -> str | None:
    if formal_charge is None:
        return None
    charge = int(formal_charge)
    symbol = str(element).strip().upper()
    if charge == 1:
        return symbol
    if charge > 1:
        return f"{symbol}{charge}"
    return None


def _ensure_amber_ion_forcefield_loaded(water_model: str | None) -> str:
    normalized = "tip3p" if water_model is None else str(water_model).strip().lower()
    if normalized not in {"tip3p", "spce", "tip4pew", "opc"}:
        normalized = "tip3p"
    if normalized in _ION_FORCEFIELD_LOADED:
        return normalized
    if normalized == "tip3p":
        from ..forcefield.amber import tip3p  # noqa: F401
    elif normalized == "spce":
        from ..forcefield.amber import spce  # noqa: F401
    elif normalized == "tip4pew":
        from ..forcefield.amber import tip4pew  # noqa: F401
    else:
        from ..forcefield.amber import opc  # noqa: F401
    _ION_FORCEFIELD_LOADED.add(normalized)
    return normalized


def _single_atom_template_mol2_text(template_name: str, atom_name: str, atom_type: str, charge: float) -> str:
    return (
        "@<TRIPOS>MOLECULE\n"
        f"{template_name}\n"
        "1 0 1\n"
        "SMALL\n"
        "USER_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        f"1 {atom_name} 0.000000 0.000000 0.000000 {atom_type} 1 {template_name} {charge:.6f}\n"
        "@<TRIPOS>BOND\n"
    )


def _ensure_single_atom_metal_templates(request) -> list[str]:
    from .. import get_template_molecule, has_template, register_residue_templates_from_mol2_text

    refresh_molecule_indices(request.molecule)
    source_atom_to_residue = atom_to_residue_map(request.molecule)
    registered: list[str] = []
    _ensure_amber_ion_forcefield_loaded(request.water_model)
    for info in request.ion_info:
        atom = request.molecule.atoms[info.atom_id]
        residue = request.molecule.residues[source_atom_to_residue[info.atom_id]]
        if len(residue.atoms) != 1:
            continue
        template_name = residue.name
        canonical = _canonical_single_atom_ion_name(info.element, info.formal_charge)
        canonical_atom_type = None
        if canonical and has_template(canonical):
            template = get_template_molecule(canonical)
            refresh_molecule_indices(template)
            template_atom = template.atoms[0]
            atom.type = template_atom.type
            atom.mass = template_atom.mass
            canonical_atom_type = getattr(template_atom.type, "name", str(template_atom.type))
            if template_name != canonical:
                ResidueType.set_type(template_name, ResidueType.get_type(canonical))
                set_real_global_variable(template_name, ResidueType.get_type(canonical))
        elif canonical:
            raise ValueError(
                f"No packaged Amber ion template is available for {info.element} with formal_charge "
                f"{info.formal_charge}. Register a residue template or ion parameters explicitly before MCPB."
            )

        if has_template(template_name):
            template = get_template_molecule(template_name)
            refresh_molecule_indices(template)
            template_atom = template.atoms[0]
            if atom.mass <= 0.0:
                atom.mass = template_atom.mass
            if canonical_atom_type is None:
                canonical_atom_type = getattr(template_atom.type, "name", str(template_atom.type))
            continue

        charge = float(info.formal_charge if info.formal_charge is not None else atom.charge)
        atom_type = canonical_atom_type or getattr(atom.type, "name", str(atom.type)) or atom.element or info.element
        register_residue_templates_from_mol2_text(
            _single_atom_template_mol2_text(template_name, atom.name, atom_type, charge)
        )
        template = get_template_molecule(template_name)
        refresh_molecule_indices(template)
        atom.mass = template.atoms[0].mass
        registered.append(template_name)
    return registered


def _patch_parent_molecule_connectivity(request) -> list[tuple[int, int]]:
    connect_records: list[tuple[int, int]] = []
    for atom1, atom2 in request.bonded_pairs:
        record = (atom1, atom2) if atom1 < atom2 else (atom2, atom1)
        if get_residue_link_by_atom_ids(request.molecule, record[0], record[1]) is None:
            add_residue_link_by_atom_ids(request.molecule, record[0], record[1])
        connect_records.append(record)
    return connect_records


def MCPB(
    molecule,
    ion_ids,
    ion_info=None,
    *,
    method="seminario",
    model="bonded",
    cutoff=2.8,
    bonded_pairs=None,
    additional_residue_ids=None,
    charge_mode="resp",
    qm_backend="pyscf",
    basis=None,
    scale_factor=1.0,
    frcmod_files=None,
    gaff=None,
    force_field=None,
    water_model=None,
):
    """Initial MCPB workflow scaffold for Xponge-origin."""

    request = normalize_request(
        molecule,
        ion_ids,
        ion_info,
        method=method,
        model=model,
        cutoff=cutoff,
        bonded_pairs=bonded_pairs,
        additional_residue_ids=additional_residue_ids,
        charge_mode=charge_mode,
        qm_backend=qm_backend,
        basis=basis,
        scale_factor=scale_factor,
        frcmod_files=frcmod_files,
        gaff=gaff,
        force_field=force_field,
        water_model=water_model,
    )
    request, selection, warnings = validate_and_select_environment(request)
    registered_metal_templates = _ensure_single_atom_metal_templates(request)
    small_model, large_model = build_small_and_large_models(request, selection)
    connect_records = _patch_parent_molecule_connectivity(request)
    updated_charge_atoms: list[int] = []
    charge_refit_summary = None
    frcmod_path = None
    pending_requirements = []
    if request.charge_mode == "resp":
        if request.basis is None:
            pending_requirements.append("local RESP charge refit pending: provide basis to execute MCPB refit")
        else:
            charge_refit_summary = run_local_charge_refit(request, large_model)
            updated_charge_atoms = list(charge_refit_summary["updated_atom_ids"])
    else:
        pending_requirements.append("fixed-charge MCPB mode selected: local RESP refit skipped")
    if request.method == "blank":
        register_blank_parameters(request, selection)
        frcmod_path = write_blank_frcmod_artifact(request, selection)
    elif request.method == "empirical":
        register_empirical_parameters(request, selection)
        frcmod_path = write_empirical_frcmod_artifact(request, selection)
    elif request.method == "seminario":
        if request.basis is None:
            pending_requirements.append("seminario MCPB frcmod generation pending: provide basis to execute Hessian path")
        else:
            seminario_summary = register_seminario_parameters(request, selection, small_model)
            frcmod_path = write_seminario_frcmod_artifact(seminario_summary=seminario_summary)
    else:
        pending_requirements.append(f"MCPB method '{request.method}' force-constant generation not implemented yet")
    sponge_ready, sponge_error = audit_sponge_ready(request.molecule)
    if sponge_error is not None:
        pending_requirements.append(f"Save_SPONGE_Input audit failed: {sponge_error}")
    return MCPBResult(
        molecule=request.molecule,
        request=request,
        selection=selection,
        small_model=small_model,
        large_model=large_model,
        frcmod_path=frcmod_path,
        connect_records=connect_records,
        updated_charge_atoms=updated_charge_atoms,
        registered_metal_templates=registered_metal_templates,
        sponge_ready=sponge_ready,
        pending_requirements=pending_requirements,
        warnings=warnings,
        metadata={
            "method": request.method,
            "qm_backend": request.qm_backend,
            "charge_mode": request.charge_mode,
            "model": request.model,
            "small_model_residue_ids": tuple(sorted(small_model.residue_id_map)),
            "large_model_residue_ids": tuple(sorted(large_model.residue_id_map)),
            "charge_refit": charge_refit_summary,
            "frcmod_generated": frcmod_path is not None,
            "sponge_audit_error": sponge_error,
        },
    )
