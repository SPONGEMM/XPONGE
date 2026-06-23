"""Validation and environment selection helpers for MCPB."""

from __future__ import annotations

from collections.abc import Sequence

from ._compat import atom_count, atom_to_residue_map, refresh_molecule_indices, residue_count
from .models import MCPBIonInfo, MCPBRequest, MCPBSelection

PERIODIC_TABLE_SYMBOLS = (
    "H", "He",
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
    "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)
PERIODIC_TABLE_SET = set(PERIODIC_TABLE_SYMBOLS)


def normalize_element_symbol(element: str) -> str:
    text = str(element).strip()
    if not text:
        raise ValueError("element symbol must not be empty")
    if len(text) == 1:
        return text.upper()
    return text[0].upper() + text[1:].lower()


def infer_element_symbol(*candidates: str) -> str:
    for candidate in candidates:
        text = "".join(ch for ch in str(candidate).strip() if ch.isalpha())
        if not text:
            continue
        direct = normalize_element_symbol(text)
        if direct in PERIODIC_TABLE_SET:
            return direct
        if len(text) >= 2:
            two = normalize_element_symbol(text[:2])
            if two in PERIODIC_TABLE_SET:
                return two
        one = normalize_element_symbol(text[:1])
        if one in PERIODIC_TABLE_SET:
            return one
    raise ValueError(f"unsupported element symbol: {candidates[0] if candidates else ''}")


def _coerce_atom_id(value: int, atom_total: int, label: str) -> int:
    atom_id = int(value)
    if atom_id < 0 or atom_id >= atom_total:
        raise IndexError(f"{label} atom id out of range: {atom_id}")
    return atom_id


def _coerce_residue_id(value: int, residue_total: int) -> int:
    residue_id = int(value)
    if residue_id < 0 or residue_id >= residue_total:
        raise IndexError(f"additional residue id out of range: {residue_id}")
    return residue_id


def _normalize_ion_infos(request: MCPBRequest) -> tuple[MCPBIonInfo, ...]:
    provided = {int(info.atom_id): info for info in request.ion_info}
    refresh_molecule_indices(request.molecule)
    atom_to_residue = atom_to_residue_map(request.molecule)
    normalized: list[MCPBIonInfo] = []
    for atom_id in request.ion_ids:
        atom = request.molecule.atoms[atom_id]
        residue = request.molecule.residues[atom_to_residue[atom_id]]
        info = provided.get(atom_id)
        if info is None:
            element = infer_element_symbol(getattr(atom, "element", ""), atom.name, residue.name)
            info = MCPBIonInfo(
                atom_id=atom_id,
                element=element,
                formal_charge=None,
                spin=None,
                resname=residue.name,
                atom_name=atom.name,
            )
        else:
            element = normalize_element_symbol(info.element)
            info = MCPBIonInfo(
                atom_id=atom_id,
                element=element,
                formal_charge=info.formal_charge,
                spin=info.spin,
                resname=info.resname or residue.name,
                atom_name=info.atom_name or atom.name,
                metadata=dict(info.metadata),
            )
        if info.element not in PERIODIC_TABLE_SET:
            raise ValueError(f"unsupported element symbol: {info.element}")
        normalized.append(info)
    extras = sorted(set(provided) - set(request.ion_ids))
    if extras:
        raise ValueError(f"ion_info contains atom ids not listed in ion_ids: {extras}")
    return tuple(normalized)


def normalize_request(
    molecule,
    ion_ids: Sequence[int],
    ion_info: Sequence[MCPBIonInfo | dict] | None,
    *,
    method: str = "seminario",
    model: str = "bonded",
    cutoff: float = 2.8,
    bonded_pairs: Sequence[Sequence[int]] | None = None,
    additional_residue_ids: Sequence[int] | None = None,
    charge_mode: str = "resp",
    qm_backend: str = "pyscf",
    basis: str | None = None,
    scale_factor: float = 1.0,
    frcmod_files: Sequence[str] | None = None,
    gaff: int | None = None,
    force_field: str | None = None,
    water_model: str | None = None,
) -> MCPBRequest:
    if molecule is None or not hasattr(molecule, "residues"):
        raise TypeError("Xponge.MCPB expects a Molecule-like object as the first argument")
    if not ion_ids:
        raise ValueError("ion_ids must not be empty")
    if model not in {"bonded", "nonbonded"}:
        raise ValueError(f"unsupported MCPB model: {model}")
    if charge_mode not in {"resp", "fixed"}:
        raise ValueError(f"unsupported MCPB charge_mode: {charge_mode}")
    if cutoff <= 0:
        raise ValueError("cutoff must be positive")
    atom_total = atom_count(molecule)
    residue_total = residue_count(molecule)
    normalized_ion_ids = tuple(dict.fromkeys(_coerce_atom_id(atom_id, atom_total, "ion_ids") for atom_id in ion_ids))
    normalized_bonded_pairs: list[tuple[int, int]] = []
    if bonded_pairs is not None:
        for pair in bonded_pairs:
            if len(pair) != 2:
                raise ValueError(f"bonded pair must have length 2: {pair}")
            atom1 = _coerce_atom_id(pair[0], atom_total, "bonded_pairs")
            atom2 = _coerce_atom_id(pair[1], atom_total, "bonded_pairs")
            if atom1 == atom2:
                raise ValueError(f"bonded pair must connect two distinct atoms: {pair}")
            normalized = (atom1, atom2) if atom1 < atom2 else (atom2, atom1)
            if normalized not in normalized_bonded_pairs:
                normalized_bonded_pairs.append(normalized)
    if model == "bonded" and not normalized_bonded_pairs:
        raise ValueError("bonded MCPB mode currently requires explicit bonded_pairs")
    normalized_additional_residue_ids = tuple(
        dict.fromkeys(_coerce_residue_id(residue_id, residue_total) for residue_id in (additional_residue_ids or ()))
    )
    normalized_infos: list[MCPBIonInfo] = []
    for entry in ion_info or ():
        if isinstance(entry, MCPBIonInfo):
            normalized_infos.append(entry)
        elif isinstance(entry, dict):
            normalized_infos.append(MCPBIonInfo(**entry))
        else:
            raise TypeError("ion_info entries must be MCPBIonInfo instances or dicts")
    return MCPBRequest(
        molecule=molecule,
        ion_ids=normalized_ion_ids,
        ion_info=tuple(normalized_infos),
        method=str(method),
        model=str(model),
        cutoff=float(cutoff),
        bonded_pairs=tuple(normalized_bonded_pairs),
        additional_residue_ids=normalized_additional_residue_ids,
        charge_mode=str(charge_mode),
        qm_backend=str(qm_backend),
        basis=None if basis is None else str(basis),
        scale_factor=float(scale_factor),
        frcmod_files=tuple(str(path) for path in (frcmod_files or ())),
        gaff=None if gaff is None else int(gaff),
        force_field=None if force_field is None else str(force_field),
        water_model=None if water_model is None else str(water_model),
    )


def validate_and_select_environment(request: MCPBRequest) -> tuple[MCPBRequest, MCPBSelection, list[str]]:
    ion_info = _normalize_ion_infos(request)
    warnings: list[str] = []
    ion_ids = set(request.ion_ids)
    refresh_molecule_indices(request.molecule)
    atom_to_residue = atom_to_residue_map(request.molecule)
    coordinating_atom_ids: set[int] = set()
    selected_residue_ids: set[int] = set(request.additional_residue_ids)
    for info in ion_info:
        atom = request.molecule.atoms[info.atom_id]
        residue = request.molecule.residues[atom_to_residue[info.atom_id]]
        imported = infer_element_symbol(getattr(atom, "element", ""), atom.name, residue.name)
        if imported != info.element:
            raise ValueError(
                f"ion_info element {info.element} does not match imported atom element {imported} for atom {info.atom_id}"
            )
        selected_residue_ids.add(atom_to_residue[info.atom_id])
        if info.formal_charge is None:
            warnings.append(f"atom {info.atom_id} has no explicit formal_charge in ion_info")
        if info.spin is None:
            warnings.append(f"atom {info.atom_id} has no explicit spin in ion_info")
    for atom1, atom2 in request.bonded_pairs:
        if atom1 in ion_ids and atom2 in ion_ids:
            raise ValueError(f"bonded pair connects two ion atoms; expected metal-ligand pair: {(atom1, atom2)}")
        if atom1 not in ion_ids and atom2 not in ion_ids:
            raise ValueError(f"bonded pair must include one selected ion atom: {(atom1, atom2)}")
        neighbor = atom2 if atom1 in ion_ids else atom1
        coordinating_atom_ids.add(neighbor)
        selected_residue_ids.add(atom_to_residue[neighbor])
    normalized_request = MCPBRequest(
        molecule=request.molecule,
        ion_ids=request.ion_ids,
        ion_info=ion_info,
        method=request.method,
        model=request.model,
        cutoff=request.cutoff,
        bonded_pairs=request.bonded_pairs,
        additional_residue_ids=request.additional_residue_ids,
        charge_mode=request.charge_mode,
        qm_backend=request.qm_backend,
        basis=request.basis,
        scale_factor=request.scale_factor,
        frcmod_files=request.frcmod_files,
        gaff=request.gaff,
        force_field=request.force_field,
        water_model=request.water_model,
    )
    selection = MCPBSelection(
        ion_atom_ids=request.ion_ids,
        coordinating_atom_ids=tuple(sorted(coordinating_atom_ids)),
        selected_residue_ids=tuple(sorted(selected_residue_ids)),
        bonded_pairs=request.bonded_pairs,
    )
    return normalized_request, selection, warnings
