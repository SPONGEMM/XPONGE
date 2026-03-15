"""
Helpers to align alternating GAFF/GAFF2 atom-type orientation with AmberTools.
"""
from ...helper import AtomType


def _iter_unique_bonds(assign):
    for atom_i, bonded in assign.bonds.items():
        for atom_j, bond_order in bonded.items():
            if atom_i < atom_j:
                yield atom_i, atom_j, bond_order


def _iter_sequence_bonds(assign):
    seen = set()
    for atom_i, atom_j in getattr(assign, "_bond_sequence", []):
        if atom_j not in assign.bonds.get(atom_i, {}):
            continue
        seen.add((atom_i, atom_j))
        yield atom_i, atom_j, assign.bonds[atom_i][atom_j]
    for atom_i, bonded in assign.bonds.items():
        for atom_j, bond_order in bonded.items():
            if atom_i < atom_j and (atom_i, atom_j) not in seen:
                yield atom_i, atom_j, bond_order


def _normalize_to_primary(assign, primary_to_secondary):
    secondary_to_primary = {value: key for key, value in primary_to_secondary.items()}
    for atom_i, atom_type in assign.atom_types.items():
        primary = secondary_to_primary.get(atom_type.name)
        if primary is not None:
            assign.atom_types[atom_i] = AtomType.get_type(primary)


def _orient_and_convert(assign, primary_to_secondary):
    _normalize_to_primary(assign, primary_to_secondary)
    candidate_names = set(primary_to_secondary.keys())
    bond_iterator = _iter_sequence_bonds
    if "ce" not in candidate_names or not any(atom_type.name == "ce" for atom_type in assign.atom_types.values()):
        bond_iterator = _iter_unique_bonds
    atom_sign = [0 for _ in range(assign.atom_numbers)]
    is_candidate = [False for _ in range(assign.atom_numbers)]
    has_candidate = False
    for atom_i, atom_type in assign.atom_types.items():
        if atom_type.name in candidate_names:
            is_candidate[atom_i] = True
            if not has_candidate:
                atom_sign[atom_i] = 1
                has_candidate = True
    if not has_candidate:
        return

    # Propagate orientation signs along bonds in the same spirit as AmberTools atomtype.c.
    changed = True
    while changed:
        changed = False
        seeded = False
        for atom_i, atom_j, bond_order in bond_iterator(assign):
            if not (is_candidate[atom_i] and is_candidate[atom_j]):
                continue
            if not seeded and atom_sign[atom_i] == 0 and atom_sign[atom_j] == 0:
                atom_sign[atom_i] = 1
                changed = True
                seeded = True
            if atom_sign[atom_i] == 0 and atom_sign[atom_j] != 0:
                atom_sign[atom_i] = atom_sign[atom_j] if bond_order == 1 else -atom_sign[atom_j]
                changed = True
            if atom_sign[atom_j] == 0 and atom_sign[atom_i] != 0:
                atom_sign[atom_j] = atom_sign[atom_i] if bond_order == 1 else -atom_sign[atom_i]
                changed = True

    for atom_i, atom_type in assign.atom_types.items():
        if atom_sign[atom_i] != -1:
            continue
        secondary = primary_to_secondary.get(atom_type.name)
        if secondary is not None:
            assign.atom_types[atom_i] = AtomType.get_type(secondary)


def apply_amber_alternating_type_adjustment(assign):
    """
    AmberTools-compatible orientation pass for alternating type families.
    """
    _demote_specific_nc_nd_to_n2(assign)
    _orient_and_convert(
        assign,
        {
            "cc": "cd",
            "ce": "cf",
            "cg": "ch",
            "pc": "pd",
            "pe": "pf",
            "nc": "nd",
            "ne": "nf",
        },
    )
    _orient_and_convert(assign, {"cp": "cq"})


def _demote_specific_nc_nd_to_n2(assign):
    """
    Narrow correction for AR3 five-membered non-fused nitrogens where Amber keeps n2.
    """
    for atom_i, atom_type in assign.atom_types.items():
        if atom_type.name not in {"nc", "nd"}:
            continue
        marker = assign.atom_marker[atom_i]
        if "AR3" not in marker:
            continue
        if marker.get("RG5") != 1 or "RG9" in marker:
            continue
        has_single_n3_neighbor = any(
            bond_order == 1 and assign.Atom_Judge(bonded_atom, "N3")
            for bonded_atom, bond_order in assign.bonds[atom_i].items()
        )
        if has_single_n3_neighbor:
            assign.atom_types[atom_i] = AtomType.get_type("n2")
