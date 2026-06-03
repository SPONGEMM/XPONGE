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


def _orient_and_convert(assign, primary_to_secondary, same_bond_orders, opposite_bond_orders):
    _normalize_to_primary(assign, primary_to_secondary)
    candidate_names = set(primary_to_secondary.keys())
    atom_sign = [0 for _ in range(assign.atom_numbers)]
    is_candidate = [False for _ in range(assign.atom_numbers)]
    has_candidate = False
    num = 0
    for atom_i, atom_type in assign.atom_types.items():
        if atom_type.name in candidate_names:
            is_candidate[atom_i] = True
            if not has_candidate:
                atom_sign[atom_i] = 1
                has_candidate = True
            num += 1
    if not has_candidate:
        return
    num -= 1

    # Mirror AmberTools atomtype.c: orient one combined alternating-family graph.
    while num > 0:
        num -= 1
        flag = 0
        for atom_i, atom_j, bond_order in _iter_sequence_bonds(assign):
            if not (is_candidate[atom_i] and is_candidate[atom_j]):
                continue
            if flag == 0 and atom_sign[atom_i] == 0 and atom_sign[atom_j] == 0:
                atom_sign[atom_i] = 1
            if atom_sign[atom_i] == 0 and atom_sign[atom_j] != 0:
                flag = 1
                if bond_order in same_bond_orders:
                    atom_sign[atom_i] = atom_sign[atom_j]
                elif bond_order in opposite_bond_orders:
                    atom_sign[atom_i] = -atom_sign[atom_j]
            if atom_sign[atom_j] == 0 and atom_sign[atom_i] != 0:
                flag = 1
                if bond_order in same_bond_orders:
                    atom_sign[atom_j] = atom_sign[atom_i]
                elif bond_order in opposite_bond_orders:
                    atom_sign[atom_j] = -atom_sign[atom_i]

    for atom_i, atom_type in assign.atom_types.items():
        if atom_sign[atom_i] != -1:
            continue
        secondary = primary_to_secondary.get(atom_type.name)
        if secondary is not None:
            assign.atom_types[atom_i] = AtomType.get_type(secondary)


def _adjust_cp_cq(assign):
    _normalize_to_primary(assign, {"cp": "cq"})
    atom_sign = [0 for _ in range(assign.atom_numbers)]
    is_candidate = [False for _ in range(assign.atom_numbers)]
    index = 0
    num = 0
    for atom_i, atom_type in assign.atom_types.items():
        if atom_type.name == "cp":
            is_candidate[atom_i] = True
            if index == 0:
                atom_sign[atom_i] = 1
                index = 1
            num += 1
    if num == 0:
        return
    num -= 1

    while num > 0:
        num -= 1
        for atom_i, atom_j, bond_order in _iter_sequence_bonds(assign):
            if not (is_candidate[atom_i] and is_candidate[atom_j]):
                continue
            same_phase = bond_order == 1 and "AB" not in assign.bond_marker[atom_i][atom_j]
            if atom_sign[atom_i] == 0 and atom_sign[atom_j] != 0:
                atom_sign[atom_i] = atom_sign[atom_j] if same_phase else -atom_sign[atom_j]
            if atom_sign[atom_j] == 0 and atom_sign[atom_i] != 0:
                atom_sign[atom_j] = atom_sign[atom_i] if same_phase else -atom_sign[atom_i]

    for atom_i, atom_type in assign.atom_types.items():
        if atom_sign[atom_i] == -1 and atom_type.name == "cp":
            assign.atom_types[atom_i] = AtomType.get_type("cq")


def apply_amber_alternating_type_adjustment(assign):
    """
    AmberTools-compatible orientation pass for alternating type families.
    """
    _demote_specific_nc_nd_to_n2(assign)
    _demote_sequence_sensitive_ne_nf_to_n2(assign)
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
        same_bond_orders={1, 7},
        opposite_bond_orders={2, 3, 8},
    )
    _adjust_cp_cq(assign)


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


def _demote_sequence_sensitive_ne_nf_to_n2(assign):
    """
    Mirror AmberTools bondtype's direction-sensitive neutral nitroso fallback.

    When a terminal oxo bond is recorded as O->N in the input bond list, Amber's
    bondtype stage can mark that bond as type 6 before atom typing, which keeps
    the nitrogen out of the ne/nf family and leaves it as n2. Match the final
    AmberTools result here without perturbing the broader ne/nf logic.
    """
    sequence_bonds = set(getattr(assign, "_bond_sequence", []))
    for atom_i, atom_type in assign.atom_types.items():
        if atom_type.name not in {"ne", "nf"}:
            continue
        terminal_double_heteros = []
        has_other_terminal_hetero = False
        for bonded_atom, bond_order in assign.bonds[atom_i].items():
            if assign.atoms[bonded_atom] not in {"O", "S"} or len(assign.bonds[bonded_atom]) != 1:
                continue
            if bond_order == 2:
                terminal_double_heteros.append(bonded_atom)
            else:
                has_other_terminal_hetero = True
        if len(terminal_double_heteros) != 1 or has_other_terminal_hetero:
            continue
        hetero = terminal_double_heteros[0]
        if (hetero, atom_i) in sequence_bonds:
            assign.atom_types[atom_i] = AtomType.get_type("n2")
