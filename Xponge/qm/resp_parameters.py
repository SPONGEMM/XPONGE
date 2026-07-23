"""Default RESP electrostatic-potential parameters."""

from __future__ import annotations

from dataclasses import dataclass


DEFAULT_MK_RADIUS = 1.7
BASIS_ID_TO_FAMILY = {
    1: "6-31G*",
    2: "CEP-31G",
    3: "SDD",
}


@dataclass(frozen=True, slots=True)
class RespParameter:
    atomic_number: int
    symbol: str
    vdw_radius: float
    mk_radius: float
    applied: int
    basis_set: int


_PARAMETER_ROWS = (
    (1, "H", 1.20, 1.14, 0, 1),
    (2, "He", 1.40, 1.33, 1, 1),
    (3, "Li", 1.82, 1.73, 1, 1),
    (4, "Be", 0.00, 0.00, 1, 1),
    (5, "B", 0.00, 0.00, 1, 1),
    (6, "C", 1.70, 1.61, 0, 1),
    (7, "N", 1.55, 1.47, 0, 1),
    (8, "O", 1.52, 1.44, 0, 1),
    (9, "F", 1.47, 1.40, 0, 1),
    (10, "Ne", 1.54, 1.46, 1, 1),
    (11, "Na", 2.27, 2.16, 1, 1),
    (12, "Mg", 1.73, 1.64, 1, 1),
    (13, "Al", 0.00, 0.00, 1, 1),
    (14, "Si", 2.10, 1.99, 1, 1),
    (15, "P", 1.80, 1.71, 0, 1),
    (16, "S", 1.80, 1.71, 0, 1),
    (17, "Cl", 1.75, 1.66, 0, 1),
    (18, "Ar", 1.88, 1.79, 1, 1),
    (19, "K", 2.75, 2.61, 1, 1),
    (20, "Ca", 0.00, 0.00, 1, 1),
    (21, "Sc", 0.00, 0.00, 1, 1),
    (22, "Ti", 0.00, 0.00, 1, 1),
    (23, "V", 0.00, 0.00, 1, 1),
    (24, "Cr", 0.00, 0.00, 1, 1),
    (25, "Mn", 0.00, 0.00, 1, 1),
    (26, "Fe", 0.00, 0.00, 1, 1),
    (27, "Co", 0.00, 0.00, 1, 1),
    (28, "Ni", 1.63, 1.55, 1, 1),
    (29, "Cu", 1.40, 1.33, 1, 1),
    (30, "Zn", 1.39, 1.32, 0, 1),
    (31, "Ga", 1.87, 1.78, 1, 1),
    (32, "Ge", 0.00, 0.00, 1, 1),
    (33, "As", 1.85, 1.76, 1, 1),
    (34, "Se", 1.90, 1.80, 1, 1),
    (35, "Br", 1.85, 1.76, 0, 1),
    (36, "Kr", 2.02, 1.92, 1, 1),
    (37, "Rb", 0.00, 0.00, 1, 2),
    (38, "Sr", 0.00, 0.00, 1, 2),
    (39, "Y", 0.00, 0.00, 1, 2),
    (40, "Zr", 0.00, 0.00, 1, 2),
    (41, "Nb", 0.00, 0.00, 1, 2),
    (42, "Mo", 0.00, 0.00, 1, 2),
    (43, "Tc", 0.00, 0.00, 1, 2),
    (44, "Ru", 0.00, 0.00, 1, 2),
    (45, "Rh", 0.00, 0.00, 1, 2),
    (46, "Pd", 1.63, 1.55, 1, 2),
    (47, "Ag", 1.72, 1.63, 1, 2),
    (48, "Cd", 1.58, 1.50, 1, 2),
    (49, "In", 1.93, 1.83, 1, 2),
    (50, "Sn", 2.17, 2.06, 1, 2),
    (51, "Sb", 0.00, 0.00, 1, 2),
    (52, "Te", 2.06, 1.96, 1, 2),
    (53, "I", 1.98, 1.99, 1, 2),
    (54, "Xe", 2.16, 2.05, 1, 2),
    (55, "Cs", 0.00, 0.00, 1, 2),
    (56, "Ba", 0.00, 0.00, 1, 2),
    (57, "La", 0.00, 0.00, 1, 2),
    (58, "Ce", 0.00, 0.00, 1, 2),
    (59, "Pr", 0.00, 0.00, 1, 2),
    (60, "Nd", 0.00, 0.00, 1, 2),
    (61, "Pm", 0.00, 0.00, 1, 2),
    (62, "Sm", 0.00, 0.00, 1, 2),
    (63, "Eu", 0.00, 0.00, 1, 2),
    (64, "Gd", 0.00, 0.00, 1, 2),
    (65, "Tb", 0.00, 0.00, 1, 2),
    (66, "Dy", 0.00, 0.00, 1, 2),
    (67, "Ho", 0.00, 0.00, 1, 2),
    (68, "Er", 0.00, 0.00, 1, 2),
    (69, "Tm", 0.00, 0.00, 1, 2),
    (70, "Yb", 0.00, 0.00, 1, 2),
    (71, "Lu", 0.00, 0.00, 1, 2),
    (72, "Hf", 0.00, 0.00, 1, 2),
    (73, "Ta", 0.00, 0.00, 1, 2),
    (74, "W", 0.00, 0.00, 1, 2),
    (75, "Re", 0.00, 0.00, 1, 2),
    (76, "Os", 0.00, 0.00, 1, 2),
    (77, "Ir", 0.00, 0.00, 1, 2),
    (78, "Pt", 1.75, 1.66, 1, 2),
    (79, "Au", 1.66, 1.58, 1, 2),
    (80, "Hg", 1.55, 1.47, 1, 2),
    (81, "Tl", 1.96, 1.86, 1, 2),
    (82, "Pb", 2.02, 1.92, 1, 2),
    (83, "Bi", 0.00, 0.00, 1, 2),
    (84, "Po", 0.00, 0.00, 1, 2),
    (85, "At", 0.00, 0.00, 1, 2),
    (86, "Rn", 0.00, 0.00, 1, 2),
    (87, "Fr", 0.00, 0.00, 1, 3),
    (88, "Ra", 0.00, 0.00, 1, 3),
    (89, "Ac", 0.00, 0.00, 1, 3),
    (90, "Th", 0.00, 0.00, 1, 3),
    (91, "Pa", 0.00, 0.00, 1, 3),
    (92, "U", 1.86, 1.77, 1, 3),
    (93, "Np", 0.00, 0.00, 1, 3),
    (94, "Pu", 0.00, 0.00, 1, 3),
    (95, "Am", 0.00, 0.00, 1, 3),
    (96, "Cm", 0.00, 0.00, 1, 3),
    (97, "Bk", 0.00, 0.00, 1, 3),
    (98, "Cf", 0.00, 0.00, 1, 3),
    (99, "Es", 0.00, 0.00, 1, 3),
    (100, "Fm", 0.00, 0.00, 1, 3),
    (101, "Md", 0.00, 0.00, 1, 3),
    (102, "No", 0.00, 0.00, 1, 3),
    (103, "Lr", 0.00, 0.00, 1, 3),
    (104, "Rf", 0.00, 0.00, 1, 3),
    (105, "Db", 0.00, 0.00, 1, 3),
    (106, "Sg", 0.00, 0.00, 1, 3),
    (107, "Bh", 0.00, 0.00, 1, 3),
    (108, "Hs", 0.00, 0.00, 1, 3),
    (109, "Mt", 0.00, 0.00, 1, 3),
    (110, "Ds", 0.00, 0.00, 1, 3),
    (111, "Rg", 0.00, 0.00, 1, 3),
    (112, "Cn", 0.00, 0.00, 1, 3),
    (113, "Nh", 0.00, 0.00, 1, 3),
    (114, "Fl", 0.00, 0.00, 1, 3),
    (115, "Mc", 0.00, 0.00, 1, 3),
    (116, "Lv", 0.00, 0.00, 1, 3),
)

_LEGACY_SYMBOLS = {
    "Uub": "Cn",
    "Uut": "Nh",
    "Uuq": "Fl",
    "Uup": "Mc",
    "Uuh": "Lv",
}

RESP_PARAMETERS_BY_SYMBOL = {row[1]: RespParameter(*row) for row in _PARAMETER_ROWS}


def normalize_element_symbol(symbol: str) -> str:
    stripped = str(symbol).strip()
    if not stripped:
        raise ValueError("Element symbol must not be empty")
    normalized = stripped[0].upper() + stripped[1:].lower()
    return _LEGACY_SYMBOLS.get(normalized, normalized)


def get_resp_parameter(symbol: str) -> RespParameter:
    normalized = normalize_element_symbol(symbol)
    try:
        return RESP_PARAMETERS_BY_SYMBOL[normalized]
    except KeyError as exc:
        raise ValueError(f"No default RESP parameters are available for element {symbol!r}") from exc


def get_resp_mk_radius(symbol: str) -> float:
    parameter = get_resp_parameter(symbol)
    if parameter.mk_radius > 0:
        return float(parameter.mk_radius)
    return float(DEFAULT_MK_RADIUS)


def get_resp_radius_overrides(symbols) -> dict[str, float]:
    return {normalize_element_symbol(symbol): get_resp_mk_radius(symbol) for symbol in set(symbols)}


def select_resp_basis_family(symbols) -> str:
    parameters = tuple(get_resp_parameter(symbol) for symbol in symbols)
    basis_id = max(parameter.basis_set for parameter in parameters)
    # The historical table assigns first-row transition metals to the same
    # numeric family as organic elements.  That would silently choose an
    # all-electron 6-31G* setup for Fe/Zn.  Use the explicit SDD/ECP family for
    # the K-Zn block instead; mixed models then resolve one coherent family.
    if any(19 <= parameter.atomic_number <= 30 for parameter in parameters):
        basis_id = max(basis_id, 3)
    try:
        return BASIS_ID_TO_FAMILY[basis_id]
    except KeyError as exc:
        raise ValueError(f"Unknown RESP basis family id {basis_id}") from exc
