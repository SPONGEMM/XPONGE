"""Basis/ECP resolver for the default RESP ESP setup."""

from __future__ import annotations

from dataclasses import dataclass

from .resp_parameters import get_resp_parameter, normalize_element_symbol, select_resp_basis_family


@dataclass(frozen=True, slots=True)
class ResolvedRespBasis:
    label: str
    basis: str | dict[str, str]
    ecp: dict[str, str] | None
    cart: bool
    references: tuple[str, ...] = ()


def _normalize_basis_family(basis_family: str) -> str:
    value = str(basis_family).strip().lower().replace("_", "-")
    aliases = {
        "6-31g*": "6-31G*",
        "6-31g(d)": "6-31G*",
        "631g*": "6-31G*",
        "cep-31g": "CEP-31G",
        "cep31g": "CEP-31G",
        "sdd": "SDD",
    }
    return aliases.get(value, str(basis_family).strip())


def resolve_resp_basis(basis_family: str, elements) -> ResolvedRespBasis:
    family = _normalize_basis_family(basis_family)
    normalized_elements = {normalize_element_symbol(element) for element in elements}
    if family == "6-31G*":
        return ResolvedRespBasis(
            label="6-31G*",
            basis="6-31g*",
            ecp=None,
            cart=True,
            references=("HariharanPople1973_631Gstar", "Francl1982_631Gstar_second_row"),
        )
    if family == "CEP-31G":
        return _resolve_cep31g(normalized_elements)
    if family == "SDD":
        return _resolve_sdd(normalized_elements)
    return ResolvedRespBasis(label=str(basis_family), basis=str(basis_family), ecp=None, cart=False)


def resolve_default_resp_basis(elements) -> ResolvedRespBasis:
    return resolve_resp_basis(select_resp_basis_family(elements), elements)


def _resolve_cep31g(elements: set[str]) -> ResolvedRespBasis:
    basis: dict[str, str] = {}
    ecp: dict[str, str] = {}
    for symbol in sorted(elements, key=lambda item: get_resp_parameter(item).atomic_number):
        z = get_resp_parameter(symbol).atomic_number
        if z <= 2:
            basis[symbol] = "dz"
        else:
            basis[symbol] = "sbkjc"
            ecp[symbol] = "sbkjc"
    return ResolvedRespBasis(
        label="CEP-31G",
        basis=basis,
        ecp=ecp or None,
        cart=True,
        references=("StevensBaschKrauss1984_CEP", "StevensKraussBaschJasien1992_CEP"),
    )


def _resolve_sdd(elements: set[str]) -> ResolvedRespBasis:
    basis: dict[str, str] = {}
    ecp: dict[str, str] = {}
    for symbol in sorted(elements, key=lambda item: get_resp_parameter(item).atomic_number):
        z = get_resp_parameter(symbol).atomic_number
        if 1 <= z <= 10 or 13 <= z <= 17:
            basis[symbol] = "dz"
        elif z in {11, 12, 18}:
            basis[symbol] = "6-31g"
        elif 19 <= z <= 30 or 37 <= z <= 48:
            basis[symbol] = "stuttgart_rsc"
            ecp[symbol] = "stuttgart_rsc"
        elif 31 <= z <= 36 or 49 <= z <= 53:
            basis[symbol] = "stuttgart"
            ecp[symbol] = "stuttgart"
        else:
            raise ValueError(f"SDD RESP basis mapping is not defined for element {symbol}")
    return ResolvedRespBasis(
        label="SDD",
        basis=basis,
        ecp=ecp or None,
        cart=True,
        references=("DunningHay1977_DZ", "Andrae1990_StuttgartDresden"),
    )
