"""Protect the process-wide Amber registry from incompatible base force fields."""


class ForceFieldFamilyConflictError(RuntimeError):
    """Raised when two alternative force fields target the same global family."""


_ACTIVE_FORCEFIELDS = {}


def activate_forcefield_family(family, forcefield):
    """Activate *forcefield* for *family*, rejecting incompatible replacements."""
    active = _ACTIVE_FORCEFIELDS.get(family)
    if active is None:
        _ACTIVE_FORCEFIELDS[family] = forcefield
        return forcefield
    if active != forcefield:
        raise ForceFieldFamilyConflictError(
            f"Cannot activate {forcefield}: force-field family {family!r} is already using {active}. "
            "Start a fresh Python process to select another base force field."
        )
    return active


def get_active_forcefield(family):
    """Return the active force field for *family*, or ``None`` when unset."""
    return _ACTIVE_FORCEFIELDS.get(family)


def require_forcefield_family(family, allowed):
    """Require an already active family member from *allowed*."""
    allowed = frozenset(allowed)
    active = get_active_forcefield(family)
    if active not in allowed:
        choices = ", ".join(sorted(allowed))
        if active is None:
            detail = "no base force field is active"
        else:
            detail = f"{active} is active"
        raise RuntimeError(
            f"Force-field extension requires family {family!r} to use one of [{choices}], but {detail}."
        )
    return active
