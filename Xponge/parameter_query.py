"""Read-only queries against the live Amber bonded-parameter registries."""
from itertools import product


def amber_parameter_exists(section, atom_types):
    """Check an ordered type tuple, including X wildcards and type case.

    Improper permutations are supplied by the caller, as in XpongeCPP. This
    query does not load a force field, create parameters, or change topology.
    """
    from .forcefield.base import angle_base, bond_base, dihedral_base

    specifications = {
        "BOND": (bond_base.BondType, 2),
        "ANGLE": (angle_base.AngleType, 3),
        "DIHE": (dihedral_base.ProperType, 4),
        "IMPROPER": (dihedral_base.ImproperType, 4),
    }
    if section not in specifications:
        raise ValueError("unsupported Amber section or atom-type count")
    cls, width = specifications[section]
    if len(atom_types) != width:
        raise ValueError("unsupported Amber section or atom-type count")
    registry = cls.get_all_types()
    orders = [tuple(atom_types)]
    if section != "IMPROPER":
        orders.append(tuple(reversed(atom_types)))
    return any("-".join(candidate) in registry
               for ordered in orders
               for candidate in product(*[(name, "X") for name in ordered]))
