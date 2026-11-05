"""Validated explicit SPONGE topology subset, in angstrom/radian/kcal/mol.

Periodic torsions are ordered energy terms, including periodic impropers;
they are never re-inferred from a chemical graph. See docs/explicit-topology.md.
"""
import math
import operator
import hashlib


def _index(value, count):
    if isinstance(value, bool):
        raise ValueError("Boolean atom index")
    try:
        value = operator.index(value)
    except TypeError as exc:
        raise ValueError("Expected an integer atom index") from exc
    if not 0 <= value < count:
        raise ValueError("Atom index out of range")
    return value


def _number(value):
    if isinstance(value, (bool, str)):
        raise ValueError("Expected a real parameter")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Nonfinite parameter")
    return value


def _validate(topology, count):
    widths = {"bonds": (2, 4), "angles": (3, 5),
              "dihedrals": (4, 7), "nb14s": (2, 4)}
    if set(topology) != set(widths) | {"exclusions"}:
        raise ValueError("Expected bonds, angles, dihedrals, nb14s and exclusions")
    result = {}
    for key, (natom, width) in widths.items():
        result[key] = []
        for row in topology[key]:
            if len(row) != width:
                raise ValueError("Invalid " + key + " row width")
            ids = [_index(v, count) for v in row[:natom]]
            if len(set(ids)) != natom:
                raise ValueError("Repeated atom in interaction")
            params = [_number(v) for v in row[natom:]]
            if key == "dihedrals":
                params[0] = _index(row[4], 2**31)
                if params[0] < 1:
                    raise ValueError("Nonpositive periodicity")
            elif params[0] < 0:
                raise ValueError("Negative force constant/scaling")
            if key == "bonds" and params[1] <= 0:
                raise ValueError("Nonpositive bond length")
            if key == "angles" and not 0 <= params[1] <= math.pi:
                raise ValueError("Angle outside [0, pi]")
            if key == "nb14s" and params[1] < 0:
                raise ValueError("Negative electrostatic scaling")
            result[key].append(ids + params)
    rows = topology["exclusions"]
    if len(rows) != count:
        raise ValueError("Exclusions must have one row per atom")
    exclusions = []
    for i, row in enumerate(rows):
        ids = [_index(v, count) for v in row]
        if i in ids or len(set(ids)) != len(ids):
            raise ValueError("Self or duplicate exclusion")
        exclusions.append(set(ids))
    for i, row in enumerate(exclusions):
        if any(i not in exclusions[j] for j in row):
            raise ValueError("Exclusions must be symmetric")
    result["exclusions"] = [sorted(row) for row in exclusions]
    return result


def _lj_identity(epsilon, rmin):
    epsilon, rmin = _number(epsilon), _number(rmin)
    if epsilon < 0 or rmin < 0:
        raise ValueError("Negative LJ parameter")
    name = "EXPLICIT_LJ_" + hashlib.sha256(
        (epsilon.hex() + ":" + rmin.hex()).encode()).hexdigest()
    return name, epsilon, rmin


def get_explicit_topology(molecule):
    """Return the five materialized topology tables (other force groups remain separate)."""
    from .build import build_bonded_force
    from .forcefield.base.exclude_base import Exclude
    build_bonded_force(molecule)
    molecule.get_atoms()
    indices = {atom: i for i, atom in enumerate(molecule.atoms)}
    result = {k: [] for k in ("bonds", "angles", "dihedrals", "nb14s")}
    for kind, key, names in (("bond", "bonds", ("k", "b")),
                             ("angle", "angles", ("k", "b")),
                             ("nb14", "nb14s", ("kLJ", "kee"))):
        for term in molecule.bonded_forces.get(kind, []):
            result[key].append([*[indices[a] for a in term.atoms],
                                *[getattr(term, name) for name in names]])
    for term in molecule.bonded_forces.get("dihedral", []):
        for n, k, phase in zip(term.periodicitys, term.ks, term.phi0s):
            result["dihedrals"].append([*[indices[a] for a in term.atoms], n, k, phase])
    for term in molecule.bonded_forces.get("improper", []):
        result["dihedrals"].append([*[indices[a] for a in term.atoms],
                                    term.periodicity, term.k, term.phi0])
    if Exclude.current is None:
        raise ValueError("Load a force field with an exclusion policy first")
    excluded = Exclude.current.get_excluded_atoms(molecule)
    result["exclusions"] = [sorted(indices[a] for a in excluded[atom]) for atom in molecule.atoms]
    return result


def set_explicit_topology(molecule, topology):
    """Set materialized energy terms without changing global bonded type registries.

    The supplied exclusions must include exclusions implied by the current
    linked-atom graph. This API does not change chemical connectivity.
    """
    from .forcefield.base import bond_base, angle_base, dihedral_base, nb14_base
    from .forcefield.base.exclude_base import Exclude
    atoms = [a for residue in molecule.residues for a in residue.atoms]
    top = _validate(topology, len(atoms))
    if Exclude.current is None:
        raise ValueError("Load a force field with an exclusion policy first")
    indices = {a: i for i, a in enumerate(atoms)}
    for i, atom in enumerate(atoms):
        implied = set(atom.linked_atoms.get("v", ()))
        for distance in range(2, Exclude.current.n + 1):
            implied.update(atom.linked_atoms.get(distance, ()))
        if not {indices[a] for a in implied} <= set(top["exclusions"][i]):
            raise ValueError("Cannot remove graph-implied exclusions")
    forces = dict(molecule.bonded_forces)
    for key in ("bond", "angle", "dihedral", "improper", "nb14"):
        forces[key] = []
    for key, kind, cls, names, width in (
        ("bonds", "bond", bond_base.BondType, ("k", "b"), 2),
        ("angles", "angle", angle_base.AngleType, ("k", "b"), 3),
        ("nb14s", "nb14", nb14_base.NB14Type, ("kLJ", "kee"), 2),
    ):
        for row in top[key]:
            term = cls.entity([atoms[i] for i in row[:width]], cls.get_type("UNKNOWNS"))
            for name, value in zip(names, row[width:]):
                setattr(term, name, value)
            forces[kind].append(term)
    # SPONGE uses the same ordered periodic potential for proper and improper
    # terms. Store one entity per row; never infer/improper-permute its order.
    for row in top["dihedrals"]:
        cls = dihedral_base.ProperType
        term = cls.entity([atoms[i] for i in row[:4]], cls.get_type("UNKNOWNS"))
        term.multiple_numbers = 1
        term.periodicitys, term.ks, term.phi0s = [row[4]], [row[5]], [row[6]]
        forces["dihedral"].append(term)
    # No mutation until validation and entity construction have succeeded.
    molecule.atoms, molecule.atom_index = atoms, indices
    molecule.bonded_forces = forces
    for i, atom in enumerate(atoms):
        atom.linked_atoms["extra_excluded_atoms"] = {atoms[j] for j in top["exclusions"][i]}
    molecule.built = True
    molecule._fep_materialized_topology = True


def set_atom_lj(molecule, atom_index, epsilon, rmin):
    """Assign LJ epsilon and Rmin/2 using immutable content-addressed type names."""
    from .forcefield.base.lj_base import LJType
    if not molecule.built or not getattr(molecule, "_fep_materialized_topology", False):
        raise ValueError("Set explicit topology before assigning atom LJ")
    atom_index = _index(atom_index, len(molecule.atoms))
    name, epsilon, rmin = _lj_identity(epsilon, rmin)
    key = name + "-" + name
    if key in LJType.get_all_types():
        old = LJType.get_type(key)
        if old.epsilon != epsilon or old.rmin != rmin:
            raise ValueError("Explicit LJ type collision")
    else:
        LJType(name=key, epsilon=epsilon, rmin=rmin)
    molecule.atoms[atom_index].LJtype = name
    return name
