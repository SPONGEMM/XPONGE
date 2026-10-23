"""AMBER periodic improper policy amber-deterministic-v1.

Preserve original parameter slots. Resolve equally matching peripheral atoms
by (actual type, residue position, unique local atom name), never object hashes
or export atom IDs. This intentionally does not emulate LEaP's unit-cache loss
of sOrder. Already materialized topologies and non-AMBER types are not reordered.
"""
from itertools import permutations

IMPROPER_ORDERING_POLICY = "amber-deterministic-v1"
_order = 0


def load_parameters(force_type, text):
    """Keep each periodicity and the original slots through alias-based updates."""
    global _order
    lines = [line for line in text.splitlines()
             if line.strip() and not line.lstrip().startswith(("!", "#", ";", "/"))]
    if not lines:
        return
    header = lines[0]
    for line in lines[1:]:
        existing = force_type._types_different_name.get(line.split()[0])
        records = getattr(existing, "_amber_improper_records", {}).copy()
        force_type.New_From_String(header + "\n" + line)
        pattern = tuple(line.split()[0].split("-"))
        parameter = force_type._types_different_name[line.split()[0]]
        _order += 1
        records[parameter.periodicity] = (pattern, parameter.k, parameter.phi0,
                                         parameter.periodicity, _order)
        parameter._amber_improper_records = records
        # A later generic/custom definition must not inherit AMBER metadata.
        parameter._parameter_update_callback = lambda p=parameter: p.__dict__.pop("_amber_improper_records", None)
    force_type._parameter_matcher = lambda candidates, owner: match_parameters(force_type, candidates, owner)


def match_parameters(force_type, candidates, owner):
    parameters = [p for p in force_type._types.values()
                  if getattr(p, "_amber_improper_records", None)]
    if not parameters:
        return None  # Preserve the ordinary path for other force fields.
    first = candidates[0]
    center = first[2]
    neighbors = [first[i] for i in (0, 1, 3)]
    residue_positions = {residue: i for i, residue in enumerate(getattr(owner, "residues", []))}

    def identity(atom):
        return (residue_positions.get(getattr(atom, "residue", None), 0), atom.name)

    options = []
    for outer in permutations(neighbors):
        atoms = [outer[0], outer[1], center, outer[2]]
        types = tuple(atom.type.name for atom in atoms)
        for parameter in parameters:
            for pattern, k, phase, n, order in parameter._amber_improper_records.values():
                if pattern[2] != types[2]:
                    continue
                if all(p == "X" or p == t for p, t in zip(pattern, types)):
                    key = tuple((atom.type.name, *identity(atom)) for atom in atoms)
                    options.append((pattern.count("X"), -order, key, n, parameter, atoms, k, phase))
    if not options:
        # Untagged custom/other-force-field definitions keep their legacy path.
        return None
    specificity = min(option[0] for option in options)
    # A generic/custom exact override must not be hidden by an older tagged
    # AMBER wildcard still present elsewhere in the same global registry.
    for parameter in force_type._types.values():
        if getattr(parameter, "_amber_improper_records", None):
            continue
        pattern = parameter.name.split("-")
        if len(pattern) != 4 or pattern.count("X") > specificity:
            continue
        for outer in permutations(neighbors):
            atoms = [outer[0], outer[1], center, outer[2]]
            if all(p == "X" or p == a.type.name for p, a in zip(pattern, atoms)):
                return None
    identities = [identity(atom) for atom in first]
    if len(set(identities)) != len(identities):
        raise ValueError("AMBER improper ordering requires unique local atom names")
    best = {}
    for option in options:
        if option[0] != specificity:
            continue
        n = option[3]
        if n not in best or option[:3] < best[n][:3]:
            best[n] = option
    return [(parameter, atoms, {"k": k, "phi0": phase, "periodicity": n})
            for _, _, _, n, parameter, atoms, k, phase in sorted(best.values(), key=lambda v: v[3])]
