"""
This **module** gives the basic functions for fep calculations
"""
import re
import copy
import math
from itertools import combinations
from ...helper import source, set_global_alternative_names, Xdict, Xopen, Xprint, kabsch
from ..base import lj_base, exclude_base, bond_base, angle_base, dihedral_base, nb14_extra_base, soft_bond_base

source("....")

lj_type = lj_base.LJType

lj_type.New_From_String("""name A  B
ZERO_LJ_ATOM-ZERO_LJ_ATOM  0  0
""")

AtomType.Add_Property({"LJtypeB": str})
AtomType.Add_Property({"subsys": int})


def _find_common_forces(forcetype, aforces, bforces, mol_b2mol_a):
    """

    :param forcetype:
    :param aforces:
    :param bforces:
    :param mol_b2mol_a:
    :return:
    """
    # Match entities one-to-one. Multiple terms on the same atom tuple are
    # additive; never reuse a B entity for several A entities or overwrite it.
    candidates = {}
    consumed = set()
    for force in bforces:
        permutations = ([force.atoms, force.atoms[::-1]] if forcetype.get_class_name() == "improper"
                        else forcetype.Same_Force(force.atoms))
        keys = {tuple(mol_b2mol_a[atom] for atom in atoms) for atoms in permutations}
        for key in keys:
            candidates.setdefault(key, []).append(force)
    toret = []
    for force in aforces:
        match = next((other for other in candidates.get(tuple(force.atoms), ())
                      if other not in consumed), None)
        if match is not None:
            consumed.add(match)
        toret.append([force, match])
    toret.extend([None, force] for force in bforces if force not in consumed)
    return toret


FEP_BONDED_FORCE_MERGE_RULE = Xdict()


# pylint: disable=unused-argument
def _nb14_extra_merge_rule(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_, mol_r2mol_a, mol_a2mol_r,
                           mol_r2mol_b, mol_b2mol_r, lj_weights=None):
    """

    :param mol_r:
    :param mol_a:
    :param mol_b:
    :param forcetype:
    :param rforces:
    :param bforces:
    :param lambda_:
    :param mol_r2mol_a:
    :param mol_a2mol_r:
    :param mol_r2mol_b:
    :param mol_b2mol_r:
    :return:
    """
    coefficients = {}
    lj_weights = lj_weights or {}
    for forces, weight, atom_map in ((mol_a.bonded_forces.get("nb14_extra", []), 1 - lambda_, mol_a2mol_r),
                                      (bforces, lambda_, mol_b2mol_r)):
        for force in forces:
            atoms = tuple(sorted((atom_map[atom] for atom in force.atoms), key=mol_r.atom_index.get))
            values = coefficients.setdefault(atoms, [0., 0., 0., 0.])
            ee = force.nb14_ee_factor
            if ee is None:
                ee = force.kee * force.atoms[0].charge * force.atoms[1].charge
            lj_weight = lj_weights.get(force, weight)
            values[0] += lj_weight * force.A
            values[1] += lj_weight * force.B
            values[2] += weight * ee
            values[3] += weight * force.kee
    result = []
    for atoms, (a, b, ee, fallback) in coefficients.items():
        product = atoms[0].charge * atoms[1].charge
        if product == 0. and ee != 0.:
            raise ValueError("FEP NB14 Coulomb coefficient cannot be represented at zero charge product")
        kee = fallback if product == 0. else ee / product
        if not math.isfinite(kee):
            raise ValueError("non-finite FEP NB14 charge compensation")
        force = nb14_extra_base.NB14Type.entity(list(atoms), nb14_extra_base.NB14Type.get_type("UNKNOWNS"))
        force.A, force.B, force.kee, force.nb14_ee_factor = a, b, kee, ee
        result.append(force)
        atoms[0].Extra_Exclude_Atom(atoms[1])
    mol_r.bonded_forces["nb14_extra"] = result


FEP_BONDED_FORCE_MERGE_RULE["nb14_extra"] = {"lambda_name": "dihedral", "merge_function": _nb14_extra_merge_rule}


def _bond_base_merge_rule(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_,
                          mol_r2mol_a, mol_a2mol_r, mol_r2mol_b, mol_b2mol_r):
    """

    :param mol_r:
    :param mol_a:
    :param mol_b:
    :param forcetype:
    :param rforces:
    :param bforces:
    :param lambda_:
    :param mol_r2mol_a:
    :param mol_a2mol_r:
    :param mol_r2mol_b:
    :param mol_b2mol_r:
    :return:
    """
    forcepair = _find_common_forces(forcetype, rforces, bforces, mol_b2mol_r)
    for f_r, f_b in forcepair:
        if f_b is None:
            mol_r.bonded_forces["bond"].remove(f_r)
            f_r.from_AorB = 0
            mol_r.Add_Bonded_Force(f_r, "bond_soft")
        elif f_r is None:
            f_r = bond_base.BondType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type, f_b.name)
            f_r.k = f_b.k
            f_r.b = f_b.b
            f_r.from_AorB = 1
            mol_r.Add_Bonded_Force(f_r, "bond_soft")
        elif abs(f_r.b - f_b.b) < 1e-5:
            f_r.k = f_r.k * (1 - lambda_) + f_b.k * lambda_
        else:
            f_r2 = bond_base.BondType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type, f_b.name)
            f_r2.k = f_b.k * lambda_
            f_r2.b = f_b.b
            f_r.k *= (1 - lambda_)
            mol_r.Add_Bonded_Force(f_r2)


FEP_BONDED_FORCE_MERGE_RULE["bond"] = {"lambda_name": "bond", "merge_function": _bond_base_merge_rule}


def _angle_merge_rule(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_, mol_r2mol_a, mol_a2mol_r, mol_r2mol_b,
                      mol_b2mol_r):
    """

    :param mol_r:
    :param mol_a:
    :param mol_b:
    :param forcetype:
    :param rforces:
    :param bforces:
    :param lambda_:
    :param mol_r2mol_a:
    :param mol_a2mol_r:
    :param mol_r2mol_b:
    :param mol_b2mol_r:
    :return:
    """
    forcepair = _find_common_forces(forcetype, rforces, bforces, mol_b2mol_r)
    for f_r, f_b in forcepair:
        if f_b is None:
            f_r.k *= 1 - lambda_
        elif f_r is None:
            f_r = angle_base.AngleType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type, f_b.name)
            f_r.k *= lambda_
            f_r.b = f_b.b
            mol_r.Add_Bonded_Force(f_r)
        elif abs(f_r.b - f_b.b) < 1e-5:
            f_r.k = f_r.k * (1 - lambda_) + f_b.k * lambda_
        else:
            f_r2 = angle_base.AngleType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type, f_b.name)
            f_r2.k = f_b.k * lambda_
            f_r2.b = f_b.b
            f_r.k *= (1 - lambda_)
            mol_r.Add_Bonded_Force(f_r2)


FEP_BONDED_FORCE_MERGE_RULE["angle"] = {"lambda_name": "angle", "merge_function": _angle_merge_rule}


def _dihedral_merge_rule(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_, mol_r2mol_a, mol_a2mol_r,
                         mol_r2mol_b, mol_b2mol_r):
    """

    :param mol_r:
    :param mol_a:
    :param mol_b:
    :param forcetype:
    :param rforces:
    :param bforces:
    :param lambda_:
    :param mol_r2mol_a:
    :param mol_a2mol_r:
    :param mol_r2mol_b:
    :param mol_b2mol_r:
    :return:
    """
    forcepair = _find_common_forces(forcetype, rforces, bforces, mol_b2mol_r)
    for f_r, f_b in forcepair:
        if f_b is None:
            for i in range(f_r.multiple_numbers):
                f_r.ks[i] *= 1 - lambda_
        elif f_r is None:
            f_r2 = dihedral_base.ProperType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type, f_b.name)
            f_r2.multiple_numbers = f_b.multiple_numbers
            f_r2.ks = [f_b.ks[i] * lambda_ for i in range(f_b.multiple_numbers)]
            f_r2.phi0s = list(f_b.phi0s[:f_b.multiple_numbers])
            f_r2.periodicitys = list(f_b.periodicitys[:f_b.multiple_numbers])
            mol_r.Add_Bonded_Force(f_r2)
        else:
            sameforce = f_r.multiple_numbers == f_b.multiple_numbers
            if sameforce:
                for i in range(f_r.multiple_numbers):
                    # Coefficients are interpolated by index below. A matching
                    # set with a different order must use the superposition
                    # branch, not pair coefficients of different frequencies.
                    if (f_r.periodicitys[i] != f_b.periodicitys[i] or
                            abs(f_b.phi0s[i] - f_r.phi0s[i]) > 1e-5):
                        sameforce = False
                        break
            if sameforce:
                for i in range(f_r.multiple_numbers):
                    f_r.ks[i] = f_r.ks[i] * (1 - lambda_) + f_b.ks[i] * lambda_
            else:
                for i in range(f_r.multiple_numbers):
                    f_r.ks[i] *= (1 - lambda_)
                f_r2 = dihedral_base.ProperType.entity(list(map(lambda x: mol_b2mol_r[x], f_b.atoms)), f_b.type,
                                                       f_b.name)
                f_r2.multiple_numbers = f_b.multiple_numbers
                f_r2.ks = [f_b.ks[i] * lambda_ for i in range(f_b.multiple_numbers)]
                f_r2.phi0s = list(f_b.phi0s[:f_b.multiple_numbers])
                f_r2.periodicitys = list(f_b.periodicitys[:f_b.multiple_numbers])
                mol_r.Add_Bonded_Force(f_r2)


FEP_BONDED_FORCE_MERGE_RULE["dihedral"] = {"lambda_name": "dihedral", "merge_function": _dihedral_merge_rule}


def _improper_merge_rule(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_, mol_r2mol_a, mol_a2mol_r,
                         mol_r2mol_b, mol_b2mol_r):
    """

    :param mol_r:
    :param mol_a:
    :param mol_b:
    :param forcetype:
    :param rforces:
    :param bforces:
    :param lambda_:
    :param mol_r2mol_a:
    :param mol_a2mol_r:
    :param mol_r2mol_b:
    :param mol_b2mol_r:
    :return:
    """
    # Preserve each ordered Fourier term. Pairing only by atom tuple loses
    # multiplicities when several periodicities/phases share the same atoms.
    # Whole reversal is safe; other equivalent permutations may remain separate
    # weighted terms. Never reselect a physical endpoint on union atom IDs.
    result = {}
    indices = mol_r.atom_index
    for forces, weight, mapping in ((list(rforces), 1 - lambda_, None),
                                    (bforces, lambda_, mol_b2mol_r)):
        for source_force in forces:
            atoms = list(source_force.atoms) if mapping is None else [mapping[a] for a in source_force.atoms]
            ids = tuple(indices[a] for a in atoms)
            key = (min(ids, ids[::-1]), source_force.periodicity, source_force.phi0)
            if key not in result:
                force = dihedral_base.ImproperType.entity(atoms, source_force.type, source_force.name)
                force.periodicity, force.phi0 = source_force.periodicity, source_force.phi0
                force.k = 0.
                result[key] = force
            result[key].k += weight * source_force.k
    mol_r.bonded_forces["improper"] = [result[key] for key in sorted(result)]


FEP_BONDED_FORCE_MERGE_RULE["improper"] = {"lambda_name": "dihedral", "merge_function": _improper_merge_rule}


def save_hard_core_lj():
    """
    This **function** is used to save hard core lj

    :return: None
    """
    Molecule.Set_Save_SPONGE_Input("LJ")(lj_base.write_lj)
    from ..base import nb14_base
    Molecule._save_functions.pop("LJ_soft_core", None)
    Molecule._save_functions.pop("subsys_division", None)
    Molecule.Set_Save_SPONGE_Input("nb14")(nb14_base.write_nb14)
    Molecule.Set_Save_SPONGE_Input("nb14_extra")(nb14_extra_base.write_nb14)


def save_soft_core_lj():
    """
    This **function** is used to save soft core lj

    :return: None
    """
    Molecule._save_functions.pop("LJ", None)
    Molecule.Set_Save_SPONGE_Input("nb14")(lambda molecule: None)
    Molecule.Set_Save_SPONGE_Input("nb14_extra")(
        lambda molecule: nb14_extra_base.write_nb14(molecule, include_regular=True))

    def write_subsys_division(self):
        towrite = "%d\n" % len(self.atoms)
        for atom in self.atoms:
            if getattr(atom, "subsys", None) is None:
                towrite += "%d\n" % 0
            else:
                towrite += "%d\n" % atom.subsys
        return towrite

    Molecule.Set_Save_SPONGE_Input("subsys_division")(write_subsys_division)

    def write_soft_lj(self):
        """

        :param self:
        :return:
        """
        lj_types = []
        lj_typemap = Xdict()
        lj_typeb = []
        lj_typemapb = Xdict()
        for atom in self.atoms:
            if atom.LJtype not in lj_typemap.keys():
                lj_typemap[atom.LJtype] = len(lj_types)
                lj_types.append(atom.LJtype)
            if atom.LJtypeB is None:
                atom.LJtypeB = atom.LJtype
            if atom.LJtypeB not in lj_typemapb.keys():
                lj_typemapb[atom.LJtypeB] = len(lj_typeb)
                lj_typeb.append(atom.LJtypeB)

        # pylint: disable=protected-access
        lj_as, lj_bs = lj_base._find_ab_lj(lj_types)
        lj_asb, lj_bsb = lj_base._find_ab_lj(lj_typeb)
        checks = lj_base._get_checks(lj_types, lj_as, lj_bs)
        same_type = lj_base._judge_same_type(lj_types, checks)
        real_lj_types = lj_base._get_real_lj(lj_types, same_type)
        real_as, real_bs = lj_base._find_ab_lj(real_lj_types, False)

        checks = lj_base._get_checks(lj_typeb, lj_asb, lj_bsb)
        same_typeb = lj_base._judge_same_type(lj_typeb, checks)
        real_lj_typesb = lj_base._get_real_lj(lj_typeb, same_typeb)
        real_asb, real_bsb = lj_base._find_ab_lj(real_lj_typesb, False)

        towrite = "%d %d %d\n\n" % (len(self.atoms), len(real_lj_types), len(real_lj_typesb))
        count = 0
        for i in range(len(real_lj_types)):
            for _ in range(i + 1):
                towrite += "%16.7e" % real_as[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_lj_types)):
            for _ in range(i + 1):
                towrite += "%16.7e" % real_bs[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_lj_typesb)):
            for _ in range(i + 1):
                towrite += "%16.7e" % real_asb[count] + " "
                count += 1
            towrite += "\n"
        towrite += "\n"

        count = 0
        for i in range(len(real_lj_typesb)):
            for _ in range(i + 1):
                towrite += "%16.7e" % real_bsb[count] + " "
                count += 1
            towrite += "\n"

        towrite += "\n"
        towrite += "\n".join(
            ["%d %d" % (same_type[lj_typemap[atom.LJtype]], same_typeb[lj_typemapb[atom.LJtypeB]]) for atom in
             self.atoms])
        return towrite

    Molecule.Set_Save_SPONGE_Input("LJ_soft_core")(write_soft_lj)


def intramolecule_nb_to_nb14(mol_a, perturbating_residues):
    """
    This **function** convert the non bonded interactions to nb14 interactions within the molecule

    :param mol_a: the Molecule instance
    :param perturbating_residues: the residue(s) to be perturbed
    :return:
    """
    if isinstance(perturbating_residues, Residue):
        perturbating_residues = [perturbating_residues]
    build.Build_Bonded_Force(mol_a)
    a_exclude = exclude_base.Exclude.current.Get_Excluded_Atoms(mol_a)
    for residue1 in perturbating_residues:
        for atom_a1 in residue1.atoms:
            for residue2 in perturbating_residues:
                for atom_a2 in residue2.atoms:
                    if atom_a1 == atom_a2:
                        continue
                    if atom_a2 not in a_exclude[atom_a1]:
                        temp_a, temp_b = nb14_extra_base.get_nb14_extra_lj(atom_a1, atom_a2)
                        new_force = nb14_extra_base.NB14Type.entity([atom_a1, atom_a2],
                                                                    nb14_extra_base.NB14Type.get_type("UNKNOWNS"))
                        new_force.A = temp_a
                        new_force.B = temp_b
                        new_force.kee = 1
                        mol_a.Add_Bonded_Force(new_force)

                        atom_a1.Extra_Exclude_Atom(atom_a2)
                        a_exclude[atom_a1].add(atom_a2)
                        a_exclude[atom_a2].add(atom_a1)


def get_free_molecule(mol_a, perturbating_residues, intra_fep=False):
    """
    This **function** makes the molecule to be "free", having no interaction with other molecules

    :param mol_a: the Molecule instance
    :param perturbating_residues: the residues to be perturbed
    :param intra_fep: whether clear intramolecular non bonded interactions
    :return: a new Molecule instance
    """
    if isinstance(perturbating_residues, Residue):
        perturbating_residues = [perturbating_residues]
    build.Build_Bonded_Force(mol_a)
    nb14_extra_base.NB14_To_NB14EXTRA(mol_a)
    mol_b = mol_a.deepcopy()
    mol_a2mol_b = Xdict()
    for i, atom_a in enumerate(mol_a.atoms):
        mol_a2mol_b[atom_a] = mol_b.atoms[i]

    for residue in perturbating_residues:
        for atom_a in residue.atoms:
            atom = mol_a2mol_b[atom_a]
            atom.charge = 0
            atom.LJtype = "ZERO_LJ_ATOM"
            atom.subsys = 1
            atom_a.subsys = 1

    return mol_b


def _correct_residueb_coordinates(residue_a, residue_b, matchmap):
    """

    :param residue_a:
    :param residue_b:
    :param matchmap:
    :return:
    """
    template_positions = []
    certified_positions = []
    residue_type_a = residue_a.type
    for i, atom in enumerate(residue_b.atoms):
        if i in matchmap.keys():
            temp_atom = residue_a.name2atom(residue_type_a.atoms[matchmap[i]].name)
            template_positions.append([temp_atom.x, temp_atom.y, temp_atom.z])
            certified_positions.append([atom.x, atom.y, atom.z])

    rotation, center1, center2 = kabsch(template_positions, certified_positions)

    for atom in residue_b.atoms:
        temp_position = np.array([getattr(atom, xyz) for xyz in "xyz"])
        tnew_position = np.dot(rotation, (temp_position - center2)) + center1
        atom.x = tnew_position[0]
        atom.y = tnew_position[1]
        atom.z = tnew_position[2]


def _get_extra_atoms_and_rbmap(restype_ab, residue_type_a, residue_type_b, residue_a,
                               forcopy, matchmap, match_a, match_b):
    """
        get the extra atoms of the residue and the map from the mixed residue to residue b
    """
    extra_a = []
    extra_b = []
    rbmap = {value: key for key, value in matchmap.items()}
    for i, atom in enumerate(residue_type_a.atoms):
        atom.copied[forcopy].x = residue_a.name2atom(atom.name).x
        atom.copied[forcopy].y = residue_a.name2atom(atom.name).y
        atom.copied[forcopy].z = residue_a.name2atom(atom.name).z
        if i not in match_a:
            extra_a.append(atom.copied[forcopy])
            extra_a[-1].subsys = 1

    for i, atom in enumerate(residue_type_b.atoms):
        if i not in match_b:
            rbmap[len(restype_ab.atoms)] = i
            restype_ab.Add_Atom(atom.name + "R2", atom.type, atom.x, atom.y, atom.z)
            atom.copied[forcopy] = restype_ab.atoms[-1]
            atom.copied[forcopy].contents = dict(atom.contents.items())
            atom.copied[forcopy].name = atom.name + "R2"
            extra_b.append(atom.copied[forcopy])
            extra_b[-1].subsys = 2
        else:
            residue_type_b.atoms[i].copied[forcopy] = restype_ab.atoms[matchmap[i]]
    return extra_a, extra_b, rbmap


def _link_restypeb_atoms(residue_type_b, forcopy, matchmap):
    """

    :param residue_type_b:
    :param forcopy:
    :param matchmap:
    :return:
    """
    for atom in residue_type_b.atoms:
        for key, _ in atom.copied[forcopy].linked_atoms.items():
            for aton in atom.linked_atoms.get(key, []):
                if not (residue_type_b.atom2index(aton) in matchmap.keys()
                        and residue_type_b.atom2index(atom) in matchmap.keys()):
                    atom.copied[forcopy].Link_Atom(key, aton.copied[forcopy])


def _make_unique_fep_restype_name(base_name, nonce):
    """Return a ResidueType registry-safe temporary name for FEP mixed topologies."""
    candidate = f"{base_name}__FEP_TMP_{nonce}"
    suffix = 0
    while candidate in ResidueType._types:
        suffix += 1
        candidate = f"{base_name}__FEP_TMP_{nonce}_{suffix}"
    return candidate


def _get_residue_ab(residue_type_a, residue_type_b, residue_a, restype_name,
                    forcopy, matchmap, match_a, match_b):
    """
        get the mixed residue type
    """
    restype_ab = residue_type_a.deepcopy(restype_name, forcopy)

    extra_a, extra_b, rbmap = _get_extra_atoms_and_rbmap(restype_ab, residue_type_a, residue_type_b, residue_a, forcopy,
                                                         matchmap, match_a, match_b)

    for atomi in extra_a:
        for atomj in extra_b:
            atomi.Extra_Exclude_Atom(atomj)

    for atom, connect_set in residue_type_b.connectivity.items():
        for aton in connect_set:
            restype_ab.Add_Connectivity(atom.copied[forcopy], aton.copied[forcopy])

    for bond_entities in residue_type_b.bonded_forces.values():
        for bond_entity in bond_entities:
            tocopy = False
            for atom in bond_entity.atoms:
                if residue_type_b.atom2index(atom) not in matchmap.keys():
                    tocopy = True
                    break
            if tocopy:
                restype_ab.Add_Bonded_Force(bond_entity.deepcopy(forcopy))

    _link_restypeb_atoms(residue_type_b, forcopy, matchmap)

    return restype_ab, rbmap


def _copy_atom_contents(target_atom, source_atom):
    """Copy the atom typing and parameters from one atom to another."""
    target_atom.type = source_atom.type
    for key, value in source_atom.contents.items():
        if key != "name":
            target_atom.contents[key] = value
    target_atom.charge = source_atom.charge
    target_atom.LJtype = source_atom.LJtype


def _set_zero_lj_atom(atom, subsys):
    """Hide one endpoint atom by clearing its non-bonded state."""
    atom.LJtype = "ZERO_LJ_ATOM"
    atom.charge = 0
    atom.contents["subsys"] = subsys


def _project_bonded_force(force, atom_map):
    """Copy a bonded force entity onto another residue type via an atom map."""
    new_force = force.type.entity([atom_map[atom] for atom in force.atoms], force.type, force.name)
    new_force.contents = copy.deepcopy(force.contents)
    return new_force


def _project_physical_endpoint(target, sources, auxiliary_sources=()):
    """Project physical terms plus unchanged opposite-branch auxiliary terms.

    A dual-topology dummy loses nonbonded interactions, not its bonded scaffold.
    Copy only source entities touching an exclusive atom; never re-parameterize
    the union graph (which could invent cross-branch angles or impropers).
    """
    target.get_atoms()
    target.bonded_forces = {}
    for atom in target.atoms:
        atom.linked_atoms = Xdict({i: set() for i in range(2, GlobalSetting.farthest_bonded_force + 1)})
        atom.linked_atoms["extra_excluded_atoms"] = set()
    for source, atom_map in sources:
        for key, forces in source.bonded_forces.items():
            for force in forces:
                if all(atom in atom_map for atom in force.atoms):
                    target.Add_Bonded_Force(_project_bonded_force(force, atom_map), key)
        for atom, mapped in atom_map.items():
            for distance, neighbors in atom.linked_atoms.items():
                mapped.linked_atoms.setdefault(distance, set()).update(
                    atom_map[other] for other in neighbors if other in atom_map)
    for source, atom_map, exclusive in auxiliary_sources:
        for key, forces in source.bonded_forces.items():
            for force in forces:
                if any(atom in exclusive for atom in force.atoms):
                    target.Add_Bonded_Force(_project_bonded_force(force, atom_map), key)
        for atom, mapped in atom_map.items():
            for distance, neighbors in atom.linked_atoms.items():
                mapped.linked_atoms.setdefault(distance, set()).update(
                    atom_map[other] for other in neighbors
                    if other in atom_map and (atom in exclusive or other in exclusive))
    target.built = True
    target._fep_materialized_topology = True


def _nhydrogen(rdmol, j):
    """get the number of hydrogens of the j-th atom in an rdkit rdmol instance"""
    return rdmol.GetAtomWithIdx(j).GetTotalNumHs(includeNeighbors=True)


def merge_dual_topology(mol, residue_a, residue_b, assign_a, assign_b,
                        tmcs=60, image_path=None, similarity_limit=0, imcs=None):
    """
    This **function** perturbs a residue in the molecule into another type in the dual topology way

    :param mol: the Molecule instance
    :param residue_a: the Residue which needs to be perturbed
    :param residue_b: the Residue which is perturbed to
    :param assign_a: the Assign instance corresponding to ``residue_a``
    :param assign_b: the Assign instance corresponding to ``residue_b``
    :param tmcs: the max time to find the max common structure
    :param image_path: the path to save the mcs image
    :param similarity_limit: the limitation of the similarity. The similarity is calculated by \
the tanimoto coefficient of the max common structure.
    :param imcs: the input file to store the mcs information
    :return: two molecules in the initial and final lambda stat respectively, and the matchmap
    """
    residue_index = mol.residues.index(residue_a)
    mol = mol.deepcopy()
    residue_a = mol.residues[residue_index]
    # A materialized single-residue Molecule carries externally assigned
    # parameters that cannot be recovered from its ResidueType registry.
    replacement_molecule = None
    if isinstance(residue_b, Molecule):
        if len(residue_b.residues) != 1:
            raise ValueError("FEP replacement Molecule must contain exactly one residue")
        replacement_molecule = residue_b.deepcopy()
        build.Build_Bonded_Force(replacement_molecule)
        residue_b = replacement_molecule.residues[0]
    replacement_input = residue_b
    replacement_type = residue_b.type if isinstance(residue_b, Residue) else residue_b
    replacement_public_name = replacement_type.name
    replacement_name = _make_unique_fep_restype_name(replacement_public_name, time.time_ns())
    residue_b = replacement_type.deepcopy(replacement_name)
    if isinstance(replacement_input, Residue):
        for atom in replacement_input.atoms:
            _copy_atom_contents(residue_b.name2atom(atom.name), atom)
    build.Build_Bonded_Force(mol)
    build.Build_Bonded_Force(residue_b)
    target_atoms = set(residue_a.atoms)
    for forces in mol.bonded_forces.values():
        for force in forces:
            membership = [atom in target_atoms for atom in force.atoms]
            if any(membership) and not all(membership):
                raise ValueError("coupled FEP target requires explicit physical B environment topology")
    for atom in mol.atoms:
        if any((atom in target_atoms) != (other in target_atoms) for other in atom.extra_excluded_atoms):
            raise ValueError("coupled FEP exclusions require explicit physical B environment topology")
    if (mol.bonded_forces.get("nb14_extra") or residue_b.bonded_forces.get("nb14_extra")
            or (replacement_molecule is not None and replacement_molecule.bonded_forces.get("nb14_extra"))):
        raise ValueError("dual topology of pre-existing explicit pairs is not supported; convert after dual assembly")

    from ...helper.rdkit import assign_to_rdmol
    from rdkit.Chem import rdFMCS as MCS
    from rdkit.Chem import Draw, AllChem, RemoveHs, SanitizeMol, rdmolops

    rdmol_a = assign_to_rdmol(assign_a)
    rdmol_b = assign_to_rdmol(assign_b)
    if imcs is None:
        Xprint("FINDING MAXIMUM COMMON SUBSTRUCTURE\n")
        flags = rdmolops.SanitizeFlags
        rdmol_a_no_h = RemoveHs(rdmol_a, sanitize=False)
        SanitizeMol(rdmol_a_no_h, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
        rdmol_b_no_h = RemoveHs(rdmol_b, sanitize=False)
        SanitizeMol(rdmol_b_no_h, flags.SANITIZE_ALL & ~flags.SANITIZE_PROPERTIES)
        result = MCS.FindMCS([rdmol_a_no_h, rdmol_b_no_h],
                             completeRingsOnly=True,
                             bondCompare=MCS.BondCompare.CompareOrderExact,
                             timeout=tmcs)
        rdmol_mcs = result.queryMol
        if rdmol_mcs is None:
            raise OSError("No common substructure found")
        match_a = rdmol_a.GetSubstructMatch(rdmol_mcs)
        a_hydros = {i: _nhydrogen(rdmol_a, j) for i, j in enumerate(match_a)}
        match_b = sorted(rdmol_b.GetSubstructMatches(rdmol_mcs, uniquify=False),
                         key=lambda x: sum(min(_nhydrogen(rdmol_b, j), a_hydros[i]) for i,j in enumerate(x)))[-1]
        matchmap = {match_b[j]: match_a[j] for j in range(len(match_a))}
        extra_map = {}
        for j, i in matchmap.items():
            a = rdmol_a.GetAtomWithIdx(i)
            a_hydrogens = []
            for ai in a.GetNeighbors():
                if ai.GetAtomicNum() == 1:
                    a_hydrogens.append(ai.GetIdx())
            b = rdmol_b.GetAtomWithIdx(j)
            for bi in b.GetNeighbors():
                if not a_hydrogens:
                    break
                if bi.GetAtomicNum() == 1:
                    extra_map[bi.GetIdx()] = a_hydrogens.pop()
        matchmap.update(extra_map)
        match_a, match_b = zip(*[[value, key] for key, value in matchmap.items()])
        tanimoto = len(match_a)
        tanimoto /= len(assign_a.atoms) + len(assign_b.atoms) - tanimoto
        Xprint(f"similarity: {tanimoto}")
        if image_path:
            draw_a = assign_to_rdmol(assign_a)
            draw_b = assign_to_rdmol(assign_b)
            AllChem.Compute2DCoords(draw_a)
            AllChem.Compute2DCoords(draw_b)
            for atom in draw_a.GetAtoms():
                atom.SetProp("atomLabel", atom.GetSymbol())
            for atom in draw_b.GetAtoms():
                atom.SetProp("atomLabel", atom.GetSymbol())
            img = Draw.MolsToGridImage([draw_a, draw_b],
                                       molsPerRow=1,
                                       subImgSize=(1200, 600),
                                       highlightAtomLists=[match_a, match_b])
            img.save(image_path)
            f = Xopen(os.path.splitext(image_path)[0] + ".txt", "w")
            f.write(f"similarity: {tanimoto}\n")
            f.write(f"mcs atoms in r1: {match_a}\n\n")
            f.write(f"mcs atoms in r2: {match_b}\n\n")
            f.close()
    else:
        match_a = None
        match_b = None
        with open(imcs) as f:
            for line in f:
                if line.startswith("mcs atoms in"):
                    result = re.search(r"r(\d): (.*)", line)
                    if result is None or len(result.groups()) != 2:
                        raise ValueError(f"the format of the mcs file {imcs} is not right")
                    id_ = int(result.group(1))
                    atoms = tuple(int(i) for i in result.group(2).replace('(', "")
                                                            .replace(')', "")
                                                            .split(",") if i.strip())
                    if id_ == 1:
                        match_a = atoms
                    elif id_ == 2:
                        match_b = atoms
                    else:
                        raise ValueError(f"the format of the mcs file {imcs} is not right")
        if match_a is None or match_b is None or len(match_a) != len(match_b):
            raise ValueError(f"the format of the mcs file {imcs} is not right")
        matchmap = {match_b[j]: match_a[j] for j in range(len(match_a))}
        tanimoto = len(match_a)
        tanimoto /= len(assign_a.atoms) + len(assign_b.atoms) - tanimoto
        Xprint(f"similarity: {tanimoto}")
    if tanimoto <= similarity_limit:
        Xprint(f"similarity (={tanimoto}) should be greater than {similarity_limit}", "ERROR")
        sys.exit(2)
    Xprint("ALIGNING TOPOLOGY AND COORDINATE\n")

    residue_type_a = residue_a.type

    if isinstance(residue_b, Residue):
        residue_type_b = residue_b.type
    elif isinstance(residue_b, ResidueType):
        residue_type_b = residue_b
    else:
        raise TypeError

    _correct_residueb_coordinates(residue_a, residue_type_b, matchmap)

    forcopy = hash(str(time.time()))
    canonical_ab_name = residue_type_a.name + "_" + replacement_public_name
    canonical_ba_name = replacement_public_name + "_" + residue_type_a.name
    nonce = hex(forcopy & 0xffffffffffffffff)[2:]
    restype_ab_name = _make_unique_fep_restype_name(canonical_ab_name, nonce)
    restype_ab, rbmap = _get_residue_ab(residue_type_a, residue_type_b, residue_a, restype_ab_name,
                                        forcopy, matchmap, match_a, match_b)

    restype_ba_name = _make_unique_fep_restype_name(canonical_ba_name, nonce)
    restype_ba = restype_ab.deepcopy(restype_ba_name)

    for i, atomi in enumerate(restype_ab.atoms):
        if i < len(residue_type_a.atoms):
            _copy_atom_contents(atomi, residue_a.name2atom(residue_type_a.atoms[i].name))
        else:
            _set_zero_lj_atom(restype_ab.atoms[i], 2)

        if i in rbmap:
            _copy_atom_contents(restype_ba.atoms[i], residue_type_b.atoms[rbmap[i]])
        else:
            _set_zero_lj_atom(restype_ba.atoms[i], 1)
        subsys = 2 if i >= len(residue_type_a.atoms) else int(i not in rbmap)
        restype_ab.atoms[i].contents["subsys"] = subsys
        restype_ba.atoms[i].contents["subsys"] = subsys
    mol_a = Molecule(mol.name + "A")
    mol_b = Molecule(mol.name + "B")

    ri = 0
    res2index = Xdict()
    for i, res in enumerate(mol.residues):
        res2index[res] = i
        if res == residue_a:
            mol_a.Add_Residue(restype_ab)
            mol_b.Add_Residue(restype_ba)
            mol_a.residues[-1].name = canonical_ab_name
            mol_b.residues[-1].name = canonical_ba_name
            ri = i
        else:
            mol_a.Add_Residue(res.deepcopy())
            mol_b.Add_Residue(res.deepcopy())

    for reslink in mol.residue_links:
        atom1 = reslink.atom1
        atom2 = reslink.atom2
        new_atom1_a = mol_a.residues[res2index[atom1.residue]].name2atom(atom1.name)
        new_atom1_b = mol_b.residues[res2index[atom1.residue]].name2atom(atom1.name)
        new_atom2_a = mol_a.residues[res2index[atom2.residue]].name2atom(atom2.name)
        new_atom2_b = mol_b.residues[res2index[atom2.residue]].name2atom(atom2.name)
        if atom1.residue == residue_a:
            new_atom1_a = mol_a.residues[ri].name2atom(new_atom1_a.name)
            new_atom1_b = mol_b.residues[ri].name2atom(new_atom1_b.name)
        if atom2.residue == residue_a:
            new_atom2_a = mol_a.residues[ri].name2atom(new_atom2_a.name)
            new_atom2_b = mol_b.residues[ri].name2atom(new_atom2_b.name)
        mol_a.add_residue_link(new_atom1_a, new_atom2_a)
        mol_b.add_residue_link(new_atom1_b, new_atom2_b)
    a_map = {atom: mol_a.residues[index].name2atom(atom.name)
             for index, residue in enumerate(mol.residues) for atom in residue.atoms}
    environment_map = {atom: mol_b.residues[index].name2atom(atom.name)
                       for index, residue in enumerate(mol.residues) if index != ri for atom in residue.atoms}
    physical_b = (replacement_molecule if replacement_molecule is not None
                  else Molecule.cast(residue_type_b, deepcopy=True))
    build.Build_Bonded_Force(physical_b)
    b_map = {physical_b.atoms[b_index]: mol_b.residues[ri].atoms[a_index]
             for a_index, b_index in rbmap.items()}
    a_in_b = {atom: mol_b.residues[index].name2atom(atom.name)
              for index, residue in enumerate(mol.residues) for atom in residue.atoms}
    b_in_a = {physical_b.atoms[b_index]: mol_a.residues[ri].atoms[a_index]
              for a_index, b_index in rbmap.items()}
    exclusive_a = {residue_a.atoms[i] for i in range(len(residue_a.atoms)) if i not in match_a}
    exclusive_b = {physical_b.atoms[i] for i in range(len(physical_b.atoms)) if i not in match_b}
    _project_physical_endpoint(mol_a, [(mol, a_map)], [(physical_b, b_in_a, exclusive_b)])
    _project_physical_endpoint(mol_b, [(mol, environment_map), (physical_b, b_map)],
                               [(mol, a_in_b, exclusive_a)])
    for endpoint in (mol_a, mol_b):
        a_only = [endpoint.residues[ri].atoms[i] for i in range(len(residue_a.atoms)) if i not in match_a]
        b_only = [endpoint.residues[ri].atoms[i] for i, j in rbmap.items() if j not in match_b]
        for atom_a in a_only:
            for atom_b in b_only:
                atom_a.extra_excluded_atoms.add(atom_b)
                atom_b.extra_excluded_atoms.add(atom_a)
    for endpoint in (mol_a, mol_b):
        endpoint.box_length = copy.deepcopy(mol.box_length)
        endpoint.box_origin = copy.deepcopy(mol.box_origin)
        endpoint.box_angle = copy.deepcopy(mol.box_angle)
    return mol_a, mol_b, matchmap


SUPPORTS_BRANCH_NB14_LJ_WEIGHTS = True


def core_pair_transitions(endpoint_a, endpoint_b, residue_index):
    """List changed nonbonded categories using the physical union endpoint topology."""
    modes = []
    core_indices = []
    for endpoint in (endpoint_a, endpoint_b):
        if not endpoint.built:
            endpoint = endpoint.deepcopy()
            build.Build_Bonded_Force(endpoint)
        core = tuple(endpoint.atom_index[atom] for atom in endpoint.residues[residue_index].atoms
                     if not getattr(atom, "subsys", None))
        core_indices.append(core)
        excluded = exclude_base.Exclude.current.Get_Excluded_Atoms(endpoint)
        native, explicit = (
            {tuple(sorted(endpoint.atom_index[atom] for atom in force.atoms))
             for force in endpoint.bonded_forces.get(key, ())}
            for key in ("nb14", "nb14_extra")
        )
        state = {}
        for i, j in combinations(core, 2):
            pair = i, j
            first, second = endpoint.atoms[i], endpoint.atoms[j]
            if pair in native:
                mode = "native14+explicit" if pair in explicit else "native14"
            elif pair in explicit:
                mode = "explicit"
            elif second in excluded[first] or first in excluded[second]:
                mode = "excluded"
            else:
                mode = "ordinary"
            state[pair] = mode
        modes.append(state)
    if (len(endpoint_a.atoms) != len(endpoint_b.atoms) or core_indices[0] != core_indices[1]
            or any(endpoint_a.atoms[i].name != endpoint_b.atoms[i].name for i in core_indices[0])):
        raise ValueError("FEP core pair classification requires aligned dual topology endpoints")
    return [(i, j, mode_a, modes[1][i, j]) for (i, j), mode_a in modes[0].items()
            if mode_a != modes[1][i, j]]


def _native_nb14_lj_weights(endpoint, branch, suffix, specific_lambda):
    """Convert native pairs and assign independent weights to eligible LJ terms."""
    explicit = set(endpoint.bonded_forces.get("nb14_extra", ()))
    nb14_extra_base.nb14_to_nb14_extra(endpoint)
    junction = specific_lambda.get(f"nb14_junction_lj_{suffix}")
    internal = specific_lambda.get(f"nb14_internal_lj_{suffix}")
    if junction is None and internal is None:
        return {}
    for value in (junction, internal):
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError("Native 1-4 LJ weights must be finite and in [0, 1]")
    residues = {atom: residue for residue in endpoint.residues for atom in residue.atoms}
    weights = {}
    for force in endpoint.bonded_forces.get("nb14_extra", ()):
        if force in explicit:
            continue
        atom_a, atom_b = force.atoms
        if residues[atom_a] is not residues[atom_b] or any(
            atom.LJtype == "ZERO_LJ_ATOM" for atom in force.atoms
        ):
            continue
        groups = tuple(getattr(atom, "subsys", None) or 0 for atom in force.atoms)
        if junction is not None and set(groups) == {0, branch}:
            weights[force] = junction
        if internal is not None and groups == (branch, branch):
            weights[force] = internal
    return weights


def merge_force_field(mol_a, mol_b, default_lambda, specific_lambda=None, intra_fep=False):
    """
    This **function** merges one molecule in two different force fields (two Molecule instances) into one

    :param mol_a: the Molecule instance in the initial lambda stat
    :param mol_b: the Molecule instance in the final lambda stat
    :param default_lambda: the lambda to scale the force if no ``specific_lambda`` is set for the force
    :param specific_lambda: a dict to map the force to its special scale factor lambda
    :param intra_fep: whether clear intramolecular non bonded interactions
    :return: the Molecule instance merged
    """
    if specific_lambda is None:
        specific_lambda = Xdict()
    # Building force caches and converting NB14 representations must not mutate
    # either caller-owned endpoint. This also makes repeated lambda generation
    # independent of call order.
    endpoints = []
    for endpoint in (mol_a, mol_b):
        duplicate = endpoint.deepcopy()
        for original_residue, copied_residue in zip(endpoint.residues, duplicate.residues):
            copied_residue.name = original_residue.name
        endpoints.append(duplicate)
    mol_a, mol_b = endpoints
    build.Build_Bonded_Force(mol_a)
    build.Build_Bonded_Force(mol_b)

    lj_weights = _native_nb14_lj_weights(mol_a, 1, "a", specific_lambda)
    lj_weights.update(_native_nb14_lj_weights(mol_b, 2, "b", specific_lambda))

    assert len(mol_a.atoms) == len(mol_b.atoms)
    mol_a2mol_b = Xdict()
    for i, atom_a in enumerate(mol_a.atoms):
        mol_a2mol_b[atom_a] = mol_b.atoms[i]

    mol_r = mol_a.deepcopy()
    for i, residue in enumerate(mol_r.residues):
        residue.name = mol_a.residues[i].name

    mol_r2mol_a = Xdict()
    mol_r2mol_b = Xdict()
    for i, atom_ret in enumerate(mol_r.atoms):
        mol_r2mol_a[atom_ret] = mol_a.atoms[i]
        mol_r2mol_b[atom_ret] = mol_b.atoms[i]
    mol_b2mol_r = {value: key for key, value in mol_r2mol_b.items()}
    mol_a2mol_r = {value: key for key, value in mol_r2mol_a.items()}

    charge_lambda = specific_lambda.get("charge", default_lambda)
    mass_lambda = specific_lambda.get("mass", default_lambda)

    for i, atom in enumerate(mol_r.atoms):
        atom.charge = mol_r2mol_a[atom].charge * (1 - charge_lambda) + mol_r2mol_b[atom].charge * charge_lambda
        atom.mass = mol_r2mol_a[atom].mass * (1 - mass_lambda) + mol_r2mol_b[atom].mass * mass_lambda
        atom.contents["LJtypeB"] = mol_r2mol_b[atom].LJtype
        for distance, neighbors in mol_r2mol_b[atom].linked_atoms.items():
            atom.linked_atoms.setdefault(distance, set()).update(mol_b2mol_r[other] for other in neighbors)

    # Deepcopy retains populated force groups only. Include B-only groups so
    # their native merge rule still runs when A has no term of that kind.
    for forcename in mol_b.bonded_forces:
        mol_r.bonded_forces.setdefault(forcename, [])
    for forcename, rforces in list(mol_r.bonded_forces.items()):
        if forcename in FEP_BONDED_FORCE_MERGE_RULE.keys():
            temp_lambda = specific_lambda.get(FEP_BONDED_FORCE_MERGE_RULE[forcename].get("lambda_name"),
                                              default_lambda)
            bforces = mol_b.bonded_forces.get(forcename, [])
            if forcename == "nb14_extra":
                _nb14_extra_merge_rule(
                    mol_r, mol_a, mol_b, GlobalSetting.BondedForcesMap[forcename],
                    rforces, bforces, temp_lambda, mol_r2mol_a, mol_a2mol_r,
                    mol_r2mol_b, mol_b2mol_r, lj_weights=lj_weights,
                )
                continue
            FEP_BONDED_FORCE_MERGE_RULE[forcename].get("merge_function")(mol_r, mol_a, mol_b,
                                                                         GlobalSetting.BondedForcesMap[forcename],
                                                                         rforces, bforces, temp_lambda, mol_r2mol_a,
                                                                         mol_a2mol_r, mol_r2mol_b, mol_b2mol_r)
        elif rforces or mol_b.bonded_forces.get(forcename):
            raise NotImplementedError(forcename + " is not supported for FEP to merge force field yet.")

    mol_r._fep_materialized_topology = True
    return mol_r


set_global_alternative_names()
