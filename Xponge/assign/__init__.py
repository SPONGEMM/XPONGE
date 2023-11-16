"""
This **package** is used to assign the properties for atoms, residues and molecules
"""
#pylint: disable=cyclic-import
import heapq
import io
import re
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Iterable
from itertools import groupby
import numpy as np
from ..helper import AtomType, ResidueType, Residue, Xopen, Xdict, set_real_global_variable, \
    set_global_alternative_names, Guess_Element_From_Mass, Xprint, \
    get_basis_vectors_from_length_and_angle


__all__ = ["Assign", "get_assignment_from_pdb", "get_assignment_from_mol2", "get_assignment_from_pubchem",
           "get_assignment_from_residuetype", "get_assignment_from_xyz", "get_assignment_from_smiles",
           "get_assignment_from_cif"]


class AssignRule:
    """
    This **class** is to be the rule to determine the atom type for one atom

    :param name: the name of the rule
    """
    all = Xdict(not_found_message="AssignRule {} not found. Did you import the proper force field?")

    def __init__(self, name, pure_string=False, pre_action=None, post_action=None):
        self.name = name
        AssignRule.all[name] = self
        self.rules = OrderedDict()
        self.priority = Xdict()
        self.built = False
        self.pure_string = pure_string
        self.pre_action = pre_action
        self.post_action = post_action

    def add_rule(self, atomtype, priority=0):
        """
        This **function** is used as a **decorator** to add the atom type - judge function

        :param atomtype: a string or an AtomType instance
        :param priority: if more than one judge function returns True,  \
the atom type with higher priority will be chosen. \
If the priority levels of the functions are the same, the atom type which is added first will be chosen.
        :return: a **function**, which wraps a judge function (receiving the Assign instance and the atom index \
and giving True or False as a result)
        """
        if isinstance(atomtype, str):
            if not self.pure_string:
                atomtype = AtomType.get_type(atomtype)
        elif not isinstance(atomtype, AtomType):
            raise TypeError("atomtype should be a string or AtomType")

        def wrapper(rule_function):
            self.rules[atomtype] = rule_function
            self.priority[atomtype] = -priority

        return wrapper

    def set_pre_action(self, function):
        self.pre_action = function

    def set_post_action(self, function):
        self.post_action = function

class _RING():
    """
    This **class** is used to help with the ring assignment.
    """
    def __init__(self, atom_list):
        min_index = np.argmin(atom_list)
        self.atoms = atom_list[min_index:] + atom_list[:min_index]
        reverse_list = self.atoms[::-1]
        reverse_list = reverse_list[-1:] + reverse_list[:-1]
        if reverse_list[1] < self.atoms[1]:
            self.atoms = reverse_list
        self.tohash = "-".join(["%d" % atom for atom in self.atoms])
        self.is_pure_aromatic_ring = None
        self.is_pure_aliphatic_ring = None
        self.is_planar_ring = None
        self.out_plane_double_bond = None
        self.is_aromatic = False

    def __repr__(self):
        return self.tohash

    def __hash__(self):
        return hash(self.tohash)

    def __eq__(self, other):
        return isinstance(other, _RING) and self.tohash == other.tohash

    def __len__(self):
        return len(self.atoms)

    @staticmethod
    def add_rings_basic_marker(assign, rings):
        """

        :param assign:
        :param rings:
        :return:
        """
        for ring in rings:
            for atom in ring.atoms:
                assign.Add_Atom_Marker(atom, "RG")
                assign.Add_Atom_Marker(atom, "RG%d" % len(ring.atoms))

    @staticmethod
    def check_rings_type(assign, rings):
        """

        :param assign:
        :param rings:
        :return:
        """
        for ring in rings:
            ring.check_pure_aromatic(assign)
            ring.check_pure_aliphatic_and_planar(assign)
            ring.check_out_plane_double_bond(assign)

            if not ring.is_pure_aromatic_ring:
                for atom in ring.atoms:
                    if ring.is_pure_aliphatic_ring:
                        assign.Add_Atom_Marker(atom, "AR5")
                    elif ring.is_planar_ring:
                        if ring.out_plane_double_bond:
                            assign.Add_Atom_Marker(atom, "AR3")
                        else:
                            assign.Add_Atom_Marker(atom, "AR2")
                    else:
                        assign.Add_Atom_Marker(atom, "AR4")

    @staticmethod
    def get_rings(assign):
        """

        :param assign:
        :return:
        """
        current_path = []
        current_path_sons = Xdict()
        current_work = []
        have_found_rings = set([])
        for atom0 in range(len(assign.atoms)):
            current_path.append(atom0)
            current_work.extend([[atom, atom0] for atom in assign.bonds[atom0].keys()])
            current_path_sons[atom0] = len(assign.bonds[atom0])
            current_path_father = []
            while current_path:
                work_atom, from_atom = current_work.pop()
                current_path.append(work_atom)
                current_path_father.append(from_atom)
                bond_atom = []
                for atom in assign.bonds[work_atom].keys():
                    if atom != from_atom:
                        try:
                            index = current_path.index(atom)
                            have_found_rings.add(_RING(current_path[index:]))
                        except ValueError:
                            bond_atom.append([atom, work_atom])

                if len(current_path) < 9:
                    current_path_sons[work_atom] = len(bond_atom)
                    current_work.extend(bond_atom)

                else:
                    current_path_sons[work_atom] = 0

                for atom in current_path[::-1]:
                    if current_path_sons[atom] == 0:
                        pop_atom = current_path.pop()
                        current_path_sons.pop(pop_atom)
                        if current_path_father:
                            father_atom = current_path_father.pop()
                            current_path_sons[father_atom] -= 1
        return have_found_rings

    def get_3_neighbors(self):
        """

        :return:
        """
        for i, atom in enumerate(self.atoms):
            yield self.atoms[i - 2], self.atoms[i - 1], atom

    def check_pure_aromatic(self, assign):
        """

        :param assign:
        :return:
        """
        if len(self.atoms) == 6:
            self.is_pure_aromatic_ring = True
            for atom in self.atoms:
                if not assign.Atom_Judge(atom, "C3") and not assign.Atom_Judge(atom, "N2") and not assign.Atom_Judge(
                        atom, "N3"):
                    self.is_pure_aromatic_ring = False
                    break
                if assign.Atom_Judge(atom, "N3"):
                    temp = 0
                    for bonded_atom, bond_order in assign.bonds[atom].items():
                        temp += bond_order
                    if temp == 3:
                        self.is_pure_aromatic_ring = False
                        break
                for bonded_atom, bond_order in assign.bonds[atom].items():
                    if bond_order == 2 and "RG" not in assign.atom_marker[bonded_atom].keys():
                        self.is_pure_aromatic_ring = False
                        break
                if not self.is_pure_aromatic_ring:
                    break
        else:
            self.is_pure_aromatic_ring = False

    def check_pure_aliphatic_and_planar(self, assign):
        """

        :param assign:
        :return:
        """
        self.is_pure_aliphatic_ring = True
        self.is_planar_ring = True
        for atom in self.atoms:
            if self.is_pure_aromatic_ring:
                assign.Add_Atom_Marker(atom, "AR1")
                for i in range(6):
                    assign.Add_Bond_Marker(self.atoms[i - 1], self.atoms[i], "AB")
            if not assign.Atom_Judge(atom, "C4"):
                self.is_pure_aliphatic_ring = False
            if (not assign.Atom_Judge(atom, "C3") and not assign.Atom_Judge(atom, "N2")
                    and not assign.Atom_Judge(atom, "N3") and not assign.Atom_Judge(atom, "O2")
                    and not assign.Atom_Judge(atom, "S2") and not assign.Atom_Judge(atom, "P2")
                    and not assign.Atom_Judge(atom, "P3")):
                self.is_planar_ring = False

    def check_out_plane_double_bond(self, assign):
        """

        :param assign:
        :return:
        """
        self.out_plane_double_bond = False
        for atom in self.atoms:
            for bonded_atom, order in assign.bonds[atom].items():
                if assign.atoms[bonded_atom] != "C" and order == 2 and bonded_atom not in self.atoms:
                    self.out_plane_double_bond = True

    def check_aromatic(self, assign):
        """

        :param assign:
        :return:
        """
        if len(self) < 4:
            self.is_aromatic = False
            return False
        pi_electron = 0
        for atom0, atom1, atom2 in self.get_3_neighbors():
            degree = len(assign.bonds[atom1])
            valence = sum(assign.bonds[atom1].values())
            charge = assign.formal_charge[atom1]
            if assign.atoms[atom1] == "C":
                if charge == 0:
                    if degree == 3:
                        outside_atom = next(iter(set(assign.bonds[atom1]) - {atom0, atom2}))
                        if assign.atoms[outside_atom] == "C" or assign.bonds[atom1][outside_atom] != 2:
                            pi_electron += 1
                    else:
                        break
                elif charge == 1 and valence == 3:
                    if degree == 2:
                        pi_electron += 1
                    elif degree != 3:
                        break
                elif charge == -1 and valence == 3:
                    if degree == 3:
                        pi_electron += 2
                    elif degree == 2:
                        pi_electron += 1
                    else:
                        break
                else:
                    break
            elif assign.atoms[atom1] in ("P", "N"):
                if charge == 0:
                    if valence == 3:
                        if degree == 3:
                            pi_electron += 2
                        elif degree == 2:
                            pi_electron += 1
                    elif valence == 5:
                        outside_atom = next(iter(set(assign.bonds[atom1]) - {atom0, atom2}))
                        if assign.atoms[outside_atom] == "O":
                            pi_electron += 1
                        else:
                            pi_electron += 2
                    else:
                        break
                elif charge == 1 and valence == 4 and degree == 3:
                    pi_electron += 1
                elif charge == -1 and valence == 2 and degree == 2:
                    pi_electron += 2
                else:
                    break
            elif assign.atoms[atom1] == "O":
                if charge == 0 and valence == 2 and degree == 2:
                    pi_electron += 2
                elif charge == 1 and valence == 3 and degree == 2:
                    pi_electron += 1
                else:
                    break
            elif assign.atoms[atom1] == "S":
                if charge == 0:
                    if valence == 2 and degree == 2:
                        pi_electron += 2
                    else:
                        outside_atom = next(iter(set(assign.bonds[atom1]) - {atom0, atom2}))
                        if  degree == 3 and valence == 4 and assign.atoms[outside_atom] == "O":
                            pi_electron += 2
                        else:
                            break
                elif charge == 1:
                    if valence == 3 and degree == 2:
                        pi_electron += 1
                    else:
                        outside_atom = next(iter(set(assign.bonds[atom1]) - {atom0, atom2}))
                        if  degree == 3 and valence == 3 and assign.atoms[outside_atom] == "O":
                            pi_electron += 2
                        else:
                            break
                else:
                    break
        else:
            if pi_electron % 4 == 2:
                for atom in self.atoms:
                    assign.add_atom_marker(atom, "AR0")
                self.is_aromatic = True
                return True
        self.is_aromatic = False
        return False


class Assign():
    """
    This **class** is used to assign properties for atoms, which is called an "assignment"

    :param name: the name of the molecule
    """
    XX = set("CNOPS")
    XA = set("OS")
    XB = set("NP")
    XC = set(["F", "Cl", "Br", "I"])
    XD = set("SP")
    XE = set(["N", "O", "F", "Cl", "Br", "S", "I"])

    CONNECTIVITY_RADII = {"H": 0.35, "C": 0.73, "N": 0.66, "O": 0.69, "F": 0.68,
                          "P": 1.04, "S": 0.96, "Cl": 0.95, "Br": 1.08, "I": 1.26}

    def __init__(self, name="ASN"):
        self.name = name
        self.atom_numbers = 0
        self.atoms = []
        self.names = []
        self.formal_charge = []
        self.element_details = []
        self.coordinate = None
        self.charge = None
        self._built = False
        self.kekulized = False
        self.atom_types = Xdict()
        self.atom_marker = Xdict()
        self.bonds = Xdict()
        self.bond_marker = Xdict()

    @property
    def built(self):
        return self._built

    @built.setter
    def built(self, value):
        if value is True:
            self._built = True
        else:
            self._built = False
            self.kekulized = False

    def add_index_to_name(self):
        """
        This **function** renames the atoms by adding the index to the element name

        :return: None
        """
        count = Xdict()
        for i in range(self.atom_numbers):
            atom_name = self.names[i] if self.names[i] else self.atoms[i]
            if atom_name not in count:
                count[atom_name] = 0
            else:
                count[atom_name] += 1
                atom_name += f"{count[atom_name]}"
            self.names[i] = atom_name

    def atom_judge(self, atom, string):
        """
        This **function** helps judge whether the atom belongs to the mask. For example, "O2" means an oxygen atom \
connected to two other atoms, "N4" means a nitrogen atom connected to four other atoms.

        :param atom: the index of the atom
        :param string: a string mask  of a list of string masks.
        :return:
        """
        assert isinstance(string, (Iterable, str))
        if isinstance(string, str):
            todo = [string]
        else:
            todo = string
        judge = False
        for s in todo:
            element, links = [''.join(list(g)) for k, g in groupby(s, key=lambda x: x.isdigit())]
            if self.atoms[atom] == element and len(self.bonds[atom]) == int(links):
                judge = True
                break
        return judge

    def add_atom(self, element, x, y, z, name="", charge=0.0):
        """
        This **function** adds an atom to the Assign

        :param element: the chemical symbol for the element. "O" - oxygen, "H" - hydrogen for example.
        :param x: the x coordinate
        :param y: the y coordinate
        :param z: the z coordinate
        :param name: the name of the atom
        :param charge: the charge of the atom
        :return: None
        """
        if "." in element:
            element, element_detail = element.split(".")
            element_detail = "." + element_detail
        else:
            element_detail = ""
        self.built = False
        self.element_details.append(element_detail)
        self.atoms.append(element)
        self.bonds[self.atom_numbers] = Xdict()
        self.bond_marker[self.atom_numbers] = Xdict()
        self.atom_marker[self.atom_numbers] = Xdict()
        self.atom_types[self.atom_numbers] = None
        self.atom_numbers += 1
        self.names.append(name)
        if self.coordinate is None:
            self.coordinate = np.array([[float(x), float(y), float(z)]])
        else:
            self.coordinate = np.vstack((self.coordinate, np.array([x, y, z])))
        if self.charge is None:
            self.charge = np.array([charge])
        else:
            self.charge = np.hstack((self.charge, np.array([charge])))
        self.formal_charge.append(0)

    def add_atom_marker(self, atom, marker):
        """
        This **function** adds a marker to an atom

        :param atom: the atom index
        :param marker: the marker
        :return: None
        """
        if marker in self.atom_marker[atom].keys():
            self.atom_marker[atom][marker] += 1
        else:
            self.atom_marker[atom][marker] = 1

    def add_bond(self, atom1, atom2, order=-1):
        """
        This **function** adds a bond to two atoms

        :param atom1: the index of the first atom
        :param atom2: the index of the the second atom
        :param order: the bond order
        :return: None
        """
        self.built = False
        self.bonds[atom1][atom2] = order
        self.bond_marker[atom1][atom2] = set([])
        self.bonds[atom2][atom1] = order
        self.bond_marker[atom2][atom1] = set([])

    def add_bond_marker(self, atom1, atom2, marker, only1=False):
        """
        This **function** adds a marker to a bond

        :param atom1: the index of the first atom
        :param atom2: the index of the the second atom
        :param marker: the marker
        :param only1: only add the marker to the atom1 - atom2 bond instead of the atom2 - atom1 bond
        :return: None
        """
        self.bond_marker[atom1][atom2].add(marker)
        if marker in self.atom_marker[atom1]:
            self.atom_marker[atom1][marker] += 1
        else:
            self.atom_marker[atom1][marker] = 1
        if not only1:
            self.bond_marker[atom2][atom1].add(marker)
            if marker in self.atom_marker[atom2]:
                self.atom_marker[atom2][marker] += 1
            else:
                self.atom_marker[atom2][marker] = 1

    def delete_atom(self, atom):
        """
        This **function** deletes the atom

        :param atom: the index of the atom to delete
        :return: None
        """
        if atom < 0 or atom >= self.atom_numbers:
            raise ValueError(f"the index of the atom to delete should be in the range of\
0 ~ {self.atom_numbers}, but {atom} got")
        self.atom_numbers -= 1
        self.atoms.pop(atom)
        self.names.pop(atom)
        self.formal_charge.pop(atom)
        self.element_details.pop(atom)
        self.coordinate = np.delete(self.coordinate, atom, 0)
        self.charge = np.delete(self.charge, atom, 0)
        self.built = False
        #move forward keys of Xdict
        for btom in range(atom + 1, self.atom_numbers + 1):
            self.bonds[btom - 1] = self.bonds.pop(btom)
            self.atom_types[btom - 1] = self.atom_types.pop(btom)
        for btom, bond in self.bonds.items():
            self.bonds[btom] = Xdict({key - 1 if key > btom else key : value
                for key, value in bond.items() if key != btom})

    def delete_bond(self, atom1, atom2):
        """
        This **function** deletes the bond between two atoms

        :param atom1: the index of the first atom
        :param atom2: the index of the the second atom
        :return: None
        """
        self.built = False
        self.bonds[atom1].pop(atom2, None)
        self.bonds[atom2].pop(atom1, None)

    def determine_equal_atoms(self):
        """
        This **function** dertermines the chemical equalvalent atoms

        .. NOTE::

            The pyckage RDKit is needed for this **function**

        :return: a list of equalvalent atom index lists
        """
        from ..helper.rdkit import Find_Equal_Atoms
        return Find_Equal_Atoms(self)

    def determine_ring_and_bond_type(self):
        """
        This **function** determine the ring and the bond type

        :return: None
        """
        if not self.check_connectivity():
            Xprint("The atoms in the assignment are not all in a connected graph", "ERROR")
        for atom in range(len(self.atoms)):
            self.atom_marker[atom].clear()
            dlo = 0
            noto = 0
            for atom2, order in self.bonds[atom].items():
                self.bond_marker[atom][atom2].clear()
                if self.Atom_Judge(atom2, "O1"):
                    dlo += 1
                else:
                    noto += 1
            if dlo >= 1 >= noto:
                for atom2, order in self.bonds[atom].items():
                    if self.Atom_Judge(atom2, "O1"):
                        self.Add_Bond_Marker(atom, atom2, "DLB")
            for atom2, order in self.bonds[atom].items():
                if "DLB" in self.bond_marker[atom][atom2]:
                    self.Add_Bond_Marker(atom, atom2, "DL", True)
                    self.Add_Bond_Marker(atom, atom2, "sb", True)
                elif order == 1:
                    self.Add_Bond_Marker(atom, atom2, "sb", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "SB", True)
                elif order == 2:
                    self.Add_Bond_Marker(atom, atom2, "db", True)
                    if "AB" not in self.bond_marker[atom][atom2]:
                        self.Add_Bond_Marker(atom, atom2, "DB", True)
                else:
                    self.Add_Bond_Marker(atom, atom2, "tb", True)
        self.rings = _RING.get_rings(self)
        _RING.add_rings_basic_marker(self, self.rings)
        _RING.check_rings_type(self, self.rings)
        self.built = True

    def determine_atom_type(self, rule):
        """
        This **function** determines the atom type.

        :param rule: a string or an AssignRule instance
        :return: if the attribute "pure_string" of the rule is False, the atom types will be saved inplace and return
None, else return the atom types.
        """
        if not self.built:
            self.determine_ring_and_bond_type()
        if isinstance(rule, str):
            rule = AssignRule.all[rule]
        if not rule.built:
            rule.rules = OrderedDict(sorted(rule.rules.items(),
                                            key=lambda x: rule.priority[x[0]]))
            rule.built = True
        if rule.pure_string:
            backup = deepcopy(self.atom_types)
        if rule.pre_action:
            rule.pre_action(self)
        for i in range(len(self.atoms)):
            for atom_type, type_rule in rule.rules.items():
                if type_rule(i, self):
                    self.atom_types[i] = atom_type
                    break
            else:
                raise KeyError("No atom type found for assignment %s of atom #%d" % (self.name, i))
        if rule.post_action:
            rule.post_action(self)
        if rule.pure_string:
            backup, self.atom_types = self.atom_types, backup
            return backup
        return None

    def determine_connectivity(self, simple_cutoff=None, tolerance=1.00):
        """
        This **function** determines the connectivity based on atomic distances

        :param simple_cutoff: the distance cutoff to determine whether the two atoms are connected. If None (default),
the rule described in the reference (J. Wang et al., J. Am. Chem. Soc, 2001) will be used.
        :param tolerance: the tolerance factor for the default method
        :return: None
        """
        if simple_cutoff is None:
            for i, ci in enumerate(self.coordinate):
                for j in range(i + 1, self.atom_numbers):
                    dij = np.linalg.norm(np.array(self.coordinate[j]) - np.array(ci))
                    rij = self.CONNECTIVITY_RADII.get(self.atoms[i], 1.25) + \
                         self.CONNECTIVITY_RADII.get(self.atoms[j], 1.25)
                    if dij <= 1.5:
                        factor = 1 - 0.15
                    elif dij <= 1.9:
                        factor = 1 - 0.11
                    elif dij <= 2.05:
                        factor = 1 - 0.09
                    else:
                        factor = 1 - 0.08
                    factor /= tolerance
                    if factor < rij / dij < 2:
                        self.add_bond(i, j, -1)
        else:
            for i, ci in enumerate(self.coordinate):
                for j in range(i + 1, self.atom_numbers):
                    dij = np.linalg.norm(np.array(self.coordinate[j]) - np.array(ci))
                    if dij < simple_cutoff:
                        self.add_bond(i, j, -1)

    def determine_bond_order(self, max_step=2000, max_stat=20000, penalty_scores=None,
                             total_charge=None, extra_criteria=None):
        """
        This **function** determines the bond order based on connectivities

        :param max_step: the max iterative step
        :param max_stat: the max iterative stat
        :param penalty_scores: the penalty scores for every atom. \
This should be a list of ordered dicts, and every ordered dict stores the valence-penalty pairs for every atom, \
and it is sorted by the penalty scores. If None(default), \
a set of penalty scores described in the reference (J. Wang et al., J. Mol. Graph. Model., 2006) will be used.
        :param total_charge: the total charge of the molecule
        :param extra_criteria: a function as the extra convergence criteria. \
The function will receive the assignment as input, and give True or False as output.
        :return: ReasonedBool, True for success, False for failure.
        """
        from .bond_order import BondOrderAssignment
        bo_assign = BondOrderAssignment(penalty_scores, max_step, max_stat, self, total_charge, extra_criteria)
        return bo_assign.main()

    def check_connectivity(self):
        """
            This **function** checks whether all atoms are connected in one graph

            :return: True or False
        """
        visited = set()
        stack = [next(iter(self.bonds))]
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                stack.extend(self.bonds.get(node, {}))
        return len(visited) == self.atom_numbers


    def to_residuetype(self, name, charge=None):
        """
        This **function** converts the Assign instance to the ResidueType instance

        :param name: the name of the ResidueType instance
        :param charge: the charge of atoms. If set to None, internal charge will be used
        :return: the ResidueType instance
        """
        temp = ResidueType(name=name)
        if not charge:
            if self.charge is None:
                charge = np.zeros(self.atom_numbers)
            else:
                charge = self.charge
        self.add_index_to_name()
        for i in range(self.atom_numbers):
            temp.Add_Atom(self.names[i], self.atom_types[i], x=self.coordinate[i][0],
                          y=self.coordinate[i][1], z=self.coordinate[i][2])
            temp.atoms[-1].charge = charge[i]
        for i, bondi in self.bonds.items():
            for j in bondi.keys():
                temp.Add_Connectivity(temp.atoms[i], temp.atoms[j])
        set_real_global_variable(name, temp)
        return temp

    def calculate_charge(self, method, **parameters):
        """
        This **function** calculates the partial charge for every atom.
        the method "RESP" to calculate the partial charge is not available on Windows.

        :param method: the method to calculate the charge
        :param parameters: the parameters to calculate the charge
        :return: None
        """
        method = method.upper()
        if method == "RESP":
            from . import resp
            self.charge = resp.RESP_Fit(self, basis=parameters.get("basis", "6-31g*"), opt=parameters.get("opt", False),
                                        charge=parameters.get("charge", int(round(sum(self.formal_charge)))),
                                        spin=parameters.get("spin", 0),
                                        extra_equivalence=parameters.get("extra_equivalence", []),
                                        grid_density=parameters.get("grid_density", 6),
                                        grid_cell_layer=parameters.get("grid_cell_layer", 4),
                                        a1=parameters.get("a1", 0.0005),
                                        a2=parameters.get("a2", 0.001), two_stage=parameters.get("two_stage", True),
                                        only_esp=parameters.get("only_esp", False),
                                        radius=parameters.get("radius", None))
        elif method == "GASTEIGER":
            from . import gasteiger
            self.charge = gasteiger.Gasteiger(self)
        elif method == "TPACM4":
            from . import tpacm4
            self.charge = tpacm4.tpacm4(self, charge=parameters.get("charge", int(round(sum(self.formal_charge)))))
        else:
            raise ValueError("methods should be one of the following: 'RESP', 'GASTEIGER', 'TPACM4' (case-insensitive)")

    def kekulize(self):
        """
        This **function** kekulizes the structure.
        The marker "AR0" will be added to the aromatic atom and \
the marker "ar" will be added to the aromatic bond.

        :return: None
        """
        if not self.built:
            self.determine_ring_and_bond_type()
        if not self.kekulized:
            for ring in self.rings:
                if ring.check_aromatic(self):
                    for atom0, atom1, _ in ring.get_3_neighbors():
                        self.add_bond_marker(atom0, atom1, "ar")
            self.kekulized = True

    def uff_optimize(self):
        """
        This **function** uses rdkit and uff to optimize the structure

        :return None:
        """
        from rdkit.Chem import AllChem
        from ..helper.rdkit import assign_to_rdmol
        rdmol = assign_to_rdmol(self)
        AllChem.UFFOptimizeMolecule(rdmol)
        for i in range(self.atom_numbers):
            self.coordinate[i] = rdmol.GetConformer().GetAtomPosition(i)

    def save_as_pdb(self, filename):
        """
        This **function** saves the instance as a pdb file

        :param filename: the name of the output file
        :return: None
        """
        if not isinstance(filename, str):
            raise TypeError("filename needed to save an assignment to a pdb file")
        towrite = towrite = "REMARK   Generated By Xponge (Assignment)\n"
        count = Xdict()
        for i in range(self.atom_numbers):
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d" % count[self.atoms[i]]
                self.names[i] = atom_name
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
                self.names[i] = atom_name
        for i, atom in enumerate(self.atoms):
            towrite += "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%22s%2s\n" % (i + 1,
                                                                                 self.names[i],
                                                                                 self.name,
                                                                                 " ", 1,
                                                                                 self.coordinate[i][0],
                                                                                 self.coordinate[i][1],
                                                                                 self.coordinate[i][2], " ", atom)

        for i in range(self.atom_numbers):
            bonded_atoms = list(self.bonds[i].keys())
            bonded_atoms.sort()
            bonded_atoms = [bonded_atoms[i:i + 4] for i in range(0, len(bonded_atoms), 4)]
            if bonded_atoms:
                for atoms in bonded_atoms:
                    towrite += "CONECT %4d" % (i + 1)
                    for atom in atoms:
                        towrite += " %4d" % (atom + 1)
                    towrite += "\n"

        f = Xopen(filename, "w")
        f.write(towrite)
        f.close()

    def save_as_mol2(self, filename):
        """
        This **function** saves the instance as a mol2 file

        :param filename: the name of the output file
        :return: None
        """
        if not isinstance(filename, str):
            raise TypeError("filename needed to save an assignment as a mol2 file")
        import Xponge.forcefield.sybyl #pylint: disable=unused-import
        atom_types = self.determine_atom_type("sybyl")
        for atom, atype in enumerate(atom_types):
            self.element_details[atom] = "." + atype.split(".")[1] if "." in atype else ""
        bonds = []
        for i in range(self.atom_numbers):
            for j, order in self.bonds[i].items():
                if i < j:
                    if "ar" in self.bond_marker[i][j]:
                        bonds.append("%6d %6d ar\n" % (i + 1, j + 1))
                    elif "am" in self.bond_marker[i][j]:
                        bonds.append("%6d %6d am\n" % (i + 1, j + 1))
                    elif order == -1:
                        bonds.append("%6d %6d un\n" % (i + 1, j + 1))
                    else:
                        bonds.append("%6d %6d %1d\n" % (i + 1, j + 1, order))
        bonds.sort(key=lambda x: list(map(int, x.split()[:2])))
        count = Xdict()
        for i in range(self.atom_numbers):
            if self.names[i]:
                atom_name = self.names[i]
            elif self.atoms[i] in count.keys():
                atom_name = self.atoms[i] + "%d" % count[self.atoms[i]]
                self.names[i] = atom_name
                count[self.atoms[i]] += 1
            else:
                count[self.atoms[i]] = 1
                atom_name = self.atoms[i]
                self.names[i] = atom_name
        towrite = "@<TRIPOS>MOLECULE\n%s\n %d %d 1 0 1\nSMALL\nUSER_CHARGES\n" % (
            self.name, self.atom_numbers, len(bonds))
        towrite += "@<TRIPOS>ATOM\n"
        for i, atom in enumerate(self.atoms):
            towrite += "%6d %4s %8.4f %8.4f %8.4f   %-8s %5d %8s %10.6f\n" % (
                i + 1, self.names[i], self.coordinate[i][0], self.coordinate[i][1], self.coordinate[i][2],
                atom + self.element_details[i], 1, self.name, self.charge[i])
        charged_atom = list(filter(lambda i: i[1] != 0, enumerate(self.formal_charge)))
        if charged_atom:
            towrite += "@<TRIPOS>UNITY_ATOM_ATTR\n"
        for i, charge in charged_atom:
            towrite += f"{i + 1} 1\ncharge {charge}\n"
        towrite += "@<TRIPOS>BOND\n"
        for i, bond in enumerate(bonds):
            towrite += "%6d %s" % (i + 1, bond)
        towrite += "@<TRIPOS>SUBSTRUCTURE\n"
        towrite += "%5d %8s %6d ****               0 ****  **** \n" % (1, self.name, 1)

        f = Xopen(filename, "w")
        f.write(towrite)
        f.close()

    def set_ph(self, ph):
        """
            This **function** sets the pH value, and adds or deletes the related hydrogens.

            :param ph: the pH value
            :return: the sum of final formal charge
        """
        from .phmodel import PHModelAssignment
        phma = PHModelAssignment(self, ph)
        return phma.main()

def get_assignment_from_pubchem(parameter, keyword):
    """
    This **function** gets an Assign instance from PubChem

    usage example::

        a1 = Get_Assignment_From_PubChem("ethane", "name")
        a2 = Get_Assignment_From_PubChem("CC", "smiles")

    :param parameter: the parameter to search on PubChem
    :param keyword: the keyword to search on PubChem
    :return: the Assign instance
    """
    import pubchempy as pcp
    cs = pcp.get_compounds(parameter, keyword, record_type='3d')
    if not cs:
        raise pcp.NotFoundError
    if len(cs) == 1:
        assign = Assign()
        c = cs[0]
        for atom in c.atoms:
            assign.Add_Atom(atom.element, atom.x, atom.y, atom.z)
        for bond in c.bonds:
            assign.Add_Bond(bond.aid1 - 1, bond.aid2 - 1, bond.order)
        assign.Determine_Ring_And_Bond_Type()
        return assign
    raise NotImplementedError


def get_assignment_from_smiles(smiles):
    """
    This **function** gets an Assign instance from smiles

    usage example::

        a1 = get_assignment_from_smiles("CC", "smiles")

    :param smiles: the smiles to get
    :return: the Assign instance
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from ..helper.rdkit import rdmol_to_assign
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return rdmol_to_assign(mol)


def get_assignment_from_pdb(file, only_residue="", bond_tolerance=1.0, total_charge=None):
    """
    This **function** gets an Assign instance from a pdb file

    :param file: the name of the input file or an instance of io.IOBase
    :param only_residue: only get the residue with the name same as ``only_residue``
    :param bond_tolerance: the parameter to determine the atomic connections. \
The larger tolerance, the easier to set a bond between two atoms
    :param total_charge: the total charge of the molecule used when aligning bond orders. \
If None is given, the total charge will not be checked
    :return: the Assign instance
    """
    assign = Assign()
    index_atom_map = Xdict()
    has_conect = False
    if not isinstance(file, io.IOBase):
        filename = file
        file = open(file)
    else:
        filename = "in-memory string"
    with file as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                if only_residue:
                    if resname != only_residue:
                        continue
                assign.name = resname
                index = int(line[6:11])
                index_atom_map[index] = assign.atom_numbers
                atom_name = line[12:16].strip()
                element = line[76:78].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                assign.Add_Atom(element, x, y, z, atom_name)
            if line.startswith("CONECT"):
                has_conect = True
                atom = int(line[6:11])
                if atom not in index_atom_map.keys():
                    continue
                for bonded_atom_i in range(11, 31, 5):
                    try:
                        temp = line[bonded_atom_i:bonded_atom_i + 5]
                        bonded_atom = int(temp)
                    except ValueError:
                        break
                    if bonded_atom in index_atom_map.keys():
                        assign.Add_Bond(index_atom_map[atom], index_atom_map[int(bonded_atom)])
    if assign.atom_numbers == 0:
        raise OSError(f"The file {filename} is not a pdb file")
    if not has_conect:
        assign.determine_connectivity(tolerance=bond_tolerance)
    success = assign.Determine_Bond_Order(total_charge=total_charge)
    if not success:
        Xprint(f"The connectivity or the bond orders in {filename} are not reasonable", "WARNING")
    return assign


def get_assignment_from_residuetype(restype):
    """
    This **function** gets an Assign instance from a ResidueType instance

    :param restype: the ResidueType instance
    :return: the Assign instance
    """
    if not isinstance(restype, (ResidueType, Residue)):
        raise OSError(f"{restype} is not a ResidueType instance")
    assign = Assign()
    for atom in restype.atoms:
        assign.Add_Atom(Guess_Element_From_Mass(atom.mass),
                        atom.x, atom.y, atom.z,
                        atom.name, atom.charge)
    for atom in restype.atoms:
        i = restype.atom2index(atom)
        for atomb in restype.connectivity[atom]:
            j = restype.atom2index(atomb)
            if i < j:
                assign.Add_Bond(i, j)
    total_charge = int(round(sum(assign.charge)))
    success = assign.Determine_Bond_Order(total_charge=total_charge)
    if not success:
        Xprint(f"The connectivity, the bond orders or the charges of the ResidueType {restype.name} are not reasonable",
            "WARNING")
    return assign


def get_assignment_from_xyz(file, bond_tolerance=1.0, total_charge=None):
    """
    This **function** gets an Assign instance from a xyz file

    :param file: the name of the input file or an instance of io.IOBase
    :param bond_tolerance: the parameter to determine the atomic connections. \
The larger tolerance, the easier to set a bond between two atoms
    :param total_charge: the total charge of the molecule used when aligning bond orders. \
If None is given, the total charge will not be checked
    :return: the Assign instance
    """
    assign = None
    if not isinstance(file, io.IOBase):
        filename = file
        file = open(file)
    else:
        filename = "in-memory string"
    with file as f:
        atom_numbers = int(f.readline())
        assign = Assign()
        assign.name = f.readline().strip()
        for _ in range(atom_numbers):
            atom_name, x, y, z = f.readline().split()
            assign.Add_Atom(atom_name, float(x), float(y), float(z))
        assign.determine_connectivity(tolerance=bond_tolerance)
        success = assign.Determine_Bond_Order(total_charge=total_charge)
        if not success:
            Xprint(f"The connectivity or the bond orders in {filename} are not reasonable", "WARNING")
    if assign is None:
        raise OSError(f"The file {filename} is not a xyz file")
    return assign


def get_assignment_from_mol2(file, total_charge=None):
    """
    This **function** gets an Assign instance from a mol2 file

    :param file: the name of the input file or an instance of io.IOBase
    :param total_charge: the total charge of the molecule used when aligning bond orders. \
If "sum" is given, the sum of the partial charges will be used; \
If None is given, the total charge will not be checked
    :return: the Assign instance
    """
    if not isinstance(file, io.IOBase):
        filename = file
        file = open(file)
    else:
        filename = "in-memory string"
    with file as f:
        assign = None
        flag = None
        subflag = None
        real_index = Xdict(not_found_message="Atom #{} not found")
        for line in f:
            if not line.strip():
                continue
            if line.startswith("@<TRIPOS>"):
                flag = line[9:].strip()
            elif flag == "MOLECULE":
                if subflag is None:
                    assign = Assign(line.strip())
                    subflag = "0"
            # 处理原子信息
            elif flag == "ATOM":
                words = line.split()
                atom_name = words[1]
                real_index[words[0]] = assign.atom_numbers
                element = words[5]
                x = float(words[2])
                y = float(words[3])
                z = float(words[4])
                charge = float(words[8])
                assign.Add_Atom(element, x, y, z, atom_name, charge)
            elif flag == "UNITY_ATOM_ATTR":
                words = line.split()
                atom = real_index[words[0]]
                for _ in range(int(words[1])):
                    line = f.readline()
                    words = line.split()
                    attr = words[0]
                    if attr == "charge":
                        assign.formal_charge[atom] = int(words[1])
                    else:
                        raise NotImplementedError("Unknown UNITY_ATOM_ATTR flag")
            elif flag == "BOND":
                words = line.split()
                if words[1] not in real_index or words[2] not in real_index:
                    continue
                atom1 = real_index[words[1]]
                atom2 = real_index[words[2]]
                if words[3] in "123456789":
                    assign.Add_Bond(atom1 , atom2, int(words[3]))
                elif words[3] == "ar":
                    assign.Add_Bond(atom1, atom2, -1)
                    assign.add_bond_marker(atom1, atom2, "mol2_ar")
                elif words[3] == "am":
                    assign.Add_Bond(atom1, atom2, -1)
                    assign.add_bond_marker(atom1, atom2, "mol2_am")
                elif words[3] == "un":
                    assign.Add_Bond(atom1, atom2, -1)
                else:
                    raise NotImplementedError(f"No implemented method to process bond #{words[0]} type {words[3]}")
    if total_charge == "sum":
        total_charge = int(round(sum(assign.charge)))
    count_h = assign.atoms.count("H")
    if count_h < len(assign.atoms) // 4:
        Xprint(f"The number of hydrogen atoms in {filename} is {count_h}, \
which is less than a quarter of the number of atoms. \
Xponge.Assign is designed for molecules with explicit hydrogen atoms.", "WARNING")
    success = assign.Determine_Bond_Order(total_charge=total_charge)
    if not success:
        for bond in assign.bonds.values():
            for j in bond:
                bond[j] = -1
        success = assign.Determine_Bond_Order(total_charge=total_charge)
        if not success:
            Xprint(f"The bond orders in {filename} are not reasonable", "WARNING")
        else:
            Xprint(f"The bond orders in {filename} are not reasonable and have been modified", "WARNING")
    if assign is None:
        raise OSError(f"The file {filename} is not a mol2 file")

    sum_of_formal_charge = sum(assign.formal_charge)
    sum_of_partial_charge = int(round(sum(assign.charge)))
    if  sum_of_formal_charge != sum_of_partial_charge:
        Xprint(f"For {filename}, the sum of formal charges ({sum_of_formal_charge}) \
!= the sum of partial charges ({sum_of_partial_charge})", "WARNING")
    return assign

def _cif_find_box_information(key, contents, filename):
    """
        Read one line of the box information from the CIF file
    """
    pattern = "_" + key + r"\s+(\d+\.\d+)"
    match = re.search(pattern, contents)
    if not match:
        raise ValueError(f"There is no {key} found in {filename}")
    return float(match.group(1))

def _cif_get_loop_with_keyword(key, contents):
    """
        Read the loop_ block with the provided keyword
    """
    if key not in contents:
        return None
    pattern = r"^loop_\s*\n((_[^\n]*\n)*_" + key + r"\s*\n(_[^\n]*\n)*)(([^_][^\n]*\n)*)"
    match = re.search(pattern, contents, flags=re.MULTILINE | re.DOTALL)
    keys = match.group(1).split()
    values = Xdict({key : [] for key in keys})
    valuelines = match.group(4).replace("loop_", "").strip().split("\n")
    for value in valuelines:
        for i, v in enumerate(value.split()):
            values[keys[i]].append(v)
    return len(valuelines), values

def _get_cif_float(string, hint, filename):
    """
        convert a string in cif to float
    """
    string = string.strip()
    if string in (".", "?"):
        raise ValueError(f"{hint} in {filename} is {string}, which is not right")
    if "(" in string:
        return float(string.split("(")[0])
    return float(string)

def _parse_cif_symops(symops, lattice_info):
    """
        parse the symmetry operations in the file as basis_position
    """
    if set(symops) - set("+-,xyz0123456789\n /"):
        raise ValueError("the symmetry operator strings can only be simple math expression of x, y, z")
    symops = symops.replace("x", "1").replace("y", "1").replace("z", "1")
    lattice_info["basis_position"] = [[eval(op) for op in symop.split(",")] #pylint:disable=eval-used
                                        for symop in symops.split("\n")]

def get_assignment_from_cif(file, total_charge=0):
    """
    This **function** gets an Assign instance and a preprocessed lattice information from a cif file

    :param file: the name of the input file or an instance of io.IOBase
    :param total_charge: the total charge of the molecule used when aligning bond orders. 0 for default.
    :return: the Assign instance and a dict which stores the preprocessed lattice information
    """
    assign = None
    lattice_info = Xdict({"scale": 1, "style": "custom"})
    if not isinstance(file, io.IOBase):
        filename = file
        file = open(file)
    else:
        filename = "in-memory string"
    with file as f:
        contents = f.read()
        matches = re.findall(r"^data_.+$", contents, flags=re.MULTILINE)
        if len(matches) == 0:
            raise ValueError(f"There is no data block found in {filename}")
        if len(matches) > 1:
            raise NotImplementedError("The support for CIF files with more than one data block is not implemented now")
        assign = Assign(name=matches[0][5:])
        symops = re.findall(r"(_symmetry_equiv_pos_as_xyz|_space_group_symop_operation_xyz)\s+(.+?)(?!\_)\n_\S+",
            contents, flags=re.DOTALL)
        if symops:
            _parse_cif_symops(symops[0][1], lattice_info)
        la = _cif_find_box_information("cell_length_a", contents, filename)
        lb = _cif_find_box_information("cell_length_b", contents, filename)
        lc = _cif_find_box_information("cell_length_c", contents, filename)
        lattice_info["cell_length"] = [la, lb, lc]
        alpha = _cif_find_box_information("cell_angle_alpha", contents, filename)
        beta = _cif_find_box_information("cell_angle_beta", contents, filename)
        gamma = _cif_find_box_information("cell_angle_gamma", contents, filename)
        lattice_info["cell_angle"] = [alpha, beta, gamma]
        basis = get_basis_vectors_from_length_and_angle(la, lb, lc, alpha, beta, gamma)
        n_atom, atom_info = _cif_get_loop_with_keyword("atom_site_type_symbol", contents)
        if not atom_info:
            raise ValueError(f"There is no atom information found in {filename}")
        name_map = Xdict(not_found_message="No atom named {}, but found in bond information")
        for i in range(n_atom):
            if "_atom_site_type_symbol" not in atom_info:
                raise ValueError(f"There is no atom element information found in {filename}")
            element = atom_info["_atom_site_type_symbol"][i]
            if "_atom_site_label" not in atom_info:
                raise ValueError(f"There is no atom name information found in {filename}")
            name = atom_info["_atom_site_label"][i]
            name_map[name] = i
            x = None
            y = None
            z = None
            if "_atom_site_Cartn_x" in atom_info:
                if "_atom_site_Cartn_y" not in atom_info or "_atom_site_Cartn_z" not in atom_info:
                    raise ValueError(f"There is no correct atom position found in {filename}")
                x = _get_cif_float(atom_info["_atom_site_Cartn_x"][i], "_atom_site_Cartn_x", filename)
                y = _get_cif_float(atom_info["_atom_site_Cartn_y"][i], "_atom_site_Cartn_y", filename)
                z = _get_cif_float(atom_info["_atom_site_Cartn_z"][i], "_atom_site_Cartn_z", filename)
            if "_atom_site_fract_x" in atom_info:
                if "_atom_site_fract_y" not in atom_info or "_atom_site_fract_z" not in atom_info:
                    raise ValueError(f"There is no correct atom position found in {filename}")
                fx = _get_cif_float(atom_info["_atom_site_fract_x"][i], "_atom_site_fract_x", filename)
                fy = _get_cif_float(atom_info["_atom_site_fract_y"][i], "_atom_site_fract_y", filename)
                fz = _get_cif_float(atom_info["_atom_site_fract_z"][i], "_atom_site_fract_z", filename)
                x = fx * basis[0][0] + fy * basis[1][0] + fz * basis[2][0]
                y = fx * basis[0][1] + fy * basis[1][1] + fz * basis[2][1]
                z = fx * basis[0][2] + fy * basis[1][2] + fz * basis[2][2]
            if x is None or y is None or z is None:
                raise ValueError(f"There is no atom position found in {filename}")
            assign.add_atom(element, x, y, z, name=name)
        n_bond, bond_info = _cif_get_loop_with_keyword("geom_bond_atom_site_label_1", contents)
        if bond_info:
            if "_geom_bond_atom_site_label_2" not in bond_info:
                raise ValueError(f"There is no correct bond information found in {filename}")
            for i in range(n_bond):
                atom1 = name_map[bond_info["_geom_bond_atom_site_label_1"][i]]
                atom2 = name_map[bond_info["_geom_bond_atom_site_label_2"][i]]
                assign.add_bond(atom1, atom2, -1)
    success = assign.Determine_Bond_Order(total_charge=total_charge)
    if not success:
        Xprint(f"The bond orders in {filename} are not reasonable", "WARNING")
    if assign is None:
        raise OSError(f"The file {filename} is not a mol2 file")
    return assign, lattice_info

set_global_alternative_names()
