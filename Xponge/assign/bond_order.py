"""
This **module** helps to assign bond orders
"""
from itertools import product, combinations
# pylint: disable=cyclic-import
from . import AssignRule, Xdict, OrderedDict, set_attribute_alternative_names, deepcopy, np
from ..helper import ReasonedBool
 
bo = AssignRule("bo", pure_string=True)

@bo.add_rule("X")
def _(i, assign):
    return assign.atoms[i] in ("H", "F", "Cl", "Br", "I")


@bo.add_rule("Cn1")
def _(i, assign):
    return assign.Atom_Judge(i, "C1") and [assign.atoms[key] for key in assign.bonds[i].keys()][0] == "N"


@bo.add_rule("Cx1")
def _(i, assign):
    return assign.Atom_Judge(i, "C1")

#The reference defined this type, but I think it is useless
#@bo.add_rule("Co2")
#def _(i, assign):
#    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
#    return assign.Atom_Judge(i, "C3") and bonded_os == 2


@bo.add_rule("C")
def _(i, assign):
    return assign.atoms[i] == "C"


@bo.add_rule("Nnn1")
def _(i, assign):
    return assign.Atom_Judge(i, "N1") and [assign.atoms[key] for key in assign.bonds[i].keys()][0] == "N"


@bo.add_rule("Nx1")
def _(i, assign):
    return assign.Atom_Judge(i, "N1")


@bo.add_rule("Nnn2")
def _(i, assign):
    return assign.Atom_Judge(i, "N2") and "N1" in {assign.atoms[key] + str(len(assign.bonds[key])) for key in assign.bonds[i]}


@bo.add_rule("Nx2")
def _(i, assign):
    return assign.Atom_Judge(i, "N2")


@bo.add_rule("No2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "N3") and bonded_os == 2


@bo.add_rule("No1")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "N3") and bonded_os == 1


@bo.add_rule("Nx3")
def _(i, assign):
    return assign.Atom_Judge(i, "N3")


@bo.add_rule("Nx4")
def _(i, assign):
    return assign.Atom_Judge(i, "N4")


@bo.add_rule("On1")
def _(i, assign):
    return assign.Atom_Judge(i, "O1") and [assign.atoms[key] for key in assign.bonds[i].keys()].count("N") == 1


@bo.add_rule("Ox1")
def _(i, assign):
    return assign.Atom_Judge(i, "O1")


@bo.add_rule("Ox2")
def _(i, assign):
    return assign.Atom_Judge(i, "O2")


@bo.add_rule("Px1")
def _(i, assign):
    return assign.Atom_Judge(i, "P1")


@bo.add_rule("Px2")
def _(i, assign):
    return assign.Atom_Judge(i, "P2")


@bo.add_rule("Px3")
def _(i, assign):
    return assign.Atom_Judge(i, "P3")


@bo.add_rule("Po3")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "P4") and bonded_os == 3



@bo.add_rule("Po2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "P4") and bonded_os == 2


@bo.add_rule("Px4")
def _(i, assign):
    return assign.Atom_Judge(i, "P4")


@bo.add_rule("Sn1")
def _(i, assign):
    return assign.Atom_Judge(i, "S1") and [assign.atoms[key] for key in assign.bonds[i].keys()].count("N") == 1


@bo.add_rule("Sx1")
def _(i, assign):
    return assign.Atom_Judge(i, "S1")


@bo.add_rule("Sx2")
def _(i, assign):
    return assign.Atom_Judge(i, "S2")


@bo.add_rule("Sx3")
def _(i, assign):
    return assign.Atom_Judge(i, "S3")


@bo.add_rule("So3")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "S4") and bonded_os >= 3


@bo.add_rule("So2")
def _(i, assign):
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "S4") and bonded_os == 2


@bo.add_rule("Sx4")
def _(i, assign):
    return assign.Atom_Judge(i, "S4")


class BondOrderAssignment:
    """
    This **class** includes the functions to assign bond orders

    :param original_penalties: the original penalties dict
    :param max_stat: the max valence stats to iterate
    :param assign: the father Assignment instance
    :param debug: whether to print the debug information
    :param total_charge: the total charge of the molecule
    :param extra_criteria: a function as the extra convergence criteria. The function will receive the assignment as input, and give True or False as output.
    """
    atomic_valence = Xdict({
        "X": OrderedDict({1: 0}),
        "Cn1": OrderedDict({3: 0, 4: 1}),
        "Cx1": OrderedDict({3: 1, 4: 0}),
        "C": OrderedDict({2: 64, 3: 32, 4: 0}),
        "Nnn1": OrderedDict({3: 0, 2: 0}),
        "Nx1": OrderedDict({2: 3, 3: 0, 4: 32}),
        "Nnn2": OrderedDict({3: 1, 4: 0}),
        "Nx2": OrderedDict({2: 4, 3: 0, 4: 32}),
        "No2": OrderedDict({3: 64, 4: 32, 5: 0}),
        "No1": OrderedDict({3: 1, 4: 0}),
        "Nx3": OrderedDict({2: 32, 3: 0, 4: 1, 5: 2}),
        "Nx4": OrderedDict({3: 64, 4: 0, 5: 64}),
        "On1": OrderedDict({1: 0, 0: 1}),
        "Ox1": OrderedDict({2: 0, 1: 1, 3: 64}),
        "Ox2": OrderedDict({1: 32, 2: 0, 3:64}),
        "Px1": OrderedDict({2: 2, 3: 0, 4: 32}),
        "Px2": OrderedDict({2: 4, 3: 0, 4: 2}),
        "Px3": OrderedDict({3: 32, 4: 0, 5: 1, 6: 2}),
        "Po2": OrderedDict({5: 32, 6: 0, 7: 32}),
        "Po3": OrderedDict({6: 32, 7: 0}),
        "Px4": OrderedDict({3: 64, 4: 1, 5: 0, 6: 32}),
        "Sn1": OrderedDict({1: 0, 2: 1}),
        "Sx1": OrderedDict({1: 2, 2: 0, 3: 64}),
        "Sx2": OrderedDict({1: 32, 2: 0, 3: 32}),
        "Sx3": OrderedDict({3: 1, 4: 0, 5: 2, 6: 2}),
        "So2": OrderedDict({6: 0, 7: 32}),
        "So3": OrderedDict({6: 32, 7: 0}),
        "Sx4": OrderedDict({4: 4, 5: 2, 6:0}),
    })
    atomic_formal_valence = Xdict({
        "H": {1: 0},
        "Cl": {1: 0},
        "Br": {1: 0},
        "F": {1: 0},
        "I": {1: 0},
        "C": {4: 0, 3: 1, 5: -1},
        "N": {3: 0, 4: 1, 2: -1, 5: 0},
        "O": {1: -1, 2: 0, 3: 1},
        "P": {3: 0, 4: 1, 2: -1, 5: 0, 7: 0},
        "S": {1: -1, 2: 0, 3: 1, 4: 0, 6: 0},
    })
    failure = ReasonedBool(False, "the calculation can not converge")
    def __init__(self, original_penalties, max_step, max_stat, assign, debug=False, total_charge=0, extra_criteria=None):
        self.prepare_success = True
        self.extra_criteria = extra_criteria
        if original_penalties is None:
            try:
                original_penalties = [self.atomic_valence[type_]
                                      for type_ in assign.determine_atom_type("bo").values()]
            except KeyError as e:
                if "No atom type found for assignment" in e.args[0]:
                    n = e.args[0].split("#")[1]
                    self.prepare_success = ReasonedBool(False, f"the valence of atom #{n} can not be found in the default table")
                else:
                    raise e
        if self.prepare_success:
            self.debug = debug
            self.max_step = max_step
            self.max_stat = max_stat
            self.assign = assign
            self.original_penalties = original_penalties
            self.uc = [set(filter(lambda aj: self.assign.bonds[ai][aj] == -1, self.assign.bonds[ai]))
                         for ai in range(self.assign.atom_numbers)]
            self.connected = [sum(filter(lambda x: x > 0, self.assign.bonds[ai].values()))
                         for ai in range(self.assign.atom_numbers)]
            self.valence_best = [next(iter(pi)) - self.connected[i] if self.uc[i] else 0 for i, pi in enumerate(self.original_penalties)]
            self.valence = deepcopy(self.valence_best)
            self.penalties = Xdict(not_found_method=lambda x: [])
            self._get_penalties(original_penalties)
            self.stat_position = 0
            self.current_stat = 0
            self.points = []
            self.cached = {}
            self.total_charge = total_charge
        set_attribute_alternative_names(self)

    def _get_penalties(self, original_penalties):
        """

        :param original_penalties:
        :return:
        """
        for atom, pi in enumerate(original_penalties):
            for valence, penalty in pi.items():
                self.penalties[penalty].append((atom, penalty, valence))

    def _preprocess_penalties(self, n):
        """
        Get the atomic valence combinations when penality == n
        """
        if n in self.cached:
            return self.cached[n]
        if n == 1:
            toret = [[point] for point in self.penalties[n]]
        else:
            toret = []
            have_added = set()
            for i in range(1, n // 2 + 1):
                r1 = self.cached[i]
                r2 = self.cached[n - i]
                for ri, rj in product(r1, r2):
                    if {r0[0] for r0 in ri} & {r0[0] for r0 in rj}:
                        continue
                    rij = ri + rj
                    rij.sort()
                    rij_checkstring = "+".join([f"{r0[0]}-{r0[1]}-{r0[2]}" for r0 in rij])
                    if rij_checkstring not in have_added:
                        have_added.add(rij_checkstring)
                        toret.append(rij)
            toret.extend([[point] for point in self.penalties.get(n, [])])
        self.cached[n] = toret
        return toret

    def _get_next_valence(self):
        """
        Get the next atomic valence combination to search
        """
        if self.current_stat < self.max_stat:
            if self.stat_position < len(self.points):
                self.valence = deepcopy(self.valence_best)
                has_negative_value = False
                for point in self.points[self.stat_position]:
                    self.valence[point[0]] = point[2] - self.connected[point[0]]
                    if self.valence[point[0]] < 0:
                        has_negative_value = True
                        break
                self.stat_position += 1
                if has_negative_value:
                    return True
            else:
                self.current_stat += 1
                self.stat_position = 0
                self.points = self._preprocess_penalties(self.current_stat)
                if self.debug:
                    print("stat=", self.current_stat)
                    print("points=\n", self.points)
                return True
        return False

    def _check_formal_charge(self, bonds):
        """
        check whether the formal charge is right
        """
        success = True
        formal_charge_iter = None
        total_charge = 0
        c3_atoms = []
        for atom, bond in bonds.items():
            valence = sum(bond.values())
            if self.assign.atoms[atom] == "C" and valence == 3:
                c3_atoms.append(atom)
            formal_charge = self.atomic_formal_valence[self.assign.atoms[atom]].get(valence, None)
            if formal_charge is None:
                success = self.failure
                break
            total_charge += formal_charge
            self.assign.formal_charge[atom] = formal_charge
        if self.total_charge is not None and self.total_charge != total_charge:
            delta_charge = total_charge - self.total_charge
            need_to_change = delta_charge // 2
            if delta_charge % 2 == 1:
                success = self.failure
            elif need_to_change == len(c3_atoms):
                for atom in c3_atoms:
                    self.assign.formal_charge[atom] = -1
            elif 0 < need_to_change < len(c3_atoms):
                formal_charge_iter = combinations(c3_atoms, need_to_change)
                success = self.failure
            else:
                success = self.failure
        return success, formal_charge_iter, c3_atoms

    def _assign_bond_order_one_try(self):
        """
        try to assign bond orders
        """
        bonds = deepcopy(self.assign.bonds)
        uc = deepcopy(self.uc)
        valence = deepcopy(self.valence)
        valence_backup = []
        bonds_backup = []
        uc_backup = []
        atom_guessed = set()
        guess_bonds = []
        success = False
        determined = False
        if self.debug:
            print(f"The initial undetermined valence and undetermined conected atoms for every atom:\n{valence}\n{uc}\n\n")
        while not determined and not success:
            index_sort = np.argsort([len(uci) or float("inf") for uci in uc]).tolist()
            no_basic_rule = True
            for i in index_sort:
                if len(uc[i]) == valence[i] and valence[i] > 0:
                    while uc[i]:
                        j = uc[i].pop()
                        bonds[i][j] = 1
                        bonds[j][i] = 1
                        uc[j].remove(i)
                        valence[j] -= 1
                        if self.debug:
                            print(f"try to assign the order of the bond between {i} {j} to {1}")
                    valence[i] = 0
                    no_basic_rule = False
                    if self.debug:
                        print(f"{valence}\n{uc}\n\n")
                elif len(uc[i]) == 1 and valence[i] > 0:
                    j = uc[i].pop()
                    uc[j].remove(i)
                    bonds[j][i] = valence[i]
                    bonds[i][j] = valence[i]
                    valence[j] -= valence[i]
                    valence[i] -= valence[i]
                    no_basic_rule = False
                    if self.debug:
                        print(f"try to assign the order of the bond between {i} {j} to {bonds[j][i]}\n{valence}\n{uc}\n\n")
            success = True
            for i in range(self.assign.atom_numbers):
                if valence[i] != 0 or len(uc[i]) != 0:
                    success = self.failure
                if (len(uc[i]) == 0 and valence[i] != 0) or (valence[i] == 0 and len(uc[i]) != 0) or valence[i] < 0:
                    success = self.failure
                    determined = True
            if not success and determined and guess_bonds:
                determined = False
                i, j = guess_bonds[-1]
                trial = bonds[i][j]
                if trial == 3:
                    guess_bonds.pop()
                    if not guess_bonds:
                        break
                    uc = uc_backup.pop()
                    uc[i].add(j)
                    uc[j].add(i)
                    bonds = bonds_backup.pop()
                    valence = valence_backup.pop()
                    continue
                uc = deepcopy(uc_backup[-1])
                bonds = deepcopy(bonds_backup[-1])
                valence = deepcopy(valence_backup[-1])
                trial += 1
                bonds[i][j] = trial
                bonds[j][i] = trial
                valence[i] -= trial
                valence[j] -= trial
                if self.debug:
                    print(f"guessing the order of the bond between {i} {j} is {trial}\n{valence}\n{uc}\n\n")
            if not success and not determined and no_basic_rule:
                index_sort = np.argsort([len(uci) if i not in atom_guessed else -1 for i, uci in enumerate(uc)]).tolist()
                i = index_sort[-1]
                if i in atom_guessed or not uc[i]:
                    break
                j = uc[i].pop()
                atom_guessed.add(j)
                atom_guessed.add(i)
                uc[j].remove(i)
                guess_bonds.append([i, j])
                uc_backup.append(deepcopy(uc))
                valence_backup.append(deepcopy(valence))
                valence[i] -= 1
                valence[j] -= 1
                bonds_backup.append(deepcopy(bonds))
                bonds[i][j] = 1
                bonds[j][i] = 1
                if self.debug:
                    print(f"guessing the order of the bond between {i} {j} is {1}\n{valence}\n{uc}\n\n")

        return success, bonds

    def _assign_formal_charge_one_try(self, c3_atoms, formal_charge_iter):
        """
        try to assign formal charges
        """
        for atom in c3_atoms:
            self.assign.formal_charge[atom] = 1
        try:
            for atom in next(formal_charge_iter):
                self.assign.formal_charge[atom] = -1
            return True, formal_charge_iter
        except StopIteration:
            return self.failure, None

    def main(self):
        """
        This **function** is the main function to do the bond order assignment

        :return: True for success, False for failure
        """
        if not self.prepare_success:
            return self.prepare_success
        count = 0
        success = self.failure
        formal_charge_iter = None
        while count < self.max_step and self.current_stat < self.max_stat and not success:
            if formal_charge_iter is None:
                success, bonds = self._assign_bond_order_one_try()
                if success:
                    success, formal_charge_iter, c3_atoms = self._check_formal_charge(bonds)
            else:
                success, formal_charge_iter = self._assign_formal_charge_one_try(c3_atoms, formal_charge_iter)
            if success and self.extra_criteria is not None:
                self.assign.bonds, bonds = bonds, self.assign.bonds
                success = self.extra_criteria(self.assign)
                self.assign.bonds, bonds = bonds, self.assign.bonds
            if not success and formal_charge_iter is None:
                count += 1
                while self._get_next_valence():
                    pass
                if self.debug and self.current_stat != self.max_stat:
                    print("-"*20, self.points[self.stat_position - 1])
        if success:
            self.assign.bonds = bonds
        return success


for key, value in BondOrderAssignment.atomic_valence.items():
    BondOrderAssignment.atomic_valence[key] = OrderedDict(sorted(value.items(), key=lambda t: t[1]))
