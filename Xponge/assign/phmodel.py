"""
This **module** helps to assign hydrogens according to pH
"""
from . import AssignRule, Xdict, np
#pylint: disable=unused-argument, missing-function-docstring
phmodel = AssignRule("phmodel", pure_string=True)

@phmodel.add_rule("B-phenol")
def _(i, assign):
    return assign.atoms[i] == "O" and \
           assign.formal_charge[i] == -1 and \
           "AR0" in assign.atom_marker[next(iter(assign.bonds[i]))]


@phmodel.add_rule("A-phenol")
def _(i, assign):
    if assign.atoms[i] != "H":
        return False
    j = next(iter(assign.bonds[i]))
    if assign.atoms[j] != "O":
        return False
    for k in assign.bonds[j]:
        if j != k:
            break
    else:
        return False
    return "AR0" in assign.atom_marker[k]


@phmodel.add_rule("B-carboxylic")
def _(i, assign):
    if not (assign.atoms[i] == "O" and assign.formal_charge[i] == -1):
        return False
    j = next(iter(assign.bonds[i]))
    for k, bond_order in assign.bonds[j].items():
        if bond_order == 2 and assign.atoms[k] == "O":
            return True
    return False


@phmodel.add_rule("A-carboxylic")
def _(i, assign):
    if assign.atoms[i] != "H":
        return False
    j = next(iter(assign.bonds[i]))
    if assign.atoms[j] != "O":
        return False
    for k in assign.bonds[j]:
        if j != k:
            break
    else:
        return False
    if not assign.atom_judge(k, "C3"):
        return False
    for l, bond_order in assign.bonds[k].items():
        if bond_order == 2 and assign.atoms[l] == "O":
            return True
    return False


@phmodel.add_rule("B-alcohol")
def _(i, assign):
    return assign.atoms[i] == "O" and assign.formal_charge[i] == -1


@phmodel.add_rule("A-alcohol")
def _(i, assign):
    return assign.atoms[i] == "H" and assign.atoms[next(iter(assign.bonds[i]))] == "O"


@phmodel.add_rule("N")
def _(i, assign):
    return True


class PHModelAssignment:
    """
    This **class** includes the functions to assign the pH model

    :param assign: the father Assignment instance
    :param ph: the pH value
    """
    pka_dict = Xdict({
        "carboxylic": 4.0,
        "phenol": 10.0,
        "alcohol": 15.9,
    })
    def __init__(self, assign, ph):
        self.assign = assign
        self.ph = ph

    def _get_hydrogen_position(self, i):
        """
        Get the position for the hydrogen to add
        """
        x, y, z = 0, 0, 0
        for j in self.assign.bonds[i]:
            displacement = self.assign.coordinate[i] - self.assign.coordinate[j]
            displacement = self.assign.coordinate[i] + displacement / np.linalg.norm(displacement)
            x += displacement[0]
            y += displacement[1]
            z += displacement[2]
        return x, y, z

    def main(self):
        """
        This **function** is the main function to do the pH model assignment

        :return: the sum of the formal charge now
        """
        phtypes = self.assign.determine_atom_type("phmodel")
        to_add = []
        to_delete = []
        for i, atype in phtypes.items():
            if atype != "N":
                a_or_b, name = atype.split("-")
                pka = self.pka_dict[name]
                if a_or_b == "A" and pka < self.ph:
                    to_delete.append(i)
                elif a_or_b == "B" and pka > self.ph:
                    to_add.append(i)
        for i in to_add:
            x, y, z = self._get_hydrogen_position(i)
            self.assign.add_atom("H", x=x, y=y, z=z)
            self.assign.add_bond(self.assign.atom_numbers - 1, i, 1)

        to_delete.sort(reverse=True)
        for i in to_delete:
            self.assign.delete_atom(i)
        self.assign.determine_bond_order()
        return int(round(sum(self.assign.formal_charge)))

set_global_alternative_names()
