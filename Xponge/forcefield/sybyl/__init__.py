"""
This **module** defines the atom types for SYBYL
"""
from ... import assign

atom_types = assign.AssignRule("sybyl", pure_string=True)

@atom_types.set_pre_action
def _(assign):
    assign.Kekulize()

@atom_types.set_post_action
def _(assign):
    replace_dict = {"N.Ccat": "N.pl3", "N.no2":  "N.pl3"}
    assign.atom_types = [replace_dict.get(t, t) for t in assign.atom_types.values()]

@atom_types.Add_Rule("C.3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "C4") or (assign.Atom_Judge(i, "C3") and assign.formal_charge[i] == -1)

@atom_types.Add_Rule("C.ar")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.atoms[i] == "C" and "AR0" in assign.atom_marker[i]

@atom_types.Add_Rule("C.cat")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.Atom_Judge(i, "C3"):
        for j in assign.bonds[i]:
            if assign.atoms[j] != "N" or "RG" in assign.atom_marker[j]:
                return False
        return True
    return False

@atom_types.Add_Rule("C.2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "C3") or (assign.Atom_Judge(i, "C2") and assign.formal_charge[i] == -1)

@atom_types.Add_Rule("C.1")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, ("C2", "C1"))

@atom_types.Add_Rule("N.ar")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.atoms[i] == "N" and "AR0" in assign.atom_marker[i]

@atom_types.Add_Rule("N.am")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.atoms[i] != "N":
        return False
    for j in assign.bonds[i]:
        if assign.atoms[j] == "C":
            for k, order in assign.bonds[j].items():
                if order == 2 and assign.atoms[k] in ("O", "S"):
                    assign.add_bond_marker(i, j, "am")
                    return True
    return False

@atom_types.Add_Rule("N.4")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "N4")

@atom_types.Add_Rule("N.1")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return (assign.Atom_Judge(i, "N1") and assign.formal_charge[i] != -1) or (assign.Atom_Judge(i, "N2") and sum(assign.bonds[i].values())  >= 4)

@atom_types.Add_Rule("N.Ccat")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.atoms[i] != "N" or "RG" in assign.atom_marker[i]:
        return False
    for j in assign.bonds[i]:
        if "AR0" not in assign.atom_marker[j] and atom_types.rules["C.cat"](j, assign):
            return True
    return False

@atom_types.Add_Rule("N.no2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    bonded_os = sum(assign.Atom_Judge(key, "O1") or assign.Atom_Judge(key, "S1") for key in assign.bonds[i].keys())
    return assign.Atom_Judge(i, "N3") and bonded_os == 2

@atom_types.Add_Rule("N.2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.atoms[i] != "N":
        return False
    for j, bond in assign.bonds[i].items():
        if bond == 2:
            return True
    return False

@atom_types.Add_Rule("N.pl3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.atoms[i] != "N":
        return False
    for j in assign.bonds[i]:
        for k, order in assign.bonds[j].items():
            if order > 1:
                return True
    return False

@atom_types.Add_Rule("N.3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.atoms[i] == "N"

@atom_types.Add_Rule("O.co2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if not assign.atom_judge(i, "O1"):
        return False
    if assign.formal_charge[i] == -1:
        return True
    j = next(iter(assign.bonds[i]))
    if not assign.atom_judge(j, "C3"):
        return False
    for k in assign.bonds[j]:
        if k != i and assign.atom_judge(k, "O1"):
            return True
    return False

@atom_types.Add_Rule("O.2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    if assign.atoms[i] != "O":
        return False
    if len(assign.bonds[i]) == 1 or "AR0" in assign.atom_marker[i]:
        return True
    for j in assign.bonds[i]:
        if "AR0" in assign.atom_marker[j]:
            return False
        for k, order in assign.bonds[j].items():
            if order == 2 and assign.atoms[j] == "C" and assign.atoms[k] in ("O", "S"):
                return False
    for j in assign.bonds[i]:
        for order in assign.bonds[j].values():
            if order > 1:
                return True
    return False

@atom_types.Add_Rule("O.3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "O2")

@atom_types.Add_Rule("S.O")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "S3")

@atom_types.Add_Rule("S.O2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "S4")

@atom_types.Add_Rule("S.2")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "S1") or (assign.Atom_Judge(i, "S2") 
                                          and "AR0" in assign.atom_marker[i].keys())

@atom_types.Add_Rule("S.3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.Atom_Judge(i, "S2")

@atom_types.Add_Rule("P.3")
def _(i, assign):
    """

    :param i:
    :param assign:
    :return:
    """
    return assign.atoms[i] == "P"

def _new_rule(element):
    @atom_types.Add_Rule(element)
    def _(i, assign):
        """

        :param i:
        :param assign:
        :return:
        """
        return assign.atoms[i] == element

for i in ["C", "O", "N", "H", "F", "Cl", "Br", "I"]:
    _new_rule(i)
