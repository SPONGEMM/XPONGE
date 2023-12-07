"""
This **package** is used to calculate the TPACM4 (Transferable Partial Atomic Charge Model upto 4 bonds) charge.
"""
from copy import deepcopy
from pathlib import Path
from itertools import permutations
from .. import Xprint, Xdict, AssignRule
from ...helper import set_global_alternative_names


__all__ = ["tpacm4"]

#pylint: disable=unused-argument
def _init():
    """initialize the package"""
    file = Path(__file__).resolve().parent
    with open(file/"ATOMTYPE.dat") as fa, open(file/"CHARGE.dat") as fc:
        type_charge_mapper_ = Xdict({a: float(c) for a, c in zip(fa.read().split(), fc.read().split())})
    string_mapper_ = Xdict({"Cl": "B", "Br": "B"})
    string_mapper_.update({element: element for element in "CNOFPSH"})
    tpacm4_helper = AssignRule("_tpacm4", pure_string=True)
    @tpacm4_helper.add_rule("NH4")
    def _(i, assign):
        return assign.atoms[i] == "N" and sum(assign.bonds[i].values()) == 4

    @tpacm4_helper.add_rule("EST")
    def _(i, assign):
        if assign.atoms[i] != "O":
            return False
        for j, jorder in assign.bonds[i].items():
            if not assign.atom_judge(j, "C3") or jorder != 1:
                continue
            for k, korder in assign.bonds[j].items():
                if assign.atoms[k] == "O" and korder == 2:
                    return True
        return False

    @tpacm4_helper.add_rule("ACM")
    def _(i, assign):
        if assign.atoms[i] != "N":
            return False
        for j, jorder in assign.bonds[i].items():
            if not assign.atom_judge(j, "C3") or jorder != 1:
                continue
            for k, korder in assign.bonds[j].items():
                if assign.atoms[k] == "O" and korder == 2:
                    return True
        return False

    @tpacm4_helper.add_rule("")
    def _(i, assign):
        return True

    return type_charge_mapper_, string_mapper_

def _ring_process(ring, atom_type_alls, assign, extra_strings):
    """process the extra strings for a ring"""
    atom_numbers = len(ring.atoms)
    if atom_numbers < 6:
        for atom in ring.atoms:
            extra_strings[atom].append(f"{atom_numbers}MR")
    if ring.is_aromatic and len(ring.atoms) <= 6:
        for i, atom in enumerate(ring.atoms):
            ring_name = None
            if assign.atoms[atom] == "N":
                ring_name = "PY"
            elif assign.atoms[atom] == "O":
                ring_name = "FU"
            elif assign.atoms[atom] == "S":
                ring_name = "TF"
            if ring_name is not None:
                btom = ring.atoms[i - 1]
                for ctom in assign.bonds[btom]:
                    if "ar" not in assign.bond_marker[btom][ctom]:
                        extra_strings[ctom].append(f"o{ring_name}")
                btom = ring.atoms[(i + 1) % len(ring.atoms)]
                for ctom in assign.bonds[btom]:
                    if "ar" not in assign.bond_marker[btom][ctom]:
                        extra_strings[ctom].append(f"o{ring_name}")
                btom = ring.atoms[i - 2]
                for ctom in assign.bonds[btom]:
                    if "ar" not in assign.bond_marker[btom][ctom]:
                        extra_strings[ctom].append(f"m{ring_name}")
                btom = ring.atoms[(i + 2) % len(ring.atoms)]
                for ctom in assign.bonds[btom]:
                    if "ar" not in assign.bond_marker[btom][ctom]:
                        extra_strings[ctom].append(f"m{ring_name}")
                if len(ring.atoms) == 6:
                    btom = ring.atoms[i - 3]
                    extra_strings[btom].append(f"P{ring_name}")

_extra_string_sort = ["NH4", "3MR", "4MR", "5MR", "OEW",
                      "CC4", "OED", "CO2", "CN2", "CN3",
                      "OCO", "oPY", "mPY", "PPY", "oFU",
                      "mFU", "oTF", "mTF", "OXX"]
_extra_string_sort = Xdict({j:i for i, j in enumerate(_extra_string_sort)})
def _find_extra_string(atom_type_alls, assign):
    """find extra string for every atom"""
    extra_strings = Xdict({i:[] for i in range(assign.atom_numbers)})
    for ring in assign.rings:
        _ring_process(ring, atom_type_alls, assign, extra_strings)
    for i in range(assign.atom_numbers):
        if "3MR" in extra_strings[i] or "4MR" in extra_strings[i] or "5MR" in extra_strings[i]:
            continue
        if "NH4" in atom_type_alls[i] or "NO2" in atom_type_alls[i] or \
           "CN3" in atom_type_alls[i] or "SO2" in atom_type_alls[i]:
            for j in assign.bonds[i]:
                if assign.atoms[j] != "C":
                    continue
                for k in assign.bonds[j]:
                    if k != i and assign.atoms[k] == "C":
                        extra_strings[k].append("OEW")
        elif (("OH1" in atom_type_alls[i] or "OC1" in atom_type_alls[i]) and assign.atom_judge(i, "O2")) or \
             (("NH1" in atom_type_alls[i] or "NC1" in atom_type_alls[i]) and assign.atom_judge(i, "N3")):
            for j in assign.bonds[i]:
                if assign.atoms[j] != "C" or "CO2" in atom_type_alls[j]:
                    continue
                for k in assign.bonds[j]:
                    if k != i and assign.atoms[k] == "C":
                        extra_strings[k].append("OED")
        elif assign.atoms[i] == "F":
            for j in assign.bonds[i]:
                if assign.atoms[j] != "C":
                    continue
                for k in assign.bonds[j]:
                    if k != i and assign.atoms[k] == "C":
                        extra_strings[k].append("OXX")
        elif ("SH1" in atom_type_alls[i] or "SC1" in atom_type_alls[i]):
            for j in assign.bonds[i]:
                if assign.atoms[j] != "C":
                    continue
                for k in assign.bonds[j]:
                    if k != i and assign.atoms[k] == "C":
                        extra_strings[k].append("OSH")
    for i in range(assign.atom_numbers):
        if assign.atoms[i] == "H":
            j = next(iter(assign.bonds[i]))
            if "OED" in extra_strings[j]:
                extra_strings[i].append("OED")
            if "OEW" in extra_strings[j]:
                extra_strings[i].append("OEW")
            for k in assign.bonds[j]:
                if "CO2" in atom_type_alls[k]:
                    extra_strings[i].append("CO2")
                elif "CN2" in atom_type_alls[k]:
                    extra_strings[i].append("CN2")
                elif "NC3" in atom_type_alls[k] or "CN3" in atom_type_alls[k]:
                    extra_strings[i].append("CN3")
                elif "NH4" in atom_type_alls[k]:
                    extra_strings[i].append("NH4")
        extra_strings[i].sort(key=lambda x: _extra_string_sort.get(x, 999))
    extra_strings = Xdict({i:"".join(extra_strings[i]) for i in range(assign.atom_numbers)})
    return extra_strings

def _get_type_distance(s1, s2):
    """ get the distance between two type strings """
    if len(s1) != len(s2):
        return 0
    ans = 0
    for i, s1i in enumerate(s1):
        if i % 3 != 2 and s1i != s2[i]:
            return 0
        if i % 3 == 2 and s1i != s2[i]:
            ans += 1
    return ans


def tpacm4(assign, charge):
    """
    This **function** is used to calculate the Gasteiger charge of an assignment

    :param assign: the Assign instance
    :param charge: the total charge of the Assignment
    :return: a list of charges
    """
    assign.kekulize()
    type_suffix = assign.determine_atom_type("_tpacm4")
    atom_type_alls = []
    charges = []
    for atom in range(assign.atom_numbers):
        atom_type = []
        for btom, order in assign.bonds[atom].items():
            if "ar" in assign.bond_marker[atom][btom]:
                order = 4
            atom_type.append(f"{string_mapper[assign.atoms[atom]]}{string_mapper[assign.atoms[btom]]}{order}")
        found = False
        for one_pos in permutations(atom_type):
            temp_atom_type = "".join(one_pos) + type_suffix[atom]
            if temp_atom_type in type_charge_mapper:
                atom_type_alls.append(temp_atom_type)
                charges.append(type_charge_mapper[temp_atom_type])
                found = True
                break
        total_length = len(atom_type)
        for dis in range(0, total_length):
            if found:
                break
            for one_pos in permutations(atom_type):
                if found:
                    break
                temp_atom_type = "".join(one_pos)
                for type_, charge_ in type_charge_mapper.items():
                    if _get_type_distance(type_, temp_atom_type) == dis:
                        charges.append(charge_)
                        atom_type_alls.append(type_)
                        found = True
                        break
        if not found:
            charges.append(0)
            atom_type_alls.append("XXX")
            Xprint(f"The pattern of atom #{atom} can not be found. The molecule {assign.name} is not in the training dataset of TPACM4 charge model", "WARNING")
    Xprint(f"The simple patterns are \n{atom_type_alls}", "DEBUG")
    extra_strings = _find_extra_string(atom_type_alls, assign)
    Xprint(f"The extra strings are \n{extra_strings}", "DEBUG")
    for atom in range(assign.atom_numbers):
        neighbors = []
        for btom in assign.bonds[atom]:
            neighbors.append(atom_type_alls[btom])
        neighbors.sort(key=lambda x: x.replace("H", "Z"))
        for one_pos in permutations(neighbors):
            temp_atom_type = atom_type_alls[atom] + "".join(one_pos) + extra_strings[atom]
            if temp_atom_type in type_charge_mapper:
                charges[atom] = type_charge_mapper[temp_atom_type]
                Xprint(f"The full pattern of atom #{atom} [{atom_type_alls[atom]}] is {temp_atom_type} [c={charges[atom]}]", "DEBUG")
                break
        else:
            Xprint(f"The full pattern of atom #{atom} [{atom_type_alls[atom]}] can not be found [c={charges[atom]}]", "DEBUG")
    delta = (sum(charges) - charge) / assign.atom_numbers
    for i in range(assign.atom_numbers):
        charges[i] -= delta
    return charges

set_global_alternative_names()
type_charge_mapper, string_mapper = _init()
Xprint("""Reference for TPACM4 charge
  Mukherjee, G., Patra, N., Barua, P., Jayaram, B. 
    A Fast Empirical GAFF Compatible Partial Atomic Charge Assignment Scheme for Modeling Interactions of Small Molecules with Biomolecular Targets
    Journal of Computational Chemistry, 2010 32(5) 893-907
    DOI: 10.1002/jcc.21671""")
