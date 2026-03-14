"""
This **module** is used to load and read
"""
import os
import re
import io

import numpy as np

from .helper import Molecule, Residue, ResidueType, AtomType, GlobalSetting, Xdict, Xopen, \
    set_real_global_variable, set_global_alternative_names, Xprint
from .helper.math import get_length_angle_from_basis_vectors


##########################################################################
# General Format
##########################################################################
def _mol2_atom(line, current_residue_index, current_residue, ignore_atom_type, temp, current_molecule,
               atom_residue_map):
    """
        parse atom lines in a mol2 file
    """
    words = line.split()
    if current_residue_index is None or int(words[6]) != current_residue_index:
        current_residue_index = int(words[6])
        if words[7] not in ResidueType.get_all_types():
            set_real_global_variable(words[7], ResidueType(name=words[7]))
            temp = True
        else:
            temp = False
        if current_residue:
            current_molecule.Add_Residue(current_residue)
        current_residue = Residue(ResidueType.get_type(words[7]))
    if ignore_atom_type:
        temp_atom_type = AtomType.get_type("UNKNOWN")
    else:
        temp_atom_type = AtomType.get_type(words[5])

    if temp:
        current_residue.type.Add_Atom(words[1], temp_atom_type, *words[2:5])
        current_residue.type.atoms[-1].Update(**{"charge[e]": float(words[8])})
    current_residue.Add_Atom(words[1], temp_atom_type, *words[2:5])
    current_residue.atoms[-1].Update(**{"charge[e]": float(words[8])})
    atom_residue_map[words[0]] = [words[1], current_residue, current_residue_index, temp, current_residue.atoms[-1],
                                  current_residue.type.atoms[-1]]
    return current_residue_index, current_residue, temp


def _mol2_bond(line, current_molecule, atom_residue_map):
    """
        parse bond lines in a mol2 file
    """
    words = line.split()
    if atom_residue_map[words[1]][1] == atom_residue_map[words[2]][1]:
        if atom_residue_map[words[1]][3]:
            atom_residue_map[words[1]][1].type.Add_Connectivity(atom_residue_map[words[1]][0],
                                                                atom_residue_map[words[2]][0])
        atom_residue_map[words[1]][1].Add_Connectivity(atom_residue_map[words[1]][0], atom_residue_map[words[2]][0])
    else:
        current_molecule.Add_Residue_Link(atom_residue_map[words[1]][4], atom_residue_map[words[2]][4])
        index_diff = atom_residue_map[words[1]][2] - atom_residue_map[words[2]][2]
        if abs(index_diff) == 1:
            if atom_residue_map[words[1]][3]:
                if index_diff < 0:
                    atom_residue_map[words[1]][1].type.tail = atom_residue_map[words[1]][0]
                else:
                    atom_residue_map[words[1]][1].type.head = atom_residue_map[words[1]][0]
            if atom_residue_map[words[2]][3]:
                if index_diff < 0:
                    atom_residue_map[words[2]][1].type.head = atom_residue_map[words[2]][0]
                else:
                    atom_residue_map[words[2]][1].type.tail = atom_residue_map[words[2]][0]


def _mol2_template_atom(line, current_residue_index, temp, current_residue, atom_residue_map):
    """
        parse atom lines in a mol2 file as a template
    """
    words = line.split()
    if current_residue_index is None or int(words[6]) != current_residue_index:
        current_residue_index = int(words[6])
        if words[7] not in ResidueType.get_all_types():
            set_real_global_variable(words[7], ResidueType(name=words[7]))
            temp = True
        else:
            temp = False
        current_residue = ResidueType.get_type(words[7])

    temp_atom_type = AtomType.get_type(words[5])

    if temp:
        current_residue.Add_Atom(words[1], temp_atom_type, *words[2:5])
        current_residue.atoms[-1].Update(**{"charge[e]": float(words[8])})
    atom_residue_map[words[0]] = [words[1], current_residue, current_residue_index, temp]
    return current_residue_index, current_residue, temp


def _mol2_template_bond(line, atom_residue_map):
    """
        parse bond lines in a mol2 file as a template
    """
    words = line.split()
    atom1_index, atom2_index = words[1], words[2]
    atom_name, residue_type, residue_index, is_new = 0, 1, 2, 3
    atom1_info, atom2_info = atom_residue_map[atom1_index], atom_residue_map[atom2_index]
    if atom1_info[residue_type] == atom2_info[residue_type] and atom1_info[is_new]:
        atom1_info[residue_type].Add_Connectivity(atom1_info[atom_name], atom2_info[atom_name])
    else:
        index_diff = atom1_info[residue_index] - atom2_info[residue_index]
        if abs(index_diff) == 1:
            if atom1_info[is_new]:
                if index_diff < 0:
                    atom1_info[residue_type].tail = atom1_info[atom_name]
                else:
                    atom1_info[residue_type].head = atom1_info[atom_name]
            if atom2_info[is_new]:
                if index_diff < 0:
                    atom2_info[residue_type].head = atom2_info[atom_name]
                else:
                    atom2_info[residue_type].tail = atom2_info[atom_name]


def _mol2_template(file):
    """
    This is used to make loading more efficient when importing the force fields
    """
    current_residue_index = None
    current_residue = None
    atom_residue_map = Xdict()
    temp = None
    with file as f:
        flag = None
        for line in f:
            if not line.strip():
                continue
            if line.startswith("@<TRIPOS>"):
                flag = line[9:].strip()
            elif flag == "ATOM":
                current_residue_index, current_residue, temp = _mol2_template_atom(line, current_residue_index, temp,
                                                                                   current_residue, atom_residue_map)
            elif flag == "BOND":
                _mol2_template_bond(line, atom_residue_map)


def load_mol2(file, ignore_atom_type=False, as_template=False):
    """
    This **function** is used to load a mol2 file

    :param file: the name of the input file or an instance of io.IOBase
    :param ignore_atom_type: ignore the atom types in the mol2 file
    :param as_template: only read the mol2 file as some residue types and no molecule will created
        **New From 1.2.6.8**
    :return: a Molecule instance if as_template is False
    """
    if not isinstance(file, io.IOBase):
        file = open(file)
    if as_template:
        return _mol2_template(file)
    with file as f:
        # 存储读的时候的临时信息，key是编号
        # value是list：原子名(0)、residue(1)、residue编号(2)、是否是新的residue type(3)、该原子(4)、residue type的最新原子(5)
        atom_residue_map = {}
        flag = None
        nline = 0
        current_molecule = None
        current_residue = None
        current_residue_index = None
        temp = None
        for line in f:
            if line.strip():
                nline += 1
            else:
                continue

            if line.startswith("@<TRIPOS>"):
                if flag == "ATOM":
                    current_molecule.Add_Residue(current_residue)
                flag = line[9:].strip()
                nline = 0
            elif flag == "MOLECULE":
                if nline == 1:
                    current_molecule = Molecule(line.strip())
            elif flag == "ATOM":
                current_residue_index, current_residue, temp = _mol2_atom(line, current_residue_index, current_residue,
                                                                          ignore_atom_type, temp, current_molecule,
                                                                          atom_residue_map)
            elif flag == "BOND":
                _mol2_bond(line, current_molecule, atom_residue_map)
    return current_molecule


def _pdb_ssbond_before(chain, residue_type_map, ssbonds):
    """
        change the residue name of the SS bonds to CYX
    """
    for ssbond in ssbonds:
        res_a_index=chain[ssbond[15]][int(ssbond[17:21])]
        if residue_type_map[res_a_index] in GlobalSetting.PDBResidueNameMap["head"].values():
            residue_type_map[res_a_index] = "NCYX"
        elif residue_type_map[res_a_index] in GlobalSetting.PDBResidueNameMap["tail"].values():
            residue_type_map[res_a_index] = "CCYX"
        else:
            residue_type_map[res_a_index] = "CYX"
            
        res_b_index= chain[ssbond[29]][int(ssbond[31:35])]
        if residue_type_map[res_b_index] in GlobalSetting.PDBResidueNameMap["head"].values():
            residue_type_map[res_b_index] = "NCYX"
        elif residue_type_map[res_b_index] in GlobalSetting.PDBResidueNameMap["tail"].values():
            residue_type_map[res_b_index] = "CCYX"
        else:
            residue_type_map[res_b_index] = "CYX"


def _pdb_ssbond_after(chain, ssbonds, molecule):
    """
        connect the SS bonds
    """
    for ssbond in ssbonds:
        res_a_pdb = _pdb_parse_atom_serial(ssbond[17:21])
        res_b_pdb = _pdb_parse_atom_serial(ssbond[31:35])
        if res_a_pdb is None or res_b_pdb is None:
            continue
        res_a_index = chain[ssbond[15]][res_a_pdb]
        res_b_index = chain[ssbond[29]][res_b_pdb]
        if res_a_index > res_b_index:
            res_a_index, res_b_index = res_b_index, res_a_index
        res_a = molecule.residues[res_a_index]
        res_b = molecule.residues[res_b_index]
        molecule.Add_Residue_Link(res_a.name2atom(res_a.type.connect_atoms["ssbond"]),
                                  res_b.name2atom(res_b.type.connect_atoms["ssbond"]))


def _pdb_link_after(chain, links, molecule):
    """
        connect residue links
    """
    for link in links:
        chain_a, chain_b = link[21], link[51]
        res_a_pdb = _pdb_parse_atom_serial(link[22:26])
        res_b_pdb = _pdb_parse_atom_serial(link[52:56])
        if res_a_pdb is None or res_b_pdb is None:
            continue
        if chain_a not in chain or chain_b not in chain or \
                res_a_pdb not in chain[chain_a] or res_b_pdb not in chain[chain_b]:
            Xprint(f"The link between {chain_a}{res_a_pdb} and {chain_b}{res_b_pdb} is not valid", "WARNING")
            continue
        res_a_index = chain[link[21]][res_a_pdb]
        res_b_index = chain[link[51]][res_b_pdb]
        if res_a_index > res_b_index:
            res_a_index, res_b_index = res_b_index, res_a_index
        res_a = molecule.residues[res_a_index]
        res_b = molecule.residues[res_b_index]
        molecule.Add_Residue_Link(res_a.name2atom(link[12:16].strip()),
                                  res_b.name2atom(link[42:46].strip()))


def _pdb_add_atom(current_residue, atomname, x, y, z,
                  ignore_hydrogen, ignore_unknown_name, atom_map, atom_index):
    """
        add the atom to the residue
    """
    if ignore_hydrogen and (atomname.startswith("H") or
                            (len(atomname) > 1 and atomname[0] in "123" and atomname[1] == "H")):
        return False
    try:
        current_residue.Add_Atom(atomname, x=x, y=y, z=z)
        if atom_index not in atom_map:
            atom_map[atom_index] = current_residue.atoms[-1]
        return True
    except KeyError as ke:
        if ignore_unknown_name:
            return False
        raise ke


def _pdb_add_residue(f, molecule, position_need, residue_type_map, ignore_hydrogen, ignore_unknown_name):
    """
        add a residue to the molecule
    """
    atom_map = Xdict()
    current_residue_count = -1
    current_residue_index = None
    current_insertion_code = None
    current_residue = None
    current_resname = None
    links = []
    f.seek(0)
    for line in f:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            extra = line[16]
            insertion_code = line[26]
            resindex = _pdb_parse_atom_serial(line[22:26])
            if resindex is None:
                continue
            resname = line[17:20].strip()
            atomname = line[12:16].strip()
            atom_index = _pdb_parse_atom_serial(line[6:11])
            if atom_index is None:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if current_residue_index is None or current_residue_index != resindex or \
                    current_resname != resname or current_insertion_code != insertion_code:
                current_residue_count += 1
                if current_residue:
                    molecule.Add_Residue(current_residue)
                    if current_residue_index is not None and current_residue.type.tail and ResidueType.get_type(
                            residue_type_map[current_residue_count]).head:
                        links.append(len(molecule.residues))
                current_residue = Residue(ResidueType.get_type(residue_type_map[current_residue_count]))
                current_residue_index = resindex
                current_insertion_code = insertion_code
                current_resname = resname
            if extra not in (" ", position_need):
                continue
            _pdb_add_atom(current_residue, atomname, x, y, z,
                          ignore_hydrogen, ignore_unknown_name, atom_map, atom_index)
        elif line.startswith("TER"):
            current_residue_index = None
            current_resname = None

    if current_residue:
        molecule.Add_Residue(current_residue)
    for count in links:
        res_a = molecule.residues[count]
        res_b = molecule.residues[count - 1]
        molecule.Add_Residue_Link(res_a.name2atom(res_a.type.head), res_b.name2atom(res_b.type.tail))
    return atom_map


def _pdb_judge_histone(judge_histone, residue_type_map, current_histone_information):
    """
        judge the protonation state of the histone residue
    """
    if judge_histone and residue_type_map and residue_type_map[-1] in GlobalSetting.HISMap["HIS"].keys():
        if current_histone_information["DeltaH"]:
            if current_histone_information["EpsilonH"]:
                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIP"]
            else:
                residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HID"]
        else:
            residue_type_map[-1] = GlobalSetting.HISMap["HIS"][residue_type_map[-1]]["HIE"]
        current_histone_information["DeltaH"] = False
        current_histone_information["EpsilonH"] = False


def _pdb_read_sequences(line, sequences):
    """
        read the sequence
    """
    chain_id = line[11]
    if chain_id not in sequences:
        sequences[chain_id] = []
    sequence = sequences[chain_id]
    sequence.extend(line[19:].split())


def _pdb_find_missing_residues(mol, sequences, chain, residue_type_map):
    """
        find missing residues
    """
    if sequences:
        original_residue_names = [GlobalSetting.PDBResidueNameMap["save"].get(name, name) for name in residue_type_map]
        for chain_id, sequence in sequences.items():
            pdb_index_to_mol_index = chain[chain_id]
            name_index_to_mol_index = Xdict()
            pdb_index = list(pdb_index_to_mol_index.keys())
            pdb_index.sort()
            names_length = pdb_index[-1] - pdb_index[0] + 1
            if pdb_index[-1] > 0 > pdb_index[0] and 0 not in pdb_index_to_mol_index:
                names_length -= 1
            names_with_none = [None for _ in range(names_length)]
            for i in pdb_index:
                index = i - pdb_index[0]
                if pdb_index[0] < 0 < i and 0 not in pdb_index_to_mol_index:
                    index -= 1
                names_with_none[index] = original_residue_names[pdb_index_to_mol_index[i]]
                name_index_to_mol_index[index] = pdb_index_to_mol_index[i]
            offset = 0
            for offset in range(len(sequence) - len(names_with_none) + 1):
                if all(b in (a, None) or
                       (a == "HIS" and b in (None, "HID", "HIE", "HIP"))
                       for a, b in zip(sequence[offset:], names_with_none)):
                    break
            if offset:
                mol.set_missing_residues_info(None, name_index_to_mol_index[0], sequence[:offset])
            tail_offset = len(sequence) - offset - len(names_with_none)
            if tail_offset:
                tail_index = name_index_to_mol_index[len(names_with_none) - 1]
                mol.set_missing_residues_info(tail_index, None, sequence[-tail_offset:])
            find_none_index = None
            for i, name in enumerate(names_with_none):
                if find_none_index is None and name is None:
                    find_none_index = i
                elif find_none_index is not None and name is not None:
                    start = name_index_to_mol_index[find_none_index - 1]
                    end = name_index_to_mol_index[i]
                    seq_start = find_none_index + offset
                    seq_end = i + offset
                    mol.set_missing_residues_info(start, end, sequence[seq_start:seq_end])
                    find_none_index = None


_HY36_DIGITS_UPPER = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_HY36_DIGITS_LOWER = _HY36_DIGITS_UPPER.lower()
_HY36_UPPER_VALUES = {digit: value for value, digit in enumerate(_HY36_DIGITS_UPPER)}
_HY36_LOWER_VALUES = {digit: value for value, digit in enumerate(_HY36_DIGITS_LOWER)}


def _hy36_decode_pure(digits_values, text):
    """
        decode string using digits-values mapping
    """
    result = 0
    base = len(digits_values)
    for char in text:
        result = result * base + digits_values[char]
    return result


def _pdb_hybrid36_decode(width, field):
    """
        decode hybrid-36 encoded integer field (e.g. PDB atom serial)
    """
    if len(field) != width:
        raise ValueError("invalid number literal.")
    first = field[0]
    if first in ("-", " ") or first.isdigit():
        try:
            return int(field)
        except ValueError:
            if field == " " * width:
                return 0
    elif first in _HY36_UPPER_VALUES:
        try:
            return _hy36_decode_pure(_HY36_UPPER_VALUES, field) - 10 * 36 ** (width - 1) + 10 ** width
        except KeyError:
            pass
    elif first in _HY36_LOWER_VALUES:
        try:
            return _hy36_decode_pure(_HY36_LOWER_VALUES, field) + 16 * 36 ** (width - 1) + 10 ** width
        except KeyError:
            pass
    raise ValueError("invalid number literal.")


def _pdb_parse_atom_serial(field):
    """
        parse PDB atom serial with possible hybrid-36 encoding
    """
    text = field.strip()
    if not text:
        return None
    try:
        if text[0] in " +-0123456789":
            return int(text)
        return _pdb_hybrid36_decode(len(field), field)
    except ValueError as exc:
        Xprint(f"Invalid atom serial {field!r}: {exc}", "WARNING")
        return None


def _pdb_connects(mol, connects, atom_map):
    """
        process the CONECT part of the pdb file
    """
    for connect in connects:
        atom = _pdb_parse_atom_serial(connect[6:11])
        if atom is None:
            continue
        if atom not in atom_map:
            Xprint(f"CONECT of {atom} is not valid \
because {atom} is not found in the file", "WARNING")
            continue
        atom_ = atom_map[atom]
        connect = connect[11:]
        for i in range(0, len(connect), 5):
            to_connect = connect[i:i+5]
            if not to_connect.strip():
                break
            to_connect = _pdb_parse_atom_serial(to_connect)
            if to_connect is None:
                continue
            if to_connect not in atom_map:
                Xprint(f"CONECT of {atom} to {to_connect} is not valid \
because {to_connect} is not found in the file", "WARNING")
                continue
            to_connect_ = atom_map[to_connect]
            if mol.get_residue_link(atom_, to_connect_):
                continue
            if atom_.residue == to_connect_.residue:
                Xprint(f"CONECT of {atom} to {to_connect} is not valid \
because {to_connect} and {atom} is in one residue", "WARNING")
                continue
            mol.add_residue_link(atom_, to_connect_)


def load_pdb(file, judge_histone=True, position_need="A", ignore_hydrogen=False,
             ignore_unknown_name=False, ignore_seqres=True, ignore_conect=True, read_cryst1=True):
    """
    This **function** is used to load a pdb file

    :param file: the name of the input file or an instance of io.IOBase
    :param judge_histone: judge the protonized state of the histone residues
    :param position_need: the position character to read
    :param ignore_hydrogen: do not read the atom with a name beginning with "H" or "[123]H"
    :param ignore_unknown_name: do not read the atom with an unknown name **New From 1.2.6.4**
    :param ignore_seqres: do not read the SEQRES lines **New From 1.2.6.7**
    :param ignore_conect: do not read the CONECT lines **New From 1.2.6.7**
    :param read_cryst1: read the CRYST1 line to set box length and angle
    :return: a Molecule instance
    """
    if not isinstance(file, io.IOBase):
        filename = file
        file = open(file)
    else:
        filename = "in-memory-string"
    molecule = Molecule(os.path.splitext(os.path.basename(filename))[0])
    chain = Xdict(not_found_method=lambda key: Xdict())
    sequences = Xdict()
    ssbonds = []
    links = []
    connects = []
    residue_type_map = []
    insertion_count = 0
    current_residue_count = -1
    current_insertion_code = None
    current_residue_index = None
    current_resname = None
    chain_id_processed = set()
    current_histone_information = {"DeltaH": False, "EpsilonH": False}
    cryst1 = None
    with file as f:
        for line in f:
            if read_cryst1 and line.startswith("CRYST1"):
                try:
                    a = float(line[6:15])
                    b = float(line[15:24])
                    c = float(line[24:33])
                    alpha = float(line[33:40])
                    beta = float(line[40:47])
                    gamma = float(line[47:54])
                    cryst1 = ([a, b, c], [alpha, beta, gamma])
                except ValueError:
                    Xprint(f"Invalid CRYST1 line: {line.strip()}", "WARNING")
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resindex = _pdb_parse_atom_serial(line[22:26])
                if resindex is None:
                    continue
                resindex += insertion_count
                insertion_code = line[26]
                chain_id = line[21]
                if not chain_id.strip() or chain_id in chain_id_processed:
                    chain_id = 0
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                if current_residue_index is None:
                    _pdb_judge_histone(judge_histone, residue_type_map, current_histone_information)
                    current_residue_count += 1
                    current_resname = resname
                    current_insertion_code = insertion_code
                    if resname in GlobalSetting.PDBResidueNameMap["head"].keys():
                        resname = GlobalSetting.PDBResidueNameMap["head"][resname]
                    residue_type_map.append(resname)
                    current_residue_index = resindex
                    chain[chain_id][resindex] = current_residue_count
                elif (current_residue_index != resindex or current_insertion_code != insertion_code) or \
                        current_resname != resname:
                    _pdb_judge_histone(judge_histone, residue_type_map, current_histone_information)
                    current_residue_count += 1
                    current_resname = resname
                    current_insertion_code = insertion_code
                    if insertion_code != " ":
                        insertion_count += 1
                        resindex += 1
                    current_residue_index = resindex
                    chain[chain_id][resindex] = current_residue_count
                    residue_type_map.append(resname)
                if judge_histone and resname in GlobalSetting.HISMap["HIS"].keys():
                    if atomname == GlobalSetting.HISMap["DeltaH"]:
                        current_histone_information["DeltaH"] = True
                    elif atomname == GlobalSetting.HISMap["EpsilonH"]:
                        current_histone_information["EpsilonH"] = True

            elif line.startswith("TER"):
                current_residue_index = None
                current_resname = None
                current_insertion_code = None
                insertion_count = 0
                chain_id_processed.add(chain_id)
                if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
                    residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]
                _pdb_judge_histone(judge_histone, residue_type_map, current_histone_information)

            elif line.startswith("SSBOND"):
                ssbonds.append(line)

            elif line.startswith("LINK"):
                links.append(line)

            elif line.startswith("SEQRES") and not ignore_seqres:
                _pdb_read_sequences(line, sequences)

            elif line.startswith("CONECT") and not ignore_conect:
                connects.append(line)

        current_residue_index = None
        if residue_type_map[-1] in GlobalSetting.PDBResidueNameMap["tail"].keys():
            residue_type_map[-1] = GlobalSetting.PDBResidueNameMap["tail"][residue_type_map[-1]]

        _pdb_ssbond_before(chain, residue_type_map, ssbonds)
        atom_map = _pdb_add_residue(file, molecule, position_need,
                                    residue_type_map, ignore_hydrogen, ignore_unknown_name)
        _pdb_ssbond_after(chain, ssbonds, molecule)
        _pdb_link_after(chain, links, molecule)
        _pdb_find_missing_residues(molecule, sequences, chain, residue_type_map)
        _pdb_connects(molecule, connects, atom_map)
        if cryst1 is not None:
            molecule.box_length = cryst1[0]
            molecule.box_angle = cryst1[1]

    return molecule


##########################################################################
# SPONGE Format
##########################################################################
def load_coordinate(filename, mol=None):
    """
    This **function** is used to read the SPONGE coordinate in file

    :param filename: the coordinate file to load
    :param mol: the molecule or residue to load the coordinate into
    :return: two numpy arrays, representing the coordinates and the box information respectively
    """
    with Xopen(filename, "r") as f:
        atom_numbers = int(f.readline().split()[0])
        crd = np.zeros((atom_numbers, 3), dtype=np.float32)
        box = np.zeros(6, dtype=np.float32)
        for i in range(atom_numbers):
            crd[i][:] = np.array([float(x) for x in f.readline().split()])
        box = np.array(f.readline().split(), dtype=np.float32)
    if mol:
        for i, atom in enumerate(mol.get_atoms()):
            atom.x = crd[i][0]
            atom.y = crd[i][1]
            atom.z = crd[i][2]
        if isinstance(mol, Molecule):
            mol.box_length = box[:3]
    return crd, box


##########################################################################
# amber Format
##########################################################################
def _frcmod_nb14(line, atoms):
    nb14ee = re.findall(r"SCEE=[\d+\.]+", line)
    nb14lj = re.findall(r"SCNB=[\d+\.]+", line)
    if not nb14ee and not nb14lj:
        return ""
    nb14ee = 1.0 / float(nb14ee[0].split("=")[1]) if nb14ee else 1.0 / 1.2
    nb14lj = 1.0 / float(nb14lj[0].split("=")[1]) if nb14lj else 1.0 / 2.0
    return f"{atoms[0]}-{atoms[3]} {nb14ee} {nb14lj}\n"


def _frcmod_cmap(line, cmap, temp_cmp, cmap_flag):
    """

    :param line:
    :param cmap:
    :param temp_cmp:
    :param cmap_flag:
    :return:
    """
    if line.startswith("%FLAG"):
        if "CMAP_COUNT" in line:
            if temp_cmp:
                for res in temp_cmp["residues"]:
                    cmap[f"C-N-{res}@XC-C-N"] = {"resolution": temp_cmp["info"]["resolution"],
                                                 "parameters": temp_cmp["info"]["parameters"]}
            temp_cmp = {"residues": [], "info": {"resolution": 24, "count": int(line.split()[-1]), "parameters": []}}
            cmap_flag = "CMAP_COUNT"
        elif "CMAP_RESOLUTION" in line:
            temp_cmp["info"]["resolution"] = int(line.split()[-1])
            cmap_flag = "CMAP_RESOLUTION"
        elif "CMAP_RESLIST" in line:
            cmap_flag = "CMAP_RESLIST"
        elif "CMAP_TITLE" in line:
            cmap_flag = "CMAP_TITLE"
        elif "CMAP_PARAMETER" in line:
            cmap_flag = "CMAP_PARAMETER"
    elif cmap_flag == "CMAP_RESLIST":
        temp_cmp["residues"].extend(line.split())
    elif cmap_flag == "CMAP_PARAMETER":
        temp_cmp["info"]["parameters"].extend([float(x) for x in line.split()])
    return temp_cmp, cmap_flag


def _frcmod_atoms_words(line, n, last_atoms=None):
    """

    :param line:
    :param n:
    :return:
    """
    if line[0] != " ":
        return [word.strip() for word in line[:n].split("-")], line[n:].split()
    return [last_atoms, line[n:].split()]


def load_frcmod(filename, nbtype="RE"):
    """
    This **function** is used to load a frcmod file

    :param filename: the name of the file to load
    :param nbtype: the non-bonded interaction recording type in the frcmod file.
    :return: a list of strings, including atoms, bonds, angles, propers, impropers, ljs, cmap information respectively
    """
    with open(filename) as f:
        f.readline()
        flag = None
        atom_types = {}  # 元素符号和质量
        bonds = ["name  k[kcal/mol·A^-2]    b[A]\n"]
        angles = ["name  k[kcal/mol·rad^-2]    b[degree]\n"]
        propers = ["name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"]
        reset = 1
        impropers = ["name  k[kcal/mol]    phi0[degree]    periodicity\n"]
        cmap = {}
        cmap_flag = None
        temp_cmp = {"residues": []}
        if nbtype == "SK":
            raise NotImplementedError
        if nbtype == "AC":
            ljs = ["name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"]
        elif nbtype == "RE":
            ljs = ["name rmin[A]   epsilon[kcal/mol]\n"]

        for line in f:
            if not line.strip():
                continue
            words = line.split()
            if flag != "CMAP" and len(words) == 1:
                flag = line.strip()
            elif flag[:4] == "MASS":
                atom_types[words[0]] = words[1]
            elif flag[:4] == "BOND":
                atoms, words = _frcmod_atoms_words(line, 5)
                bonds.append("-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n")
            elif flag[:4] == "ANGL":
                atoms, words = _frcmod_atoms_words(line, 8)
                angles.append("-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n")
            elif flag[:4] == "DIHE":
                atoms, words = _frcmod_atoms_words(line, 11)
                propers.append(
                    "-".join(atoms) + "\t" + str(float(words[1]) / int(words[0])) + "\t" + words[2] + "\t" + str(
                        abs(int(float(words[3])))) + "\t" + str(reset) + "\n"
                )
                if int(float(words[3])) < 0:
                    reset = 0
                else:
                    reset = 1
            elif flag[:4] == "IMPR":
                atoms, words = _frcmod_atoms_words(line, 11)
                impropers.append(
                    "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\t" + str(int(float(words[2]))) + "\n"
                )
            elif flag[:4] == "NONB":
                words = line.split()
                ljs.append(words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n")
            elif flag[:4] == "CMAP":
                temp_cmp, cmap_flag = _frcmod_cmap(line, cmap, temp_cmp, cmap_flag)

    for res in temp_cmp["residues"]:
        cmap[f"C-N-{res}@XC-C-N"] = {"resolution": temp_cmp["info"]["resolution"],
                                     "parameters": temp_cmp["info"]["parameters"]}
    atoms = ["name  mass  LJtype\n"]
    for atom, mass in atom_types.items():
        atoms.append(atom + "\t" + mass + "\t" + atom + "\n")
    toret = [
        "".join(atoms),
        "".join(bonds),
        "".join(angles),
        "".join(propers),
        "".join(impropers),
        "".join(ljs),
        cmap,
    ]
    return toret


def _parmdat_read_harmonic_bonds(f, bonds, n):
    """

    :param f:
    :param bonds:
    :param n:
    :return:
    """
    for line in f:
        if not line.strip():
            break
        atoms, words = _frcmod_atoms_words(line, n)
        bonds.append("-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\n")
    return bonds


def load_parmdat(filename):
    """
    This **function** is used to load a parmdat file

    :param filename: the name of the file to load
    :return: a list of strings, including atoms, bonds, angles, propers, impropers, ljs, nb14s information respectively
    """
    with open(filename) as f:
        f.readline()
        # 读原子
        atom_types = {}  # 元素符号和质量
        lj_types = Xdict()  # 元素符号和LJ类型
        for line in f:
            if not line.strip():
                break

            words = line.split()
            atom_types[words[0]] = words[1]
            lj_types[words[0]] = words[0]
        f.readline()
        # 读键长
        bonds = ["name  k[kcal/mol·A^-2]    b[A]\n"]
        bonds = _parmdat_read_harmonic_bonds(f, bonds, 5)

        # 读键角
        angles = ["name  k[kcal/mol·rad^-2]    b[degree]\n"]
        angles = _parmdat_read_harmonic_bonds(f, angles, 8)

        # 读恰当二面角
        reset = 1
        propers = ["name  k[kcal/mol]    phi0[degree]    periodicity    reset\n"]
        nb14s = ["name    kLJ     kee\n"]
        atoms = None
        for line in f:
            if not line.strip():
                break
            last_atoms = atoms
            atoms, words = _frcmod_atoms_words(line, 11, last_atoms)
            nb14s.append(_frcmod_nb14(line, atoms))
            propers.append(
                "-".join(atoms) + "\t" + str(float(words[1]) / int(words[0])) + "\t" + words[2] + "\t" + str(
                    abs(int(float(words[3])))) + "\t" + str(reset) + "\n"
            )
            if int(float(words[3])) < 0:
                reset = 0
            else:
                reset = 1

        # 读非恰当二面角
        impropers = ["name  k[kcal/mol]    phi0[degree]    periodicity\n"]
        for line in f:
            if not line.strip():
                break
            atoms, words = _frcmod_atoms_words(line, 11)
            impropers.append(
                "-".join(atoms) + "\t" + words[0] + "\t" + words[1] + "\t" + str(int(float(words[2]))) + "\n"
            )

        # 跳过水的信息
        f.readline()
        f.readline()

        # 读LJ种类
        for line in f:
            if not line.strip():
                break
            atoms = line.split()
            atom0 = atoms.pop(0)
            for atom in atoms:
                lj_types[atom] = atom0

        # 读LJ信息
        word = f.readline().split()[1]
        if word == "SK":
            raise NotImplementedError
        if word == "AC":
            ljs = ["name A[kcal/mol·A^-12]   B[kcal/mol·A^-6]\n"]
        elif word == "RE":
            ljs = ["name rmin[A]   epsilon[kcal/mol]\n"]

        for line in f:
            if not line.strip():
                break
            words = line.split()
            ljs.append(words[0] + "-" + words[0] + "\t" + words[1] + "\t" + words[2] + "\n")

    atoms = ["name  mass  LJtype\n"]
    for atom, mass in atom_types.items():
        atoms.append(atom + "\t" + mass + "\t" + lj_types[atom] + "\n")
    toret = [
        "".join(atoms),
        "".join(bonds),
        "".join(angles),
        "".join(propers),
        "".join(impropers),
        "".join(ljs),
        "".join(nb14s),
    ]
    return toret


def load_rst7(filename, mol=None):
    """
    This **function** is used to load a rst7 file

    :param filename: the name of the file to load
    :param mol: the molecule to load the coordinates
    :return: a tuple including coordinates and box information
    """
    crds = []
    with open(filename) as f:
        f.readline()
        words = f.readline().split()
        atom_numbers = int(words[0])
        line = ""
        for line in f:
            words = line.split()
            while words and len(crds) < atom_numbers * 3:
                crds.append(float(words.pop(0)))
        box = [float(i) for i in line.split()]
    crds = np.array(crds).reshape((-1, 3))
    mol.box_length = box[:3]
    count = -1
    for residue in mol.residues:
        for atom in residue.atoms:
            count += 1
            atom.x = crds[count][0]
            atom.y = crds[count][1]
            atom.z = crds[count][2]

    return crds, box


##########################################################################
# GROMACS Format
##########################################################################
class GromacsTopologyIterator():
    """
    This **class** is used to read a GROMACS topology file

    usage example::

        f = GromacsTopologyIterator("example.itp")
        for line in f:
            print(line)

    :param filename: the name of the file to read
    :param macros: the macros used to read the Gromacs topology file

    """

    def __init__(self, filename=None, macros=None):
        self.files = []
        self.filenames = []
        self.stack = []
        self.flag = ""
        self.macro_define_stat = []

        if macros:
            self.defined_macros = macros
        else:
            self.defined_macros = {}
        if filename:
            self._add_iterator_file(filename)

    def __iter__(self):
        self.flag = ""
        self.macro_define_stat = []
        self.stack = []
        return self

    def __next__(self):
        while self.files:
            f = self.files[-1]
            line = f.readline()
            if line:
                line = self._line_preprocess(line)
                if line[0] == "#":
                    line = self._line_define(line)
                elif self.macro_define_stat and not self.macro_define_stat[-1]:
                    self.stack.append("false line")
                    line = next(self)
                    self.stack.pop()
                elif "[" in line and "]" in line:
                    self.flag = line.strip()[1:-1].strip()
                    self.stack.append(f"flag line: {self.flag}")
                    line = next(self)
                    self.stack.pop()
                for macro, tobecome in self.defined_macros.items():
                    line = line.replace(macro, tobecome)
                return line

            f.close()
            self.files.pop()
            self.filenames.pop()

        raise StopIteration

    def _add_iterator_file(self, filename):
        """

        :param filename:
        :return:
        """
        if self.files:
            filename = os.path.abspath(os.path.join(os.path.dirname(self.filenames[-1]), filename.replace('"', '')))
        else:
            filename = os.path.abspath(filename.replace('"', ''))

        if not os.path.exists(filename):
            base_name = os.path.basename(filename)
            for include_dir in GlobalSetting.GMXIncludePaths:
                candidate = os.path.abspath(os.path.join(include_dir, base_name))
                if os.path.exists(candidate):
                    filename = candidate
                    break

        f = Xopen(filename, "r")
        self.files.append(f)
        self.filenames.append(filename)

    def _line_preprocess(self, line):
        """

        :param line:
        :return:
        """
        line = line.strip()
        comment = line.find(";")
        if comment >= 0:
            line = line[:comment]
        while line and line[-1] == "\\":
            self.stack.append("extra line")
            line = line[:-1] + " " + next(self).strip()
            self.stack.pop()
        if not line:
            self.stack.append("blank line")
            line = next(self)
            self.stack.pop()
        return line

    def _line_define(self, line):
        """

        :param line:
        :return:
        """
        words = line.split()
        need_new_line = True
        if words[0] == "#ifdef":
            macro = words[1]
            if self.macro_define_stat and not self.macro_define_stat[-1]:
                self.macro_define_stat.append(False)
            elif macro in self.defined_macros.keys():
                self.macro_define_stat.append(True)
            else:
                self.macro_define_stat.append(False)
        elif words[0] == "#ifndef":
            macro = words[1]
            if self.macro_define_stat and not self.macro_define_stat[-1]:
                self.macro_define_stat.append(False)
            elif macro not in self.defined_macros.keys():
                self.macro_define_stat.append(True)
            else:
                self.macro_define_stat.append(False)
        elif words[0] == "#else":
            if len(self.macro_define_stat) <= 1 or self.macro_define_stat[-2]:
                self.macro_define_stat[-1] = not self.macro_define_stat[-1]
        elif words[0] == "#endif":
            self.macro_define_stat.pop()
        elif self.macro_define_stat and not self.macro_define_stat[-1]:
            self.stack.append(f"false define: {len(self.macro_define_stat)}")
            line = next(self)
            self.stack.pop()
            need_new_line = False
        elif words[0] == "#define":
            if len(words) > 2:
                self.defined_macros[words[1]] = line[line.find(words[1]) + len(words[1]):].strip()
            else:#elif len(words) > 1:
                self.defined_macros[words[1]] = ""
        elif words[0] == "#include":
            self._add_iterator_file(words[1])
        elif words[0] == "#undef":
            self.defined_macros.pop(words[1])
        elif words[0] == "#error":
            raise AssertionError(line)
        if need_new_line:
            self.stack.append(f"process define: {line}")
            line = next(self)
            self.stack.pop()
        return line


def _ffitp_dihedrals(line, buffers):
    """

    :param line:
    :param output:
    :return:
    """
    words = line.split()
    func = words[4]
    if func == "1":
        buffers["dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + " 0\n")
    elif func == "2":
        temp = [words[1], words[2], words[0], words[3]]
        temp2 = [words[1], words[2], words[3], words[0]]
        if words[0][0] == "O":
            temp = temp2
        elif words[3] == "C":
            temp = temp2
        elif words[3] == "CN3T":
            temp = temp2
        buffers["impropers"].append(
            "-".join(temp) + " {b} {k}".format(b=float(words[5]), k=float(words[6]) / 2) + "\n"
        )
    elif func == "3":
        buffers["RB_dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + "\n")
    elif func == "4":
        buffers["periodic_impropers"].append("-".join(words[:4]) + " " + " ".join(words[5:]) + "\n")
    elif func == "9":
        for i in range(5, len(words), 20):
            buffers["dihedrals"].append("-".join(words[:4]) + " " + " ".join(words[i:i + 3]) + " 0\n")
    else:
        raise NotImplementedError(f"Unsupported dihedral function type {func} for line:\n{line}")


def load_ffitp(filename, macros=None):
    """
    This **function** is used to load a fftip file

    .. ATTENTION::

        This is used to read a force field itp (ffitp) file instead of a simple itp file for a molecule.

    :param filename: the name of the file to load
    :param macros: the macros used to read the Gromacs topology file
    :return: a dict, which stores the name of the forcefield term - the corresponding information mapping
    """
    iterator = GromacsTopologyIterator(filename, macros)
    output = Xdict()
    buffers = {
        "nb14": ["name  kLJ  kee\n"],
        "atomtypes": ["name mass charge[e] LJtype\n"],
        "bonds": ["name b[nm] k[kJ/mol·nm^-2]\n"],
        "angles": ["name b[degree] k[kJ/mol·rad^-2]\n"],
        "Urey-Bradley": ["name b[degree] k[kJ/mol·rad^-2] r13[nm] kUB[kJ/mol·nm^-2]\n"],
        "dihedrals": ["name phi0[degree] k[kJ/mol] periodicity  reset\n"],
        "periodic_impropers": ["name phi0[degree] k[kJ/mol] periodicity\n"],
        "impropers": ["name phi0[degree] k[kJ/mol·rad^-2]\n"],
        "RB_dihedrals": ["name c0[kJ/mol] c1[kJ/mol] c2[kJ/mol] c3[kJ/mol] c4[kJ/mol] c5[kJ/mol]\n"],
    }
    output["cmaps"] = Xdict()
    output["bond_type_names"] = Xdict(not_found_message="The bond type of {} can not be found")
    for line in iterator:
        if iterator.flag == "":
            continue
        if iterator.flag == "defaults":
            words = line.split()
            assert int(words[0]) == 1, "SPONGE Only supports Lennard-Jones now"
            if int(words[1]) == 1:
                buffers["LJ"] = ["name A[kJ/mol·nm^6] B[kJ/mol·nm^12]\n"]
                buffers["nb14_extra"] = ["name A[kJ/mol·nm^6] B[kJ/mol·nm^12] kee\n"]
            else:
                buffers["LJ"] = ["name sigma[nm] epsilon[kJ/mol] \n"]
                buffers["nb14_extra"] = ["name sigma[nm] epsilon[kJ/mol] kee\n"]
            if len(words) > 3:
                fudge_lj = float(words[3])
            else:
                fudge_lj = 1
            if len(words) > 3:
                fudge_qq = float(words[4])
            else:
                fudge_qq = 1
            if len(words) > 2 and words[2] == "yes":
                buffers["nb14"].append("X-X {fudgeLJ} {fudgeQQ}\n".format(fudgeLJ=fudge_lj, fudgeQQ=fudge_qq))

        elif iterator.flag == "atomtypes":
            words = line.split()
            offset = 0
            if len(words) == 8:
                output["bond_type_names"][words[0]] = words.pop(1)
            elif len(words) == 6:
                offset = 1
            buffers["atomtypes"].append(
                "{type} {mass} {charge} {type}\n".format(
                    type=words[0], mass=float(words[2 - offset]), charge=float(words[3 - offset])
                )
            )
            buffers["LJ"].append(
                "{type}-{type} {V} {W}\n".format(type=words[0], V=float(words[5 - offset]), W=float(words[6 - offset]))
            )
        elif iterator.flag == "pairtypes":
            words = line.split()
            if len(words) <= 3:
                buffers["nb14"].append(
                    "{atom1}-{atom2} {kLJ} {kee}\n".format(
                        atom1=words[0], atom2=words[1], kLJ=fudge_lj, kee=fudge_qq
                    )
                )
            elif words[2] == "1":
                buffers["nb14_extra"].append(
                    "{atom1}-{atom2} {V} {W} {kee}\n".format(
                        atom1=words[0], atom2=words[1], V=float(words[3]), W=float(words[4]), kee=fudge_qq
                    )
                )
                buffers["nb14"].append("{atom1}-{atom2} 0 0\n".format(atom1=words[0], atom2=words[1]))
            elif words[2] == "2":
                raise NotImplementedError
        elif iterator.flag == "bondtypes":
            words = line.split()
            func = words[2]
            if func == "1":
                buffers["bonds"].append(
                    "{atom1}-{atom2} {b} {k}\n".format(
                        atom1=words[0], atom2=words[1], b=float(words[3]), k=float(words[4]) / 2
                    )
                )
            else:
                raise NotImplementedError
        elif iterator.flag == "angletypes":
            words = line.split()
            func = words[3]
            if func == "1":
                buffers["angles"].append(
                    "-".join(words[:3]) + " {b} {k}".format(b=float(words[4]), k=float(words[5]) / 2) + "\n"
                )
            elif func == "5":
                buffers["Urey-Bradley"].append(
                    "-".join(words[:3]) + " {b} {k} {b2} {k2}".format(
                        b=float(words[4]), k=float(words[5]) / 2, b2=float(words[6]), k2=float(words[7]) / 2
                    ) + "\n"
                )
            else:
                raise NotImplementedError
        elif iterator.flag == "dihedraltypes":
            _ffitp_dihedrals(line, buffers)
        elif iterator.flag == "cmaptypes":
            words = line.split()
            output["cmaps"]["-".join(words[:5])] = {"resolution": int(words[7]),
                                                    "parameters": list(map(lambda x: float(x) / 4.184, words[8:]))}
        elif iterator.flag == "nonbond_params":
            words = line.split()
            buffers["LJ"].append(
                "{type1}-{type2} {V} {W}\n".format(type1=words[0], type2=words[1], V=float(words[3]), W=float(words[4]))
            )
    for key, value in buffers.items():
        output[key] = "".join(value)
    return output

def _molitp_find_tail_residue(filename, macros, water_replace):
    """ find the index of the tail residue for every molecule """
    heads = {}
    tails = {}
    system_molecules = set()
    current = None
    iterator = GromacsTopologyIterator(filename, macros)
    for line in iterator:
        if iterator.flag == "moleculetype":
            current = line.split()[0]
            tails[current] = -999999
        elif iterator.flag == "atoms":
            resnr = int(line.split()[2])
            if resnr > tails[current]:
                tails[current] = resnr
            if current not in heads:
                heads[current] = resnr
        elif iterator.flag == "molecules":
            molname = line.split()[0]
            if not (molname == "SOL" and water_replace):
                system_molecules.add(molname)
    return heads, tails, system_molecules

def _molitp_parse_atoms(line, current_mol, current_stat, heads, tails, head_prefix, tail_prefix):
    """ parse one line for atoms of a molitp file """
    words = line.split()
    nr = int(words[0])
    resnr = int(words[2])
    resname = words[3]
    is_head = tails[current_mol.name] != heads[current_mol.name] and resnr == heads[current_mol.name]
    is_tail = tails[current_mol.name] != heads[current_mol.name] and resnr == tails[current_mol.name]
    if is_head:
        resname = head_prefix + resname
    if is_tail:
        resname = tail_prefix + resname
    atom_type = words[1]
    atom_name = words[4]
    if resnr != current_stat["resnr"]:
        current_stat["resnr"] = resnr
        current_stat["new_res_type"] = False
        if resname not in ResidueType.get_all_types():
            set_real_global_variable(resname, ResidueType(name=resname))
            current_stat["new_res_type"] = True
        current_mol.add_residue(Residue(ResidueType.get_type(resname)))
        if current_stat["new_res_type"]:
            current_stat[current_mol.residues[-1]] = True
    current_residue = current_mol.residues[-1]
    atom_type_name = atom_type if isinstance(atom_type, str) else atom_type.name
    expected_charge = float(words[6])
    existing_atom = current_residue.type._name2atom.get(atom_name)
    if existing_atom is not None:
        existing_charge = existing_atom.contents.get("charge")
        if existing_charge is None:
            existing_charge = existing_atom.contents.get("charge[e]")
        mismatch = (
            existing_atom.type.name != atom_type_name or
            (existing_charge is not None and abs(existing_charge - expected_charge) > 1e-6)
        )
        if mismatch:
            if is_head or is_tail:
                alt_base = resname
            else:
                Xprint(
                    f"ResidueType mismatch for {current_residue.type.name}.{atom_name} outside head/tail; "
                    f"auto-creating a new residue type.",
                    "WARNING",
                )
                alt_base = f"{current_residue.type.name}__{atom_name}_{atom_type_name}"
            compatible_type = None
            for name in ResidueType.get_all_types().keys():
                if name == alt_base or name.startswith(alt_base + "_"):
                    candidate = ResidueType.get_type(name)
                    candidate_atom = candidate._name2atom.get(atom_name)
                    if candidate_atom is None:
                        continue
                    if candidate_atom.type.name != atom_type_name:
                        continue
                    candidate_charge = candidate_atom.contents.get("charge")
                    if candidate_charge is None:
                        candidate_charge = candidate_atom.contents.get("charge[e]")
                    if candidate_charge is not None and abs(candidate_charge - expected_charge) > 1e-6:
                        continue
                    compatible_type = candidate
                    break
            if compatible_type is not None:
                current_residue.set_type(compatible_type, add_missing_atoms=False)
                current_stat[current_residue] = True
                existing_atom = current_residue.type._name2atom.get(atom_name)
            else:
                alt_resname = alt_base
                suffix = 1
                while alt_resname in ResidueType.get_all_types():
                    alt_resname = f"{alt_base}_{suffix}"
                    suffix += 1
                new_res_type = current_residue.type.deepcopy(alt_resname)
                set_real_global_variable(alt_resname, new_res_type)
                conflict_atom = new_res_type.name2atom(atom_name)
                conflict_atom.type = AtomType.get_type(atom_type_name)
                conflict_atom.Update(**{"charge[e]": expected_charge})
                current_residue.set_type(new_res_type, add_missing_atoms=False)
                current_stat[current_residue] = True
                existing_atom = current_residue.type._name2atom.get(atom_name)
    if current_stat["new_res_type"] or existing_atom is None:
        current_residue.type.Add_Atom(atom_name, atom_type, 0, 0, 0)
        current_residue.type.atoms[-1].Update(**{"charge[e]": expected_charge})
        current_stat[current_residue] = True
    current_residue.Add_Atom(atom_name, atom_type, 0, 0, 0)
    current_residue.atoms[-1].Update(**{"charge[e]": float(words[6])})
    current_stat[nr] = current_residue.atoms[-1]

def _molitp_parse_bonds(line, current_mol, current_stat):
    """ parse one line for bonds of a molitp file """
    words = line.split()
    func = int(words[2])
    atom1 = current_stat[int(words[0])]
    atom2 = current_stat[int(words[1])]
    if atom1.residue == atom2.residue:
        if current_stat[atom1.residue]:
            atom1.residue.type.add_connectivity(atom1.name, atom2.name)
        atom1.residue.add_connectivity(atom1.name, atom2.name)
    else:
        current_mol.add_residue_link(atom1, atom2)
    parser = GlobalSetting.get_gmx_bonded_type_parser("bond", func)
    parser(words, current_mol, current_stat)

def _molitp_build_system(current_stat, system):
    """ Directly build the system from the molecules """
    offset = 0
    system.get_atoms()
    for i in range(len(current_stat)):
        mol, n = current_stat[i]
        mol.built = True
        mol.get_atoms()
        for j in range(n):
            for force_type, forces in mol.bonded_forces.items():
                for force in forces:
                    new_atoms = [system.atoms[offset + mol.atom_index[atom]] for atom in force.atoms]
                    new_force = force.type.entity(new_atoms, force.type.getType("UNKNOWNS"))
                    new_force.contents = {**force.contents}
                    system.add_bonded_force(new_force)
            offset += len(mol.atoms)
    system.built = True

def load_molitp(filename, water_replace=True, head_prefix="N", tail_prefix="C", macros=None):
    """
    This **function** is used to load a molitp file

    .. ATTENTION::

        This is used to read a itp file for molecules (molitp or top) file instead of a force field file.

    :param filename: the name of the file to load
    :param water_replace: whether to change the SOL molecule in GMX to the water in Xponge. True as default.
    :param head_prefix: a string, the prefix will be added to the name of the first residue of each molecule
    :param tail_prefix: a string, the prefix will be added to the name of the last residue of each molecule
    :param macros: the macros used to read the Gromacs topology file
    :return: 1. an Xponge.Molecule representing the systema. None if not define
             2. an Xponge.Xdict, which maps the names of the molecules to Xponge.Molecule
    """
    sys = None
    mols = Xdict(not_found_message="{} is not a defined molecule in the system")
    current_mol = None
    current_stat = Xdict(not_found_method=lambda key: None)
    heads, tails, system_molecules = _molitp_find_tail_residue(filename, macros, water_replace)
    iterator = GromacsTopologyIterator(filename, macros)
    for line in iterator:
        Xprint(f"file={iterator.filenames[-1]}", "DEBUG")
        Xprint(f"flag={iterator.flag}", "DEBUG")
        Xprint(f"line={line}", "DEBUG")
        if iterator.flag == "moleculetype":
            name = line.split()[0]
            current_stat.clear()
            if water_replace and name == "SOL":
                mols[name] = ResidueType.get_type("WAT")
                current_stat["skip"] = True
            elif name in system_molecules:
                current_mol = Molecule(name=name)
                mols[name] = current_mol
            else:
                current_stat["skip"] = True
        elif iterator.flag == "atoms" and not current_stat["skip"]:
            _molitp_parse_atoms(line, current_mol, current_stat, heads, tails, head_prefix, tail_prefix)
        elif iterator.flag == "bonds" and not current_stat["skip"]:
            _molitp_parse_bonds(line, current_mol, current_stat)
        elif iterator.flag == "system":
            sys = Molecule(name=line.strip())
            current_stat.clear()
        elif iterator.flag == "molecules":
            words = line.split()
            mol = mols[words[0]]
            sys += mol * int(words[1])
            current_stat[len(current_stat)] = [mol, int(words[1])]
    _molitp_build_system(current_stat, sys)
    return sys, mols


##########################################################################
# CHARMM PSF Format
##########################################################################
_PSF_SECTION_RE = re.compile(r"^\s*(\d+)\s+!([A-Za-z0-9]+)")


def _psf_read_ints(fileobj, count):
    """Read ``count`` integers from a PSF section."""
    ints = []
    while len(ints) < count:
        line = fileobj.readline()
        if not line:
            raise EOFError("Unexpected EOF while reading PSF section")
        ints.extend([int(i) for i in line.split()])
    return ints


def _psf_get_atom_type(atom_type_name):
    """Get AtomType instance, creating a minimal one when missing."""
    try:
        return AtomType.get_type(atom_type_name)
    except KeyError:
        AtomType.New_From_String(f"name\n{atom_type_name}")
        return AtomType.get_type(atom_type_name)


def _psf_set_charge(atom, value):
    """Set charge if the AtomType exposes a charge field."""
    if "charge[e]" in atom.contents:
        atom.Update(**{"charge[e]": value})
    elif "charge" in atom.contents:
        atom.Update(**{"charge": value})


def _psf_get_charge(atom):
    """Get charge if present on the Atom instance."""
    if "charge[e]" in atom.contents:
        return atom.contents.get("charge[e]")
    if "charge" in atom.contents:
        return atom.contents.get("charge")
    return None


def _psf_parse_atom(line, current_mol, current_stat, residue_index_map, residue_newtype_map):
    """Parse one atom record in PSF and update molecule/residue state."""
    words = line.split()
    if len(words) < 8:
        raise ValueError(f"Bad PSF atom line: {line.rstrip()}")
    atom_index = int(words[0])
    segid = words[1]
    resnr = int(words[2])
    resname = words[3]
    atom_name = words[4]
    atom_type_name = words[5]
    charge = float(words[6])

    res_key = (segid, resnr, resname)
    if current_stat.get("res_key") != res_key:
        current_stat["res_key"] = res_key
        current_stat["new_res_type"] = False
        if resname not in ResidueType.get_all_types():
            set_real_global_variable(resname, ResidueType(name=resname))
            current_stat["new_res_type"] = True
        current_mol.add_residue(Residue(ResidueType.get_type(resname)))
        current_residue = current_mol.residues[-1]
        residue_index_map[current_residue] = len(current_mol.residues) - 1
        residue_newtype_map[current_residue] = current_stat["new_res_type"]
        current_stat[current_residue] = current_stat["new_res_type"]
        current_stat["resnr"] = resnr
    else:
        current_residue = current_mol.residues[-1]

    atom_type = _psf_get_atom_type(atom_type_name)
    expected_charge = charge
    existing_atom = current_residue.type._name2atom.get(atom_name)
    if existing_atom is not None:
        existing_charge = existing_atom.contents.get("charge")
        if existing_charge is None:
            existing_charge = existing_atom.contents.get("charge[e]")
        mismatch = (
            existing_atom.type.name != atom_type.name or
            (existing_charge is not None and abs(existing_charge - expected_charge) > 1e-6)
        )
        if mismatch:
            alt_base = current_residue.type.name
            compatible_type = None
            for name in ResidueType.get_all_types().keys():
                if name == alt_base or name.startswith(alt_base + "_"):
                    candidate = ResidueType.get_type(name)
                    candidate_atom = candidate._name2atom.get(atom_name)
                    if candidate_atom is None:
                        continue
                    if candidate_atom.type.name != atom_type.name:
                        continue
                    candidate_charge = candidate_atom.contents.get("charge")
                    if candidate_charge is None:
                        candidate_charge = candidate_atom.contents.get("charge[e]")
                    if candidate_charge is not None and abs(candidate_charge - expected_charge) > 1e-6:
                        continue
                    compatible_type = candidate
                    break
            if compatible_type is not None:
                current_residue.set_type(compatible_type, add_missing_atoms=False)
                residue_newtype_map[current_residue] = True
                current_stat[current_residue] = True
                existing_atom = current_residue.type._name2atom.get(atom_name)
            else:
                alt_resname = alt_base
                suffix = 1
                while alt_resname in ResidueType.get_all_types():
                    alt_resname = f"{alt_base}_{suffix}"
                    suffix += 1
                new_res_type = current_residue.type.deepcopy(alt_resname)
                set_real_global_variable(alt_resname, new_res_type)
                conflict_atom = new_res_type.name2atom(atom_name)
                conflict_atom.type = atom_type
                _psf_set_charge(conflict_atom, expected_charge)
                current_residue.set_type(new_res_type, add_missing_atoms=False)
                residue_newtype_map[current_residue] = True
                current_stat[current_residue] = True
                existing_atom = current_residue.type._name2atom.get(atom_name)

    if current_stat.get(current_residue) or existing_atom is None:
        current_residue.type.Add_Atom(atom_name, atom_type, 0, 0, 0)
        _psf_set_charge(current_residue.type.atoms[-1], expected_charge)
        current_stat[current_residue] = True
    current_residue.Add_Atom(atom_name, atom_type, 0, 0, 0)
    _psf_set_charge(current_residue.atoms[-1], expected_charge)
    current_stat[atom_index] = current_residue.atoms[-1]


def _psf_add_bond(atom1, atom2, current_mol, residue_index_map, residue_newtype_map):
    """Add bond connectivity from PSF indices."""
    if atom1.residue == atom2.residue:
        if residue_newtype_map.get(atom1.residue, False):
            atom1.residue.type.add_connectivity(atom1.name, atom2.name)
        atom1.residue.add_connectivity(atom1, atom2)
    else:
        current_mol.add_residue_link(atom1, atom2)
        index_diff = residue_index_map[atom1.residue] - residue_index_map[atom2.residue]
        if abs(index_diff) == 1:
            if residue_newtype_map.get(atom1.residue, False):
                if index_diff < 0:
                    atom1.residue.type.tail = atom1.name
                else:
                    atom1.residue.type.head = atom1.name
            if residue_newtype_map.get(atom2.residue, False):
                if index_diff < 0:
                    atom2.residue.type.head = atom2.name
                else:
                    atom2.residue.type.tail = atom2.name


def _psf_split_by_connectivity(current_mol):
    """Split molecule into connected components using atom connectivity."""
    adjacency = Xdict(not_found_method=lambda key: set())
    for residue in current_mol.residues:
        for atom, partners in residue.connectivity.items():
            for partner in partners:
                adjacency[atom].add(partner)
                adjacency[partner].add(atom)
    for link in current_mol.residue_links:
        adjacency[link.atom1].add(link.atom2)
        adjacency[link.atom2].add(link.atom1)

    components = {}
    visited = set()
    comp_id = 0
    for residue in current_mol.residues:
        for atom in residue.atoms:
            if atom in visited:
                continue
            comp_id += 1
            stack = [atom]
            components[atom] = comp_id
            visited.add(atom)
            while stack:
                cur = stack.pop()
                for nxt in adjacency[cur]:
                    if nxt in visited:
                        continue
                    visited.add(nxt)
                    components[nxt] = comp_id
                    stack.append(nxt)

    mols = Xdict()
    comp_to_mol = {}
    comp_to_atom_map = {}
    for residue in current_mol.residues:
        if not residue.atoms:
            continue
        comp = components.get(residue.atoms[0], None)
        if comp is None:
            continue
        if comp not in comp_to_mol:
            new_name = f"{current_mol.name}_{comp}"
            new_mol = Molecule(name=new_name)
            comp_to_mol[comp] = new_mol
            mols[new_name] = new_mol
            comp_to_atom_map[comp] = {}
        new_mol = comp_to_mol[comp]
        new_res = Residue(residue.type)
        for atom in residue.atoms:
            new_res.Add_Atom(atom.name, atom.type, atom.x, atom.y, atom.z)
            charge = _psf_get_charge(atom)
            if charge is not None:
                _psf_set_charge(new_res.atoms[-1], charge)
            comp_to_atom_map[comp][atom] = new_res.atoms[-1]
        for atom, partners in residue.connectivity.items():
            for partner in partners:
                new_res.Add_Connectivity(atom.name, partner.name)
        new_mol.Add_Residue(new_res)

    for link in current_mol.residue_links:
        comp = components.get(link.atom1, None)
        if comp is None:
            continue
        atom_map = comp_to_atom_map[comp]
        if link.atom1 in atom_map and link.atom2 in atom_map:
            comp_to_mol[comp].Add_Residue_Link(atom_map[link.atom1], atom_map[link.atom2])

    return mols


def _psf_split_residues_by_connectivity(current_mol):
    """Split residues containing multiple disconnected components (e.g., repeated resnr water)."""
    linked_atoms = set()
    for link in current_mol.residue_links:
        linked_atoms.add(link.atom1)
        linked_atoms.add(link.atom2)

    new_residues = []
    for residue in current_mol.residues:
        if not residue.atoms:
            new_residues.append(residue)
            continue
        if any(atom in linked_atoms for atom in residue.atoms):
            new_residues.append(residue)
            continue
        bond_count = sum(len(partners) for partners in residue.connectivity.values())
        if bond_count == 0:
            new_residues.append(residue)
            continue

        visited = set()
        components = []
        for atom in residue.atoms:
            if atom in visited:
                continue
            stack = [atom]
            visited.add(atom)
            comp = []
            while stack:
                cur = stack.pop()
                comp.append(cur)
                for nxt in residue.connectivity.get(cur, []):
                    if nxt in visited:
                        continue
                    visited.add(nxt)
                    stack.append(nxt)
            components.append(comp)

        if len(components) <= 1:
            new_residues.append(residue)
            continue

        for comp in components:
            new_res = Residue(residue.type)
            old_to_new = {}
            for atom in comp:
                new_res.Add_Atom(atom.name, atom.type, atom.x, atom.y, atom.z)
                new_atom = new_res.atoms[-1]
                new_atom.contents = dict(atom.contents.items())
                old_to_new[atom] = new_atom
            for atom in comp:
                for partner in residue.connectivity.get(atom, []):
                    if partner in old_to_new:
                        new_res.Add_Connectivity(old_to_new[atom].name, old_to_new[partner].name)
            new_residues.append(new_res)

    current_mol.residues = new_residues


def load_molpsf(filename, split_by="connectivity"):
    """
    This **function** is used to load a PSF file (topology only).

    :param filename: the name of the PSF file
    :param split_by: whether to split molecules ("connectivity" or None)
    :return: 1. a Molecule representing the system
             2. an Xdict mapping molecule names to Molecule
    """
    if not isinstance(filename, io.IOBase):
        fileobj = open(filename, "r")
    else:
        fileobj = filename

    with fileobj as f:
        header = f.readline()
        if not header.strip().startswith("PSF"):
            raise ValueError("Not a PSF file")

        mol_name = os.path.splitext(os.path.basename(getattr(f, "name", "psf")))[0] or "psf"
        current_mol = Molecule(name=mol_name)
        mols = Xdict()
        mols[mol_name] = current_mol

        current_stat = Xdict(not_found_method=lambda key: None)
        residue_index_map = {}
        residue_newtype_map = {}

        line = f.readline()
        while line:
            if not line.strip():
                line = f.readline()
                continue
            match = _PSF_SECTION_RE.match(line)
            if not match:
                line = f.readline()
                continue
            count = int(match.group(1))
            section = match.group(2).upper()

            if section == "NTITLE":
                for _ in range(count):
                    f.readline()
            elif section == "NATOM":
                for _ in range(count):
                    atom_line = f.readline()
                    if not atom_line:
                        raise EOFError("Unexpected EOF while reading NATOM section")
                    _psf_parse_atom(atom_line, current_mol, current_stat,
                                    residue_index_map, residue_newtype_map)
            elif section == "NBOND":
                bonds = _psf_read_ints(f, count * 2)
                for i in range(0, len(bonds), 2):
                    atom1 = current_stat.get(bonds[i])
                    atom2 = current_stat.get(bonds[i + 1])
                    if atom1 is None or atom2 is None:
                        continue
                    _psf_add_bond(atom1, atom2, current_mol,
                                  residue_index_map, residue_newtype_map)

            line = f.readline()

    _psf_split_residues_by_connectivity(current_mol)
    if split_by == "connectivity":
        mols = _psf_split_by_connectivity(current_mol)
    return current_mol, mols


def load_gro(filename, mol=None, read_box_angle=True):
    """
    This **function** is used to read the GROMACS coordinate file

    :param filename: the coordinate file to load
    :param mol: the molecule or residue to load the coordinate into
    :param read_box_angle: read the box angle information when provided
    :return: two numpy arrays, representing the coordinates and the box information respectively
    """
    with Xopen(filename, "r") as f:
        f.readline()
        atom_numbers = int(f.readline().split()[0])
        crd = np.zeros((atom_numbers, 3), dtype=np.float32)
        for i in range(atom_numbers):
            line = f.readline()
            crd[i][0] = float(line[20:28]) * 10
            crd[i][1] = float(line[28:36]) * 10
            crd[i][2] = float(line[36:44]) * 10
        line = f.readline().split()
        if len(line) == 3:
            box_length = [float(line[0]) * 10, float(line[1]) * 10, float(line[2]) * 10]
            box_angle = [90.0, 90.0, 90.0]
        elif len(line) == 9:
            vals = list(map(float, line))
            # GROMACS triclinic box: v1x v2y v3z v1y v1z v2x v2z v3x v3y
            v1 = [vals[0], vals[3], vals[4]]
            v2 = [vals[5], vals[1], vals[6]]
            v3 = [vals[7], vals[8], vals[2]]
            box_length, box_angle = get_length_angle_from_basis_vectors(v1, v2, v3, angle_in_degree=True)
            box_length = [i * 10 for i in box_length]
        else:
            raise NotImplementedError("Unsupported GRO box format")
        if not read_box_angle:
            box_angle = [90.0, 90.0, 90.0]
        box = np.array([box_length[0], box_length[1], box_length[2],
                        box_angle[0], box_angle[1], box_angle[2]], dtype=float)
    if mol:
        for i, atom in enumerate(mol.get_atoms()):
            atom.x = crd[i][0]
            atom.y = crd[i][1]
            atom.z = crd[i][2]
        if isinstance(mol, Molecule):
            mol.box_length = box[:3].tolist()
            mol.box_angle = box[3:].tolist()
    return crd, box

set_global_alternative_names()
