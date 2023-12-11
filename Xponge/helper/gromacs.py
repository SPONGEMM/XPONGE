"""
    This **module** contains functions to interact with gromacs
"""
import os
from .. import Molecule, set_global_alternative_names, get_assignment_from_mol2
from . import Xdict, Xprint

def _check_new_residue_when_sorting(mol, mol_res_id,
  sort_map, cri, ri, rname):
    """ check whether it is a new residue """
    if ri != cri:
        cri = ri
        if mol_res_id >= 0:
            residue = mol.residues[mol_res_id]
            restype = residue.type
            if "O1P" in sort_map:
                sort_map["OP1"] = sort_map["O1P"]
            if "O2P" in sort_map:
                sort_map["OP2"] = sort_map["O2P"]
            if "C5M" in sort_map and restype.name.startswith("DT"):
                sort_map["C7"] = sort_map["C5M"]
            if len(restype.name) == 4 and restype.name[0] == "C" and "OT1" in sort_map and "OT2" in sort_map:
                sort_map["O"] = sort_map.pop("OT1")
                sort_map["OXT"] = sort_map.pop("OT2")
            if len(restype.name) == 4 and restype.name[0] == "C" and "O1" in sort_map and "O2" in sort_map:
                sort_map["O"] = sort_map.pop("O1")
                sort_map["OXT"] = sort_map.pop("O2")
            for atom in residue.atoms:
                if atom.name.startswith("H") or atom.name in ("1H", "2H", "3H"):
                    atom_connected = next(iter(restype.connectivity[restype.name2atom(atom.name)]))
                    find_name = atom_connected.name if atom_connected.name in sort_map else atom_connected.name[:-1]
                    if atom_connected.name not in sort_map:
                        Xprint(f"{atom_connected.name} not found in sort_map", "WARNING")
                    sort_map[atom.name] = sort_map[find_name] + 0.5
            residue.atoms.sort(key=lambda atom: sort_map[atom.name] if atom.name in sort_map else sort_map[atom.name[:-1]])
        sort_map.clear()
        sort_map.not_found_message = "The name of the atom '{}' in the %d-th residue in \
Xponge.Molecule can not be found in the gro file (%5d%-5s)"%(mol_res_id + 1, ri, rname)
        mol_res_id += 1
    return mol_res_id, cri

def sort_atoms_by_gro(mol, gro):
    """
       This **function** sorts the atoms in a Xponge.Molecule according to the index in a gro file

       :param mol: a Xponge.Molecule
       :param gro: the gro file
    """
    if not isinstance(mol, Molecule):
        raise TypeError("The input for sorting should be an Xponge.Molecule")
    mol.get_atoms()
    current_residue_index = 0
    mol_res_id = -1
    current_sort_map = Xdict()
    with open(gro) as f:
        for li, line in enumerate(f):
            if li == 1:
                atom_numbers = int(line)
                if atom_numbers != len(mol.atoms):
                    raise ValueError("The number of atoms is not equal \
in the gro file and in the Xponge.Molecule instance")
            elif 1 < li < len(mol.atoms) + 2:
                residue_index = int(line[:5])
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                if residue_name in ("WAT", "H2O", "HOH", "SOL"):
                    break
                mol_res_id, current_residue_index = _check_new_residue_when_sorting(mol, mol_res_id,
                   current_sort_map, current_residue_index, residue_index, residue_name)
                if atom_name.startswith("H") or atom_name in ("1H", "2H", "3H"):
                    continue
                current_sort_map[atom_name] = len(current_sort_map)
    _check_new_residue_when_sorting(mol, mol_res_id,
       current_sort_map, current_residue_index, current_residue_index + 1, "")

def _output_mol2(prefix, flag, ans, atoms, bond_numbers):
    """ output 1 rtp residue to 1 mol2 file """
    if flag not in ("bondedtypes", ""):
        ans = ans%(len(atoms), bond_numbers)
        atoms.clear()
        with open(os.path.join(prefix, f"{flag}.mol2"), "w") as f:
            f.write(ans)
    return 0


def rtp_to_mol2(rtp, prefix):
    """
        This **function** convert a rtp file into several mol2 files

        :param rtp: the name of the rtp file
        :param prefix: the output prefix for the mol files
    """
    with open(rtp) as f:
        flag = ""
        small_flag = ""
        ans = ""
        bond_numbers = 0
        atoms = {}
        for line in f:
            line = line.split(";")[0]
            if not line.strip():
                continue
            if line.strip().startswith("["):
                small_flag = line.split("[")[1].split("]")[0].strip()
                if small_flag not in ("bondedtypes", "atoms", "bonds", "impropers", "cmap"):
                    bond_numbers = _output_mol2(prefix, flag, ans, atoms, bond_numbers)
                    flag = small_flag
                    ans = f"""@<TRIPOS>MOLECULE
{flag}
 %d  %d  1  0  1 
SMALL
USER_CHARGES
@<TRIPOS>ATOM\n"""
                if small_flag == "bonds":
                    ans += "@<TRIPOS>BOND\n"
            elif flag == "bondedtypes":
                continue
            elif small_flag in ("impropers",):
                continue
            elif small_flag == "atoms":
                words = line.split()
                atoms[words[0]] = len(atoms) + 1
                ans += f"{len(atoms)} {words[0]} 0 0 0 {words[1]} 1 {flag} {words[2]}\n"
            elif small_flag == "bonds":
                words = line.split()
                if {words[0][0], words[1][0]} & {"-", "+"}:
                    continue
                bond_numbers += 1
                ans += f"{bond_numbers} {atoms[words[0]]} {atoms[words[1]]} 1\n"
        _output_mol2(prefix, flag, ans, atoms, bond_numbers)

def read_tdb(tdb, rule, mol2_in, mol2_out, newname):
    """
        This **function** reads a tdb file and according to the rule modify the mol2 file

        :param tdb: the name of the tdb file
        :param rule: the rule to use
        :param mol2_in: the input mol2 file
        :param mol2_out: the output mol2 file
        :param newname: the new name of the residue
    """
    assign = get_assignment_from_mol2(mol2_in)
    atom_names = {atom: i for i, atom in enumerate(assign.names)}
    with open(tdb) as f:
        flag = ""
        small_flag = ""
        for line in f:
            line = line.split(";")[0]
            if not line.strip():
                continue
            if line.strip().startswith("["):
                small_flag = line.split("[")[1].split("]")[0].strip()
                if small_flag.lower() not in ("None", "delete", "replace", "add"):
                    flag = small_flag
                small_flag = small_flag.lower()
            elif flag != rule:
                continue
            elif small_flag == "delete":
                if line.strip() in atom_names:
                    i = atom_names.pop(line.strip())
                    assign.delete_atom(i)
                    for name, value in atom_names.items():
                        if value > i:
                            atom_names[name] -= 1
            elif small_flag == "replace":
                words = line.split()
                if words[0] in atom_names:
                    i = atom_names.pop(words[0])
                    if len(words) == 5:
                        assign.names[i] = words[1]
                        assign.atoms[i] = words[2]
                        assign.charge[i] = float(words[4])
                        atom_names[words[1]] = i
                    else:
                        assign.atoms[i] = words[1]
                        assign.charge[i] = float(words[3])
                        atom_names[words[0]] = i
            elif small_flag == "add":
                words = line.split()
                words2 = f.readline().split()
                if words[3] in atom_names:
                    i = atom_names[words[3]]
                    for j in range(int(words[0])):
                        if words[2] + str(j + 1) not in atom_names:
                            assign.add_atom(words2[0], 0, 0, 0, words[2] + str(j + 1), float(words2[2]))
                            atom_names[words[2] + str(j + 1)] = len(assign.atoms) - 1
                            assign.add_bond(len(assign.atoms) - 1, i, -1)

    assign.determine_ring_and_bond_type()
    assign.name = newname
    assign.save_as_mol2(mol2_out, atomtype=None)

#pylint: disable=line-too-long, pointless-string-statement
"""
===========================================
Convert protein templates in GROMACS to SPONGE
Using charmm36-jul2022 as example
===========================================
#0. read ffitp files
#-----------------------
#1. get good residues
#------------------------
import os
from pathlib import Path
import Xponge
import Xponge.forcefield.amber.ff14sb
from Xponge.helper.gromacs import rtp_to_mol2, read_tdb
DIR = Path("/usr/local/gromacs/share/gromacs/top/charmm36-jul2022.ff")
rtp_to_mol2(DIR / "aminoacids.rtp", "./")
#name map comes from .r2b file
namemap = {"CYX": "CYS2", "HID": "HSD", "HIE": "HSE", "HIP": "HSP"}
alls = "ACE ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL NME"
for i in  alls.split():
    ii = namemap.get(i, i)
    res = Xponge.ResidueType.get_type(i)
    Xponge.save_mol2(res, f"amber_{i}.mol2")
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {ii}.mol2 \
-tformat gaff_mol2 -tfile amber_{i}.mol2 -oformat gaff_mol2 -ofile good_{i}.mol2 -ff Xponge.forcefield.amber.ff14sb  Xponge.forcefield.charmm.charmm36 -cpcrd -ores {i}")
    if i in ("ACE", "NME"):
        continue
#2. get good 5 residues
#------------------------
    j = "N" + i
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "aminoacids.n.tdb", "NH3+", f"{ii}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.ff14sb  Xponge.forcefield.charmm.charmm36 -cpcrd -ores {j}")
#3. get good 3 residues
#------------------------
    j = "C" + i
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "aminoacids.c.tdb", "COO-", f"{ii}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.ff14sb  Xponge.forcefield.charmm.charmm36 -cpcrd -ores {j}")

#4. load all templates
#---------------------------
Xponge.ResidueType.clear_type()
import Xponge.forcefield.charmm.charmm36
templates = None
for i in alls.split():
    for j in " NC":
        k = j + i
        k = k.strip()
        t = load_mol2(f"good_{k}.mol2")
        if templates is not None:
            templates |= t
        else:
            templates = t
        if i in ("ACE", "NME"):
            break

Xponge.save_mol2(templates, "dna.mol2")

#5. define the linkage
#----------------------- 
"""

"""
===========================================
Convert DNA templates in GROMACS to SPONGE
Using charmm27 as example
===========================================
#0. read ffitp files
#-----------------------
#1. get good DA DT DC DG
#------------------------
import os
from pathlib import Path
import Xponge
import Xponge.forcefield.amber.bsc1
from Xponge.helper.gromacs import rtp_to_mol2, read_tdb
DIR = Path("/usr/local/gromacs/share/gromacs/top/charmm27.ff")
rtp_to_mol2(DIR / "dna.rtp", "./")
for i in "ATCG":
    i = "D" + i
    res = Xponge.ResidueType.get_type(i)
    Xponge.save_mol2(res, f"amber_{i}.mol2")
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {i}.mol2 \
-tformat gaff_mol2 -tfile amber_{i}.mol2 -oformat gaff_mol2 -ofile good_{i}.mol2 -ff Xponge.forcefield.amber.bsc1  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {i}")
#2. get good 5 residues
#------------------------
    j = i + "5"
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "dna.n.tdb", "5'", f"{i}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.bsc1  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {j}")
#3. get good 3 residues
#------------------------
    j = i + "3"
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "dna.c.tdb", "3'", f"{i}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.bsc1  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {j}")

#4. load all templates
#---------------------------
Xponge.ResidueType.clear_type()
import Xponge.forcefield.charmm.charmm27
templates = None
for i in "ATCG":
    for j in " 35":
        k = "D" + i + j
        k = k.strip()
        t = load_mol2(f"good_{k}.mol2")
        if templates is not None:
            templates |= t
        else:
            templates = t
Xponge.save_mol2(templates, "dna.mol2")

#5. define the linkage
#----------------------- 
"""

"""
===========================================
Convert RNA templates in GROMACS to SPONGE
Using charmm27 as example
===========================================
#0. read ffitp files
#-----------------------
#1. get good A U C G
#------------------------
import os
from pathlib import Path
import Xponge
import Xponge.forcefield.amber.ol3
from Xponge.helper.gromacs import rtp_to_mol2, read_tdb
DIR = Path("/usr/local/gromacs/share/gromacs/top/charmm27.ff")
rtp_to_mol2(DIR / "rna.rtp", "./")
for i in "AUCG":
    res = Xponge.ResidueType.get_type(i)
    Xponge.save_mol2(res, f"amber_{i}.mol2")
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile R{i}.mol2 \
-tformat gaff_mol2 -tfile amber_{i}.mol2 -oformat gaff_mol2 -ofile good_{i}.mol2 -ff Xponge.forcefield.amber.ol3  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {i}")
#2. get good 5 residues
#------------------------
    j = i + "5"
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "rna.n.tdb", "5'", f"{i}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.ol3  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {j}")
#3. get good 3 residues
#------------------------
    j = i + "3"
    res = Xponge.ResidueType.get_type(j)
    Xponge.save_mol2(res, f"amber_{j}.mol2")
    read_tdb(DIR / "rna.c.tdb", "3'", f"{i}.mol2", f"{j}.mol2", j)
    os.system(f"Xponge name2name -fformat gaff_mol2 -ffile {j}.mol2 \
-tformat gaff_mol2 -tfile amber_{j}.mol2 -oformat gaff_mol2 -ofile good_{j}.mol2 -ff Xponge.forcefield.amber.ol3  Xponge.forcefield.charmm.charmm27 -cpcrd -ores {j}")

#4. load all templates
#---------------------------
Xponge.ResidueType.clear_type()
import Xponge.forcefield.charmm.charmm27
templates = None
for i in "ATCG":
    for j in " 35":
        k = "D" + i + j
        k = k.strip()
        t = load_mol2(f"good_{k}.mol2")
        if templates is not None:
            templates |= t
        else:
            templates = t
Xponge.save_mol2(templates, "dna.mol2")

#5. define the linkage
#----------------------- 
"""

set_global_alternative_names()
