"""
This **module** implements the terminal commands
"""
import os
import shutil
import sys

from ..helper import source, GlobalSetting, Xopen, Xprint
from ..mdrun import run
from .unittests import mytest

def converter(args):
    """
    This **function** converts the format of coordinate file

    :param args: arguments from argparse
    :return: None
    """
    from ..analysis import md_analysis as xmda #pylint:disable=unused-import
    import MDAnalysis as mda

    if args.c:
        if args.cf == "guess":
            u = mda.Universe(args.p, args.c)
        elif args.cf == "sponge_crd":
            u = mda.Universe(args.p, args.c, format=xmda.SpongeCoordinateReader)
        elif args.cf == "sponge_traj":
            dirname, basename = os.path.split(args.c)
            if basename == "mdcrd.dat":
                box = "mdbox.txt"
            else:
                box = basename.replace(".dat", ".box")
            box = os.path.join(dirname, box)
            if box and os.path.exists(box):
                u = mda.Universe(args.p, args.c, box=box, format=xmda.SpongeTrajectoryReader)
            else:
                u = mda.Universe(args.p, args.c, format=xmda.SpongeTrajectoryReader)
    else:
        u = mda.Universe(args.p)

    if args.of == "sponge_crd":
        with xmda.SpongeCoordinateWriter(args.o) as w:
            w.write(u)
    elif args.of == "sponge_traj":
        with xmda.SpongeTrajectoryWriter(args.o) as w:
            for _ in u.trajectory:
                w.write(u)
    else:
        with mda.Writer(args.o, n_atoms=len(u.coord.positions)) as w:
            for _ in u.trajectory:
                w.write(u)


def maskgen(args):
    """
    This **function** uses VMD to generate mask

    :param args: arguments from argparse
    :return: None
    """
    import MDAnalysis as mda
    import Xponge.analysis.md_analysis as xmda #pylint: disable=unused-import
    if not os.path.exists(args.p):
        raise FileNotFoundError(f"can not find {args.p}")
    kwargs = {}
    if args.pf is not None:
        kwargs["topology_format"] = args.pf
    if args.cf is not None:
        kwargs["format"] = args.cf
    if args.c is not None:
        if not os.path.exists(args.c):
            raise FileNotFoundError(f"can not find {args.c}")
        u = mda.Universe(args.p, args.c, **kwargs)
    else:
        u = mda.Universe(args.p, **kwargs)
    id2index_map = {atom.id: str(i) for i, atom in enumerate(u.atoms)}
    atoms = u.select_atoms(args.s)
    with open(args.o, "w") as f:
        f.write("\n".join([id2index_map[atom.id] for atom in atoms]))
    if args.oc is not None:
        with open(args.oc, "w") as f:
            f.write("\n".join([f"{atom.position[0]} {atom.position[1]} {atom.position[2]}" for atom in atoms]))


def exgen(args):
    """
    This **function** reads the SPONGE input files for bonded interactions and generate a exclude file

    :param args: arguments from argparse
    :return: None
    """
    partners = [set([]) for i in range(args.n)]

    def exclude_2_atoms(words):
        i, j = int(words[0]), int(words[1])
        partners[i].add(j)
        partners[j].add(i)

    def exclude_3_atoms(words):
        i, k = int(words[0]), int(words[2])
        partners[i].add(k)
        partners[k].add(i)

    def exclude_4_atoms(words):
        i, l = int(words[0]), int(words[3])
        partners[i].add(l)
        partners[l].add(i)

    for bond in args.bond:
        with open(bond) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_2_atoms(words)

    for angle in args.angle:
        with open(angle) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_3_atoms(words)
    for dihedral in args.dihedral:
        with open(dihedral) as f:
            f.readline()
            for line in f:
                words = line.split()
                exclude_4_atoms(words)

    for virtual in args.virtual:
        with open(virtual) as f:
            for line in f:
                words = line.split()
                t = int(words[0])
                if t == 0:
                    exclude_2_atoms(words[1:])
                elif t == 1:
                    exclude_3_atoms(words[1:])
                elif t in (2, 3):
                    exclude_4_atoms(words[1:])
                else:
                    raise TypeError("virtual atom type wrong: are you sure this is a SPONGE virtual atom file?")

    for exclude in args.exclude:
        with open(exclude) as f:
            f.readline()
            count = 0
            for line in f:
                words = line.split()
                t = set(words[1:])
                partners[count] = partners[count].union(t)
                count += 1

    total = 0
    towrite = "{} {}\n"
    for i, p in enumerate(partners):
        newp = []
        for pi in p:
            if pi > i:
                newp.append(pi)
        towrite += "%d " % len(newp)
        towrite += ("{} " * len(newp)).format(*newp) + "\n"
        total += len(newp)
        towrite = towrite.format(args.n, total)

    f = Xopen(args.o, "w")
    f.write(towrite)
    f.close()


def name2name(args):
    """
    This **function** change the atom names from one file to another file

    :param args: arguments from argparse
    :return: None
    """
    from rdkit.Chem import rdFMCS
    source("..")
    rdktool = source("..helper.rdkit")

    if args.to_format == "mol2":
        to_ = assign.Get_Assignment_From_Mol2(args.to_file, total_charge="sum")
    elif args.to_format == "gaff_mol2":
        source("..forcefield.amber.gaff")
        to_ = load_mol2(args.to_file).residues[0]
        to_ = assign.Get_Assignment_From_ResidueType(to_)
    elif args.to_format == "pdb":
        to_ = assign.Get_Assignment_From_PDB(args.to_file,
                                             only_residue=args.to_residue)

    if args.from_format == "mol2":
        from_ = assign.Get_Assignment_From_Mol2(args.from_file, total_charge="sum")
    elif args.from_format == "gaff_mol2":
        source("..forcefield.amber.gaff")
        from_ = load_mol2(args.from_file).residues[0]
        from_ = assign.Get_Assignment_From_ResidueType(from_)
    elif args.from_format == "pdb":
        from_ = assign.Get_Assignment_From_PDB(args.from_file,
                                               only_residue=args.from_residue)

    from_.add_index_to_name()
    rdmol_a = rdktool.assign_to_rdmol(to_, ignore_bond_type=True)
    rdmol_b = rdktool.assign_to_rdmol(from_, ignore_bond_type=True)

    result = rdFMCS.FindMCS([rdmol_b, rdmol_a], timeout=args.tmcs)

    match_a = rdmol_a.GetSubstructMatch(result.queryMol)
    match_b = rdmol_b.GetSubstructMatch(result.queryMol)
    matchmap = {from_.names[match_b[j]]: to_.names[match_a[j]] for j in range(len(match_a))}

    from_.names = [matchmap.get(name, name) for name in from_.names]
    from_.name = args.out_residue

    if args.out_format == "mol2":
        from_.Save_As_Mol2(args.out_file)
    elif args.out_format == "pdb":
        from_.Save_As_PDB(args.out_file)
    elif args.out_format == "mcs_pdb":
        towrite = towrite = "REMARK   Generated By Xponge (Max Common Structure)\n"
        for i, atom in enumerate(from_.atoms):
            if i in match_b:
                towrite += "ATOM  %5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%17s%2s\n" % (i + 1, from_.names[i],
                                                                                     from_.name, " ", 1,
                                                                                     from_.coordinate[i][0],
                                                                                     from_.coordinate[i][1],
                                                                                     from_.coordinate[i][2], " ", atom)
        f = Xopen(args.out_file, "w")
        f.write(towrite)
        f.close()


def mol2rfe(args):
    """
    This **function** helps with the relative free energy calculation

    :param args: arguments from argparse
    :return: None
    """
    from .mol2rfe import _mol2rfe_build, _mol2rfe_min, _mol2rfe_pre_equilibrium, \
 _mol2rfe_equilibrium, _mol2rfe_analysis

    source("..")
    source("..forcefield.special.fep")
    source("..forcefield.special.min")

    if args.fl is not None:
        l = list(np.loadtxt(args.fl))
        args.l = sorted(list(set(l)))
    elif args.l is not None:
        args.l = sorted(list(set(args.l)))

    if args.l is not None:
        for li in args.l:
            if li > 1 or li < 0:
                Xprint(f"There is a weird lambda == {li}", "WARNING")
        if args.l[-1] != 1:
            Xprint(f"The largest lambda {args.l[-1]} != 1", "WARNING")
        if args.l[0] != 0:
            Xprint(f"The smallest lambda {args.l[0]} != 0", "WARNING")
        args.nl = len(args.l) - 1
    else:
        args.l = np.linspace(0.0, 1.0, args.nl + 1)

    if not args.ff:
        source("..forcefield.amber.gaff")
        source("..forcefield.amber.ff14sb")
        source("..forcefield.amber.tip3p")
    else:
        import_python_script(args.ff)

    if not args.do:
        args.do = [["build", "min", "pre_equilibrium", "equilibrium", "analysis"]]
    args.do = args.do[0]
    if "debug" in args.do:
        Debug()
    from_res_type_ = load_mol2(args.r1).residues[0].type
    from_ = assign.Get_Assignment_From_ResidueType(from_res_type_)
    if not args.ff:
        parmchk2_gaff(args.r1, args.temp + "_TMP1.frcmod")

    to_res = load_mol2(args.r2).residues[0]
    to_ = assign.Get_Assignment_From_ResidueType(to_res.type)
    if not args.ff:
        parmchk2_gaff(args.r2, args.temp + "_TMP2.frcmod")

    for extrai, mol2file in enumerate(args.r0):
        load_mol2(mol2file)
        parmchk2_gaff(mol2file, f"{args.temp}_TMP{3+extrai}.frcmod")

    rmol = load_pdb(args.pdb)

    merged_from, merged_to, matchmap = Merge_Dual_Topology(rmol, rmol.residues[args.ri],
                                                    to_res, from_, to_,
                                                    args.tmcs, f"{args.fmcs}",
                                                    args.lmcs, args.imcs)

    if args.dohmr:
        H_Mass_Repartition(merged_from)
        H_Mass_Repartition(merged_to)

    _mol2rfe_build(args, merged_from, merged_to)

    _mol2rfe_min(args, [0])

    _mol2rfe_pre_equilibrium(args, [0])

    _mol2rfe_min(args, range(1, args.nl + 1))

    _mol2rfe_pre_equilibrium(args, range(1, args.nl + 1))

    _mol2rfe_equilibrium(args)

    _mol2rfe_analysis(args, merged_from, merged_to, matchmap, from_, to_)

def mm_gbsa(args):
    """
    This **function** helps with the MM/GBSA calculation

    :param args: arguments from argparse
    :return: None
    """
    from .mmgbsa import _mmgbsa_build, _mmgbsa_min, _mmgbsa_pre_equilibrium, \
_mmgbsa_equilibrium, _mmgbsa_analysis

    if not args.do:
        args.do = [["build", "min", "pre_equilibrium", "equilibrium", "analysis"]]
    args.do = args.do[0]
    if "debug" in args.do:
        Debug()

    _mmgbsa_build(args)

    _mmgbsa_min(args)

    _mmgbsa_pre_equilibrium(args)

    _mmgbsa_equilibrium(args)

    _mmgbsa_analysis(args)
