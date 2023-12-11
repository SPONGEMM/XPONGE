"""
    This **module** includes unittests of the Xponge.forcefield.charmm36
"""
import os
from pathlib import Path

__all__ = ["test_protein", "test_rna", "test_dna"]

def _get_energies(filename):
    """ get the energies from the log file """
    import numpy as np
    gmx_mdout = {}
    keywords = []
    with open(filename) as f:
        t = f.read()
    t = t.split("   Energies (kJ/mol)")[1:]
    stop = False
    for tt in t:
        count = 0
        if stop:
            break
        for line in tt.split("\n"):
            if line.strip().startswith("<=="):
                stop = True
                break
            start = 0
            word = line[start:start+15].strip()
            while word:
                start += 15
                try:
                    word = float(word)
                except ValueError:
                    pass
                if isinstance(word, float):
                    gmx_mdout[keywords[count]].append(word)
                    count += 1
                else:
                    if word not in gmx_mdout:
                        gmx_mdout[word] = []
                        keywords.append(word)
                word = line[start:start+15].strip()
    for key, value in gmx_mdout.items():
        gmx_mdout[key] = np.array(value)
    return gmx_mdout

def _check_one_energy(gmx_mdout, term_name, sponge_out):
    """check one energy term"""
    import matplotlib.pyplot as plt
    import numpy as np
    unit = "kcal/mol"
    gmx_mdout /= 4.184
    if abs(np.mean(sponge_out)) > 1000:
        unit = "Mcal/mol"
        gmx_mdout /= 1000
        sponge_out /= 1000
    k, b = np.polyfit(gmx_mdout, sponge_out, 1)
    r = np.corrcoef(gmx_mdout, sponge_out)[0][1]
    plt.plot([np.min(gmx_mdout), np.max(gmx_mdout)],
             [k * np.min(gmx_mdout) + b, k * np.max(gmx_mdout) + b],
             label=f"y={k:.3f}x{b:+.3f},r={r:.3f}")
    plt.plot(gmx_mdout, sponge_out, "o")
    plt.xlabel(f"Result from GROMACS [{unit}]")
    plt.ylabel(f"Result from SPONGE [{unit}]")
    plt.legend()
    plt.savefig(f"{term_name}.png")
    plt.clf()


def test_protein():
    """
        test the single point energy for residues of protein
    """
    import Xponge
    import Xponge.forcefield.charmm.charmm36
    import Xponge.forcefield.charmm.tip3p_charmm
    from Xponge.helper.gromacs import Sort_Atoms_By_Gro
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run

    s = "ALA ARG ASN ASP CYS GLN GLU GLY ILE LEU " + \
        "LYS MET PHE PRO SER THR TRP TYR VAL"

    mol = Xponge.ResidueType.get_type("NALA")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("CALA")

    Path("protein").mkdir(exist_ok=True)
    Xponge.save_pdb(mol, "protein/protein.pdb")
    assert os.system("cd protein && gmx pdb2gmx -f protein.pdb -ff charmm36-jul2022 \
-water tip3p -ignh > pdb2gmx.log 2>&1") == 0

    Sort_Atoms_By_Gro(mol, "protein/conf.gro")
    Xponge.save_sponge_input(mol, "protein/protein")

    assert run("SPONGE -mode minimization -default_in_file_prefix protein/protein -step_limit 2000 \
-rst protein/min > protein/min.log") == 0
    Xponge.load_coordinate("protein/min_coordinate.txt", mol)
    Xponge.save_gro(mol, "protein/min.gro")

    with open("protein/protein.mdp", "w") as f:
        f.write("""integrator = sd
dt = 2e-3
constraint_algorithm = lincs
constraints = h-bonds
nsteps = 1000
nstxout = 1
nstlog = 1
coulombtype = PME
vdw-modifier = None
ref_t = 300
tau_t = 1
nstlist = 1
tc_grps = system
DispCorr = Ener
""")
    assert os.system("cd protein && gmx grompp -f protein.mdp -c min.gro -p topol.top \
-o run.tpr -maxwarn 1 > runpp.log 2>&1") == 0
    assert os.system("cd protein && gmx mdrun -deffnm run > mdrun.log 2>&1") == 0

    assert os.system("cd protein && Xponge converter -p protein_mass.txt -c run.trr -o run.dat > convert.log") == 0
    assert run("SPONGE -mode rerun -default_in_file_prefix protein/protein -crd protein/run.dat \
-mdout protein/mdout.txt -box protein/run.box > protein/rerun.log") == 0

    gmx_mdout = _get_energies("protein/run.log")
    t = MdoutReader("protein/mdout.txt")

    Path("protein_plot").mkdir(exist_ok=True)
    _check_one_energy(gmx_mdout["Potential"], "protein_plot/Potential", t.potential)
    _check_one_energy(gmx_mdout["Bond"], "protein_plot/Bond", t.bond)
    _check_one_energy(gmx_mdout["U-B"], "protein_plot/Urey_Bradly", t.urey_bradley)
    _check_one_energy(gmx_mdout["Proper Dih."], "protein_plot/Dihedral", t.dihedral)
    _check_one_energy(gmx_mdout["Improper Dih."], "protein_plot/Improper", t.improper_dihedral)
    _check_one_energy(gmx_mdout["CMAP Dih."], "protein_plot/Cmap", t.cmap)
    _check_one_energy(gmx_mdout["LJ-14"], "protein_plot/1-4 NB", t.nb14_LJ)
    _check_one_energy(gmx_mdout["Coulomb-14"], "protein_plot/1-4 EEL", t.nb14_EE)
    _check_one_energy(gmx_mdout["LJ (SR)"] + gmx_mdout["Disper. corr."], "protein_plot/LJ", t.LJ)
    _check_one_energy(gmx_mdout["Coulomb (SR)"] + gmx_mdout["Coul. recip."], "protein_plot/PME", t.PME)

def test_rna():
    """
        test the single point energy for residues of RNA
    """
    import Xponge
    import Xponge.forcefield.charmm.charmm36
    from Xponge.helper.gromacs import Sort_Atoms_By_Gro
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run

    s = "A U C G C U A C A U G"

    mol = Xponge.ResidueType.get_type("A")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("C")

    Path("rna").mkdir(exist_ok=True)
    Xponge.save_pdb(mol, "rna/rna.pdb")
    with open("rna/rna.pdb") as f:
        contents = f.read()
    contents = contents.replace("OP1", "O1P").replace("OP2", "O2P").replace(" C7", "C5M")
    with open("rna/rna.pdb", "w") as f:
        f.write(contents)
    with open("rna/ter.in", "w") as f:
        f.write("4\n6\n")
    assert os.system("cd rna && gmx pdb2gmx -f rna.pdb -ff charmm36-jul2022 -water tip3p -ignh -ter \
< ter.in > pdb2gmx.log 2>&1") == 0

    mol.residues[0].set_type("A5")
    mol.residues[12].set_type("C3")
    mol.add_missing_atoms()

    Sort_Atoms_By_Gro(mol, "rna/conf.gro")
    Xponge.save_sponge_input(mol, "rna/rna")

    assert run("SPONGE -mode minimization -default_in_file_prefix rna/rna -step_limit 2000 \
-rst rna/min > rna/min.log") == 0
    Xponge.load_coordinate("rna/min_coordinate.txt", mol)
    Xponge.save_gro(mol, "rna/min.gro")

    with open("rna/rna.mdp", "w") as f:
        f.write("""integrator = sd
dt = 2e-3
constraint_algorithm = lincs
constraints = h-bonds
nsteps = 1000
nstxout = 1
nstlog = 1
coulombtype = PME
vdw-modifier = None
ref_t = 300
tau_t = 1
nstlist = 1
tc_grps = system
DispCorr = Ener
""")
    assert os.system("cd rna && gmx grompp -f rna.mdp -c min.gro -p topol.top \
-o run.tpr -maxwarn 2 > runpp.log 2>&1") == 0
    assert os.system("cd rna && gmx mdrun -deffnm run > mdrun.log 2>&1") == 0

    assert os.system("cd rna && Xponge converter -p rna_mass.txt -c run.trr -o run.dat > convert.log") == 0
    assert run("SPONGE -mode rerun -default_in_file_prefix rna/rna -crd rna/run.dat \
-mdout rna/mdout.txt -box rna/run.box > rna/rerun.log") == 0

    gmx_mdout = _get_energies("rna/run.log")
    t = MdoutReader("rna/mdout.txt")

    Path("rna_plot").mkdir(exist_ok=True)
    _check_one_energy(gmx_mdout["Potential"], "rna_plot/Potential", t.potential)
    _check_one_energy(gmx_mdout["Bond"], "rna_plot/Bond", t.bond)
    _check_one_energy(gmx_mdout["U-B"], "rna_plot/Urey_Bradly", t.urey_bradley)
    _check_one_energy(gmx_mdout["Proper Dih."], "rna_plot/Dihedral", t.dihedral)
    _check_one_energy(gmx_mdout["Improper Dih."], "rna_plot/Improper", t.improper_dihedral)
    _check_one_energy(gmx_mdout["LJ-14"], "rna_plot/1-4 NB", t.nb14_LJ)
    _check_one_energy(gmx_mdout["Coulomb-14"], "rna_plot/1-4 EEL", t.nb14_EE)
    _check_one_energy(gmx_mdout["LJ (SR)"] + gmx_mdout["Disper. corr."], "rna_plot/LJ", t.LJ)
    _check_one_energy(gmx_mdout["Coulomb (SR)"] + gmx_mdout["Coul. recip."], "rna_plot/PME", t.PME)


def test_dna():
    """
        test the single point energy for residues of DNA
    """
    import Xponge
    import Xponge.forcefield.charmm.charmm36
    from Xponge.helper.gromacs import Sort_Atoms_By_Gro
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run

    s = "DA DT DC DG DC DT DA"

    mol = Xponge.ResidueType.get_type("DA")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("DC")

    Path("dna").mkdir(exist_ok=True)
    Xponge.save_pdb(mol, "dna/dna.pdb")
    with open("dna/dna.pdb") as f:
        contents = f.read()
    contents = contents.replace("OP1", "O1P").replace("OP2", "O2P").replace(" C7", "C5M")
    with open("dna/dna.pdb", "w") as f:
        f.write(contents)
    with open("dna/ter.in", "w") as f:
        f.write("4\n6\n")
    assert os.system("cd dna && gmx pdb2gmx -f dna.pdb -ff charmm36-jul2022 -water tip3p -ignh -ter \
< ter.in > pdb2gmx.log 2>&1") == 0

    mol.residues[0].set_type("DA5")
    mol.residues[8].set_type("DC3")
    mol.add_missing_atoms()

    Sort_Atoms_By_Gro(mol, "dna/conf.gro")
    Xponge.save_sponge_input(mol, "dna/dna")

    assert run("SPONGE -mode minimization -default_in_file_prefix dna/dna -step_limit 2000 \
-rst dna/min > dna/min.log") == 0
    Xponge.load_coordinate("dna/min_coordinate.txt", mol)
    Xponge.save_gro(mol, "dna/min.gro")

    with open("dna/dna.mdp", "w") as f:
        f.write("""integrator = sd
dt = 2e-3
constraint_algorithm = lincs
constraints = h-bonds
nsteps = 1000
nstxout = 1
nstlog = 1
coulombtype = PME
vdw-modifier = None
ref_t = 300
tau_t = 1
nstlist = 1
tc_grps = system
DispCorr = Ener
""")
    assert os.system("cd dna && gmx grompp -f dna.mdp -c min.gro -p topol.top \
-o run.tpr -maxwarn 2 > runpp.log 2>&1") == 0
    assert os.system("cd dna && gmx mdrun -deffnm run > mdrun.log 2>&1") == 0

    assert os.system("cd dna && Xponge converter -p dna_mass.txt -c run.trr -o run.dat > convert.log") == 0
    assert run("SPONGE -mode rerun -default_in_file_prefix dna/dna -crd dna/run.dat \
-mdout dna/mdout.txt -box dna/run.box > dna/rerun.log") == 0

    gmx_mdout = _get_energies("dna/run.log")
    t = MdoutReader("dna/mdout.txt")

    Path("dna_plot").mkdir(exist_ok=True)
    _check_one_energy(gmx_mdout["Potential"], "dna_plot/Potential", t.potential)
    _check_one_energy(gmx_mdout["Bond"], "dna_plot/Bond", t.bond)
    _check_one_energy(gmx_mdout["U-B"], "dna_plot/Urey_Bradly", t.urey_bradley)
    _check_one_energy(gmx_mdout["Proper Dih."], "dna_plot/Dihedral", t.dihedral)
    _check_one_energy(gmx_mdout["Improper Dih."], "dna_plot/Improper", t.improper_dihedral)
    _check_one_energy(gmx_mdout["LJ-14"], "dna_plot/1-4 NB", t.nb14_LJ)
    _check_one_energy(gmx_mdout["Coulomb-14"], "dna_plot/1-4 EEL", t.nb14_EE)
    _check_one_energy(gmx_mdout["LJ (SR)"] + gmx_mdout["Disper. corr."], "dna_plot/LJ", t.LJ)
    _check_one_energy(gmx_mdout["Coulomb (SR)"] + gmx_mdout["Coul. recip."], "dna_plot/PME", t.PME)
