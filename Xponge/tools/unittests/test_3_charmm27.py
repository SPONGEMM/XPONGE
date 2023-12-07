"""
    This **module** includes unittests of the Xponge.forcefield.charmm27
"""
import os

__all__ = ["test_protein"]

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
    import Xponge.forcefield.charmm.charmm27
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
    Xponge.add_solvent_box(mol, Xponge.ResidueType.get_type("WAT"), 10, tolerance=3)
    Xponge.save_pdb(mol, "protein.pdb")
    assert os.system("gmx pdb2gmx -f protein.pdb -ff charmm27 -water tips3p -ignh > pdb2gmx.log 2> pdb2gmx.log") == 0

    Sort_Atoms_By_Gro(mol, "conf.gro")
    Xponge.save_sponge_input(mol, "protein")

    assert run("SPONGE -mode minimization -default_in_file_prefix protein -step_limit 2000 -rst min > min.log") == 0
    Xponge.load_coordinate("min_coordinate.txt", mol)
    Xponge.save_gro(mol, "min.gro")

    with open("protein.mdp", "w") as f:
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
    assert os.system("gmx grompp -f protein.mdp -c min.gro -p topol.top \
-o run.tpr -maxwarn 1 > runpp.log 2> runpp.log") == 0
    assert os.system("gmx mdrun -deffnm run > mdrun.log 2> mdrun.log") == 0

    assert os.system("Xponge converter -p protein_mass.txt -c run.trr -o run.dat > convert.log") == 0
    assert run("SPONGE -mode rerun -default_in_file_prefix protein -crd run.dat -box run.box > rerun.log") == 0

    gmx_mdout = _get_energies("run.log")
    t = MdoutReader("mdout.txt")

    _check_one_energy(gmx_mdout["Potential"], "Potential", t.potential)
    _check_one_energy(gmx_mdout["Bond"], "Bond", t.bond)
    _check_one_energy(gmx_mdout["U-B"], "Urey_Bradly", t.urey_bradley)
    _check_one_energy(gmx_mdout["Proper Dih."], "Dihedral", t.dihedral)
    _check_one_energy(gmx_mdout["Improper Dih."], "Improper", t.improper_dihedral)
    _check_one_energy(gmx_mdout["CMAP Dih."], "Cmap", t.cmap)
    _check_one_energy(gmx_mdout["LJ-14"], "1-4 NB", t.nb14_LJ)
    _check_one_energy(gmx_mdout["Coulomb-14"], "1-4 EEL", t.nb14_EE)
    _check_one_energy(gmx_mdout["LJ (SR)"] + gmx_mdout["Disper. corr."], "LJ", t.LJ)
    _check_one_energy(gmx_mdout["Coulomb (SR)"] + gmx_mdout["Coul. recip."], "PME", t.PME)
