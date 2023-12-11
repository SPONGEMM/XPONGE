"""
    This **module** includes unittests of the Xponge.forcefield.opls.oplsaam
"""
import os
from pathlib import Path

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
            while word and start <= 60:
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
    import Xponge.forcefield.opls.oplsaam
    from Xponge.helper.gromacs import Sort_Atoms_By_Gro
    from Xponge.analysis import MdoutReader
    from Xponge.mdrun import run

    s = "ALA ARG ASN ASP CYS GLN GLU GLY ILE LEU " + \
        "LYS MET PHE PRO SER THR TRP TYR VAL"

    mol = Xponge.ResidueType.get_type("ALA")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("ALA")

    Path("protein").mkdir(exist_ok=True)
    Xponge.save_pdb(mol, "protein/protein.pdb")
    assert os.system("cd protein && gmx pdb2gmx -f protein.pdb -ff oplsaam -water tip4p -ignh \
> pdb2gmx.log 2>&1") == 0

    mol.residues[0].set_type("NALA")
    mol.residues[-1].set_type("CALA")
    mol.add_missing_atoms()
    mol.box_length = [25, 25, 80]
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
    _check_one_energy(gmx_mdout["Angle"], "protein_plot/Angle", t.angle)
    _check_one_energy(gmx_mdout["Ryckaert-Bell."], "protein_plot/Ryckaert_Bell", t.Ryckaert_Bellemans)
    _check_one_energy(gmx_mdout["Coulomb-14"], "protein_plot/1-4 EEL", t.nb14_EE)
    _check_one_energy(gmx_mdout["LJ (SR)"] + gmx_mdout["Disper. corr."], "protein_plot/LJ", t.LJ)
    _check_one_energy(gmx_mdout["Coulomb (SR)"] + gmx_mdout["Coul. recip."], "protein_plot/PME", t.PME)
    _check_one_energy(gmx_mdout["Periodic Improp"], "protein_plot/Improper", t.dihedral)
    _check_one_energy(gmx_mdout["er Dih."], "protein_plot/1-4 NB", t.nb14_LJ)
