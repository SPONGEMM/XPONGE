"""
    This **module** gives the unit tests of umbrella sampling
"""

__all__ = ["test_meta1d", "test_meta2d"]

def test_meta1d():
    """ test 1-dimensional metadynamics """
    import numpy as np
    import matplotlib.pyplot as plt
    import Xponge
    from Xponge.forcefield.base.angle_base import AngleType
    import Xponge.forcefield.amber.tip3p
    from Xponge.helper.cv import CVSystem
    from Xponge.mdrun import run
    from Xponge.analysis import MdoutReader
    from scipy.stats import gaussian_kde
    from pathlib import Path

    Path("meta1d").mkdir(exist_ok=True)
    min_command = "SPONGE -mode minimization -step_limit 10000 -default_in_file_prefix meta1d/test \
                   -default_out_file_prefix meta1d/min \
                   -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000 > temp.out"

    run_command = "SPONGE -mode NVT -dt 1e-3 -step_limit 1000000 -default_in_file_prefix meta1d/test \
                   -cv_in_file meta1d/cv.txt -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000 \
                   -thermostat andersen_thermostat -coordinate_in_file meta1d/min_coordinate.txt \
                   -default_out_file_prefix meta1d/out -write_information_interval 100 > temp.out"

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    AngleType.New_From_String("""name       k   b
                                 HW-HW-HW   50  1.7""")
    assign.atom_types = [hw, hw, hw, hw]
    tes = assign.to_residuetype("TES")
    mol = Xponge.save_sponge_input(tes, "meta1d/test")
    with open("meta1d/test_dihedral.txt", "w") as f:
        f.write("""2
2 0 1 3 2 15 0.4
2 0 1 3 3 5 -0.6
""")

    cv = CVSystem(mol)
    cv.add_cv_dihedral("torsion", mol.atoms[2], mol.atoms[0], mol.atoms[1], mol.atoms[3])
    cv.meta1d("torsion", CV_grid=6285, CV_minimal=-3.142, CV_maximum=3.142, welltemp_factor=50, height=1, CV_sigma=0.5)
    cv.print("torsion")
    cv.output("meta1d/cv.txt")

    assert run(min_command) == 0
    assert run(run_command) == 0
    t = MdoutReader("meta1d/out.out")
    bias = t.metad
    t = t.torsion
    kt = -8.314 * 300 / 4184
    w = np.exp(-bias/kt)
    t = np.concatenate((t, t + np.pi * 2, t - np.pi * 2))
    w = np.concatenate((w, w, w))
    kernel = gaussian_kde(t, weights=w, bw_method=0.01)
    positions = np.linspace(-np.pi, np.pi, 300)
    result = kernel(positions)
    result = -8.314 * 300 / 4184 * np.log(result)
    result -= min(result)
    theory = 0.3 * np.cos(2 * positions - 0.4) + 0.1 * np.cos(3 * positions + 0.6)
    theory *= 50
    theory -= min(theory)
    toread = np.loadtxt("meta1d/out_metad_potential.txt", skiprows=2)
    toread = -toread
    toread -= np.min(toread)
    plt.plot(positions, result, label="simulated results")
    plt.plot(positions, theory, label="potential")
    plt.plot(np.linspace(-3.142, 3.142, 6285), toread, label="meta1d potential")
    plt.legend()
    plt.savefig("meta1d.png")
    plt.clf()

def test_meta2d():
    """ test 2-dimensional metadynamics """
    
