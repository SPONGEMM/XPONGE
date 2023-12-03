"""
    This **module** gives the unit tests of selective integrated enhancing sampling
"""

__all__ = ["test_sits"]

def test_sits():
    """ test SITS """
    import numpy as np
    import matplotlib.pyplot as plt
    import Xponge
    from Xponge.forcefield.base.angle_base import AngleType
    import Xponge.forcefield.amber.tip3p
    from Xponge.helper.cv import CVSystem
    from Xponge.mdrun import run
    from Xponge.analysis import MdoutReader
    from scipy.stats import gaussian_kde

    min_command = "SPONGE -mode minimization -step_limit 10000 -default_in_file_prefix test \
                   -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000 -rst min > temp.out"

    run_command = "SPONGE -mode NVT -dt 1e-3 -default_in_file_prefix test \
                   -sits_dihedral_in_file test_dihedral.txt \
                   -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000 \
                   -thermostat andersen_thermostat -coordinate_in_file min_coordinate.txt \
                   -cv_in_file cv.txt -write_information_interval 100"

    assign = Xponge.get_assignment_from_smiles("OO")
    hw = Xponge.AtomType.get_type("HW")
    AngleType.New_From_String("""name       k   b
                                 HW-HW-HW   50  1.7""")
    assign.atom_types = [hw, hw, hw, hw]
    tes = assign.to_residuetype("TES")
    mol = Xponge.save_sponge_input(tes, "test")
    with open("test_dihedral.txt", "w") as f:
        f.write("""2
2 0 1 3 2 15 0.4
2 0 1 3 3 5 -0.6
""")

    cv = CVSystem(mol)
    cv.add_cv_dihedral("torsion", mol.atoms[2], mol.atoms[0], mol.atoms[1], mol.atoms[3])
    cv.print("torsion")
    cv.output("cv.txt")

    assert run(min_command) == 0
    assert run(run_command + " -SITS_mode observation -SITS_atom_numbers 4 -step_limit 20000 > temp.out") == 0
    assert run(run_command + " -SITS_mode iteration  -SITS_atom_numbers 4 \
-SITS_T_low 100 -SITS_T_high 6000 -SITS_k_numbers 6000 -SITS_pe_b -0.0 -step_limit 100000 > temp.out") == 0
    assert run(run_command + " -SITS_mode production  -SITS_atom_numbers 4 \
-SITS_T_low 100 -SITS_T_high 6000 -SITS_k_numbers 6000 -SITS_pe_b -0.0 \
-SITS_nk_in_file SITS_nk_rest.txt -step_limit 980000 > temp.out") == 0

    t = MdoutReader("mdout.txt")
    bias = t.SITS_bias
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
    plt.plot(positions, result, label="simulated results")
    plt.plot(positions, theory, label="potential")
    plt.legend()
    plt.savefig("sits.png")
    plt.clf()
