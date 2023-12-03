"""
    This **module** gives the unit tests of umbrella sampling
"""

__all__ = ["test_umbrella"]

def test_umbrella():
    """ test umbrella sampling """
    import shutil
    import numpy as np
    import matplotlib.pyplot as plt
    import Xponge
    from Xponge.forcefield.base.angle_base import AngleType
    import Xponge.forcefield.amber.tip3p
    from Xponge.helper.cv import CVSystem
    from Xponge.mdrun import run
    from Xponge.analysis import wham

    min_command = "SPONGE -mode minimization -step_limit 10000 -default_in_file_prefix test \
                   -cv_in_file cv.txt -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000"

    run_command = "SPONGE -mode NVT -dt 1e-3 -step_limit 20000 -default_in_file_prefix test \
                   -cv_in_file cv.txt -cutoff 1 -skin 1 -neighbor_list_refresh_interval 100000 \
                   -thermostat andersen_thermostat -coordinate_in_file restart_coordinate.txt \
                   -write_information_interval 100"
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
    cv.restrain("torsion", weight=200, reference="res_cv_ref")
    cv.print("torsion")
    cv.output("cv.txt")

    for i, ref in enumerate(np.linspace(-np.pi,np.pi,50)):
        assert run(min_command + f" -res_cv_ref {ref} > temp.out") == 0
        assert run(run_command + f" -res_cv_ref {ref} > temp.out") == 0
        shutil.copy("mdout.txt", f"{i}.mdout")

    w = wham.WHAM(np.linspace(-np.pi, np.pi, 51), 300, 200, np.linspace(-np.pi, np.pi, 50), 2 * np.pi)
    w.get_data_from_mdout("*.mdout", "torsion")
    x, y, f = w.main()
    plt.plot(f)
    plt.savefig("free_energy_of_simulations.png")
    plt.clf()
    plt.plot(x, y, label="umbrella sampling")
    y2 = 15 * np.cos(2 * x - 0.4) + 5 * np.cos(3 * x + 0.6)
    y2 -= min(y2)
    plt.plot(x, y2, label="potential")
    plt.legend()
    plt.savefig("sampling.png")
    plt.clf()
