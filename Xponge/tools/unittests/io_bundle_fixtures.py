"""Shared fixtures for SPONGE bundled I/O integration tests."""

from pathlib import Path

import numpy as np


def write_basic_case(case_dir: Path) -> None:
    (case_dir / "coordinate.txt").write_text(
        "2\n"
        "1.0 2.0 3.0\n"
        "4.0 5.0 6.0\n"
        "10.0 20.0 30.0 90.0 90.0 90.0\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(
        "2\n"
        "0.1 0.2 0.3\n"
        "0.4 0.5 0.6\n",
        encoding="utf-8",
    )
    cv_config = (
        "print\n"
        "{\n"
        "    CV = distance\n"
        "}\n"
        "distance\n"
        "{\n"
        "    CV_type = distance\n"
        "    atom = 0 1\n"
        "}\n"
    )
    (case_dir / "cv.txt").write_text(cv_config, encoding="utf-8")
    (case_dir / "restrain.txt").write_text(cv_config, encoding="utf-8")
    (case_dir / "restrain_cv.txt").write_text(
        "restrain\n"
        "{\n"
        "    CV = distance\n"
        "    weight = 10.0\n"
        "    reference = 1.5\n"
        "}\n",
        encoding="utf-8",
    )
    (case_dir / "steer_cv.txt").write_text(
        "steer\n"
        "{\n"
        "    CV = distance\n"
        "    weight = 2.0\n"
        "}\n",
        encoding="utf-8",
    )
    (case_dir / "sits.txt").write_text(
        "SITS\n"
        "{\n"
        "    CV = distance\n"
        "    nk = 2\n"
        "}\n",
        encoding="utf-8",
    )
    (case_dir / "meta_edge.txt").write_text(
        "0.0 0.0 1.0 0.1 0.2\n"
        "1.0 0.0 2.0 0.3 0.4\n",
        encoding="utf-8",
    )
    (case_dir / "meta_potential.txt").write_text(
        "Meta potential\n"
        "0.0 1.0 0.5\n"
        "0.0 2.0 1.0\n"
        "2 2 4\n"
        "0.0 0.0 0.0 0.1 0.2 3.0\n"
        "0.5 0.0 0.0 0.3 0.4 4.0\n"
        "0.0 1.0 0.0 0.5 0.6 5.0\n"
        "0.5 1.0 0.0 0.7 0.8 6.0\n",
        encoding="utf-8",
    )
    (case_dir / "meta_scatter.txt").write_text(
        "Meta scatter\n"
        "0.0 1.0 0.5\n"
        "2 2\n"
        "0.0 0.0 1.0\n"
        "0.5 0.0 2.0\n",
        encoding="utf-8",
    )
    (case_dir / "hills.txt").write_text(
        "0.1 0.2 1.5\n"
        "0.3 0.4 1.6\n",
        encoding="utf-8",
    )
    (case_dir / "restrain.rst7").write_text(
        "restrain reference\n"
        "2 0.5\n"
        "1.0 2.0 3.0 4.0 5.0 6.0\n"
        "0.1 0.2 0.3 0.4 0.5 0.6\n"
        "10.0 20.0 30.0 90.0 90.0 90.0\n",
        encoding="utf-8",
    )
    (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("2\n0.1\n-0.1\n", encoding="utf-8")
    (case_dir / "residue.txt").write_text("2 1\n2\n", encoding="utf-8")
    (case_dir / "exclude.txt").write_text("2 1\n1 1\n0 \n", encoding="utf-8")
    (case_dir / "bond.txt").write_text("1\n0 1 100.0 1.5\n", encoding="utf-8")
    (case_dir / "angle.txt").write_text("1\n0 1 0 50.0 1.9\n", encoding="utf-8")
    (case_dir / "dihedral.txt").write_text("1\n0 1 0 1 3 2.0 3.14\n", encoding="utf-8")
    (case_dir / "improper.txt").write_text("1\n0 1 0 1 2.5 3.14\n", encoding="utf-8")
    (case_dir / "nb14_extra.txt").write_text("1\n0 1 1.0 2.0 0.5\n", encoding="utf-8")
    (case_dir / "urey_bradley.txt").write_text("1\n0 1 0 50.0 1.9 10.0 2.5\n", encoding="utf-8")
    (case_dir / "cmap.txt").write_text(
        "1 1\n"
        "2\n"
        "1.0 2.0\n"
        "3.0 4.0\n"
        "0 1 0 1 0 0\n",
        encoding="utf-8",
    )
    (case_dir / "gb.txt").write_text("2\n1.5 0.8\n1.2 0.85\n", encoding="utf-8")
    (case_dir / "virtual_atom.txt").write_text("2 1 0 1 0 0.5 0.5\n", encoding="utf-8")
    (case_dir / "lj.txt").write_text(
        "2 1\n\n"
        "1.0000000e+00\n\n"
        "2.0000000e+00\n\n"
        "0\n"
        "0\n",
        encoding="utf-8",
    )
    (case_dir / "lj_soft_core.txt").write_text(
        "2 1 1\n"
        "1.0\n"
        "2.0\n"
        "3.0\n"
        "4.0\n"
        "0 0\n"
        "0 0\n",
        encoding="utf-8",
    )
    (case_dir / "subsys_division.txt").write_text("2\n1\n2\n", encoding="utf-8")
    (case_dir / "eam.txt").write_text(
        "Synthetic Cu funcfl\n"
        "29 63.55 3.615 FCC\n"
        "2 0.5 2 0.1 5.0\n"
        "1.0 2.0\n"
        "3.0 4.0\n"
        "5.0 6.0\n",
        encoding="utf-8",
    )
    (case_dir / "eam_atom_type.txt").write_text("0\n0\n", encoding="utf-8")
    (case_dir / "sw.txt").write_text(
        "2 1\n"
        "# pair\n"
        "0 0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0\n"
        "# triple\n"
        "0 0 0 9.0 10.0 11.0\n"
        "# atom types\n"
        "0\n"
        "0\n",
        encoding="utf-8",
    )
    (case_dir / "edip.txt").write_text(
        "2 1\n"
        "# pair\n"
        "0 0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0\n"
        "# triple\n"
        "0 0 0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0\n"
        "# atom types\n"
        "0\n"
        "0\n",
        encoding="utf-8",
    )
    (case_dir / "tersoff.txt").write_text(
        "2 1\n"
        "Si\n"
        "Si Si Si 3.0 1.0 0.0 25000.0 4.3484 -0.89 0.72751 0.000000125724 2.199 340.0 1.95 0.05 3.568 1380.0\n"
        "# Atom types\n"
        "0 0\n",
        encoding="utf-8",
    )
    (case_dir / "qc_type.txt").write_text(
        "2 0 1\n"
        "0 H\n"
        "1 O\n",
        encoding="utf-8",
    )
    (case_dir / "reaxff_type.txt").write_text(
        "2\n"
        "O\n"
        "H\n",
        encoding="utf-8",
    )
    (case_dir / "reaxff.txt").write_text(
        "Synthetic ReaxFF parameter file\n"
        "2 ! Number of general parameters\n"
        "50.0 !p(boc1)\n"
        "9.5 !p(boc2)\n"
        "2 ! Nr of atoms\n"
        "atom header 1\n"
        "atom header 2\n"
        "atom header 3\n"
        "O 1.2 2.0 15.999 1.9 0.09 1.05 1.08 6.0 ! O line 1\n"
        "10.2 7.7 4.0 36.9 116.0 8.5 8.9 2.0\n"
        "0.9 1.0 60.8 20.4 3.3 0.2 0.9 0.0\n"
        "-3.6 2.7 1.0 4.0 2.9 0.0 0.0 0.0\n"
        "H 0.7 1.0 1.008 1.5 0.04 1.02 -0.1 1.0 ! H line 1\n"
        "9.3 5.0 1.0 0.0 121.1 5.3 7.4 1.0\n"
        "-0.1 0.0 62.4 1.9 3.3 0.7 1.0 0.0\n"
        "-15.7 2.1 1.0 1.0 2.8 0.0 0.0 0.0\n"
        "1 ! Nr of bonds\n"
        "bond header\n"
        "1 2 170.0 0.0 0.0 -0.5 0.0 1.0 6.0 0.7\n"
        "5.2 1.0 0.0 1.0 -0.05 6.8 0.0 0.0\n"
        "1 ! Nr of off-diagonal terms\n"
        "1 2 0.12 1.4 9.8 1.1 -1.0 -1.0\n"
        "1 ! Nr of angles\n"
        "1 2 1 70.0 25.0 3.4 0.0 0.0 0.0 3.0\n"
        "1 ! Nr of torsions\n"
        "0 1 2 0 0.0 0.1 0.02 -2.5 0.0 0.0 0.0\n"
        "1 ! Nr of hydrogen bonds\n"
        "1 2 1 1.9 -4.4 1.7 3.0\n",
        encoding="utf-8",
    )
    (case_dir / "pairwise_force.txt").write_text(
        "[[[ custom_pair ]]]\n"
        "[[ potential ]]\n"
        "E = epsilon_ij * powf(sigma_ij / r_ij, 12.0f);\n"
        "[[ parameters ]]\n"
        "float epsilon_ij, float sigma_ij\n"
        "[[ with_ele ]]\n"
        "false\n"
        "[[ end ]]\n",
        encoding="utf-8",
    )
    (case_dir / "custom_pair.txt").write_text(
        "2 1\n"
        "1.0\n"
        "2.0\n"
        "0\n"
        "0\n",
        encoding="utf-8",
    )
    (case_dir / "listed_forces.txt").write_text(
        "[[[ custom_bond ]]]\n"
        "[[ potential ]]\n"
        "E = k * (r_ab - r0) * (r_ab - r0);\n"
        "[[ parameters ]]\n"
        "int atom_a, int atom_b, float k, float r0\n"
        "[[ connected_atoms ]]\n"
        "ab\n"
        "[[ constrain_distance ]]\n"
        "r0\n"
        "[[ end ]]\n",
        encoding="utf-8",
    )
    (case_dir / "custom_bond.txt").write_text(
        "1\n"
        "0 1 10.0 1.5\n",
        encoding="utf-8",
    )
    (case_dir / "sits_nk.txt").write_text("1.0 2.0\n", encoding="utf-8")
    (case_dir / "restrain_atom_id.txt").write_text("0\n1\n", encoding="utf-8")
    (case_dir / "restrain_weight.txt").write_text("10.0 0.0 0.0\n0.0 20.0 0.0\n", encoding="utf-8")
    (case_dir / "restrain_coordinate.txt").write_text(
        "2\n"
        "1.0 2.0 3.0\n"
        "4.0 5.0 6.0\n",
        encoding="utf-8",
    )
    (case_dir / "nhc_restart.txt").write_text(
        "0.1 0.2\n"
        "0.3 0.4\n",
        encoding="utf-8",
    )
    (case_dir / "constrain.txt").write_text("1\n0 1 1.5\n", encoding="utf-8")
    (case_dir / "sits_atom.txt").write_text("0\n1\n", encoding="utf-8")
    (case_dir / "soft_walls.txt").write_text(
        "[[[ z_wall ]]]\n"
        "[[ potential ]]\n"
        "E = powf((z - 5.0f) * (z - 5.0f), -6.0f);\n"
        "[[ end ]]\n",
        encoding="utf-8",
    )
    np.asarray(
        [
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
        ],
        dtype=np.float32,
    ).tofile(case_dir / "traj.dat")
    np.asarray(
        [
            [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            [[0.2, 0.3, 0.4], [0.5, 0.6, 0.7]],
        ],
        dtype=np.float32,
    ).tofile(case_dir / "traj_vel.dat")
    (case_dir / "traj_box.dat").write_text(
        "10.0 20.0 30.0 90.0 90.0 90.0\n"
        "11.0 21.0 31.0 90.0 90.0 90.0\n",
        encoding="utf-8",
    )
    (case_dir / "mdin.spg.toml").write_text(
        'mode = "rerun"\n'
        "step_limit = 10\n"
        'coordinate_in_file = "coordinate.txt"\n'
        'velocity_in_file = "velocity.txt"\n'
        'mass_in_file = "mass.txt"\n'
        'charge_in_file = "charge.txt"\n'
        'residue_in_file = "residue.txt"\n'
        'exclude_in_file = "exclude.txt"\n'
        'bond_in_file = "bond.txt"\n'
        'angle_in_file = "angle.txt"\n'
        'dihedral_in_file = "dihedral.txt"\n'
        'improper_dihedral_in_file = "improper.txt"\n'
        'nb14_extra_in_file = "nb14_extra.txt"\n'
        'urey_bradley_in_file = "urey_bradley.txt"\n'
        'cmap_in_file = "cmap.txt"\n'
        'gb_in_file = "gb.txt"\n'
        'virtual_atom_in_file = "virtual_atom.txt"\n'
        'LJ_in_file = "lj.txt"\n'
        'LJ_soft_core_in_file = "lj_soft_core.txt"\n'
        'subsys_division_in_file = "subsys_division.txt"\n'
        'SW_in_file = "sw.txt"\n'
        'EDIP_in_file = "edip.txt"\n'
        'pairwise_force_in_file = "pairwise_force.txt"\n'
        'custom_pair_in_file = "custom_pair.txt"\n'
        'listed_forces_in_file = "listed_forces.txt"\n'
        'custom_bond_in_file = "custom_bond.txt"\n'
        'qc_type_in_file = "qc_type.txt"\n'
        'cv_in_file = "cv.txt"\n'
        'restrain_in_file = "restrain.txt"\n'
        'restrain_cv_in_file = "restrain_cv.txt"\n'
        'steer_cv_in_file = "steer_cv.txt"\n'
        'SITS_in_file = "sits.txt"\n'
        'SITS_nk_in_file = "sits_nk.txt"\n'
        'meta_edge_in_file = "meta_edge.txt"\n'
        'meta_potential_in_file = "meta_potential.txt"\n'
        'meta_scatter_in_file = "meta_scatter.txt"\n'
        'hills_in_file = "hills.txt"\n'
        'restrain_atom_id = "restrain_atom_id.txt"\n'
        'restrain_weight_in_file = "restrain_weight.txt"\n'
        'restrain_coordinate_in_file = "restrain_coordinate.txt"\n'
        'nose_hoover_chain_restart_input = "nhc_restart.txt"\n'
        'constrain_in_file = "constrain.txt"\n'
        'SITS_atom_in_file = "sits_atom.txt"\n'
        'soft_walls_in_file = "soft_walls.txt"\n'
        'crd = "traj.dat"\n'
        'box = "traj_box.dat"\n'
        'vel = "traj_vel.dat"\n'
        'output_h5_trajectory_path = "prod.spg.h5md"\n'
        "output_h5_trajectory_vds = true\n"
        'mdout = "mdout.txt"\n'
        '[EAM]\n'
        'in_file = "eam.txt"\n'
        'atom_type_in_file = "eam_atom_type.txt"\n'
        '[TERSOFF]\n'
        'in_file = "tersoff.txt"\n'
        '[REAXFF]\n'
        'in_file = "reaxff.txt"\n'
        'type_in_file = "reaxff_type.txt"\n',
        encoding="utf-8",
    )
