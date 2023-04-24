"""
This **module** gives the helper functions to do mol2rfe
"""
import sys
import os
import shutil

from .. import Molecule, load, build
from ..analysis import MdoutReader
from ..helper import Xopen, Xprint
from ..mdrun import run


__all__ = ["_mol2rfe_build", "_mol2rfe_min", "_mol2rfe_pre_equilibrium",
           "_mol2rfe_equilibrium", "_mol2rfe_analysis"]

# pylint: disable=unused-argument
def _do_nothing(mol_r, mol_a, mol_b, forcetype, rforces, bforces, lambda_, mol_r2mol_a, mol_a2mol_r,
                mol_r2mol_b, mol_b2mol_r):
    """do nothing to save time"""

def _mol2rfe_build(args, merged_from, merged_to):
    """

    :param args:
    :param merged_from:
    :param merged_to:
    :return:
    """
    from ..forcefield.special import fep

    if "build" in args.do:
        Xprint("\nBUILDING TOPOLOGY\n", "INFO")

        fep.Save_Soft_Core_LJ()
        unchanged = ["bond", "angle", "dihedral", "coordinate", "exclude",
                     "LJ_soft_core", "mass", "residue", "subsys_division"]
        for i in range(args.nl + 1):
            if os.path.exists("%d" % i):
                shutil.rmtree("%d" % i)
            os.mkdir("%d" % i)
            if i == 1:
                for j in unchanged:
                    getattr(Molecule, "_save_functions").pop(j)
                    fep.FEP_BONDED_FORCE_MERGE_RULE[j] = {"lambda_name": "_", "merge_function": _do_nothing}
            if i > 0:
                for j in unchanged:
                    shutil.copy(f"0/{args.temp}_{j}.txt", f"{i}/{args.temp}_{j}.txt")
            tt = fep.Merge_Force_Field(merged_from, merged_to, args.l[i], {"charge": args.l[i] ** args.cp})
            build.save_mol2(tt, "%d/%s.mol2" % (i, args.temp))
            build.Save_SPONGE_Input(tt, "%d/%s" % (i, args.temp))
            Xprint(f"{i} built success")

def _mol2rfe_output_path(subdir, workdir, tempname):
    """

    :param subdir:
    :param workdir:
    :param tempname:
    :return:
    """
    toadd = " -mdinfo {2}/{0}/{1}.mdinfo -mdout {2}/{0}/{1}.mdout".format(subdir, tempname, workdir)
    toadd += " -rst {2}/{0}/{1} -crd {2}/{0}/{1}.dat -box {2}/{0}/{1}.box".format(subdir, tempname, workdir)
    return toadd


def _mol2rfe_min(args, iteror):
    """

    :param args:
    :return:
    """
    if "min" in args.do:
        for i in iteror:
            if os.path.exists("%d/min" % i):
                shutil.rmtree("%d/min" % i)
            os.mkdir("%d/min" % i)
            basic = f"SPONGE -default_in_file_prefix {i}/{args.temp}"
            lambda_ = args.l[i]
            basic += f" -mode minimization -lambda_lj {lambda_} -cutoff 8 -write_information_interval 1000"
            basic += f" -neighbor_list_max_atom_in_grid_numbers 128 -neighbor_list_max_neighbor_numbers 1200"
            basic += _mol2rfe_output_path("min", i, args.temp)
            if i != 0:
                basic += f" -coordinate_in_file 0/pre_equilibrium/{args.temp}_coordinate.txt"
            dt_factor = 1e-2
            inc_rate = 1.5
            if not args.mi:
                cif = " -minimization_dynamic_dt 1"
                exit_code = run(f"{basic} {cif} -step_limit {args.min_step} \
-minimization_dt_factor {dt_factor} -minimization_dt_increasing_rate {inc_rate}")
                out = MdoutReader(f"{i}/min/{args.temp}.mdout").potential[-1]
                min_time = 0
                while (out > 0 or exit_code != 0) and min_time < 10:
                    if exit_code != 0:
                        dt_factor /= 3
                        inc_rate -= 0.049
                    else:
                        dt_factor *= 2
                        inc_rate += 0.04
                    Xprint("Minimization will be repeated to reduce the potential to 0", "WARNING")
                    min_time += 1
                    exit_code = run(f"{basic} {cif} -step_limit {args.min_step} \
-minimization_dt_factor {dt_factor} -minimization_dt_increasing_rate {inc_rate}")
                    out = MdoutReader(f"{i}/min/{args.temp}.mdout").potential[-1]
                if min_time >= 10:
                    Xprint("Minimization has been repeated for 10 times and \
the potential still can not be reduced to 0", "ERROR")
                    sys.exit(1)
            else:
                command += f" -mdin {args.mi}"
                exit_code = run(command)

def _mol2rfe_pre_equilibrium(args, iteror):
    """

    :param args:
    :return:
    """
    if "pre_equilibrium" in args.do:
        for i in iteror:
            if os.path.exists("%d/pre_equilibrium" % i):
                shutil.rmtree("%d/pre_equilibrium" % i)
            os.mkdir("%d/pre_equilibrium" % i)
            command = f"SPONGE -default_in_file_prefix {i}/{args.temp}"
            lambda_ = args.l[i]
            command += f" -lambda_lj {lambda_} -cutoff 8"
            command += f" -neighbor_list_max_atom_in_grid_numbers 128 -neighbor_list_max_neighbor_numbers 1200"
            command += _mol2rfe_output_path("pre_equilibrium", i, args.temp)
            command += f" -coordinate_in_file {i}/min/{args.temp}_coordinate.txt"
            if not args.pi:
                command += f" -mode NPT -step_limit {args.pre_equilibrium_step}"
                command += f" -dt {args.dt} -constrain_mode SHAKE"
                command += " -barostat andersen_barostat -thermostat middle_langevin"
                command += " -middle_langevin_gamma 10 -middle_langevin_velocity_max 20"
                exit_code = run(command)
            else:
                command += f" -mdin {args.pi}"
                exit_code = run(command)
            if exit_code != 0:
                Xprint(f"The pre_equilibrium of lambda {i} exited with code {exit_code}", "ERROR")
                sys.exit(exit_code)


def _mol2rfe_equilibrium(args):
    """

    :param args:
    :return:
    """
    if "equilibrium" in args.do:
        for i in range(args.nl + 1):
            if os.path.exists("%d/equilibrium" % i):
                os.system("rm -rf %d/equilibrium" % i)
            os.mkdir("%d/equilibrium" % i)
            command = f"SPONGE -default_in_file_prefix {i}/{args.temp}"
            lambda_ = args.l[i]
            command += f" -lambda_lj {lambda_} -cutoff 8 -molecule_map_output 0"
            command += f" -neighbor_list_max_atom_in_grid_numbers 128 -neighbor_list_max_neighbor_numbers 1200"
            command += _mol2rfe_output_path("equilibrium", i, args.temp)
            command += f" -coordinate_in_file {i}/pre_equilibrium/{args.temp}_coordinate.txt"
            command += f" -velocity_in_file {i}/pre_equilibrium/{args.temp}_velocity.txt"
            if not args.ei:
                command += f" -mode NPT -step_limit {args.equilibrium_step} -dt {args.dt} -constrain_mode SHAKE"
                command += " -barostat andersen_barostat -thermostat middle_langevin"
                command += " -middle_langevin_gamma 10 -middle_langevin_velocity_max 20"
                command += f" -write_information_interval 100 -write_restart_file_interval {args.equilibrium_step}"
                run(command)
            else:
                command += f" -mdin {args.pi}"
                run(command)


def _mol2rfe_analysis(args, merged_from, merged_to):
    """

    :param args:
    :param merged_from:
    :param merged_to:
    :return:
    """
    if "analysis" in args.do:
        f = Xopen("dh_dlambda.txt", "w")
        f.close()
        resname = merged_from.residues[args.ri].name
        draw_r1_mol = merged_from.deepcopy()
        draw_r2_mol = merged_to.deepcopy()
        load.load_coordinate(f"0/equilibrium/{args.temp}_coordinate.txt", draw_r1_mol)
        load.load_coordinate(f"{args.nl}/equilibrium/{args.temp}_coordinate.txt", draw_r2_mol)
        draw_r1_res = draw_r1_mol.residues[args.ri]
        draw_r2_res = draw_r2_mol.residues[args.ri]
        draw_r1_res.name = resname.split("_")[0]
        draw_r2_res.name = resname.split("_")[1]
        to_delete = []
        for atom in draw_r1_res.atoms:
            if atom.LJtype == "ZERO_LJ_ATOM":
                to_delete.append(atom)
        for atom in to_delete:
            draw_r1_res.atoms.remove(atom)
        to_delete = []
        for atom in draw_r2_res.atoms:
            if atom.LJtype == "ZERO_LJ_ATOM":
                to_delete.append(atom)
            if atom.name.endswith("R2"):
                atom.name = atom.name[:-2]
        for atom in to_delete:
            draw_r2_res.atoms.remove(atom)
        build.save_pdb(draw_r1_mol, "r1.pdb")
        build.save_pdb(draw_r2_mol, "r2.pdb")
        if args.method == "TI":
            from .ti import ti_analysis
            ti_analysis(args, merged_from)
        elif args.method == "FEP_BAR":
            raise NotImplementedError
