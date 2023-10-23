"""
This **module** gives the helper functions to do MM/GBSA
"""
import sys
import os
import shutil
import MDAnalysis as mda
import numpy as np

from .. import load, build, source, import_python_script
from ..analysis import MdoutReader
from ..analysis.md_analysis import SpongeTrajectoryReader
from ..helper import Xprint
from ..mdrun import run

__all__ = ["_mmgbsa_build", "_mmgbsa_min", "_mmgbsa_pre_equilibrium", "_mmgbsa_equilibrium",
"_mmgbsa_analysis"]

def _restrain_process(args):
    u = mda.Universe(f"run/{args.temp}_mass.txt", f"run/{args.temp}_coordinate.txt")
    urr = u.select_atoms(args.sr)
    shutil.copy(f"run/{args.temp}_coordinate.txt", f"run/{args.temp}_restrain_coordinate.txt")
    id2index_map = {atom.id: str(i) for i, atom in enumerate(u.atoms)}
    with open(f"run/{args.temp}_restrain_atom_id.txt", "w") as f:
        f.write("\n".join([id2index_map[atom.id] for atom in urr.atoms]))

def _mmgbsa_build(args):
    """ build the systems"""
    if "build" in args.do:
        if not args.ff:
            gaff = source("...forcefield.amber.gaff", False)
            source("...forcefield.amber.ff14sb", False)
            source("...forcefield.amber.tip3p", False)
            for extrai, mol2file in enumerate(args.r0):
                load.load_mol2(mol2file, as_template=True)
                gaff.parmchk2_gaff(mol2file, f"{args.temp}_TMP{extrai}.frcmod")
        else:
            import_python_script(args.ff)
        pdb = load.load_pdb(args.pdb)
        if os.path.exists("run"):
            shutil.rmtree("run")
        os.mkdir("run")
        if args.ivacuum:
            pdb.vacuum_layer = args.ivacuum
        Xprint("building system run")
        build.save_mol2(pdb, "run/%s.mol2" % (args.temp))
        build.Save_SPONGE_Input(pdb, "run/%s" % (args.temp))
        u = mda.Universe(f"run/{args.temp}.mol2")
        if args.sr:
            _restrain_process(args)
        ur1 = u.select_atoms(args.s1)
        ur2 = u.select_atoms(args.s2)
        ur1.write(f"{args.temp}_r1.mol2")
        ur2.write(f"{args.temp}_r2.mol2")
        r1 = load.load_mol2(f"{args.temp}_r1.mol2")
        r2 = load.load_mol2(f"{args.temp}_r2.mol2")
        complex_ = r1 | r2
        gb = source("...forcefield.special.gb", False)
        gb.set_gb_radius(r1)
        gb.set_gb_radius(r2)
        gb.set_gb_radius(complex_)
        if os.path.exists("complex"):
            shutil.rmtree("complex")
        os.mkdir("complex")
        build.save_mol2(complex_, "complex/%s.mol2" % (args.temp))
        Xprint("building system complex")
        build.Save_SPONGE_Input(complex_, "complex/%s" % (args.temp))
        if os.path.exists("part1"):
            shutil.rmtree("part1")
        os.mkdir("part1")
        build.save_mol2(r1, "part1/%s.mol2" % (args.temp))
        Xprint("building system part1")
        build.Save_SPONGE_Input(r1, "part1/%s" % (args.temp))
        if os.path.exists("part2"):
            shutil.rmtree("part2")
        os.mkdir("part2")
        build.save_mol2(r2, "part2/%s.mol2" % (args.temp))
        Xprint("building system part2")
        build.Save_SPONGE_Input(r2, "part2/%s" % (args.temp))

def _mmgbsa_output_path(subdir, tempname, add_crd=True):
    """give the output file commands"""
    toadd = " -mdinfo {0}/{1}.mdinfo -mdout {0}/{1}.mdout".format(subdir, tempname)
    if add_crd:
        toadd += " -rst {0}/{1} -crd {0}/{1}.dat -box {0}/{1}.box".format(subdir, tempname)
    return toadd

def _mmgbsa_min(args):
    """do minimization"""
    if "min" in args.do:
        if os.path.exists("run/min"):
            shutil.rmtree("run/min")
        os.mkdir("run/min")
        basic = f"SPONGE -default_in_file_prefix run/{args.temp} -device {args.device}"
        basic += _mmgbsa_output_path("run/min", args.temp)
        if args.sr:
            basic += f" -restrain_atom_id run/{args.temp}_restrain_atom_id.txt -restrain_weight {args.rw}"
        dt_factor = 1e-2
        inc_rate = 1.5
        if not args.mi:
            basic += " -mode minimization -cutoff 8"
            cif = " -minimization_dynamic_dt 1"
            exit_code = run(f"{basic} {cif} -step_limit {args.min_step} \
-minimization_dt_factor {dt_factor} -minimization_dt_increasing_rate {inc_rate}")
            out = MdoutReader(f"run/min/{args.temp}.mdout").potential[-1]
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
                out = MdoutReader(f"run/min/{args.temp}.mdout").potential[-1]
            if min_time >= 10:
                Xprint("Minimization has been repeated for 10 times and \
the potential still can not be reduced to 0", "ERROR")
                sys.exit(1)
        else:
            command += f" -mdin {args.mi}"
            exit_code = run(command)
            if exit_code != 0:
                Xprint(f"The minimization exited with code {exit_code}", "ERROR")
                sys.exit(exit_code)


def _mmgbsa_pre_equilibrium(args):
    """do pre_equilibrium simulation"""
    if "pre_equilibrium" in args.do:
        if os.path.exists("run/pre_equilibrium"):
            shutil.rmtree("run/pre_equilibrium")
        os.mkdir("run/pre_equilibrium")
        command = f"SPONGE -default_in_file_prefix run/{args.temp} -device {args.device}"
        command += _mmgbsa_output_path("run/pre_equilibrium", args.temp)
        if args.sr:
            command += f" -restrain_atom_id run/{args.temp}_restrain_atom_id.txt -restrain_weight {args.rw}"
        command += f" -coordinate_in_file run/min/{args.temp}_coordinate.txt"
        if not args.pi:
            command += f" -step_limit {args.pre_equilibrium_step}"
            command += " -cutoff 8"
            command += f" -dt {args.dt} -constrain_mode SHAKE"
            command += " -thermostat middle_langevin"
            if not args.nvt:
                command += " -mode NPT -barostat andersen_barostat"
            else:
                command += " -mode NVT"
            command += " -middle_langevin_gamma 10 -velocity_max 20"
            exit_code = run(command)
        else:
            command += f" -mdin {args.pi}"
            exit_code = run(command)
        if exit_code != 0:
            Xprint(f"The pre_equilibrium exited with code {exit_code}", "ERROR")
            sys.exit(exit_code)

def _mmgbsa_equilibrium(args):
    """do equilibrium production simulation"""
    if "equilibrium" in args.do:
        if os.path.exists("run/equilibrium"):
            os.system("rm -rf run/equilibrium")
        os.mkdir("run/equilibrium")
        command = f"SPONGE -default_in_file_prefix run/{args.temp} -device {args.device}"
        command += _mmgbsa_output_path("run/equilibrium", args.temp)
        command += f" -coordinate_in_file run/pre_equilibrium/{args.temp}_coordinate.txt"
        command += f" -velocity_in_file run/pre_equilibrium/{args.temp}_velocity.txt"
        if args.sr:
            command += f" -restrain_atom_id run/{args.temp}_restrain_atom_id.txt -restrain_weight {args.rw}"
        if not args.ei:
            command += f" -step_limit {args.equilibrium_step} -dt {args.dt} -constrain_mode SHAKE  -cutoff 8"
            command += " -thermostat middle_langevin"
            if not args.nvt:
                command += " -mode NPT -barostat andersen_barostat"
            else:
                command += " -mode NVT"
            command += " -middle_langevin_gamma 10 -velocity_max 20"
            command += f" -write_information_interval 100 -write_restart_file_interval {args.equilibrium_step}"
            exit_code = run(command)
        else:
            command += f" -mdin {args.pi}"
            exit_code = run(command)
        if exit_code != 0:
            Xprint(f"The equilibrium exited with code {exit_code}", "ERROR")
            sys.exit(exit_code)

def _mmgbsa_analysis(args):
    """do analysis"""
    if "analysis" in args.do:
        u = mda.Universe(f"run/{args.temp}_mass.txt", f"run/equilibrium/{args.temp}.dat",
                         box=f"run/equilibrium/{args.temp}.box", format=SpongeTrajectoryReader)
        r1 = u.select_atoms(args.s1)
        r2 = u.select_atoms(args.s2)
        complex_ = r1 + r2
        r1.write(f"part1/{args.temp}.dat", "SPONGE_TRAJ", frames="all")
        r2.write(f"part2/{args.temp}.dat", "SPONGE_TRAJ", frames="all")
        complex_.write(f"complex/{args.temp}.dat", "SPONGE_TRAJ", frames="all")

        for i in ["complex", "part1", "part2"]:
            exit_code = run(f"SPONGE_NOPBC -mode rerun -crd {i}/TMP.dat -device {args.device} \
-box {i}/TMP.box -default_in_file_prefix {i}/TMP -mdout {i}/TMP.mdout -mdinfo {i}/TMP.mdinfo \
-cutoff 999")
            if exit_code != 0:
                Xprint(f"The analysis of {i} exited with code {exit_code}", "ERROR")
                sys.exit(exit_code)
        with open("free_energy.txt", "w") as f:
            complex_ene = MdoutReader("complex/TMP.mdout")
            r1_ene = MdoutReader("part1/TMP.mdout")
            r2_ene = MdoutReader("part2/TMP.mdout")
            delta_ene = complex_ene.potential - r1_ene.potential - r2_ene.potential
            delta_gb = complex_ene.gb - r1_ene.gb - r2_ene.gb
            delta_lj = complex_ene.LJ - r1_ene.LJ - r2_ene.LJ
            delta_ee = complex_ene.Coulomb - r1_ene.Coulomb - r2_ene.Coulomb
            f.write("total\t\tgb\t\tLJ\t\tCoulomb\n")
            f.write(f"{np.mean(delta_ene):.2f} +- {np.std(delta_ene):.2f}\t")
            f.write(f"{np.mean(delta_gb):.2f} +- {np.std(delta_gb):.2f}\t")
            f.write(f"{np.mean(delta_lj):.2f} +- {np.std(delta_lj):.2f}\t")
            f.write(f"{np.mean(delta_ee):.2f} +- {np.std(delta_ee):.2f}\t")
