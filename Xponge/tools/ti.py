"""
This **module** contains the functions for ti analysis
"""
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from gamda.special import TI
from ..mdrun import run
from ..helper import Xopen
from ..analysis import MdoutReader

def _rerun_ith_traj_to_obtain_dh_dlambda(args, i):
    """rerun the i-th trajecotry to obtain the dH/dlambda"""
    if os.path.exists("%d/ti" % i):
        shutil.rmtree("%d/ti" % i)
    os.mkdir("%d/ti" % i)
    inprefix = f"{i}/{args.temp}"
    command = f"SPONGE_TI -mass_in_file {i}/{args.temp}_mass.txt -LJ_soft_core_in_file {inprefix}_LJ_soft_core.txt"
    command += " -exclude_in_file {0}_exclude.txt -charge_in_file {0}_charge.txt".format(inprefix)
    command += f" -chargeA_in_file 0/{args.temp}_charge.txt"
    command += f" -chargeB_in_file {args.nl}/{args.temp}_charge.txt"
    lambda_ = args.l[i]
    command += f" -lambda_lj {lambda_}"
    command += f" -subsys_division_in_file {inprefix}_subsys_division.txt  -charge_perturbated {args.cp}"
    inprefix = f"{i}/ti/{args.temp}"
    command += f" -mdinfo {inprefix}.mdinfo -mdout {inprefix}.mdout"
    inprefix = f"{i}/equilibrium/{args.temp}"
    command += f" -crd {inprefix}.dat -box {inprefix}.box"
    if not args.ai:
        command += " -neighbor_list_max_atom_in_grid_numbers 128"
        command += " -neighbor_list_max_neighbor_numbers 1200 -cutoff 8"
        run(command)
    else:
        command += f" -mdin {args.ai}"
        run(command)

def ti_analysis(args):
    """
    This **function** is used to do the ti analysis

    :param args: the arguments from the command line
    :return: None
    """
    frame = args.equilibrium_step // args.wi
    n_lambda = args.nl + 1
    dh_dlambda = np.zeros((n_lambda, frame), dtype=np.float32)
    weights = np.zeros((n_lambda, frame), dtype=np.float32)
    space = np.zeros(args.nl, dtype=np.float32)
    for i in range(n_lambda):
        if os.path.exists("%d/equilibrium/reweighting_factor.txt" % i):
            weight = np.loadtxt("%d/equilibrium/reweighting_factor.txt" % i, dtype=np.longdouble).reshape(-1)
        else:
            weight = np.ones(frame, dtype=float)
        weights[i][:] = weight
        if not args.nar:
            _rerun_ith_traj_to_obtain_dh_dlambda(args, i)
        dh_dlambda[i][:] = MdoutReader(f"{i}/ti/{args.temp}.mdout").dH_dlambda
        if i != args.nl:
            space[i] = (args.l[i+1] - args.l[i]) / 2
    ti = TI(dh_dlambda, weights, space, 1024)
    ti.run()
    fe = np.mean(ti.f, axis=0)
    error = np.std(ti.f, axis=0)
    f = Xopen("TI.txt", "w")
    f.write("lambda_state\tFE(i+1)-FE(i)[kcal/mol]\tFE(i+1)-FE(0)[kcal/mol]\tSigma(FE(i+1)-FE(0))[kcal/mol]\n")
    f.write("\n".join(
        [f"{i}\t\t{fe[i+1] - fe[i]: .2f}\t\t\t{fe[i+1]: .2f}\t\t\t{error[i+1]:.2f}" for i in range(args.nl)]))
    f.close()
    ans = ti.f[:, -1].get()
    kernel = gaussian_kde(ans, bw_method=1)
    x = np.linspace(np.min(ans), np.max(ans), 1024)
    y = kernel(x)
    plt.plot(x, y)
    plt.xlabel("Free Energy [kcal/mol]")
    plt.ylabel("Probability")
    plt.savefig("TI.png")
    plt.clf()
