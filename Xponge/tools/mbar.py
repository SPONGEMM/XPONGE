"""
This **module** contains the functions for mbar analysis
"""
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from gamda.special import MBAR
from scipy.stats import gaussian_kde
from .. import kb
from ..mdrun import run
from ..helper import Xopen
from ..analysis import MdoutReader

def _rerun_ith_traj_with_jth_forcefield(args, i, j):
    """rerun the i-th trajecotry using the j-th forcefield"""
    if os.path.exists("%d/mbar/%d"%(i, j)):
        shutil.rmtree("%d/mbar/%d"%(i, j))
    os.mkdir("%d/mbar/%d"%(i, j))
    lambda_ = args.l[j]
    command = f"SPONGE -mode rerun -default_in_file_prefix {j}/{args.temp} "
    command += f" -crd {i}/equilibrium/{args.temp}.dat -box {i}/equilibrium/{args.temp}.box -lambda_lj {lambda_} "
    command += f" -mdinfo {i}/mbar/{j}/{args.temp}.mdinfo -mdout {i}/mbar/{j}/{args.temp}.mdout -PME_print_detail 1"
    if not args.ai:
        command += " -cutoff 8"
        exit_code = run(command)
    else:
        command += f" -mdin {args.ai}"
        exit_code = run(command)
    assert exit_code == 0, f"Wrong for rerun Trajectory {i} using Forcefiled {j}"

def mbar_analysis(args):
    """
    This **function** is used to do the mbar analysis

    :param args: the arguments from the command line
    :return: None
    """
    beta = 4184 / 300 / 8.314
    frame = args.equilibrium_step // args.wi
    n_lambda = args.nl + 1
    weights = np.zeros((n_lambda, frame), dtype=np.float32)
    for i in range(n_lambda):
        if os.path.exists("%d/equilibrium/reweighting_factor.txt" % i):
            weight = np.loadtxt("%d/equilibrium/reweighting_factor.txt" % i, dtype=np.float32).reshape(-1)
        else:
            weight = np.ones(frame, dtype=np.float32)
        weights[i][:] = weight
        if not args.nar:
            if os.path.exists("%d/mbar" % i):
                shutil.rmtree("%d/mbar" % i)
            os.mkdir("%d/mbar" % i)
            for j in range(n_lambda):
                _rerun_ith_traj_with_jth_forcefield(args, i, j)
    enes = np.zeros((n_lambda, n_lambda, frame), dtype=np.float32)
    for i in range(n_lambda):
        for j in range(n_lambda):
            mdout = MdoutReader(f"{i}/mbar/{j}/{args.temp}.mdout")
            enes[i][j][:] = mdout.potential
    mbar = MBAR(enes, weights, 1.0 / kb / 300, 1024)
    mbar.run()
    fe = np.mean(mbar.f, axis=0)
    error = np.std(mbar.f, axis=0)
    f = Xopen("MBAR.txt", "w")
    f.write("lambda_state\tFE(i+1)-FE(i)[kcal/mol]\tFE(i+1)-FE(0)[kcal/mol]\tSigma(FE(i+1)-FE(0))[kcal/mol]\n")
    f.write("\n".join(
        [f"{i}\t\t{fe[i+1] - fe[i]: .2f}\t\t\t{fe[i+1]: .2f}\t\t\t{error[i+1]:.2f}" for i in range(args.nl)]))
    f.close()
    ans = mbar.f[:, -1].get()
    kernel = gaussian_kde(ans, bw_method=1)
    x = np.linspace(np.min(ans), np.max(ans), 1024)
    y = kernel(x)
    plt.plot(x, y)
    plt.xlabel("Free Energy [kcal/mol]")
    plt.ylabel("Probability")
    plt.savefig("MBAR.png")
    plt.clf()
