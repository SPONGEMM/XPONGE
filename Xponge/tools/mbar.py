"""
This **module** contains the functions for mbar analysis
"""
import os
import shutil
import numpy as np
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
    enes = np.zeros((n_lambda, n_lambda, frame))
    bias = np.zeros((n_lambda, frame))
    for i in range(n_lambda):
        if os.path.exists("%d/mbar" % i):
            shutil.rmtree("%d/mbar" % i)
        if os.path.exists("%d/equilibrium/reweighting_factor.txt" % i):
            weight = np.loadtxt("%d/equilibrium/reweighting_factor.txt" % i, dtype=np.longfloat).reshape(-1)
        else:
            weight = np.ones(frame, dtype=float)
        bias[i][:] = np.log(weight) / beta
        os.mkdir("%d/mbar" % i)
        for j in range(0, n_lambda):
            _rerun_ith_traj_with_jth_forcefield(args, i, j)
            mdout = MdoutReader(f"{i}/mbar/{j}/{args.temp}.mdout")
            enes[i][j][:] = mdout.potential
    bias = bias.reshape((n_lambda, 1, frame))
    enes -= np.min(enes - bias, axis=1, keepdims=True)
    f = np.zeros(n_lambda)
    last_f = np.ones(n_lambda)
    while np.any(np.abs(last_f - f) > 0.001):
        last_f = f
        sum_j = np.sum(np.exp(-beta * (enes - bias - f.reshape((1, -1, 1)))), axis=1, keepdims=True)
        f = -np.log(np.sum(np.exp(-beta * enes)/ sum_j, axis=(0, 2))) / beta
        f -= np.min(f)
    theta_i = -beta * np.diagonal(enes)
    sigma_ij = np.zeros(args.nl)
    sigma_i0 = np.zeros(args.nl)
    factor = np.array([[1, -1], [-1, 1]])
    for i in range(args.nl):
        sigma_ij[i] = np.sqrt(np.sum(factor * np.cov(theta_i[i], theta_i[i+1])))
        sigma_i0[i] = np.sqrt(np.sum(factor * np.cov(theta_i[0], theta_i[i+1])))
    fw = Xopen("free_energy.txt", "w")
    fw.write("lambda_state\tFE(i+1)-FE(i)[kcal/mol]\tFE(i+1)-FE(0)[kcal/mol]\n")
    fw.write("\n".join(
        [f"{i}\t\t{f[i+1]-f[i]:.2f} +- {sigma_ij[i]:.2f}\t\t{f[i+1]-f[0]:.2f} +- {sigma_i0[i]:.2f}"
            for i in range(args.nl)]))
    fw.close()
