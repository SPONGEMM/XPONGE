"""
This **module** contains the functions for mbar analysis
"""
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from ..mdrun import run
from ..helper import Xopen
from ..analysis import MdoutReader

def _get_observe_sits_atoms(args, merged_from, ofile):
    count = 0
    for i, r in enumerate(merged_from.residues):
        if i == args.ri:
            count = [str(i) for i in range(count, count + len(r.atoms))]
            break
        else:
            count += len(r.atoms)
    with open(ofile, "w") as f:
        f.write("\n".join(count))

def _rerun_ith_traj_with_jth_forcefield(args, i, j):
    """rerun the i-th trajecotry using the j-th forcefield"""
    if os.path.exists("%d/mbar/%d"%(i, j)):
        shutil.rmtree("%d/mbar/%d"%(i, j))
    os.mkdir("%d/mbar/%d"%(i, j))
    lambda_ = args.l[j]
    command = f"SPONGE -mode rerun -default_in_file_prefix {j}/{args.temp} "
    command += f" -crd {i}/equilibrium/{args.temp}.dat -box {i}/equilibrium/{args.temp}.box -lambda_lj {lambda_} "
    command += " -SITS_atom_in_file sits_atom_in_file.txt -SITS_cross_enhancing_factor 1 -SITS_mode observation "
    command += f" -mdinfo {i}/mbar/{j}/{args.temp}.mdinfo -mdout {i}/mbar/{j}/{args.temp}.mdout"
    if not args.ai:
        command += f" -neighbor_list_max_atom_in_grid_numbers 128 -neighbor_list_max_neighbor_numbers 1200 -cutoff 8"
        exit_code = run(command)
    else:
        command += f" -mdin {args.ai}"
        exit_code = run(command)
    assert exit_code == 0, f"Wrong for rerun Trajectory {i} using Forcefiled {j}"

def _get_neighbors_of_i(args, i, n_neighbor=1)
    """get the neighbor lambda of the i-th trajecotry""" 
    for j in range(min(0, i - n_neighbor), max(args.nl + 1, i + n_neighbor), 1):
        yield j

def mbar_analysis(args, merged_from):
    """
    This **function** is used to do the mbar analysis

    :param args: the arguments from the command line
    :return: None
    """
    _get_observe_sits_atoms(args, merged_from, "sits_atom_in_file.txt")
    frame = args.equilibrium_step // args.wi
    enes = np.zeros((args.nl, args.nl, frame))
    for i in range(args.nl + 1):
        if os.path.exists("%d/mbar" % i):
            shutil.rmtree("%d/mbar" % i)
        if os.path.exists("%d/equilibrium/reweighting_factor.txt" % i):
            weight = np.loadtxt("%d/equilibrium/reweighting_factor.txt" % i, dtype=np.float128).reshape(-1)
        else:
            weight = np.ones(frame, dtype=float)
        os.mkdir("%d/mbar" % i)
        for j in _get_neighbors_of_i(args, i):
            _rerun_ith_traj_with_jth_forcefield(args, i, j)
            enes[i][j][:] = MdoutReader(f"{i}/mbar/{j}/{args.temp}.mdout").SITS_AA_kAB
    
