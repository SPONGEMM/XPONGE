"""
This **module** contains the functions for mbar analysis
"""
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from .. import kb
from ..mdrun import run
from ..helper import Xopen
from ..analysis import MdoutReader

try:
    from gamda.special import MBAR
except ImportError:
    from pymbar import MBAR as _PyMBAR

    class _ArrayWithGet:
        """A light wrapper to mimic cupy-like `.get()` on numpy arrays."""

        def __init__(self, arr):
            self._arr = np.asarray(arr)

        def __array__(self, dtype=None):
            return np.asarray(self._arr, dtype=dtype)

        def __getitem__(self, key):
            return _ArrayWithGet(self._arr[key])

        def get(self):
            return self._arr

    class MBAR:  # pylint: disable=too-few-public-methods
        """Fallback MBAR implementation based on pymbar."""

        def __init__(self, enes, weights, beta, _n_bootstrap):
            self._enes = np.asarray(enes, dtype=np.float64)
            self._weights = np.asarray(weights, dtype=np.float64)
            self._beta = float(beta)
            self.f = None

        def run(self):
            n_lambda, _, frame = self._enes.shape
            n_total = n_lambda * frame
            u_kn = np.zeros((n_lambda, n_total), dtype=np.float64)
            n_k = np.full(n_lambda, frame, dtype=np.int64)
            # Flatten energies by trajectory blocks: (i, t) -> i * frame + t
            for i in range(n_lambda):
                start = i * frame
                end = start + frame
                w = self._weights[i]
                if w.ndim == 0:
                    w = np.full(frame, float(w), dtype=np.float64)
                # A simple reweighting-factor correction in reduced-energy space.
                # weight is expected positive.
                corr = -np.log(np.clip(w, 1e-12, None))
                for j in range(n_lambda):
                    u_kn[j, start:end] = self._beta * self._enes[i, j, :] + corr

            mbar = _PyMBAR(u_kn, n_k, verbose=False)
            delta_f = mbar.compute_free_energy_differences()["Delta_f"]
            # Store kcal/mol free energies relative to state 0 as shape (1, n_lambda).
            f0 = delta_f[0, :] / self._beta
            self.f = _ArrayWithGet(f0.reshape(1, -1))

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
