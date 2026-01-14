"""
This **package** sets the basic configuration of the Environment-Dependent Interatomic Potential (EDIP)
"""
import os
from ..base import edip_base, charge_base, mass_base
from ...helper import Xprint

def load_parameter_from_lammps(filename, folder, units, output=False):
    """
    This **function** is used to get edit force field parameters from LAMMPS

    :param filename: the name of the input file
    :param prefix: the folder of the file
    :param units: the unit of LAMMPS parameters
    :param output: whether output the converted Xponge string
    :return: None
    """
    if units not in ("real", "metal"):
        raise NotImplementedError("The units can only be 'real' or 'metal' now")
    elif units == "real":
        toread = ["name A  B  a  c alpha   beta   eta gamma   l    mu   rho   sigma   Q0  u1  u2  u3  u4\n"]
    else:
        toread = ["name A[eV] B a c alpha   beta   eta gamma   l[eV]  mu   rho   sigma   Q0  u1   u2   u3   u4\n"]
    filename = os.path.join(folder, filename)
    values = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            values.extend(line.split())
    for i in range(0, len(values), 20):
        toread.append("-".join(values[i: i + 2]) + " " + " ".join(values[i + 3: i + 20]) + "\n")
        toread.append("-".join(values[i: i + 3]) + " " + " ".join(values[i + 3: i + 20]) + "\n")
    toread = "".join(toread)
    if output:
        Xprint(toread)
    edip_base.EDIPType.New_From_String(toread)
