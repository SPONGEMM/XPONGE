"""
This **package** sets the basic configuration of Stillinger-Weber force field
"""
import os
from ..base import sw_base, charge_base, mass_base
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
        toread = "name epsilon sigma a l gamma b A B p q\n"
    else:
        toread = "name epsilon[eV] sigma a l gamma b A B p q\n"
    filename = os.path.join(folder, filename)
    values = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue
            values += line.split()
    for i in range(0, len(values), 13):
        toread += "-".join(values[i: i + 2]) + " " + " ".join(values[i + 3: i + 13]) + "\n"
        toread += "-".join(values[i: i + 3]) + " " + " ".join(values[i + 3: i + 13]) + "\n"
    if output:
        Xprint(toread)
    sw_base.SWType.New_From_String(toread)
