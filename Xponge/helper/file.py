"""
    This **module** helps to process the files of molecular modelling
"""
import io
import re
import sys
from pathlib import Path
from importlib import import_module
from . import set_global_alternative_names

__all__ = ["file_filter", "pdb_filter", "import_python_script"]

def import_python_script(path):
    if not isinstance(path, Path):
        path = Path(path)
    sys.path.append(str(path.parent))
    if path.suffix != ".py":
        raise TypeError(f"{path} should be a python script")
    import_module(path.stem)

def file_filter(infile, outfile, reg_exp, replace_dict):
    """
        This **function** finds the lines which contains any of the given regular expressions and replace some parts.

        :param infile: the input file or filename
        :param outfile: the output file or filename
        :param reg_exp: a list of regular expressions. Lines which match any regular expressions will be kept.
        :param replace_dict: a dict of regular expressions and the replacement
    """
    if not isinstance(reg_exp, list):
        raise TypeError('reg_exp should be a list of regular expressions')
    if not isinstance(replace_dict, dict):
        raise TypeError('replace_dict should be a dict of regular expressions and the replacement')
    if not isinstance(infile, io.IOBase):
        infile = open(infile, "r")
    lines = ""
    with infile as f:
        for line in infile:
            for keyword in reg_exp:
                if not isinstance(keyword, str):
                    raise TypeError('reg_exp should be a list of regular expressions')
                if re.match(keyword, line):
                    for reg, rep in replace_dict.items():
                        line = re.sub(reg, rep, line)
                    lines += line
                    break
    if not isinstance(outfile, io.IOBase):
        with open(outfile, "w") as f:
            f.write(lines)
    else:
        outfile.write(lines)


def pdb_filter(infile, outfile, heads, hetero_residues, chains=None, rename_ions=None):
    """
        This **function** finds the lines in pdb which meets the need

        :param infile: the input file or filename
        :param outfile: the output file or filename
        :param head: a list of heads which will be included
        :param hetero_residues: a list of hetero residue names which will be included
        :param chains: a list of the code for the chains you need. None for all (default).
        :param rename_ions: a dict to rename the ions
    """
    if not isinstance(heads, list):
        raise TypeError("heads should be a list")
    if not isinstance(hetero_residues, list):
        raise TypeError("hetero_residues should be a list")
    if rename_ions is None:
        rename_ions = {}
    if not isinstance(rename_ions, dict):
        raise TypeError('replace_dict should be a dict of regular expressions and the replacement')
    replace_dict = {}
    for a, b in rename_ions.items():
        if len(a) == 1:
            aname = f"{a}   | {a}  |  {a} |   {a}"
            rname = f"{a}  | {a} |  {a}"
        elif len(a) == 2:
            aname = f"{a}  | {a} |  {a}"
            rname = f" {a}|{a} "
        elif len(a) == 3:
            aname = f"{a} | {a}"
            rname = a
        else:
            raise ValueError("The ion name in a pdb file should not be longer than 3 characters")
        replace_dict["(^HETATM [ 0-9]{4} )(%s)(.)(%s)"%(aname, rname)] = r"\g<1>" + f"{b:4s}" + r"\g<3>" + f"{b:3s}"
    reg_exp = []
    for head in heads:
        if head == "ATOM" and chains is not None:
            for chain in chains:
                reg_exp.append("^ATOM.{17}%s"%chain)
        elif head == "SEQRES" and chains is not None:
            for chain in chains:
                reg_exp.append("^SEQRES.{5}%s"%chain)
        elif head == "TER" and chains is not None:
            reg_exp.append(r"^TER\s*$")
            for chain in chains:
                reg_exp.append("^TER.{18}%s"%chain)
        else:
            reg_exp.append(f"^{head}")
    for hetres in hetero_residues:
        if len(hetres) == 1:
            hetres = f"{hetres}  | {hetres} |  {hetres}"
        elif len(hetres) == 2:
            hetres = f" {hetres}|{hetres} "
        reg_exp.append("^HETATM.{11}%s"%(hetres))
    file_filter(infile, outfile, reg_exp, replace_dict)

set_global_alternative_names()
