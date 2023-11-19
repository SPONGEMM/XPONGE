"""
    This **module** helps with the definition of CVs in SPONGE
"""
from pathlib import Path
from collections import Iterable
from . import Xdict


class _CVVirtualAtom:
    """ an meta class for virtual atom in CV system"""
    def __str__(self):
        return self.name

class _CV:
   """ an meta class for CV in CV system """
   def __str__(self):
        return self.name

class _CVBias:
   """ an meta class for bias in CV system """

class _Center(_CVVirtualAtom):
    """ cv virtual atom: center """
    def __init__(self, name, atom, weight):
        self.name = name
        self.atom = atom
        self.weight = weight

    def to_string(self, folder):
        prefix = Path(folder)
        if len(self.atom) > 10:
            with open(prefix / (self.name + "_atom.txt"), "w") as f:
                f.write("\n".join(map(str, self.atom)))
            with open(prefix / (self.name + "_weight.txt"), "w") as f:
                f.write("\n".join(map(str, self.weight)))
            atom_line = f"atom_in_file = {str(prefix / (self.name + '_atom.txt'))}"
            weight_line = f"weight_in_file = {str(prefix / (self.name + '_weight.txt'))}"
        else:
            atom_line = " ".join(map(str, self.atom))
            weight_line = " ".join(map(str, self.weight))
        return f"""{self.name}
{{
    vatom_type = center
    {atom_line}
    {weight_line} 
}}
"""

class _COM(_CVVirtualAtom):
    """ cv virtual atom: center """
    def __init__(self, name, atom):
        self.name = name
        self.atom = atom

    def to_string(self, folder):
        prefix = Path(prefix)
        if len(self.atom) > 10:
            with open(prefix / (self.name + "_atom.txt"), "w") as f:
                f.write("\n".join(map(str, self.atom)))
            atom_line = f"atom_in_file = {str(prefix / (self.name + '_atom.txt'))}"
        else:
            atom_line = " ".join(map(str, self.atom))
        return f"""{self.name}
{{
    vatom_type = COM
    {atom_line}
}}
"""

class _PositionX(_CV):
    """ cv: position """
    def __init__(self, name, atom):
        self.name = name
        self.atom = atom

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = position_x
    atom = {self.atom}
}}
"""

class _PrintCV(_CVBias):
    """ bias: print """
    def __init__(self, cv):
        if isinstance(cv, Iterable):
            self.cv = cv
        else:
            self.cv = [cv]

    def to_string(self, folder):
        return f"""print
{{
    CV = {" ".join([str(cv) for cv in self.cv])}
}}
"""

class _SteerCV(_CVBias):
    """ bias: steer """
    def __init__(self, cv, weight):
        if isinstance(cv, Iterable):
            self.cv = cv
            self.weight = weight
        else:
            self.cv = [cv]
            self.weight = [weight]

    def to_string(self, folder):
        return f"""steer
{{
    CV = {" ".join([str(cv) for cv in self.cv])}
    weight = {" ".join([str(weight) for weight in self.weight])}
}}
"""

class CVSystem:
    """
        This **class** is used to help with the definition of CVs in SPONGE

        :param molecule: the molecule to define the CVs
    """
    __slots__ = ["_u", "_id2index", "molecule", "virtual_atom", "cv", "bias", "names", "_association"]
    BIASES = {"print", "restrain", "steer", "meta1d"}
    SPONGE_NAMES = {"bond", "angle", "dihedral"}
    def __init__(self, molecule):
        self._u = None
        self._id2index = None
        self._association = Xdict(not_found_message="No name {} found in the system")
        self.molecule = molecule
        self.virtual_atom = Xdict(not_found_message="No virtual atom named {} found in the system")
        self.cv = Xdict(not_found_message="No collected variable named {} found in the system")
        self.bias = Xdict(not_found_message="No bias named {} found in the system")
        self.names = Xdict(not_found_message="No name {} found in the system")
        self.names.fromkeys(self.BIASES, [])
        self.names.fromkeys(self.SPONGE_NAMES, [])

    @property
    def u(self):
        """ the MDAnalysis.Universe of the molecule """
        if self._u is None:
            from ..analysis.md_analysis import mda, XpongeMoleculeReader
            self._u = mda.Universe(self.molecule, format=XpongeMoleculeReader)
        return self._u

    @property
    def id2index(self):
        """ the MDAnalysis.Universe of the molecule """
        if self._id2index is None:
            self._id2index = {atom.id : i for i, atom in enumerate(self.u.atoms)}
        return self._id2index

    def remove(self, name):
        """
           Remove a name from the system
        """
        return NotImplementedError

    def add_center(self, name, select, weight=None):
        """
            Add a virtual atom with the type of "center" to the system

            :param name: the name of the virtual atom
            :param select: a selection string of MDAnalysis
            :param weight: weight of the atoms, None for 1/N
            :return: None
        """
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        atom = [self.id2index[atom.id] for atom in self.u.select_atoms(select)]
        if weight is None:
            weight = [1.0 / len(atom)] * len(atom)
        self.virtual_atom[name] = _Center(name, atom, weight)
        self.names[name] = self.virtual_atom[name]
        self._association[name] = []

    def add_cv_position(self, name, atom, xyz):
        """
            Add a CV with the type of "position" to the system

            :param name: the name of the CV
            :param atom: an int or a name of virtual atom
            :param xyz: the axis of the position
            :return: None
        """
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        self.cv[name] = _PositionX(name, atom)
        self.names[name] = self.cv[name]
        self._association[name] = []
        if atom in self.virtual_atom:
            self._association[name].append(self.names[name])
        elif not isinstance(atom, int):
            raise TypeError(f"atom should be an int or a name of virtual atom, but {atom} is given")

    def print(self, name):
        """
            Add a CV to print

            :param name: the name of the CV
            :return: None
        """
        if "print" not in self.bias:
            self.bias["print"] = _PrintCV(self.cv[name])
            self.names["print"] = self.bias["print"]
        else:
            self.bias["print"].cv.append(self.cv[name])
        self._association[name].append(self.bias["print"])

    def steer(self, name, weight):
        """
            Add a CV to steer

            :param name: the name of the CV
            :param weight: the weight for steering
            :return: None
        """
        if "steer" not in self.bias:
            self.bias["steer"] = _SteerCV(self.cv[name], weight)
            self.names["steer"] = self.bias["steer"]
        else:
            self.bias["steer"].cv.append(self.cv[name])
            self.bias["steer"].weight.append(weight)

    def output(self, filename, folder="."):
        """
            Output the recorded system to a file

            :param filename: the name of the output file for cv_in_file
            :param folder: the folder of the output files, the current working folder for default
            :return: None
        """
        with open(filename, "w") as f:
            if self.virtual_atom:
                f.write("#" * 30 + "\n#definition of virtual atoms\n" + "#" * 30 + "\n")
                for virtual_atom in self.virtual_atom.values():
                    f.write(virtual_atom.to_string(folder))
            if self.cv:
                f.write("#" * 30 + "\n#definition of collected variables\n" + "#" * 20 + "\n")
                for cv in self.cv.values():
                    f.write(cv.to_string(folder))
            if self.bias:
                f.write("#" * 30 + "\n#definition of bias\n" + "#" * 30 + "\n")
                for bias in self.bias.values():
                    f.write(bias.to_string(folder))

