"""
    This **module** helps with the definition of CVs in SPONGE
"""
from pathlib import Path
from collections.abc import Iterable
import numpy as np
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
    vatom_type = center_of_mass
    {atom_line}
}}
"""

class _Position(_CV):
    """ cv: position """
    def __init__(self, name, atom, xyz, scaled):
        self.name = name
        self.atom = atom
        self.xyz = xyz
        if not scaled:
            self.scaled = ""
        else:
            self.scaled = "scaled_"
            self.period = 1

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = {self.scaled}position_{self.xyz}
    atom = {self.atom}
}}
"""

class _Boxlength(_CV):
    """ cv: boxlength """
    def __init__(self, name, xyz):
        self.name = name
        self.xyz = xyz

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = box_length_{self.xyz}
}}
"""

class _Distance(_CV):
    """ cv: distance """
    def __init__(self, name, atom1, atom2):
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = distance
    atom = {self.atom1} {self.atom2}
}}
"""

class _Displacement(_CV):
    """ cv: displacement """
    def __init__(self, name, atom1, atom2, xyz):
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.xyz = xyz

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = displacement_{self.xyz}
    atom = {self.atom1} {self.atom2}
}}
"""

class _Angle(_CV):
    """ cv: angle """
    period = 2 * np.pi
    def __init__(self, name, atom1, atom2, atom3):
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = angle
    atom = {self.atom1} {self.atom2} {self.atom3}
}}
"""

class _Dihedral(_CV):
    """ cv: dihedral """
    period = 2 * np.pi
    def __init__(self, name, atom1, atom2, atom3, atom4):
        self.name = name
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def to_string(self, folder):
        return f"""{self.name}
{{
    CV_type = dihedral
    atom = {self.atom1} {self.atom2} {self.atom3} {self.atom4}
}}
"""

class _RMSD(_CV):
    """ cv: rmsd """
    def __init__(self, name, atom, coordinate):
        self.name = name
        self.atom = atom
        self.coordinate = coordinate

    def to_string(self, folder):
        prefix = Path(folder)
        with open(prefix / (self.name + "_atom.txt"), "w") as f:
            f.write("\n".join(map(str, self.atom)))
        np.savetxt(prefix / (self.name + "_coordinate.txt"), self.coordinate)
        atom_line = f"atom_in_file = {str(prefix / (self.name + '_atom.txt'))}"
        crd_line = f"coordinate_in_file = {str(prefix / (self.name + '_coordinate.txt'))}"
        return f"""{self.name}
{{
    CV_type = rmsd
    {atom_line}
    {crd_line}
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

class _RestrainCV(_CVBias):
    """ bias: steer """
    def __init__(self, cv, weight, reference):
        if isinstance(cv, Iterable):
            self.cv = cv
            self.weight = weight
            self.reference = reference
        else:
            self.cv = [cv]
            self.weight = [weight]
            self.reference = [reference]

    def to_string(self, folder):
        need_period = False
        period_line = " "
        periods = []
        for cv in self.cv.values():
            if hasattr(cv, "period"):
                need_period = True
            periods.append(str(getattr(cv, "period", 0)))
        if need_period:
            period_line = "    " + " ".join(periods) + "\n"
        return f"""restrain
{{
    CV = {" ".join([str(cv) for cv in self.cv])}
    weight = {" ".join([str(weight) for weight in self.weight])}
    reference = {" ".join([str(reference) for reference in self.reference])}
{period_line}}}
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
        if weight == "mass":
            self.virtual_atom[name] = _COM(name, atom)
        else:
            self.virtual_atom[name] = _Center(name, atom, weight)
        self.names[name] = self.virtual_atom[name]
        self._association[name] = []

    def add_cv_position(self, name, atom, xyz, scaled):
        """
            Add a CV with the type of "position" to the system

            :param name: the name of the CV
            :param atom: an int or a name of virtual atom
            :param xyz: the axis of the position
            :param scaled: whether the posithion need to be scaled
            :return: None
        """
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        if xyz not in ("x", "y", "z"):
            raise ValueError(f"xyz should be one of 'x', 'y' or 'z', but {xyz} is given")
        self.cv[name] = _Position(name, atom, xyz, scaled)
        self.names[name] = self.cv[name]
        self._association[name] = []
        if atom in self.virtual_atom:
            self._association[name].append(self.names[name])
        elif not isinstance(atom, int):
            raise TypeError(f"atom should be an int or a name of virtual atom, but {atom} is given")


    def add_cv_distance(self, name, atom1, atom2):
        """
            Add a CV with the type of "distance" to the system

            :param name: the name of the CV
            :param atom1: an int or a name of virtual atom
            :param atom2: an int or a name of virtual atom
            :return: None
        """
        self.cv[name] = _Distance(name, atom1, atom2)
        self.names[name] = self.cv[name]
        self._association[name] = []
        for atom in [atom1, atom2]:
            if atom in self.virtual_atom:
                self._association[atom].append(self.names[name])
            elif not isinstance(atom, int):
                raise TypeError(f"atom should be an int or a name of virtual atom, but {atom1} is given")

    def add_cv_displacement(self, name, atom1, atom2, xyz):
        """
            Add a CV with the type of "displacement" to the system

            :param name: the name of the CV
            :param atom1: an int or a name of virtual atom
            :param atom2: an int or a name of virtual atom
            :param xyz: the axis of the position
            :return: None
        """
        if xyz not in ("x", "y", "z"):
            raise ValueError(f"xyz should be one of 'x', 'y' or 'z', but {xyz} is given")
        self.cv[name] = _Displacement(name, atom1, atom2, xyz)
        self.names[name] = self.cv[name]
        self._association[name] = []
        for atom in [atom1, atom2]:
            if atom in self.virtual_atom:
                self._association[atom].append(self.names[name])
            elif not isinstance(atom, int):
                raise TypeError(f"atom should be an int or a name of virtual atom, but {atom1} is given")

    def add_cv_box_length(self, name, xyz):
        """
            Add a CV with the type of "box_length" to the system

            :param name: the name of the CV
            :param xyz: the axis of the position
            :return: None
        """
        if xyz not in ("x", "y", "z"):
            raise ValueError(f"xyz should be one of 'x', 'y' or 'z', but {xyz} is given")
        self.cv[name] = _Boxlength(name, xyz)
        self.names[name] = self.cv[name]
        self._association[name] = []

    def add_cv_angle(self, name, atom1, atom2, atom3):
        """
            Add a CV with the type of "angle" to the system

            :param name: the name of the CV
            :param atom1: an int or a name of virtual atom
            :param atom2: an int or a name of virtual atom
            :param atom3: an int or a name of virtual atom
            :return: None
        """
        self.cv[name] = _Angle(name, atom1, atom2, name3)
        self.names[name] = self.cv[name]
        self._association[name] = []
        for atom in [atom1, atom2, atom3]:
            if atom in self.virtual_atom:
                self._association[atom].append(self.names[name])
            elif not isinstance(atom, int):
                raise TypeError(f"atom should be an int or a name of virtual atom, but {atom1} is given")

    def add_cv_dihedral(self, name, atom1, atom2, atom3, atom4):
        """
            Add a CV with the type of "dihedral" to the system

            :param name: the name of the CV
            :param atom1: an int or a name of virtual atom
            :param atom2: an int or a name of virtual atom
            :param atom3: an int or a name of virtual atom
            :param atom4: an int or a name of virtual atom
            :return: None
        """
        self.cv[name] = _Dihedral(name, atom1, atom2, atom3, atom4)
        self.names[name] = self.cv[name]
        self._association[name] = []
        for atom in [atom1, atom2, atom3, atom4]:
            if atom in self.virtual_atom:
                self._association[atom].append(self.names[name])
            elif not isinstance(atom, int):
                raise TypeError(f"atom should be an int or a name of virtual atom, but {atom1} is given")

    def add_cv_rmsd(self, name, select):
        """
            Add a CV with the type of "displacement" to the system

            :param name: the name of the CV
            :param select: a string of selection 
            :return: None
        """
        if name in self.names:
            raise ValueError(f"{name} has been defined in the name system")
        ag = self.u.select_atoms(select)
        atom = [self.id2index[atom.id] for atom in ag]
        coordinate = ag.positions
        self.cv[name] = _RMSD(name, atom, coordinate)
        self.names[name] = self.cv[name]
        self._association[name] = []

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

    def restrain(self, name, weight, reference):
        """
            Add a CV to restrain

            :param name: the name of the CV
            :param weight: the weight for steering
            :param reference: the reference for steering
            :return: None
        """
        if "restrain" not in self.bias:
            self.bias["restrain"] = _SteerCV(self.cv[name], weight, reference)
            self.names["restrain"] = self.bias["restrain"]
        else:
            self.bias["restrain"].cv.append(self.cv[name])
            self.bias["restrain"].weight.append(weight)
            self.bias["restrain"].reference.append(reference)

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

