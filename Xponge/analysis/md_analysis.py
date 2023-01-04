"""
This **module** gives functions and classes to use MDAnalysis to analyze the trajectories
"""
import os.path

import numpy as np
from .. import Residue, ResidueType
from ..helper import Xopen, guess_element_from_mass, set_global_alternative_names, set_attribute_alternative_name
try:
    import MDAnalysis as mda
    from MDAnalysis.coordinates import base
    from MDAnalysis.lib import util
    from MDAnalysis.topology.base import TopologyReaderBase
    from MDAnalysis.core import topologyattrs
    from MDAnalysis.core.topology import Topology
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "'MDAnalysis' package needed. Maybe you need 'pip install MDAnalysis'") from exc


class SpongeNoneReader(base.ReaderBase):
    def __init__(self, _, n_atoms, **kwargs):
        super().__init__(_, **kwargs)
        self._n_atoms = n_atoms

    @property
    def n_atoms(self):
        return self._n_atoms

    @property
    def n_frames(self):
        return 0

    def close(self):
        """
        fake close function for api
        """
        return


class SpongeInputReader(TopologyReaderBase):
    """
    This **class** is used to read the SPONGE input to mdanalysis

    Create the following attributes:
        Masses
        Charges
        Bonds
        Angles

    Guesses the following attributes:
        Atomnames
        Atomtypes
        Elements
    """
    #pylint: disable=unused-argument
    def parse(self, **kwargs):
        """
        This **function** reads the file and returns the structure

        :param kwargs: keyword arguments
        :return: MDAnalysis Topology object
        """
        attrs = [topologyattrs.Segids(np.array(['SYSTEM'], dtype=object))]
        nres = 1
        if os.path.exists(self.filename + "_mass.txt"):
            with util.openany(self.filename + "_mass.txt") as fm:
                natoms = int(fm.readline())
                masses = [float(line.strip()) for line in fm]
                atom_names = [guess_element_from_mass(mass) for mass in masses]
                attrs.append(topologyattrs.Masses(masses))
                attrs.append(topologyattrs.Atomnames(atom_names, guessed=True))
                attrs.append(topologyattrs.Atomtypes(atom_names, guessed=True))
                attrs.append(topologyattrs.Elements(atom_names, guessed=True))
        if os.path.exists(self.filename + "_charge.txt"):
            with util.openany(self.filename + "_charge.txt") as fm:
                natoms = int(fm.readline())
                charges = [float(line.strip()) / 18.2223 for line in fm]
                attrs.append(topologyattrs.Charges(charges))
        resid = np.zeros(natoms, dtype=np.int32)
        nres = 1
        if os.path.exists(self.filename + "_residue.txt"):
            with util.openany(self.filename + "_residue.txt") as fm:
                natoms, nres = fm.readline().split()
                natoms, nres = int(natoms), int(nres)
                resid = np.zeros(natoms, dtype=np.int32)
                count = 0
                for i, line in enumerate(fm):
                    res_length = int(line.strip())
                    resid[count:count + res_length] = i
                    count += res_length
        attrs.append(topologyattrs.Resids(np.arange(nres) + 1))
        attrs.append(topologyattrs.Atomids(np.arange(natoms) + 1))
        attrs.append(topologyattrs.Resnums(np.arange(nres) + 1))
        if os.path.exists(self.filename + "_bond.txt"):
            with util.openany(self.filename + "_bond.txt") as fm:
                fm.readline()
                bonds = [[int(words) for words in line.split()[:2]] for line in fm]
                attrs.append(topologyattrs.Bonds(bonds))
        self._n_atoms = natoms
        return Topology(natoms, nres, 1, attrs, resid, None)


class XpongeResidueReader(base.ReaderBase):
    """
    This **class** is used to read the Xponge Residue or ResidueType to mdanalysis

    Create the following attributes:
        Atomnames
        Atomtypes
        Masses
        Charges
        Bonds

    Guesses the following attributes:
        Elements
    """
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts.positions = [[getattr(atom, i) for i in "xyz"] for atom in self.filename.atoms]

    @property
    def n_atoms(self):
        return len(self.filename.atoms)

    @property
    def n_frames(self):
        return 1

    #pylint: disable=unused-argument
    def parse(self, **kwargs):
        """
        This **function** reads the file and returns the structure

        :param kwargs: keyword arguments
        :return: MDAnalysis Topology object
        """
        attrs = [topologyattrs.Segids(np.array(['SYSTEM'], dtype=object))]
        residue = self.filename
        natoms = len(residue.atoms)
        nres = 1
        masses = [atom.mass for atom in residue.atoms]
        elements = [guess_element_from_mass(mass) for mass in masses]
        attrs.append(topologyattrs.Masses(masses))
        attrs.append(topologyattrs.Elements(elements, guessed=True))
        attrs.append(topologyattrs.Atomnames([atom.name for atom in residue.atoms]))
        attrs.append(topologyattrs.Atomtypes([atom.type.name for atom in residue.atoms]))
        attrs.append(topologyattrs.Charges([atom.charge for atom in residue.atoms]))
        attrs.append(topologyattrs.Resids(np.arange(nres) + 1))
        attrs.append(topologyattrs.Atomids(np.arange(natoms) + 1))
        attrs.append(topologyattrs.Resnums(np.arange(nres) + 1))
        attrs.append(topologyattrs.Resnames([residue.name]))
        resid = np.zeros(natoms, dtype=np.int32)
        if isinstance(residue, ResidueType):
            bonds = [[residue.atom2index(ai), residue.atom2index(aj)] for ai, bondi in residue.connectivity.items() for aj in bondi]
        else:
            bonds = [[residue.name2index(residue.type.atom2name(ai)), residue.name2index(residue.type.atom2name(aj))] for ai, bondi in residue.type.connectivity.items() for aj in bondi]
        attrs.append(topologyattrs.Bonds(bonds))
        return Topology(natoms, nres, 1, attrs, resid, None)


class XpongeMoleculeReader(base.ReaderBase):
    """
    This **class** is used to read the Xponge Molecule to mdanalysis

    Create the following attributes:
        Atomnames
        Atomtypes
        Masses
        Charges
        Bonds

    Guesses the following attributes:
        Elements
    """
    def __init__(self, filename, **kwargs):
        filename.atoms = [atom for residue in filename.residues for atom in residue.atoms]
        filename.atom2index = {atom: i for i, atom in enumerate(filename.atoms)}
        super().__init__(filename, **kwargs)
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts.positions = [[getattr(atom, i) for i in "xyz"] for atom in self.filename.atoms]

    @property
    def n_atoms(self):
        return len(self.filename.atoms)

    @property
    def n_frames(self):
        return 1

    #pylint: disable=unused-argument
    def parse(self, **kwargs):
        """
        This **function** reads the file and returns the structure

        :param kwargs: keyword arguments
        :return: MDAnalysis Topology object
        """
        attrs = [topologyattrs.Segids(np.array(['SYSTEM'], dtype=object))]
        molecule = self.filename
        natoms = len(molecule.atoms)
        nres = len(molecule.residues)
        masses = [atom.mass for atom in molecule.atoms]
        elements = [guess_element_from_mass(mass) for mass in masses]
        attrs.append(topologyattrs.Masses(masses))
        attrs.append(topologyattrs.Elements(elements, guessed=True))
        attrs.append(topologyattrs.Atomnames([atom.name for atom in molecule.atoms]))
        attrs.append(topologyattrs.Atomtypes([atom.type.name for atom in molecule.atoms]))
        attrs.append(topologyattrs.Charges([atom.charge for atom in molecule.atoms]))
        attrs.append(topologyattrs.Resids(np.arange(nres) + 1))
        attrs.append(topologyattrs.Atomids(np.arange(natoms) + 1))
        attrs.append(topologyattrs.Resnums(np.arange(nres) + 1))
        attrs.append(topologyattrs.Resnames([residue.name for residue in molecule.residues]))
        resid = np.zeros(natoms, dtype=np.int32)
        count = 0
        for i, res in enumerate(molecule.residues):
            resid[count:count + len(res.atoms)] = i
            count += len(res.atoms)
        bonds = [[molecule.atom2index[self._t2a(residue, ai)], molecule.atom2index[self._t2a(residue, aj)]] for residue in molecule.residues for ai, bondi in residue.type.connectivity.items() for aj in bondi]
        bonds.extend([[molecule.atom2index[rl.atom1], molecule.atom2index[rl.atom2]] for rl in molecule.residue_links])
        
        attrs.append(topologyattrs.Bonds(bonds))
        return Topology(natoms, nres, 1, attrs, resid, None)

    @staticmethod
    def _t2a(residue, atom):
        return residue.name2atom(residue.type.atom2name(atom))

class SpongeTrajectoryReader(base.ReaderBase):
    """
    This **class** is the interface to MDAnalysis.

    :param dat_file_name: the name of the SPONGE dat trajectory file
    :param box: the name of the box file or a list of 6 ``int`` or ``float`` \
representing the 3 box lengths and 3 box angles.
    :param n_atoms: the number of atoms
    """
    def __init__(self, dat_file_name, n_atoms, **kwargs):
        super().__init__(dat_file_name, **kwargs)
        box = kwargs.get("box", None)
        if isinstance(box, str):
            self.boxname = box
            self.box = None
            self._get_box_offset()
        else:
            self.boxname = None
            self.box = box
        self._n_atoms = n_atoms
        self._n_frames = os.path.getsize(dat_file_name) // 12 // self.n_atoms
        self.trajfile = None
        self.boxfile = None
        self.ts = base.Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()

    @property
    def n_frames(self):
        """
        The total number of frames in the trajectory file
        """
        return self._n_frames

    @property
    def n_atoms(self):
        """
        The total number of atoms in the trajectory file
        """
        return self._n_atoms

    @classmethod
    def with_arguments(cls, **kwargs):
        """
        This **function** binds the arguments to the reader to initialize
        **New From 1.2.6.8**

        :param kwargs: the arguments
        :return: a subclass of SpongeTrajectoryReader
        """
        class SpongeTrajectoryReaderWithArguments(cls):
            def __init__(self, dat_file_name, n_atoms, **kwargs_):
                kwargs_.update(kwargs)
                super().__init__(dat_file_name, n_atoms, **kwargs_)
        return SpongeTrajectoryReaderWithArguments

    def close(self):
        """
        Close all the opened file

        :return: None
        """
        if self.trajfile is not None:
            self.trajfile.close()
            self.trajfile = None
        if self.boxfile is not None:
            self.boxfile.close()
            self.boxfile = None

    def open_trajectory(self):
        """
        Open the trajectory file

        :return: trajectory file and box file
        """
        self.trajfile = util.anyopen(self.filename, "rb")
        if self.box is None:
            self.boxfile = util.anyopen(self.boxname)
        ts = self.ts
        ts.frame = -1
        return self.trajfile, self.boxfile

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        """

        :param frame:
        :return:
        """
        if self.trajfile is None:
            self.open_trajectory()
        if self.boxfile is not None:
            self.boxfile.seek(self._offsets[frame])
        self.trajfile.seek(self.n_atoms * 12 * frame)
        self.ts.frame = frame - 1
        return self._read_next_timestep()

    def _read_next_timestep(self):
        """

        :return:
        """
        ts = self.ts
        if self.trajfile is None:
            self.open_trajectory()
        t = self.trajfile.read(12 * self.n_atoms)
        if not t:
            raise EOFError
        setattr(ts, "_pos", np.frombuffer(t, dtype=np.float32).reshape(self.n_atoms, 3))

        if self.box is not None:
            ts.dimensions = self.box
        else:
            ts.dimensions = list(map(float, self.boxfile.readline().split()))
        ts.frame += 1
        return ts

    def _get_box_offset(self):
        """

        :return:
        """
        self._offsets = [0]
        with util.openany(self.boxname) as f:
            line = f.readline()
            while line:
                self._offsets.append(f.tell())
                line = f.readline()
        self._offsets.pop()


class SpongeTrajectoryWriter():
    """
    This **class** is used to write the SPONGE trajectory (xxx.dat and xxx.box)

    usage example::

        import Xponge.analysis.md_analysis as xmda
        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PDB, XTC

        u = mda.Universe(PDB, XTC)

        with xmda.SpongeTrajectoryWriter("mda_test") as W:
            for ts in u.trajectory:
                W.write(u)

    :param prefix: the prefix of the output files
    :param write_box: whether to write the box file **New From 1.2.7.0**
    """
    def __init__(self, prefix, write_box=True):
        self.write_box = write_box
        self.datname = prefix + ".dat"
        self.boxname = prefix + ".box"
        self.datfile = None
        self.boxfile = None
        set_attribute_alternative_name(self, self.open)
        set_attribute_alternative_name(self, self.close)
        set_attribute_alternative_name(self, self.write)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        """
        This **function** opens the trajectory files

        :return: None
        """
        self.datfile = Xopen(self.datname, "wb")
        self.boxfile = Xopen(self.boxname, "w")

    def close(self):
        """
        This **function** closes the trajectory files

        :return: None
        """
        self.datfile.close()
        self.boxfile.close()

    def write(self, u):
        """
        This **function** writes the coordinates of the Universe to the output files

        :param u: an MDAnalysis.Universe instance
        :return: None
        """
        if isinstance(u, mda.Universe):
            ts = u.coord
        elif isinstance(u, mda.AtomGroup):
            ts = u.ts
        else:
            raise TypeError(f"u should be Universe or AtomGroup, but {type(u)} got")
        self.datfile.write(ts.positions.astype(np.float32).tobytes())
        if self.self.write_box:
            self.boxfile.write(" ".join([f"{i}" for i in ts.dimensions]) + "\n")


class SpongeCoordinateReader(base.ReaderBase):
    """
    This **class** is the interface to MDAnalysis.

    :param file_name: the name of the SPONGE coordinate trajectory file
    """
    def __init__(self, file_name, n_atoms, **kwargs):
        super().__init__(file_name, **kwargs)
        self._n_atoms = n_atoms
        self._n_frames = 1
        self.file = None
        self.start = 0
        self.ts = base.Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()

    @property
    def n_frames(self):
        """
        The total number of frames in the trajectory file
        """
        return self._n_frames

    @property
    def n_atoms(self):
        """
        The total number of atoms in the trajectory file
        """
        return self._n_atoms

    def close(self):
        """
        Close all the opened file

        :return: None
        """
        if self.file is not None:
            self.file.close()
            self.file = None

    def open_file(self):
        """
        Open the coordinate file

        :return: trajectory file and box file
        """
        self.file = util.anyopen(self.filename, "r")
        self.file.readline()
        self.start = self.file.tell()

    def _reopen(self):
        self.close()
        self.open_file()

    def _read_frame(self, frame):
        """

        :param frame:
        :return:
        """
        if self.file is None:
            self.open_file()
        self.file.seek(self.start)
        self.ts.frame = frame - 1
        return self._read_next_timestep()

    def _read_next_timestep(self):
        """

        :return:
        """
        ts = self.ts
        if self.file is None:
            self.open_file()
        if self.file.tell() != self.start:
            raise EOFError
        t = np.loadtxt(self.file, max_rows=self.n_atoms)
        setattr(ts, "_pos", t)
        box = list(map(float, self.file.readline().split()))
        ts.dimensions = box
        ts.frame += 1
        return ts


class SpongeCoordinateWriter():
    """
    This **class** is used to write the SPONGE coordinate file

    usage example::

        import Xponge.analysis.md_analysis as xmda
        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PDB, XTC

        u = mda.Universe(PDB, XTC)

        with xmda.SpongeCoordinateWriter("mda_test") as W:
            for ts in u.trajectory:
                W.write(u)

    :param file_name: the name of the output file
    :param n_atoms: the total number of atoms this Timestep describes
    """
    def __init__(self, file_name, n_atoms=None):
        self.filename = file_name
        self.file = None
        self.n_atoms = n_atoms
        set_attribute_alternative_name(self, self.open)
        set_attribute_alternative_name(self, self.close)
        set_attribute_alternative_name(self, self.write)

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        """
        This **function** opens the coordinate file

        :return: None
        """
        self.file = Xopen(self.filename, "w")

    def close(self):
        """
        This **function** closes the coordinate file

        :return: None
        """
        self.file.close()

    def write(self, u):
        """
        This **function** writes the coordinates of the Universe to the output files

        :param u: an MDAnalysis.Universe instance
        :return: None
        """
        if isinstance(u, mda.Universe):
            ts = u.coord
        elif isinstance(u, mda.AtomGroup):
            ts = u.ts
        else:
            raise TypeError(f"u should be Universe or AtomGroup, but {type(u)} got")
        if self.n_atoms is None:
            self.n_atoms = len(ts.positions)
        towrite = f"{self.n_atoms}\n"
        for crd in ts.positions[:self.n_atoms]:
            towrite += f"{crd[0]} {crd[1]} {crd[2]}\n"
        if ts.dimensions:
            towrite += " ".join([f"{i}" for i in ts.dimensions]) + "\n"
        else:
            towrite += "999 999 999 90 90 90\n"
        self.file.write(towrite)


set_global_alternative_names()
