"""
This **module** implements the SASA calculation
"""
from functools import partial
import numpy as np
from .. import Xprint
from ..helper.math import get_fibonacci_grid
from . import md_analysis as xmda
mda = xmda.mda
class SASA:
    """
        This **class** implements the SASA calculation

        :param u: the instance of the MDAnalysis universe
        :param r_probe: radius of the probe
        :param n_points:  resolution of the surface of each atom
        :param radii_dict: dict of atomic radii to update the default RADII dict
        :param r_atoms: if specify, the atomic radii will directly use the value here \
instead of looking up in radii_dict.
    """
    RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52,
             "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75,
             "Br": 1.85, "I": 1.98}
    def __init__(self, u, r_probe=1.4, n_points=100, radii_dict=None, r_atoms=None):
        real_radii_dict = {}
        real_radii_dict.update(self.RADII)
        if radii_dict is not None:
            if not isinstance(radii_dict, dict):
                raise TypeError(f"radii_dict should be a dict, but {type(radii_dict)} got")
            real_radii_dict.update(radii_dict)
        if r_atoms is not None:
            if not isinstance(r_atoms, (int, float)):
                raise TypeError(f"r_atoms should be a number, but {type(r_atoms)} got")
            self.r_atoms = [r_atoms for _ in u.atoms]
        else:
            if hasattr(u.atoms, "elements"):
                self.r_atoms = [real_radii_dict.get(atom.element, 2) for atom in u.atoms]
            elif hasattr(u.atoms, "names"):
                guess = mda.topology.guessers.guess_atom_element
                self.r_atoms = [real_radii_dict.get(guess(atom.name), 2) for atom in u.atoms]
            else:
                raise ValueError("No name or element information in the universe")
        self.r_atoms = np.array(self.r_atoms)
        self.r = self.r_atoms + r_probe
        self.r_max = np.max(self.r) * 2
        self._bad_points = get_fibonacci_grid(n_points, 0, 1)
        self.u = u
        self.atom_map = {atom.index: i for i, atom in enumerate(u.atoms)}
        self.surface_area = None
        self.surface = None

    def main(self, need_area=True, need_surface=True):
        """
            This **function** does the real calculation.

            :param need_area: whether keep the area. If True, the surface will stored in self.surface_area
            :param need_surface: whether need the surface. If True, the surface will stored in self.surface
        """
        self.surface_area = None
        self.surface = None
        if not need_area and not need_surface:
            Xprint("Both need_area and need_surface are False, the calculation stops", "WARNING")
            return
        u = self.u
        n_points = len(self._bad_points)
        n_atoms = len(self.u.atoms)

        if need_area:
            area = np.zeros(n_atoms, dtype=np.int64)
        if need_surface:
            surface = []
        alls = set(range(n_points))
        if hasattr(u, "trajectory") and hasattr(u.trajectory.ts, "dimensions"):
            search_neighbor = partial(mda.lib.distances.capped_distance, box=u.trajectory.ts.dimensions)
        elif hasattr(u, "ts") and hasattr(u.ts, "dimensions"):
            search_neighbor = partial(mda.lib.distances.capped_distance, box=u.ts.dimensions)
        else:
            search_neighbor = mda.lib.distances.capped_distance
        for i, r_i in enumerate(self.r):
            real_points = np.array(self._bad_points, copy=True) * r_i + u.atoms[i].position
            real_alls = alls.copy()
            pairs, distances = search_neighbor(u.atoms.positions, u.atoms[i].position, self.r_max)
            for j, d in zip(pairs[:,0], distances):
                if i == j:
                    continue
                r_j = self.r[j]
                if d < (r_i + r_j):
                    bad_pairs = search_neighbor(real_points, u.atoms[j].position, r_j, return_distances=False)
                    real_alls -= set(bad_pairs[:,0])
            if need_area:
                area[i] = len(real_alls)
            if need_surface:
                surface.append(real_points[list(real_alls)])

        if need_surface:
            self.surface = np.vstack(surface)
        if need_area:
            self.surface_area = 4 * np.pi * self.r * self.r * area / n_points

    def write_surface_xyz(self, filename, headlines=True):
        """
            This **function** writes the surface information to a file

            :param filename: the name of the output file
            :param headlines: whether the headlines and fake elements are written to the file
        """
        if self.surface is None:
            raise ValueError("You should call main(need_surface=True) to get the surface information first")
        if headlines:
            n = len(self.surface)
            header = f"{n}\n{n}"
            fmt = 'H %3.8e %3.8e %3.8e'
        else:
            header = ""
            fmt = '%3.8e %3.8e %3.8e'
        np.savetxt(filename, self.surface, header=header, fmt=fmt, comments='')

    def get_sasa_result(self, select="all"):
        """
            This **function** gets the sasa result

            :param select: the sasa  sum of the selected atoms
        """
        if self.surface_area is None:
            raise ValueError("You should call main(need_area=True) to get the surface information first")
        atoms = self.u.select_atoms(select)
        return np.sum([self.surface_area[self.atom_map[atom.index]] for atom in atoms])

Xprint("""Reference for sasa:
  Shrake, A; Rupley, JA.
    Environment and exposure to solvent of protein atoms. Lysozyme and insulin
    Journal of Molecular Biology 1973 79 (2) 351-364
    DOI: 10.1016/0022-2836(73)90011-9""")
