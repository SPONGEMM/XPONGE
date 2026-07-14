from io import StringIO

import Xponge
import Xponge.forcefield.amber.ff14sb  # noqa: F401


PDB_TEXT = """\
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.450  20.000  30.000  1.00  0.00           C
ATOM      3  C   ALA A   1      12.000  21.420  30.000  1.00  0.00           C
ATOM      4  O   ALA A   1      11.300  22.360  30.000  1.00  0.00           O
END
"""


def test_explicit_periodic_box_translates_export_without_wrapping(tmp_path):
    molecule = Xponge.load_pdb(StringIO(PDB_TEXT))
    molecule.ignore_missing_atoms = True
    molecule.set_periodic_box(origin=[11.0, 21.0, 31.0], lengths=[20.0, 21.0, 22.0])

    Xponge.save_sponge_input(molecule, prefix="explicit", dirname=str(tmp_path))

    coordinate_lines = (tmp_path / "explicit_coordinate.txt").read_text().splitlines()
    assert coordinate_lines[1] == "-1.000000 -1.000000 -1.000000"
    assert coordinate_lines[-1] == "20.000000 21.000000 22.000000 90.000000 90.000000 90.000000"
    assert [molecule.atoms[0].x, molecule.atoms[0].y, molecule.atoms[0].z] == [10.0, 20.0, 30.0]
