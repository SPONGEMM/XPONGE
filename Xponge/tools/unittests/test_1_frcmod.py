"""Regression tests for Amber frcmod parsing."""

import pytest

from Xponge.load import load_frcmod
from Xponge.forcefield import amber
from Xponge.forcefield.base import nb14_base


def test_frcmod_dihedral_continuation_and_nb14(tmp_path):
    frcmod = tmp_path / "continuation.frcmod"
    frcmod.write_text(
        """continuation test
DIHE
A -B -C -D    1  0.6  0.0  -3.0  SCEE=1.0 SCNB=1.0
              1  0.2  180.0  1.0  SCEE=1.0 SCNB=1.0

NONB
A  1.5  0.1
""",
        encoding="utf-8",
    )

    parsed = load_frcmod(frcmod, include_nb14=True)
    propers = parsed[3].splitlines()
    nb14s = parsed[6].splitlines()

    assert propers[1].split() == ["A-B-C-D", "0.6", "0.0", "3", "1"]
    assert propers[2].split() == ["A-B-C-D", "0.2", "180.0", "1", "0"]
    assert nb14s[1].split() == ["A-D", "1.0", "1.0"]
    assert nb14s[2].split() == ["A-D", "1.0", "1.0"]
    assert len(load_frcmod(frcmod)) == 7


def test_frcmod_registration_includes_nb14_scaling(tmp_path):
    frcmod = tmp_path / "nb14.frcmod"
    frcmod.write_text(
        """nb14 registration test
DIHE
qA-qB-qC-qD    1  0.6  0.0  3.0  SCEE=1.0 SCNB=1.0
""",
        encoding="utf-8",
    )

    amber.load_parameters_from_frcmod(frcmod, prefix=False)

    nb14 = nb14_base.NB14Type.get_type("qA-qD")
    assert nb14.kLJ == pytest.approx(1.0)
    assert nb14.kee == pytest.approx(1.0)


def test_frcmod_rejects_orphan_dihedral_continuation(tmp_path):
    frcmod = tmp_path / "orphan.frcmod"
    frcmod.write_text(
        """orphan continuation test
DIHE
              1  0.2  180.0  1.0
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="no preceding atom types"):
        load_frcmod(frcmod)
