import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
def _run_python(code, *args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        [sys.executable, *args],
        input=code,
        text=True,
        capture_output=True,
        cwd=ROOT,
        env=env,
        check=True,
    )


def test_cli_version_and_plain_import_have_no_reference_side_effect():
    version = _run_python("", "-m", "Xponge", "-v")
    assert version.stdout.strip() == "1.6b4"
    assert version.stderr == ""

    imported = _run_python("import Xponge\n")
    assert imported.stdout == ""
    assert imported.stderr == ""


def test_cli_help_and_subcommand_help_remain_available():
    top_level = _run_python("", "-m", "Xponge", "--help")
    subcommand = _run_python("", "-m", "Xponge", "test", "--help")
    assert "show the version of Xponge" in top_level.stdout
    assert "forcefield_loading" not in top_level.stderr
    assert "--purpose" in subcommand.stdout


def test_resp_references_remain_explicitly_available():
    result = _run_python(
        "from Xponge.assign import resp\n"
        "assert 'DOI: 10.1021/j100142a004' in resp.RESP_REFERENCE_TEXT\n"
        "resp.print_references()\n"
    )
    assert (result.stdout + result.stderr).count("Reference for resp.py:") == 1


def test_frcmod_scee_and_scnb_are_not_swapped(tmp_path):
    from Xponge.load import load_frcmod

    frcmod = tmp_path / "nb14.frcmod"
    frcmod.write_text(
        """NB14 scale regression
DIHE
qA-qB-qC-qD  1  0.2  0.0  1.0  SCEE=1.2 SCNB=2.0
qE-qF-qG-qH  1  0.2  0.0  1.0  SCEE=1.0 SCNB=1.0
"""
    )
    parsed = load_frcmod(frcmod, include_nb14=True)
    nb14 = parsed[6]
    assert "qA-qB-qC-qD 0.5 0.8333333333333334" in nb14
    assert "qE-qF-qG-qH 1.0 1.0" in nb14
