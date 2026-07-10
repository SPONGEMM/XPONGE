import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _run(code):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
    )


@pytest.mark.parametrize(
    ("first", "second", "family", "active", "requested"),
    [
        ("ff14sb", "ff19sb", "protein", "ff14sb", "ff19sb"),
        ("ff19sb", "ff14sb", "protein", "ff19sb", "ff14sb"),
        ("gaff", "gaff2", "small_molecule", "gaff", "gaff2"),
        ("gaff2", "gaff", "small_molecule", "gaff2", "gaff"),
    ],
)
def test_incompatible_amber_forcefield_families_fail_before_replacement(
    first, second, family, active, requested
):
    result = _run(
        "import importlib\n"
        "from Xponge.forcefield.amber._forcefield_family import (\n"
        "    ForceFieldFamilyConflictError, get_active_forcefield\n"
        ")\n"
        f"importlib.import_module('Xponge.forcefield.amber.{first}')\n"
        "try:\n"
        f"    importlib.import_module('Xponge.forcefield.amber.{second}')\n"
        "except ForceFieldFamilyConflictError as exc:\n"
        f"    assert {requested!r} in str(exc)\n"
        "else:\n"
        "    raise AssertionError('expected a force-field family conflict')\n"
        f"assert get_active_forcefield({family!r}) == {active!r}\n"
    )
    assert result.returncode == 0, result.stderr


def test_cross_family_amber_forcefields_can_be_combined():
    result = _run(
        "import Xponge.forcefield.amber.ff14sb\n"
        "import Xponge.forcefield.amber.gaff2\n"
        "import Xponge.forcefield.amber.lipid17\n"
        "from Xponge.forcefield.amber._forcefield_family import get_active_forcefield\n"
        "assert get_active_forcefield('protein') == 'ff14sb'\n"
        "assert get_active_forcefield('small_molecule') == 'gaff2'\n"
        "assert get_active_forcefield('lipid') == 'lipid17'\n"
    )
    assert result.returncode == 0, result.stderr


def test_family_activation_is_idempotent_and_dependencies_are_explicit():
    from Xponge.forcefield.amber._forcefield_family import (
        activate_forcefield_family,
        require_forcefield_family,
    )

    assert activate_forcefield_family("test-only", "one") == "one"
    assert activate_forcefield_family("test-only", "one") == "one"
    assert require_forcefield_family("test-only", {"one", "two"}) == "one"

