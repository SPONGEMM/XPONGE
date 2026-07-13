"""Tests for direct XPONGE object to bundled input saving."""

from pathlib import Path
import subprocess

import numpy as np
import pytest

import Xponge
from Xponge.io_bundle import convert_bundle_to_legacy
from Xponge.io_bundle.contracts import classify_sponge_serializer_key
from Xponge.io_bundle.errors import BundleCapabilityError, BundlePathError
from Xponge.io_bundle.topology_parsers import parse_topology_file


def _alanine_dipeptide():
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.get_peptide_from_sequence("AA")


def _alanine_residue_type():
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.ResidueType.get_type("ALA")


def _write_saver_mdin(root: Path, prefix: str, *, mode: str = "nvt") -> None:
    (root / "mdin.bundled.spg.toml").write_text(
        f'mode = "{mode}"\n'
        f'input_h5_topology_path = "{prefix}_topology.spgt.h5"\n'
        f'input_h5_protocol_path = "{prefix}_protocol.spgp.h5"\n'
        f'input_h5_restart_path = "{prefix}_restart.spgr.h5"\n'
        'input_h5_restart_load = "structural"\n',
        encoding="utf-8",
    )


def test_save_sponge_input_bundle_writes_loadable_artifacts(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()

    returned = Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)

    assert returned is molecule
    topology = tmp_path / "system_topology.spgt.h5"
    protocol = tmp_path / "system_protocol.spgp.h5"
    restart = tmp_path / "system_restart.spgr.h5"
    assert topology.is_file()
    assert protocol.is_file()
    assert restart.is_file()
    with h5py.File(topology, "r") as handle:
        assert handle["/atoms/mass"].shape == (len(molecule.atoms),)
        assert handle["/atoms/charge"].shape == (len(molecule.atoms),)
        assert handle["/parameters/xponge/atoms/name"].shape == (len(molecule.atoms),)
        assert handle["/parameters/xponge/residues/name"].shape == (len(molecule.residues),)
        assert handle["/topology/atom_count"][()] == len(molecule.atoms)
    with h5py.File(protocol, "r") as handle:
        assert handle["/schema/name"][()].decode() == "sponge.protocol.h5"
    with h5py.File(restart, "r") as handle:
        assert handle["/particles/all/position/value"].shape == (1, len(molecule.atoms), 3)
        assert handle["/particles/all/box/edges/value"].shape == (1, 3, 3)


def test_bundled_saver_matches_legacy_saver_semantics(tmp_path):
    molecule = _alanine_dipeptide()
    legacy_dir = tmp_path / "legacy"
    bundle_dir = tmp_path / "bundled"
    restored_dir = tmp_path / "restored"
    legacy_dir.mkdir()
    bundle_dir.mkdir()

    Xponge.save_sponge_input(molecule, "system", legacy_dir)
    molecule.box_length = None
    Xponge.save_sponge_input_bundle(molecule, "system", bundle_dir)
    _write_saver_mdin(bundle_dir, "system")
    convert_bundle_to_legacy(bundle_dir, restored_dir, prefix="system")

    original_files = {path.name: path for path in legacy_dir.glob("system_*.txt")}
    restored_files = {path.name: path for path in restored_dir.glob("system_*.txt")}
    assert restored_files.keys() == original_files.keys()
    for filename, original_path in original_files.items():
        serializer_key = filename.removeprefix("system_").removesuffix(".txt")
        classification, legacy_key = classify_sponge_serializer_key(serializer_key)
        restored_path = restored_files[filename]
        if classification == "metadata":
            assert restored_path.read_text(encoding="utf-8") == original_path.read_text(encoding="utf-8")
        elif serializer_key == "coordinate":
            original_coordinate, original_box = Xponge.load_coordinate(str(original_path))
            restored_coordinate, restored_box = Xponge.load_coordinate(str(restored_path))
            assert np.allclose(restored_coordinate, original_coordinate)
            assert np.allclose(restored_box, original_box)
        else:
            original = parse_topology_file(legacy_key, original_path)
            restored = parse_topology_file(legacy_key, restored_path)
            _assert_dataset_lists_equal(original, restored)


@pytest.mark.parametrize("source_kind", ("residue", "residue_type"))
def test_save_sponge_input_bundle_accepts_residue_inputs(tmp_path, source_kind):
    residue_type = _alanine_residue_type()
    if source_kind == "residue_type":
        source = residue_type
    else:
        source = Xponge.Residue(residue_type)
        for atom in residue_type.atoms:
            source.Add_Atom(atom)

    molecule = Xponge.save_sponge_input_bundle(source, source_kind, tmp_path)

    assert isinstance(molecule, Xponge.Molecule)
    assert (tmp_path / f"{source_kind}_topology.spgt.h5").is_file()
    assert (tmp_path / f"{source_kind}_restart.spgr.h5").is_file()


def test_save_sponge_input_bundle_rejects_unknown_nonempty_serializer(tmp_path):
    molecule = _alanine_dipeptide()

    @Xponge.Molecule.Set_Save_SPONGE_Input("future_force")
    def write_future_force(_):
        return "1\n"

    try:
        with pytest.raises(BundleCapabilityError, match="unclassified serializer 'future_force'"):
            Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    finally:
        Xponge.Molecule.Del_Save_SPONGE_Input("future_force")


def test_save_sponge_input_bundle_rejects_escaping_prefix(tmp_path):
    with pytest.raises(BundlePathError, match="escapes output directory"):
        Xponge.save_sponge_input_bundle(_alanine_dipeptide(), "../escape", tmp_path)


def test_saver_bundle_without_velocity_starts_optional_sponge(tmp_path):
    sponge = Path("/home/youmans/sidereus/SPONGE/build-dev-cuda13/SPONGE")
    if not (sponge.is_file() and sponge.stat().st_mode & 0o111):
        pytest.skip("SPONGE executable is unavailable")

    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    _write_saver_mdin(tmp_path, "system", mode="minimization")

    result = subprocess.run(
        [
            str(sponge),
            "-mdin",
            str(tmp_path / "mdin.bundled.spg.toml"),
            "-step_limit",
            "0",
        ],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert "All velocity will be set to 0" in result.stdout


def test_reversed_saver_bundle_starts_optional_sponge(tmp_path):
    sponge = Path("/home/youmans/sidereus/SPONGE/build-dev-cuda13/SPONGE")
    if not (sponge.is_file() and sponge.stat().st_mode & 0o111):
        pytest.skip("SPONGE executable is unavailable")

    bundle_dir = tmp_path / "bundle"
    legacy_dir = tmp_path / "legacy"
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    Xponge.save_sponge_input_bundle(molecule, "system", bundle_dir)
    _write_saver_mdin(bundle_dir, "system", mode="minimization")
    manifest = convert_bundle_to_legacy(bundle_dir, legacy_dir, prefix="system")

    result = subprocess.run(
        [
            str(sponge),
            "-mdin",
            str(manifest.generated_mdin),
            "-step_limit",
            "0",
        ],
        cwd=legacy_dir,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)


def test_save_sponge_input_bundle_restart_passes_optional_sponge_probe(tmp_path):
    probe = Path("/home/youmans/sidereus/SPONGE/build-dev-cuda13/tests/h5_bundle/h5_restart_read_probe")
    if not (probe.is_file() and probe.stat().st_mode & 0o111):
        pytest.skip("SPONGE restart H5 probe is unavailable")
    molecule = _alanine_dipeptide()
    Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    restart_path = tmp_path / "system_restart.spgr.h5"
    result = subprocess.run(
        [str(probe), str(restart_path), str(len(molecule.atoms))],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)


def _assert_dataset_lists_equal(original, restored):
    original_by_path = {dataset.path: dataset.data for dataset in original}
    restored_by_path = {dataset.path: dataset.data for dataset in restored}
    assert restored_by_path.keys() == original_by_path.keys()
    for path, expected in original_by_path.items():
        actual = restored_by_path[path]
        if np.asarray(expected).dtype.kind in {"f", "c"}:
            assert np.allclose(actual, expected, equal_nan=True), path
        else:
            assert np.array_equal(actual, expected), path
