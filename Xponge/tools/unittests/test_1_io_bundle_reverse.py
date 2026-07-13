"""Integration tests for bundled input to direct/legacy conversion."""

from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from Xponge.io_bundle import convert_bundle_to_legacy, convert_legacy_to_bundle
from Xponge.io_bundle.errors import (
    BundleConflictError,
    BundlePathError,
    BundleSchemaError,
    BundleValidationError,
)
from Xponge.io_bundle.state_parsers import parse_protocol_or_restart_file
from Xponge.io_bundle.topology_parsers import (
    listed_force_schemas,
    pairwise_force_schema,
    parse_listed_force_data_file,
    parse_pairwise_force_data_file,
    parse_topology_file,
)
from Xponge.io_bundle.trajectory_exporters import edges_to_box
from Xponge.io_bundle.trajectory_parsers import box_to_edges, parse_trajectory_file
from Xponge.tools.unittests.io_bundle_fixtures import write_basic_case


def _make_bundle(tmp_path: Path) -> Path:
    case_dir = tmp_path / "case"
    converted_dir = tmp_path / "converted"
    case_dir.mkdir()
    write_basic_case(case_dir)
    convert_legacy_to_bundle(case_dir, converted_dir)
    return converted_dir / "bundle"


def test_bundle_to_legacy_full_case(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"

    manifest = convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system")
    entries = {entry.contract_id: entry for entry in manifest.entries}

    assert entries["topology.mass"].status == "typed_exported"
    assert entries["topology.charge"].status == "typed_exported"
    assert entries["topology.residue"].status == "typed_exported"
    assert entries["topology.bond"].status == "typed_exported"
    assert entries["topology.angle"].status == "typed_exported"
    assert entries["topology.dihedral"].status == "typed_exported"
    assert entries["topology.improper_dihedral"].status == "typed_exported"
    assert entries["topology.LJ"].status == "typed_exported"
    assert entries["restart.coordinate"].status == "typed_exported"
    assert entries["restart.velocity"].status == "typed_exported"
    assert entries["trajectory.crd"].status == "typed_exported"
    assert entries["trajectory.box"].status == "typed_exported"
    assert entries["trajectory.vel"].status == "typed_exported"
    assert entries["topology.cmap"].status == "typed_exported"
    assert entries["topology.pairwise_force.data.custom_pair"].status == "typed_exported"
    assert entries["topology.listed_forces.data.custom_bond"].status == "typed_exported"

    mdin_text = (output_dir / "mdin.legacy.spg.toml").read_text(encoding="utf-8")
    assert "input_h5_topology_path" not in mdin_text
    assert "input_h5_protocol_path" not in mdin_text
    assert "input_h5_restart_path" not in mdin_text
    assert 'mass_in_file = "system_mass.txt"' in mdin_text
    assert 'coordinate_in_file = "system_coordinate.txt"' in mdin_text
    assert 'improper_dihedral_in_file = "system_improper_dihedral.txt"' in mdin_text
    assert 'custom_pair_in_file = "system_custom_pair.txt"' in mdin_text
    assert 'crd = "system_crd.dat"' in mdin_text

    mass = parse_topology_file("mass_in_file", output_dir / "system_mass.txt")
    mass_by_path = {dataset.path: dataset.data for dataset in mass}
    assert np.allclose(mass_by_path["/atoms/mass"], np.asarray([1.0, 16.0], dtype=np.float32))

    bond = parse_topology_file("bond_in_file", output_dir / "system_bond.txt")
    bond_by_path = {dataset.path: dataset.data for dataset in bond}
    assert np.array_equal(bond_by_path["/forcefield/bond/atoms"], np.asarray([[0, 1]], dtype=np.int32))
    assert np.allclose(bond_by_path["/forcefield/bond/k"], np.asarray([100.0], dtype=np.float32))

    rerun = parse_trajectory_file("crd", output_dir / "system_crd.dat", atom_count=2)
    rerun_by_path = {dataset.path: dataset.data for dataset in rerun}
    assert rerun_by_path["/particles/all/position/value"].shape == (2, 2, 3)
    assert (output_dir / "manifest.bundle_to_legacy.json").is_file()


def test_bundle_to_legacy_uses_typed_data_instead_of_stale_sidecar(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "topology.spgt.h5", "a") as handle:
        handle["/atoms/mass"][...] = np.asarray([2.0, 18.0], dtype=np.float32)

    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="updated")
    datasets = parse_topology_file("mass_in_file", output_dir / "updated_mass.txt")
    by_path = {dataset.path: dataset.data for dataset in datasets}
    assert np.allclose(by_path["/atoms/mass"], np.asarray([2.0, 18.0], dtype=np.float32))


def test_bundle_to_legacy_improper_uses_mutated_typed_pk_without_sidecar(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "topology.spgt.h5", "a") as handle:
        sidecar_keys = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in handle["/parameters/sponge/files/legacy_sidecars/key"][()]
        ]
        assert "improper_dihedral_in_file" not in sidecar_keys
        assert "improper_in_file" not in sidecar_keys
        handle["/forcefield/improper/pk"][...] = np.asarray([7.25], dtype=np.float32)

    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="updated")
    datasets = parse_topology_file(
        "improper_dihedral_in_file", output_dir / "updated_improper_dihedral.txt"
    )
    by_path = {dataset.path: dataset.data for dataset in datasets}
    assert np.allclose(by_path["/forcefield/improper/pk"], np.asarray([7.25], dtype=np.float32))


@pytest.mark.parametrize(
    ("key", "source_name", "exported_name"),
    (
        ("improper_dihedral_in_file", "improper.txt", "system_improper_dihedral.txt"),
        ("nb14_extra_in_file", "nb14_extra.txt", "system_nb14_extra.txt"),
        ("urey_bradley_in_file", "urey_bradley.txt", "system_urey_bradley.txt"),
        ("cmap_in_file", "cmap.txt", "system_cmap.txt"),
        ("gb_in_file", "gb.txt", "system_gb.txt"),
        ("virtual_atom_in_file", "virtual_atom.txt", "system_virtual_atom.txt"),
        ("LJ_soft_core_in_file", "lj_soft_core.txt", "system_LJ_soft_core.txt"),
        ("subsys_division_in_file", "subsys_division.txt", "system_subsys_division.txt"),
        ("EAM_atom_type_in_file", "eam_atom_type.txt", "system_EAM_atom_type.txt"),
        ("EAM_in_file", "eam.txt", "system_EAM.txt"),
        ("SW_in_file", "sw.txt", "system_SW.txt"),
        ("EDIP_in_file", "edip.txt", "system_EDIP.txt"),
        ("TERSOFF_in_file", "tersoff.txt", "system_TERSOFF.txt"),
        ("qc_type_in_file", "qc_type.txt", "system_qc_type.txt"),
        ("REAXFF_in_file", "reaxff.txt", "system_REAXFF.txt"),
        ("REAXFF_type_in_file", "reaxff_type.txt", "system_REAXFF_type.txt"),
        ("pairwise_force_in_file", "pairwise_force.txt", "system_pairwise_force.txt"),
        ("listed_forces_in_file", "listed_forces.txt", "system_listed_forces.txt"),
    ),
)
def test_topology_typed_exporters_roundtrip(tmp_path, key, source_name, exported_name):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system")

    original = parse_topology_file(key, tmp_path / "case" / source_name)
    restored = parse_topology_file(key, output_dir / exported_name)
    original_by_path = {dataset.path: dataset.data for dataset in original}
    restored_by_path = {dataset.path: dataset.data for dataset in restored}
    assert restored_by_path.keys() == original_by_path.keys()
    for path, expected in original_by_path.items():
        actual = restored_by_path[path]
        if np.asarray(expected).dtype.kind in {"f", "c"}:
            assert np.allclose(actual, expected, equal_nan=True), path
        else:
            assert np.array_equal(actual, expected), path


def test_custom_force_dynamic_exporters_roundtrip(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system")

    original_pair_schema = pairwise_force_schema(tmp_path / "case" / "pairwise_force.txt")
    restored_pair_schema = pairwise_force_schema(output_dir / "system_pairwise_force.txt")
    assert restored_pair_schema == original_pair_schema
    force_name, parameter_types, parameter_names, ij_count = restored_pair_schema
    original_pair = parse_pairwise_force_data_file(
        force_name,
        tmp_path / "case" / "custom_pair.txt",
        parameter_types=parameter_types,
        parameter_names=parameter_names,
        ij_parameter_count=ij_count,
    )
    restored_pair = parse_pairwise_force_data_file(
        force_name,
        output_dir / "system_custom_pair.txt",
        parameter_types=parameter_types,
        parameter_names=parameter_names,
        ij_parameter_count=ij_count,
    )
    _assert_dataset_lists_equal(original_pair, restored_pair)

    original_listed_schema = listed_force_schemas(tmp_path / "case" / "listed_forces.txt")
    restored_listed_schema = listed_force_schemas(output_dir / "system_listed_forces.txt")
    assert restored_listed_schema == original_listed_schema
    force_name, parameter_types, parameter_names = restored_listed_schema[0]
    original_listed = parse_listed_force_data_file(
        force_name,
        tmp_path / "case" / "custom_bond.txt",
        parameter_types=parameter_types,
        parameter_names=parameter_names,
    )
    restored_listed = parse_listed_force_data_file(
        force_name,
        output_dir / "system_custom_bond.txt",
        parameter_types=parameter_types,
        parameter_names=parameter_names,
    )
    _assert_dataset_lists_equal(original_listed, restored_listed)


@pytest.mark.parametrize(
    ("key", "source_name", "exported_name", "ignored_paths"),
    (
        ("cv_in_file", "cv.txt", "system_cv.txt", ()),
        ("restrain_in_file", "restrain.txt", "system_restrain.txt", ()),
        ("restrain_cv_in_file", "restrain_cv.txt", "system_restrain_cv.txt", ()),
        ("steer_cv_in_file", "steer_cv.txt", "system_steer_cv.txt", ()),
        ("SITS_in_file", "sits.txt", "system_SITS.txt", ()),
        ("SITS_nk_in_file", "sits_nk.txt", "system_SITS_nk.txt", ()),
        ("restrain_atom_id", "restrain_atom_id.txt", "system_restrain_atom_id.txt", ()),
        ("restrain_weight_in_file", "restrain_weight.txt", "system_restrain_weight.txt", ()),
        (
            "restrain_coordinate_in_file",
            "restrain_coordinate.txt",
            "system_restrain_coordinate.txt",
            (),
        ),
        (
            "nose_hoover_chain_restart_input",
            "nhc_restart.txt",
            "system_nose_hoover_chain_restart_input.txt",
            (),
        ),
        ("constrain_in_file", "constrain.txt", "system_constrain.txt", ()),
        ("SITS_atom_in_file", "sits_atom.txt", "system_SITS_atom.txt", ()),
        ("soft_walls_in_file", "soft_walls.txt", "system_soft_walls.txt", ()),
        ("meta_edge_in_file", "meta_edge.txt", "system_meta_edge.txt", ()),
        (
            "meta_potential_in_file",
            "meta_potential.txt",
            "system_meta_potential.txt",
            ("/parameters/restart/bias/meta/default/potential_export",),
        ),
        ("meta_scatter_in_file", "meta_scatter.txt", "system_meta_scatter.txt", ()),
        (
            "hills_in_file",
            "hills.txt",
            "system_hills.txt",
            ("/parameters/restart/bias/meta/default/hills",),
        ),
    ),
)
def test_state_typed_exporters_roundtrip(tmp_path, key, source_name, exported_name, ignored_paths):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system")

    original = parse_protocol_or_restart_file(key, tmp_path / "case" / source_name)
    restored = parse_protocol_or_restart_file(key, output_dir / exported_name)
    _assert_dataset_lists_equal(original, restored, ignored_paths=set(ignored_paths))


def test_bundle_to_legacy_dry_run_does_not_create_output(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"

    manifest = convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system", dry_run=True)

    assert manifest.generated_mdin == str(output_dir / "mdin.legacy.spg.toml")
    assert not output_dir.exists()


def test_full_legacy_bundle_legacy_bundle_semantic_roundtrip(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_a = _make_bundle(tmp_path)
    legacy_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_a, legacy_dir, prefix="roundtrip")
    converted_c = tmp_path / "converted_c"
    convert_legacy_to_bundle(
        legacy_dir,
        converted_c,
        mdin="mdin.legacy.spg.toml",
    )
    bundle_c = converted_c / "bundle"

    paths_by_file = {
        "topology.spgt.h5": (
            "/atoms/mass",
            "/atoms/charge",
            "/residues/atom_offset",
            "/forcefield/bond/atoms",
            "/forcefield/improper/pk",
            "/forcefield/improper/phi0",
            "/forcefield/cmap/grid_value",
            "/forcefield/lj_soft_core/pair_AA",
            "/manybody/eam/embed/raw_ev",
            "/manybody/sw/triple/parameters",
            "/manybody/tersoff/entry/parameters_raw",
            "/manybody/reaxff/parameters/general/value",
            "/forcefield/custom_force/pairwise/data/custom_pair/parameter/value",
            "/forcefield/custom_force/listed/data/custom_bond/parameter/value",
        ),
        "protocol.spgp.h5": (
            "/cv/config/section/name",
            "/restraint/default/weight",
            "/meta/default/grid/coordinate",
        ),
        "restart.spgr.h5": (
            "/particles/all/position/value",
            "/particles/all/velocity/value",
            "/parameters/restart/bias/sits/SITS/nk",
            "/parameters/restart/bias/meta/default/potential/value",
        ),
        "trajectory.spg.h5md": (
            "/particles/all/position/value",
            "/particles/all/velocity/value",
            "/particles/all/box/edges/value",
        ),
    }
    for filename, paths in paths_by_file.items():
        with h5py.File(bundle_a / filename, "r") as left, h5py.File(bundle_c / filename, "r") as right:
            for path in paths:
                expected = np.asarray(left[path][...])
                actual = np.asarray(right[path][...])
                assert actual.shape == expected.shape, (filename, path)
                if expected.dtype.kind in {"f", "c"}:
                    assert np.allclose(actual, expected, equal_nan=True), (filename, path)
                else:
                    assert np.array_equal(actual, expected), (filename, path)


def test_bundle_to_legacy_refuses_existing_targets(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy"
    output_dir.mkdir()
    (output_dir / "system_mass.txt").write_text("user data\n", encoding="utf-8")

    with pytest.raises(BundleConflictError, match="already exist"):
        convert_bundle_to_legacy(bundle_dir, output_dir, prefix="system")
    assert (output_dir / "system_mass.txt").read_text(encoding="utf-8") == "user data\n"


def test_bundle_to_legacy_rejects_schema_mismatch(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "topology.spgt.h5", "a") as handle:
        del handle["/schema/version"]
        handle.create_dataset("/schema/version", data="future.v2", dtype=h5py.string_dtype("utf-8"))

    with pytest.raises(BundleSchemaError, match="schema version"):
        convert_bundle_to_legacy(bundle_dir, tmp_path / "legacy")


def test_bundle_to_legacy_rejects_sidecar_path_traversal(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "topology.spgt.h5", "a") as handle:
        paths = handle["/parameters/sponge/files/legacy_sidecars/path"]
        paths[0] = "../../outside.txt"

    with pytest.raises(BundlePathError, match="escapes bundle root"):
        convert_bundle_to_legacy(bundle_dir, tmp_path / "legacy")


def test_bundle_to_legacy_rejects_missing_declared_sidecar(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "topology.spgt.h5", "r") as handle:
        raw_path = handle["/parameters/sponge/files/legacy_sidecars/path"][0]
        sidecar = raw_path.decode("utf-8") if isinstance(raw_path, bytes) else str(raw_path)
    (bundle_dir / sidecar).unlink()

    with pytest.raises(BundleValidationError, match="legacy sidecar .* does not exist"):
        convert_bundle_to_legacy(bundle_dir, tmp_path / "legacy")


def test_bundle_to_legacy_uses_compatibility_payload_without_sidecar(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    topology = bundle_dir / "topology.spgt.h5"
    with h5py.File(topology, "a") as handle:
        atom_count = int(handle["/topology/atom_count"][()])
        sidecar_key_path = "/parameters/sponge/files/legacy_sidecars/key"
        sidecar_value_path = "/parameters/sponge/files/legacy_sidecars/path"
        sidecar_keys = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in handle[sidecar_key_path][()]
        ]
        sidecar_values = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in handle[sidecar_value_path][()]
        ]
        retained_sidecars = [
            (key, value)
            for key, value in zip(sidecar_keys, sidecar_values)
            if key != "mass_in_file"
        ]
        del handle[sidecar_key_path]
        del handle[sidecar_value_path]
        dtype = h5py.string_dtype("utf-8")
        handle.create_dataset(
            sidecar_key_path,
            data=[key for key, _ in retained_sidecars],
            dtype=dtype,
        )
        handle.create_dataset(
            sidecar_value_path,
            data=[value for _, value in retained_sidecars],
            dtype=dtype,
        )
        del handle["/atoms/mass"]
        handle.create_dataset(
            "/compatibility/legacy_import/mass_in_file/raw_text",
            data=f"{atom_count}\n" + "\n".join("1.0" for _ in range(atom_count)) + "\n",
            dtype=dtype,
        )

    output_dir = tmp_path / "compatibility_legacy"
    manifest = convert_bundle_to_legacy(bundle_dir, output_dir, prefix="compat")
    entries = {entry.contract_id: entry for entry in manifest.entries}
    assert entries["topology.mass"].status == "compatibility_restored"
    datasets = parse_topology_file("mass_in_file", output_dir / "compat_mass.txt")
    assert np.allclose(datasets[0].data, 1.0)


def test_bundle_to_legacy_validates_restart_load_policy(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "restart.spgr.h5", "a") as handle:
        del handle["/parameters/restart/thermostat/nose_hoover_chain"]
    mdin_path = bundle_dir / "mdin.bundled.spg.toml"
    mdin_path.write_text(
        mdin_path.read_text(encoding="utf-8").replace(
            'input_h5_restart_load = "full"',
            'input_h5_restart_load = "dynamic"',
        ),
        encoding="utf-8",
    )

    with pytest.raises(BundleValidationError, match="no corresponding restart state"):
        convert_bundle_to_legacy(bundle_dir, tmp_path / "legacy")

    manifest = convert_bundle_to_legacy(
        bundle_dir,
        tmp_path / "legacy_non_strict",
        strict=False,
        dry_run=True,
    )
    assert any("no corresponding restart state" in warning for warning in manifest.warnings)


def test_bundle_to_legacy_selects_configured_particle_stream(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    with h5py.File(bundle_dir / "trajectory.spg.h5md", "a") as handle:
        handle.move("/particles/all", "/particles/custom")
    mdin_path = bundle_dir / "mdin.bundled.spg.toml"
    mdin_path.write_text(
        mdin_path.read_text(encoding="utf-8").replace(
            'input_h5_trajectory_particle_stream = "all"',
            'input_h5_trajectory_particle_stream = "custom"',
        ),
        encoding="utf-8",
    )

    output_dir = tmp_path / "legacy"
    convert_bundle_to_legacy(bundle_dir, output_dir, prefix="custom")
    datasets = parse_trajectory_file("crd", output_dir / "custom_crd.dat", atom_count=2)
    assert datasets[0].data.shape == (2, 2, 3)


def test_bundle_to_legacy_rejects_mismatched_trajectory_frames(tmp_path):
    h5py = pytest.importorskip("h5py")
    bundle_dir = _make_bundle(tmp_path)
    trajectory = bundle_dir / "trajectory.spg.h5md"
    with h5py.File(trajectory, "a") as handle:
        values = handle["/particles/all/velocity/value"][...][:1]
        del handle["/particles/all/velocity/value"]
        handle.create_dataset("/particles/all/velocity/value", data=values)

    with pytest.raises(BundleValidationError, match="mismatched frame counts"):
        convert_bundle_to_legacy(bundle_dir, tmp_path / "legacy")


def test_edges_to_box_roundtrip_for_triclinic_cell():
    box = np.asarray([10.0, 11.0, 12.0, 75.0, 80.0, 70.0], dtype=np.float32)
    edges = box_to_edges(box)
    restored = edges_to_box(edges)
    assert np.allclose(restored, box, atol=1e-5)


def test_bundle_to_legacy_cli(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy_cli"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "Xponge",
            "bundle-to-legacy",
            str(bundle_dir),
            "-o",
            str(output_dir),
            "--prefix",
            "cli",
        ],
        cwd=Path(__file__).resolve().parents[3],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert "typed exported entries:" in result.stdout
    assert (output_dir / "cli_mass.txt").is_file()
    assert (output_dir / "mdin.legacy.spg.toml").is_file()


def test_bundle_to_legacy_cli_dry_run(tmp_path):
    bundle_dir = _make_bundle(tmp_path)
    output_dir = tmp_path / "legacy_cli_dry"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "Xponge",
            "bundle-to-legacy",
            str(bundle_dir),
            "-o",
            str(output_dir),
            "--dry-run",
        ],
        cwd=Path(__file__).resolve().parents[3],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (result.stdout, result.stderr)
    assert "generated mdin:" in result.stdout
    assert not output_dir.exists()


def _assert_dataset_lists_equal(original, restored, *, ignored_paths=None):
    ignored_paths = ignored_paths or set()
    original_by_path = {dataset.path: dataset.data for dataset in original}
    restored_by_path = {dataset.path: dataset.data for dataset in restored}
    for path in ignored_paths:
        original_by_path.pop(path, None)
        restored_by_path.pop(path, None)
    assert restored_by_path.keys() == original_by_path.keys()
    for path, expected in original_by_path.items():
        actual = restored_by_path[path]
        if np.asarray(expected).dtype.kind in {"f", "c"}:
            assert np.allclose(actual, expected, equal_nan=True), path
        else:
            assert np.array_equal(actual, expected), path
