"""Tests for direct XPONGE object to bundled input saving."""

import os
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pytest

import Xponge
from Xponge.io_bundle import convert_bundle_to_legacy
from Xponge.io_bundle.contracts import classify_sponge_serializer_key
from Xponge.io_bundle.errors import (
    BundleCapabilityError,
    BundlePathError,
    BundleValidationError,
)
from Xponge.io_bundle.topology_parsers import parse_topology_file


def _optional_executable(env_var: str, command: str | None = None) -> Path:
    configured = os.environ.get(env_var)
    discovered = configured or (shutil.which(command) if command else None)
    if not discovered:
        pytest.skip(f"set {env_var} to run this optional SPONGE integration test")
    executable = Path(discovered).expanduser()
    if not (executable.is_file() and executable.stat().st_mode & 0o111):
        pytest.skip(f"{env_var} does not point to an executable file: {executable}")
    return executable


def _alanine_dipeptide():
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.get_peptide_from_sequence("AA")


def _alanine_residue_type():
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    return Xponge.ResidueType.get_type("ALA")


def _write_saver_mdin(
    root: Path,
    prefix: str,
    *,
    mode: str = "nvt",
    restart_load: str = "structural",
) -> None:
    (root / "mdin.bundled.spg.toml").write_text(
        f'mode = "{mode}"\n'
        f'input_h5_topology_path = "{prefix}_topology.spgt.h5"\n'
        f'input_h5_protocol_path = "{prefix}_protocol.spgp.h5"\n'
        f'input_h5_restart_path = "{prefix}_restart.spgr.h5"\n'
        f'input_h5_restart_load = "{restart_load}"\n',
        encoding="utf-8",
    )


def _add_native_custom_force_fixture(topology_path: Path) -> None:
    h5py = pytest.importorskip("h5py")
    with h5py.File(topology_path, "a") as handle:
        string_dtype = h5py.string_dtype(encoding="utf-8")
        pair = handle.require_group("/forcefield/custom_force/pairwise")
        pair.create_dataset("name", data="custom_pair", dtype=string_dtype)
        pair.create_dataset(
            "potential",
            data="E = epsilon_ij * powf(sigma_ij / r_ij, 12.0f);",
            dtype=string_dtype,
        )
        parameters = pair.require_group("parameters")
        parameters.create_dataset(
            "type", data=np.asarray(["float", "float"], dtype=object), dtype=string_dtype
        )
        parameters.create_dataset(
            "name",
            data=np.asarray(["epsilon_ij", "sigma_ij"], dtype=object),
            dtype=string_dtype,
        )
        parameters.create_dataset("ij_count", data=np.int64(2))
        pair.create_dataset("with_ele", data=np.bool_(False))
        pair.create_dataset("electrostatic_potential", data="", dtype=string_dtype)
        pair_data = pair.require_group("data/custom_pair")
        atom_count = handle["/topology/atom_count"][()]
        pair_data.create_dataset("atom_count", data=np.int64(atom_count))
        pair_data.create_dataset("type_count", data=np.int32(1))
        pair_data.create_dataset("pair_count", data=np.int64(1))
        pair_data.require_group("parameter").create_dataset(
            "value", data=np.asarray([[0.01], [1.0]], dtype=np.float32)
        )
        pair_data.create_dataset(
            "atom_type", data=np.zeros(atom_count, dtype=np.int32)
        )

        listed = handle.require_group("/forcefield/custom_force/listed")
        listed.create_dataset(
            "name", data=np.asarray(["custom_bond"], dtype=object), dtype=string_dtype
        )
        listed.create_dataset(
            "potential",
            data=np.asarray(
                ["E = k * (r_ab - r0) * (r_ab - r0);"], dtype=object
            ),
            dtype=string_dtype,
        )
        listed.create_dataset(
            "connected_atoms", data=np.asarray(["ab"], dtype=object), dtype=string_dtype
        )
        listed.create_dataset(
            "constrain_distance", data=np.asarray(["r0"], dtype=object), dtype=string_dtype
        )
        listed_parameters = listed.require_group("parameters")
        listed_parameters.create_dataset(
            "type",
            data=np.asarray(["int", "int", "float", "float"], dtype=object),
            dtype=string_dtype,
        )
        listed_parameters.create_dataset(
            "name",
            data=np.asarray(["atom_a", "atom_b", "k", "r0"], dtype=object),
            dtype=string_dtype,
        )
        listed_parameters.create_dataset("offset", data=np.asarray([0, 4], dtype=np.int64))
        listed_data = listed.require_group("data/custom_bond")
        listed_data.create_dataset("item_count", data=np.int64(1))
        listed_data.require_group("parameter").create_dataset(
            "value", data=np.asarray([[0.0, 1.0, 1.0, 1.5]], dtype=np.float32)
        )

        bond_soft = handle.require_group("/forcefield/bond_soft")
        bond_soft.create_dataset("atoms", data=np.asarray([[0, 1]], dtype=np.int32))
        bond_soft.create_dataset("k", data=np.asarray([1.0], dtype=np.float32))
        bond_soft.create_dataset("r0", data=np.asarray([1.5], dtype=np.float32))
        bond_soft.create_dataset("from_a_or_b", data=np.asarray([0], dtype=np.int32))


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
    with h5py.File(topology, "r") as topology_handle, h5py.File(
        protocol, "r"
    ) as protocol_handle, h5py.File(restart, "r") as handle:
        assert protocol_handle["/schema/name"][()].decode() == "sponge.protocol.h5"
        assert handle["/particles/all/position/value"].shape == (1, len(molecule.atoms), 3)
        assert handle["/particles/all/box/edges/value"].shape == (1, 3, 3)
        assert handle["/run/topology_hash"].asstr()[()] == topology_handle[
            "/topology/topology_hash"
        ].asstr()[()]
        assert handle["/run/atom_order_hash"].asstr()[()] == topology_handle[
            "/topology/atom_order_hash"
        ].asstr()[()]
        assert handle["/run/producer_protocol_hash"].asstr()[()] == protocol_handle[
            "/identity/content_hash"
        ].asstr()[()]
        assert handle["/run/state_hash"].asstr()[()].startswith("sha256:")


def test_save_sponge_input_wrapper_defaults_to_raw(tmp_path):
    molecule = _alanine_dipeptide()

    returned = Xponge.save_sponge_input(molecule, "system", tmp_path)

    assert returned is molecule
    assert (tmp_path / "system_coordinate.txt").is_file()
    assert not (tmp_path / "system_topology.spgt.h5").exists()


def test_save_sponge_input_wrapper_accepts_explicit_raw(tmp_path):
    molecule = _alanine_dipeptide()

    returned = Xponge.save_sponge_input(molecule, "system", tmp_path, format="raw")

    assert returned is molecule
    assert (tmp_path / "system_coordinate.txt").is_file()


def test_save_sponge_input_wrapper_dispatches_bundle(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()

    returned = Xponge.save_sponge_input(molecule, "system", tmp_path, format="bundle")

    assert returned is molecule
    assert (tmp_path / "system_topology.spgt.h5").is_file()
    assert (tmp_path / "system_protocol.spgp.h5").is_file()
    with h5py.File(tmp_path / "system_restart.spgr.h5", "r") as handle:
        assert handle["/particles/all/position/value"].shape == (1, len(molecule.atoms), 3)


def test_save_sponge_input_bundle_serializes_native_protocol_objects(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()
    protocol = Xponge.SpongeProtocol(
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="distance",
                type="distance",
                atom_indices=(0, 1),
                sigma=(0.5,),
                period=(0.0,),
            ),
        ),
        distance_constraints=(
            Xponge.ProtocolDistanceConstraints(
                atoms=((0, 1),),
                r0=(1.5,),
            ),
        ),
        positional_restraints=(
            Xponge.ProtocolPositionalRestraint(
                name="backbone",
                atom_indices=(0, 1),
                reference_coordinates=tuple(
                    (atom.x, atom.y, atom.z) for atom in molecule.atoms
                ),
                single_weight_default=2.0,
                refcoord_scaling_default="none",
                calc_virial_default=True,
            ),
        ),
        cv_restraints=(
            Xponge.ProtocolCVRestraint(
                name="umbrella",
                cv_refs=("distance",),
                weight=(1.0,),
                reference=(1.5,),
                start_step=(0,),
            ),
        ),
        metadynamics=(
            Xponge.ProtocolMetadynamics(
                name="native_bias",
                cv_refs=("distance",),
                grid_min=(0.0,),
                grid_max=(10.0,),
                grid_count=(64,),
                potential_update_interval_default=20,
                method_flags={"sink": 1},
            ),
        ),
        steering=Xponge.ProtocolSteering(
            cv_refs=("distance",),
            weight=(0.25,),
        ),
    )

    returned = Xponge.save_sponge_input(
        molecule,
        "system",
        tmp_path,
        format="bundle",
        protocol=protocol,
    )

    assert returned is molecule
    with h5py.File(tmp_path / "system_protocol.spgp.h5", "r") as handle:
        assert handle["/cv/distance/type"][()].decode() == "distance"
        assert np.array_equal(handle["/cv/distance/atom_indices"][:], [0, 1])
        assert np.allclose(handle["/constraint/default/pairs/r0"][:], [1.5])
        assert handle["/restraint/backbone/type"][()].decode() == "harmonic_positional"
        assert handle["/restraint/umbrella/type"][()].decode() == "cv_harmonic"
        assert handle["/meta/native_bias/cv_refs"][:].astype(str).tolist() == ["distance"]
        assert handle["/meta/native_bias/method_flags/sink"][()] == 1
        assert handle["/steer/cv_refs"][:].astype(str).tolist() == ["distance"]
        assert np.allclose(handle["/steer/weight"][:], [0.25])
        assert handle["/protocol/cv_count"][()] == 1
        assert handle["/protocol/restraint_count"][()] == 2
        assert handle["/protocol/enhanced_sampling/method"][()].decode() == "metadynamics"
        assert "/parameters/sponge/files/legacy_sidecars" not in handle
    with h5py.File(tmp_path / "system_restart.spgr.h5", "r") as handle:
        reference = handle[
            "/parameters/restart/references/restraint/backbone/coordinate"
        ][:]
        assert reference.shape == (len(molecule.atoms), 3)


def test_native_protocol_rejects_missing_cv_reference(tmp_path):
    protocol = Xponge.SpongeProtocol(
        cv_restraints=(
            Xponge.ProtocolCVRestraint(
                name="umbrella",
                cv_refs=("missing",),
                weight=(1.0,),
                reference=(1.0,),
            ),
        ),
    )

    with pytest.raises(BundleValidationError, match="missing or disabled CVs"):
        Xponge.save_sponge_input_bundle(
            _alanine_dipeptide(),
            "system",
            tmp_path,
            protocol=protocol,
        )


def test_save_sponge_input_bundle_serializes_native_sits(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()
    protocol = Xponge.SpongeProtocol(
        sits=Xponge.ProtocolSITS(
            mode="production",
            atom_indices=(0, 1),
            temperature_ladder=(280.0, 320.0),
            fb_interval=2,
            record_interval=3,
            update_interval=10,
            initial_nk=(1.0, 4.0),
            initial_log_norm=(0.0, 0.0),
            initial_log_nk=(0.0, np.log(4.0)),
        )
    )

    Xponge.save_sponge_input_bundle(
        molecule,
        "system",
        tmp_path,
        protocol=protocol,
    )

    with h5py.File(tmp_path / "system_protocol.spgp.h5", "r") as handle:
        assert handle["/sits/method/mode"][()].decode() == "production"
        assert handle["/sits/method/k_numbers"][()] == 2
        assert np.allclose(handle["/sits/method/temperature_ladder"][:], [280, 320])
        assert np.array_equal(handle["/sits/atom_indices"][:], [0, 1])
        assert handle["/protocol/enhanced_sampling/method"][()].decode() == "SITS"
    with h5py.File(tmp_path / "system_restart.spgr.h5", "r") as handle:
        root = "/parameters/restart/bias/sits/SITS"
        assert handle[root + "/schema_version"][()] == 1
        assert np.allclose(handle[root + "/nk"][:], [1.0, 4.0])
        assert np.allclose(handle[root + "/log_norm"][:], [0.0, 0.0])
        assert np.allclose(handle[root + "/log_nk"][:], [0.0, np.log(4.0)])


def test_native_sits_requires_restart_state_for_production(tmp_path):
    protocol = Xponge.SpongeProtocol(
        sits=Xponge.ProtocolSITS(
            mode="production",
            atom_numbers_policy="ALL",
            temperature_ladder=(280.0, 320.0),
        )
    )

    with pytest.raises(BundleValidationError, match="requires initial_nk"):
        Xponge.save_sponge_input_bundle(
            _alanine_dipeptide(),
            "system",
            tmp_path,
            protocol=protocol,
        )


def test_save_sponge_input_bundle_serializes_native_soft_walls(tmp_path):
    h5py = pytest.importorskip("h5py")
    protocol = Xponge.SpongeProtocol(
        soft_walls=(
            Xponge.ProtocolSoftWall(
                name="z_wall",
                potential="E = 0.5f * z * z;",
            ),
        )
    )

    Xponge.save_sponge_input_bundle(
        _alanine_dipeptide(),
        "system",
        tmp_path,
        protocol=protocol,
    )

    with h5py.File(tmp_path / "system_protocol.spgp.h5", "r") as handle:
        assert handle["/wall/soft/count"][()] == 1
        assert handle["/wall/soft/name"][:].astype(str).tolist() == ["z_wall"]
        assert handle["/wall/soft/potential"][:].astype(str).tolist() == [
            "E = 0.5f * z * z;"
        ]
        assert handle["/protocol/enhanced_sampling/method"][()].decode() == "soft_walls"


def test_save_sponge_input_bundle_serializes_native_hard_wall(tmp_path):
    h5py = pytest.importorskip("h5py")
    protocol = Xponge.SpongeProtocol(
        hard_wall=Xponge.ProtocolHardWall(
            bounds_low=(None, None, 1.5),
            bounds_high=(None, None, 12.0),
            allow_npt=True,
        )
    )

    Xponge.save_sponge_input_bundle(
        _alanine_dipeptide(),
        "system",
        tmp_path,
        protocol=protocol,
    )

    with h5py.File(tmp_path / "system_protocol.spgp.h5", "r") as handle:
        low = handle["/wall/hard/bounds_low"][:]
        high = handle["/wall/hard/bounds_high"][:]
        assert np.isneginf(low[:2]).all()
        assert np.isposinf(high[:2]).all()
        assert low[2] == pytest.approx(1.5)
        assert high[2] == pytest.approx(12.0)
        assert handle["/wall/hard/allow_npt"][()] == 1
        assert handle["/wall/hard"].attrs["allow_npt"] == 1


def test_native_hard_wall_requires_a_valid_finite_interval(tmp_path):
    with pytest.raises(BundleValidationError, match="at least one finite bound"):
        Xponge.save_sponge_input_bundle(
            _alanine_dipeptide(),
            "system",
            tmp_path,
            protocol=Xponge.SpongeProtocol(hard_wall=Xponge.ProtocolHardWall()),
        )

    with pytest.raises(BundleValidationError, match="z_low must be smaller"):
        Xponge.save_sponge_input_bundle(
            _alanine_dipeptide(),
            "system",
            tmp_path,
            protocol=Xponge.SpongeProtocol(
                hard_wall=Xponge.ProtocolHardWall(
                    bounds_low=(None, None, 2.0),
                    bounds_high=(None, None, 1.0),
                )
            ),
        )


def test_native_positional_restraint_requires_full_system_reference(tmp_path):
    molecule = _alanine_dipeptide()
    protocol = Xponge.SpongeProtocol(
        positional_restraints=(
            Xponge.ProtocolPositionalRestraint(
                name="backbone",
                atom_indices=(0, 1),
                reference_coordinates=((0.0, 0.0, 0.0), (1.5, 0.0, 0.0)),
                single_weight_default=2.0,
            ),
        ),
    )

    with pytest.raises(BundleValidationError, match="reference coordinates must have shape"):
        Xponge.save_sponge_input_bundle(
            molecule,
            "system",
            tmp_path,
            protocol=protocol,
        )


def test_save_sponge_input_rejects_protocol_in_raw_mode(tmp_path):
    with pytest.raises(ValueError, match="require format='bundle'"):
        Xponge.save_sponge_input(
            _alanine_dipeptide(),
            "system",
            tmp_path,
            format="raw",
            protocol=Xponge.SpongeProtocol(),
        )


def test_save_sponge_input_wrapper_rejects_unknown_format(tmp_path):
    with pytest.raises(ValueError, match="must be 'raw' or 'bundle'"):
        Xponge.save_sponge_input(
            _alanine_dipeptide(), "system", tmp_path, format="future"
        )


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


def test_save_sponge_input_bundle_supports_soft_bond(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()
    previous = Xponge.Molecule._save_functions.get("bond_soft")

    @Xponge.Molecule.Set_Save_SPONGE_Input("bond_soft")
    def write_soft_bond(_):
        return "1\n0 1 12.5 1.25 1"

    try:
        Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    finally:
        Xponge.Molecule.Del_Save_SPONGE_Input("bond_soft")
        if previous is not None:
            Xponge.Molecule.Set_Save_SPONGE_Input("bond_soft")(previous)

    with h5py.File(tmp_path / "system_topology.spgt.h5", "r") as handle:
        assert np.array_equal(handle["/forcefield/bond_soft/atoms"][:], [[0, 1]])
        assert np.allclose(handle["/forcefield/bond_soft/k"][:], [12.5])
        assert np.allclose(handle["/forcefield/bond_soft/r0"][:], [1.25])
        assert np.array_equal(handle["/forcefield/bond_soft/from_a_or_b"][:], [1])

    _write_saver_mdin(tmp_path, "system")
    restored_dir = tmp_path / "restored"
    convert_bundle_to_legacy(tmp_path, restored_dir, prefix="system")
    restored = parse_topology_file(
        "bond_soft_in_file", restored_dir / "system_bond_soft.txt"
    )
    restored_by_path = {dataset.path: dataset.data for dataset in restored}
    assert np.array_equal(restored_by_path["/forcefield/bond_soft/atoms"], [[0, 1]])
    assert np.array_equal(restored_by_path["/forcefield/bond_soft/from_a_or_b"], [1])


def test_save_sponge_input_bundle_supports_ryckaert_bellemans(tmp_path):
    h5py = pytest.importorskip("h5py")
    molecule = _alanine_dipeptide()
    assert classify_sponge_serializer_key("Ryckaert_Bellemans") == ("listed_force", None)
    previous_descriptor = Xponge.Molecule._save_functions.get("listed_forces")
    previous_data = Xponge.Molecule._save_functions.get("Ryckaert_Bellemans")

    @Xponge.Molecule.Set_Save_SPONGE_Input("listed_forces")
    def write_listed_forces(_):
        return (
            "[[[ Ryckaert_Bellemans ]]]\n"
            "[[ parameters ]]\n"
            "int atom_a, int atom_b, int atom_c, int atom_d, "
            "float c0, float c1, float c2, float c3, float c4, float c5\n"
            "[[ potential ]]\n"
            "E = c0 + c1;\n"
            "[[ end ]]\n"
        )

    @Xponge.Molecule.Set_Save_SPONGE_Input("Ryckaert_Bellemans")
    def write_ryckaert_bellemans(_):
        return "1\n0 1 2 3 0.1 0.2 0.3 0.4 0.5 0.6"

    try:
        Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    finally:
        for key, previous in (
            ("listed_forces", previous_descriptor),
            ("Ryckaert_Bellemans", previous_data),
        ):
            Xponge.Molecule.Del_Save_SPONGE_Input(key)
            if previous is not None:
                Xponge.Molecule.Set_Save_SPONGE_Input(key)(previous)

    root = "/forcefield/custom_force/listed/data/Ryckaert_Bellemans"
    with h5py.File(tmp_path / "system_topology.spgt.h5", "r") as handle:
        assert handle[f"{root}/item_count"][()] == 1
        assert np.array_equal(handle[f"{root}/parameter/int_value"][0, :4], [0, 1, 2, 3])
        assert np.allclose(handle[f"{root}/parameter/float_value"][0, 4:], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

    _write_saver_mdin(tmp_path, "system")
    restored_dir = tmp_path / "restored"
    convert_bundle_to_legacy(tmp_path, restored_dir, prefix="system")
    restored_lines = (restored_dir / "system_Ryckaert_Bellemans.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    assert restored_lines[0] == "1"
    assert np.allclose(
        [float(value) for value in restored_lines[1].split()],
        [0, 1, 2, 3, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    )


def test_save_sponge_input_bundle_rejects_escaping_prefix(tmp_path):
    with pytest.raises(BundlePathError, match="escapes output directory"):
        Xponge.save_sponge_input_bundle(_alanine_dipeptide(), "../escape", tmp_path)


def test_saver_bundle_without_velocity_starts_optional_sponge(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")

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


def test_native_protocol_bundle_starts_optional_sponge_without_sidecars(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    protocol = Xponge.SpongeProtocol(
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="distance",
                type="distance",
                atom_indices=(0, 1),
            ),
        ),
        cv_restraints=(
            Xponge.ProtocolCVRestraint(
                name="umbrella",
                cv_refs=("distance",),
                weight=(1.0,),
                reference=(1.5,),
            ),
        ),
        positional_restraints=(
            Xponge.ProtocolPositionalRestraint(
                name="backbone",
                atom_indices=(0, 1),
                reference_coordinates=tuple(
                    (atom.x, atom.y, atom.z) for atom in molecule.atoms
                ),
                single_weight_default=2.0,
                refcoord_scaling_default="none",
            ),
        ),
        steering=Xponge.ProtocolSteering(
            cv_refs=("distance",),
            weight=(0.25,),
        ),
    )
    Xponge.save_sponge_input_bundle(
        molecule,
        "system",
        tmp_path,
        protocol=protocol,
    )
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
    assert "END INITIALIZING STEER CV" in result.stdout
    assert "END INITIALIZING NATIVE H5 RESTRAIN" in result.stdout
    assert not (tmp_path / ".sponge_h5_restart_protocol").exists()
    assert not (tmp_path / ".sponge_h5_native_protocol" / "restraint").exists()


def test_native_metadynamics_bundle_starts_optional_sponge_without_cv_text(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    protocol = Xponge.SpongeProtocol(
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="distance",
                type="distance",
                atom_indices=(0, 1),
                sigma=(0.5,),
                period=(0.0,),
            ),
        ),
        metadynamics=(
            Xponge.ProtocolMetadynamics(
                name="native_bias",
                cv_refs=("distance",),
                grid_min=(0.0,),
                grid_max=(10.0,),
                grid_count=(8,),
                potential_update_interval_default=1000,
            ),
        ),
    )
    Xponge.save_sponge_input_bundle(
        molecule,
        "system",
        tmp_path,
        protocol=protocol,
    )
    _write_saver_mdin(tmp_path, "system", mode="nvt")
    with (tmp_path / "mdin.bundled.spg.toml").open("a", encoding="utf-8") as handle:
        handle.write(
            'dt = 0.001\n'
            'thermostat = "nose_hoover_chain"\n'
            'target_temperature = 300\n'
            '\n[nose_hoover_chain]\n'
            'length = 2\n'
        )

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
    assert not (tmp_path / ".sponge_h5_native_protocol" / "cv.txt").exists()


def test_native_sits_bundle_starts_optional_sponge_without_sits_text(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    protocol = Xponge.SpongeProtocol(
        sits=Xponge.ProtocolSITS(
            mode="production",
            atom_indices=(0, 1),
            temperature_ladder=(280.0, 320.0),
            update_interval=100,
            initial_nk=(1.0, 4.0),
        )
    )
    Xponge.save_sponge_input_bundle(
        molecule,
        "system",
        tmp_path,
        protocol=protocol,
    )
    _write_saver_mdin(
        tmp_path,
        "system",
        mode="nvt",
        restart_load="protocol",
    )
    with (tmp_path / "mdin.bundled.spg.toml").open("a", encoding="utf-8") as handle:
        handle.write(
            'dt = 0.001\n'
            'thermostat = "nose_hoover_chain"\n'
            'target_temperature = 300\n'
            '\n[nose_hoover_chain]\n'
            'length = 2\n'
        )

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
    assert "START INITIALIZING SITS" in result.stdout
    assert not (tmp_path / ".sponge_h5_native_protocol" / "sits.txt").exists()
    assert not (tmp_path / ".sponge_h5_native_protocol" / "sits_atom.txt").exists()


def test_native_soft_wall_bundle_starts_optional_sponge_without_text(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    protocol = Xponge.SpongeProtocol(
        soft_walls=(
            Xponge.ProtocolSoftWall(
                name="z_wall",
                potential="E = 0.5f * z * z;",
            ),
        )
    )
    Xponge.save_sponge_input_bundle(
        molecule,
        "system",
        tmp_path,
        protocol=protocol,
    )
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
    assert "START INITIALIZING SOFT WALLS FROM NATIVE H5" in result.stdout
    assert not (tmp_path / ".sponge_h5_native_protocol" / "soft_walls.txt").exists()


def test_native_custom_forces_start_optional_sponge_without_text(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")
    molecule = _alanine_dipeptide()
    molecule.box_length = [50.0, 50.0, 50.0]
    Xponge.save_sponge_input_bundle(molecule, "system", tmp_path)
    _add_native_custom_force_fixture(tmp_path / "system_topology.spgt.h5")
    _write_saver_mdin(tmp_path, "system", mode="minimization")
    with (tmp_path / "mdin.bundled.spg.toml").open("a", encoding="utf-8") as handle:
        handle.write("lambda_bond = 0.5\nsoft_bond_alpha = 0.2\n")

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
    assert "START INITIALIZING PAIRWISE FORCE FROM NATIVE H5" in result.stdout
    assert "START INITIALIZING LISTED FORCES FROM NATIVE H5" in result.stdout
    assert not (tmp_path / ".sponge_h5_native_custom_force").exists()


def test_reversed_saver_bundle_starts_optional_sponge(tmp_path):
    sponge = _optional_executable("SPONGE_EXECUTABLE", "SPONGE")

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
    probe = _optional_executable("SPONGE_H5_RESTART_PROBE")
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
