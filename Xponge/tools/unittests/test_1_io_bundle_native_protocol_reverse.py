"""Native bias conversion and typed-only legacy import regression tests."""
import h5py
import numpy as np
import pytest
import Xponge as Xponge
from Xponge.io_bundle import convert_bundle_to_legacy, convert_legacy_to_bundle
from Xponge.io_bundle.errors import BundleExportError
from Xponge.io_bundle.legacy_case import parse_mdin_text
from Xponge.io_bundle.state_parsers import _read_braced_sections


def _bundle(tmp_path, *, disabled=False, policy=None):
    import Xponge.forcefield.amber.ff14sb
    molecule = Xponge.get_peptide_from_sequence("AA")
    if hasattr(molecule, "set_box_padding"):
        molecule.set_box_padding(20.0)
    else:
        molecule.box_length = [50., 50., 50.]
    atoms = list(molecule.atoms)
    ref = tuple((float(a.x), float(a.y), float(a.z)) for a in atoms)
    protocol = Xponge.SpongeProtocol(
        collective_variables=(Xponge.ProtocolCollectiveVariable(name="distance_cv", type="distance", atom_indices=(0, 1)),),
        positional_restraints=(Xponge.ProtocolPositionalRestraint(
            name="anchor", atom_indices=(1, 0), reference_coordinates=ref,
            weight=((1., 2., 3.), (4., 5., 6.)), refcoord_scaling_default="all",
            calc_virial_default=False, enabled=not disabled),),
        cv_restraints=(Xponge.ProtocolCVRestraint(name="bias_a", cv_refs=("distance_cv",), weight=(2.,), reference=(1.5,), enabled=not disabled),
                       Xponge.ProtocolCVRestraint(name="bias_b", cv_refs=("distance_cv",), weight=(3.,), reference=(2.5,),
                           period=(6.,), start_step=(1,), max_step=(2,), reduce_step=(3,), stop_step=(4,), enabled=not disabled)),
        steering=Xponge.ProtocolSteering(cv_refs=("distance_cv",), weight=(0.25,), enabled=not disabled),
        sits=Xponge.ProtocolSITS(mode="production", atom_indices=(0, 1) if policy is None else (), atom_numbers_policy=policy,
            temperature_ladder=(300., 450.), initial_nk=(1., 2.), nk_rest=True,
            nk_fix=True, fb_interval=2, record_interval=3, update_interval=4,
            cross_enhance_factor=0.75, enabled=not disabled),
    )
    root = tmp_path / "bundle"
    Xponge.save_sponge_input_bundle(molecule, "system", root, protocol=protocol)
    (root / "mdin.bundled.spg.toml").write_text(
        'mode = "nve"\nstep_limit = 0\ncutoff = 4\nskin = 1\n'
        'input_h5_topology_path = "system_topology.spgt.h5"\n'
        'input_h5_protocol_path = "system_protocol.spgp.h5"\n'
        'input_h5_restart_path = "system_restart.spgr.h5"\n'
        'input_h5_restart_load = "structural"\n')
    return root, ref


def _sections(path):
    return {name: dict(items) for name, items in _read_braced_sections(path)}


@pytest.mark.parametrize("policy", [None, "ALL", 2])
def test_native_protocol_roundtrip_without_sidecars(tmp_path, policy):
    root, ref = _bundle(tmp_path, policy=policy)
    legacy = tmp_path / "legacy"
    convert_bundle_to_legacy(root, legacy, prefix="system")
    commands = parse_mdin_text((legacy / "mdin.legacy.spg.toml").read_text())
    assert commands["restrain_refcoord_scaling"] == "all"
    assert commands["restrain_calc_virial"] == "0"
    assert commands["SITS_mode"] == "production"
    assert commands["SITS_T"] == "300/450"
    assert commands["SITS_nk_rest"] == "1"
    assert commands["SITS_k_numbers"] == "2"
    assert "SITS_in_file" not in commands
    if policy is not None:
        assert commands["SITS_atom_numbers"] == str(policy)
    np.testing.assert_array_equal(np.loadtxt(legacy / commands["restrain_atom_id"], dtype=int), [1, 0])
    np.testing.assert_allclose(np.loadtxt(legacy / commands["restrain_coordinate_in_file"], skiprows=1), ref, rtol=1e-6)
    bias = _sections(legacy / commands["restrain_cv_in_file"])["restrain"]
    assert bias["CV"] == "distance_cv distance_cv"
    assert bias["period"] == "0 6"
    assert bias["start_step"] == "0 1"
    assert bias["stop_step"] == "0 4"
    assert _sections(legacy / commands["steer_cv_in_file"])["steer"] == {"CV": "distance_cv", "weight": "0.25"}
    converted = tmp_path / "converted"
    convert_legacy_to_bundle(legacy, converted, mdin="mdin.legacy.spg.toml")
    assert not (converted / "bundle/legacy_sidecars").exists()
    for filename in ("topology.spgt.h5", "protocol.spgp.h5", "restart.spgr.h5"):
        with h5py.File(converted / "bundle" / filename) as h5:
            assert "/parameters/sponge/files/legacy_sidecars" not in h5
            assert not list(h5.get("/parameters/restart/protocol_sidecars", {}))
    restored = tmp_path / "restored"
    convert_bundle_to_legacy(converted / "bundle", restored, prefix="system")
    restored_commands = parse_mdin_text((restored / "mdin.legacy.spg.toml").read_text())
    for key in ("SITS_mode", "SITS_T", "restrain_refcoord_scaling", "restrain_calc_virial"):
        assert restored_commands[key] == commands[key]
    assert _sections(restored / commands["steer_cv_in_file"]) == _sections(legacy / commands["steer_cv_in_file"])


def test_disabled_native_protocol_is_not_reactivated(tmp_path):
    root, _ = _bundle(tmp_path, disabled=True)
    convert_bundle_to_legacy(root, tmp_path / "legacy", prefix="system")
    commands = parse_mdin_text((tmp_path / "legacy/mdin.legacy.spg.toml").read_text())
    assert not any(key.startswith(("restrain_", "steer_", "SITS_")) for key in commands)


@pytest.mark.parametrize("problem,error", [("reference", "shape"), ("weight", "shape"), ("cv", "missing or disabled"), ("log_norm", "cannot restore")])
def test_native_protocol_malformed_fails_before_writes(tmp_path, problem, error):
    root, _ = _bundle(tmp_path)
    if problem in {"reference", "log_norm"}:
        with h5py.File(root / "system_restart.spgr.h5", "a") as h5:
            if problem == "reference":
                path = "/parameters/restart/references/restraint/anchor/coordinate"
                del h5[path]; h5[path] = np.zeros((1, 3), np.float32)
            else:
                h5["/parameters/restart/bias/sits/SITS/log_norm"] = np.zeros(2, np.float32)
    else:
        with h5py.File(root / "system_protocol.spgp.h5", "a") as h5:
            path = "/steer/weight" if problem == "weight" else "/steer/cv_refs"
            del h5[path]
            h5[path] = np.zeros(2, np.float32) if problem == "weight" else np.asarray(["unknown"], dtype=h5py.string_dtype())
    with pytest.raises(BundleExportError, match=error):
        convert_bundle_to_legacy(root, tmp_path / "legacy", prefix="system")
    assert not (tmp_path / "legacy").exists()


def test_native_defaults_respect_explicit_mdin(tmp_path):
    root, _ = _bundle(tmp_path)
    with (root / "mdin.bundled.spg.toml").open("a") as out:
        out.write('\n[restrain]\ncalc_virial = true\nsingle_weight = 7\n[SITS]\nfb_interval = 10\n')
    convert_bundle_to_legacy(root, tmp_path / "legacy", prefix="system")
    commands = parse_mdin_text((tmp_path / "legacy/mdin.legacy.spg.toml").read_text())
    assert commands["restrain_calc_virial"] == "true"
    assert commands["restrain_single_weight"] == "7"
    assert commands["SITS_fb_interval"] == "10"
    assert commands["SITS_mode"] == "production"


def test_protocol_conversion_runs_without_sidecars(tmp_path):
    import os
    import subprocess
    executable = os.environ.get("SPONGE_EXECUTABLE")
    if not executable:
        pytest.skip("set SPONGE_EXECUTABLE for native/legacy/reimport runtime parity")
    root, _ = _bundle(tmp_path)
    legacy = tmp_path / "legacy"
    convert_bundle_to_legacy(root, legacy, prefix="system")
    converted = tmp_path / "converted"
    convert_legacy_to_bundle(legacy, converted, mdin="mdin.legacy.spg.toml")
    assert not (converted / "bundle/legacy_sidecars").exists()
    energies = []
    for case, mdin in ((root, "mdin.bundled.spg.toml"), (legacy, "mdin.legacy.spg.toml"), (converted / "bundle", "mdin.bundled.spg.toml")):
        import json
        commands = parse_mdin_text((case / mdin).read_text())
        commands.update(mode="nve", step_limit="1", dt="0.000001", print_zeroth_frame="1",
                        write_information_interval="1", write_mdout_interval="1", mdout="energy.out",
                        write_trajectory_interval="0", write_restart_file_interval="100")
        if case != legacy:
            commands["input_h5_restart_load"] = "protocol"
        (case / "runtime.spg.toml").write_text("\n".join(f"{key} = {json.dumps(value)}" for key, value in commands.items()) + "\n")
        command = [executable, "-mdin", "runtime.spg.toml"]
        result = subprocess.run(command, cwd=case, env=dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"),
                                text=True, capture_output=True, timeout=120)
        (case / "runtime.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr
        assert "END INITIALIZING STEER CV" in result.stdout
        rows = [line.split() for line in (case / "energy.out").read_text().splitlines() if line.strip()]
        assert len(rows) >= 2
        values = dict(zip(rows[0], map(float, rows[-1])))
        assert all(np.isfinite(value) for value in values.values())
        energies.append(values)
    common = set(energies[0]) & set(energies[1]) & set(energies[2])
    for key in common - {"time", "step"}:
        np.testing.assert_allclose([energy[key] for energy in energies], energies[0][key], rtol=2e-4, atol=2e-3, err_msg=key)
