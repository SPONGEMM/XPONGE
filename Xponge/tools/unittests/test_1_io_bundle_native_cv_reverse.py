"""Native CV export and legacy round-trip coverage."""

import h5py
import numpy as np
import pytest

import Xponge as Xponge
from Xponge.io_bundle import convert_bundle_to_legacy, convert_legacy_to_bundle
from Xponge.io_bundle.errors import BundleExportError
from Xponge.io_bundle.state_parsers import _read_braced_sections


REFERENCE = ((1.125, 2.25, 3.5), (4.75, 5.125, 6.25))


def _bundle(tmp_path, *, virtual=False):
    import Xponge.forcefield.amber.ff14sb  # noqa: F401

    molecule = Xponge.get_peptide_from_sequence("AA")
    selection = {"atom_refs": ("center", 0)} if virtual else {"atom_indices": (1, 0)}
    protocol = Xponge.SpongeProtocol(
        virtual_atoms=(Xponge.ProtocolVirtualAtom(
            name="center", type="center", atom_indices=(1, 0), weight=(0.25, 0.75),
        ),) if virtual else (),
        collective_variables=(
            Xponge.ProtocolCollectiveVariable(
                name="rmsd_cv", type="rmsd", reference_coordinates=REFERENCE,
                rotate=False, sigma=(0.5,), period=(0.0,), **selection,
            ),
            Xponge.ProtocolCollectiveVariable(
                name="distance_cv", type="distance", atom_indices=(0, 1),
            ),
            Xponge.ProtocolCollectiveVariable(
                name="combined", type="combination",
                parameters={"CV": ("rmsd_cv", "distance_cv")}, function="x + y",
            ),
            Xponge.ProtocolCollectiveVariable(
                name="disabled", type="distance", atom_indices=(0, 1), enabled=False,
            ),
        ),
    )
    root = tmp_path / "bundle"
    Xponge.save_sponge_input_bundle(molecule, "system", root, protocol=protocol)
    (root / "mdin.bundled.spg.toml").write_text(
        'mode = "nvt"\n'
        'input_h5_topology_path = "system_topology.spgt.h5"\n'
        'input_h5_protocol_path = "system_protocol.spgp.h5"\n'
        'input_h5_restart_path = "system_restart.spgr.h5"\n'
        'input_h5_restart_load = "structural"\n'
    )
    return root


def _sections(path):
    return {name: dict(items) for name, items in _read_braced_sections(path)}


@pytest.mark.parametrize("virtual", [False, True])
@pytest.mark.parametrize("source", ["inline", "restart", "both"])
def test_native_cv_reference_round_trip(tmp_path, virtual, source):
    root = _bundle(tmp_path, virtual=virtual)
    if source in {"restart", "both"}:
        with h5py.File(root / "system_restart.spgr.h5", "a") as handle:
            handle.create_dataset("/parameters/restart/references/cv/rmsd_cv/coordinate",
                                  data=np.asarray(REFERENCE, dtype=np.float32))
    if source == "restart":
        with h5py.File(root / "system_protocol.spgp.h5", "a") as handle:
            del handle["/cv/rmsd_cv/coordinate"]
    legacy = tmp_path / "legacy"
    manifest = convert_bundle_to_legacy(root, legacy, prefix="system")
    sections = _sections(legacy / "system_cv.txt")
    rmsd = sections["rmsd_cv"]
    assert rmsd["CV_type"] == "rmsd"
    assert rmsd["atom"] == ("center 0" if virtual else "1 0")
    np.testing.assert_array_equal(
        np.fromstring(rmsd["coordinate"], sep=" ").reshape(2, 3), REFERENCE,
    )
    assert rmsd["rotate"] == "0"
    assert float(rmsd["sigma"]) == 0.5
    assert float(rmsd["period"]) == 0.0
    assert "coordinate_in_file" not in rmsd
    assert sections["combined"]["CV"] == "rmsd_cv distance_cv"
    assert sections["combined"]["function"] == "x + y"
    assert "disabled" not in sections
    if virtual:
        assert sections["center"] == {
            "vatom_type": "center", "atom": "1 0", "weight": "0.25 0.75",
        }
    assert 'cv_in_file = "system_cv.txt"' in (legacy / "mdin.legacy.spg.toml").read_text()
    assert any(entry.contract_id == "protocol.cv" and entry.status == "typed_exported"
               for entry in manifest.entries)

    # The forward converter represents legacy CVs as /cv/config; exporting
    # that bundle again must retain all reference values and atom ordering.
    converted = tmp_path / "converted"
    convert_legacy_to_bundle(legacy, converted, mdin="mdin.legacy.spg.toml")
    restored = tmp_path / "restored"
    convert_bundle_to_legacy(converted / "bundle", restored, prefix="system")
    assert _sections(restored / "system_cv.txt") == sections


@pytest.mark.parametrize("problem,error", [
    ("missing", "requires atoms and an RMSD reference"),
    ("rows", "must have shape"),
    ("columns", "must have shape"),
    ("nan", "finite"),
    ("overflow", "finite"),
    ("conflict", "conflicts"),
    ("wrong_type", "only supported for rmsd"),
    ("atom_selection", "mutually exclusive"),
    ("dimension", "dimension must be 1"),
    ("virtual_weight", "weight does not match"),
    ("disabled_virtual", "invalid atom reference"),
])
def test_invalid_native_cv_export_publishes_nothing(tmp_path, problem, error):
    root = _bundle(tmp_path, virtual=problem in {"virtual_weight", "disabled_virtual"})
    with h5py.File(root / "system_protocol.spgp.h5", "a") as handle:
        cv = handle["/cv/rmsd_cv"]
        if problem in {"missing", "rows", "columns", "nan", "overflow"}:
            del cv["coordinate"]
            if problem != "missing":
                value = np.asarray(REFERENCE)
                if problem == "rows":
                    value = value[:1]
                elif problem == "columns":
                    value = value[:, :2]
                else:
                    value[0, 0] = float("nan") if problem == "nan" else 1e100
                cv.create_dataset("coordinate", data=value)
        elif problem == "wrong_type":
            cv["type"][()] = "distance"
        elif problem == "atom_selection":
            cv.create_dataset("atom_refs", data=np.asarray(["1", "0"], dtype=h5py.string_dtype()))
        elif problem == "dimension":
            cv["dimension"][()] = 2
        elif problem == "virtual_weight":
            del handle["/cv/virtual_atom/center/weight"]
            handle.create_dataset("/cv/virtual_atom/center/weight", data=[1.0])
        elif problem == "disabled_virtual":
            handle["/cv/virtual_atom/center/enabled_default"][()] = 0
    if problem == "conflict":
        with h5py.File(root / "system_restart.spgr.h5", "a") as handle:
            handle.create_dataset("/parameters/restart/references/cv/rmsd_cv/coordinate",
                                  data=np.zeros((2, 3), dtype=np.float32))
    output = tmp_path / "legacy"
    with pytest.raises(BundleExportError, match=error):
        convert_bundle_to_legacy(root, output, prefix="system")
    assert not output.exists()


@pytest.mark.parametrize("conflict", [False, True])
def test_native_and_legacy_cv_fields_merge_or_reject_conflict(tmp_path, conflict):
    root = _bundle(tmp_path)
    with h5py.File(root / "system_protocol.spgp.h5", "a") as handle:
        config = handle.require_group("/cv/config")
        dtype = h5py.string_dtype()
        config.create_dataset("section/name", data=["rmsd_cv", "extra"], dtype=dtype)
        config.create_dataset("section/key_offset", data=[0, 1, 3])
        config.create_dataset("key", data=["rotate", "CV_type", "atom"], dtype=dtype)
        config.create_dataset("value", data=["1" if conflict else "0", "distance", "0 1"], dtype=dtype)
    output = tmp_path / "legacy"
    if conflict:
        with pytest.raises(BundleExportError, match="conflicts"):
            convert_bundle_to_legacy(root, output, prefix="system")
        assert not output.exists()
    else:
        convert_bundle_to_legacy(root, output, prefix="system")
        sections = _sections(output / "system_cv.txt")
        assert sections["extra"] == {"CV_type": "distance", "atom": "0 1"}
        assert "coordinate" in sections["rmsd_cv"]


def test_all_disabled_native_cvs_do_not_emit_empty_legacy_config(tmp_path):
    root = _bundle(tmp_path)
    with h5py.File(root / "system_protocol.spgp.h5", "a") as handle:
        for cv in handle["/cv"].values():
            cv["enabled_default"][()] = 0
    output = tmp_path / "legacy"
    convert_bundle_to_legacy(root, output, prefix="system")
    assert not (output / "system_cv.txt").exists()
    assert "cv_in_file" not in (output / "mdin.legacy.spg.toml").read_text()
