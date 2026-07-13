from __future__ import annotations

import h5py
import numpy as np
import pytest


mda = pytest.importorskip("MDAnalysis")

import Xponge.analysis.md_analysis as xmda  # noqa: E402
from Xponge.io_bundle.errors import (  # noqa: E402
    BundleTrajectoryError,
    UnverifiedBundlePairError,
)
from Xponge.io_bundle.bundle_builder import (  # noqa: E402
    BundleBuilder,
    BundleMetadata,
    BundlePaths,
)
from Xponge.tools.traj_analysis import load_Sponge_trajectory  # noqa: E402


def _text(group, name: str, value: str):
    return group.create_dataset(name, data=np.asarray(value, dtype=h5py.string_dtype("utf-8")))


def _write_topology(path, *, topology_hash="top-hash", atom_order_hash="order-hash"):
    with h5py.File(path, "w") as handle:
        schema = handle.create_group("schema")
        _text(schema, "name", "sponge.topology.h5")
        _text(schema, "version", "xponge.legacy_to_bundle.v1")

        topology = handle.create_group("topology")
        topology.create_dataset("atom_count", data=np.int64(2))
        _text(topology, "topology_hash", topology_hash)
        _text(topology, "atom_order_hash", atom_order_hash)

        atoms = handle.create_group("atoms")
        atoms.create_dataset("mass", data=np.asarray([12.011, 1.008]))
        charge = atoms.create_dataset("charge", data=np.asarray([18.2223, -18.2223]))
        charge.attrs["unit"] = "Amber"
        atoms.create_dataset("residue_index", data=np.asarray([0, 0], dtype=np.int32))
        handle.create_group("residues").create_dataset(
            "atom_offset", data=np.asarray([0, 2], dtype=np.int64)
        )

        parameters = handle.create_group("parameters")
        xponge = parameters.create_group("xponge")
        atom_parameters = xponge.create_group("atoms")
        string_dtype = h5py.string_dtype("utf-8")
        atom_parameters.create_dataset("name", data=np.asarray(["C", "H"], dtype=string_dtype))
        atom_parameters.create_dataset("type_name", data=np.asarray(["CT", "HC"], dtype=string_dtype))
        xponge.create_group("residues").create_dataset(
            "name", data=np.asarray(["MOL"], dtype=string_dtype)
        )
        handle.create_group("forcefield").create_group("bond").create_dataset(
            "atoms", data=np.asarray([[0, 1]], dtype=np.int32)
        )


def _write_bundle_trajectory(
    path,
    *,
    topology_hash="top-hash",
    atom_order_hash="order-hash",
    position_atoms=2,
):
    with h5py.File(path, "w") as handle:
        h5md = handle.create_group("h5md")
        h5md.attrs["version"] = np.asarray([1, 1], dtype=np.int32)
        creator = h5md.create_group("creator")
        creator.attrs["name"] = "XPONGE"
        creator.attrs["version"] = "test"

        sponge = handle.create_group("parameters").create_group("sponge")
        schema = sponge.create_group("schema")
        _text(schema, "name", "sponge.output.h5md")
        _text(schema, "version", "xponge.legacy_to_bundle.v1")
        output = sponge.create_group("output")
        _text(output, "status", "finalized")
        output.create_dataset("frame_count", data=np.asarray([2], dtype=np.int64))
        output.create_dataset(
            "particle_streams",
            data=np.asarray(["all"], dtype=h5py.string_dtype("utf-8")),
        )
        compatibility = sponge.create_group("topology_compatibility")
        _text(compatibility, "topology_hash", topology_hash)
        _text(compatibility, "atom_order_hash", atom_order_hash)

        stream = handle.create_group("particles").create_group("all")
        step = stream.create_dataset("step", data=np.asarray([10, 20], dtype=np.int64))
        time = stream.create_dataset("time", data=np.asarray([500.0, 1000.0]))
        time.attrs["unit"] = "fs"

        positions = np.zeros((2, position_atoms, 3), dtype=np.float32)
        positions[1, :, 0] = 0.2
        position = stream.create_group("position")
        position_value = position.create_dataset("value", data=positions)
        position_value.attrs["unit"] = "nm"
        position["step"] = step
        position["time"] = time

        velocity = stream.create_group("velocity")
        velocity_value = velocity.create_dataset(
            "value", data=np.full((2, position_atoms, 3), 0.1, dtype=np.float32)
        )
        velocity_value.attrs["unit"] = "nm ps-1"
        velocity["step"] = step
        velocity["time"] = time

        force = stream.create_group("force")
        force_value = force.create_dataset(
            "value", data=np.ones((2, position_atoms, 3), dtype=np.float32)
        )
        force_value.attrs["unit"] = "kcal mol-1 Angstrom-1"
        force["step"] = step
        force["time"] = time

        box = stream.create_group("box")
        box.attrs["dimension"] = 3
        box.attrs["boundary"] = np.asarray(["periodic"] * 3, dtype="S8")
        edges = box.create_group("edges")
        edge_value = edges.create_dataset(
            "value", data=np.asarray([np.eye(3), np.eye(3)], dtype=np.float32)
        )
        edge_value.attrs["unit"] = "nm"
        edges["step"] = step
        edges["time"] = time

        observable = handle.create_group("observables").create_group("all").create_group("energy")
        observable.create_dataset("value", data=np.asarray([-1.0, -2.0]))


def _write_legacy_trajectory(path, *, numbered=False):
    with h5py.File(path, "w") as handle:
        h5md = handle.create_group("h5md")
        h5md.attrs["version"] = np.asarray([1, 1], dtype=np.int32)
        creator = h5md.create_group("creator")
        creator.attrs["name"] = "SPONGE"
        creator.attrs["version"] = "legacy"
        particles = handle.create_group("particles")
        names = ["trajectory0", "trajectory1"] if numbered else ["trajectory"]
        for walker, name in enumerate(names):
            value = np.full((2, 2, 3), walker + 1, dtype=np.float32)
            particles.create_group(name).create_group("position").create_dataset("value", data=value)


def test_registered_formats_load_bundle_through_universe(tmp_path):
    topology = tmp_path / "topology.spgt.h5"
    trajectory = tmp_path / "trajectory.spg.h5md"
    _write_topology(topology)
    _write_bundle_trajectory(trajectory)

    universe = mda.Universe(
        topology,
        trajectory,
        topology_format="SPONGE_TOPOLOGY_H5",
        format="SPONGE_H5MD",
        particle_stream="all",
    )

    assert universe.trajectory.__class__ is xmda.SpongeH5MDReader
    assert universe.atoms.names.tolist() == ["C", "H"]
    assert universe.residues.resnames.tolist() == ["MOL"]
    assert universe.atoms.charges == pytest.approx([1.0, -1.0], rel=1e-5)
    assert len(universe.bonds) == 1
    assert len(universe.trajectory) == 2
    universe.trajectory[1]
    assert universe.atoms.positions[:, 0] == pytest.approx([2.0, 2.0])
    assert universe.trajectory.ts.velocities[0, 0] == pytest.approx(1.0)
    assert universe.trajectory.ts.forces[0, 0] == pytest.approx(4.184)
    assert universe.trajectory.ts.time == pytest.approx(1.0)
    assert universe.trajectory.ts.dt == pytest.approx(0.5)
    assert universe.trajectory.ts.data["step"] == 20
    assert universe.trajectory.ts.dimensions == pytest.approx([10, 10, 10, 90, 90, 90])
    assert universe.trajectory.ts.data["observables/energy"] == pytest.approx(-2.0)

    auto_detected = mda.Universe(topology, trajectory, particle_stream="all")
    assert auto_detected.trajectory.__class__ is xmda.SpongeH5MDReader

    cli_loaded = load_Sponge_trajectory(topology, trajectory, box=None)
    assert cli_loaded.trajectory.__class__ is xmda.SpongeH5MDReader

    class_selected = mda.Universe(
        topology,
        trajectory,
        topology_format=xmda.BundleTopologyParser,
        format=xmda.SpongeH5MDReader,
        particle_stream="all",
    )
    assert class_selected.atoms.positions == pytest.approx(universe.trajectory[0].positions)

    expected = np.asarray([frame.positions.copy() for frame in universe.trajectory])
    universe.transfer_to_memory()
    assert np.asarray([frame.positions for frame in universe.trajectory]) == pytest.approx(expected)


def test_wrapper_validates_hashes_before_using_registered_formats(tmp_path):
    topology = tmp_path / "topology.spgt.h5"
    trajectory = tmp_path / "trajectory.spg.h5md"
    _write_topology(topology)
    _write_bundle_trajectory(trajectory, topology_hash="different")

    with pytest.raises(UnverifiedBundlePairError, match="topology_hash"):
        xmda.load_bundle_universe(topology, trajectory)


def test_unified_reader_keeps_legacy_single_and_numbered_walkers(tmp_path):
    topology = tmp_path / "topology.spgt.h5"
    single = tmp_path / "single.h5md"
    numbered = tmp_path / "numbered.h5md"
    _write_topology(topology)
    _write_legacy_trajectory(single)
    _write_legacy_trajectory(numbered, numbered=True)

    single_universe = mda.Universe(
        topology,
        single,
        topology_format="SPONGE_TOPOLOGY_H5",
        format="SPONGE_H5MD",
    )
    assert single_universe.trajectory.layout == "legacy"
    assert single_universe.atoms.positions == pytest.approx(np.ones((2, 3)))
    historical_reader = xmda.SPONGEH5MDReader(single, n_atoms=2)
    assert historical_reader.layout == "legacy"
    historical_reader.close()

    walker_universe = mda.Universe(
        topology,
        numbered,
        topology_format="SPONGE_TOPOLOGY_H5",
        format="SPONGE_H5MD",
        walker=1,
    )
    assert walker_universe.trajectory.layout == "legacy"
    assert walker_universe.atoms.positions == pytest.approx(np.full((2, 3), 2.0))

    with pytest.raises(ValueError, match="walker indices"):
        mda.Universe(
            topology,
            numbered,
            topology_format="SPONGE_TOPOLOGY_H5",
            format="SPONGE_H5MD",
            walker=3,
        )


def test_reader_rejects_trajectory_atom_count_mismatch(tmp_path):
    topology = tmp_path / "topology.spgt.h5"
    trajectory = tmp_path / "trajectory.spg.h5md"
    _write_topology(topology)
    _write_bundle_trajectory(trajectory, position_atoms=3)

    with pytest.raises(BundleTrajectoryError, match="topology has 2"):
        mda.Universe(
            topology,
            trajectory,
            topology_format="SPONGE_TOPOLOGY_H5",
            format="SPONGE_H5MD",
        )


def test_historical_reader_name_is_supported_alias():
    assert xmda.SPONGEH5MDReader is xmda.SpongeH5MDReader
    assert xmda.register_mdanalysis_formats() == (
        "SPONGE_TOPOLOGY_H5",
        "SPONGE_H5MD",
    )


def test_bundle_builder_writes_mdanalysis_compatibility_metadata(tmp_path):
    paths = BundlePaths.canonical(tmp_path)
    builder = BundleBuilder(paths)
    builder.add_dataset("topology.spgt.h5", "/atoms/mass", np.asarray([12.011, 1.008]))
    builder.add_dataset(
        "topology.spgt.h5", "/atoms/charge", np.asarray([18.2223, -18.2223])
    )
    builder.add_dataset(
        "trajectory.spg.h5md",
        "/particles/all/position/value",
        np.zeros((1, 2, 3), dtype=np.float32),
    )
    builder.add_dataset(
        "trajectory.spg.h5md",
        "/particles/all/box/edges/value",
        np.asarray([np.eye(3)], dtype=np.float32),
    )
    builder.add_dataset(
        "trajectory.spg.h5md", "/particles/all/step", np.asarray([0], dtype=np.int64)
    )
    builder.add_dataset(
        "trajectory.spg.h5md", "/particles/all/time", np.asarray([0.0])
    )
    metadata = BundleMetadata(
        atom_count=2,
        topology_hash="top-hash",
        atom_order_hash="order-hash",
        forcefield_hash="force-hash",
        protocol_hash="protocol-hash",
    )
    builder.finalize({"topology.spgt.h5", "trajectory.spg.h5md"}, metadata)

    with h5py.File(paths.topology, "r") as topology:
        assert topology["/atoms/charge"].attrs["unit"] == "Amber"
    with h5py.File(paths.trajectory, "r") as trajectory:
        assert (
            trajectory["/parameters/sponge/topology_compatibility/topology_hash"].asstr()[()]
            == "top-hash"
        )
        assert (
            trajectory["/parameters/sponge/topology_compatibility/atom_order_hash"].asstr()[()]
            == "order-hash"
        )
        assert trajectory["/particles/all/box"].attrs["dimension"] == 3

    universe = xmda.load_bundle_universe(paths.topology, paths.trajectory)
    assert universe.atoms.n_atoms == 2
    assert universe.trajectory.n_frames == 1


def test_unified_reader_delegates_third_party_native_h5md(tmp_path):
    topology = tmp_path / "topology.spgt.h5"
    native = tmp_path / "native.h5md"
    _write_topology(topology)
    source = mda.Universe.empty(2, trajectory=True)
    source.atoms.positions = np.asarray([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    with mda.coordinates.H5MD.H5MDWriter(
        native,
        n_atoms=2,
        positions=True,
        velocities=False,
        forces=False,
        lengthunit="Angstrom",
        timeunit="ps",
    ) as writer:
        writer.write(source)

    universe = mda.Universe(
        topology,
        native,
        topology_format="SPONGE_TOPOLOGY_H5",
        format="SPONGE_H5MD",
    )
    assert universe.trajectory.layout == "h5md"
    assert universe.atoms.positions == pytest.approx(
        np.asarray([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    )
    assert not xmda.SpongeH5MDReader._format_hint(native)
    native_universe = mda.Universe(
        topology,
        native,
        topology_format="SPONGE_TOPOLOGY_H5",
    )
    assert native_universe.trajectory.__class__ is mda.coordinates.H5MD.H5MDReader
