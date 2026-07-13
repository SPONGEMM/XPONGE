"""
Unit tests for SPONGE legacy-to-bundle conversion scaffolding.
"""

from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from Xponge.io_bundle import convert_legacy_outputs_to_bundle, convert_legacy_to_bundle
from Xponge.io_bundle import converter as converter_module
from Xponge.io_bundle.legacy_case import parse_mdin_text, render_mdin_without_keys
from Xponge.io_bundle.output_parsers import parse_legacy_output_file
from Xponge.io_bundle.state_parsers import parse_protocol_or_restart_file
from Xponge.io_bundle.topology_parsers import (
    listed_force_schemas,
    pairwise_force_schema,
    parse_listed_force_data_file,
    parse_pairwise_force_data_file,
    parse_topology_file,
)
from Xponge.io_bundle.topology_exporters import export_improper
from Xponge.io_bundle.trajectory_parsers import parse_trajectory_file
from Xponge.tools.unittests.io_bundle_fixtures import write_basic_case


__all__ = [
    "test_parse_mdin_text",
    "test_improper_parser_uses_sponge_native_pk_schema",
    "test_improper_exporter_prefers_pk_and_reads_legacy_k",
    "test_typed_topology_parsers",
    "test_typed_state_parsers",
    "test_typed_trajectory_parsers",
    "test_legacy_to_bundle_dry_run_manifest",
    "test_legacy_to_bundle_default_legacy_output_plan",
    "test_legacy_to_bundle_h5_output_disables_default_legacy_outputs",
    "test_legacy_to_bundle_normal_mode_preserves_legacy_outputs",
    "test_legacy_to_bundle_bundled_mdin_sidecars_and_overrides",
    "test_legacy_to_bundle_static_default_prefix_generates_required_protocol",
    "test_legacy_to_bundle_dynamic_custom_force_default_prefix",
    "test_legacy_to_bundle_writes_typed_topology",
    "test_legacy_to_bundle_improper_alias_is_typed_only",
    "test_legacy_to_bundle_outputs_pass_sponge_h5_probes",
    "test_legacy_io_contract_matrix_dry_run",
    "test_legacy_output_mdout_parser",
    "test_legacy_output_mdout_name_alignment",
    "test_legacy_output_sits_and_meta_parsers",
    "test_legacy_output_qc_reaxff_parsers",
    "test_legacy_output_restart_parsers",
    "test_legacy_outputs_to_h5md_bundle",
    "test_legacy_outputs_diagnostic_bundle",
    "test_legacy_outputs_restart_bundle",
    "test_legacy_outputs_to_vds_bundle",
    "test_legacy_outputs_observable_path_splits_from_trajectory_bundle",
    "test_legacy_outputs_observable_only_bundle",
    "test_legacy_output_vds_mode_matrix_writes",
    "test_legacy_to_bundle_cli_dry_run",
    "test_legacy_outputs_to_bundle_cli_dry_run",
]


def _h5_string(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if hasattr(value, "decode"):
        return value.decode("utf-8")
    return str(value)


def _h5_string_list(dataset):
    return [_h5_string(value) for value in dataset[...]]


def test_parse_mdin_text():
    commands = parse_mdin_text(
        'mode = "npt"\n'
        'coordinate_in_file = "coord with spaces.txt"\n'
        "step_limit = 10 # ignored\n"
    )
    assert commands["mode"] == "npt"
    assert commands["coordinate_in_file"] == "coord with spaces.txt"
    assert commands["step_limit"] == "10"
    commands = parse_mdin_text(
        "[EAM]\n"
        'in_file = "eam.txt"\n'
        'atom_type_in_file = "eam_type.txt"\n'
        "[REAXFF]\n"
        'type_in_file = "reaxff_type.txt"\n'
        "[restrain]\n"
        'atom_id = "restrain_atom_id.txt"\n'
    )
    assert commands["EAM_in_file"] == "eam.txt"
    assert commands["EAM_atom_type_in_file"] == "eam_type.txt"
    assert commands["REAXFF_type_in_file"] == "reaxff_type.txt"
    assert commands["restrain_atom_id"] == "restrain_atom_id.txt"
    assert "in_file" not in commands
    rendered = render_mdin_without_keys(
        'mode = "nve"\n'
        '[EAM]\n'
        'in_file = "eam.txt"\n'
        'atom_type_in_file = "type.txt"\n'
        '[custom]\n'
        'keep = "value.txt"\n',
        {"EAM_in_file", "EAM_atom_type_in_file"},
        ['input_h5_topology_path = "topology.spgt.h5"'],
    )
    assert "[EAM]" not in rendered
    assert '[custom]\nkeep = "value.txt"' in rendered
    assert 'input_h5_topology_path = "topology.spgt.h5"' in rendered


def test_improper_parser_uses_sponge_native_pk_schema():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "improper.txt"
        path.write_text("1\n0 1 2 3 10.0 0.2\n", encoding="utf-8")

        for key in ("improper_dihedral_in_file", "improper_in_file"):
            datasets = {
                dataset.path: dataset.data
                for dataset in parse_topology_file(key, path)
            }
            assert np.array_equal(
                datasets["/forcefield/improper/atoms"],
                np.asarray([[0, 1, 2, 3]], dtype=np.int32),
            )
            assert np.allclose(
                datasets["/forcefield/improper/pk"],
                np.asarray([10.0], dtype=np.float32),
            )
            assert np.allclose(
                datasets["/forcefield/improper/phi0"],
                np.asarray([0.2], dtype=np.float32),
            )
            assert "/forcefield/improper/k" not in datasets


def test_improper_exporter_prefers_pk_and_reads_legacy_k():
    class Contract:
        bundle_file = "topology.spgt.h5"
        contract_id = "input.topology.improper"
        legacy_keys = ("improper_dihedral_in_file",)

    class Reader:
        def __init__(self, parameter_path):
            self.parameter_path = parameter_path
            self.read_paths = []

        def contains(self, bundle_file, dataset_path):
            assert bundle_file == Contract.bundle_file
            return dataset_path == self.parameter_path

        def read(self, bundle_file, dataset_path):
            assert bundle_file == Contract.bundle_file
            self.read_paths.append(dataset_path)
            values = {
                "/forcefield/improper/atoms": np.asarray(
                    [[0, 1, 2, 3]], dtype=np.int32
                ),
                self.parameter_path: np.asarray([10.0], dtype=np.float32),
                "/forcefield/improper/phi0": np.asarray([0.2], dtype=np.float32),
            }
            return values[dataset_path]

    for parameter_path in ("/forcefield/improper/pk", "/forcefield/improper/k"):
        reader = Reader(parameter_path)
        payloads = export_improper(Contract(), reader, None)
        assert len(payloads) == 1
        assert payloads[0].key == "improper_dihedral_in_file"
        assert payloads[0].data == "1\n0 1 2 3 10 0.200000003\n"
        assert parameter_path in reader.read_paths
        unused_path = (
            "/forcefield/improper/k"
            if parameter_path.endswith("/pk")
            else "/forcefield/improper/pk"
        )
        assert unused_path not in reader.read_paths


def test_typed_topology_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        case_dir.mkdir()
        write_basic_case(case_dir)
        (case_dir / "eam_setfl.txt").write_text(
            "comment 1\n"
            "comment 2\n"
            "comment 3\n"
            "2 Ni Cu\n"
            "2 0.5 2 0.1 5.0\n"
            "28 58.69 3.52 FCC\n"
            "1.0 2.0\n"
            "3.0 4.0\n"
            "29 63.55 3.61 FCC\n"
            "5.0 6.0\n"
            "7.0 8.0\n"
            "9.0 10.0\n"
            "11.0 12.0\n"
            "13.0 14.0\n",
            encoding="utf-8",
        )

        mass = {dataset.path: dataset.data for dataset in parse_topology_file("mass_in_file", case_dir / "mass.txt")}
        bond = {dataset.path: dataset.data for dataset in parse_topology_file("bond_in_file", case_dir / "bond.txt")}
        residue = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("residue_in_file", case_dir / "residue.txt")
        }
        cmap = {dataset.path: dataset.data for dataset in parse_topology_file("cmap_in_file", case_dir / "cmap.txt")}
        virtual_atom = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("virtual_atom_in_file", case_dir / "virtual_atom.txt")
        }
        lj_soft_core = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("LJ_soft_core_in_file", case_dir / "lj_soft_core.txt")
        }
        eam = {dataset.path: dataset.data for dataset in parse_topology_file("EAM_in_file", case_dir / "eam.txt")}
        eam_setfl = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("EAM_in_file", case_dir / "eam_setfl.txt")
        }
        eam_atom_type = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("EAM_atom_type_in_file", case_dir / "eam_atom_type.txt")
        }
        sw = {dataset.path: dataset.data for dataset in parse_topology_file("SW_in_file", case_dir / "sw.txt")}
        edip = {dataset.path: dataset.data for dataset in parse_topology_file("EDIP_in_file", case_dir / "edip.txt")}
        tersoff = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("TERSOFF_in_file", case_dir / "tersoff.txt")
        }
        qc_type = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("qc_type_in_file", case_dir / "qc_type.txt")
        }
        reaxff_type = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("REAXFF_type_in_file", case_dir / "reaxff_type.txt")
        }
        reaxff = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("REAXFF_in_file", case_dir / "reaxff.txt")
        }
        pairwise_force = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("pairwise_force_in_file", case_dir / "pairwise_force.txt")
        }
        pairwise_name, pairwise_types, pairwise_names, pairwise_ij_count = pairwise_force_schema(
            case_dir / "pairwise_force.txt"
        )
        pairwise_data = {
            dataset.path: dataset.data
            for dataset in parse_pairwise_force_data_file(
                pairwise_name,
                case_dir / "custom_pair.txt",
                parameter_types=pairwise_types,
                parameter_names=pairwise_names,
                ij_parameter_count=pairwise_ij_count,
            )
        }
        listed_forces = {
            dataset.path: dataset.data
            for dataset in parse_topology_file("listed_forces_in_file", case_dir / "listed_forces.txt")
        }
        listed_name, listed_types, listed_names = listed_force_schemas(case_dir / "listed_forces.txt")[0]
        listed_data = {
            dataset.path: dataset.data
            for dataset in parse_listed_force_data_file(
                listed_name,
                case_dir / "custom_bond.txt",
                parameter_types=listed_types,
                parameter_names=listed_names,
            )
        }

        assert np.allclose(mass["/atoms/mass"], np.asarray([1.0, 16.0], dtype=np.float32))
        assert np.array_equal(bond["/forcefield/bond/atoms"], np.asarray([[0, 1]], dtype=np.int32))
        assert np.allclose(bond["/forcefield/bond/r0"], np.asarray([1.5], dtype=np.float32))
        assert np.array_equal(residue["/residues/atom_offset"], np.asarray([0, 2], dtype=np.int64))
        assert np.array_equal(cmap["/forcefield/cmap/atoms"], np.asarray([[0, 1, 0, 1, 0]], dtype=np.int32))
        assert np.allclose(cmap["/forcefield/cmap/grid_value"], np.asarray([1.0, 2.0, 3.0, 4.0], dtype=np.float32))
        assert np.array_equal(virtual_atom["/forcefield/virtual_atom/from"], np.asarray([0, 1, 0], dtype=np.int32))
        assert np.allclose(virtual_atom["/forcefield/virtual_atom/parameter"], np.asarray([0.5, 0.5], dtype=np.float32))
        assert np.array_equal(lj_soft_core["/forcefield/lj_soft_core/atom_type_A"], np.asarray([0, 0], dtype=np.int32))
        assert np.allclose(lj_soft_core["/forcefield/lj_soft_core/pair_BB"], np.asarray([4.0], dtype=np.float32))
        assert eam["/manybody/eam/format"] == "funcfl"
        assert np.array_equal(eam["/manybody/eam/atomic_number"], np.asarray([29], dtype=np.int32))
        assert np.allclose(eam["/manybody/eam/embed/value"], np.asarray([[23.060548, 46.121096]], dtype=np.float32))
        assert np.allclose(eam["/manybody/eam/electron_density/value"], np.asarray([[5.0, 6.0]], dtype=np.float32))
        assert eam_setfl["/manybody/eam/format"] == "setfl"
        assert np.array_equal(eam_setfl["/manybody/eam/type_name"], np.asarray(["Ni", "Cu"], dtype=object))
        assert np.allclose(eam_setfl["/manybody/eam/embed/raw_ev"][1], np.asarray([5.0, 6.0], dtype=np.float32))
        assert np.allclose(
            eam_setfl["/manybody/eam/pair_potential/value"][1, 0, 1],
            12.0 / 0.1 * 23.060548,
        )
        assert np.allclose(
            eam_setfl["/manybody/eam/pair_potential/value"][0, 1],
            eam_setfl["/manybody/eam/pair_potential/value"][1, 0],
        )
        assert np.array_equal(eam_atom_type["/manybody/eam/atom_type"], np.asarray([0, 1], dtype=np.int32))
        assert np.array_equal(sw["/manybody/sw/pair/type"], np.asarray([[0, 0]], dtype=np.int32))
        assert np.allclose(sw["/manybody/sw/triple/parameters"], np.asarray([[9.0, 10.0, 11.0]], dtype=np.float32))
        assert np.array_equal(edip["/manybody/edip/atom_type"], np.asarray([0, 0], dtype=np.int32))
        assert np.allclose(
            edip["/manybody/edip/triple/parameters"],
            np.asarray([[9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]], dtype=np.float32),
        )
        assert np.array_equal(tersoff["/manybody/tersoff/type_name"], np.asarray(["Si"], dtype=object))
        assert np.array_equal(tersoff["/manybody/tersoff/entry/type"], np.asarray([[0, 0, 0]], dtype=np.int32))
        assert np.allclose(tersoff["/manybody/tersoff/entry/parameters_raw"][0, 13], 1380.0)
        assert np.allclose(tersoff["/manybody/tersoff/entry/parameters"][0, 13], 1380.0 * 23.0605480)
        assert np.array_equal(tersoff["/manybody/tersoff/map"], np.asarray([0], dtype=np.int32))
        assert np.array_equal(tersoff["/manybody/tersoff/atom_type"], np.asarray([0, 0], dtype=np.int32))
        assert np.array_equal(qc_type["/qc/type/atom_index"], np.asarray([0, 1], dtype=np.int32))
        assert np.array_equal(qc_type["/qc/type/symbol"], np.asarray(["H", "O"], dtype=object))
        assert qc_type["/qc/type/charge"] == 0
        assert qc_type["/qc/type/multiplicity"] == 1
        assert np.array_equal(reaxff_type["/manybody/reaxff/type/name"], np.asarray(["O", "H"], dtype=object))
        assert reaxff_type["/manybody/reaxff/type/count"] == 2
        assert reaxff["/manybody/reaxff/parameters/general/count"] == 2
        assert np.allclose(reaxff["/manybody/reaxff/parameters/general/value"], np.asarray([50.0, 9.5], dtype=np.float32))
        assert np.array_equal(
            reaxff["/manybody/reaxff/parameters/atom/type_name"],
            np.asarray(["O", "H"], dtype=object),
        )
        assert np.array_equal(reaxff["/manybody/reaxff/parameters/atom/value_offset"], np.asarray([0, 32, 64], dtype=np.int64))
        assert np.array_equal(reaxff["/manybody/reaxff/parameters/bond/type"], np.asarray([[1, 2]], dtype=np.int32))
        assert np.allclose(reaxff["/manybody/reaxff/parameters/bond/value"][:2], np.asarray([170.0, 0.0], dtype=np.float32))
        assert np.array_equal(reaxff["/manybody/reaxff/parameters/angle/type"], np.asarray([[1, 2, 1]], dtype=np.int32))
        assert np.allclose(reaxff["/manybody/reaxff/parameters/hydrogen_bond/value"], np.asarray([[1.9, -4.4, 1.7, 3.0]], dtype=np.float32))
        assert pairwise_force["/forcefield/custom_force/pairwise/name"] == "custom_pair"
        assert np.array_equal(
            pairwise_force["/forcefield/custom_force/pairwise/parameters/name"],
            np.asarray(["epsilon_ij", "sigma_ij"], dtype=object),
        )
        assert pairwise_force["/forcefield/custom_force/pairwise/parameters/ij_count"] == 2
        assert pairwise_force["/forcefield/custom_force/pairwise/with_ele"] == np.asarray(False, dtype=np.bool_)
        assert np.allclose(
            pairwise_data["/forcefield/custom_force/pairwise/data/custom_pair/parameter/value"],
            np.asarray([[1.0], [2.0]], dtype=np.float32),
        )
        assert np.array_equal(
            pairwise_data["/forcefield/custom_force/pairwise/data/custom_pair/atom_type"],
            np.asarray([0, 0], dtype=np.int32),
        )
        assert np.allclose(
            pairwise_data["/forcefield/custom_force/pairwise/data/custom_pair/parameter/float_value"],
            np.asarray([[1.0], [2.0]], dtype=np.float32),
        )
        assert np.array_equal(
            listed_forces["/forcefield/custom_force/listed/name"],
            np.asarray(["custom_bond"], dtype=object),
        )
        assert np.array_equal(
            listed_forces["/forcefield/custom_force/listed/parameters/name"],
            np.asarray(["atom_a", "atom_b", "k", "r0"], dtype=object),
        )
        assert np.array_equal(listed_forces["/forcefield/custom_force/listed/parameters/is_atom"], np.asarray([1, 1, 0, 0], dtype=np.int32))
        assert np.array_equal(listed_forces["/forcefield/custom_force/listed/connected_atoms"], np.asarray(["ab"], dtype=object))
        assert np.allclose(
            listed_data["/forcefield/custom_force/listed/data/custom_bond/parameter/value"],
            np.asarray([[0.0, 1.0, 10.0, 1.5]], dtype=np.float32),
        )
        assert np.array_equal(
            listed_data["/forcefield/custom_force/listed/data/custom_bond/parameter/int_value"],
            np.asarray([[0, 1, 0, 0]], dtype=np.int32),
        )


def test_typed_state_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        case_dir.mkdir()
        write_basic_case(case_dir)

        sits_nk = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("SITS_nk_in_file", case_dir / "sits_nk.txt")
        }
        atom_id = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("restrain_atom_id", case_dir / "restrain_atom_id.txt")
        }
        weight = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("restrain_weight_in_file", case_dir / "restrain_weight.txt")
        }
        coordinate = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file(
                "restrain_coordinate_in_file", case_dir / "restrain_coordinate.txt"
            )
        }
        nhc_restart = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file(
                "nose_hoover_chain_restart_input", case_dir / "nhc_restart.txt"
            )
        }
        constrain = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("constrain_in_file", case_dir / "constrain.txt")
        }
        sits_atom = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("SITS_atom_in_file", case_dir / "sits_atom.txt")
        }
        soft_walls = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("soft_walls_in_file", case_dir / "soft_walls.txt")
        }
        cv = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("cv_in_file", case_dir / "cv.txt")
        }
        restrain = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("restrain_in_file", case_dir / "restrain.txt")
        }
        restrain_cv = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("restrain_cv_in_file", case_dir / "restrain_cv.txt")
        }
        steer_cv = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("steer_cv_in_file", case_dir / "steer_cv.txt")
        }
        sits = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("SITS_in_file", case_dir / "sits.txt")
        }
        meta_edge = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("meta_edge_in_file", case_dir / "meta_edge.txt")
        }
        meta_potential = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("meta_potential_in_file", case_dir / "meta_potential.txt")
        }
        meta_scatter = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("meta_scatter_in_file", case_dir / "meta_scatter.txt")
        }
        hills = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("hills_in_file", case_dir / "hills.txt")
        }
        amber_rst7 = {
            dataset.path: dataset.data
            for dataset in parse_protocol_or_restart_file("restrain_amber_rst7", case_dir / "restrain.rst7")
        }

        assert np.allclose(sits_nk["/parameters/restart/bias/sits/SITS/nk"], np.asarray([1.0, 2.0], dtype=np.float32))
        assert np.array_equal(cv["/cv/config/section/name"], np.asarray(["print", "distance"], dtype=object))
        assert np.array_equal(cv["/cv/config/section/key_offset"], np.asarray([0, 1, 3], dtype=np.int64))
        assert np.array_equal(cv["/cv/config/key"], np.asarray(["CV", "CV_type", "atom"], dtype=object))
        assert np.array_equal(restrain["/restraint/config/section/name"], np.asarray(["print", "distance"], dtype=object))
        assert np.array_equal(restrain_cv["/restraint/cv/config/key"], np.asarray(["CV", "weight", "reference"], dtype=object))
        assert np.array_equal(steer_cv["/steer/config/value"], np.asarray(["distance", "2.0"], dtype=object))
        assert np.array_equal(sits["/sits/config/key"], np.asarray(["CV", "nk"], dtype=object))
        assert np.array_equal(meta_edge["/meta/default/grid/coordinate"], np.asarray([[0.0, 0.0], [1.0, 0.0]], dtype=np.float32))
        assert np.allclose(meta_edge["/meta/default/grid/force"], np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32))
        assert np.array_equal(meta_potential["/parameters/restart/bias/meta/default/potential/grid"], np.asarray([2, 2], dtype=np.int32))
        assert np.allclose(meta_potential["/parameters/restart/bias/meta/default/potential/value"], np.asarray([3.0, 4.0, 5.0, 6.0], dtype=np.float32))
        assert meta_potential["/parameters/restart/bias/meta/default/potential_export"].startswith("Meta potential")
        assert np.array_equal(meta_scatter["/parameters/restart/bias/meta/default/scatter/grid"], np.asarray([2], dtype=np.int32))
        assert np.allclose(meta_scatter["/parameters/restart/bias/meta/default/scatter/coordinate"], np.asarray([[0.0], [0.5]], dtype=np.float32))
        assert hills["/parameters/restart/bias/meta/default/hills"] == "0.1 0.2 1.5\n0.3 0.4 1.6\n"
        assert np.allclose(hills["/parameters/restart/bias/meta/default/hills_typed/value"], np.asarray([[0.1, 0.2, 1.5], [0.3, 0.4, 1.6]], dtype=np.float32))
        assert np.allclose(
            amber_rst7["/parameters/restart/references/restraint/default/coordinate"],
            np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
        )
        assert np.allclose(
            amber_rst7["/parameters/restart/references/restraint/default/velocity"],
            np.asarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype=np.float32),
        )
        assert np.array_equal(atom_id["/restraint/default/atom_indices"], np.asarray([0, 1], dtype=np.int32))
        assert np.allclose(
            weight["/restraint/default/weight"],
            np.asarray([[10.0, 0.0, 0.0], [0.0, 20.0, 0.0]], dtype=np.float32),
        )
        assert np.allclose(
            coordinate["/parameters/restart/references/restraint/default/coordinate"],
            np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
        )
        assert np.allclose(
            nhc_restart["/parameters/restart/thermostat/nose_hoover_chain"],
            np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
        )
        assert np.array_equal(constrain["/constraint/default/pairs/atoms"], np.asarray([[0, 1]], dtype=np.int32))
        assert np.allclose(constrain["/constraint/default/pairs/r0"], np.asarray([1.5], dtype=np.float32))
        assert np.array_equal(sits_atom["/sits/atom_indices"], np.asarray([0, 1], dtype=np.int32))
        assert np.array_equal(soft_walls["/wall/soft/name"], np.asarray(["z_wall"], dtype=object))
        assert soft_walls["/wall/soft/potential"][0] == "E = powf((z - 5.0f) * (z - 5.0f), -6.0f);"


def test_typed_trajectory_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        case_dir.mkdir()
        write_basic_case(case_dir)

        crd = {
            dataset.path: dataset.data
            for dataset in parse_trajectory_file("crd", case_dir / "traj.dat", atom_count=2)
        }
        box = {
            dataset.path: dataset.data
            for dataset in parse_trajectory_file("box", case_dir / "traj_box.dat", atom_count=2)
        }
        vel = {
            dataset.path: dataset.data
            for dataset in parse_trajectory_file("vel", case_dir / "traj_vel.dat", atom_count=2)
        }

        assert np.allclose(crd["/particles/all/position/value"][1, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))
        assert np.array_equal(crd["/particles/all/step"], np.asarray([0, 1], dtype=np.int64))
        assert np.allclose(box["/particles/all/box/edges/value"][1], np.diag([11.0, 21.0, 31.0]).astype(np.float32))
        assert np.allclose(vel["/particles/all/velocity/value"][0, 0], np.asarray([0.1, 0.2, 0.3], dtype=np.float32))


def test_legacy_to_bundle_dry_run_manifest():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        write_basic_case(case_dir)

        manifest = convert_legacy_to_bundle(case_dir, output_dir, dry_run=True)
        payload = manifest.to_dict()
        entries = {entry["contract_id"]: entry for entry in payload["entries"]}

        assert payload["schema"] == "xponge.legacy_to_bundle.manifest"
        assert payload["mode"] == "rerun"
        assert entries["restart.coordinate"]["status"] == "converted"
        assert entries["restart.velocity"]["status"] == "converted"
        assert entries["restart.box"]["status"] == "converted"
        assert entries["protocol.cv"]["status"] == "typed_converted"
        assert entries["topology.mass"]["status"] == "typed_converted"
        assert entries["topology.charge"]["status"] == "typed_converted"
        assert entries["topology.residue"]["status"] == "typed_converted"
        assert entries["topology.bond"]["status"] == "typed_converted"
        assert entries["topology.angle"]["status"] == "typed_converted"
        assert entries["topology.dihedral"]["status"] == "typed_converted"
        assert entries["topology.improper_dihedral"]["status"] == "typed_converted"
        assert entries["topology.nb14_extra"]["status"] == "typed_converted"
        assert entries["topology.urey_bradley"]["status"] == "typed_converted"
        assert entries["topology.cmap"]["status"] == "typed_converted"
        assert entries["topology.gb"]["status"] == "typed_converted"
        assert entries["topology.virtual_atom"]["status"] == "typed_converted"
        assert entries["topology.LJ"]["status"] == "typed_converted"
        assert entries["topology.LJ_soft_core"]["status"] == "typed_converted"
        assert entries["topology.subsys_division"]["status"] == "typed_converted"
        assert entries["topology.EAM"]["status"] == "typed_converted"
        assert entries["topology.EAM_atom_type"]["status"] == "typed_converted"
        assert entries["topology.SW"]["status"] == "typed_converted"
        assert entries["topology.EDIP"]["status"] == "typed_converted"
        assert entries["topology.pairwise_force"]["status"] == "typed_converted"
        assert entries["topology.pairwise_force_data.custom_pair"]["status"] == "typed_converted"
        assert entries["topology.listed_forces"]["status"] == "typed_converted"
        assert entries["topology.listed_force_data.custom_bond"]["status"] == "typed_converted"
        assert entries["topology.TERSOFF"]["status"] == "typed_converted"
        assert entries["topology.qc_type"]["status"] == "typed_converted"
        assert entries["topology.REAXFF"]["status"] == "typed_converted"
        assert entries["topology.REAXFF_type"]["status"] == "typed_converted"
        assert entries["restart.SITS_nk"]["status"] == "typed_converted"
        assert entries["protocol.restrain"]["status"] == "typed_converted"
        assert entries["protocol.restrain_atom_id"]["status"] == "typed_converted"
        assert entries["protocol.restrain_weight"]["status"] == "typed_converted"
        assert entries["protocol.restrain_cv"]["status"] == "typed_converted"
        assert entries["protocol.steer_cv"]["status"] == "typed_converted"
        assert entries["restart.restrain_coordinate"]["status"] == "typed_converted"
        assert entries["restart.nose_hoover_chain_restart_input"]["status"] == "typed_converted"
        assert entries["protocol.constrain"]["status"] == "typed_converted"
        assert entries["protocol.SITS"]["status"] == "typed_converted"
        assert entries["protocol.SITS_atom"]["status"] == "typed_converted"
        assert entries["protocol.meta_edge"]["status"] == "typed_converted"
        assert entries["restart.meta_potential"]["status"] == "typed_converted"
        assert entries["restart.meta_scatter"]["status"] == "typed_converted"
        assert entries["restart.hills"]["status"] == "typed_converted"
        assert entries["protocol.soft_walls"]["status"] == "typed_converted"
        assert entries["run_mdin.step_limit"]["status"] == "preserved_in_mdin"
        assert entries["output.h5.output_h5_trajectory_path"]["status"] == "output_plan_preserved"
        assert entries["output.legacy_sidecar.mdout"]["status"] == "legacy_output_sidecar_preserved"
        assert entries["trajectory.crd"]["status"] == "typed_converted"
        assert entries["trajectory.box"]["status"] == "typed_converted"
        assert entries["trajectory.vel"]["status"] == "typed_converted"
        assert "rerun.trajectory" not in entries


def test_legacy_to_bundle_default_legacy_output_plan():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'mode = "nvt"\n'
            'default_out_file_prefix = "prod"\n'
            "step_limit = 1\n",
            encoding="utf-8",
        )

        manifest = convert_legacy_to_bundle(case_dir, output_dir, dry_run=True)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}

        assert entries["output.legacy_sidecar.mdout"]["status"] == "legacy_output_sidecar_default"
        assert entries["output.legacy_sidecar.mdout"]["source_path"] == "prod.out"
        assert entries["output.legacy_sidecar.mdinfo"]["source_path"] == "prod.info"
        assert entries["output.legacy_sidecar.crd"]["source_path"] == "prod.dat"
        assert entries["output.legacy_sidecar.box"]["source_path"] == "prod.box"
        assert entries["output.legacy_sidecar.rst"]["source_path"] == "prod"
        assert "output.legacy_sidecar.vel" not in entries
        assert "output.legacy_sidecar.frc" not in entries
        assert entries["run_mdin.default_out_file_prefix"]["status"] == "preserved_in_mdin"


def test_legacy_to_bundle_h5_output_disables_default_legacy_outputs():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'mode = "nvt"\n'
            'default_out_file_prefix = "prod"\n'
            'output_h5_trajectory_path = "prod.spg.h5md"\n'
            'mdout = "keep_mdout.txt"\n'
            "step_limit = 1\n",
            encoding="utf-8",
        )

        manifest = convert_legacy_to_bundle(case_dir, output_dir, dry_run=True)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}

        assert entries["output.h5.output_h5_trajectory_path"]["status"] == "output_plan_preserved"
        assert entries["output.legacy_sidecar.mdout"]["status"] == "legacy_output_sidecar_preserved"
        assert entries["run_mdin.default_out_file_prefix"]["status"] == "preserved_in_mdin"
        assert "output.legacy_sidecar.mdinfo" not in entries
        assert "output.legacy_sidecar.crd" not in entries
        assert "output.legacy_sidecar.box" not in entries
        assert "output.legacy_sidecar.rst" not in entries


def test_legacy_to_bundle_normal_mode_preserves_legacy_outputs():
    originals = (
        converter_module.ensure_dataset,
        converter_module.ensure_group,
        converter_module.ensure_hard_link,
        converter_module.set_attrs,
        converter_module.write_string,
        converter_module.write_bytes,
        converter_module.write_string_array,
    )
    converter_module.ensure_dataset = lambda *args, **kwargs: None
    converter_module.ensure_group = lambda *args, **kwargs: None
    converter_module.ensure_hard_link = lambda *args, **kwargs: None
    converter_module.set_attrs = lambda *args, **kwargs: None
    converter_module.write_string = lambda *args, **kwargs: None
    converter_module.write_bytes = lambda *args, **kwargs: None
    converter_module.write_string_array = lambda *args, **kwargs: None
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "case"
            output_dir = Path(tmpdir) / "out"
            case_dir.mkdir()
            write_basic_case(case_dir)
            mdin_path = case_dir / "mdin.spg.toml"
            mdin_path.write_text(
                mdin_path.read_text(encoding="utf-8").replace('mode = "rerun"', 'mode = "nvt"'),
                encoding="utf-8",
            )

            manifest = converter_module.convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
            entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
            bundled_mdin = (output_dir / "bundle" / "mdin.bundled.spg.toml").read_text(encoding="utf-8")

            assert manifest.mode == "nvt"
            assert entries["restart.coordinate"]["status"] == "converted"
            assert entries["restart.velocity"]["status"] == "converted"
            assert entries["restart.box"]["status"] == "converted"
            assert entries["output.legacy_sidecar.crd"]["status"] == "legacy_output_sidecar_preserved"
            assert entries["output.legacy_sidecar.box"]["status"] == "legacy_output_sidecar_preserved"
            assert entries["output.legacy_sidecar.vel"]["status"] == "legacy_output_sidecar_preserved"
            assert "trajectory.crd" not in entries
            assert "trajectory.box" not in entries
            assert "trajectory.vel" not in entries
            assert 'input_h5_restart_path = "restart.spgr.h5"' in bundled_mdin
            assert 'input_h5_restart_load = "full"' in bundled_mdin
            assert "input_h5_trajectory_path" not in bundled_mdin
            assert 'crd = "traj.dat"' in bundled_mdin
            assert 'box = "traj_box.dat"' in bundled_mdin
            assert 'vel = "traj_vel.dat"' in bundled_mdin
            assert not (output_dir / "bundle" / "trajectory.spg.h5md").exists()
    finally:
        (
            converter_module.ensure_dataset,
            converter_module.ensure_group,
            converter_module.ensure_hard_link,
            converter_module.set_attrs,
            converter_module.write_string,
            converter_module.write_bytes,
            converter_module.write_string_array,
        ) = originals


def test_legacy_to_bundle_bundled_mdin_sidecars_and_overrides():
    writes = []

    def fake_ensure_dataset(file_path, dataset_path, data):
        writes.append(("dataset", Path(file_path).name, dataset_path, np.asarray(data).shape))

    def fake_ensure_group(file_path, group_path):
        writes.append(("group", Path(file_path).name, group_path))

    def fake_ensure_hard_link(file_path, target_path, link_path):
        writes.append(("hard_link", Path(file_path).name, target_path, link_path))

    def fake_set_attrs(file_path, object_path, attrs):
        writes.append(("attrs", Path(file_path).name, object_path, dict(attrs)))

    def fake_write_string(file_path, dataset_path, value):
        writes.append(("string", Path(file_path).name, dataset_path, value))

    def fake_write_bytes(file_path, dataset_path, value):
        writes.append(("bytes", Path(file_path).name, dataset_path, bytes(value)))

    def fake_write_string_array(file_path, dataset_path, values):
        writes.append(("string_array", Path(file_path).name, dataset_path, tuple(values)))

    originals = (
        converter_module.ensure_dataset,
        converter_module.ensure_group,
        converter_module.ensure_hard_link,
        converter_module.set_attrs,
        converter_module.write_string,
        converter_module.write_bytes,
        converter_module.write_string_array,
    )
    converter_module.ensure_dataset = fake_ensure_dataset
    converter_module.ensure_group = fake_ensure_group
    converter_module.ensure_hard_link = fake_ensure_hard_link
    converter_module.set_attrs = fake_set_attrs
    converter_module.write_string = fake_write_string
    converter_module.write_bytes = fake_write_bytes
    converter_module.write_string_array = fake_write_string_array
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "case"
            output_dir = Path(tmpdir) / "out"
            case_dir.mkdir()
            write_basic_case(case_dir)
            mdin_path = case_dir / "mdin.spg.toml"
            mdin_path.write_text(
                'input_h5_topology_path = "stale_topology.spgt.h5"\n'
                'input_h5_protocol_path = "stale_protocol.spgp.h5"\n'
                'input_h5_restart_path = "stale_restart.spgr.h5"\n'
                'input_h5_restart_load = "custom"\n'
                'input_h5_trajectory_path = "stale_trajectory.spg.h5md"\n'
                'input_h5_trajectory_particle_stream = "stale_stream"\n'
                + mdin_path.read_text(encoding="utf-8"),
                encoding="utf-8",
            )

            manifest = converter_module.convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
            entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
            bundled_mdin = (output_dir / "bundle" / "mdin.bundled.spg.toml").read_text(encoding="utf-8")

            assert entries["run_mdin.input_h5_topology_path"]["status"] == "legacy_input_override_replaced"
            assert entries["run_mdin.input_h5_protocol_path"]["status"] == "legacy_input_override_replaced"
            assert entries["run_mdin.input_h5_restart_path"]["status"] == "legacy_input_override_replaced"
            assert entries["run_mdin.input_h5_restart_load"]["status"] == "legacy_input_override_replaced"
            assert entries["run_mdin.input_h5_trajectory_path"]["status"] == "legacy_input_override_replaced"
            assert entries["run_mdin.input_h5_trajectory_particle_stream"]["status"] == "legacy_input_override_replaced"
            assert entries["restart.protocol_sidecar.cv_in_file"]["status"] == "sidecar_embedded"
            assert entries["restart.protocol_sidecar.meta_potential_in_file"]["status"] == "sidecar_embedded"
            assert "stale_topology.spgt.h5" not in bundled_mdin
            assert "stale_protocol.spgp.h5" not in bundled_mdin
            assert "stale_restart.spgr.h5" not in bundled_mdin
            assert 'input_h5_restart_load = "custom"' not in bundled_mdin
            assert "stale_trajectory.spg.h5md" not in bundled_mdin
            assert "stale_stream" not in bundled_mdin
            assert 'input_h5_topology_path = "topology.spgt.h5"' in bundled_mdin
            assert 'input_h5_protocol_path = "protocol.spgp.h5"' in bundled_mdin
            assert 'input_h5_restart_path = "restart.spgr.h5"' in bundled_mdin
            assert 'input_h5_restart_load = "full"' in bundled_mdin
            assert 'input_h5_trajectory_path = "trajectory.spg.h5md"' in bundled_mdin
            assert 'input_h5_trajectory_particle_stream = "all"' in bundled_mdin
            assert "mass_in_file" not in bundled_mdin
            assert "cmap_in_file" not in bundled_mdin
            assert "improper_dihedral_in_file" not in bundled_mdin
            assert "LJ_soft_core_in_file" not in bundled_mdin
            assert "subsys_division_in_file" not in bundled_mdin
            assert "EAM_in_file" not in bundled_mdin
            assert "EAM_atom_type_in_file" not in bundled_mdin
            assert "[EAM]" not in bundled_mdin
            assert "SW_in_file" not in bundled_mdin
            assert "EDIP_in_file" not in bundled_mdin
            assert "pairwise_force_in_file" not in bundled_mdin
            assert "listed_forces_in_file" not in bundled_mdin
            assert 'custom_pair_in_file = "custom_pair.txt"' not in bundled_mdin
            assert 'custom_bond_in_file = "custom_bond.txt"' not in bundled_mdin
            assert 'custom_pair_in_file = "legacy_sidecars/custom_pair_in_file/custom_pair.txt"' in bundled_mdin
            assert 'custom_bond_in_file = "legacy_sidecars/custom_bond_in_file/custom_bond.txt"' in bundled_mdin
            assert "TERSOFF_in_file" not in bundled_mdin
            assert "[TERSOFF]" not in bundled_mdin
            assert "qc_type_in_file" not in bundled_mdin
            assert "REAXFF_in_file" not in bundled_mdin
            assert "REAXFF_type_in_file" not in bundled_mdin
            assert "[REAXFF]" not in bundled_mdin
            assert "SITS_nk_in_file" not in bundled_mdin
            assert "restrain_in_file" not in bundled_mdin
            assert "restrain_cv_in_file" not in bundled_mdin
            assert "steer_cv_in_file" not in bundled_mdin
            assert "nose_hoover_chain_restart_input" not in bundled_mdin
            assert "SITS_in_file" not in bundled_mdin
            assert "SITS_atom_in_file" not in bundled_mdin
            assert "meta_edge_in_file" not in bundled_mdin
            assert "meta_potential_in_file" not in bundled_mdin
            assert "meta_scatter_in_file" not in bundled_mdin
            assert "hills_in_file" not in bundled_mdin
            assert "constrain_in_file" not in bundled_mdin
            assert "soft_walls_in_file" not in bundled_mdin
            assert "restrain_atom_id" not in bundled_mdin
            assert 'crd = "traj.dat"' not in bundled_mdin
            assert 'box = "traj_box.dat"' not in bundled_mdin
            assert 'vel = "traj_vel.dat"' not in bundled_mdin
            assert 'mdout = "mdout.txt"' in bundled_mdin
            assert 'output_h5_trajectory_path = "prod.spg.h5md"' in bundled_mdin

            bundle_dir = output_dir / "bundle"
            assert (bundle_dir / "legacy_sidecars" / "cv_in_file" / "cv.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "restrain_in_file" / "restrain.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "restrain_cv_in_file" / "restrain_cv.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "steer_cv_in_file" / "steer_cv.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "mass_in_file" / "mass.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "cmap_in_file" / "cmap.txt").exists()
            assert not (
                bundle_dir / "legacy_sidecars" / "improper_dihedral_in_file" / "improper.txt"
            ).exists()
            assert (bundle_dir / "legacy_sidecars" / "LJ_soft_core_in_file" / "lj_soft_core.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "EAM_in_file" / "eam.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "EAM_atom_type_in_file" / "eam_atom_type.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "SW_in_file" / "sw.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "EDIP_in_file" / "edip.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "pairwise_force_in_file" / "pairwise_force.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "custom_pair_in_file" / "custom_pair.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "listed_forces_in_file" / "listed_forces.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "custom_bond_in_file" / "custom_bond.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "TERSOFF_in_file" / "tersoff.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "qc_type_in_file" / "qc_type.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "REAXFF_in_file" / "reaxff.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "REAXFF_type_in_file" / "reaxff_type.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "SITS_nk_in_file" / "sits_nk.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "SITS_in_file" / "sits.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "meta_edge_in_file" / "meta_edge.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "meta_potential_in_file" / "meta_potential.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "meta_scatter_in_file" / "meta_scatter.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "hills_in_file" / "hills.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "nose_hoover_chain_restart_input" / "nhc_restart.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "SITS_atom_in_file" / "sits_atom.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "constrain_in_file" / "constrain.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "soft_walls_in_file" / "soft_walls.txt").exists()
            assert (bundle_dir / "legacy_sidecars" / "restrain_atom_id" / "restrain_atom_id.txt").exists()
            assert not (bundle_dir / "legacy_sidecars" / "crd" / "traj.dat").exists()
            assert not (bundle_dir / "legacy_sidecars" / "box" / "traj_box.dat").exists()
            assert not (bundle_dir / "legacy_sidecars" / "vel" / "traj_vel.dat").exists()

            string_arrays = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "string_array"
            }
            strings = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "string"
            }
            datasets = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "dataset"
            }
            hard_links = {
                (item[1], item[3]): item[2]
                for item in writes
                if item[0] == "hard_link"
            }
            attrs = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "attrs"
            }
            assert strings[("topology.spgt.h5", "/schema/name")] == "sponge.topology.h5"
            assert strings[("topology.spgt.h5", "/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("topology.spgt.h5", "/parameters/sponge/schema/name")] == "sponge.topology.h5"
            assert strings[("topology.spgt.h5", "/parameters/sponge/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("topology.spgt.h5", "/topology/atom_order_hash")].startswith("sha256:")
            assert strings[("topology.spgt.h5", "/topology/topology_hash")].startswith("sha256:")
            assert strings[("topology.spgt.h5", "/topology/forcefield_hash")].startswith("sha256:")
            assert datasets[("topology.spgt.h5", "/topology/atom_count")] == ()
            assert datasets[("topology.spgt.h5", "/forcefield/improper/pk")] == (1,)
            assert ("topology.spgt.h5", "/forcefield/improper/k") not in datasets
            assert datasets[
                ("topology.spgt.h5", "/forcefield/custom_force/pairwise/data/custom_pair/parameter/value")
            ] == (2, 1)
            assert datasets[
                ("topology.spgt.h5", "/forcefield/custom_force/listed/data/custom_bond/parameter/value")
            ] == (1, 4)
            assert strings[("protocol.spgp.h5", "/protocol/topology_compatibility/topology_hash")] == strings[
                ("topology.spgt.h5", "/topology/topology_hash")
            ]
            assert strings[("protocol.spgp.h5", "/schema/name")] == "sponge.protocol.h5"
            assert strings[("protocol.spgp.h5", "/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("protocol.spgp.h5", "/parameters/sponge/schema/name")] == "sponge.protocol.h5"
            assert strings[("protocol.spgp.h5", "/parameters/sponge/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("protocol.spgp.h5", "/identity/content_hash")].startswith("sha256:")
            assert datasets[("protocol.spgp.h5", "/protocol/cv_count")] == ()
            assert datasets[("protocol.spgp.h5", "/protocol/restraint_count")] == ()
            assert strings[("protocol.spgp.h5", "/protocol/enhanced_sampling/method")] == (
                "SITS,metadynamics,steer,soft_walls"
            )
            assert strings[("restart.spgr.h5", "/parameters/sponge/schema/name")] == "sponge.restart.h5"
            assert strings[("restart.spgr.h5", "/parameters/sponge/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("restart.spgr.h5", "/parameters/sponge/output/status")] == "finalized"
            assert datasets[("restart.spgr.h5", "/parameters/sponge/output/frame_count")] == (1,)
            assert datasets[("restart.spgr.h5", "/parameters/sponge/output/last_complete_step")] == (1,)
            assert datasets[("restart.spgr.h5", "/parameters/sponge/output/last_complete_time")] == (1,)
            assert attrs[("restart.spgr.h5", "/h5md")]["version"].tolist() == [1, 1]
            assert attrs[("restart.spgr.h5", "/particles/all/position/value")]["unit"] == "Angstrom"
            assert attrs[("restart.spgr.h5", "/particles/all/velocity/value")]["unit"] == "Angstrom ps-1"
            assert attrs[("restart.spgr.h5", "/particles/all/box/edges/value")]["unit"] == "Angstrom"
            assert hard_links[("restart.spgr.h5", "/particles/all/position/step")] == "/particles/all/step"
            assert hard_links[("restart.spgr.h5", "/particles/all/position/time")] == "/particles/all/time"
            assert hard_links[("restart.spgr.h5", "/particles/all/velocity/step")] == "/particles/all/step"
            assert hard_links[("restart.spgr.h5", "/particles/all/box/edges/step")] == "/particles/all/step"
            assert strings[("restart.spgr.h5", "/parameters/restart/protocol_sidecars/cv_in_file")].startswith("print")
            assert strings[("restart.spgr.h5", "/parameters/restart/protocol_sidecars/meta_potential_in_file")].startswith(
                "Meta potential"
            )
            assert strings[("trajectory.spg.h5md", "/parameters/sponge/schema/name")] == "sponge.output.h5md"
            assert strings[("trajectory.spg.h5md", "/parameters/sponge/schema/version")] == "xponge.legacy_to_bundle.v1"
            assert strings[("trajectory.spg.h5md", "/parameters/sponge/output/mode")] == "single"
            assert strings[("trajectory.spg.h5md", "/parameters/sponge/output/status")] == "finalized"
            assert datasets[("trajectory.spg.h5md", "/parameters/sponge/output/frame_count")] == (1,)
            assert datasets[("trajectory.spg.h5md", "/parameters/sponge/output/last_complete_step")] == (1,)
            assert datasets[("trajectory.spg.h5md", "/parameters/sponge/output/last_complete_time")] == (1,)
            assert string_arrays[("trajectory.spg.h5md", "/parameters/sponge/output/particle_streams")] == ("all",)
            assert attrs[("trajectory.spg.h5md", "/h5md")]["version"].tolist() == [1, 1]
            assert attrs[("trajectory.spg.h5md", "/particles/all/position/value")]["unit"] == "Angstrom"
            assert attrs[("trajectory.spg.h5md", "/particles/all/velocity/value")]["unit"] == "Angstrom ps-1"
            assert attrs[("trajectory.spg.h5md", "/particles/all/box/edges/value")]["unit"] == "Angstrom"
            assert hard_links[("trajectory.spg.h5md", "/particles/all/position/step")] == "/particles/all/step"
            assert hard_links[("trajectory.spg.h5md", "/particles/all/position/time")] == "/particles/all/time"
            assert hard_links[("trajectory.spg.h5md", "/particles/all/velocity/step")] == "/particles/all/step"
            assert hard_links[("trajectory.spg.h5md", "/particles/all/box/edges/step")] == "/particles/all/step"
            assert "mass_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "cmap_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "improper_dihedral_in_file" not in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "LJ_soft_core_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "EAM_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "EAM_atom_type_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "SW_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "EDIP_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "pairwise_force_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "listed_forces_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "custom_pair_in_file" not in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "custom_bond_in_file" not in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "TERSOFF_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "qc_type_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "REAXFF_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "REAXFF_type_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "cv_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "restrain_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "restrain_cv_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "steer_cv_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "constrain_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "SITS_atom_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "SITS_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "meta_edge_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "soft_walls_in_file" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "restrain_atom_id" in string_arrays[
                ("protocol.spgp.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "SITS_nk_in_file" in string_arrays[
                ("restart.spgr.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "nose_hoover_chain_restart_input" in string_arrays[
                ("restart.spgr.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "meta_potential_in_file" in string_arrays[
                ("restart.spgr.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "meta_scatter_in_file" in string_arrays[
                ("restart.spgr.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert "hills_in_file" in string_arrays[
                ("restart.spgr.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
    finally:
        (
            converter_module.ensure_dataset,
            converter_module.ensure_group,
            converter_module.ensure_hard_link,
            converter_module.set_attrs,
            converter_module.write_string,
            converter_module.write_bytes,
            converter_module.write_string_array,
        ) = originals


def test_legacy_to_bundle_static_default_prefix_generates_required_protocol():
    writes = []

    def fake_ensure_dataset(file_path, dataset_path, data):
        writes.append(("dataset", Path(file_path).name, dataset_path, np.asarray(data).shape))

    def fake_write_string(file_path, dataset_path, value):
        writes.append(("string", Path(file_path).name, dataset_path, value))

    def fake_write_string_array(file_path, dataset_path, values):
        writes.append(("string_array", Path(file_path).name, dataset_path, tuple(values)))

    originals = (
        converter_module.ensure_dataset,
        converter_module.ensure_group,
        converter_module.ensure_hard_link,
        converter_module.set_attrs,
        converter_module.write_string,
        converter_module.write_bytes,
        converter_module.write_string_array,
    )
    converter_module.ensure_dataset = fake_ensure_dataset
    converter_module.ensure_group = lambda *args, **kwargs: None
    converter_module.ensure_hard_link = lambda *args, **kwargs: None
    converter_module.set_attrs = lambda *args, **kwargs: None
    converter_module.write_string = fake_write_string
    converter_module.write_bytes = lambda *args, **kwargs: None
    converter_module.write_string_array = fake_write_string_array
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "case"
            output_dir = Path(tmpdir) / "out"
            case_dir.mkdir()
            (case_dir / "legacy_coordinate.txt").write_text(
                "2\n"
                "1.0 2.0 3.0\n"
                "4.0 5.0 6.0\n"
                "10.0 20.0 30.0 90.0 90.0 90.0\n",
                encoding="utf-8",
            )
            (case_dir / "legacy_velocity.txt").write_text(
                "2\n"
                "0.1 0.2 0.3\n"
                "0.4 0.5 0.6\n",
                encoding="utf-8",
            )
            (case_dir / "legacy_mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
            (case_dir / "legacy_charge.txt").write_text("2\n0.1\n-0.1\n", encoding="utf-8")
            (case_dir / "mdin.spg.toml").write_text(
                'mode = "nvt"\n'
                'default_in_file_prefix = "legacy"\n'
                "step_limit = 1\n",
                encoding="utf-8",
            )

            manifest = converter_module.convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
            entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
            bundled_mdin = (output_dir / "bundle" / "mdin.bundled.spg.toml").read_text(encoding="utf-8")
            strings = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "string"
            }
            string_arrays = {
                (item[1], item[2]): item[3]
                for item in writes
                if item[0] == "string_array"
            }

            assert entries["restart.coordinate"]["source_path"].endswith("legacy_coordinate.txt")
            assert entries["restart.velocity"]["source_path"].endswith("legacy_velocity.txt")
            assert entries["topology.mass"]["source_path"].endswith("legacy_mass.txt")
            assert entries["topology.charge"]["source_path"].endswith("legacy_charge.txt")
            assert entries["protocol.required_empty"]["status"] == "generated"
            assert entries["run_mdin.default_in_file_prefix"]["status"] == "legacy_input_override_replaced"
            assert 'default_in_file_prefix = "legacy"' not in bundled_mdin
            assert 'input_h5_topology_path = "topology.spgt.h5"' in bundled_mdin
            assert 'input_h5_protocol_path = "protocol.spgp.h5"' in bundled_mdin
            assert 'input_h5_restart_path = "restart.spgr.h5"' in bundled_mdin
            assert "mass_in_file" in string_arrays[
                ("topology.spgt.h5", "/parameters/sponge/files/legacy_sidecars/key")
            ]
            assert strings[("protocol.spgp.h5", "/protocol/topology_compatibility/topology_hash")] == strings[
                ("topology.spgt.h5", "/topology/topology_hash")
            ]
    finally:
        (
            converter_module.ensure_dataset,
            converter_module.ensure_group,
            converter_module.ensure_hard_link,
            converter_module.set_attrs,
            converter_module.write_string,
            converter_module.write_bytes,
            converter_module.write_string_array,
        ) = originals


def test_legacy_to_bundle_dynamic_custom_force_default_prefix():
    writes = []

    def fake_ensure_dataset(file_path, dataset_path, data):
        writes.append((Path(file_path).name, dataset_path, np.asarray(data).shape))

    originals = (
        converter_module.ensure_dataset,
        converter_module.ensure_group,
        converter_module.ensure_hard_link,
        converter_module.set_attrs,
        converter_module.write_string,
        converter_module.write_bytes,
        converter_module.write_string_array,
    )
    converter_module.ensure_dataset = fake_ensure_dataset
    converter_module.ensure_group = lambda *args, **kwargs: None
    converter_module.ensure_hard_link = lambda *args, **kwargs: None
    converter_module.set_attrs = lambda *args, **kwargs: None
    converter_module.write_string = lambda *args, **kwargs: None
    converter_module.write_bytes = lambda *args, **kwargs: None
    converter_module.write_string_array = lambda *args, **kwargs: None
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            case_dir = Path(tmpdir) / "case"
            output_dir = Path(tmpdir) / "out"
            case_dir.mkdir()
            write_basic_case(case_dir)
            (case_dir / "legacy_custom_pair.txt").write_text(
                (case_dir / "custom_pair.txt").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            (case_dir / "legacy_custom_bond.txt").write_text(
                (case_dir / "custom_bond.txt").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
            mdin_path = case_dir / "mdin.spg.toml"
            mdin_text = mdin_path.read_text(encoding="utf-8")
            mdin_text = mdin_text.replace('custom_pair_in_file = "custom_pair.txt"\n', "")
            mdin_text = mdin_text.replace('custom_bond_in_file = "custom_bond.txt"\n', "")
            mdin_path.write_text('default_in_file_prefix = "legacy"\n' + mdin_text, encoding="utf-8")

            manifest = converter_module.convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
            entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
            bundled_mdin = (output_dir / "bundle" / "mdin.bundled.spg.toml").read_text(encoding="utf-8")

            assert entries["topology.pairwise_force_data.custom_pair"]["source_path"].endswith(
                "legacy_custom_pair.txt"
            )
            assert entries["topology.listed_force_data.custom_bond"]["source_path"].endswith(
                "legacy_custom_bond.txt"
            )
            assert 'default_in_file_prefix = "legacy"' not in bundled_mdin
            assert 'custom_pair_in_file = "legacy_sidecars/custom_pair_in_file/legacy_custom_pair.txt"' in bundled_mdin
            assert 'custom_bond_in_file = "legacy_sidecars/custom_bond_in_file/legacy_custom_bond.txt"' in bundled_mdin
            assert (output_dir / "bundle" / "legacy_sidecars" / "custom_pair_in_file" / "legacy_custom_pair.txt").exists()
            assert (output_dir / "bundle" / "legacy_sidecars" / "custom_bond_in_file" / "legacy_custom_bond.txt").exists()
            assert ("topology.spgt.h5", "/forcefield/custom_force/pairwise/data/custom_pair/parameter/value", (2, 1)) in writes
            assert ("topology.spgt.h5", "/forcefield/custom_force/listed/data/custom_bond/parameter/value", (1, 4)) in writes
    finally:
        (
            converter_module.ensure_dataset,
            converter_module.ensure_group,
            converter_module.ensure_hard_link,
            converter_module.set_attrs,
            converter_module.write_string,
            converter_module.write_bytes,
            converter_module.write_string_array,
        ) = originals


def test_legacy_to_bundle_writes_typed_topology():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        write_basic_case(case_dir)

        manifest = convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["topology.mass"]["status"] == "typed_converted"
        assert entries["topology.cmap"]["status"] == "typed_converted"
        assert entries["topology.improper_dihedral"]["status"] == "typed_converted"
        assert entries["topology.LJ_soft_core"]["status"] == "typed_converted"
        assert entries["restart.SITS_nk"]["status"] == "typed_converted"
        assert entries["protocol.cv"]["status"] == "typed_converted"
        assert entries["protocol.restrain"]["status"] == "typed_converted"
        assert entries["protocol.restrain_cv"]["status"] == "typed_converted"
        assert entries["protocol.steer_cv"]["status"] == "typed_converted"
        assert entries["protocol.SITS"]["status"] == "typed_converted"
        assert entries["protocol.meta_edge"]["status"] == "typed_converted"
        assert entries["restart.meta_potential"]["status"] == "typed_converted"
        assert entries["restart.meta_scatter"]["status"] == "typed_converted"
        assert entries["restart.hills"]["status"] == "typed_converted"
        assert entries["restart.nose_hoover_chain_restart_input"]["status"] == "typed_converted"
        assert entries["protocol.restrain_atom_id"]["status"] == "typed_converted"
        assert entries["protocol.constrain"]["status"] == "typed_converted"
        assert entries["protocol.soft_walls"]["status"] == "typed_converted"
        assert entries["topology.SW"]["status"] == "typed_converted"
        assert entries["topology.EDIP"]["status"] == "typed_converted"
        assert entries["topology.pairwise_force"]["status"] == "typed_converted"
        assert entries["topology.listed_forces"]["status"] == "typed_converted"
        assert entries["topology.TERSOFF"]["status"] == "typed_converted"
        assert entries["topology.qc_type"]["status"] == "typed_converted"
        assert entries["topology.EAM"]["status"] == "typed_converted"
        assert entries["topology.EAM_atom_type"]["status"] == "typed_converted"
        assert entries["topology.REAXFF"]["status"] == "typed_converted"
        assert entries["topology.REAXFF_type"]["status"] == "typed_converted"
        assert entries["trajectory.crd"]["status"] == "typed_converted"

        with h5py.File(output_dir / "bundle" / "topology.spgt.h5", "r") as handle:
            assert _h5_string(handle["/schema/name"][()]) == "sponge.topology.h5"
            assert _h5_string(handle["/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.topology.h5"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/topology/atom_order_hash"][()]).startswith("sha256:")
            assert _h5_string(handle["/topology/topology_hash"][()]).startswith("sha256:")
            assert _h5_string(handle["/topology/forcefield_hash"][()]).startswith("sha256:")
            assert np.allclose(handle["/atoms/mass"][...], np.asarray([1.0, 16.0], dtype=np.float32))
            assert np.allclose(handle["/atoms/charge"][...], np.asarray([0.1, -0.1], dtype=np.float32))
            assert np.array_equal(handle["/residues/atom_offset"][...], np.asarray([0, 2], dtype=np.int64))
            assert np.array_equal(handle["/atoms/residue_index"][...], np.asarray([0, 0], dtype=np.int32))
            assert np.array_equal(handle["/forcefield/bond/atoms"][...], np.asarray([[0, 1]], dtype=np.int32))
            assert np.allclose(handle["/forcefield/bond/k"][...], np.asarray([100.0], dtype=np.float32))
            assert np.array_equal(
                handle["/forcefield/improper/atoms"][...],
                np.asarray([[0, 1, 0, 1]], dtype=np.int32),
            )
            assert np.allclose(handle["/forcefield/improper/pk"][...], np.asarray([2.5], dtype=np.float32))
            assert "/forcefield/improper/k" not in handle
            sidecar_keys = _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"])
            assert "improper_dihedral_in_file" not in sidecar_keys
            assert np.array_equal(handle["/forcefield/cmap/atoms"][...], np.asarray([[0, 1, 0, 1, 0]], dtype=np.int32))
            assert np.allclose(handle["/forcefield/gb/params"][...], np.asarray([[1.5, 0.8], [1.2, 0.85]], dtype=np.float32))
            assert np.array_equal(handle["/forcefield/virtual_atom/type"][...], np.asarray([2], dtype=np.int32))
            assert np.array_equal(handle["/forcefield/lj/type"][...], np.asarray([0, 0], dtype=np.int32))
            assert np.array_equal(handle["/forcefield/lj_soft_core/atom_type_B"][...], np.asarray([0, 0], dtype=np.int32))
            assert np.array_equal(handle["/forcefield/subsys_division"][...], np.asarray([1, 2], dtype=np.int32))
            assert _h5_string(handle["/manybody/eam/format"][()]) == "funcfl"
            assert np.array_equal(handle["/manybody/eam/atomic_number"][...], np.asarray([29], dtype=np.int32))
            assert np.allclose(handle["/manybody/eam/embed/value"][...], np.asarray([[23.060548, 46.121096]], dtype=np.float32))
            assert np.allclose(handle["/manybody/eam/electron_density/value"][...], np.asarray([[5.0, 6.0]], dtype=np.float32))
            assert np.array_equal(handle["/manybody/eam/atom_type"][...], np.asarray([0, 1], dtype=np.int32))
            assert np.array_equal(handle["/manybody/sw/pair/type"][...], np.asarray([[0, 0]], dtype=np.int32))
            assert np.allclose(
                handle["/manybody/sw/triple/parameters"][...],
                np.asarray([[9.0, 10.0, 11.0]], dtype=np.float32),
            )
            assert np.array_equal(handle["/manybody/edip/atom_type"][...], np.asarray([0, 0], dtype=np.int32))
            assert np.allclose(
                handle["/manybody/edip/triple/parameters"][...],
                np.asarray([[9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]], dtype=np.float32),
            )
            assert _h5_string(handle["/forcefield/custom_force/pairwise/name"][()]) == "custom_pair"
            assert _h5_string_list(handle["/forcefield/custom_force/pairwise/parameters/name"]) == [
                "epsilon_ij",
                "sigma_ij",
            ]
            assert handle["/forcefield/custom_force/pairwise/parameters/ij_count"][()] == 2
            assert handle["/forcefield/custom_force/pairwise/with_ele"][()] == np.bool_(False)
            assert np.allclose(
                handle["/forcefield/custom_force/pairwise/data/custom_pair/parameter/value"][...],
                np.asarray([[1.0], [2.0]], dtype=np.float32),
            )
            assert np.array_equal(
                handle["/forcefield/custom_force/pairwise/data/custom_pair/atom_type"][...],
                np.asarray([0, 0], dtype=np.int32),
            )
            assert _h5_string_list(handle["/forcefield/custom_force/listed/name"]) == ["custom_bond"]
            assert _h5_string_list(handle["/forcefield/custom_force/listed/parameters/name"]) == [
                "atom_a",
                "atom_b",
                "k",
                "r0",
            ]
            assert np.array_equal(
                handle["/forcefield/custom_force/listed/parameters/is_atom"][...],
                np.asarray([1, 1, 0, 0], dtype=np.int32),
            )
            assert _h5_string_list(handle["/forcefield/custom_force/listed/connected_atoms"]) == ["ab"]
            assert np.allclose(
                handle["/forcefield/custom_force/listed/data/custom_bond/parameter/value"][...],
                np.asarray([[0.0, 1.0, 10.0, 1.5]], dtype=np.float32),
            )
            assert _h5_string_list(handle["/manybody/tersoff/type_name"]) == ["Si"]
            assert np.array_equal(handle["/manybody/tersoff/entry/type"][...], np.asarray([[0, 0, 0]], dtype=np.int32))
            assert np.allclose(handle["/manybody/tersoff/entry/parameters_raw"][0, 13], 1380.0)
            assert np.allclose(handle["/manybody/tersoff/entry/parameters"][0, 13], 1380.0 * 23.0605480)
            assert np.array_equal(handle["/manybody/tersoff/map"][...], np.asarray([0], dtype=np.int32))
            assert np.array_equal(handle["/manybody/tersoff/atom_type"][...], np.asarray([0, 0], dtype=np.int32))
            assert np.array_equal(handle["/qc/type/atom_index"][...], np.asarray([0, 1], dtype=np.int32))
            assert _h5_string_list(handle["/qc/type/symbol"]) == ["H", "O"]
            assert handle["/qc/type/charge"][()] == 0
            assert handle["/qc/type/multiplicity"][()] == 1
            assert handle["/manybody/reaxff/parameters/general/count"][()] == 2
            assert np.allclose(
                handle["/manybody/reaxff/parameters/general/value"][...],
                np.asarray([50.0, 9.5], dtype=np.float32),
            )
            assert _h5_string_list(handle["/manybody/reaxff/parameters/atom/type_name"]) == ["O", "H"]
            assert np.array_equal(
                handle["/manybody/reaxff/parameters/atom/value_offset"][...],
                np.asarray([0, 32, 64], dtype=np.int64),
            )
            assert np.array_equal(handle["/manybody/reaxff/parameters/bond/type"][...], np.asarray([[1, 2]], dtype=np.int32))
            assert np.array_equal(handle["/manybody/reaxff/parameters/torsion/type"][...], np.asarray([[0, 1, 2, 0]], dtype=np.int32))
            assert np.allclose(
                handle["/manybody/reaxff/parameters/hydrogen_bond/value"][...],
                np.asarray([[1.9, -4.4, 1.7, 3.0]], dtype=np.float32),
            )
            assert _h5_string_list(handle["/manybody/reaxff/type/name"]) == ["O", "H"]
            assert handle["/manybody/reaxff/type/count"][()] == 2

        with h5py.File(output_dir / "bundle" / "protocol.spgp.h5", "r") as handle:
            assert _h5_string(handle["/schema/name"][()]) == "sponge.protocol.h5"
            assert _h5_string(handle["/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.protocol.h5"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert handle["/protocol/cv_count"][()] > 0
            assert handle["/protocol/restraint_count"][()] > 0
            assert _h5_string(handle["/protocol/enhanced_sampling/method"][()]) == (
                "SITS,metadynamics,steer,soft_walls"
            )
            assert _h5_string_list(handle["/cv/config/section/name"]) == ["print", "distance"]
            assert _h5_string_list(handle["/cv/config/key"]) == ["CV", "CV_type", "atom"]
            assert _h5_string_list(handle["/restraint/config/section/name"]) == ["print", "distance"]
            assert _h5_string_list(handle["/restraint/cv/config/key"]) == ["CV", "weight", "reference"]
            assert _h5_string_list(handle["/steer/config/value"]) == ["distance", "2.0"]
            assert _h5_string_list(handle["/sits/config/key"]) == ["CV", "nk"]
            assert np.array_equal(
                handle["/meta/default/grid/coordinate"][...],
                np.asarray([[0.0, 0.0], [1.0, 0.0]], dtype=np.float32),
            )
            assert np.allclose(
                handle["/meta/default/grid/force"][...],
                np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
            )
            assert np.array_equal(handle["/restraint/default/atom_indices"][...], np.asarray([0, 1], dtype=np.int32))
            assert np.allclose(
                handle["/restraint/default/weight"][...],
                np.asarray([[10.0, 0.0, 0.0], [0.0, 20.0, 0.0]], dtype=np.float32),
            )
            assert np.array_equal(handle["/constraint/default/pairs/atoms"][...], np.asarray([[0, 1]], dtype=np.int32))
            assert np.allclose(handle["/constraint/default/pairs/r0"][...], np.asarray([1.5], dtype=np.float32))
            assert np.array_equal(handle["/sits/atom_indices"][...], np.asarray([0, 1], dtype=np.int32))
            assert _h5_string_list(handle["/wall/soft/name"]) == ["z_wall"]
            assert _h5_string_list(handle["/wall/soft/potential"]) == [
                "E = powf((z - 5.0f) * (z - 5.0f), -6.0f);"
            ]
            assert _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"])

        with h5py.File(output_dir / "bundle" / "restart.spgr.h5", "r") as handle:
            assert np.allclose(handle["/parameters/restart/bias/sits/SITS/nk"][...], np.asarray([1.0, 2.0], dtype=np.float32))
            assert np.allclose(
                handle["/parameters/restart/references/restraint/default/coordinate"][...],
                np.asarray([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32),
            )
            assert np.allclose(
                handle["/parameters/restart/thermostat/nose_hoover_chain"][...],
                np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
            )
            assert np.array_equal(
                handle["/parameters/restart/bias/meta/default/potential/grid"][...],
                np.asarray([2, 2], dtype=np.int32),
            )
            assert np.allclose(
                handle["/parameters/restart/bias/meta/default/potential/value"][...],
                np.asarray([3.0, 4.0, 5.0, 6.0], dtype=np.float32),
            )
            assert _h5_string(handle["/parameters/restart/bias/meta/default/potential_export"][()]).startswith(
                "Meta potential"
            )
            assert np.array_equal(
                handle["/parameters/restart/bias/meta/default/scatter/grid"][...],
                np.asarray([2], dtype=np.int32),
            )
            assert np.allclose(
                handle["/parameters/restart/bias/meta/default/hills_typed/value"][...],
                np.asarray([[0.1, 0.2, 1.5], [0.3, 0.4, 1.6]], dtype=np.float32),
            )
            assert _h5_string(handle["/parameters/restart/bias/meta/default/hills"][()]) == "0.1 0.2 1.5\n0.3 0.4 1.6\n"
            assert _h5_string(handle["/parameters/restart/protocol_sidecars/cv_in_file"][()]).startswith("print")
            assert _h5_string(handle["/parameters/restart/protocol_sidecars/meta_potential_in_file"][()]).startswith(
                "Meta potential"
            )
            assert np.array_equal(handle["/h5md"].attrs["version"], np.asarray([1, 1], dtype=np.int32))
            assert _h5_string(handle["/h5md/creator"].attrs["name"]) == "XPONGE"
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.restart.h5"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/status"][()]) == "finalized"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([1], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([0], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([0.0]))
            assert np.array_equal(handle["/run/current_step"][...], np.asarray([0], dtype=np.int64))
            assert np.allclose(handle["/run/current_time"][...], np.asarray([0.0]))
            assert _h5_string(handle["/run/state_type"][()]) == "restart"
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert _h5_string(handle["/particles/all/position/value"].attrs["unit"]) == "Angstrom"
            assert _h5_string(handle["/particles/all/velocity/value"].attrs["unit"]) == "Angstrom ps-1"
            assert _h5_string(handle["/particles/all/box/edges/value"].attrs["unit"]) == "Angstrom"
            assert _h5_string(handle["/particles/all/time"].attrs["unit"]) == "ps"
            assert h5py.h5o.get_info(handle["/particles/all/position/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/position/time"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/time"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/velocity/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/box/edges/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr

        with h5py.File(output_dir / "bundle" / "trajectory.spg.h5md", "r") as handle:
            assert np.array_equal(handle["/h5md"].attrs["version"], np.asarray([1, 1], dtype=np.int32))
            assert _h5_string(handle["/h5md/creator"].attrs["name"]) == "XPONGE"
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.output.h5md"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "single"
            assert _h5_string(handle["/parameters/sponge/output/status"][()]) == "finalized"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([1.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert _h5_string(handle["/particles/all/position/value"].attrs["unit"]) == "Angstrom"
            assert _h5_string(handle["/particles/all/velocity/value"].attrs["unit"]) == "Angstrom ps-1"
            assert _h5_string(handle["/particles/all/box/edges/value"].attrs["unit"]) == "Angstrom"
            assert _h5_string(handle["/particles/all/time"].attrs["unit"]) == "ps"
            assert h5py.h5o.get_info(handle["/particles/all/position/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/position/time"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/time"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/velocity/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/particles/all/box/edges/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert np.allclose(handle["/particles/all/position/value"][1, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))
            assert np.array_equal(handle["/particles/all/step"][...], np.asarray([0, 1], dtype=np.int64))
            assert np.allclose(handle["/particles/all/box/edges/value"][1], np.diag([11.0, 21.0, 31.0]).astype(np.float32))
            assert np.allclose(handle["/particles/all/velocity/value"][0, 0], np.asarray([0.1, 0.2, 0.3], dtype=np.float32))


def test_legacy_to_bundle_improper_alias_is_typed_only():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        write_basic_case(case_dir)
        mdin_path = case_dir / "mdin.spg.toml"
        mdin_path.write_text(
            mdin_path.read_text(encoding="utf-8").replace(
                "improper_dihedral_in_file", "improper_in_file"
            ),
            encoding="utf-8",
        )

        manifest = convert_legacy_to_bundle(case_dir, output_dir)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["topology.improper"]["status"] == "typed_converted"

        bundle_dir = output_dir / "bundle"
        bundled_mdin = (bundle_dir / "mdin.bundled.spg.toml").read_text(encoding="utf-8")
        assert "improper_in_file" not in bundled_mdin
        assert not (bundle_dir / "legacy_sidecars" / "improper_in_file" / "improper.txt").exists()
        with h5py.File(bundle_dir / "topology.spgt.h5", "r") as handle:
            assert np.allclose(handle["/forcefield/improper/pk"][...], np.asarray([2.5], dtype=np.float32))
            sidecar_keys = _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"])
            assert "improper_in_file" not in sidecar_keys


def test_legacy_to_bundle_outputs_pass_sponge_h5_probes():
    try:
        import h5py
    except ImportError:
        return

    h5py  # Keep the import local while making the dependency explicit.
    probe_dir = Path("/home/youmans/sidereus/SPONGE/build-dev-cuda13/tests/h5_bundle")
    restart_probe = probe_dir / "h5_restart_read_probe"
    trajectory_probe = probe_dir / "h5_trajectory_read_probe"
    if not (
        restart_probe.is_file()
        and restart_probe.stat().st_mode & 0o111
        and trajectory_probe.is_file()
        and trajectory_probe.stat().st_mode & 0o111
    ):
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        write_basic_case(case_dir)

        convert_legacy_to_bundle(case_dir, output_dir, dry_run=False)
        bundle_dir = output_dir / "bundle"
        restart_result = subprocess.run(
            [str(restart_probe), str(bundle_dir / "restart.spgr.h5"), "2"],
            check=False,
            capture_output=True,
            text=True,
        )
        assert restart_result.returncode == 0, (
            restart_result.stdout,
            restart_result.stderr,
        )
        trajectory_result = subprocess.run(
            [str(trajectory_probe), str(bundle_dir / "trajectory.spg.h5md"), "all", "2"],
            check=False,
            capture_output=True,
            text=True,
        )
        assert trajectory_result.returncode == 0, (
            trajectory_result.stdout,
            trajectory_result.stderr,
        )


def test_legacy_io_contract_matrix_dry_run():
    for mode in ("normal", "rerun"):
        for legacy_input in (False, True):
            for legacy_output in (False, True):
                for vds in (False, True):
                    with tempfile.TemporaryDirectory() as tmpdir:
                        case_dir = Path(tmpdir) / "case"
                        output_dir = Path(tmpdir) / "out"
                        case_dir.mkdir()
                        mdin_lines = [f'mode = "{mode}"', "step_limit = 1"]
                        if legacy_input:
                            (case_dir / "coordinate.txt").write_text(
                                "2\n"
                                "1.0 2.0 3.0\n"
                                "4.0 5.0 6.0\n"
                                "10.0 20.0 30.0 90.0 90.0 90.0\n",
                                encoding="utf-8",
                            )
                            (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
                            mdin_lines.extend(
                                [
                                    'input_h5_topology_path = "old_topology.spgt.h5"',
                                    'input_h5_protocol_path = "old_protocol.spgp.h5"',
                                    'input_h5_restart_path = "old_restart.spgr.h5"',
                                    'coordinate_in_file = "coordinate.txt"',
                                    'mass_in_file = "mass.txt"',
                                ]
                            )
                        else:
                            mdin_lines.extend(
                                [
                                    'input_h5_topology_path = "existing_topology.spgt.h5"',
                                    'input_h5_protocol_path = "existing_protocol.spgp.h5"',
                                    'input_h5_restart_path = "existing_restart.spgr.h5"',
                                ]
                            )
                        if legacy_output:
                            mdin_lines.append('default_out_file_prefix = "prod"')
                            if vds:
                                mdin_lines.append("output_h5_trajectory_vds = true")
                        else:
                            mdin_lines.append('output_h5_trajectory_path = "prod.spg.h5md"')
                            mdin_lines.append(f"output_h5_trajectory_vds = {'true' if vds else 'false'}")
                        (case_dir / "mdin.spg.toml").write_text("\n".join(mdin_lines) + "\n", encoding="utf-8")

                        manifest = convert_legacy_to_bundle(case_dir, output_dir, dry_run=True)
                        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}

                        if legacy_input:
                            assert entries["restart.coordinate"]["status"] == "converted"
                            assert entries["restart.box"]["status"] == "converted"
                            assert entries["topology.mass"]["status"] == "typed_converted"
                            assert entries["protocol.required_empty"]["status"] == "generated"
                            assert entries["run_mdin.input_h5_topology_path"]["status"] == "legacy_input_override_replaced"
                        else:
                            assert entries["run_mdin.input_h5_topology_path"]["status"] == "preserved_in_mdin"
                            assert entries["run_mdin.input_h5_protocol_path"]["status"] == "preserved_in_mdin"
                            assert entries["run_mdin.input_h5_restart_path"]["status"] == "preserved_in_mdin"
                            assert "restart.coordinate" not in entries

                        if legacy_output:
                            assert entries["output.legacy_sidecar.mdout"]["status"] == "legacy_output_sidecar_default"
                            assert entries["output.legacy_sidecar.mdout"]["source_path"] == "prod.out"
                            assert "output.h5.output_h5_trajectory_path" not in entries
                            if vds:
                                assert entries["output.h5.output_h5_trajectory_vds"]["status"] == "output_plan_preserved"
                        else:
                            assert entries["output.h5.output_h5_trajectory_path"]["status"] == "output_plan_preserved"
                            assert entries["output.h5.output_h5_trajectory_vds"]["status"] == "output_plan_preserved"
                            assert "output.legacy_sidecar.mdout" not in entries


def test_legacy_output_mdout_parser():
    with tempfile.TemporaryDirectory() as tmpdir:
        mdout_path = Path(tmpdir) / "prod.out"
        mdout_path.write_text(
            "step time temperature potential meta rbias rct QC QC_S_sq REAXFF_BOND custom-force\n"
            "0 0.0 300.0 -12.5 1.0 2.0 3.0 -0.4 0.75 4.5 1.2\n"
            "10 0.02 301.0 -12.0 1.5 2.5 3.5 -0.5 0.76 4.6 1.4\n",
            encoding="utf-8",
        )

        datasets = {dataset.path: dataset.data for dataset in parse_legacy_output_file("mdout", mdout_path)}

        assert np.array_equal(datasets["/observables/all/step"], np.asarray([0, 10], dtype=np.int64))
        assert np.allclose(datasets["/observables/all/time"], np.asarray([0.0, 0.02], dtype=np.float64))
        assert np.allclose(datasets["/observables/all/temperature/value"], np.asarray([300.0, 301.0]))
        assert np.allclose(datasets["/observables/all/metadynamics/meta/value"], np.asarray([1.0, 1.5]))
        assert np.allclose(datasets["/observables/all/metadynamics/rbias/value"], np.asarray([2.0, 2.5]))
        assert np.allclose(datasets["/observables/all/metadynamics/rct/value"], np.asarray([3.0, 3.5]))
        assert np.allclose(datasets["/observables/all/qc/energy/value"], np.asarray([-0.4, -0.5]))
        assert np.allclose(datasets["/observables/all/qc/spin_square/value"], np.asarray([0.75, 0.76]))
        assert np.allclose(datasets["/observables/all/reaxff/REAXFF_BOND/value"], np.asarray([4.5, 4.6]))
        assert np.allclose(datasets["/observables/all/custom_force/value"], np.asarray([1.2, 1.4]))
        assert list(datasets["/parameters/sponge/mdout/columns/hdf5_name"]) == [
            "temperature",
            "potential",
            "metadynamics/meta",
            "metadynamics/rbias",
            "metadynamics/rct",
            "qc/energy",
            "qc/spin_square",
            "reaxff/REAXFF_BOND",
            "custom_force",
        ]
        assert list(datasets["/parameters/sponge/mdout/columns/category"]) == [
            "core",
            "core",
            "metadynamics",
            "metadynamics",
            "metadynamics",
            "qc",
            "qc",
            "reaxff",
            "custom",
        ]


def test_legacy_output_mdout_name_alignment():
    with tempfile.TemporaryDirectory() as tmpdir:
        mdout_path = Path(tmpdir) / "prod.out"
        mdout_path.write_text(
            "step time A-B A/B A_B A_B 1 TEMP(K) mdout_step mdout_time step time reaxff REAXFF REAXFF_BOND\n"
            "0 0.0 1 2 3 4 5 6 7 8 9 10 11 12 13\n"
            "1 0.1 14 15 16 17 18 19 20 21 22 23 24 25 26\n",
            encoding="utf-8",
        )

        datasets = {dataset.path: dataset.data for dataset in parse_legacy_output_file("mdout", mdout_path)}

        assert list(datasets["/parameters/sponge/mdout/columns/hdf5_name"]) == [
            "A_B",
            "A_B_1",
            "A_B_2",
            "A_B_3",
            "_1",
            "TEMP_K_",
            "mdout_step",
            "mdout_time",
            "mdout_step_1",
            "mdout_time_1",
            "reaxff",
            "reaxff/REAXFF",
            "reaxff/REAXFF_BOND",
        ]
        assert list(datasets["/parameters/sponge/mdout/columns/category"])[10:] == [
            "custom",
            "reaxff",
            "reaxff",
        ]
        assert np.allclose(datasets["/observables/all/A_B_1/value"], np.asarray([2.0, 15.0]))
        assert np.allclose(datasets["/observables/all/_1/value"], np.asarray([5.0, 18.0]))
        assert np.allclose(datasets["/observables/all/mdout_step_1/value"], np.asarray([9.0, 22.0]))
        assert np.allclose(datasets["/observables/all/mdout_time_1/value"], np.asarray([10.0, 23.0]))
        assert np.allclose(datasets["/observables/all/reaxff/value"], np.asarray([11.0, 24.0]))
        assert np.allclose(datasets["/observables/all/reaxff/REAXFF/value"], np.asarray([12.0, 25.0]))
        assert np.allclose(datasets["/observables/all/reaxff/REAXFF_BOND/value"], np.asarray([13.0, 26.0]))


def test_legacy_output_sits_and_meta_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        sits_traj_path = tmp_path / "sits_nk.dat"
        sits_rest_path = tmp_path / "sits_nk.rest"
        meta_path = tmp_path / "Meta_Potential.txt"
        direct_path = tmp_path / "Meta_directly.txt"
        hills_path = tmp_path / "myhill.log"
        history_path = tmp_path / "history.log"
        edge_path = tmp_path / "sumhill.log"
        np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32).tofile(sits_traj_path)
        sits_rest_path.write_text("3.0 4.0\n", encoding="utf-8")
        meta_path.write_text(
            "Meta potential\n"
            "0.0 1.0 0.5\n"
            "0.0 2.0 1.0\n"
            "2 2 4\n"
            "0.0 0.0 0.1 0.2 3.0\n"
            "0.5 0.0 0.3 0.4 4.0\n"
            "0.0 1.0 0.5 0.6 5.0\n"
            "0.5 1.0 0.7 0.8 6.0\n",
            encoding="utf-8",
        )
        direct_path.write_text("DIRECT\n", encoding="utf-8")
        hills_path.write_text("HILLS\n", encoding="utf-8")
        history_path.write_text("HISTORY\n", encoding="utf-8")
        edge_path.write_text("EDGE\n", encoding="utf-8")

        sits_traj = {
            dataset.path: dataset.data
            for dataset in parse_legacy_output_file("SITS_nk_traj_file", sits_traj_path, sits_nk_width=2)
        }
        sits_rest = {
            dataset.path: dataset.data
            for dataset in parse_legacy_output_file("SITS_nk_rest_file", sits_rest_path)
        }
        meta = {dataset.path: dataset.data for dataset in parse_legacy_output_file("meta_potential_out_file", meta_path)}
        direct = {dataset.path: dataset.data for dataset in parse_legacy_output_file("meta_direct_export", direct_path)}
        hills = {dataset.path: dataset.data for dataset in parse_legacy_output_file("meta_hills_log", hills_path)}
        history = {dataset.path: dataset.data for dataset in parse_legacy_output_file("meta_history_log", history_path)}
        edge = {dataset.path: dataset.data for dataset in parse_legacy_output_file("meta_edge_log", edge_path)}

        assert np.allclose(
            sits_traj["/observables/all/sits/SITS/nk/value"],
            np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        )
        assert np.array_equal(sits_traj["/observables/all/sits/SITS/nk/step"], np.asarray([0, 1], dtype=np.int64))
        assert np.allclose(
            sits_rest["/parameters/sponge/restart_exports/sits/SITS/nk/value"],
            np.asarray([3.0, 4.0], dtype=np.float32),
        )
        assert np.array_equal(
            meta["/parameters/sponge/metadynamics/default/potential_export/grid/count"],
            np.asarray([2, 2], dtype=np.int32),
        )
        assert np.allclose(
            meta["/parameters/sponge/metadynamics/default/potential_export/coordinate"],
            np.asarray([[0.0, 0.0], [0.5, 0.0], [0.0, 1.0], [0.5, 1.0]], dtype=np.float32),
        )
        assert np.allclose(
            meta["/parameters/sponge/metadynamics/default/potential_export/force"],
            np.asarray([[0.1, 0.2], [0.3, 0.4], [0.5, 0.6], [0.7, 0.8]], dtype=np.float32),
        )
        assert np.allclose(
            meta["/parameters/sponge/metadynamics/default/potential_export/value"],
            np.asarray([3.0, 4.0, 5.0, 6.0], dtype=np.float32),
        )
        assert direct["/parameters/sponge/metadynamics/default/direct_export"] == "DIRECT\n"
        assert hills["/parameters/sponge/metadynamics/default/hills"] == "HILLS\n"
        assert history["/parameters/sponge/metadynamics/default/history"] == "HISTORY\n"
        assert edge["/parameters/sponge/metadynamics/default/edge"] == "EDGE\n"


def test_legacy_output_qc_reaxff_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        qc_path = tmp_path / "qc_scf.out"
        eeq_path = tmp_path / "eeq_charges.txt"
        qc_path.write_text("SCF iteration 1\nEnergy -75.0\n", encoding="utf-8")
        eeq_path.write_text("1 -0.123456\n2 0.123456\n", encoding="utf-8")

        qc = {dataset.path: dataset.data for dataset in parse_legacy_output_file("qc_scf_output", qc_path)}
        eeq = {dataset.path: dataset.data for dataset in parse_legacy_output_file("reaxff_eeq_charges", eeq_path)}

        assert qc["/parameters/sponge/qc/scf_output"] == "SCF iteration 1\nEnergy -75.0\n"
        assert np.array_equal(
            eeq["/parameters/sponge/reaxff/eeq_charges/atom_index"],
            np.asarray([1, 2], dtype=np.int32),
        )
        assert np.allclose(
            eeq["/parameters/sponge/reaxff/eeq_charges/value"],
            np.asarray([-0.123456, 0.123456], dtype=np.float32),
        )
        assert eeq["/parameters/sponge/reaxff/eeq_charges/raw_text"] == "1 -0.123456\n2 0.123456\n"


def test_legacy_output_restart_parsers():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        rst_path = tmp_path / "prod.rst7"
        nhc_path = tmp_path / "nhc_restart.txt"
        rst_path.write_text(
            "restart state\n"
            "2 0.5\n"
            "1.0 2.0 3.0 4.0 5.0 6.0\n"
            "0.1 0.2 0.3 0.4 0.5 0.6\n"
            "10.0 20.0 30.0 90.0 90.0 90.0\n",
            encoding="utf-8",
        )
        nhc_path.write_text("0.1 0.2\n0.3 0.4\n", encoding="utf-8")

        rst = {dataset.path: dataset.data for dataset in parse_legacy_output_file("rst", rst_path)}
        nhc = {
            dataset.path: dataset.data
            for dataset in parse_legacy_output_file("nose_hoover_chain_restart_output", nhc_path)
        }

        assert np.array_equal(rst["/particles/all/step"], np.asarray([0], dtype=np.int64))
        assert np.allclose(rst["/particles/all/time"], np.asarray([0.5], dtype=np.float64))
        assert np.allclose(
            rst["/particles/all/position/value"],
            np.asarray([[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]], dtype=np.float32),
        )
        assert np.allclose(
            rst["/particles/all/velocity/value"],
            np.asarray([[[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]], dtype=np.float32),
        )
        assert np.allclose(rst["/particles/all/box/edges/value"][0], np.diag([10.0, 20.0, 30.0]).astype(np.float32))
        assert rst["/parameters/restart/source/title"] == "restart state"
        assert np.allclose(
            nhc["/parameters/restart/thermostat/nose_hoover_chain"],
            np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
        )


def test_legacy_outputs_to_h5md_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "legacy_outputs.spg.h5md"
        case_dir.mkdir()
        (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
        (case_dir / "mdin.spg.toml").write_text(
            'mass_in_file = "mass.txt"\n'
            'mdinfo = "prod.info"\n'
            'mdout = "prod.out"\n'
            'crd = "prod.dat"\n'
            'box = "prod.box"\n'
            'vel = "prod.vel"\n'
            'frc = "prod.frc"\n'
            'nose_hoover_chain_crd = "nhc.crd"\n',
            encoding="utf-8",
        )
        (case_dir / "prod.info").write_text("SPONGE mdinfo\n", encoding="utf-8")
        (case_dir / "prod.out").write_text(
            "step time temperature potential pressure\n"
            "0 0.0 300.0 -10.0 1.0\n"
            "5 0.01 301.0 -9.5 1.1\n",
            encoding="utf-8",
        )
        np.asarray(
            [
                [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
            ],
            dtype=np.float32,
        ).tofile(case_dir / "prod.dat")
        np.asarray(
            [
                [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
                [[0.2, 0.3, 0.4], [0.5, 0.6, 0.7]],
            ],
            dtype=np.float32,
        ).tofile(case_dir / "prod.vel")
        np.asarray(
            [
                [[-0.1, -0.2, -0.3], [-0.4, -0.5, -0.6]],
                [[-0.2, -0.3, -0.4], [-0.5, -0.6, -0.7]],
            ],
            dtype=np.float32,
        ).tofile(case_dir / "prod.frc")
        (case_dir / "prod.box").write_text(
            "10.0 20.0 30.0 90.0 90.0 90.0\n"
            "11.0 21.0 31.0 90.0 90.0 90.0\n",
            encoding="utf-8",
        )
        (case_dir / "nhc.crd").write_text("0.1 0.2\n0.3 0.4\n", encoding="utf-8")

        manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.legacy_import.mdout"]["status"] == "typed_converted"
        assert entries["output.legacy_import.crd"]["status"] == "typed_converted"
        assert entries["output.legacy_import.frc"]["status"] == "typed_converted"
        assert entries["output.legacy_import.nose_hoover_chain_crd"]["status"] == "typed_converted"
        assert entries["output.h5.trajectory"]["status"] == "generated"
        assert entries["output.h5.trajectory"]["bundle_file"] == "legacy_outputs.spg.h5md"

        with h5py.File(output_path, "r") as handle:
            assert np.array_equal(handle["/h5md"].attrs["version"], np.asarray([1, 1], dtype=np.int32))
            assert _h5_string(handle["/h5md/creator"].attrs["name"]) == "XPONGE"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/status"][()]) == "finalized"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([1.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert _h5_string_list(handle["/parameters/sponge/output/observable_streams"]) == ["all"]
            assert _h5_string(handle["/parameters/sponge/log/mdinfo_text"][()]) == "SPONGE mdinfo\n"
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars["mdinfo"] == "prod.info"
            assert sidecars["mdout"] == "prod.out"
            assert sidecars["crd"] == "prod.dat"
            assert sidecars["frc"] == "prod.frc"
            assert np.allclose(handle["/particles/all/position/value"][1, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))
            assert np.allclose(handle["/particles/all/force/value"][0, 0], np.asarray([-0.1, -0.2, -0.3], dtype=np.float32))
            assert np.allclose(handle["/particles/all/box/edges/value"][1], np.diag([11.0, 21.0, 31.0]).astype(np.float32))
            assert np.array_equal(handle["/particles/all/step"][...], np.asarray([0, 1], dtype=np.int64))
            assert _h5_string(handle["/particles/all/force/value"].attrs["unit"]) == "kcal mol-1 Angstrom-1"
            assert h5py.h5o.get_info(handle["/particles/all/force/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert np.array_equal(handle["/observables/all/step"][...], np.asarray([0, 5], dtype=np.int64))
            assert np.allclose(handle["/observables/all/temperature/value"][...], np.asarray([300.0, 301.0]))
            assert np.allclose(
                handle["/observables/all/thermostat/nose_hoover_chain/coordinate/value"][...],
                np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64),
            )
            assert h5py.h5o.get_info(handle["/observables/all/temperature/step"].id).addr == h5py.h5o.get_info(
                handle["/observables/all/step"].id
            ).addr
            assert _h5_string_list(handle["/parameters/sponge/mdout/columns/hdf5_name"]) == [
                "temperature",
                "potential",
                "pressure",
            ]


def test_legacy_outputs_restart_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "legacy_outputs.spg.h5md"
        restart_path = Path(tmpdir) / "custom_restart.spgr.h5"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'rst = "prod.rst7"\n'
            'nose_hoover_chain_restart_output = "nhc_restart.txt"\n'
            f'output_h5_restart_path = "{restart_path}"\n',
            encoding="utf-8",
        )
        (case_dir / "prod.rst7").write_text(
            "restart state\n"
            "2 0.5\n"
            "1.0 2.0 3.0 4.0 5.0 6.0\n"
            "0.1 0.2 0.3 0.4 0.5 0.6\n"
            "10.0 20.0 30.0 90.0 90.0 90.0\n",
            encoding="utf-8",
        )
        (case_dir / "nhc_restart.txt").write_text("0.1 0.2\n0.3 0.4\n", encoding="utf-8")

        manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.legacy_import.rst"]["status"] == "typed_converted"
        assert entries["output.legacy_import.rst"]["bundle_file"] == "custom_restart.spgr.h5"
        assert entries["output.legacy_import.nose_hoover_chain_restart_output"]["status"] == "typed_converted"
        assert entries["output.h5.restart"]["status"] == "generated"
        assert restart_path.exists()
        assert not output_path.exists()

        with h5py.File(restart_path, "r") as handle:
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.restart.h5"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/status"][()]) == "finalized"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([1], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([0], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([0.5]))
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars == {
                "rst": "prod.rst7",
                "nose_hoover_chain_restart_output": "nhc_restart.txt",
            }
            assert np.array_equal(handle["/particles/all/step"][...], np.asarray([0], dtype=np.int64))
            assert np.allclose(handle["/particles/all/time"][...], np.asarray([0.5], dtype=np.float64))
            assert np.allclose(
                handle["/particles/all/position/value"][0, 1],
                np.asarray([4.0, 5.0, 6.0], dtype=np.float32),
            )
            assert np.allclose(
                handle["/particles/all/velocity/value"][0, 0],
                np.asarray([0.1, 0.2, 0.3], dtype=np.float32),
            )
            assert np.allclose(
                handle["/parameters/restart/thermostat/nose_hoover_chain"][...],
                np.asarray([[0.1, 0.2], [0.3, 0.4]], dtype=np.float32),
            )
            assert h5py.h5o.get_info(handle["/particles/all/position/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert _h5_string(handle["/particles/all/position/value"].attrs["unit"]) == "Angstrom"


def test_legacy_outputs_diagnostic_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "diagnostic.spg.h5md"
        shard_dir = Path(tmpdir) / "diagnostic.spg.shards"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'qc_scf_output = "qc_scf.out"\n'
            "output_h5_trajectory_vds = true\n",
            encoding="utf-8",
        )
        (case_dir / "qc_scf.out").write_text("SCF iteration 1\nEnergy -75.0\n", encoding="utf-8")
        (case_dir / "eeq_charges.txt").write_text("1 -0.123456\n2 0.123456\n", encoding="utf-8")

        manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.h5.observable"]["status"] == "generated"
        assert entries["output.h5.observable"]["bundle_file"] == "diagnostic.spg.h5md"
        assert entries["output.legacy_import.qc_scf_output"]["status"] == "typed_converted"
        assert entries["output.legacy_import.reaxff_eeq_charges"]["status"] == "typed_converted"
        assert "output.h5.vds" not in entries
        assert not shard_dir.exists()

        with h5py.File(output_path, "r") as handle:
            assert "/particles" not in handle
            assert "/observables" not in handle
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "legacy_import"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([0], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([-1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([0.0]))
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars["qc_scf_output"] == "qc_scf.out"
            assert sidecars["reaxff_eeq_charges"] == "eeq_charges.txt"
            assert _h5_string(handle["/parameters/sponge/qc/scf_output"][()]) == "SCF iteration 1\nEnergy -75.0\n"
            assert np.array_equal(
                handle["/parameters/sponge/reaxff/eeq_charges/atom_index"][...],
                np.asarray([1, 2], dtype=np.int32),
            )
            assert np.allclose(
                handle["/parameters/sponge/reaxff/eeq_charges/value"][...],
                np.asarray([-0.123456, 0.123456], dtype=np.float32),
            )


def test_legacy_outputs_to_vds_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "prod.spg.h5md"
        shard_dir = Path(tmpdir) / "prod.spg.shards"
        case_dir.mkdir()
        (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
        (case_dir / "mdin.spg.toml").write_text(
            'mass_in_file = "mass.txt"\n'
            'mdout = "prod.out"\n'
            'crd = "prod.dat"\n'
            'box = "prod.box"\n'
            'nose_hoover_chain_crd = "nhc.crd"\n'
            "output_h5_trajectory_vds = true\n"
            "output_h5_trajectory_chunk_size = 1\n",
            encoding="utf-8",
        )
        (case_dir / "prod.out").write_text(
            "step time temperature meta rbias rct QC REAXFF_BOND\n"
            "0 0.0 300.0 1.0 2.0 3.0 -0.5 4.5\n"
            "5 0.01 301.0 1.5 2.5 3.5 -0.6 4.6\n",
            encoding="utf-8",
        )
        (case_dir / "nhc.crd").write_text("0.1 0.2\n0.3 0.4\n", encoding="utf-8")
        np.asarray(
            [
                [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
            ],
            dtype=np.float32,
        ).tofile(case_dir / "prod.dat")
        (case_dir / "prod.box").write_text(
            "10.0 20.0 30.0 90.0 90.0 90.0\n"
            "11.0 21.0 31.0 90.0 90.0 90.0\n",
            encoding="utf-8",
        )

        manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.h5.trajectory"]["status"] == "generated"
        assert entries["output.h5.vds"]["status"] == "generated"
        assert (shard_dir / "segment_000000.spg.h5md").exists()
        assert (shard_dir / "segment_000001.spg.h5md").exists()

        with h5py.File(output_path, "r") as handle:
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "vds"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([1.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert _h5_string_list(handle["/parameters/sponge/output/observable_streams"]) == ["all"]
            assert handle["/particles/all/position/value"].is_virtual
            assert handle["/observables/all/temperature/value"].is_virtual
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars["mdout"] == "prod.out"
            assert sidecars["crd"] == "prod.dat"
            assert sidecars["box"] == "prod.box"
            assert np.allclose(handle["/particles/all/position/value"][1, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))
            assert np.allclose(handle["/observables/all/temperature/value"][...], np.asarray([300.0, 301.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/shard_manifest/path"]) == [
                "prod.spg.shards/segment_000000.spg.h5md",
                "prod.spg.shards/segment_000001.spg.h5md",
            ]
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/frame_start"][...],
                np.asarray([0, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/observable_frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/nhc_frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/metadynamics_scalar_frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/qc_frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/reaxff_frame_count"][...],
                np.asarray([1, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/step_start"][...],
                np.asarray([0, 1], dtype=np.int64),
            )
            assert np.array_equal(
                handle["/parameters/sponge/output/shard_manifest/step_end"][...],
                np.asarray([0, 1], dtype=np.int64),
            )
            assert np.allclose(
                handle["/parameters/sponge/output/shard_manifest/time_start"][...],
                np.asarray([0.0, 1.0]),
            )
            assert np.allclose(
                handle["/parameters/sponge/output/shard_manifest/time_end"][...],
                np.asarray([0.0, 1.0]),
            )
            assert h5py.h5o.get_info(handle["/particles/all/position/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr

        with h5py.File(shard_dir / "segment_000001.spg.h5md", "r") as handle:
            assert np.array_equal(handle["/h5md"].attrs["version"], np.asarray([1, 1], dtype=np.int32))
            assert _h5_string(handle["/parameters/sponge/schema/name"][()]) == "sponge.output.h5md"
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/status"][()]) == "finalized"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([1], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([1.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert _h5_string_list(handle["/parameters/sponge/output/observable_streams"]) == ["all"]
            assert handle["/parameters/sponge/shard/index"][()] == 1
            assert handle["/parameters/sponge/shard/frame_start"][()] == 1
            assert handle["/parameters/sponge/shard/frame_count"][()] == 1
            assert handle["/parameters/sponge/shard/last_complete_step"][()] == 1
            assert np.isclose(handle["/parameters/sponge/shard/last_complete_time"][()], 1.0)
            assert np.allclose(handle["/particles/all/position/value"][0, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))
            assert handle["/particles/all/position/value"].attrs["unit"] == "Angstrom"
            assert handle["/particles/all/time"].attrs["unit"] == "ps"
            assert h5py.h5o.get_info(handle["/particles/all/position/step"].id).addr == h5py.h5o.get_info(
                handle["/particles/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/observables/all/temperature/step"].id).addr == h5py.h5o.get_info(
                handle["/observables/all/step"].id
            ).addr
            assert _h5_string_list(handle["/parameters/sponge/mdout/columns/hdf5_name"]) == [
                "temperature",
                "metadynamics/meta",
                "metadynamics/rbias",
                "metadynamics/rct",
                "qc/energy",
                "reaxff/REAXFF_BOND",
            ]


def test_legacy_outputs_observable_path_splits_from_trajectory_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        trajectory_path = Path(tmpdir) / "prod.spg.h5md"
        observable_path = Path(tmpdir) / "prod.obs.spg.h5md"
        shard_dir = Path(tmpdir) / "prod.spg.shards"
        case_dir.mkdir()
        (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
        (case_dir / "mdin.spg.toml").write_text(
            f'output_h5_observable_path = "{observable_path}"\n'
            "output_h5_trajectory_vds = true\n"
            "output_h5_trajectory_chunk_size = 1\n"
            'mass_in_file = "mass.txt"\n'
            'mdinfo = "prod.info"\n'
            'mdout = "prod.out"\n'
            'crd = "prod.dat"\n'
            'qc_scf_output = "qc_scf.out"\n',
            encoding="utf-8",
        )
        (case_dir / "prod.info").write_text("SPONGE mdinfo\n", encoding="utf-8")
        (case_dir / "prod.out").write_text(
            "step time temperature QC\n"
            "0 0.0 300.0 -0.5\n"
            "5 0.01 301.0 -0.6\n",
            encoding="utf-8",
        )
        (case_dir / "qc_scf.out").write_text("SCF iteration 1\nEnergy -75.0\n", encoding="utf-8")
        np.asarray(
            [
                [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
            ],
            dtype=np.float32,
        ).tofile(case_dir / "prod.dat")

        manifest = convert_legacy_outputs_to_bundle(case_dir, trajectory_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.legacy_import.crd"]["bundle_file"] == "prod.spg.h5md"
        assert entries["output.legacy_import.mdout"]["bundle_file"] == "prod.obs.spg.h5md"
        assert entries["output.legacy_import.mdinfo"]["bundle_file"] == "prod.obs.spg.h5md"
        assert entries["output.legacy_import.qc_scf_output"]["bundle_file"] == "prod.obs.spg.h5md"
        assert entries["output.h5.trajectory"]["bundle_file"] == "prod.spg.h5md"
        assert entries["output.h5.observable"]["bundle_file"] == "prod.obs.spg.h5md"
        assert entries["output.h5.vds"]["bundle_file"] == "prod.spg.h5md"
        assert (shard_dir / "segment_000000.spg.h5md").exists()
        assert observable_path.exists()

        with h5py.File(trajectory_path, "r") as handle:
            assert "/observables" not in handle
            assert "/parameters/sponge/log" not in handle
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "vds"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([1], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([1.0]))
            assert _h5_string_list(handle["/parameters/sponge/output/particle_streams"]) == ["all"]
            assert "/parameters/sponge/output/observable_streams" not in handle
            assert _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]) == ["crd"]
            assert _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]) == ["prod.dat"]
            assert handle["/particles/all/position/value"].is_virtual
            assert np.allclose(handle["/particles/all/position/value"][1, 1], np.asarray([4.5, 5.5, 6.5], dtype=np.float32))

        with h5py.File(observable_path, "r") as handle:
            assert "/particles" not in handle
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "legacy_import"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([5], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([0.01]))
            assert "/parameters/sponge/output/particle_streams" not in handle
            assert _h5_string_list(handle["/parameters/sponge/output/observable_streams"]) == ["all"]
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars == {
                "mdinfo": "prod.info",
                "mdout": "prod.out",
                "qc_scf_output": "qc_scf.out",
            }
            assert _h5_string(handle["/parameters/sponge/log/mdinfo_text"][()]) == "SPONGE mdinfo\n"
            assert _h5_string(handle["/parameters/sponge/qc/scf_output"][()]) == "SCF iteration 1\nEnergy -75.0\n"
            assert np.allclose(handle["/observables/all/temperature/value"][...], np.asarray([300.0, 301.0]))
            assert np.allclose(handle["/observables/all/qc/energy/value"][...], np.asarray([-0.5, -0.6]))
            assert h5py.h5o.get_info(handle["/observables/all/qc/energy/step"].id).addr == h5py.h5o.get_info(
                handle["/observables/all/step"].id
            ).addr


def test_legacy_outputs_observable_only_bundle():
    try:
        import h5py
    except ImportError:
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "observable_only.spg.h5md"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'default_out_file_prefix = "prod"\n'
            'mdout = "prod.out"\n'
            'SITS_nk_traj_file = "sits_nk.dat"\n'
            'SITS_nk_rest_file = "sits_nk.rest"\n',
            encoding="utf-8",
        )
        (case_dir / "prod.out").write_text(
            "step time temperature meta\n"
            "0 0.0 300.0 1.0\n"
            "5 0.01 301.0 1.5\n",
            encoding="utf-8",
        )
        np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32).tofile(case_dir / "sits_nk.dat")
        (case_dir / "sits_nk.rest").write_text("3.0 4.0\n", encoding="utf-8")
        (case_dir / "prod_Meta_Potential.txt").write_text(
            "Meta potential\n"
            "0.0 1.0 0.5\n"
            "0.0 2.0 1.0\n"
            "2 2 4\n"
            "0.0 0.0 0.1 0.2 3.0\n"
            "0.5 0.0 0.3 0.4 4.0\n"
            "0.0 1.0 0.5 0.6 5.0\n"
            "0.5 1.0 0.7 0.8 6.0\n",
            encoding="utf-8",
        )
        (case_dir / "Meta_directly.txt").write_text("DIRECT\n", encoding="utf-8")
        (case_dir / "myhill.log").write_text("HILLS\n", encoding="utf-8")
        (case_dir / "history.log").write_text("HISTORY\n", encoding="utf-8")
        (case_dir / "sumhill.log").write_text("EDGE\n", encoding="utf-8")

        manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
        entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
        assert entries["output.h5.observable"]["status"] == "generated"
        assert entries["output.h5.observable"]["bundle_file"] == "observable_only.spg.h5md"
        assert entries["output.legacy_import.SITS_nk_traj_file"]["status"] == "typed_converted"
        assert entries["output.legacy_import.SITS_nk_rest_file"]["status"] == "typed_converted"
        assert entries["output.legacy_import.meta_potential_out_file"]["status"] == "typed_converted"
        assert entries["output.legacy_import.meta_potential_out_file"]["source_path"].endswith(
            "prod_Meta_Potential.txt"
        )
        assert entries["output.legacy_import.meta_direct_export"]["status"] == "typed_converted"
        assert entries["output.legacy_import.meta_hills_log"]["status"] == "typed_converted"
        assert entries["output.legacy_import.meta_history_log"]["status"] == "typed_converted"
        assert entries["output.legacy_import.meta_edge_log"]["status"] == "typed_converted"

        with h5py.File(output_path, "r") as handle:
            assert "/particles" not in handle
            assert np.array_equal(handle["/h5md"].attrs["version"], np.asarray([1, 1], dtype=np.int32))
            assert _h5_string(handle["/parameters/sponge/schema/version"][()]) == "xponge.legacy_to_bundle.v1"
            assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == "legacy_import"
            assert np.array_equal(handle["/parameters/sponge/output/frame_count"][...], np.asarray([2], dtype=np.int64))
            assert np.array_equal(
                handle["/parameters/sponge/output/last_complete_step"][...],
                np.asarray([5], dtype=np.int64),
            )
            assert np.allclose(handle["/parameters/sponge/output/last_complete_time"][...], np.asarray([0.01]))
            assert _h5_string_list(handle["/parameters/sponge/output/observable_streams"]) == ["all"]
            sidecars = dict(
                zip(
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/key"]),
                    _h5_string_list(handle["/parameters/sponge/files/legacy_sidecars/path"]),
                )
            )
            assert sidecars["mdout"] == "prod.out"
            assert sidecars["SITS_nk_traj_file"] == "sits_nk.dat"
            assert sidecars["meta_potential_out_file"] == "prod_Meta_Potential.txt"
            assert np.array_equal(handle["/observables/all/step"][...], np.asarray([0, 5], dtype=np.int64))
            assert np.allclose(handle["/observables/all/temperature/value"][...], np.asarray([300.0, 301.0]))
            assert np.allclose(handle["/observables/all/metadynamics/meta/value"][...], np.asarray([1.0, 1.5]))
            assert np.allclose(
                handle["/observables/all/sits/SITS/nk/value"][...],
                np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
            )
            assert np.allclose(
                handle["/parameters/sponge/restart_exports/sits/SITS/nk/value"][...],
                np.asarray([3.0, 4.0], dtype=np.float32),
            )
            assert np.allclose(
                handle["/parameters/sponge/metadynamics/default/potential_export/value"][...],
                np.asarray([3.0, 4.0, 5.0, 6.0], dtype=np.float32),
            )
            assert _h5_string(handle["/parameters/sponge/metadynamics/default/direct_export"][()]) == "DIRECT\n"
            assert _h5_string(handle["/parameters/sponge/metadynamics/default/hills"][()]) == "HILLS\n"
            assert _h5_string(handle["/parameters/sponge/metadynamics/default/history"][()]) == "HISTORY\n"
            assert _h5_string(handle["/parameters/sponge/metadynamics/default/edge"][()]) == "EDGE\n"
            assert h5py.h5o.get_info(handle["/observables/all/temperature/step"].id).addr == h5py.h5o.get_info(
                handle["/observables/all/step"].id
            ).addr
            assert h5py.h5o.get_info(handle["/observables/all/metadynamics/meta/step"].id).addr == h5py.h5o.get_info(
                handle["/observables/all/step"].id
            ).addr


def test_legacy_output_vds_mode_matrix_writes():
    try:
        import h5py
    except ImportError:
        return

    for mode in ("normal", "rerun"):
        for vds in (False, True):
            with tempfile.TemporaryDirectory() as tmpdir:
                case_dir = Path(tmpdir) / "case"
                output_path = Path(tmpdir) / f"{mode}_{'vds' if vds else 'single'}.spg.h5md"
                case_dir.mkdir()
                (case_dir / "mass.txt").write_text("2\n1.0\n16.0\n", encoding="utf-8")
                (case_dir / "mdin.spg.toml").write_text(
                    f'mode = "{mode}"\n'
                    'mass_in_file = "mass.txt"\n'
                    'mdout = "prod.out"\n'
                    'crd = "prod.dat"\n'
                    'box = "prod.box"\n'
                    f"output_h5_trajectory_vds = {'true' if vds else 'false'}\n"
                    "output_h5_trajectory_chunk_size = 1\n",
                    encoding="utf-8",
                )
                (case_dir / "prod.out").write_text(
                    "step time temperature\n"
                    "0 0.0 300.0\n"
                    "5 0.01 301.0\n",
                    encoding="utf-8",
                )
                np.asarray(
                    [
                        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                        [[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]],
                    ],
                    dtype=np.float32,
                ).tofile(case_dir / "prod.dat")
                (case_dir / "prod.box").write_text(
                    "10.0 20.0 30.0 90.0 90.0 90.0\n"
                    "11.0 21.0 31.0 90.0 90.0 90.0\n",
                    encoding="utf-8",
                )

                manifest = convert_legacy_outputs_to_bundle(case_dir, output_path)
                assert manifest.mode == mode
                entries = {entry["contract_id"]: entry for entry in manifest.to_dict()["entries"]}
                assert entries["output.legacy_import.mdout"]["status"] == "typed_converted"
                assert entries["output.legacy_import.crd"]["status"] == "typed_converted"

                with h5py.File(output_path, "r") as handle:
                    assert np.allclose(
                        handle["/particles/all/position/value"][1, 1],
                        np.asarray([4.5, 5.5, 6.5], dtype=np.float32),
                    )
                    assert np.allclose(handle["/observables/all/temperature/value"][...], np.asarray([300.0, 301.0]))
                    assert _h5_string(handle["/parameters/sponge/output/mode"][()]) == ("vds" if vds else "legacy_import")
                    assert handle["/particles/all/position/value"].is_virtual == vds
                    assert handle["/observables/all/temperature/value"].is_virtual == vds
                    if vds:
                        assert entries["output.h5.vds"]["status"] == "generated"
                        assert _h5_string_list(handle["/parameters/sponge/output/shard_manifest/status"]) == [
                            "complete",
                            "complete",
                        ]
                    else:
                        assert "output.h5.vds" not in entries
                        assert "/parameters/sponge/output/shard_manifest" not in handle


def test_legacy_to_bundle_cli_dry_run():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_dir = Path(tmpdir) / "out"
        case_dir.mkdir()
        write_basic_case(case_dir)

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "Xponge",
                "legacy-to-bundle",
                str(case_dir),
                "-o",
                str(output_dir),
                "--dry-run",
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, result.stderr
        assert "converted entries:" in result.stdout
        assert "typed converted entries:" in result.stdout
        typed_line = next(line for line in result.stdout.splitlines() if line.startswith("typed converted entries:"))
        assert int(typed_line.split(":", 1)[1].strip()) > 0
        assert "compatibility imported entries:" not in result.stdout
        assert "unsupported entries:" not in result.stdout


def test_legacy_outputs_to_bundle_cli_dry_run():
    with tempfile.TemporaryDirectory() as tmpdir:
        case_dir = Path(tmpdir) / "case"
        output_path = Path(tmpdir) / "legacy_outputs.spg.h5md"
        case_dir.mkdir()
        (case_dir / "mdin.spg.toml").write_text(
            'mdout = "prod.out"\n'
            'mdinfo = "prod.info"\n',
            encoding="utf-8",
        )
        (case_dir / "prod.out").write_text(
            "step time temperature\n"
            "0 0.0 300.0\n"
            "5 0.01 301.0\n",
            encoding="utf-8",
        )
        (case_dir / "prod.info").write_text("SPONGE mdinfo\n", encoding="utf-8")

        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "Xponge",
                "legacy-outputs-to-bundle",
                str(case_dir),
                "-o",
                str(output_path),
                "--dry-run",
            ],
            check=False,
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0, result.stderr
        assert "typed converted output entries: 2" in result.stdout
        assert "missing legacy output entries:" not in result.stdout
        assert not output_path.exists()
