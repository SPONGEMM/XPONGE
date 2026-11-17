import numpy as np
import pytest

from Xponge.io_bundle.topology_parsers import parse_topology_file


@pytest.mark.parametrize("key,group,width,param,shape", [
    ("bond_in_file", "bond", 2, "k", (0,)),
    ("angle_in_file", "angle", 3, "k", (0,)),
    ("dihedral_in_file", "dihedral", 4, "k", (0,)),
    ("improper_dihedral_in_file", "improper", 4, "pk", (0,)),
    ("nb14_in_file", "nb14", 2, "params", (0, 2)),
    ("nb14_extra_in_file", "nb14_extra", 2, "params", (0, 3)),
])
def test_zero_count_retains_typed_matrix_axes(tmp_path, key, group, width, param, shape):
    path = tmp_path / "empty.txt"
    path.write_text("0\n")
    data = {item.path: item.data for item in parse_topology_file(key, path)}
    root = "/forcefield/" + group
    assert data[root + "/atoms"].shape == (0, width)
    assert data[root + "/atoms"].dtype == np.int32
    assert data[root + "/" + param].shape == shape
    assert int(data[root + "/count"]) == 0
