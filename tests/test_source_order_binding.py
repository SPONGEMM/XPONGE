import importlib
import json
import warnings
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import Xponge
from Xponge.analysis import md_analysis as xmda
from Xponge.io_bundle.errors import BundleTopologyError, BundleTrajectoryError
from Xponge.io_bundle.source_order import (
    bind_bundle_source_order, make_source_order_binding, read_source_order_binding,
    validate_source_order_binding, validate_trajectory_source_order,
)


def _text(handle, path, value):
    if path in handle:
        del handle[path]
    handle.create_dataset(path, data=value, dtype=h5py.string_dtype('utf-8'))


def _bundle(tmp_path):
    paths = SimpleNamespace(topology=tmp_path/'top.h5', restart=tmp_path/'rst.h5')
    with h5py.File(paths.topology, 'w') as f:
        f['/topology/atom_count'] = 2
        _text(f, '/topology/atom_order_hash', 'native-order')
        _text(f, '/topology/topology_hash', 'native-topology')
    with h5py.File(paths.restart, 'w') as f:
        _text(f, '/run/atom_order_hash', 'native-order')
    return paths


def _trajectory(path, binding, stream='all'):
    with h5py.File(path, 'w') as f:
        f.require_group('/h5md').attrs['version'] = [1, 1]
        f.require_group('/h5md/creator').attrs.update(name='SPONGE', version='test')
        _text(f, '/parameters/sponge/schema/name', 'sponge.output.h5md')
        _text(f, '/parameters/sponge/schema/version', 'sponge.output.v2')
        _text(f, '/parameters/sponge/output/status', 'finalized')
        for key in ('atom_order_hash', 'topology_hash'):
            _text(f, '/parameters/sponge/topology_compatibility/'+key, binding[key])
        p = f.require_group('/particles/'+stream)
        p['step'] = [0, 1]
        p['time'] = [0., 1.]
        p['time'].attrs['unit'] = 'ps'
        position = p.create_group('position')
        position['value'] = np.zeros((2, 2, 3), dtype=np.float32)
        position['value'].attrs['unit'] = 'Angstrom'
        position['step'], position['time'] = p['step'], p['time']
        box = p.create_group('box')
        box.attrs.update(dimension=3, boundary=np.asarray(['periodic']*3, dtype='S8'))
        edges = box.create_group('edges')
        edges['value'] = np.asarray([np.eye(3)*10]*2)
        edges['value'].attrs['unit'] = 'Angstrom'
        edges['step'], edges['time'] = p['step'], p['time']


def _cif_pair(tmp_path):
    import hashlib
    atoms = [dict(simulation_index=i, external_id=f'atom:{i}', canonical_atom_id=i+1,
                  simulation_residue_index=0, simulation_residue_id='res:0') for i in range(2)]
    digest = hashlib.sha256(json.dumps({'atom_mapping': atoms}, sort_keys=True,
                                      separators=(',', ':')).encode()).hexdigest()
    cif = tmp_path/'system_trajectory_topology.cif'
    tags = ['id', 'type_symbol', 'label_atom_id', 'label_comp_id', 'label_asym_id',
            'label_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    mapping_tags = ['simulation_index', 'atom_site_id', 'canonical_atom_id', 'external_id',
                    'simulation_residue_index', 'simulation_residue_id', 'mapping_hash']
    cif.write_text('data_test\nloop_\n' + '\n'.join('_atom_site.'+tag for tag in tags)
                   + '\n1 C C1 MOL A 1 0 0 0\n2 C C2 MOL A 1 1 0 0\n#\nloop_\n'
                   + '\n'.join('_mokda_trajectory_mapping.'+tag for tag in mapping_tags)
                   + '\n' + '\n'.join(f'{i} {i+1} {i+1} atom:{i} 0 res:0 {digest}' for i in range(2))+'\n')
    binding = make_source_order_binding(['atom:0', 'atom:1'], 'native-order', 'native-topology')
    document = dict(schema='sponge-atom-order-mapping', schema_version=1,
                    mapping_hash=digest, atoms=atoms, topology_binding=binding)
    sidecar = tmp_path/'system_atom_order_mapping.json'
    sidecar.write_text(json.dumps(document))
    trajectory = tmp_path/'trajectory.h5md'
    _trajectory(trajectory, binding)
    return cif, sidecar, trajectory, document


def test_binding_changes_for_identity_permutation():
    a = make_source_order_binding(['C1', 'C2'], 'same-physical-properties', 'same-topology')
    b = make_source_order_binding(['C2', 'C1'], 'same-physical-properties', 'same-topology')
    assert a['atom_order_hash'] != b['atom_order_hash']
    assert validate_source_order_binding(a, ['C1', 'C2']) == a
    with pytest.raises(ValueError, match='identity sequence'):
        validate_source_order_binding(a, ['C2', 'C1'])


@pytest.mark.parametrize('ids', [[], [''], ['a', 'a'], ['a', None]])
def test_binding_rejects_invalid_identities(ids):
    with pytest.raises(ValueError):
        make_source_order_binding(ids, 'order', 'topology')


def test_staged_binding_updates_restart_and_checks_metadata(tmp_path):
    paths = _bundle(tmp_path)
    binding = bind_bundle_source_order(paths, ['a', 'b'])
    assert read_source_order_binding(paths.topology, ['a', 'b']) == binding
    with h5py.File(paths.restart) as f:
        assert f['/run/atom_order_hash'].asstr()[()] == binding['atom_order_hash']
    with pytest.raises(ValueError):
        read_source_order_binding(paths.topology, ['b', 'a'])
    with h5py.File(paths.topology, 'r+') as f:
        _text(f, '/topology/atom_order_hash', 'wrong')
    with pytest.raises(ValueError, match='atom_order_hash'):
        read_source_order_binding(paths.topology, ['a', 'b'])


def test_cif_h5md_bound_without_native_topology(tmp_path):
    cif, _, trajectory, _ = _cif_pair(tmp_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        universe = xmda.load_cif_h5md_universe(cif, trajectory)
    try:
        assert universe.cif_metadata['trajectory_order_validation']['verified']
        assert not any('CIF/H5MD compatibility' in str(w.message) for w in caught)
    finally:
        universe.trajectory.close()


@pytest.mark.parametrize('field', ['atom_order_hash', 'topology_hash'])
def test_cif_rejects_same_count_wrong_trajectory(tmp_path, field):
    cif, _, trajectory, document = _cif_pair(tmp_path)
    with h5py.File(trajectory, 'r+') as f:
        _text(f, '/parameters/sponge/topology_compatibility/'+field, 'other-system')
    with pytest.raises(BundleTrajectoryError, match=field):
        xmda.load_cif_h5md_universe(cif, trajectory, strict=False)


@pytest.mark.parametrize('change', ['permuted', 'version', 'algorithm', 'null', 'count'])
def test_cif_rejects_invalid_sidecar_binding(tmp_path, change):
    cif, sidecar, trajectory, document = _cif_pair(tmp_path)
    if change == 'permuted':
        document['topology_binding'] = make_source_order_binding(['atom:1','atom:0'], 'native-order', 'native-topology')
    elif change == 'null':
        document['topology_binding'] = None
    else:
        key = {'version': 'schema_version', 'algorithm': 'hash_algorithm', 'count': 'atom_count'}[change]
        document['topology_binding'][key] = 'invalid'
    sidecar.write_text(json.dumps(document))
    with pytest.raises(BundleTopologyError, match='topology_binding'):
        xmda.load_cif_h5md_universe(cif, trajectory)


@pytest.mark.parametrize('legacy', ['no_sidecar', 'old_sidecar', 'missing_hash'])
def test_legacy_inputs_still_load_unverified(tmp_path, legacy):
    cif, sidecar, trajectory, document = _cif_pair(tmp_path)
    if legacy == 'no_sidecar':
        sidecar.unlink()
    elif legacy == 'old_sidecar':
        del document['topology_binding']
        sidecar.write_text(json.dumps(document))
    else:
        with h5py.File(trajectory, 'r+') as f:
            del f['/parameters/sponge/topology_compatibility/atom_order_hash']
    with pytest.warns(RuntimeWarning, match='CIF/H5MD compatibility'):
        universe = xmda.load_cif_h5md_universe(cif, trajectory)
    assert not universe.cif_metadata['trajectory_order_validation']['verified']
    universe.trajectory.close()


def test_global_hash_does_not_verify_custom_particle_stream(tmp_path):
    binding = make_source_order_binding(['a', 'b'], 'order', 'topology')
    trajectory = tmp_path/'custom.h5md'
    _trajectory(trajectory, binding, stream='subset')
    assert not validate_trajectory_source_order(trajectory, binding)['verified']


def test_real_saver_binds_actual_serialization_order(tmp_path):
    importlib.import_module('Xponge.forcefield.amber.ff14sb')
    molecule = Xponge.get_peptide_from_sequence('AA')
    if hasattr(molecule, 'get_atoms'):
        molecule.get_atoms()
    ids = [f'atom:{i}' for i, _ in enumerate(molecule.atoms)]
    bindings = []
    for prefix, identities in [('first', ids), ('second', list(reversed(ids)))]:
        _, mapping = Xponge.save_sponge_input(molecule, prefix, dirname=str(tmp_path),
                                             format='bundle', source_atom_ids=identities,
                                             return_mapping=True)
        serialized = [row['source_atom_id'] for row in mapping]
        topology = tmp_path/f'{prefix}_topology.spgt.h5'
        binding = read_source_order_binding(topology, serialized)
        assert binding is not None
        with h5py.File(tmp_path/f'{prefix}_restart.spgr.h5') as f:
            assert f['/run/atom_order_hash'].asstr()[()] == binding['atom_order_hash']
        bindings.append(binding)
    assert bindings[0]['base_atom_order_hash'] == bindings[1]['base_atom_order_hash']
    assert bindings[0]['atom_order_hash'] != bindings[1]['atom_order_hash']
