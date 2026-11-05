"""Same public-contract and energy tests run in isolated provider processes."""
import copy
import math
import os
from io import StringIO
from uuid import uuid4

import pytest

if os.environ.get('EXPLICIT_PROVIDER') == 'cpp':
    import XpongeCPP as x
    from XpongeCPP.forcefield.special import fep
else:
    import Xponge as x
    import Xponge.forcefield.amber.gaff
    from Xponge.forcefield.special import fep


def molecule():
    if x.__name__ == 'XpongeCPP':
        label = 'EXPLICITTEST'
        x.register_amber_lj_parameter(label, label, .1, 1.8)
        x.register_amber_bond_parameter(label, label, 100., 1.5)
        x.register_amber_angle_parameter([label]*3, 50., 1.9)
        x.register_amber_proper_dihedral_parameter([label]*4, 1, 2., .3, reset=True)
        return x.load_mol2(StringIO('''@<TRIPOS>MOLECULE
LIG
4 3 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C0 0.0 0.0 0.0 EXPLICITTEST 1 LIG 0.2
2 C1 1.0 1.0 0.2 EXPLICITTEST 1 LIG -0.2
3 C2 2.0 0.0 0.4 EXPLICITTEST 1 LIG 0.3
4 C3 3.0 1.0 0.6 EXPLICITTEST 1 LIG -0.3
@<TRIPOS>BOND
1 1 2 1
2 2 3 1
3 3 4 1
'''))
    residue = x.ResidueType(name='EXPLICIT_' + uuid4().hex)
    for i in range(4):
        residue.Add_Atom('C'+str(i), 'c3', float(i), float(i % 2), .2*i)
        residue.atoms[-1].charge = (.2, -.2, .3, -.3)[i]
    for i in range(3):
        residue.Add_Connectivity(residue.atoms[i], residue.atoms[i+1])
    mol = x.Molecule('explicit')
    mol.Add_Residue(residue)
    x.build_bonded_force(mol)
    return mol


def parameterized(scale=1.):
    mol = molecule()
    top = x.get_explicit_topology(mol)
    top['bonds'][0][2:] = [123. * scale, 1.4]
    top['dihedrals'] = [[0, 1, 2, 3, 1, 2.*scale, .3],
                        [0, 1, 2, 3, 1, 3.*scale, .3],
                        [2, 0, 1, 3, 2, 7.*scale, .7]]
    top['nb14s'] = [[0, 3, .5, .5]]
    x.set_explicit_topology(mol, top)
    for i in range(4):
        x.set_atom_lj(mol, i, .13 + .01*i, 1.6 + .1*i)
    return mol, top


def test_roundtrip_copy_export_and_environment(tmp_path):
    mol, top = parameterized()
    assert x.get_explicit_topology(mol) == top
    assert x.get_explicit_topology(mol.deepcopy()) == top
    other, _ = parameterized(2.)
    # Disjoint assembly must offset all atom references, without retyping.
    merged = mol | other
    mt = x.get_explicit_topology(merged)
    assert mt['dihedrals'] == top['dihedrals'] + [
        [*[i+4 for i in row[:4]], *row[4:]]
        for row in x.get_explicit_topology(other)['dihedrals']]
    x.save_sponge_input(mol, prefix='explicit', dirname=str(tmp_path))
    rows = (tmp_path/'explicit_dihedral.txt').read_text().splitlines()
    assert int(rows[0]) == 3
    assert sorted(float(r.split()[5]) for r in rows[1:] if r.strip()) == [2., 3., 7.]
    assert '0 3' in (tmp_path/'explicit_nb14.txt').read_text()
    assert x.get_explicit_topology(mol) == top


@pytest.mark.parametrize('bad', ['missing', 'nan', 'index', 'fractional', 'negative', 'self', 'asymmetric', 'periodicity'])
def test_invalid_set_is_transactional(bad):
    mol, original = parameterized()
    top = copy.deepcopy(original)
    if bad == 'missing': del top['bonds']
    elif bad == 'nan': top['bonds'][0][2] = float('nan')
    elif bad == 'index': top['angles'][0][0] = 40
    elif bad == 'fractional': top['angles'][0][0] = .5
    elif bad == 'negative': top['bonds'][0][2] = -1
    elif bad == 'self': top['bonds'][0][1] = top['bonds'][0][0]
    elif bad == 'asymmetric': top['exclusions'][0].remove(1)
    elif bad == 'periodicity': top['dihedrals'][0][4] = 1.5
    with pytest.raises(ValueError): x.set_explicit_topology(mol, top)
    assert x.get_explicit_topology(mol) == original


def test_lj_does_not_mutate_other_molecule(tmp_path):
    first, _ = parameterized()
    second = first.deepcopy()
    x.save_sponge_input(first, prefix='before', dirname=str(tmp_path))
    x.set_atom_lj(second, 0, .9, 2.8)
    x.save_sponge_input(first, prefix='after', dirname=str(tmp_path))
    assert (tmp_path/'before_LJ.txt').read_bytes() == (tmp_path/'after_LJ.txt').read_bytes()
    x.save_sponge_input(second, prefix='changed', dirname=str(tmp_path))
    assert (tmp_path/'before_LJ.txt').read_bytes() != (tmp_path/'changed_LJ.txt').read_bytes()


def test_fep_preserves_periodic_energy_and_endpoint_ownership(tmp_path):
    a, ta = parameterized()
    b, tb = parameterized(3.)
    for lam in (0., .5, 1.):
        merged = fep.Merge_Force_Field(a, b, lam)
        x.save_sponge_input(merged, prefix='merged', dirname=str(tmp_path))
        rows = [list(map(float, line.split())) for line in
                (tmp_path/'merged_dihedral.txt').read_text().splitlines()[1:] if line.strip()]
        # Evaluate each distinct ordered torsion separately, including a
        # non-proper atom ordering. Reversal has the same periodic angle.
        for ids, phi in (((0,1,2,3), .67), ((2,0,1,3), -1.2)):
            selected = [r for r in rows if tuple(r[:4]) in (ids, ids[::-1])]
            energy = sum(r[5]*(1+math.cos(r[4]*phi-r[6])) for r in selected)
            expected = (1+2*lam)*sum(r[5]*(1+math.cos(r[4]*phi-r[6]))
                                     for r in ta['dihedrals'] if tuple(r[:4]) == ids)
            assert energy == pytest.approx(expected, abs=1e-5)
        assert x.get_explicit_topology(a) == ta
        assert x.get_explicit_topology(b) == tb


def test_origin_rebuild_of_invalidated_materialization_is_rejected():
    if x.__name__ == 'XpongeCPP':
        pytest.skip('Python object model invalidation guard')
    mol, _ = parameterized()
    mol.built = False
    with pytest.raises(ValueError, match='invalidated'):
        x.build_bonded_force(mol)


@pytest.mark.parametrize('shared', [2, 4])
def test_dual_topology_keeps_explicit_replacement_parameters(tmp_path, shared):
    a, ta = parameterized()
    b, tb = parameterized(3.)
    if x.__name__ == 'XpongeCPP':
        left, right, _ = fep.merge_dual_topology(a, 0, b, {i:i for i in range(shared)})
    else:
        from Xponge.assign import get_assignment_from_residuetype
        mapping = tuple(range(shared))
        path = tmp_path/'mcs.txt'
        path.write_text(f'mcs atoms in r1: {mapping}\nmcs atoms in r2: {mapping}\n')
        left, right, _ = fep.merge_dual_topology(
            a, a.residues[0], b,
            get_assignment_from_residuetype(a.residues[0]),
            get_assignment_from_residuetype(b.residues[0]), imcs=str(path))
    lt, rt = x.get_explicit_topology(left), x.get_explicit_topology(right)
    # With exclusive branches the opposite endpoint's bonded scaffold remains.
    assert sorted(r[5] for r in lt['dihedrals']) == sorted([2.,3.,7.] + ([6.,9.,21.] if shared == 2 else []))
    assert sorted(r[5] for r in rt['dihedrals']) == sorted([6.,9.,21.] + ([2.,3.,7.] if shared == 2 else []))
    assert x.get_explicit_topology(a) == ta
    assert x.get_explicit_topology(b) == tb
    for lam in (0., .5, 1.):
        merged = fep.Merge_Force_Field(left, right, lam)
        x.save_sponge_input(merged, prefix=f'dual{lam}', dirname=str(tmp_path))
        assert (tmp_path/f'dual{lam}_dihedral.txt').is_file()
