"""FEP parameter ownership, endpoint identity, and lambda-order regressions."""
import copy
import math
from uuid import uuid4

import pytest
import Xponge as x
import Xponge.forcefield.amber.gaff
from Xponge.forcefield.base import dihedral_base
from Xponge.forcefield.special import fep


def _molecule():
    residue = x.ResidueType(name='ISOLATION_' + uuid4().hex)
    for i in range(4):
        residue.Add_Atom('C' + str(i), 'c3', float(i), float(i % 2), .2 * i)
        residue.atoms[-1].charge = 0.
    for i in range(3):
        residue.Add_Connectivity(residue.atoms[i], residue.atoms[i + 1])
    mol = x.Molecule('isolation')
    mol.Add_Residue(residue)
    x.build_bonded_force(mol)
    return mol


def _torsion(mol, values):
    term = dihedral_base.ProperType.entity(mol.atoms[:4], dihedral_base.ProperType.get_type('c3-c3-c3-c3'))
    term.multiple_numbers = len(values)
    term.ks = [v[0] for v in values]
    term.phi0s = [v[1] for v in values]
    term.periodicitys = [v[2] for v in values]
    mol.bonded_forces['dihedral'] = [term]
    return term


def _fingerprint(mol):
    return copy.deepcopy({kind: [term.contents for term in terms] for kind, terms in mol.bonded_forces.items()})


def _energy(mol, phi):
    return sum(term.ks[i] * (1 + math.cos(term.periodicitys[i] * phi - term.phi0s[i]))
               for term in mol.bonded_forces.get('dihedral', []) for i in range(term.multiple_numbers))


def test_force_entities_own_template_and_copy_parameter_lists():
    mol = _molecule()
    template = dihedral_base.ProperType.get_type('c3-c3-c3-c3')
    original = copy.deepcopy(template.contents)
    first = dihedral_base.ProperType.entity(mol.atoms[:4], template)
    second = dihedral_base.ProperType.entity(mol.atoms[:4], template)
    assert first.ks is not template.ks and first.ks is not second.ks
    first.ks[0] += 123.
    assert template.contents == original
    _torsion(mol, [(2., 0., 1), (3., math.pi, 2)])
    duplicate = mol.deepcopy()
    a, b = mol.bonded_forces['dihedral'][0], duplicate.bonded_forces['dihedral'][0]
    assert b.type is a.type  # Registry/type identity is intentionally retained.
    assert all(atom in duplicate.atoms for atom in b.atoms)
    assert all(atom not in mol.atoms for atom in b.atoms)
    assert b.ks is not a.ks and b.phi0s is not a.phi0s and b.periodicitys is not a.periodicitys
    b.ks[0] += 12.
    assert a.ks == [2., 3.]


@pytest.mark.parametrize('b_values', [
    [(4., 0., 1), (7., math.pi, 2)],
    [(7., math.pi, 2), (4., 0., 1)],
    [(4., .7, 3)],
])
def test_merge_preserves_endpoints_and_is_independent_of_lambda_order(b_values):
    a, b = _molecule(), _molecule()
    _torsion(a, [(2., 0., 1), (3., math.pi, 2)])
    _torsion(b, b_values)
    before_a, before_b = _fingerprint(a), _fingerprint(b)
    stored = {}
    for lam in (0., .5, 1., .25, 1., .5, 0.):
        result = fep.Merge_Force_Field(a, b, lam, {'charge': lam ** 2})
        assert _fingerprint(a) == before_a
        assert _fingerprint(b) == before_b
        values = [_energy(result, phi) for phi in (-2.1, -.4, 0., .8, 2.9)]
        expected = [(1-lam)*_energy(a, phi)+lam*_energy(b, phi) for phi in (-2.1, -.4, 0., .8, 2.9)]
        assert values == pytest.approx(expected)
        for term in result.bonded_forces['dihedral']:
            assert len(term.ks) == len(term.phi0s) == len(term.periodicitys) == term.multiple_numbers
        if lam in stored:
            assert values == pytest.approx(stored[lam])
        stored[lam] = values


def test_differing_bond_lengths_preserve_harmonic_endpoints():
    a, b = _molecule(), _molecule()
    a.bonded_forces['bond'][0].k = 100.
    a.bonded_forces['bond'][0].b = 1.3
    b.bonded_forces['bond'][0].k = 200.
    b.bonded_forces['bond'][0].b = 1.7
    atoms = {atom:i for i,atom in enumerate(a.atoms)}
    target = {atoms[atom] for atom in a.bonded_forces['bond'][0].atoms}
    def energy(mol, distance):
        indices = {atom:i for i,atom in enumerate(mol.atoms)}
        return sum(t.k*(distance-t.b)**2 for t in mol.bonded_forces['bond']
                   if {indices[atom] for atom in t.atoms} == target)
    for lam in (0., .5, 1.):
        result = fep.Merge_Force_Field(a, b, lam)
        for distance in (1., 1.3, 1.7, 2.):
            assert energy(result, distance) == pytest.approx((1-lam)*100*(distance-1.3)**2+lam*200*(distance-1.7)**2)


@pytest.mark.parametrize('missing_endpoint', ['a', 'b'])
def test_endpoint_only_torsions_do_not_append_template_terms(missing_endpoint):
    a, b = _molecule(), _molecule()
    _torsion(a, [(2., .3, 1)])
    _torsion(b, [(4., .7, 3)])
    (a if missing_endpoint == 'a' else b).bonded_forces['dihedral'] = []
    before_a, before_b = _fingerprint(a), _fingerprint(b)
    for lam in (1., .5, 0.):
        result = fep.Merge_Force_Field(a, b, lam)
        for term in result.bonded_forces['dihedral']:
            assert len(term.ks) == len(term.phi0s) == len(term.periodicitys) == term.multiple_numbers
        for phi in (-2., .2, 1.):
            assert _energy(result, phi) == pytest.approx((1-lam)*_energy(a, phi)+lam*_energy(b, phi))
        assert _fingerprint(a) == before_a and _fingerprint(b) == before_b


def test_normal_build_and_serialization_are_unchanged_by_copy(tmp_path):
    mol = _molecule()
    before = _fingerprint(mol)
    duplicate = mol.deepcopy()
    x.save_sponge_input(mol, prefix='original', dirname=str(tmp_path))
    x.save_sponge_input(duplicate, prefix='duplicate', dirname=str(tmp_path))
    assert _fingerprint(mol) == before
    original = {p.name.removeprefix('original_'): p.read_bytes() for p in tmp_path.glob('original_*.txt')}
    copied = {p.name.removeprefix('duplicate_'): p.read_bytes() for p in tmp_path.glob('duplicate_*.txt')}
    assert original == copied
