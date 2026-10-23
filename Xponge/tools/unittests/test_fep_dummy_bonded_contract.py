"""Ordinary dual-topology branches keep auxiliary bonded terms at both endpoints."""
import pytest
from .test_9_fep import _load_self_contained_merge


@pytest.mark.parametrize('lam', [0., .25, .5, .75, 1.])
def test_non_mcs_bonded_terms_do_not_fade(lam):
    from Xponge.forcefield.special import fep
    _, left, right, _, _ = _load_self_contained_merge()
    branch = lambda force: any(atom.name.endswith('R2') for atom in force.atoms)
    for key in ('bond', 'angle', 'dihedral', 'improper'):
        a = [f for f in left.bonded_forces.get(key, []) if branch(f)]
        b = [f for f in right.bonded_forces.get(key, []) if branch(f)]
        assert len(a) == len(b), key
    mixed = fep.merge_force_field(left, right, lam)
    assert not mixed.bonded_forces.get('bond_soft')
    def signature(force):
        atoms = tuple(atom.name for atom in force.atoms)
        parameters = []
        for field in ('k', 'b', 'phi0', 'periodicity', 'ks', 'phi0s', 'periodicitys'):
            value = getattr(force, field, None)
            if value is not None:
                parameters.append((field, tuple(value) if isinstance(value, list) else value))
        return min(atoms, atoms[::-1]), tuple(parameters)
    for key in ('bond', 'angle', 'dihedral', 'improper'):
        expected = sorted(signature(f) for f in right.bonded_forces.get(key, []) if branch(f))
        assert sorted(signature(f) for f in mixed.bonded_forces.get(key, []) if branch(f)) == expected


def test_two_exclusive_branches_do_not_create_cross_forces(tmp_path):
    import Xponge as x
    from Xponge.forcefield.special import fep
    from .test_9_fep import _build_smiles_residue

    # Two terminal heavy substituents share only the first two heavy atoms.
    # Freeze the MCS rather than allowing a new search to remove this boundary.
    aa, _, a = _build_smiles_residue('CCO', 'BRANCH_A')
    bb, _, b = _build_smiles_residue('CCN', 'BRANCH_B')
    mapping = tmp_path / 'mcs.in'
    mapping.write_text('mcs atoms in r1: (0, 1)\nmcs atoms in r2: (0, 1)\n')
    left, right, _ = fep.merge_dual_topology(a, a.residues[0], b.residues[0], aa, bb, imcs=str(mapping))
    a_only = {atom.name for atom in a.atoms[2:]}
    for lam in (0., .25, .5, .75, 1.):
        mixed = fep.merge_force_field(left, right, lam)
        assert not mixed.bonded_forces.get('bond_soft')
        for kind, forces in mixed.bonded_forces.items():
            for force in forces:
                atom_names = {atom.name for atom in force.atoms}
                assert not (atom_names & a_only and any(n.endswith('R2') for n in atom_names)), kind
        for atom in mixed.atoms:
            if atom.name in a_only:
                excluded = {other.name for other in atom.extra_excluded_atoms}
                assert {other.name for other in mixed.atoms if other.name.endswith('R2')} <= excluded
