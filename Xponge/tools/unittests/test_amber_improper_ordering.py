"""Deterministic parameter slots and term-wise FEP; no simulation required."""
from types import SimpleNamespace
import math
import numpy as np
import pytest

from Xponge.forcefield.amber import improper_ordering as policy
from Xponge.forcefield.base.dihedral_base import ImproperType
from Xponge.forcefield.special.fep import _improper_merge_rule


@pytest.fixture
def registry(monkeypatch):
    monkeypatch.setattr(ImproperType, '_types', type(ImproperType._types)())
    monkeypatch.setattr(ImproperType, '_types_different_name', type(ImproperType._types_different_name)())
    return ImproperType


def atom(name, typ):
    class Atom:
        pass
    a = Atom();a.name = name;a.type = SimpleNamespace(name=typ)
    return a


def register(registry, rows):
    policy.load_parameters(registry, 'name k phi0 periodicity\n' + '\n'.join(rows))


def test_parameter_slots_and_named_ties_ignore_candidate_and_owner_order(registry):
    register(registry, ['ca-ca-ca-c3 1.1 3.141592653589793 2'])
    a,b,c,d = [atom(n,t) for n,t in [('A','ca'),('B','ca'),('C','ca'),('D','c3')]]
    for candidate in ([d,b,c,a],[b,a,c,d],[a,d,c,b]):
        for order in ([a,b,c,d],[d,c,b,a]):
            selected = policy.match_parameters(registry, [candidate], SimpleNamespace(atoms=order))
            assert [v.name for v in selected[0][1]] == ['A','B','C','D']


def test_override_updates_original_slots_and_preserves_other_periodicities(registry):
    register(registry, ['ca-ca-ca-c3 1.1 0.37 2', 'ca-ca-ca-c3 0.4 0.2 3'])
    register(registry, ['c3-ca-ca-ca 2.2 0.73 2'])
    a,b,c,d = [atom(n,t) for n,t in [('A','ca'),('B','ca'),('C','ca'),('D','c3')]]
    selected = policy.match_parameters(registry, [[a,b,c,d]], SimpleNamespace())
    assert [(r[2]['periodicity'],r[2]['k']) for r in selected] == [(2,2.2),(3,.4)]
    assert [v.name for v in selected[0][1]] == ['D','A','C','B']
    assert [v.name for v in selected[1][1]] == ['A','B','C','D']


def test_specificity_and_wildcard_identity_are_explicit(registry):
    register(registry, ['X-X-ca-X 9 0 2', 'X-X-ca-c3 1.1 0 2'])
    a,b,c,d = [atom(n,t) for n,t in [('B','ca'),('A','ca'),('C','ca'),('D','c3')]]
    selected = policy.match_parameters(registry, [[d,a,c,b]], SimpleNamespace())
    assert selected[0][2]['k'] == 1.1
    assert [v.name for v in selected[0][1]] == ['A','B','C','D']


def test_untagged_types_keep_legacy_path_and_ambiguous_names_reject(registry):
    registry.New_From_String('name k phi0 periodicity\nX-X-ca-c3 1.1 0 2')
    atoms=[atom(n,t) for n,t in [('A','ca'),('A','ca'),('C','ca'),('D','c3')]]
    assert policy.match_parameters(registry,[atoms],SimpleNamespace()) is None
    register(registry,['X-X-ca-c3 1.1 0 2'])
    with pytest.raises(ValueError,match='unique local atom names'):
        policy.match_parameters(registry,[atoms],SimpleNamespace())


def test_generic_override_does_not_inherit_amber_slot_metadata(registry):
    register(registry,['X-X-ca-X 9.0 0 2', 'ca-ca-ca-c3 1.1 0 2'])
    atoms=[atom(n,t) for n,t in [('A','ca'),('B','ca'),('C','ca'),('D','c3')]]
    registry.New_From_String('name k phi0 periodicity\nca-ca-ca-c3 7.7 0.4 2')
    assert policy.match_parameters(registry,[atoms],SimpleNamespace()) is None


def energy(forces, xyz, lookup):
    value=0.
    for f in forces:
        p=xyz[[lookup[a] for a in f.atoms]]
        axis=p[2]-p[1];axis/=np.linalg.norm(axis)
        u=p[0]-p[1];v=p[3]-p[2]
        u-=np.dot(u,axis)*axis;v-=np.dot(v,axis)*axis
        phi=np.arctan2(np.dot(np.cross(axis,u),v),np.dot(u,v))
        value+=f.k*(1+np.cos(f.periodicity*phi-f.phi0))
    return value


@pytest.mark.parametrize('lam',[0.,.25,.5,.75,1.])
def test_fep_preserves_multiple_terms_phases_and_non_equivalent_orders(registry,lam):
    register(registry,['X-X-ca-c3 1.1 0 2'])
    param=next(iter(registry._types.values()))
    atoms=[atom(n,t) for n,t in [('A','ca'),('B','ca'),('C','ca'),('D','c3')]]
    indices={a:i for i,a in enumerate(atoms)}
    def force(ids,n,k,p):
        f=registry.entity([atoms[i] for i in ids],param)
        f.periodicity=n;f.k=k;f.phi0=p
        return f
    a=[force((0,1,2,3),2,1.1,.37),force((0,1,2,3),3,.4,.2)]
    b=[force((0,1,2,3),2,2.2,.73),force((0,3,2,1),3,.8,.4)]
    before=[(f.k,f.phi0,f.periodicity,list(f.atoms)) for f in a+b]
    result=SimpleNamespace(atom_index=indices,bonded_forces={})
    _improper_merge_rule(result,None,None,registry,a,b,lam,None,None,None,{v:v for v in atoms})
    mixed=result.bonded_forces['improper']
    xyz=np.array([[1.,.2,.3],[0,1.,0],[0,0,0],[-1.,-1.,.7]])
    def expected(x):return (1-lam)*energy(a,x,indices)+lam*energy(b,x,indices)
    assert energy(mixed,xyz,indices)==pytest.approx(expected(xyz),abs=1e-12)
    for i in range(4):
        for j in range(3):
            left=xyz.copy();right=xyz.copy();left[i,j]-=1e-5;right[i,j]+=1e-5
            actual=(energy(mixed,right,indices)-energy(mixed,left,indices))/2e-5
            target=(expected(right)-expected(left))/2e-5
            assert actual==pytest.approx(target,abs=1e-8)
    assert before==[(f.k,f.phi0,f.periodicity,list(f.atoms)) for f in a+b]
