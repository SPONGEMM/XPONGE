"""Physical endpoint and serialization regressions shared with the CPP repair."""
import copy
from uuid import uuid4

import pytest
import Xponge as x
import Xponge.forcefield.amber.gaff
from Xponge.forcefield.base import dihedral_base
from Xponge.forcefield.special import fep


def molecule(size=4, star=False):
    rt = x.ResidueType(name="PHYSICAL_" + uuid4().hex)
    for i in range(size):
        rt.Add_Atom("C" + str(i), "c3", float(i), float(i % 2), .2 * i)
        rt.atoms[-1].charge = 0.
    for i in range(1, size):
        rt.Add_Connectivity(rt.atoms[0 if star else i - 1], rt.atoms[i])
    m = x.Molecule("physical")
    m.Add_Residue(rt)
    x.build_bonded_force(m)
    return m


def dual(a, b, mapping, tmp_path):
    indices_b = sorted(mapping)
    path = tmp_path / "mcs.txt"
    path.write_text("mcs atoms in r1: " + str(tuple(mapping[i] for i in indices_b)) +
                    "\nmcs atoms in r2: " + str(tuple(indices_b)) + "\n")
    def assignment(m):
        result = x.Assign()
        rt = m.residues[0].type
        indices = {atom: i for i, atom in enumerate(rt.atoms)}
        for atom in rt.atoms:
            result.add_atom("C", atom.x, atom.y, atom.z, atom.name)
        for atom, neighbors in rt.connectivity.items():
            for other in neighbors:
                if indices[atom] < indices[other]:
                    result.add_bond(indices[atom], indices[other], 1)
        return result
    aa, ab = assignment(a), assignment(b)
    return fep.merge_dual_topology(a, a.residues[0], b.residues[0], aa, ab, imcs=str(path))


@pytest.mark.parametrize("lam", [0., .25, 1.])
def test_one_sided_bonds_only_enter_soft_table(tmp_path, lam):
    a, b = molecule(2), molecule(2)
    left, right, _ = dual(a, b, {0: 0}, tmp_path)
    mixed = fep.merge_force_field(left, right, lam)
    assert not mixed.bonded_forces["bond"]
    terms = mixed.bonded_forces["bond_soft"]
    assert len(terms) == 2
    assert sorted(t.from_AorB for t in terms) == [0, 1]
    assert all(t.k == a.bonded_forces["bond"][0].k for t in terms)
    x.save_sponge_input(mixed, prefix="test", dirname=str(tmp_path))
    assert len((tmp_path / "test_bond_soft.txt").read_text().splitlines()) == 3
    duplicate = mixed.deepcopy()
    assert not duplicate.bonded_forces.get("bond")
    assert len(duplicate.bonded_forces["bond_soft"]) == 2


@pytest.mark.parametrize("mapping", [{0: 0, 1: 3, 2: 1, 3: 2}, {0: 2, 1: 0, 2: 3, 3: 1}])
def test_physical_improper_keeps_ordered_atom_identity(tmp_path, mapping):
    a, b = molecule(star=True), molecule(star=True)
    term = dihedral_base.ImproperType.entity([b.atoms[i] for i in (1, 2, 0, 3)],
                                            dihedral_base.ImproperType.get_type("UNKNOWNS"))
    term.k, term.phi0, term.periodicity = 1.1, 3.141592653589793, 2
    # Use a physical template entity so building a replacement copies this term.
    rt = b.residues[0].type
    template_term = dihedral_base.ImproperType.entity([rt.atoms[i] for i in (1, 2, 0, 3)], term.type)
    template_term.contents = copy.deepcopy(term.contents)
    rt.bonded_forces["improper"] = [template_term]
    b.bonded_forces["improper"] = [term]
    before = copy.deepcopy(term.contents)
    left, right, _ = dual(a, b, mapping, tmp_path)
    mixed = fep.merge_force_field(left, right, 1.)
    active = [t for t in mixed.bonded_forces["improper"] if abs(t.k) > 1e-12]
    assert len(active) == 1
    assert tuple(atom.name for atom in active[0].atoms) == tuple("C" + str(mapping[i]) for i in (1, 2, 0, 3))
    assert active[0].contents == before
    assert term.contents == before


def test_nonzero_pair_at_zero_interpolated_charge_is_rejected():
    a, b = molecule(), molecule()
    a.atoms[0].charge = a.atoms[-1].charge = 1.
    b.atoms[0].charge = b.atoms[-1].charge = -1.
    assert a.bonded_forces["nb14"]
    with pytest.raises(ValueError, match="zero charge product"):
        fep.merge_force_field(a, b, .5)


def test_materialized_append_keeps_explicit_parameters(tmp_path):
    a, b = molecule(), molecule()
    a.bonded_forces["bond"][0].k = 12345.
    b.bonded_forces["bond"][0].k = 54321.
    merged = fep.merge_force_field(a, b, .25)
    before = sorted(t.k for t in merged.bonded_forces["bond"])
    extra = molecule(2)
    joined = merged | extra
    x.save_sponge_input(joined, prefix="joined", dirname=str(tmp_path))
    assert sorted(t.k for t in joined.bonded_forces["bond"]) == sorted(before + [extra.bonded_forces["bond"][0].k])
    assert sorted(t.k for t in merged.bonded_forces["bond"]) == before


def test_residue_force_parameters_do_not_alias_template():
    m = molecule()
    force = m.residues[0].bonded_forces["bond"][0]
    original = m.residues[0].type.bonded_forces["bond"][0].k
    force.k = 98765.
    assert m.residues[0].type.bonded_forces["bond"][0].k == original


def test_softcore_export_converts_nb14_without_mutating_normal_topology(tmp_path):
    m = molecule()
    before = copy.deepcopy([term.contents for term in m.bonded_forces["nb14"]])
    fep.save_soft_core_lj()
    try:
        x.save_sponge_input(m, prefix="soft", dirname=str(tmp_path))
        assert (tmp_path / "soft_LJ_soft_core.txt").is_file()
        assert not (tmp_path / "soft_nb14.txt").exists()
        assert (tmp_path / "soft_nb14_extra.txt").is_file()
        assert [term.contents for term in m.bonded_forces["nb14"]] == before
    finally:
        fep.save_hard_core_lj()
    x.save_sponge_input(m, prefix="hard", dirname=str(tmp_path))
    assert (tmp_path / "hard_LJ.txt").is_file()
    assert (tmp_path / "hard_nb14.txt").is_file()
