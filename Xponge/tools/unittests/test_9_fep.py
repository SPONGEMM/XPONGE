"""
    test the workflow of FEP
"""

__all__ = ["test_uncovalent", "test_covalent"]

import pytest
from uuid import uuid4


def _find_force_by_indices(forces, atom2index, indices, same_force):
    target = tuple(indices)
    for force in forces:
        if callable(atom2index):
            force_indices = [atom2index(atom) for atom in force.atoms]
        else:
            force_indices = [atom2index[atom] for atom in force.atoms]
        for permutation in same_force(force_indices):
            if tuple(permutation) == target:
                return force
    raise AssertionError(f"force with indices {indices} not found")


def _unique_name(prefix):
    return f"{prefix}_{uuid4().hex[:8]}"


def _build_smiles_residue(smiles, prefix):
    import Xponge
    from Xponge import build
    import Xponge.forcefield.amber.gaff

    assign_ = Xponge.get_assignment_from_smiles(smiles)
    assign_.determine_atom_type("gaff")
    assign_.calculate_charge("tpacm4")
    residue_type = assign_.to_residuetype(_unique_name(prefix))
    build.build_bonded_force(residue_type)
    mol = Xponge.Molecule(_unique_name(f"{prefix}_mol"))
    mol.add_residue(residue_type)
    return assign_, residue_type, mol


def _build_self_contained_merge_inputs(zero_b_extra_charge=False):
    import Xponge
    from Xponge import build

    from_assign, from_res_type, from_mol = _build_smiles_residue("CCO", "FEP_A")
    to_assign, to_res_type, to_mol = _build_smiles_residue("CCON", "FEP_B")
    if zero_b_extra_charge:
        to_res_type.atoms[3].charge = 0.0
        build.build_bonded_force(to_res_type)
        to_mol = Xponge.Molecule(_unique_name("FEP_B_mol"))
        to_mol.add_residue(to_res_type)
    return Xponge, from_mol, to_mol, from_assign, to_assign


def _load_self_contained_merge(zero_b_extra_charge=False):
    from Xponge.forcefield.special import fep

    Xponge, from_mol, to_mol, from_assign, to_assign = _build_self_contained_merge_inputs(zero_b_extra_charge)

    merged_from, merged_to, matchmap = fep.Merge_Dual_Topology(
        from_mol, from_mol.residues[0], to_mol.residues[0], from_assign, to_assign)
    return Xponge, merged_from, merged_to, to_mol.residues[0], matchmap


def _find_b_only_heavy_nb14_pair(forces):
    for force in forces:
        if sum(atom.mass > 2 for atom in force.atoms) != 2:
            continue
        if any(atom.name.endswith("R2") for atom in force.atoms):
            return force
    raise AssertionError("B-only heavy-heavy nb14 pair not found")

def test_uncovalent():
    """ test whether the program is able to run """
    import sys
    from subprocess import Popen, PIPE
    import Xponge
    from Xponge import GlobalSetting, Xprint
    import Xponge.forcefield.amber.tip3p
    import Xponge.forcefield.amber.gaff
    if getattr(GlobalSetting, "purpose", None) == "academic":
        return

    t = Xponge.get_assignment_from_smiles("C")
    t.determine_atom_type("gaff")
    t.calculate_charge("tpacm4")
    a = t.to_residuetype("A")
    wat = Xponge.ResidueType.get_type("WAT")
    t = Xponge.add_solvent_box(a, wat, 10)
    Xponge.save_pdb(t, "test.pdb")
    Xponge.save_mol2(a)
    t = Xponge.get_assignment_from_smiles("CC")
    t.determine_atom_type("gaff")
    t.calculate_charge("tpacm4")
    Xponge.save_mol2(t.to_residuetype("B"))
    with Popen([sys.executable, "-m", "Xponge", "mol2rfe", "-do", "build", "-pdb", "test.pdb", "-r1", "A.mol2",
                "-r2", "B.mol2", "-nl", "1", "-p1step", "5000", "-estep", "5000"],
               stdout=PIPE, stderr=PIPE, stdin=PIPE) as p:
        outs, hints = p.communicate()
        Xprint(hints.decode("utf-8"))
        Xprint(outs.decode("utf-8"))
        assert p.returncode == 0

def test_covalent():
    """ test whether the program is able to run """
    import sys
    from subprocess import Popen, PIPE
    import Xponge
    from Xponge import GlobalSetting, Xprint
    import Xponge.forcefield.amber.tip3p
    import Xponge.forcefield.amber.ff14sb
    if getattr(GlobalSetting, "purpose", None) == "academic":
        return
    ace = Xponge.ResidueType.get_type("ACE")
    wat = Xponge.ResidueType.get_type("WAT")
    ala = Xponge.ResidueType.get_type("ALA")
    gly = Xponge.ResidueType.get_type("GLY")
    nme = Xponge.ResidueType.get_type("NME")
    t = ace + ala + nme
    t = Xponge.add_solvent_box(t, wat, 10)
    Xponge.save_pdb(t, "test.pdb")
    Xponge.save_mol2(ala, "A.mol2")
    Xponge.save_mol2(gly, "B.mol2")
    with Popen([sys.executable, "-m", "Xponge", "mol2rfe", "-do", "build", "-pdb", "test.pdb", "-r1", "A.mol2",
                "-r2", "B.mol2", "-nl", "1", "-p1step", "5000", "-estep", "5000", "-ri", "1"],
               stdout=PIPE, stderr=PIPE, stdin=PIPE) as p:
        outs, hints = p.communicate()
        Xprint(hints.decode("utf-8"))
        Xprint(outs.decode("utf-8"))
        assert p.returncode == 0


def test_merge_dual_topology_switches_common_forces_to_b_endpoint():
    from Xponge.forcefield.base import bond_base, angle_base

    _, merged_from, merged_to, to_res, _ = _load_self_contained_merge()

    res_from = merged_from.residues[0].type
    res_to = merged_to.residues[0].type
    bond_from = _find_force_by_indices(res_from.bonded_forces["bond"], res_from.atom2index, (1, 2),
                                       bond_base.BondType.Same_Force)
    bond_to = _find_force_by_indices(res_to.bonded_forces["bond"], res_to.atom2index, (1, 2),
                                     bond_base.BondType.Same_Force)
    angle_to = _find_force_by_indices(res_to.bonded_forces["angle"], res_to.atom2index, (0, 1, 2),
                                      angle_base.AngleType.Same_Force)
    expected_bond = _find_force_by_indices(
        to_res.type.bonded_forces["bond"], to_res.type.atom2index,
        (1, 2),
        bond_base.BondType.Same_Force)
    expected_angle = _find_force_by_indices(
        to_res.type.bonded_forces["angle"], to_res.type.atom2index,
        (0, 1, 2),
        angle_base.AngleType.Same_Force)

    assert res_to.atoms[2].type.name == to_res.type.atoms[2].type.name
    assert bond_from.type.name != bond_to.type.name
    assert bond_to.type.name == expected_bond.type.name
    assert bond_to.k == pytest.approx(expected_bond.k)
    assert bond_to.b == pytest.approx(expected_bond.b)
    assert angle_to.type.name == expected_angle.type.name
    assert angle_to.k == pytest.approx(expected_angle.k)
    assert angle_to.b == pytest.approx(expected_angle.b)


def test_b_only_nb14_uses_zero_charge_endpoint():
    from Xponge.forcefield.special import fep

    _, merged_from, merged_to, _, _ = _load_self_contained_merge(zero_b_extra_charge=True)

    merged_1 = fep.Merge_Force_Field(merged_from, merged_to, 1.0, {"charge": 1.0})
    merged_02 = fep.Merge_Force_Field(merged_from, merged_to, 0.2, {"charge": 0.2})

    pair_1 = _find_b_only_heavy_nb14_pair(merged_1.bonded_forces["nb14_extra"])
    pair_02 = _find_b_only_heavy_nb14_pair(merged_02.bonded_forces["nb14_extra"])

    assert pair_1.A > 0
    assert pair_1.B > 0
    assert pair_1.kee == pytest.approx(0.0)
    assert pair_02.kee == pytest.approx(0.0)
    assert pair_02.A == pytest.approx(pair_1.A * 0.2)
    assert pair_02.B == pytest.approx(pair_1.B * 0.2)
    assert any(atom.name.endswith("R2") and atom.charge == pytest.approx(0.0) for atom in pair_1.atoms)


def test_merge_dual_topology_reuses_process_without_restype_name_collision():
    from Xponge.forcefield.special import fep

    _, from_mol, to_mol, from_assign, to_assign = _build_self_contained_merge_inputs()
    canonical_ab_name = f"{from_mol.residues[0].type.name}_{to_mol.residues[0].type.name}"
    canonical_ba_name = f"{to_mol.residues[0].type.name}_{from_mol.residues[0].type.name}"

    first_from, first_to, _ = fep.Merge_Dual_Topology(
        from_mol, from_mol.residues[0], to_mol.residues[0], from_assign, to_assign)
    second_from, second_to, _ = fep.Merge_Dual_Topology(
        from_mol, from_mol.residues[0], to_mol.residues[0], from_assign, to_assign)

    assert first_from.residues[0].type.name != second_from.residues[0].type.name
    assert first_to.residues[0].type.name != second_to.residues[0].type.name
    assert "__FEP_TMP_" in first_from.residues[0].type.name
    assert "__FEP_TMP_" in first_to.residues[0].type.name
    assert first_from.residues[0].name == canonical_ab_name
    assert first_to.residues[0].name == canonical_ba_name
    assert second_from.residues[0].name == canonical_ab_name
    assert second_to.residues[0].name == canonical_ba_name


def test_merge_force_field_preserves_canonical_residue_names():
    from Xponge.forcefield.special import fep

    _, from_mol, to_mol, from_assign, to_assign = _build_self_contained_merge_inputs()
    canonical_ab_name = f"{from_mol.residues[0].type.name}_{to_mol.residues[0].type.name}"

    merged_from, merged_to, _ = fep.Merge_Dual_Topology(
        from_mol, from_mol.residues[0], to_mol.residues[0], from_assign, to_assign)
    merged_lambda = fep.Merge_Force_Field(merged_from, merged_to, 0.5, {"charge": 0.5})

    assert merged_lambda.residues[0].name == canonical_ab_name
    assert "__FEP_TMP_" not in merged_lambda.residues[0].name
