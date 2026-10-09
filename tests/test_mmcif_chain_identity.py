from io import StringIO
import pytest
import Xponge as xp
import Xponge.forcefield.amber.ff14sb  # noqa: F401


def structure(chains=("ML1", "ML2", "ML61"), sequences=(1, 2)):
    columns = "group_PDB id type_symbol label_atom_id auth_atom_id label_comp_id auth_comp_id label_asym_id auth_asym_id label_seq_id auth_seq_id pdbx_PDB_ins_code label_alt_id Cartn_x Cartn_y Cartn_z pdbx_PDB_model_num".split()
    rows = ["data_chains", "loop_"] + ["_atom_site." + c for c in columns]
    serial = 0
    for i, chain in enumerate(chains):
        for seq in sequences:
            for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
                serial += 1
                rows.append(f"ATOM {serial} {element} {name} {name} ALA ALA L{i} {chain} {seq} {seq} ? . {serial}.0 0 0 1")
    return "\n".join(rows) + "\n"


def linked_indices(mol):
    result = set()
    for link in mol.residue_links:
        if isinstance(link, (tuple, list)):
            result.add(tuple(sorted(link)))
        elif hasattr(link.atom1, "index"):
            result.add(tuple(sorted((link.atom1.index, link.atom2.index))))
        else:
            atoms = [atom for residue in mol.residues for atom in residue.atoms]
            result.add(tuple(sorted(atoms.index(atom) for atom in (link.atom1, link.atom2))))
    return result


def test_full_chain_identity_and_no_cross_chain_links():
    mol = xp.load_mmcif(StringIO(structure()), infer_terminals=False)
    assert [r.chain_id for r in mol.residues] == ["ML1", "ML1", "ML2", "ML2", "ML61", "ML61"]
    assert linked_indices(mol) == {(2, 4), (10, 12), (18, 20)}


@pytest.mark.parametrize("selector", [
    {"chain_id": "ML2", "residue_seq": 1, "n_terminal": True},
    ("ML2", 1, "", "N"),
])
def test_terminal_selector_preserves_full_chain(selector):
    mol = xp.load_mmcif(StringIO(structure()), infer_terminals=False, terminal_residues=[selector])
    assert [r.name for r in mol.residues] == ["ALA", "ALA", "NALA", "ALA", "ALA", "ALA"]


@pytest.mark.parametrize("selector", ["ML2:1", ("ML2", 1)])
def test_unterminal_selector_preserves_full_chain(selector):
    mol = xp.load_mmcif(StringIO(structure()), unterminal_residues=[selector])
    assert [r.name for r in mol.residues] == ["NALA", "CALA", "ALA", "CALA", "NALA", "CALA"]


def test_external_link_targets_exact_chain():
    def atom(chain, seq, name):
        return dict(chain_id=chain, residue_seq=seq, residue_name="ALA", atom_name=name)
    mol = xp.load_mmcif(StringIO(structure()), infer_terminals=False, residue_links=[
        dict(atom_a=atom("ML1", 2, "C"), atom_b=atom("ML2", 1, "N"))
    ])
    assert linked_indices(mol) == {(2, 4), (10, 12), (18, 20), (6, 8)}


def test_adjacent_residues_with_same_number_remain_separate():
    mol = xp.load_mmcif(StringIO(structure(sequences=(1,))), infer_terminals=False)
    assert [r.chain_id for r in mol.residues] == ["ML1", "ML2", "ML61"]
    assert [len(r.atoms) for r in mol.residues] == [4, 4, 4]
    assert not mol.residue_links


def test_chain_identity_is_case_sensitive():
    mol = xp.load_mmcif(StringIO(structure(chains=("ML", "Ml"))), infer_terminals=False)
    assert [r.chain_id for r in mol.residues] == ["ML", "ML", "Ml", "Ml"]
    assert linked_indices(mol) == {(2, 4), (10, 12)}
