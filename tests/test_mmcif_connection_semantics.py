"""File connections are optional and only covalent semantics become force-field edges."""
from io import StringIO
import pytest
from test_explicit_residue_link_precedence import x, templates, cif, links


def load(text, **kwargs):
    templates()
    x.ResidueType.get_type('XLF').head = None
    x.ResidueType.get_type('XLF').tail = None
    return x.load_mmcif(StringIO(text), **kwargs)


@pytest.mark.parametrize('kind,expected', [('covale',4),('disulf',4),('hydrog',0),
    ('metalc',0),('nonbonded',0),('unknown',0),('?',0)])
def test_struct_conn_semantics(kind, expected):
    assert len(links(load(cif().replace(' covale ', ' '+kind+' ')))) == expected


@pytest.mark.parametrize('kind,expected', [('covalent',1),('disulfide',1),
    ('coordination',0),('nonbonded_contact',0),('unknown',0)])
def test_mokda_semantics(kind, expected):
    text = cif(explicit=False) + ('loop_\n_mokda_bond_semantic.atom_site_id_1\n'
        '_mokda_bond_semantic.atom_site_id_2\n_mokda_bond_semantic.connection_semantic\n'
        '1 5 '+kind+'\n#\n')
    assert len(links(load(text))) == expected


def test_ignore_file_connections_keeps_external_links():
    text = cif() + ('loop_\n_mokda_bond_semantic.atom_site_id_1\n'
        '_mokda_bond_semantic.atom_site_id_2\n_mokda_bond_semantic.connection_semantic\n'
        '2 6 covalent\n#\nloop_\n_mokda_edit_operation.operation_type\n'
        '_mokda_edit_operation.atom_site_id_1\n_mokda_edit_operation.atom_site_id_2\n'
        '_mokda_edit_operation.connection_semantic\nadd_bond 3 7 covalent\n#\n')
    assert links(load(text, ignore_conect=True)) == set()
    external = [{'atom_a':dict(chain_id='A', residue_seq=1, residue_name='XHB', atom_name='A0'),
                 'atom_b':dict(chain_id='A', residue_seq=5, residue_name='XLF', atom_name='A0')}]
    assert links(load(text, ignore_conect=True, residue_links=external)) == {(0,7)}


def test_noncovalent_edit_does_not_reintroduce_connection():
    text=cif(explicit=False)+('loop_\n_mokda_edit_operation.operation_type\n'
        '_mokda_edit_operation.atom_site_id_1\n_mokda_edit_operation.atom_site_id_2\n'
        '_mokda_edit_operation.connection_semantic\nadd_bond 1 5 nonbonded_contact\n#\n')
    assert links(load(text)) == set()
