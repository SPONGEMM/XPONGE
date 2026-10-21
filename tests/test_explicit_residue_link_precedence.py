"""Explicit branch connections must not acquire extra implicit head/tail edges."""
from io import StringIO

import pytest
import Xponge as x
import Xponge.forcefield.amber.gaff


def templates():
    for name,count in [('XHB',4),('XLF',1)]:
        atoms='\n'.join(f'{i+1} A{i} {i}.0 0.0 0.0 c3 1 {name} 0.0' for i in range(count))
        x.load_mol2(StringIO(f'@<TRIPOS>MOLECULE\n{name}\n{count} 0 1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n{atoms}\n@<TRIPOS>BOND\n'),as_template=True)
    x.ResidueType.get_type('XLF').head = 'A0'
    x.ResidueType.get_type('XLF').tail = 'A0'


def pdb(record='ATOM',explicit=True):
    lines=[]
    for serial,(name,residue,seq) in enumerate(
        [(f'A{i}','XHB',1) for i in range(4)]+[('A0','XLF',i+2) for i in range(4)],1
    ):
        lines.append(f'{record:<6}{serial:5d} {name:<4} {residue:>3} A{seq:4d}    {float(serial):8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C')
    if explicit:lines.extend(f'CONECT{i+1:5d}{i+5:5d}' for i in range(4))
    return '\n'.join(lines+['END',''])


def cif(record='ATOM',explicit=True):
    columns=['group_PDB','id','type_symbol','label_atom_id','label_comp_id','label_asym_id','label_seq_id','Cartn_x','Cartn_y','Cartn_z']
    rows=[]
    for serial,(name,residue,seq) in enumerate(
        [(f'A{i}','XHB',1) for i in range(4)]+[('A0','XLF',i+2) for i in range(4)],1
    ):rows.append(f'{record} {serial} C {name} {residue} A {seq} {serial}.0 0.0 0.0')
    text='data_links\nloop_\n'+'\n'.join('_atom_site.'+v for v in columns)+'\n'+'\n'.join(rows)+'\n#\n'
    if explicit:
        columns=['id','conn_type_id','ptnr1_label_asym_id','ptnr1_label_seq_id','ptnr1_label_comp_id','ptnr1_label_atom_id','ptnr2_label_asym_id','ptnr2_label_seq_id','ptnr2_label_comp_id','ptnr2_label_atom_id']
        text+='loop_\n'+'\n'.join('_struct_conn.'+v for v in columns)+'\n'
        text+='\n'.join(f'c{i} covale A 1 XHB A{i} A {i+2} XLF A0' for i in range(4))+'\n#\n'
    return text


def links(m):
    if hasattr(m, 'get_atoms'):
        m.get_atoms()
    if isinstance(m.residue_links, set):
        return {tuple(sorted((m.atom_index[link.atom1], m.atom_index[link.atom2]))) for link in m.residue_links}
    return {tuple(sorted(pair)) for pair in m.residue_links}


@pytest.mark.parametrize('record',['ATOM','HETATM'])
@pytest.mark.parametrize('fmt',['pdb','cif'])
def test_explicit_connections_prevent_spurious_branch_cycles(record,fmt):
    templates()
    m=x.load_pdb(StringIO(pdb(record)),ignore_conect=False) if fmt=='pdb' else x.load_mmcif(StringIO(cif(record)))
    m.get_atoms()
    assert len(m.atoms)==8 and len(m.residues)==5
    assert links(m)=={(0,4),(1,5),(2,6),(3,7)}


@pytest.mark.parametrize('fmt',['pdb','cif'])
def test_implicit_multiresidue_chain_without_explicit_edges_is_preserved(fmt):
    templates()
    m=x.load_pdb(StringIO(pdb(explicit=False))) if fmt=='pdb' else x.load_mmcif(StringIO(cif(explicit=False)))
    assert links(m)=={(4,5),(5,6),(6,7)}


def test_ignore_conect_keeps_existing_opt_out_semantics():
    templates()
    assert links(x.load_pdb(StringIO(pdb()),ignore_conect=True))=={(4,5),(5,6),(6,7)}


def test_explicit_ring_closure_does_not_erase_implicit_polymer_chain():
    templates()
    text=pdb(explicit=False).replace('END','CONECT    5    8\nEND')
    assert links(x.load_pdb(StringIO(text),ignore_conect=False))=={(4,5),(5,6),(6,7),(4,7)}


def test_explicit_neighbor_link_is_preserved_not_duplicated():
    templates()
    text=pdb(explicit=False).replace('END','CONECT    5    6\nCONECT    7    8\nEND')
    # Both middle anchors are explicitly linked elsewhere, so no inferred 6--7.
    assert links(x.load_pdb(StringIO(text),ignore_conect=False))=={(4,5),(6,7)}


@pytest.mark.parametrize('boundary', ['TER', 'chain'])
def test_pdb_chain_boundaries_do_not_acquire_template_links(boundary):
    templates()
    rows = pdb(explicit=False).splitlines()
    if boundary == 'TER':
        rows.insert(6, 'TER')
    else:
        rows[6:8] = [row[:21] + 'B' + row[22:] for row in rows[6:8]]
    m = x.load_pdb(StringIO('\n'.join(rows)))
    assert links(m) == {(4, 5), (6, 7)}
