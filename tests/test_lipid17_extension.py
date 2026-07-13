import json
import os
import subprocess
import sys
from io import StringIO
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "lipid17"
AMBER_DATA = ROOT / "Xponge" / "forcefield" / "amber"
TEMPLATE_NAMES = [
    "2-A", "2-B", "2-C", "CAM", "CLI", "ERG", "H2A", "H2B", "H2C",
    "P2A", "P2B", "P2C", "P2D", "P2E", "P2F", "P3-", "P3A", "P3B",
    "P3C", "P3D", "P3E", "P3F", "P3H", "PC1", "PC2", "PE1", "PE2",
    "PG1", "PG2", "PH3", "PH4", "PH5", "PI", "PI3", "PI4", "PI5",
    "SIT", "STI",
]


def _run_python(code, *args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        [sys.executable, *args],
        input=code,
        text=True,
        capture_output=True,
        cwd=ROOT,
        env=env,
        check=True,
    )


def _pdb_as_mmcif(path):
    tags = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "auth_atom_id",
        "label_comp_id", "auth_comp_id", "label_asym_id", "auth_asym_id",
        "label_seq_id", "auth_seq_id", "pdbx_PDB_ins_code", "label_alt_id",
        "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv",
        "pdbx_PDB_model_num",
    ]
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom = line[12:16].strip()
        residue = line[17:20].strip()
        sequence = int(line[22:26])
        element = line[76:78].strip() or next(char for char in atom if char.isalpha())
        rows.append(
            f"ATOM {int(line[6:11])} {element.upper()} {atom} {atom} {residue} {residue} "
            f"A A {sequence} {sequence} ? . {float(line[30:38])} {float(line[38:46])} "
            f"{float(line[46:54])} 1.0 0.0 1"
        )
    return "data_lipid17\nloop_\n" + "\n".join(f"_atom_site.{tag}" for tag in tags) + "\n" + "\n".join(rows) + "\n#\n"


def test_off_conversion_is_reproducible(tmp_path):
    mol2 = tmp_path / "lipid_ext.mol2"
    manifest = tmp_path / "lipid_ext_manifest.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "convert_amber_off_to_mol2.py"),
            str(DATA / "lipid_ext.lib"),
            str(mol2),
            str(manifest),
        ],
        check=True,
    )
    assert mol2.read_bytes() == (AMBER_DATA / "lipid_ext.mol2").read_bytes()
    assert manifest.read_bytes() == (AMBER_DATA / "lipid_ext_manifest.json").read_bytes()


def test_all_templates_and_connection_metadata_are_registered():
    import Xponge
    import Xponge.forcefield.amber.lipid17  # noqa: F401

    manifest = json.loads((AMBER_DATA / "lipid_ext_manifest.json").read_text())
    assert manifest["template_count"] == 38
    assert [entry["template"] for entry in manifest["templates"]] == TEMPLATE_NAMES

    for entry in manifest["templates"]:
        residue_type = Xponge.ResidueType.get_type(entry["template"])
        assert len(residue_type.atoms) == entry["atom_count"]
        assert abs(sum(atom.charge for atom in residue_type.atoms) - entry["total_charge"]) < 1e-7
        assert residue_type.head == entry["head_atom"]
        assert residue_type.tail == entry["tail_atom"]
        assert residue_type.head_next == entry["head_next_atom"]
        assert residue_type.tail_next == entry["tail_next_atom"]
        Xponge.build_bonded_force(residue_type)


def test_lipid17_loads_glycam_parameters_without_glycam_templates():
    result = _run_python(
        "import json\n"
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid17\n"
        "names = Xponge.ResidueType.get_all_types()\n"
        "print(json.dumps({'PI': 'PI' in names, 'ROH': 'ROH' in names, 'OME': 'OME' in names}))\n"
    )
    state = json.loads(result.stdout.splitlines()[-1])
    assert state == {"PI": True, "ROH": False, "OME": False}


def test_representative_pdbs_load_with_expected_links():
    import Xponge
    import Xponge.forcefield.amber.lipid17  # noqa: F401

    expected = {
        "AHPC": (140, ["AR", "PC", "DHA"]),
        "AHPI": (143, ["AR", "PI", "DHA"]),
        "AHPI3": (146, ["AR", "PI3", "DHA"]),
        "AHPI45H": (151, ["AR", "H2C", "DHA"]),
        "AHPI345H": (155, ["AR", "P3H", "DHA"]),
    }
    for name, (atom_count, residue_names) in expected.items():
        molecule = Xponge.load_pdb(DATA / f"{name}.pdb")
        assert sum(len(residue.atoms) for residue in molecule.residues) == atom_count
        assert [residue.name for residue in molecule.residues] == residue_names
        assert len(molecule.residue_links) == 2
        Xponge.build_bonded_force(molecule)


def test_representative_mmcif_loads_with_fragment_links():
    import Xponge
    import Xponge.forcefield.amber.lipid17  # noqa: F401

    molecule = Xponge.load_mmcif(StringIO(_pdb_as_mmcif(DATA / "AHPI3.pdb")))
    assert [residue.name for residue in molecule.residues] == ["AR", "PI3", "DHA"]
    assert sum(len(residue.atoms) for residue in molecule.residues) == 146
    assert len(molecule.residue_links) == 2
    Xponge.build_bonded_force(molecule)


def test_lysophospholipid_templates_have_exactly_one_linkable_end():
    import Xponge
    import Xponge.forcefield.amber.lipid17  # noqa: F401

    residue_type = Xponge.ResidueType.get_type
    for first, second in [
        ("PA", "PC1"), ("PC2", "PA"),
        ("PA", "PE1"), ("PE2", "PA"),
        ("PA", "PG1"), ("PG2", "PA"),
    ]:
        molecule = residue_type(first) + residue_type(second)
        assert len(molecule.residue_links) == 1
        Xponge.build_bonded_force(molecule)


def test_ff14sb_numeric_lj_values_survive_lipid17_extension_loading():
    result = _run_python(
        "import json\n"
        "import Xponge.forcefield.amber.ff14sb\n"
        "from Xponge.helper import AtomType\n"
        "from Xponge.forcefield.base.lj_base import LJType\n"
        "def values(name):\n"
        "    lj = LJType.get_type(AtomType.get_type(name).LJtype + '-' + AtomType.get_type(name).LJtype)\n"
        "    return [lj.rmin, lj.epsilon]\n"
        "before = {name: values(name) for name in ('N3', 'NT')}\n"
        "import Xponge.forcefield.amber.lipid17\n"
        "after = {name: values(name) for name in ('N3', 'NT')}\n"
        "print(json.dumps([before, after]))\n"
    )
    before, after = json.loads(result.stdout.splitlines()[-1])
    assert after == before
