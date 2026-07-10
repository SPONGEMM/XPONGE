import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "lipid21"
LIPID17_DATA = Path(__file__).resolve().parent / "data" / "lipid17"
AMBER_DATA = ROOT / "Xponge" / "forcefield" / "amber"
TEMPLATE_NAMES = [
    "AR", "CHL", "DHA", "LAL", "MY", "OL", "PA", "PC", "PE", "PGR",
    "PGS", "PH-", "PS", "SA", "SPM", "ST",
]


def _run(code, *args):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT)
    return subprocess.run(
        [sys.executable, "-c", code, *map(str, args)],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=True,
    )


def test_lipid21_off_conversion_is_reproducible(tmp_path):
    mol2 = tmp_path / "lipid21.mol2"
    manifest = tmp_path / "lipid21_manifest.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "convert_amber_off_to_mol2.py"),
            "--source-license",
            "Public Domain (AmberTools dat/leap)",
            str(DATA / "lipid21.lib"),
            str(mol2),
            str(manifest),
        ],
        check=True,
    )
    assert mol2.read_bytes() == (AMBER_DATA / "lipid21.mol2").read_bytes()
    assert manifest.read_bytes() == (AMBER_DATA / "lipid21_manifest.json").read_bytes()


def test_lipid21_packaged_data_and_special_connection_metadata():
    manifest = json.loads((AMBER_DATA / "lipid21_manifest.json").read_text())
    assert manifest["format_version"] == 2
    assert manifest["source_sha256"] == hashlib.sha256((DATA / "lipid21.lib").read_bytes()).hexdigest()
    assert manifest["source_license"] == "Public Domain (AmberTools dat/leap)"
    assert [entry["template"] for entry in manifest["templates"]] == TEMPLATE_NAMES
    entries = {entry["template"]: entry for entry in manifest["templates"]}
    assert entries["PA"]["head_link_conditions"][1]["parameter_degrees"] == -120.0
    assert entries["PA"]["tail_link_conditions"][1]["parameter_degrees"] == 120.0
    assert entries["PGS"]["head_rule_source"] == "carbonyl_head"
    assert entries["SA"]["head_next_atom"] == "C13"
    assert entries["SA"]["head_reference_atom"] == "H2R"
    assert entries["SPM"]["head_atom"] == "C11"
    assert entries["SPM"]["head_next_atom"] == "N11"
    assert entries["SPM"]["tail_atom"] == "C1"
    assert entries["SPM"]["tail_next_atom"] == "C2"


def test_lipid21_composed_connection_geometry_avoids_fragment_overlap():
    _run(
        "import numpy as np\n"
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid21\n"
        "rt = Xponge.ResidueType.get_type\n"
        "m = rt('PA') + rt('SPM') + rt('SA')\n"
        "distance = lambda a, b: np.linalg.norm(np.array([a.x-b.x, a.y-b.y, a.z-b.z]))\n"
        "assert distance(m.residues[0].name2atom('C13'), m.residues[1].name2atom('C11')) > 2.0\n"
        "assert distance(m.residues[1].name2atom('C1'), m.residues[2].name2atom('C13')) > 2.0\n"
        "m = rt('PGS') + rt('PA')\n"
        "assert distance(m.residues[0].name2atom('C21'), m.residues[1].name2atom('C13')) > 2.0\n"
    )


def test_lipid21_import_registers_base_and_extension_once():
    result = _run(
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid21\n"
        "from Xponge.forcefield.amber._lipid_ext import register_lipid_extension\n"
        "before = len(Xponge.ResidueType.get_type('PI3').head_link_conditions)\n"
        "register_lipid_extension()\n"
        "after = len(Xponge.ResidueType.get_type('PI3').head_link_conditions)\n"
        "assert before == after\n"
    )
    output = result.stdout + result.stderr
    assert output.count("Reference for Lipid21:") == 1
    assert output.count("Reference for Xponge lipid extensions:") == 1
    assert "Reference for Lipid17:" not in output


def test_all_lipid21_and_extension_templates_have_complete_parameters(tmp_path):
    _run(
        "import json, sys\n"
        "from pathlib import Path\n"
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid21\n"
        "base = json.loads(Path(sys.argv[1]).read_text())\n"
        "ext = json.loads(Path(sys.argv[2]).read_text())\n"
        "for entry in base['templates'] + ext['templates']:\n"
        "    residue = Xponge.ResidueType.get_type(entry['template'])\n"
        "    assert len(residue.atoms) == entry['atom_count']\n"
        "    assert abs(sum(atom.charge for atom in residue.atoms) - entry['total_charge']) < 1e-7\n"
        "    Xponge.build_bonded_force(residue)\n",
        AMBER_DATA / "lipid21_manifest.json",
        AMBER_DATA / "lipid_ext_manifest.json",
    )


def test_lipid21_representative_pdbs_and_extension_links(tmp_path):
    _run(
        "import sys\n"
        "from pathlib import Path\n"
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid21\n"
        "data, ext, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])\n"
        "expected = {\n"
        " 'PSM': (127, ['PA','SPM','SA']),\n"
        " 'SSM': (133, ['ST','SPM','SA']),\n"
        " 'PGS': (123, ['PA','PGS','PA']),\n"
        "}\n"
        "for name, (count, residues) in expected.items():\n"
        "    molecule = Xponge.load_pdb(data / (name + '.pdb'))\n"
        "    assert sum(len(residue.atoms) for residue in molecule.residues) == count\n"
        "    assert [residue.name for residue in molecule.residues] == residues\n"
        "    assert len(molecule.residue_links) == 2\n"
        "    Xponge.save_sponge_input(molecule, name, dirname=str(out))\n"
        "molecule = Xponge.load_mmcif(data / 'PSM.cif')\n"
        "assert [residue.name for residue in molecule.residues] == ['PA','SPM','SA']\n"
        "assert sum(len(residue.atoms) for residue in molecule.residues) == 127\n"
        "assert len(molecule.residue_links) == 2\n"
        "Xponge.save_sponge_input(molecule, 'PSM_cif', dirname=str(out))\n"
        "molecule = Xponge.load_pdb(ext / 'AHPI3.pdb')\n"
        "assert [residue.name for residue in molecule.residues] == ['AR','PI3','DHA']\n"
        "assert sum(len(residue.atoms) for residue in molecule.residues) == 146 and len(molecule.residue_links) == 2\n"
        "Xponge.save_sponge_input(molecule, 'AHPI3', dirname=str(out))\n",
        DATA,
        LIPID17_DATA,
        tmp_path,
    )


def test_lipid21_preserves_torsion_specific_nb14_scaling(tmp_path):
    _run(
        "import sys\n"
        "from pathlib import Path\n"
        "import Xponge\n"
        "import Xponge.forcefield.amber.lipid21\n"
        "Xponge.save_sponge_input(Xponge.ResidueType.get_type('PA'), 'pa', dirname=sys.argv[1])\n"
        "rows = (Path(sys.argv[1]) / 'pa_nb14.txt').read_text().splitlines()[1:]\n"
        "scales = {tuple(row.split()[2:]) for row in rows}\n"
        "assert ('0.166667', '0.833333') in scales\n"
        "assert ('0.500000', '0.833333') in scales\n",
        tmp_path,
    )
