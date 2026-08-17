# Welcome to Use Xponge!

## Introduction

``Xponge`` is a lightweight and easy to customize python package to perform pre- and post-processing of molecular simulations.

### What can Xponge do?

Xponge includes three major categories of functionality, namely, the simulation system construction, simulation data transformation and analysis, and automated workflows for complex simulations. ``Xponge`` is mainly designed for the molecular dynamics (MD) program [SPONGE](https://onlinelibrary.wiley.com/doi/epdf/10.1002/cjoc.202100456)[1], but it can also output some general format files such as mol2 and PDB, so it may help the other molecular modelling programs too.

## Installation

Xponge can be used on all operating systems (Windows/Linux/macOS), including its quantum chemistry workflows. Xponge prefers PySCF and automatically falls back to Psi4 when PySCF is not installed; these workflows report an error if neither backend is available. Until official Windows PySCF wheels are released, build and install a compatible PySCF wheel separately to use the PySCF backend.

### 1. pip install

```bash
pip install Xponge
```

### 2. source setup

- 2.1 Download or clone the source from the gitee or github repository

    The gitee repository is [here](https://gitee.com/gao_hyp_xyj_admin/xponge).
    The github repository is [here](https://github.com/xia-yijie/xponge).

        git clone http://gitee.com/gao_hyp_xyj_admin/xponge.git
        git clone http://github.com/xia-yijie/xponge.git

- 2.2 Open the directory where you download or clone the repository

- 2.3 (Optional) Configure the environment

    It is recommended to use `conda` to configure the environment. Two `yml` files named `install_requirements.yml` and `extras_requirements.yml` are provided in the repository.

    It is recommanded to use the file `install_requirements.yml` to configure the environment. The file will only install the basic dependent packages. If a `ModuleNotFoundError` is raised when you are using `Xponge`, then install the module. This allows you to avoid installing many modules that you will never use, and also makes `Xponge` more cross-platform compatible. Here are the commands to use `install_requirements.yml`.

    ```bash
    conda env create -f install_requirements.yml
    conda activate Xponge
    ```

    All the dependent packages are listed in the [dependencies](#dependencies) section. If you don't want to install the dependent packages one by one (which can be really annoying), the file `extras_requirements.yml` can help you with the environment configuration except the packages `mindspore` and `mindsponge`. The two packages should be installed according to your device (e.g. whether the backend is CPU, GPU or Huawei Ascend Processors) and can not be simply installed by conda. Here are the commands to use `extras_requirements.yml`.

    ```bash
    conda env create -f extras_requirements.yml
    conda activate Xponge
    ```

- 2.4 Run the command

    ```bash
    python setup.py install
    ```

### Installation check

There are some unit tests in ``Xponge``. You can do the basic test to check whether the installation is successful like this:

```bash
Xponge test --do base -o test --verbose 1
```

Here, ``Xponge`` can be replaced to ``python -m Xponge``, ``python3 -m Xponge`` and so on according to your settings of the environmental variables. Some files will be generated after the test is finished.

## Quickstart

Here is a simple example.

```python
import Xponge
# Import the force field you need
import Xponge.forcefield.amber.ff14sb
# Build the molecule like this
peptide = ACE + ALA + NME
# or like this
peptide2 = NALA + ALA * 10 + CALA
# or like this
peptide3 = Xponge.Get_Peptide_From_Sequence("AAAAA")
# See the documentation for more usage!
# Save them as your favorite format
Xponge.save_pdb(peptide, "ala.pdb")
Xponge.save_mol2(peptide2, "ala12.mol2")
Xponge.save_sponge_input(peptide3, "ala5")
```

### Metal force-field overlays

`Xponge.metal_assignment` extends an already loaded and normally assigned
`Molecule`; it is not a separate structure-building pipeline. Callers provide
the locked metal sites, topology mapping, electronic state and parameter
protocol, then continue with the ordinary saver:

```python
import Xponge
import Xponge.forcefield.amber.ff19sb
import Xponge.forcefield.amber.tip3p

mol = Xponge.load_mmcif("system.cif")
result = Xponge.metal_assignment.assign(
    mol,
    package=metal_assignment_package,
    atom_mapping=atom_mapping,
    operation="fit-bonded",
    water_model="tip3p",
    bonded_fit_input=bonded_fit_input,
)
Xponge.save_sponge_input(result.molecule, "system", format="bundle")
```

The API also exposes separate `parameterize` and transactional `apply`
operations. Xponge never splits residues, chooses coordination bonds or builds
caps; the caller must provide a complete topology and stable mapping.

For `manual_bonded`, `bonded_fit_input` carries explicit, hash-closed reference
bond lengths and angles plus site force constants. Xponge validates term
coverage and applies the resulting overlay without running Hessian or
Seminario and without inferring terms from distances. The potential convention
is `E = k * delta^2`; canonical units are
`kcal/mol/angstrom^2` for bonds and `kcal/mol/rad^2` for angles. This path is
intended for experimental engineering validation, not as a claim of
QM-fitted scientific accuracy.

### Amber lipid force fields

Lipid17 and Lipid21 are available as mutually exclusive base lipid force fields:

```python
import Xponge.forcefield.amber.ff14sb   # optional protein family
import Xponge.forcefield.amber.gaff2    # optional small-molecule family
import Xponge.forcefield.amber.lipid21  # or lipid17, but not both
```

Loading either lipid force field also loads the bundled PACKMOL-Memgen-derived
PI/phosphoinositide/LysoPL extension. Lipids use Amber's split residue templates
(for example `PA + SPM + SA`); automatic conversion of full names such as `PSM`
or `POPC` is not provided. The extension combines documented Lipid, GLYCAM,
phosphate, and selected GAFF2-derived terms; its references are printed on load.

The `protein` (`ff14sb`/`ff19sb`), `small_molecule` (`gaff`/`gaff2`), and
`lipid` (`lipid17`/`lipid21`) families each allow one active base force field per
Python process. Use a fresh process to compare alternatives.

Then we can see `ala12.mol2` in VMD:

![pic2](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/2.jpg)

Here is another simple example.

```python
import Xponge
import Xponge.forcefield.amber.tip3p

box = Xponge.BlockRegion(0, 0, 0, 60, 60, 60)
region_1 = Xponge.BlockRegion(0, 0, 20, 20, 20, 40)
region_2 = Xponge.BlockRegion(0, 0, 40, 20, 20, 60)
region_3 = Xponge.BlockRegion(0, 0, 0, 20, 20, 20)
region_4 = Xponge.SphereRegion(20, 10, 30, 10)
region_5 = Xponge.BlockRegion(20, 0, 20, 60, 60, 60)
region_2or3 = Xponge.UnionRegion(region_2, region_3)
region_4and5 = Xponge.IntersectRegion(region_4, region_5)
t = Xponge.Lattice("bcc", basis_molecule=CL, scale=4)
t2 = Xponge.Lattice("fcc", basis_molecule=K, scale=3)
t3 = Xponge.Lattice("sc", basis_molecule=NA, scale=3)
mol = t.Create(box, region_1)
mol = t2.create(box, region_2or3, mol)
mol = t3.create(box, region_4and5, mol)
Xponge.Save_PDB(mol, "out.pdb")
```

Then we can see `out.pdb` in VMD:

![pic1](https://gitee.com/gao_hyp_xyj_admin/xponge/raw/master/README_PICTURE/1.jpg)

## Detailed usage and API documentation

All can be seen [here](https://spongemm.cn/%E6%96%87%E6%A1%A3/Xponge%E6%96%87%E6%A1%A3/Xponge).

## SPONGE bundled inputs

XPONGE can write the topology, protocol, and structural restart HDF5 inputs
directly from a `Molecule`, `Residue`, or `ResidueType`:

```python
import Xponge
import Xponge.forcefield.amber.ff14sb

peptide = Xponge.get_peptide_from_sequence("AAAAA")
Xponge.save_sponge_input(
    peptide,
    prefix="ala5",
    dirname="inputs",
    format="bundle",
)
```

This creates:

```text
inputs/ala5_topology.spgt.h5
inputs/ala5_protocol.spgp.h5
inputs/ala5_restart.spgr.h5
```

Bind them in the SPONGE mdin with `input_h5_topology_path`,
`input_h5_protocol_path`, and `input_h5_restart_path`. The saver does not add
run policy such as mode, step limit, thermostat, or output paths.

Reusable simulation protocol objects can be attached without compiling legacy
`cv.toml` or other text inputs. Object references use stable names and are
validated before any artifact is published:

```python
protocol = Xponge.SpongeProtocol(
    collective_variables=(
        Xponge.ProtocolCollectiveVariable(
            name="distance",
            type="distance",
            atom_indices=(0, 1),
            sigma=(0.5,),
        ),
    ),
    cv_restraints=(
        Xponge.ProtocolCVRestraint(
            name="umbrella",
            cv_refs=("distance",),
            weight=(1.0,),
            reference=(1.5,),
        ),
    ),
    metadynamics=(
        Xponge.ProtocolMetadynamics(
            name="bias",
            cv_refs=("distance",),
            grid_min=(0.0,),
            grid_max=(10.0,),
            grid_count=(64,),
        ),
    ),
    steering=Xponge.ProtocolSteering(
        cv_refs=("distance",),
        weight=(0.25,),
    ),
    hard_wall=Xponge.ProtocolHardWall(
        bounds_low=(0.0, None, None),
        bounds_high=(40.0, None, None),
    ),
)

Xponge.save_sponge_input(
    peptide,
    prefix="ala5",
    dirname="inputs",
    format="bundle",
    protocol=protocol,
)
```

The native object API currently covers scalar CVs, extra distance constraints,
positional restraints, CV harmonic restraints, one enabled metadynamics
object, one static steering definition, and one typed SITS protocol. SITS
method and atom-selection fields are written to `protocol.spgp.h5`; its
versioned launch state (`nk`, optional `log_norm`, and optional `log_nk`) is
written to `restart.spgr.h5`. Named `ProtocolSoftWall` objects are stored as
typed columnar definitions and are compiled directly by SPONGE without an
intermediate configuration file. `ProtocolHardWall` uses `None` for an
unbounded axis, requires at least one finite bound, and defaults to rejecting
NPT unless `allow_npt=True` explicitly opts in. Positional and RMSD-style CV
reference coordinates are written into the structural restart artifact rather
than the reusable protocol artifact. A positional-restraint reference contains
one XYZ row for every system atom; `atom_indices` selects which rows are
restrained.

SITS is represented independently rather than compiled from a legacy section:

```python
sits_protocol = Xponge.SpongeProtocol(
    sits=Xponge.ProtocolSITS(
        mode="production",
        atom_indices=(0, 1),
        temperature_ladder=(280.0, 320.0),
        initial_nk=(1.0, 4.0),
    )
)
```

The default remains the legacy raw-text format. It can be selected explicitly
with ``format="raw"``; the raw implementation is also available directly as
``Xponge.save_sponge_input_raw``. ``Xponge.save_sponge_input_bundle`` remains
available for callers that prefer the format-specific API.

Existing direct/legacy input cases can be converted in either direction:

```bash
Xponge legacy-to-bundle CASE_DIR -m mdin.spg.toml -o CONVERTED
Xponge bundle-to-legacy CONVERTED/bundle -o RESTORED --prefix system
```

Reverse conversion treats typed HDF5 datasets as authoritative. Embedded
legacy sidecars are used only for contracts without a typed representation.
The current reader accepts the `xponge.legacy_to_bundle.v1` input schema and
rejects unknown schema names or versions in strict mode.
Use `--dry-run` to validate and inspect the conversion manifest without
writing output files. Existing targets are preserved unless `--overwrite` is
specified.

Bundled topology and trajectory files are registered as MDAnalysis formats.
Import the XPONGE adapter once, then use the standard `Universe` entry point;
no mdin, protocol, or restart file is required for analysis:

```python
import MDAnalysis as mda
import Xponge.analysis.md_analysis

universe = mda.Universe(
    "topology.spgt.h5",
    "trajectory.spg.h5md",
    topology_format="SPONGE_TOPOLOGY_H5",
    format="SPONGE_H5MD",
    particle_stream="all",
)
```

`Xponge.analysis.md_analysis.load_bundle_universe(topology, trajectory)` is a
stricter convenience wrapper that additionally verifies topology and atom
order hashes. The same `SPONGE_H5MD` reader retains support for historical
`/particles/trajectory` and numbered walker layouts through the `walker`
keyword, and delegates ordinary third-party H5MD files to MDAnalysis' native
reader.

## CLI quick start (trajectory analysis)

The `Xponge traj` subcommand provides cpptraj-like post-analysis for SPONGE trajectories.

Basic syntax:

```
Xponge traj -p TOPO -c TRAJ [-c TRAJ ...] -b BOX [-b BOX ...] [--traj-mode {separate,concat}] -o OUTDIR <analysis> [options]
```

Notes:
- `-c/--traj` and `-b/--box` are repeatable; use `-b` once for a single box file or repeat it to match each trajectory.
- `--traj-mode separate` means multiple `-c` inputs are analyzed independently, with outputs written to per-trajectory subfolders.
- `--traj-mode concat` means multiple `-c` inputs are treated as successive segments of one trajectory and analyzed together in the given order.
- In `-i/--input` command files, `traj=`/`box=` accept comma-separated lists.
- In `-i/--input` command files, `traj_mode=` accepts `separate` or `concat` and defaults to the command-line `--traj-mode` value.

Common analyses:

```
# RMSD
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out rmsd -s "backbone" --dt-ps 1

# RMSF
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out rmsf -s "backbone and name CA"

# Radius of gyration
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out rgyr -s "protein and backbone" --dt-ps 1

# Hydrogen bonds
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out hbond --between "protein,resname SOL" --dt-ps 1

# PCA
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out pca -n 2 -s "backbone"

# Free energy surface (CVs: rmsd:SELECTION, rgyr:SELECTION)
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out fes --cv1 "rmsd:backbone" --cv2 "rgyr:protein and backbone" --bins 20 --temperature 300

# Extract PDBs at specific times (ns)
Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out extract_pdb --times 0 0.1 0.2 --dt-ps 1 -s "all"

# multiple trajectories analyzed separately
Xponge traj -p top.txt -c mdcrd_1.dat -c mdcrd_2.dat -b mdbox_1.txt -b mdbox_2.txt -o out rmsd -s "backbone" --dt-ps 1 --traj-mode separate

# successive trajectory segments concatenated into one analysis
Xponge traj -p top.txt -c md_0_100.dat -c md_100_200.dat -c md_200_300.dat -b box_0_100.txt -b box_100_200.txt -b box_200_300.txt -o out rmsd -s "backbone" --dt-ps 1 --traj-mode concat
```

cpptraj-like command file:

```
# analysis.in
rmsd selection="backbone" dt_ps=1
rmsf selection="backbone and name CA"
rgyr selection="protein and backbone" dt_ps=1
hbond between="protein,resname SOL" dt_ps=1
pca n=2 selection="backbone"
fes cv1="rmsd:backbone" cv2="rgyr:protein and backbone" bins=20 temp=300
extract_pdb times=0,0.1,0.2 selection="all" dt_ps=1

# Optional per-line trajectory/box override (comma-separated for multiple)
rmsd selection="backbone" dt_ps=1 traj="mdcrd_1.dat,mdcrd_2.dat" box="mdbox_1.txt,mdbox_2.txt" traj_mode=concat

Xponge traj -p top.txt -c mdcrd.dat -b mdbox.txt -o out -i analysis.in
```

Outputs are written as PNG figures and JSON data files in `OUTDIR`.

Convert JSON outputs to CSV for tools such as gnuplot:

```bash
# convert one analysis result
Xponge json2csv -i out/rmsd_xxxxxx.json

# choose the output CSV path explicitly
Xponge json2csv -i out/free_energy_surface.json -o out/free_energy_surface.csv

# convert every JSON file in a result directory
Xponge json2csv -i out -o out_csv -r
```

Supported JSON schemas currently include line-series results with `x`/`y`, `extract_pdb` time-file tables,
free-energy surfaces, PCA outputs, and extra sidecar CSVs for fields such as RMSD statistics or raw hydrogen-bond lists.

## Contribution Guideline

If you want to contribute to the main codebase or report some issues, see [here](https://spongemm.cn/zh/%E8%B4%A1%E7%8C%AE) for the guides.

## Dependencies

`Xponge` does not depend on other packages except numpy for its basic use.

However, there are some complicated functions that depend on some other packages. If you do not install the dependent package, you can not use the related functions.

Here is the list of all packages which may be uesd:

| package name      | description                       | how to install                 |
| ------------------| --------------------------------- | ------------------------------ |
| XpongeLib         | c/c++ compiled library for Xponge | `pip install XpongeLib`        |
| pyscf [2-4]       | quantum chemistry                 | `pip install pyscf`            |
| geometric[5]      | geometry optimization             | `pip install geometric`        |
| rdkit[6]          | cheminformatics                   | `pip install rdkit`            |
| MDAnalysis[7-8]   | trajectory analysis               | `pip install MDAnalysis`       |
| matplotlib        | plot and visualization            | `pip install matplotlib`       |
| mindspore[9]      | AI framework for machine learning | See the [official website](https://www.mindspore.cn/install)|
| mindsponge[1]     | end-to-end differentiable MD      | See the [official website](https://www.mindspore.cn/mindscience/docs/en/master/mindsponge/intro_and_install.html)|

## References

[0] Y. Xia, Y. Q. Gao, *J. Open Source Softw.* (2022) DOI:[10.21105/joss.04467](https://doi.org/10.21105/joss.04467)

[1] Y.-P. Huang, et al. *Chinese J. Chem.* (2022) DOI: [10.1002/cjoc.202100456](https://doi.org/10.1002/cjoc.202100456)

[2] Q. Sun, et al. *J. Chem. Phys.* (2020) DOI: [10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

[3] Q. Sun, et al. Wiley Interdiscip. *Rev. Comput. Mol. Sci.* (2018) DOI: [10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340)

[4] Q. Sun, *J. Comp. Chem.* (2015) DOI: [10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981)

[5] L.-P. Wang, C.C. Song, *J. Chem. Phys.* (2016) DOI: [10.1063/1.4952956](https://doi.org/10.1063/1.4952956)

[6] RDKit: Open-source cheminformatics. https://www.rdkit.org

[7] R. J. Gowers, et al. Proceedings of the 15th Python in Science Conference (2016) DOI: [10.25080/majora-629e541a-00e](https://doi.org/10.25080/majora-629e541a-00e)

[8] N. Michaud-Agrawal, et al. *J. Comput. Chem.* (2011) DOI: [10.1002/jcc.21787](https://10.1002/jcc.21787)

[9] MindSpore: An Open AI Framwork. https://www.mindspore.cn/
