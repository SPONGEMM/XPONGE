"""
This **module** gives the unit tests of the functions in Xponge

`Xponge test` only gives the tests for the basic functions of Xponge

`Xponge test --do all` starts all of the tests

`Xponge test --do XXX` starts the tests for XXX

Requirements:

- conda install python==3.11
- pip install numpy
- pip install rdkit-pypi
- pip install xpongelib
- pip install netcdf4
- pip install pubchempy
- pip install MDAnalysis (for nojump)
- pip install pyscf (for RESP calculation only)
- pip install zenodo_get
- vmd (for manully result check)
- AMBER (version: 24)
- GROMACS (version:24.2)
    * charmm36-jul2022.ff
    * oplsaam.ff

The tests include:

0. base: run Xponge
    - base: run Xponge
        * hint: print this hint
        * import: import the Xponge forcefield module
        * assign: create a simple water Xponge.assign
        * molecule: create a simple methane Xponge.molecule
1. building: build systems
    - assign: Xponge.assign
        * get_assign: get an assignment
        * residuetype_to_assign: convert Xponge.ResidueType to Xponge.Assign
        * objective_assign_rule: create a rule to assign atom types (Xponge.AtomType)
        * string_assign_rule: create a rule to assign atom types (str)
        * atom_deletion: delete an atom from Xponge.assign
        * bond_deletion: delete a bond from Xponge.assign
        * equal_atom_search: search the equalvalent atoms in an Xponge.assign object
        * atom_type_determination: determine the atom type
        * uff_optimization: optimize the coordinates using UFF
        * saving: save the object
    - charge: calculate partial charge
        * tpacm4
        * resp (todo)
        * gasteiger (todo)
    - cif: process CIF files
        * cof: covalent organic framework
        * som: small organic molecule (todo)
    - helper: help to process files and Xponge objects
        * adding_of_residues: add and multiply Xponge.ResidueType and Xponge.Moleule
        * omitting_residue_type: omit atoms for a residue type
        * reseting_type: change the type of certain residue
        * pdb_filter: clean the PDB file (Xponge.pdb_filter)
    - io: input and output
        * pdb_general: load a PDB file
        * pdb_ssbond_link_and_conect: read the conect and the ssbond part of a PDB file
        * mol2_general: load a MOL2 file
    - lattice: create a lattice
    - process: run Xponge.process
        * impose: impose the bond, angle or dihedral of a molecule
        * solvent_process: add the cubic solvent box, rotate the main axis and hydrogen mass repartition
        * octbox: add the truncated octahedron solvent box (todo)
    - sasa: calculate the solvent accessible surface area
2. forcefield_loading: load forcefield
    - ffbase: load basic parameters
3. forcefield_using: use forcefield
    - bsc1
    - charmm27
    - charmm36
    - ff14sb
    - ff19sb
    - ol3
    - oplsaam
    - mW
4. MD_efficiency: test the efficiency of the simulations
    - min_efficiency: minimize a benzene molecule with a random-generated comformation
    - nvt_efficiency: simulate in NVT of 6915 atoms (waters)
    - npt_efficiency: simulate in NPT of 6915 atoms (waters)
    - large_system_efficiency (todo): simulate in NPT of 500,000 atoms
5. MD_function: test extra functions for MD
    - nopbc: simulate the system without periodic boundary condition
        * gb: simulate general born model
        * cv_run: steer MD simulation without pbc
    - rerun: rerun the trajectory
    - wall: reflect the atoms
        * soft_wall: with an extra position bias potential
        * hard_wall: reflect the velocity
    - wrap: wrap the molecules
        * periodic: there is one molecule which has bonds across box
        * unperiodic: there is no molecule which has bonds across box
    - PMC-IZ (todo)
6. MD_thermodynamics: test the thermodynamic properties
    - andersen_barostat
        * reweighting
    - berendsen_barostat
        * avarage_value
    - bussi_barostat
        * reweighting
    - monte_carlo_barostat
        * reweighting
    - andersen_thermostat
        * reweighting
    - berendsen_thermostat
        * avarage_value
    - bussi_thermotat
        * reweighting (todo)
    - nose_hoover_chain
        * reweighting
    - middle_langevin
        * reweighting
    - langevin
        * reweighting
    - rdf
7. MD_kinetics: test the kinetic properties
    - msd
    - tip4p
8. enhancing sampling:
    - metadynamics
        * 1d_sampling
        * 2d_sampling (todo)
    - umbrella
    - tmd
    - sits
        * 1d_sampling
        * 2d_sampling
9. workflow:
    - fep
        * uncovalent
        * covalent (todo)

Some results you may check it manually:

0. base
    - base
        * hint: read `test_hint.log`
1. building
    - assign:
        * saving: `vmd -m ben.pdb ben.mol2`
    - cif:
        * cof: `vmd -sponge_mass ./cof_mass.txt ./cof.dat`
          cof: `vmd -sponge_mass ./cof2_mass.txt ./cof2.dat`
        * som: `vmd -sponge_mass ./som_mass.txt ./som.dat`
    - helper:
        * pdb_filter: read and visualize `test_filter.pdb`
    - lattice:
        * all: `vmd out.pdb`
    - sasa:
        * surface: `vmd -m protein.pdb surface.xyz`
3. forcefield_using:
    - all: check the png files
    - mw: read `sponge/test.out` and `lammps/lammps.log`
4. MD_efficiency:
    - min_efficiency: read `amber_cpu.out`, `amber_gpu.out` and `sponge.log`
    - nvt_efficiency: see `mdout` and `mdinfo.txt`
    - npt_efficiency: see `mdout` and `mdinfo.txt`
    - large_system_efficiency: see `amber_nvt.out`, `amber_npt.out`, `sponge_nvt.info` and `sponge_npt.info`
5. MD_function:
    - nopbc:
        * gb: check the png files
        * cv_run: `vmd -sponge_mass cv_mass.txt ./mdcrd.dat`
    - wall:
        * hard_wall: `vmd -sponge_mass wats_mass.txt ./hard.dat`
        * soft_wall: `vmd -sponge_mass wats_mass.txt ./soft.dat`
    - wrap:
        * periodic: `vmd -sponge_mass p_mass.txt ./p.dat`
        * unperiodic: `vmd -sponge_mass up_mass.txt ./up.dat`
6. MD_thermodynamics
    - andersen_barostat: see the png file, and read `test_reweighting.log`
    - berendsen_barostat: read `test_average_value.log`
    - bussi_barostat: see the png file, and read `test_reweighting.log`
    - monte_carlo_barostat: 
    - andersen_thermostat: see the png file
    - nose_hoover_chain: see the png file
    - middle_langevin: see the png file
    - langevin: see the png file
    - rdf: see the png files
7. MD_kinetics
    - msd: see the png file
    - tip4p: see the png file
8. enhancing_sampling
    - tmd: `vmd -sponge_mass ./trp_cage_mass.txt -sponge_crd ./trp_cage_coordinate.txt ./mdcrd.dat` and compare the RMSD value (selection for VMD: protein and name C CA N O)
    - umbrella: see the png files
    - metadynamics: see the png files
    - sits: see the png files (todo)
9. workflow
    - fep (todo)
"""

import sys
import os
import pathlib
import warnings
import re
import logging
import importlib.util as iu
import unittest
from ...helper import Xdict, Xopen, Xprint, GlobalSetting, source
warnings.filterwarnings("ignore")

CATEGORY = Xdict({'0': "base",
                  '1': "building",
                  '2': "forcefield_loading",
                  '3': "forcefield_using",
                  '4': "MD_efficiency",
                  '5': "MD_function",
                  '6': "MD_thermodynamics",
                  '7': "MD_kinetics",
                  '8': "enhancing_sampling",
                  '9': "workflow",
                  '100': "application"},
                  not_found_message="{} is not a valid unittest category")

class XpongeTestRunner(unittest.TextTestRunner):
    """ the unittest wrapper of Xponge tests """
    def run(self, test):
        result = self._makeResult()
        test(result)
        if result.errors:
            for error in result.errors:
                Xprint(error[1], "ERROR")
        if result.failures:
            for error in result.failures:
                Xprint(error[1], "ERROR")
        return result

def _find_tests(todo):
    """ find all tests in the folder"""
    module_dir = os.path.dirname(os.path.abspath(__file__))
    file_list = os.listdir(module_dir)
    file_list.sort()
    tests = []
    for file_name in file_list:
        result = re.search(r"test_(\d+)_(.+)\.py", file_name)
        if result:
            file_path = os.path.join(module_dir, file_name)
            index = result.group(1)
            module_name = result.group(2)
            if todo == "all":
                tests.append([module_name, file_path, index])
            elif todo == module_name:
                spec = iu.spec_from_file_location(module_name, file_path)
                module = iu.module_from_spec(spec)
                spec.loader.exec_module(module)
                tests.append([])
                for case in module.__all__:
                    tests[0].append(getattr(module, case))
                tests.append(CATEGORY[index])
    return tests

def _check_test_file(f):
    """ check the cases in the test file """
    if not os.path.exists(f):
        raise ValueError(f"{f} does not exist")
    result = re.search(r"test_(\d+)_(.+)\.py", f)
    category = CATEGORY[result.group(1)]
    module_name = result.group(2)
    spec = iu.spec_from_file_location(module_name, f)
    module = iu.module_from_spec(spec)
    spec.loader.exec_module(module)
    tests = []
    for case in module.__all__:
        tests.append(getattr(module, case))
    return tests, module_name, category

def _run_one_test(case, verbose):
    """ Run one test"""
    for handle in GlobalSetting.logger.handlers:
        handle.setLevel("CRITICAL")
    log_file = f'{case.__name__}.log'
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.set_name("temp")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s\n%(message)s')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(verbose)
    GlobalSetting.logger.addHandler(file_handler)
    runner = XpongeTestRunner()
    result = runner.run(unittest.FunctionTestCase(case))
    for handler in GlobalSetting.logger.handlers:
        if handler.get_name() == "temp":
            GlobalSetting.logger.removeHandler(handler)
        else:
            handler.setLevel(verbose)
    return result

def _run_several_tests(tests, name, args, catogory):
    """run several tests"""
    Xprint(f"{len(tests)} test case(s) for {catogory} - {name}")
    failures = []
    errors = []
    for case in tests:
        result = _run_one_test(case, args)
        if result.failures:
            failures.append(case.__name__[5:])
        if result.errors:
            errors.append(case.__name__[5:])
    if failures:
        Xprint(f"\nFailed function(s): {', '.join(failures)}")
    if errors:
        Xprint(f"\nError function(s): {', '.join(errors)}")
    if not failures and not errors:
        Xprint("")

def mytest(args):
    """
    This **function** does the tests for Xponge

    :param args: arguments from argparse
    :return: None
    """
    GlobalSetting.logger.setLevel(args.verbose)
    GlobalSetting.purpose = args.purpose
    if args.file:
        tests, name, category = _check_test_file(args.file)
        _run_several_tests(tests, name, args.verbose, category)
    elif args.do != "all":
        tests = _find_tests(args.do)
        if not tests:
            raise ValueError(f"No test named {args.do} found")
        _run_several_tests(tests[0], args.do, args.verbose, tests[1])
    else:
        tests = _find_tests(args.do)
        Xprint(f"{len(tests)} test script(s)\n{'='*30}")
        for case, f, index in tests:
            folder = pathlib.Path(f"{CATEGORY[index]}") / f"{case}"
            folder.mkdir(exist_ok=True, parents=True)
            os.system(f"cd {folder} && {sys.argv[0]} test -f {f} -v {args.verbose} -p {args.purpose}")
            Xprint("-"*30)
