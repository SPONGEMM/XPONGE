"""
    This **module** gives the unit tests of the functions in Xponge
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
