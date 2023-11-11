"""
    This **module** gives the unit tests of the functions in Xponge
"""
import sys
import os
import pathlib
import re
import logging
import importlib.util as iu
import unittest
from ...helper import Xopen, Xprint, GlobalSetting, source

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
    tests = []
    for file_name in file_list:
        result = re.search(r"test_\d_(.+)\.py", file_name)
        if result:
            file_path = os.path.join(module_dir, file_name)
            module_name = result.group(1)
            if todo == "all":
                tests.append(module_name)
            elif todo == module_name:
                spec = iu.spec_from_file_location(module_name, file_path)
                module = iu.module_from_spec(spec)
                spec.loader.exec_module(module)
                for case in module.__all__:
                    tests.append(getattr(module, case))
    return tests

def _run_one_test(case, args):
    """ Run one test"""
    for handle in GlobalSetting.logger.handlers:
        handle.setLevel("CRITICAL")
    log_file = f'{case.__name__}.log'
    file_handler = logging.FileHandler(log_file, mode="w")
    file_handler.set_name("temp")
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s\n%(message)s')
    file_handler.setFormatter(formatter)
    file_handler.setLevel(args.verbose)
    GlobalSetting.logger.addHandler(file_handler)
    runner = XpongeTestRunner()
    result = runner.run(unittest.FunctionTestCase(case))
    for handler in GlobalSetting.logger.handlers:
        if handler.get_name() == "temp":
            GlobalSetting.logger.removeHandler(handler)
        else:
            handler.setLevel(args.verbose)
    return result

def mytest(args):
    """
    This **function** does the tests for Xponge

    :param args: arguments from argparse
    :return: None
    """
    GlobalSetting.logger.setLevel(args.verbose)
    tests = _find_tests(args.do)
    if args.do != "all":
        Xprint(f"{len(tests)} test case(s) for {args.do}")
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
    else:
        Xprint(f"{len(tests)} test script(s)\n{'='*30}")
        for case in tests:
            pathlib.Path(f"{case}").mkdir(exist_ok=True)
            os.system(f"cd {case} & {sys.argv[0]} test -d {case} -v {args.verbose}")
            Xprint("-"*30)
