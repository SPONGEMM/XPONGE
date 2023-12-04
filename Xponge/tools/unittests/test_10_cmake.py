"""
Test the cmake building system of Xponge
"""
__all__ = ["test_cmake"]
import os
import shutil
import re


# check if in the source directory
# get current working directory
def test_cmake():
    cwd = os.getcwd()
    # for test target "all", the test will create ./cmake/cmake so the source code is in ../..
    if cwd.endswith("cmake"):
        # source code in parent directory
        src_dir = cwd.replace("/cmake/cmake", "")
    # for test target "cmake", we create a new directory ./cmake for test outputs
    else:
        src_dir = cwd
        os.mkdir("cmake")
        os.chdir("cmake")
    # check if in the source directory
    if not (
        os.path.exists(os.path.join(src_dir, "CMakeLists.txt"))
        and os.path.exists(os.path.join(src_dir, "SPONGE"))
    ):
        raise ValueError("Please run this test in the source directory")

    # check if cmake is installed
    if shutil.which("cmake") is None:
        raise ValueError("cmake is not installed")

    # build the project
    os.system(f"cmake -S {src_dir} -B {src_dir}/build")
    os.system(f"cmake --build {src_dir}/build --config Release -j 12")
    # add cmake output to the path
    os.environ["PATH"] += os.pathsep + os.path.join(src_dir, "build", "SPONGE")
    # get all tests
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
            if module_name != "cmake":
                tests.append(module_name)
    # run Xponge tests
    for test in tests:
        os.mkdir(test)
        os.chdir(test)
        os.system(f"Xponge test -d {test} -v INFO")
        os.chdir("..")
