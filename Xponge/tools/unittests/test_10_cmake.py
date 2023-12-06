"""
Test the cmake building system of Xponge
"""
__all__ = ["test_cmake"]
import os
import shutil
import sys


# check if in the source directory
# get current working directory
def test_cmake():
    # parse sys.argv, get src_dir from -S or --src_dir
    src_found_flag = 0
    for i, arg in enumerate(sys.argv):
        if arg == "-S" or arg == "--src_dir":
            src_dir = sys.argv[i + 1]
            src_found_flag = 1
    if src_found_flag == 0:
        raise ValueError(
            "Source code directory must be provided with -S or --src_dir flag."
        )
    src_dir = os.path.abspath(os.path.expanduser(src_dir))
    cwd = os.getcwd()
    # check if in the source directory
    if not (
        os.path.exists(os.path.join(src_dir, "CMakeLists.txt"))
        and os.path.exists(os.path.join(src_dir, "SPONGE"))
    ):
        raise FileNotFoundError(
            "CMakeLists.txt or the child directory 'SPONGE' not found. Please provide the correct source code directory!"
        )

    # check if cmake is installed
    if shutil.which("cmake") is None:
        raise FileNotFoundError("CMake not found!")

    # build the project
    assert os.system(f"cmake -S {src_dir} -B {cwd}/build") == 0
    assert os.system(f"cmake --build {cwd}/build --config Release -j 12") == 0
