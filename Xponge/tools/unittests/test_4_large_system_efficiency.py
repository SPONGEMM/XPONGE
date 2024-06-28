"""
    This **module** gives the unit tests of the speed of SPONGE for large systems
"""

__all__ = ["test_q"]

def test_q():
    """ test the speed of NPT simulation for the Q system"""
    import os
    import Xponge
    from Xponge.mdrun import run

    assert os.system("zenodo_get 10.5281/zenodo.12200626 > download.log 2>&1") == 0
    assert os.system("pmemd.cuda -i amber_nvt.in -p Q.parm7 -c Q.rst7 -o amber_nvt.out -O > amber_nvt.log 2>&1") == 0
    assert os.system("pmemd.cuda -i amber_npt.in -p Q.parm7 -c Q.rst7 -o amber_npt.out -O > amber_npt.log 2>&1") == 0
    assert run("SPONGE -mdin sponge_npt.in -amber_parm7 Q.parm7 -amber_rst7 Q.rst7 \
-default_out_file_prefix sponge_npt > sponge_npt.log 2>&1") == 0
    assert run("SPONGE -mdin sponge_nvt.in -amber_parm7 Q.parm7 -amber_rst7 Q.rst7 \
-default_out_file_prefix sponge_nvt > sponge_nvt.log 2>&1") == 0

