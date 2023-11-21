"""
    This **module** gives the unit tests of the speed of SPONGE minimization
"""

__all__ = ["test_min"]

def test_min():
    """ test the speed of minimization """
    import os
    from Xponge.mdrun import run

    s = r"""@<TRIPOS>MOLECULE
BEN
 12 12 1 0 1
SMALL
USER_CHARGES
@<TRIPOS>ATOM
     1    C   0.0100   1.3950   0.0000   ca         1      BEN   0.000000
     2   C1  -1.2030   0.7060   0.0000   ca         1      BEN   0.000000
     3   C2  -0.0100  -1.3950   0.0000   ca         1      BEN   0.000000
     4   C3  -1.2130  -0.6880  -0.0000   ca         1      BEN   0.000000
     5   C4   1.2030  -0.7060   0.0000   ca         1      BEN   0.000000
     6   C5   1.2130   0.6880   0.0000   ca         1      BEN   0.000000
     7    H   0.0180   2.4810   0.0000   ha            1      BEN   0.000000
     8   H1  -2.1580  -1.2240   0.0000   ha            1      BEN   0.000000
     9   H2  -2.1390   1.2560   0.0000   ha            1      BEN   0.000000
    10   H3   2.1390  -1.2560  -0.0000   ha            1      BEN   0.000000
    11   H4  -0.0180  -2.4810   0.0000   ha            1      BEN   0.000000
    12   H5   2.1580   1.2240   0.0000   ha           1      BEN   0.000000
@<TRIPOS>BOND
     1      1      2 ar
     2      1      3 ar
     3      1      8 1
     4      2      4 ar
     5      2      9 1
     6      3      6 ar
     7      3     10 1
     8      4      5 ar
     9      4      7 1
    10      5      6 ar
    11      5     12 1
    12      6     11 1
@<TRIPOS>SUBSTRUCTURE
    1      BEN      1 ****               0 ****  **** 
"""
    with open("ben.mol2", "w") as f:
        f.write(s)

    with open("leaprc", "w") as f:
        f.write("""source leaprc.water.tip3p
source leaprc.gaff
t = loadmol2 ben.mol2
solvateBox t WAT 20
saveamberparm t t.parm7 t.rst7
quit""")
    with open("mdin", "w") as f:
        f.write("""test efficiency
&cntrl
  imin = 1
  maxcyc = 1000
  ncyc = 500
  ntwx = 500
  ntpr = 500
/
""")
    assert os.system("tleap > tleap.log") == 0
    assert os.system("pmemd.cuda -p t.parm7 -c t.rst7 -i mdin -O -o amber_gpu.out \
> pmemd_cuda.log 2> pmemd_cuda.log") == 0
    assert os.system("pmemd -p t.parm7 -c t.rst7 -i mdin -O -o amber_cpu.out > pmemd.log 2> pmemd.log") == 0
    assert run("SPONGE -cutoff 8 -amber_parm7 t.parm7 -amber_rst7 t.rst7 -mode minimization \
-write_information_interval 500 > sponge.log") == 0
