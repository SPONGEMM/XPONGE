"""
This **module** defines the terminal commands
"""
import argparse
from . import tools

def _mytest(subparsers):
    """

    :param subparsers:
    :return:
    """
    mytest = subparsers.add_parser("test", help="test the basic function of Xponge")
    mytest.add_argument("-p", "--purpose", metavar="programmatic", default="programmatic",
                        choices=["programmatic", "academic"],
                        help="the purpose of the unittests. a test with programmatic purpose will \
only test the program, and a test with academic purpose will do more to exam the system")
    mytest.add_argument("-v", "--verbose", metavar="INFO", default="INFO",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                        help="the verbose level for output, 0")
    mytest.add_argument("-d", "--do", metavar="todo",
                        default="base",
                        help="the unit tests need to do")
    mytest.add_argument("-f", "--file", metavar="file",
                        help="the unit test file. If set, argument -d/--do will be suppressed")
    mytest.set_defaults(func=tools.mytest)


def _converter(subparsers):
    """

    :param subparsers:
    :return:
    """
    converter = subparsers.add_parser("converter", help="convert the format of coordinate file")
    converter.add_argument("-p", required=True, metavar="TOP", help="the name of the topology file")
    converter.add_argument("-c", metavar="CRD", help="the name of the coordinate file")
    converter.add_argument("-o", required=True, metavar="OUT", help="the name of the output file")
    converter.add_argument("-cf", metavar="GUESS", default="guess", choices=["guess", "sponge_crd", "sponge_traj"],
                           help='''the format of the topology file, can be "guess", "sponge_crd" or "sponge_traj"''')
    converter.add_argument("-of", metavar="GUESS", default="guess", choices=["guess", "sponge_crd", "sponge_traj"],
                           help='''the format of the output file, can be "guess", "sponge_crd" or "sponge_traj"''')
    converter.set_defaults(func=tools.converter)


def _maskgen(subparsers):
    """

    :param subparsers:
    :return:
    """
    maskgen = subparsers.add_parser("maskgen", help="""use MDAnalysis to generate a file to record the atom
 indexes (and coordinates ) of the corresponding mask""")
    maskgen.add_argument("-p", required=True, help="the topology file")
    maskgen.add_argument("-pf", help="the topology file format")
    maskgen.add_argument("-c", help="the coordinate file")
    maskgen.add_argument("-cf", help="the coordinate file format")
    maskgen.add_argument("-b", help="the box file, which is required for the SPONGE trajectory file")
    maskgen.add_argument("-o", required=True, help="the output file")
    maskgen.add_argument("-s", required=True, help="the select string")
    maskgen.add_argument("-oc", help="the output coordinate file")
    maskgen.set_defaults(func=tools.maskgen)


def _exgen(subparsers):
    """

    :param subparsers:
    :return:
    """
    exgen = subparsers.add_parser("exgen",
                                  help='process bond-like, angle-like, dihedral-like files to get the atoms to exclude')
    exgen.add_argument('-n', type=int, required=True, help='the atom numbers')
    exgen.add_argument('-o', required=True, help='output exclude file name')
    exgen.add_argument('-b', '--bond', default=[], nargs='+', help='''bond-like input files: skip the first line,
 and there are 2 atoms in the head of following lines''')
    exgen.add_argument('-a', '--angle', default=[], nargs='+', help='''angle-like input files: skip the first line,
 and there are 3 atoms in the head of following lines''')
    exgen.add_argument('-d', '--dihedral', default=[], nargs='+', help='''dihedral-like input files:
 skip the first line, and there are 4 atoms in the head of following lines''')
    exgen.add_argument('-v', '--virtual', default=[], nargs='+', help='''virtual-atom-like input files:
 the first number indicates the virtual type''')
    exgen.add_argument('-e', '--exclude', default=[], nargs='+', help='''exclude-like input files:
 add the information of another exclude file''')
    exgen.set_defaults(func=tools.exgen)


def _name2name(subparsers):
    """

    :param subparsers:
    :return:
    """
    name2name = subparsers.add_parser("name2name",
                                      help="change the atom names of a residue from one file to another file")
    name2name.add_argument("-fformat", "-from_format", dest="from_format", choices=["mol2", "pdb", "gaff_mol2"],
                           required=True, help="the format of the file which is needed to change from")
    name2name.add_argument("-ffile", "-from_file", dest="from_file", required=True,
                           help="the name of the file which is needed to change from")
    name2name.add_argument("-fres", "-from_residue", dest="from_residue", default="",
                           help="the residue name in ffile if fformat == pdb")

    name2name.add_argument("-tformat", "-to_format", dest="to_format", choices=["mol2", "pdb", "gaff_mol2"],
                           required=True, help="the format of the file which is needed to change to")
    name2name.add_argument("-tfile", "-to_file", dest="to_file", required=True,
                           help="the name of the file which is needed to change to")
    name2name.add_argument("-tres", "-to_residue", dest="to_residue", default="",
                           help="the residue name in tfile if tformat == pdb")

    name2name.add_argument("-oformat", "-out_format", dest="out_format", choices=["mol2", "pdb", "mcs_pdb"],
                           required=True, help="the format of the output file")
    name2name.add_argument("-ofile", "-out_file", dest="out_file", required=True, help="the name of the output file")
    name2name.add_argument("-ores", "-out_residue", dest="out_residue", default="ASN",
                           help="the name of the output residue")
    name2name.add_argument("-tmcs", type=int, default=10, help="the time to find max common structure")
    name2name.set_defaults(func=tools.name2name)


def _mol2rfe(subparsers):
    """

    :param subparsers:
    :return:
    """
    mol2rfe = subparsers.add_parser("mol2rfe",
                                    help='calculate the relative binding energy of a small molecule using SPONGE')
    mol2rfe.add_argument("-do", metavar="todo", nargs="*", action="append", help="""the things need to do,
 should be one or more of 'build', 'min', 'pre_equilibrium', 'equilibrium', 'analysis', 'debug'""",
                         choices=["build", "min", "pre_equilibrium", "equilibrium", "analysis", "debug"])

    mol2rfe.add_argument("-pdb", required=True, help="the initial conformation given by the pdb file")
    mol2rfe.add_argument("-r2", "-residuetype2", required=True,
                         help="molecule mutated to by an Xponge ResidueType mol2 file")
    mol2rfe.add_argument("-r1", "-residuetype1", required=True,
                         help="molecule mutated from by an Xponge ResidueType mol2 file")
    mol2rfe.add_argument("-r0", "-residuetype0", nargs="*", default=[],
                         help="small molecules or non-standard residues that are also in the system")
    mol2rfe.add_argument("-ri", "-residue_index", type=int, metavar=0, default=0,
                         help="the residue index of the molecule to mutate")

    mol2rfe.add_argument("-fl", "-lambda_in_file", type=str, default=None,
                         help="Specify lambda values in file. \
                            Lambda numbers will be counted in this file. \
                            This option overrides the -nl and -l option.")
    mol2rfe.add_argument("-l", "-lambda", type=float, nargs="+", default=None,
                         help="Specify lambda values in command line. \
                            Lambda numbers will be counted in this file. \
                            This option overrides the -nl option and is overrided by \
                            the -fl option.")
    mol2rfe.add_argument("-nl", "-lambda_numbers", metavar=20, type=int, default=20,
                         help="the number of lambda groups - 1, default 20 for 0, 0.05, 0.10, 0.15..., 1.0")
    mol2rfe.add_argument("-cp", "-charge_power", metavar=2, type=int, default=2,
                         help="The power of charge as a function of the lambda")

    mol2rfe.add_argument("-dohmr", "-do_hydrogen_mass_repartition", action="store_true",
                         help="use the hydrogen mass repartition method")
    mol2rfe.add_argument("-ff", "-forcefield", help="Use this force field file instead of the default ff14SB and gaff")
    mol2rfe.add_argument("-mi", "-min_mdin", nargs="*", help="Use this minimization mdin file \
instead of the default one")
    mol2rfe.add_argument("-pi", "-pre_equilibrium_mdin", help="Use this pre-equilibrium mdin file \
instead of the default one")
    mol2rfe.add_argument("-ei", "-equilibrium_mdin", help="Use this equilibrium mdin file instead of the default one")
    mol2rfe.add_argument("-ai", "-analysis_mdin", help="Use this analysis mdin file instead of the default one")
    mol2rfe.add_argument("-method", default="TI", choices=["TI", "MBAR"], help="the method to analyze the free energy")
    mol2rfe.add_argument("-temp", default="TMP", metavar="TMP", help="the temporary file name prefix")

    mol2rfe.add_argument("-tmcs", default=60, type=int, metavar="60",
                         help="the timeout parameter for max common structure in unit of second")
    mol2rfe.add_argument("-fmcs", default="mcs.png", type=str, metavar="mcs.png",
                         help="the file name for max common structure image")
    mol2rfe.add_argument("-imcs", type=str,
                         help="the input max common structure file name")
    mol2rfe.add_argument("-lmcs", default=0.0, type=float, metavar="0.0",
                         help="minimum limitation of the Tanimoto coefficient of max common structure")
    mol2rfe.add_argument("-dt", default=2e-3, type=float, metavar="dt",
                         help="the dt used for simulation when mdin is not provided")
    mol2rfe.add_argument("-mstep", "-min_step", dest="min_step", default=5000, type=int,
                         metavar="5000",
                         help="the minimization step used for simulation")
    mol2rfe.add_argument("-mpy", help="the python script to do after minimization")
    mol2rfe.add_argument("-pstep", "-pre_equilibrium_step", dest="pre_equilibrium_step", default=50000, type=int,
                         metavar="50000",
                         help="the pre-equilibrium step used for simulation")
    mol2rfe.add_argument("-ppy", help="the python script to do after pre_equilibrium")
    mol2rfe.add_argument("-p1step", "-first_pre_equilibrium_step", dest="p1step", default=100000, type=int,
                         metavar="100000",
                         help="the equilibrium step used for the first lambda simulation")
    mol2rfe.add_argument("-estep", "-equilibrium_step", dest="equilibrium_step", default=500000, type=int,
                         metavar="500000",
                         help="the equilibrium step used for simulation")
    mol2rfe.add_argument("-epy", help="the python script to do after equilibrium")
    mol2rfe.add_argument("-wi", default=100, metavar="100", type=int,
                         help="the writing information interval for equilibrium simulations")
    mol2rfe.add_argument("-thermostat", default="middle_langevin",
                         metavar="middle_langevin", choices=["middle_langevin"],
                         help="the thermostat used for simulation when mdin is not provided")
    mol2rfe.add_argument("-barostat", default="andersen_barostat",
                         metavar="andersen_barostat", choices=["andersen_barostat"],
                         help="the barostat used for simulation when mdin is not provided")
    mol2rfe.add_argument("-device", default=0, type=int,
                         metavar="0",
                         help="the index of the cuda device to use")
    mol2rfe.set_defaults(func=tools.mol2rfe)


def _mm_gbsa(subparsers):
    """

    :param subparsers:
    :return:
    """
    mm_gbsa = subparsers.add_parser("mmgbsa",
                                    help='calculate the absolute binding energy using SPONGE via MM/GBSA')
    mm_gbsa.add_argument("-do", metavar="todo", nargs="*", action="append", help="""the things need to do,
 should be one or more of 'build', 'min', 'pre_equilibrium', 'equilibrium', 'analysis', 'debug'""",
                         choices=["build", "min", "pre_equilibrium", "equilibrium", "analysis", "debug"])

    mm_gbsa.add_argument("-pdb", required=True, help="the initial conformation given by the pdb file")
    mm_gbsa.add_argument("-s1", default="resid 1", metavar='"resid 1"',
                         help="the MDAnalysis selection of the first part")
    mm_gbsa.add_argument("-s2", default="not resid 1 and protein", metavar='"not resid 1 and protein"',
                         help="the MDAnalysis selection of the second part")
    mm_gbsa.add_argument("-sr", help="the MDAnalysis selection of the restrained part")
    mm_gbsa.add_argument("-rw", help="the restraint weight", type=float, default=5, metavar="5")
    mm_gbsa.add_argument("-r0", "-residuetypes", nargs="*", default=[],
                         help="small molecules or non-standard residues in the system given by the Xponge mol2 file(s)")
    mm_gbsa.add_argument("-dohmr", "-do_hydrogen_mass_repartition", action="store_true",
                         help="use the hydrogen mass repartition method")
    mm_gbsa.add_argument("-ff", "-forcefield", help="Use this force field file instead of the default ff14SB and gaff")
    mm_gbsa.add_argument("-mi", "-min_mdin", nargs="*", help="Use this minimization mdin file \
instead of the default one")
    mm_gbsa.add_argument("-pi", "-pre_equilibrium_mdin", help="Use this pre-equilibrium mdin file \
instead of the default one")
    mm_gbsa.add_argument("-ei", "-equilibrium_mdin", help="Use this equilibrium mdin file instead of the default one")
    mm_gbsa.add_argument("-ai", "-analysis_mdin", help="Use this analysis mdin file instead of the default one")
    mm_gbsa.add_argument("-temp", default="TMP", metavar="TMP", help="the temporary file name prefix")

    mm_gbsa.add_argument("-dt", default=2e-3, type=float, metavar="dt",
                         help="the dt used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-mstep", "-min_step", dest="min_step", default=5000, type=int,
                         metavar="5000",
                         help="the minimization step used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-pstep", "-pre_equilibrium_step", dest="pre_equilibrium_step", default=50000, type=int,
                         metavar="50000",
                         help="the pre-equilibrium step used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-estep", "-equilibrium_step", dest="equilibrium_step", default=500000, type=int,
                         metavar="500000",
                         help="the equilibrium step used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-astart", "-analysis_start", default=200, type=int,
                         metavar="200", help="the start frame for analysis")
    mm_gbsa.add_argument("-astop", "-analysis_stop", default=-1, type=int,
                         metavar="-1", help="the stop frame for analysis")
    mm_gbsa.add_argument("-astride", "-analysis_stride", default=500, type=int,
                         metavar="500", help="the stride frame for analysis")
    mm_gbsa.add_argument("-thermostat", default="middle_langevin",
                         metavar="middle_langevin", choices=["middle_langevin"],
                         help="the thermostat used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-barostat", default="andersen_barostat",
                         metavar="andersen_barostat", choices=["andersen_barostat"],
                         help="the barostat used for simulation when mdin is not provided")
    mm_gbsa.add_argument("-ivacuum", nargs=3, type=float,
                         metavar="x y z",
                         help="the input initial vacuum layer thickness in unit of Angstrom")
    mm_gbsa.add_argument("-device", default=0, type=int,
                         metavar="0",
                         help="the index of the cuda device to use")
    mm_gbsa.add_argument("-nvt", action="store_true", help="use NVT instead of NPT ensemble")

    mm_gbsa.set_defaults(func=tools.mm_gbsa)


def main():
    """

    :return:
    """

    parser = argparse.ArgumentParser(prog="Xponge")
    parser.add_argument("-v", "--version", action="store_true", help="show the version of Xponge")
    subparsers = parser.add_subparsers(help="subcommands",
                                       description="Tools for SPONGE. Use Xponge XXX -h for the help of tool 'XXX'.")
    _mytest(subparsers)
    _maskgen(subparsers)
    _exgen(subparsers)
    _name2name(subparsers)
    _mol2rfe(subparsers)
    _converter(subparsers)
    _mm_gbsa(subparsers)

    args = parser.parse_args()

    if args.version:
        from . import __version__
        print(__version__)
    elif hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
