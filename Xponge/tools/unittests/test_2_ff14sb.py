"""
    This **module** includes unittests of the Xponge.forcefield.amber.ff14sb
"""

__all__ = ["test_ff14sb"]

def test_ff14sb():
    """
        test the single point energy for residues
    """
    import Xponge
    import Xponge.forcefield.amber.ff14sb
    from Xponge.mdrun import run
    import os

    s = "ALA ARG ASN ASP CYS CYX GLN GLU GLY HID HIE HIP ILE LEU " + \
        "LYS MET PHE PRO SER THR TRP TYR VAL HIS"
    mol = Xponge.ResidueType.get_type("ACE")
    for res in s.split():
        mol += Xponge.ResidueType.get_type(res)
    mol += Xponge.ResidueType.get_type("NME")

    Xponge.save_sponge_input(mol, "TEST")
    assert run("SPONGE -mode NVT -default_in_file_prefix TEST " + \
               "-step_limit 5 -write_information_interval 1 " + \
               f"-thermostat middle_langevin -cutoff 8 > {os.devnull}") == 0
