"""
    test the workflow of FEP
"""

__all__ = ["test_uncovalent", "test_covalent"]

def test_uncovalent():
    """ test whether the program is able to run """
    from subprocess import Popen, PIPE
    import Xponge
    from Xponge import GlobalSetting, Xprint
    import Xponge.forcefield.amber.tip3p
    import Xponge.forcefield.amber.gaff
    if GlobalSetting.purpose == "academic":
        return

    t = Xponge.get_assignment_from_smiles("C")
    t.determine_atom_type("gaff")
    t.calculate_charge("tpacm4")
    a = t.to_residuetype("A")
    wat = Xponge.ResidueType.get_type("WAT")
    t = Xponge.add_solvent_box(a, wat, 10)
    Xponge.save_pdb(t, "test.pdb")
    Xponge.save_mol2(a)
    t = Xponge.get_assignment_from_smiles("CC")
    t.determine_atom_type("gaff")
    t.calculate_charge("tpacm4")
    Xponge.save_mol2(t.to_residuetype("B"))
    with Popen("Xponge mol2rfe -pdb test.pdb -r1 A.mol2 -r2 B.mol2 -nl 1 -p1step 5000 -estep 5000".split(),
                stdout=PIPE, stderr=PIPE, stdin=PIPE) as p:
        outs, hints = p.communicate()
        Xprint(hints.decode("utf-8"))
        Xprint(outs.decode("utf-8"))
        assert p.returncode == 0

def test_covalent():
    """ test whether the program is able to run """
    from subprocess import Popen, PIPE
    import Xponge
    from Xponge import GlobalSetting, Xprint
    import Xponge.forcefield.amber.tip3p
    import Xponge.forcefield.amber.ff14sb
    if GlobalSetting.purpose == "academic":
        return
    ace = Xponge.ResidueType.get_type("ACE")
    wat = Xponge.ResidueType.get_type("WAT")
    ala = Xponge.ResidueType.get_type("ALA")
    gly = Xponge.ResidueType.get_type("GLY")
    nme = Xponge.ResidueType.get_type("NME")
    t = ace + ala + nme
    t = Xponge.add_solvent_box(t, wat, 10)
    Xponge.save_pdb(t, "test.pdb")
    Xponge.save_mol2(ala, "A.mol2")
    Xponge.save_mol2(gly, "B.mol2")
    with Popen("Xponge mol2rfe -pdb test.pdb -r1 A.mol2 -r2 B.mol2 -nl 1 -p1step 5000 -estep 5000 -ri 1".split(),
         stdout=PIPE, stderr=PIPE, stdin=PIPE) as p:
        outs, hints = p.communicate()
        Xprint(hints.decode("utf-8"))
        Xprint(outs.decode("utf-8"))
        assert p.returncode == 0
