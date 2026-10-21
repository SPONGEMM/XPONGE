import pytest
import Xponge as xp
from Xponge.forcefield.base import angle_base, bond_base, dihedral_base


def test_live_registry_query_includes_wildcards_without_mutation():
    classes = (bond_base.BondType, angle_base.AngleType,
               dihedral_base.ProperType, dihedral_base.ImproperType)
    snapshots = [(cls, dict(cls._types), dict(cls._types_different_name)) for cls in classes]
    try:
        dihedral_base.ImproperType(name="X-X-Qcenter-Qoxygen", k=10.5, phi0=0., periodicity=2)
        assert xp.amber_parameter_exists("IMPROPER", ["Qalpha", "qn", "Qcenter", "Qoxygen"])
        assert not xp.amber_parameter_exists("IMPROPER", ["Qalpha", "qn", "qcenter", "Qoxygen"])
        bond_base.BondType(name="QA-qb", k=100., b=1.5)
        assert xp.amber_parameter_exists("BOND", ["qb", "QA"])
        assert not xp.amber_parameter_exists("BOND", ["QA", "QB"])
        angle_base.AngleType(name="QA-QB-qc", k=30., b=1.)
        assert xp.amber_parameter_exists("ANGLE", ["qc", "QB", "QA"])
        dihedral_base.ProperType(name="X-QA-QB-X", k=1., phi0=0., periodicity=3)
        assert xp.amber_parameter_exists("DIHE", ["qa", "QA", "QB", "qd"])
        before = [dict(cls._types) for cls in classes]
        assert not xp.amber_parameter_exists("IMPROPER", ["missing"] * 4)
        assert before == [dict(cls._types) for cls in classes]
        with pytest.raises(ValueError):
            xp.amber_parameter_exists("BOND", ["QA"])
    finally:
        for cls, types, aliases in snapshots:
            cls._types.clear()
            cls._types.update(types)
            cls._types_different_name.clear()
            cls._types_different_name.update(aliases)
