"""Register terminal phosphate variants with their owning Amber force field."""
from ...helper import source

source("....")
amber = source("...amber")


def register_terminal_monophosphate(family, bases):
    """Keep the HOP3 monoanion separate from the default 5′-OH mapping."""
    load_mol2(os.path.join(AMBER_DATA_DIR, "terminal_monophosphate_" + family + ".mol2"), as_template=True)
    for base in bases:
        res = ResidueType.get_type(base + "5MP")
        res.tail = "O3'"
        res.tail_next = "C3'"
        res.tail_link_conditions.append({"atoms": ["C3'", "O3'"], "parameter": 120 / 180 * np.pi})
        res.tail_link_conditions.append({"atoms": ["H3'", "C3'", "O3'"], "parameter": -54 / 180 * np.pi})
