"""Shared residue-connection configuration for Amber lipid force fields."""

import json

import numpy as np

from ...helper import ResidueType


def _append_unique(conditions, condition):
    if condition not in conditions:
        conditions.append(condition)


def configure_connection(
    residue_name,
    position,
    anchor,
    next_atom,
    angle_atoms,
    angle_degrees,
    dihedral_atoms,
    dihedral_degrees,
    length=1.5,
):
    """Configure one residue connection end without duplicating conditions."""
    residue = ResidueType.get_type(residue_name)
    setattr(residue, position, anchor)
    setattr(residue, f"{position}_next", next_atom)
    setattr(residue, f"{position}_length", length)
    conditions = getattr(residue, f"{position}_link_conditions")
    _append_unique(
        conditions,
        {"atoms": list(angle_atoms), "parameter": angle_degrees / 180 * np.pi},
    )
    _append_unique(
        conditions,
        {"atoms": list(dihedral_atoms), "parameter": dihedral_degrees / 180 * np.pi},
    )


def configure_standard_chain(residue_name):
    """Configure a Lipid17-style saturated/unsaturated chain at C12."""
    for position in ("head", "tail"):
        configure_connection(
            residue_name,
            position,
            "C12",
            "C13",
            ["H2R", "C12"],
            109.5,
            ["H2S", "H2R", "C12"],
            -120 if position == "head" else 120,
        )


def configure_standard_headgroup(residue_name):
    """Configure a two-ended glycerophospholipid headgroup."""
    configure_connection(
        residue_name, "head", "C11", "O11", ["O11", "C11"], 120,
        ["O12", "O11", "C11"], 180,
    )
    configure_connection(
        residue_name, "tail", "C21", "O21", ["O21", "C21"], 120,
        ["O22", "O21", "C21"], 180,
    )


def configure_manifest(path):
    """Apply connection metadata from a generated lipid manifest."""
    with open(path, encoding="utf-8") as manifest_file:
        manifest = json.load(manifest_file)
    for entry in manifest["templates"]:
        for position in ("head", "tail"):
            anchor = entry[f"{position}_atom"]
            if anchor is None:
                continue
            link_conditions = entry.get(f"{position}_link_conditions")
            if link_conditions:
                angle, dihedral = link_conditions
                configure_connection(
                    entry["template"],
                    position,
                    anchor,
                    entry[f"{position}_next_atom"],
                    angle["atoms"],
                    angle["parameter_degrees"],
                    dihedral["atoms"],
                    dihedral["parameter_degrees"],
                )
            else:
                next_atom = entry[f"{position}_next_atom"]
                reference = entry[f"{position}_reference_atom"]
                configure_connection(
                    entry["template"], position, anchor, next_atom,
                    [next_atom, anchor], 120,
                    [reference, next_atom, anchor], 180,
                )
    return manifest
