"""Native protocol CV objects to self-contained legacy CV sections."""

from __future__ import annotations

import re

import numpy as np

from .errors import BundleExportError
from .topology_exporters import _format_scalar, _text_vector


_PROTOCOL = "protocol.spgp.h5"
_RESTART = "restart.spgr.h5"
_RESERVED = {
    "CV_type", "atom", "coordinate", "period", "sigma", "rotate",
    "function", "min_padding", "max_padding",
}


def _name(value, path):
    if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", value):
        raise BundleExportError(f"{path} cannot be represented as a legacy CV name")
    return value


def _scalar(reader, path, default=None):
    if not reader.contains(_PROTOCOL, path) and default is not None:
        return default
    value = reader.read(_PROTOCOL, path)
    if value.shape != ():
        raise BundleExportError(f"{path} must be scalar")
    return value.item()


def _tokens(value, path):
    value = np.asarray(value)
    if value.size == 0:
        raise BundleExportError(f"{path} must not be empty")
    if value.dtype.kind in "SUO":
        if value.ndim > 1:
            raise BundleExportError(f"{path} string data must be scalar or 1D")
        text = " ".join(_text_vector(value))
    elif value.dtype.kind in "biuf":
        text = " ".join(_format_scalar(item) for item in value.reshape(-1))
    else:
        raise BundleExportError(f"{path} has unsupported data type {value.dtype}")
    if not text.strip() or any(char in text for char in "\r\n\x00{}"):
        raise BundleExportError(f"{path} cannot be represented as a legacy CV value")
    return text


def _enabled(reader, root):
    value = _scalar(reader, root + "/enabled_default", 1)
    if value not in (0, 1):
        raise BundleExportError(f"{root}/enabled_default must be 0 or 1")
    return bool(value)


def _selection(reader, root, atom_count, virtual_names):
    indices = reader.contains(_PROTOCOL, root + "/atom_indices")
    refs = reader.contains(_PROTOCOL, root + "/atom_refs")
    if indices and refs:
        raise BundleExportError(f"{root} atom_indices and atom_refs are mutually exclusive")
    if not indices and not refs:
        if reader.contains(_PROTOCOL, root + "/selection_expression"):
            raise BundleExportError(f"{root}/selection_expression must be resolved before export")
        return []
    path = root + ("/atom_indices" if indices else "/atom_refs")
    values = reader.read(_PROTOCOL, path)
    if values.ndim != 1 or (indices and values.dtype.kind not in "iu"):
        raise BundleExportError(f"{path} must be a 1D atom selection")
    result = _text_vector(values)
    for value in result:
        if refs and value in virtual_names:
            continue
        if not re.fullmatch(r"[0-9]+", value) or not 0 <= int(value) < atom_count:
            raise BundleExportError(f"{path} contains invalid atom reference {value!r}")
    if len(result) != len(set(result)):
        raise BundleExportError(f"{path} contains duplicate atom references")
    return result


def _reference(reader, root, cv_type, atom_count):
    name = root.rsplit("/", 1)[-1]
    sources = (
        (_PROTOCOL, root + "/coordinate"),
        (_RESTART, f"/parameters/restart/references/cv/{name}/coordinate"),
    )
    reference = None
    for bundle_file, path in sources:
        if not reader.contains(bundle_file, path):
            continue
        if cv_type != "rmsd":
            raise BundleExportError(f"{path} is only supported for rmsd CV objects")
        values = reader.read(bundle_file, path)
        if values.dtype.kind not in "iuf" or values.shape != (atom_count, 3):
            raise BundleExportError(f"{path} must have shape ({atom_count}, 3) and numeric values")
        # SPONGE consumes reference coordinates as float32 in both H5 paths.
        with np.errstate(over="ignore", invalid="ignore"):
            values = values.astype(np.float32)
        if not np.isfinite(values).all():
            raise BundleExportError(f"{path} must contain only finite float32 values")
        if reference is not None and not np.array_equal(reference, values):
            raise BundleExportError(f"{path} conflicts with the inline CV reference")
        reference = values
    if cv_type == "rmsd" and (atom_count == 0 or reference is None):
        raise BundleExportError(f"{root} requires atoms and an RMSD reference coordinate dataset")
    return reference


def native_cv_sections(reader):
    """Return enabled native virtual atoms and CVs in runtime section form."""
    roots = [
        (name, "/cv/" + name)
        for name in reader.list_children(_PROTOCOL, "/cv", groups_only=True)
        if name not in {"config", "virtual_atom"}
    ]
    virtual_roots = [
        (name, "/cv/virtual_atom/" + name)
        for name in reader.list_children(_PROTOCOL, "/cv/virtual_atom", groups_only=True)
    ]
    roots = [(name, root) for name, root in roots if _enabled(reader, root)]
    virtual_roots = [(name, root) for name, root in virtual_roots if _enabled(reader, root)]
    if not roots and not virtual_roots:
        return []
    atom_count = int(reader.read_scalar("topology.spgt.h5", "/topology/atom_count"))
    virtual_names = {name for name, _ in virtual_roots}
    sections = []
    for name, root in virtual_roots:
        _name(name, root)
        kind = _tokens(_scalar(reader, root + "/type"), root + "/type")
        if kind not in {"center", "center_of_mass"}:
            raise BundleExportError(f"{root} has unsupported virtual atom type {kind!r}")
        atoms = _selection(reader, root, atom_count, set())
        if not atoms:
            raise BundleExportError(f"{root} requires atom_indices")
        items = {"vatom_type": kind, "atom": " ".join(atoms)}
        if reader.contains(_PROTOCOL, root + "/weight"):
            weights = reader.read(_PROTOCOL, root + "/weight")
            if kind != "center" or weights.shape != (len(atoms),):
                raise BundleExportError(f"{root}/weight does not match the virtual atom definition")
            items["weight"] = _tokens(weights, root + "/weight")
        elif kind == "center":
            raise BundleExportError(f"{root}/weight is required for center")
        sections.append((name, items))
    for name, root in roots:
        _name(name, root)
        if name in virtual_names:
            raise BundleExportError(f"{root} conflicts with a virtual atom name")
        if _scalar(reader, root + "/dimension", 1) != 1:
            raise BundleExportError(f"{root}/dimension must be 1 for scalar CVs")
        kind = _name(_tokens(_scalar(reader, root + "/type"), root + "/type"), root + "/type")
        atoms = _selection(reader, root, atom_count, virtual_names)
        items = {"CV_type": kind}
        if atoms:
            items["atom"] = " ".join(atoms)
        for key in reader.list_children(_PROTOCOL, root + "/parameter"):
            path = root + "/parameter/" + key
            _name(key, path)
            if key in _RESERVED:
                raise BundleExportError(f"{path} uses a reserved native CV field")
            items[key] = _tokens(reader.read(_PROTOCOL, path), path)
        for key in ("period", "sigma", "rotate", "function", "min_padding", "max_padding"):
            path = root + "/" + key
            if reader.contains(_PROTOCOL, path):
                values = reader.read(_PROTOCOL, path)
                expected = (1,) if key in {"period", "sigma"} else ()
                if values.shape != expected:
                    raise BundleExportError(f"{path} must have shape {expected}")
                items[key] = _tokens(values, path)
        reference = _reference(reader, root, kind, len(atoms))
        if reference is not None:
            items["coordinate"] = _tokens(reference, root + "/coordinate")
        sections.append((name, items))
    return sections
