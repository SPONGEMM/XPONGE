"""Materialize native bias objects for the legacy SPONGE runtime."""
from __future__ import annotations

import numpy as np

from .cv_exporters import _enabled, _name, _scalar, _tokens
from .errors import BundleExportError
from .legacy_materializer import LegacyPayload
from .topology_exporters import _format_float, _text_vector

PROTOCOL = "protocol.spgp.h5"
RESTART = "restart.spgr.h5"


def _array(reader, path, shape=None, integer=False, bundle=PROTOCOL):
    values = np.asarray(reader.read(bundle, path))
    if values.dtype.kind not in ("iu" if integer else "iuf"):
        raise BundleExportError(f"{path} has invalid numeric type")
    if shape is not None and values.shape != shape:
        raise BundleExportError(f"{path} must have shape {shape}")
    if not np.isfinite(values).all():
        raise BundleExportError(f"{path} must contain finite values")
    return values


def _refs(reader, root):
    values = reader.read(PROTOCOL, root + "/cv_refs")
    if values.ndim != 1 or not values.size:
        raise BundleExportError(f"{root}/cv_refs must be a nonempty vector")
    refs = _text_vector(values)
    for ref in refs:
        _name(ref, root + "/cv_refs")
        cv = "/cv/" + ref
        if not reader.contains(PROTOCOL, cv) or not _enabled(reader, cv):
            raise BundleExportError(f"{root} references missing or disabled CV {ref!r}")
        if _scalar(reader, cv + "/dimension", 1) != 1:
            raise BundleExportError(f"{root} requires scalar CVs")
    return refs


def _text_rows(values):
    return "\n".join(" ".join(_format_float(x) for x in row) for row in values) + "\n"


def _sections_payload(key, sections):
    lines = []
    for name, values in sections:
        lines.extend([name, "{"])
        lines.extend(f"    {key} = {value}" for key, value in values.items())
        lines.append("}")
    return [LegacyPayload(key, "\n".join(lines) + "\n")] if lines else []


def config_sections(reader, root):
    root += "/config"
    if not reader.contains(PROTOCOL, root):
        return []
    names = _text_vector(reader.read(PROTOCOL, root + "/section/name"))
    offsets = _array(reader, root + "/section/key_offset", (len(names) + 1,), integer=True)
    keys = _text_vector(reader.read(PROTOCOL, root + "/key"))
    values = _text_vector(reader.read(PROTOCOL, root + "/value"))
    if offsets[0] != 0 or offsets[-1] != len(keys) or np.any(np.diff(offsets) < 0) or len(keys) != len(values):
        raise BundleExportError(f"{root} has invalid section offsets")
    return [(name, dict(zip(keys[offsets[i]:offsets[i+1]], values[offsets[i]:offsets[i+1]])))
            for i, name in enumerate(names)]


def merge_sections(legacy, native):
    sections = {}
    for name, values in legacy + native:
        target = sections.setdefault(name, {})
        for key, value in values.items():
            if key in target and target[key].split() != value.split():
                try:
                    equal = np.array_equal(np.asarray(target[key].split(), float), np.asarray(value.split(), float))
                except ValueError:
                    equal = False
                if not equal:
                    raise BundleExportError(f"section {name!r} field {key!r} conflicts with native protocol")
            target[key] = value
    return list(sections.items())


def prepare_native_protocol(reader, context):
    """Return per-key payload overrides; collect mdin defaults in context."""
    overrides = {}
    positional = []
    cv_restraints = []
    for name in reader.list_children(PROTOCOL, "/restraint", groups_only=True):
        root = "/restraint/" + name
        if not reader.contains(PROTOCOL, root + "/type"):
            continue  # compatibility datasets are handled by existing exporters
        kind = _tokens(_scalar(reader, root + "/type"), root)
        if kind not in {"harmonic_positional", "cv_harmonic"}:
            raise BundleExportError(f"{root} has unsupported restraint type {kind!r}")
        if kind == "harmonic_positional":
            for key in ("restrain_atom_id", "restrain_weight_in_file", "restrain_coordinate_in_file", "restrain_amber_rst7"):
                overrides.setdefault(key, [])
            if _enabled(reader, root):
                positional.append((name, root))
        else:
            overrides.setdefault("restrain_cv_in_file", [])
            if _enabled(reader, root):
                cv_restraints.append(root)
    if len(positional) > 1:
        raise BundleExportError("legacy runtime supports at most one enabled positional restraint")
    for name, root in positional:
        atom_count = int(reader.read_scalar("topology.spgt.h5", "/topology/atom_count"))
        indices = _array(reader, root + "/atom_indices", integer=True)
        if indices.ndim != 1 or not indices.size or np.any(indices < 0) or np.any(indices >= atom_count) or len(set(indices.tolist())) != indices.size:
            raise BundleExportError(f"{root}/atom_indices is invalid")
        overrides["restrain_atom_id"] = [LegacyPayload("restrain_atom_id", "\n".join(map(str, indices)) + "\n")]
        reference_path = f"/parameters/restart/references/restraint/{name}/coordinate"
        reference = _array(reader, reference_path, (atom_count, 3), bundle=RESTART)
        overrides["restrain_coordinate_in_file"] = [LegacyPayload("restrain_coordinate_in_file", f"{atom_count}\n" + _text_rows(reference))]
        weight_path = root + "/weight"
        if reader.contains(PROTOCOL, weight_path):
            weights = _array(reader, weight_path, (indices.size, 3))
            overrides["restrain_weight_in_file"] = [LegacyPayload("restrain_weight_in_file", _text_rows(weights))]
        elif not reader.contains(PROTOCOL, root + "/single_weight_default") and "restrain_single_weight" not in context.commands:
            raise BundleExportError(f"{root} requires weight or single_weight_default")
        for field, key in (("single_weight_default", "restrain_single_weight"),
                           ("refcoord_scaling_default", "restrain_refcoord_scaling"),
                           ("calc_virial_default", "restrain_calc_virial")):
            if reader.contains(PROTOCOL, root + "/" + field):
                context.mdin_defaults[key] = _tokens(_scalar(reader, root + "/" + field), root + "/" + field)
    if "restrain_cv_in_file" in overrides:
        arrays = {key: [] for key in ("CV", "weight", "reference", "period", "start_step", "max_step", "reduce_step", "stop_step")}
        for root in cv_restraints:
            refs = _refs(reader, root)
            arrays["CV"].extend(refs)
            for key in arrays:
                if key == "CV":
                    continue
                integer = key.endswith("step")
                path = root + ("/schedule/" if integer else "/") + key
                if reader.contains(PROTOCOL, path):
                    values = _array(reader, path, (len(refs),), integer=integer)
                elif key in {"weight", "reference"}:
                    raise BundleExportError(f"{path} is required")
                else:
                    values = np.zeros(len(refs), dtype=int if integer else float)
                arrays[key].extend(values.tolist())
        native = [("restrain", {key: _tokens(value, key) for key, value in arrays.items()})] if cv_restraints else []
        overrides["restrain_cv_in_file"] = _sections_payload("restrain_cv_in_file", merge_sections(config_sections(reader, "/restraint/cv"), native))
    if reader.contains(PROTOCOL, "/steer/cv_refs"):
        native = []
        if _enabled(reader, "/steer"):
            refs = _refs(reader, "/steer")
            weights = _array(reader, "/steer/weight", (len(refs),))
            native = [("steer", {"CV": " ".join(refs), "weight": _tokens(weights, "/steer/weight")})]
        overrides["steer_cv_in_file"] = _sections_payload("steer_cv_in_file", merge_sections(config_sections(reader, "/steer"), native))
    if reader.contains(PROTOCOL, "/sits/method"):
        overrides["SITS_in_file"] = []
        if not _enabled(reader, "/sits"):
            overrides.update({"SITS_atom_in_file": [], "SITS_nk_in_file": []})
            return overrides
        method = "/sits/method"
        mode = _tokens(_scalar(reader, method + "/mode"), method + "/mode")
        if mode not in {"observation", "iteration", "production", "empirical", "amd", "gamd"}:
            raise BundleExportError(f"invalid SITS mode {mode!r}")
        native = {"mode": mode}
        for field in ("k_numbers", "temperature_low", "temperature_high", "pe_a", "pe_b", "fb_interval", "fb_bias", "record_interval", "update_interval", "nk_rest", "nk_fix", "cross_enhance_factor"):
            if reader.contains(PROTOCOL, method + "/" + field):
                key = {"temperature_low": "T_low", "temperature_high": "T_high"}.get(field, field)
                native[key] = _tokens(_scalar(reader, method + "/" + field), field)
        if reader.contains(PROTOCOL, method + "/temperature_ladder"):
            ladder = _array(reader, method + "/temperature_ladder")
            if ladder.ndim != 1 or not ladder.size or np.any(ladder <= 0):
                raise BundleExportError("SITS temperature_ladder must be a positive vector")
            native["T"] = "/".join(_format_float(x) for x in ladder)
        if reader.contains(PROTOCOL, "/sits/atom_numbers_policy"):
            if reader.contains(PROTOCOL, "/sits/atom_indices"):
                raise BundleExportError("SITS atom selection has conflicting ownership")
            native["atom_numbers"] = _tokens(_scalar(reader, "/sits/atom_numbers_policy"), "/sits/atom_numbers_policy")
        for section, values in merge_sections(config_sections(reader, "/sits"), [("SITS", native)]):
            if section != "SITS":
                raise BundleExportError("SITS config must contain the SITS section")
            context.mdin_defaults.update({"SITS_" + key: value for key, value in values.items()})
        for field in ("log_norm", "log_nk"):
            if reader.contains(RESTART, "/parameters/restart/bias/sits/SITS/" + field):
                raise BundleExportError(f"legacy SITS input cannot restore {field}; use the H5 checkpoint to preserve restart state")
    return overrides
