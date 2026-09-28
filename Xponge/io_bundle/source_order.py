"""Versioned binding of serialized source identities to native atom order."""

from __future__ import annotations

import hashlib
import json

import h5py
import numpy as np


def _digest(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return "sha256:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()


def make_source_order_binding(source_atom_ids, base_atom_order_hash, topology_hash):
    """Bind the actual serialized identity sequence to its native order hash."""
    identities = list(source_atom_ids)
    if not identities or any(not isinstance(value, str) or not value for value in identities):
        raise ValueError("source atom IDs must be nonempty strings")
    if len(set(identities)) != len(identities):
        raise ValueError("source atom IDs must be unique")
    if not all(isinstance(value, str) and value for value in (base_atom_order_hash, topology_hash)):
        raise ValueError("source order binding requires native topology and atom-order hashes")
    source_hash = _digest({"schema": "sponge-source-atom-order-v1", "source_atom_ids": identities})
    order_hash = _digest({
        "schema": "sponge-bound-atom-order-v1",
        "base_atom_order_hash": base_atom_order_hash,
        "source_atom_order_hash": source_hash,
    })
    return {
        "schema": "sponge-source-order-binding",
        "schema_version": 1,
        "hash_algorithm": "sha256",
        "atom_count": len(identities),
        "base_atom_order_hash": base_atom_order_hash,
        "source_atom_order_hash": source_hash,
        "atom_order_hash": order_hash,
        "topology_hash": topology_hash,
    }


def validate_source_order_binding(binding, source_atom_ids):
    """Reject malformed bindings and identity sequences different from export."""
    if not isinstance(binding, dict):
        raise ValueError("source order binding must be an object")
    expected = make_source_order_binding(
        source_atom_ids, binding.get("base_atom_order_hash"), binding.get("topology_hash")
    )
    if binding != expected:
        raise ValueError("source order binding does not match the atom identity sequence or schema")
    return expected


def _text(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _write_text(handle, path, value):
    if path in handle:
        del handle[path]
    handle.create_dataset(path, data=value, dtype=h5py.string_dtype("utf-8"))


def bind_bundle_source_order(paths, source_atom_ids):
    """Bind a staged bundle before publication; keep restart lineage consistent.

    Only the exporter may call this with IDs in its actual serialization order.
    No coordinates or force-field datasets are changed.
    """
    identities = list(source_atom_ids)
    with h5py.File(paths.topology, "r+") as top, h5py.File(paths.restart, "r+") as restart:
        if "/topology/source_order_binding" in top:
            raise ValueError("bundle source atom order is already bound")
        count = int(np.asarray(top["/topology/atom_count"][()]).reshape(-1)[0])
        if len(identities) != count:
            raise ValueError("source atom IDs must cover every serialized atom")
        binding = make_source_order_binding(
            identities,
            _text(top["/topology/atom_order_hash"][()]),
            _text(top["/topology/topology_hash"][()]),
        )
        top.create_dataset("/topology/source_atom_ids", data=identities, dtype=h5py.string_dtype("utf-8"))
        _write_text(top, "/topology/source_order_binding", json.dumps(binding, sort_keys=True))
        _write_text(top, "/topology/atom_order_hash", binding["atom_order_hash"])
        _write_text(restart, "/run/atom_order_hash", binding["atom_order_hash"])
    return binding


def read_source_order_binding(topology, source_atom_ids):
    """Read an export binding, checking its stored IDs and current H5 metadata."""
    with h5py.File(topology, "r") as top:
        if "/topology/source_order_binding" not in top:
            return None
        binding = validate_source_order_binding(
            json.loads(_text(top["/topology/source_order_binding"][()])), source_atom_ids
        )
        if top["/topology/source_atom_ids"].asstr()[...].tolist() != list(source_atom_ids):
            raise ValueError("topology source atom IDs do not match the final mapping")
        for name in ("atom_order_hash", "topology_hash"):
            if _text(top[f"/topology/{name}"][()]) != binding[name]:
                raise ValueError(f"topology {name} does not match its source order binding")
        count = int(np.asarray(top["/topology/atom_count"][()]).reshape(-1)[0])
        if count != binding["atom_count"]:
            raise ValueError("topology atom count does not match its source order binding")
    return binding


def validate_trajectory_source_order(trajectory, binding, *, particle_stream=None):
    """Check available trajectory hashes; legacy missing hashes stay unverified.

    SPONGE compatibility hashes describe the global ``all`` particle stream.
    They cannot establish the identity order of an arbitrary custom stream.
    """
    if binding is None:
        return {"verified": False, "reason": "missing_source_order_binding"}
    missing = []
    with h5py.File(trajectory, "r") as handle:
        streams = list(handle.get("particles", {}))
        selected = particle_stream or (streams[0] if len(streams) == 1 else None)
        if selected != "all":
            return {"verified": False, "reason": "unbound_particle_stream"}
        for name in ("atom_order_hash", "topology_hash"):
            path = f"/parameters/sponge/topology_compatibility/{name}"
            if path not in handle or not _text(handle[path][()]):
                missing.append(name)
            elif _text(handle[path][()]) != binding[name]:
                raise ValueError(f"H5MD {name} does not match the CIF mapping source order binding")
    if missing:
        return {"verified": False, "reason": "missing_trajectory_hash", "missing": missing}
    return {"verified": True, "method": "source_order_binding_v1", **binding}
