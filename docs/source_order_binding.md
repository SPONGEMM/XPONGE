# CIF / H5MD source-order binding, version 1

The bundle saver binds source atom identities in its **final serialization
order** when `source_atom_ids` is supplied. XPONGE resolves IDs by Atom identity
after its preparation/reordering; XpongeCPP uses the native serialized
`molecule.atoms` order. Exports without source IDs retain the previous hashes.

## Hash contract

`digest(value)` is `"sha256:" + SHA256(canonical_json(value)).hexdigest()`.
Canonical JSON uses Python `json.dumps(sort_keys=True, separators=(",", ":"),
ensure_ascii=True, allow_nan=False)`, encoded as UTF-8. IDs are unique nonempty
strings. No Unicode normalization or whitespace trimming is performed.

```text
source_atom_order_hash = digest({
  "schema": "sponge-source-atom-order-v1",
  "source_atom_ids": [IDs in final simulation order]
})

atom_order_hash = digest({
  "schema": "sponge-bound-atom-order-v1",
  "base_atom_order_hash": original native atom_order_hash,
  "source_atom_order_hash": source_atom_order_hash
})
```

The original native order digest is retained as `base_atom_order_hash`. Thus
permuting distinct IDs changes the bound order digest even if the swapped atoms
have identical mass, charge and residue properties. The existing topology hash
continues to identify the native physical topology payload; the added binding
and source IDs are provenance metadata.

Before staged bundle publication, the saver writes:

- `/topology/source_atom_ids`: UTF-8 array in final simulation order.
- `/topology/source_order_binding`: UTF-8 JSON object defined below.
- `/topology/atom_order_hash`: new bound order digest.
- Restart `/run/atom_order_hash`: the same digest.

SPONGE already propagates this opaque topology order hash to
`/parameters/sponge/topology_compatibility/atom_order_hash` in output H5MD,
and restores global atom order before coordinate output. No per-frame IDs or
SPONGE format change is required.

## Mapping sidecar

Mokda checks that the serialized source IDs exactly equal the final mapping's
`external_id` sequence and that the binding agrees with current topology
metadata. It then adds an optional `topology_binding` object to the existing
`sponge-atom-order-mapping` version 1 JSON and save manifest:

```json
{
  "schema": "sponge-source-order-binding",
  "schema_version": 1,
  "hash_algorithm": "sha256",
  "atom_count": 2,
  "base_atom_order_hash": "native exporter digest",
  "source_atom_order_hash": "sha256:...",
  "atom_order_hash": "sha256:...",
  "topology_hash": "native topology digest"
}
```

This is bound during native topology export, **never by copying values from a
chosen trajectory**. Chemcore's ordered CIF and existing `mapping_hash` remain
unchanged: they prove the mapping's complete identity rows correspond to the
CIF atom order. The source-ID digest connects those rows to the native topology.

## Analysis validation

`load_cif_h5md_universe` checks:

1. CIF mapping order, IDs and digest, then equality with the sidecar.
2. Binding schema and recomputed source-ID / bound-order digests.
3. CIF/H5MD atom counts.
4. Both available H5MD `atom_order_hash` and `topology_hash` against the binding.

Conflicting available hashes or malformed bindings are rejected even with
`strict=False`. Missing sidecars, old sidecars without a binding, and older
H5MD files missing hashes remain readable with a warning and unverified status.
A custom particle stream cannot be proven by SPONGE's global hashes and is
reported as unverified. Raw exports currently have no source-order binding.

Results expose `universe.cif_metadata["trajectory_order_validation"]` with
`verified`, `method` (when verified), or `reason` (when unverified). Mokda copies
this into result `meta.trajectoryOrderValidation` and
`analysis_inputs.json.trajectory_order_validation`. This remains available when
the original topology H5 has been removed; the CIF and sidecar suffice.

This verifies export provenance and order consistency under the writer's
coordinate-order contract. It is not a signature, nor does it detect coordinates
edited after simulation while their provenance metadata is deliberately retained.
