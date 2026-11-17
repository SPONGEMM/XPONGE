## 1.7.8 · FEP endpoint topology and native 1–4 controls

`merge_force_field(..., specific_lambda=...)` accepts four optional [0, 1]
weights for native 1–4 LJ: `nb14_junction_lj_a`, `nb14_junction_lj_b`,
`nb14_internal_lj_a`, and `nb14_internal_lj_b`. Junction terms connect a
mapped atom and an atom unique to the specified side. Internal terms have
two atoms unique to that side. The controls leave 1–4 Coulomb, custom
explicit pairs, and bonded force interpolation on their existing paths.
Omitting the controls preserves the existing behavior.

This release also preserves explicit endpoint exclusions and branch identities,
and exposes common-core pair-category transitions to the unified FEP exporter.
The Python provider preserves empty force-table shapes during bundle conversion.
The two providers expose matching native 1–4 LJ control and pair-transition APIs.
Mokda owns endpoint bonded-term selection and lambda scheduling; this release
does not claim formal FEP convergence or production sampling validation.

# Xponge 1.7.7

Add public `get_explicit_topology`, `set_explicit_topology`, and `set_atom_lj`
interfaces for external materialized parameters, shared with XpongeCPP 0.2.16.
FEP dual construction accepts a single-residue Molecule as B so externally
assigned parameters survive projection. Invalidated materialized topology is
rejected instead of silently reassigning parameters from templates.

Fix repeated independent force-entity matching: each B entity is consumed at
most once. This complements earlier fixes for multiple Fourier terms inside
one dihedral entity and the separate ordered-improper path; it is not a
reversion of those fixes. Release tests cover repeated terms and endpoint energy.

The `mokda-xponge` compatibility distribution pins `Xponge==1.7.7`.
See `docs/explicit-topology.md`; no model inference, real FEP or full-engine
force validation is claimed.

# Xponge 1.7.6

Includes the deterministic AMBER improper ordering and retained dummy bonded
scaffolds described below. Installed-wheel regression checks now distinguish
ordinary dual-branch scaffolds from explicitly one-sided low-level soft bonds,
including force parameters, serialization and molecule-copy preservation.

The CDK8 parameter-provenance and cross-backend ordered-improper limitations
below remain; full-engine/MD equivalence is not claimed.

# Xponge 1.7.5 (unpublished candidate)

## deterministic AMBER impropers (2026-10-23)

AMBER-loaded periodic improper parameters preserve original slots, use stable
type/residue/local-name ties, and retain distinct periodicities through overrides.
FEP mixes ordered Fourier terms without atom-tuple-only pairing or union rematching.
Generic parameter updates invalidate AMBER-only metadata; other force-field
conventions remain on their existing path. See AMBER_IMPROPER_POLICY.md for the
explicit policy and intentional differences from ambiguous LEaP cache ordering.
No frozen inputs/trajectories are migrated; local checks are not full-engine acceptance.

## retained dummy bonded scaffolds (2026-10-23)

Dual construction now retains source bonded entities touching non-MCS atoms
in both union endpoints, including boundary angles and ordered torsions.
Dummy LJ and charge remain zero; common-only parameters still switch to the
correct physical endpoint. No force matching is rerun on the union graph.
This corrects the physical-only projection below, not ordinary modeling.
Validation is incomplete: historical CDK2 ligand-subset tables match at seven
step2 states, but CDK8 parameter provenance and one cross-backend ordered
improper discrepancy remain. No full-engine/MD equivalence is claimed.

# Xponge 1.7.4

Fix FEP endpoint parameter ownership, ordered improper terms, soft-bond
partitioning, and explicit NB14 coefficients. Physical endpoint topology and
materialized force terms survive molecule copies and disjoint assembly.

PDB and mmCIF file-declared inter-residue bonds take precedence over competing
implicit template links. mmCIF now supports `ignore_conect=True`, ignoring file
bond declarations and edits while retaining template connectivity and separately
provided residue links. Noncovalent `struct_conn`, Mokda semantic, and bond-edit
records no longer become covalent force-field connections.

Add the read-only `amber_parameter_exists` query for registered parameter
coverage, including supported Amber wildcard matching.

The `Xponge` and `mokda-xponge` compatibility distributions share version 1.7.4.
This release preserves the existing Mokda template and explicit-residue-link
handoff; no complete-connectivity interface is required.
