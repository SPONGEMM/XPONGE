# Xponge 1.7.5

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
