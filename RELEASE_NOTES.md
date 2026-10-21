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
