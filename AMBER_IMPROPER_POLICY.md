# AMBER periodic improper policy: amber-deterministic-v1

This policy applies to periodic impropers loaded through the AMBER parmdat and
frcmod loaders, including GAFF/GAFF2. It does not redefine CHARMM harmonic
impropers or OPLS ordering, and is not a claim of bitwise LEaP compatibility.

## Matching and ordering

- The third parameter slot is the central atom. Search peripheral permutations
  only to match parameter types; never assume all these periodic functions are
  mathematically equivalent.
- Prefer the most specific matching definition (fewest wildcard slots). At that
  specificity, the latest definition of each periodicity overrides earlier ones;
  distinct periodicities are retained as separate entities.
- Preserve the original slots of each selected parameter, including after a
  differently ordered alias is reloaded. Among equivalent type matches, select
  the lexicographically least sequence of `(actual type, residue position,
  unique residue-local atom name)` in parameter-slot order. Duplicate local
  names within a matching group are rejected, not resolved with object hashes.
- Names and residue identity are part of this modeling convention. Atom-vector
  reordering and exporter renumbering do not change the selected function;
  arbitrary renaming or a different template identity is not guaranteed to do so.
- Loading a generic/custom definition over an AMBER type invalidates its AMBER
  metadata. Untagged force types retain the legacy generic matching path.

## FEP and existing topology

Physical A and B terms are built separately, then projected without rematching
on union atoms. Keep every source periodicity and phase. The interpolated
potential is `(1-lambda)*U_A + lambda*U_B`. Coalesce only identical ordered
functions, including safe whole reversal. Other mathematically equivalent
permutations can remain separate weighted terms; merging them is not required
for correctness. Distinct axes must never be merged just because their unordered
atom set and scalar parameters match.

This policy is for lazy topology construction, not rewriting imported or already
materialized topology. In particular, physical endpoint projection and native
topology overrides must preserve their ordered terms.

## Compatibility and validation boundaries

LEaP is an independent reference, not an oracle for ambiguous atom permutations.
Its library-to-unit parameter copy can lose original `sOrder` metadata; we do
not reproduce occurrence/cache-dependent ordering. See OpenMM issue #220:
https://github.com/openmm/openmm/issues/220

Use fixed nonplanar geometries, energy and numerical-gradient tests, repeated
groups, atom reorderings, and intermediate-lambda tests. Cross-backend matching
does not establish full-engine correctness or MD/free-energy convergence.
Existing frozen inputs and trajectories are not migrated by this source change.
