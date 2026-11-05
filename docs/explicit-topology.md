# Explicit topology API

Both Xponge and XpongeCPP expose the following API through `import Xponge as x`:

```python
topology = x.get_explicit_topology(molecule)
topology["bonds"][0][2] = 120.0
x.set_explicit_topology(molecule, topology)
x.set_atom_lj(molecule, 0, epsilon=0.12, rmin=1.8)
```

This is a parameter interface, not a ByteFF model runner or a topology/atom
mapping generator. Call after structure assembly, with zero-based global atom
indices in residue/atom order. Assign charges and masses on the physical atoms.
Do not change connectivity after materialization without explicitly supplying
the corresponding new parameters. Disjoint materialized assembly is supported.
Xponge rejects a rebuild of an invalidated materialized molecule rather than
silently replacing external parameters with registry-derived parameters.

## Data contract

The dictionary has exactly five keys. Rows are copied on get and set:

| Key | Row | Units and potential |
| --- | --- | --- |
| bonds | i, j, k, r0 | angstrom, kcal/mol; E=k(r-r0)^2 |
| angles | i, j, k, constant, theta0 | radians, kcal/mol; E=constant(theta-theta0)^2 |
| dihedrals | i, j, k, l, periodicity, constant, phase | radians, kcal/mol; E=constant(1+cos(periodicity*phi-phase)) |
| nb14s | i, j, LJ scale, electrostatic scale | dimensionless; explicit scaled pairs |
| exclusions | one list of excluded indices per atom | symmetric, no self or duplicate indices |

Indices must be integers in range; parameters must be finite; harmonic
constants and pair scales must be nonnegative, bond lengths positive,
equilibrium angles in [0, pi], periodicities positive integers. Signed periodic
coefficients are permitted. Repeated additive energy terms are retained.
Periodic impropers use their ordered SPONGE/GROMACS atom tuple and the same
periodic energy formula. No graph-based improper permutation is performed.
Xponge stores supplied periodic terms as materialized `dihedral` entities;
chemical proper/improper classification is not part of this API.

The five tables are a **subset**, not a serialization of the whole force field.
Separate groups such as CMAP, harmonic improper, virtual sites, explicit A/B
pair potentials and soft bonds are not included; existing separate groups are
left unchanged. A five-table fingerprint must not be called a full Hamiltonian
fingerprint. Importers must reject external functions they cannot represent.

In Xponge, supplied exclusions must include those implied by the current
linked-atom graph and active exclusion policy. Explicit additions are stored
on the atoms and survive copying/projection. Load the force field first.
Neither backend infers chemical connectivity from supplied force rows.

`set_atom_lj` expects epsilon in kcal/mol and **Rmin/2** in angstrom (not sigma).
It uses immutable content-addressed LJ type names, so changing one molecule's
parameters does not overwrite another molecule's type parameters. Xponge uses
its global LJ registry with these unique names; XpongeCPP also retains local
overrides. Cross interactions follow the active force field's mixing rule.
The ByteFF adapter assumes the Amber Lorentz-Berthelot mixing convention.

## FEP compatibility

Xponge `fep.merge_dual_topology` now additionally accepts a single-residue
**Molecule** as the B replacement. Use that form for explicitly assigned
parameters: a `Residue`/`ResidueType` only carries the legacy template path.
The materialized physical B molecule is copied and projected onto the union.
XpongeCPP already accepts this form (its other arguments differ).

Copying, disjoint append, dual projection and lambda force-field interpolation
are covered by `tests/test_explicit_topology.py`. The Xponge FEP matcher now
matches repeated force entities one-to-one; multiple A entities must not share
one B entity. Complete reversal of a dihedral tuple preserves its angle;
arbitrary improper permutations are not equivalent and are retained separately.

## Validation and release boundary

Run tests in separate provider environments: both distributions own `Xponge`.
Check `Xponge.__file__` before testing. The shared API regression covers invalid
input transactional behavior, roundtrip/deepcopy, disjoint append, raw export,
LJ isolation, repeated periodic energies, full and partial MCS dual projection.
Existing FEP isolation and physical-topology regressions must also pass.

For XpongeCPP, rebuild the extension using the repository installation guide;
limit compiler/linker parallelism (for example one build job and serial LTO).
Changing Python files alone cannot add the native topology bridge. Local tests
used a rebuilt binding linked against clean-base cached objects, loaded from a
temporary provider directory, without replacing installed environments.

These are offline topology and polynomial-energy checks. No ByteFF inference,
GPU MD, real FEP, or full engine force validation was performed. The API is introduced in XpongeCPP 0.2.16 and Xponge 1.7.7.
Older installed 0.2.15/1.7.6 packages do not gain these APIs automatically.
The release CI runs installed-wheel regression tests on the new packages.
