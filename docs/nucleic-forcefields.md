# Amber nucleic-acid force fields

## Loaders

- amber.ol3: RNA; parm10 already contains OL3 (do not apply it twice).
- amber.bsc1 / amber.ol15: existing DNA alternatives.
- amber.ol24: DNA; parm10 followed by the complete official OL24 frcmod and OFF-derived templates.
OL3 automatically registers RNA HOP3-bearing 5′ monophosphate templates; OL24 registers the corresponding DNA templates. No extra module is required. Explicit selection of *5MP never changes the default 5′-OH mapping.

Internal nucleotides connect P to the preceding O3′; *5 connects only downstream, *3 only upstream, *N is an isolated nucleoside. End-residue partial-charge sums are fractional: the complete chain, rather than each residue separately, has integer charge.

## Terminal phosphate contract

The official terminal_monophosphate library is converted to names A5MP, C5MP, G5MP, U5MP, DA5MP, DC5MP, DG5MP and DT5MP. These have OP3–HOP3 and phosphate formal charge −1. RNA first-residue partial sum is −1.3081; DNA is −1.3079. A two-residue chain ending in *3 totals −2.

These Amber terminal templates retain their source atom types. In the pinned upstream revision, DNA variants use CJ/C7/C2 and therefore require OL24 parameters; they are not compatible with the bsc1/OL15 loaders alone. RNA variants use the OL3 types. No automatic state conversion is performed. A dianionic phosphate request raises AMBER_5PRIME_PHOSPHATE_STATE_UNSUPPORTED; removing HOP3 alone would produce wrong partial charges. An independently validated parameter source is required for that state.

## Provenance

The amber_sources directory contains the original OFF libraries, frcmod, leaprc and upstream license. provenance.json records upstream commit and SHA256 checksums. Generated MOL2 files preserve atom names, partial charges, atom types and connectivity; manifests record connection anchors and source identity. The converter is scripts/convert_amber_off_to_mol2.py. OFF connectivity flags represent connectivity, not chemical bond-order inference.

Upstream: https://github.com/Amber-MD/AmberClassic/tree/2fdfc7a64c63514b18340b1833341bbdafbef105/dat/leap

## References

- Zgarbová et al. Refinement of the Cornell et al. Nucleic Acids Force Field Based on Reference Quantum Chemical Calculations of Glycosidic Torsion Profiles. JCTC 2011. https://doi.org/10.1021/ct200162x
- Zgarbová et al. Refinement of the Sugar-Phosphate Backbone Torsion Beta for AMBER Force Fields Improves the Description of Z- and B-DNA. JCTC 2015. https://doi.org/10.1021/acs.jctc.5b00716
- Refinement of the Sugar Puckering Torsion Potential in the AMBER DNA Force Field. https://doi.org/10.1021/acs.jctc.4c01100

Common single-atom metal ion aliases are bound after the selected water module loads its ion parameters. Native atom matching requires a single-atom residue, a single charged template atom, and the same element; source atom name and serial are retained.
