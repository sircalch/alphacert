# Changelog

## 1.1.0 (2026-09-27)

Validated against real AlphaFold DB v6 models: haemoglobin α (P69905) and p53 (P04637). The
references are independent of AlphaCert:
- mean pLDDT and pLDDT band fractions published by the AlphaFold DB API (exact match);
- φ/ψ computed by gemmi (maximum difference 3·10⁻¹⁴°);
- the known domain architecture of p53.

### Fixed
- **Dihedral sign.** φ, ψ and ω had the opposite sign to the IUPAC convention, so every α-helix
  (φ ≈ −60°) was placed in the left-handed region. Ramachandran statistics were wrong for every
  protein; for example, p53 had 67 % "outliers".
- **Clash score.** Bonded atoms were excluded by atom index distance (fewer than 4 apart), so every
  peptide bond C(i)–N(i+1) and many 1–3 pairs counted as clashes. Haemoglobin, with a pLDDT of 98,
  scored 770. Bonds are now taken from the covalent graph and pairs up to three bonds apart are
  excluded. All atoms are checked with a KD-tree; before, only the first 1500 were checked while the
  count was divided by all atoms. N/O pairs count as clashes only below 2.5 Å, because they may be
  hydrogen bonds or salt bridges. The score is documented as a heavy-atom proxy, not the MolProbity
  clashscore.
- **Made-up Ramachandran values.** Without backbone atoms, the module returned invented values
  (95 % favoured). It now returns NaN and evaluates only the clash score.
- **Chain breaks.** φ/ψ is now computed only between residues joined by a peptide bond (same chain,
  C–N < 2 Å).
- **mmCIF.** mmCIF files, the AlphaFold DB default, were read as empty structures without any error.
  They are now parsed from the `_atom_site` loop and give the same result as the PDB file. Only the
  first model and the first alternate location are read, water is ignored, and residues keep the
  file order, so insertion codes are handled.
- **PAE domains.** The domain count came from an ad hoc windowed average and was capped at 5 (p53
  gave 5). It now uses T. Croll's greedy-modularity clustering (ChimeraX and AlphaFold DB method;
  resolution 0.5, cutoff 5 Å), discarding clusters with fewer than 30 residues or a mean pLDDT
  below 70. Haemoglobin gives 1 domain; p53 gives the DNA-binding (94–294) and tetramerisation
  (325–357) domains.
- **Complex interfaces.**
  - The interface passed if the inter-chain PAE was ≤ 10 Å *or* the ipTM was ≥ 0.6. When ipTM is
    given it is now the primary criterion, following AlphaFold's guidance (> 0.8 PASS, 0.6–0.8
    WARNING, < 0.6 FAIL).
  - The inter-chain PAE is symmetrised and covers all chain pairs, not only chain A against the rest.
- **`--pocket` residues.** The residue numbers were used as 0-based positions, which shifted the
  pocket pLDDT by one residue, or by more when numbering does not start at 1. They are now mapped
  by residue number (optionally with a chain, `A:248`), and residues that do not exist are reported.
- The generated methods text reports the installed version instead of a fixed "v1.0.0".

### Added
- `pae_domains()` and `count_heavy_atom_clashes()`.
- A `plddt` argument in `evaluate_pae_matrix` and a `domains` field in `PAEAnalysisResult`.
- `networkx` dependency.
- Tests on the real AlphaFold DB files (`tests/test_afdb_real.py`).
