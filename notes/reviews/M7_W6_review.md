# M7 W6 review (non-author; no labels read; nothing committed)
Session probe (head -c 50 on the raw CRISPR gene-effect file): BLOCKED by the label guard hook. No label files or data/locked used.
Scripts: own code in scratchpad (chk1/chk2/chk3.py; w6.py not imported). BLAS 1 thread, 1 process, ~0.6 CPU-h.

## 1. Decisions recomputed from per-unit tables (2 arms x 3 contrasts x 3 analyses = 18 rows)
- Medians, n, Wilcoxon (one-sided, zeros dropped), Holm across 3 analyses, above-permuted flag and verdict all match
  decisions.tsv. Only deviation: P1 metabolomics H_flux p 0.8571 vs 0.8574 (tables rounded to 6 digits). Immaterial.
- Verdict: CONFIRMED. All 18 "not supported". Smallest p_Holm for flux-adds is 0.38 (G* drugs); no row is positive,
  significant and above permuted. Drug set = 106 non-off-model units.

## 2. Independent refit (SVD-PCA k=20 fitted without held-out lineage, standardise, lineage one-hot, RidgeCV LOO, LOLO)
| unit | set | mine | pipeline |
| GDSC2 methotrexate | RNA | 0.206444 | 0.206444 |
| GDSC2 methotrexate | RNA+FLUX | 0.460201 | 0.460201 |
| asparagine | RNA | -0.047630 | -0.047630 |
| asparagine | RNA+FLUX | 0.318841 | 0.318841 |
- Verdict: CONFIRMED (agreement ~1e-6, well inside 1e-3). Pipeline ridge/PCA/LOLO is correct.

## 3. Structural claim (generic flux depends only on assay medium)
- Distinct states per assay (both arms): PRISM 1; GDSC 2 (RPMI, DMEM:F12); culture medium (CTRP, metabolomics,
  subtypes) 20 states from 22 medium keys; the 20 states are 20 distinct flux vectors. No line-level information.
- Of 106 primary drug units, 77 have a single state; all 77 have Delta exactly 0; all other 29 are nonzero.
- Verdict: CONFIRMED.
- README correction: all 77 constant units are PRISM; CTRP (19 units) and GDSC2 (17) are never constant. The README
  wording "PRISM, plus GDSC and CTRP units where constant" is inaccurate.
- Caveat: the big gains (methotrexate +0.25, asparagine +0.37) are medium labels (folate/thymidine, asparagine
  content of the medium), so flux acts as a medium encoding. Worth stating beside the positive-control pairs.

## 4. Controls vs DESIGN
- Permuted flux: line-to-state map shuffled within lineage, SEEDS[0:5], "above permuted" = real median above all 5.
  Matches DESIGN. CONFIRMED.
- Synthetic positive: y = tanh(2z) + noise (SNR 1), z a random projection of the flux block over states (flux-only,
  nonlinear), 20 units; recovered (median +0.42 / +0.45, p = 9.5e-7 = 2^-20, all 20 positive). CONFIRMED.
  Limitation: signal depends only on state (<= 20 states), so it shows detectability of medium-state information,
  not of line-specific flux.
- Synthetic negative: 200-column random linear map of model-gene expression added to RNA; no gain (drugs -0.007,
  p 0.97; metabolomics +0.0003, p 0.14; subtypes 0). CONFIRMED as implemented; reported in secondary.tsv only,
  with no pass/fail threshold in the decision rule.

## Overall
The "not supported" decisions (both hypotheses, both arms, all 3 analyses) stand. Minor: fix the README sentence in
sec. 3; the NaN-target cosmetic issue is already disclosed by the author.
