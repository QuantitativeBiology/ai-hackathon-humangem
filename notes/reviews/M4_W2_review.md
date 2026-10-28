# M4 W2 review (non-author), 2026-10-28
Probe (head -c 50 on the raw CRISPR gene-effect CSV): BLOCKED by the label hook (not retried). Labels only via
labels.train() and cv.lolo (fold-scoped CE_fit); no raw label file, no data/locked/. Own script:
results/w2/review/rv_w2_score.py (output rv_w2_score.json). It imports nextgem.labels/metrics/paths/w1/cv only,
not w2_score*. Attribution stratum: 466 training lines, 12 lineages, L_prob, all U genes, lineage-equal, B=2000, SEEDS[0].

## Verdicts
| # | Number | Author | Recomputed | Verdict |
|---|---|---|---|---|
| 1 | M2(f=1) vs M0 (H1) | +0.0159 [0.0146, 0.0172], 11/12 | +0.0159 [0.0146, 0.0172], 11/12 | CONFIRMED |
| 1 | M2(f=1) minus CE_fit | +0.0182 | +0.0182 [0.0167, 0.0198], 11/12 | CONFIRMED |
| 1 | M-open vs M0 | +0.0195 [0.0180, 0.0212] | +0.0195 [0.0180, 0.0212], 11/12 | CONFIRMED |
| 1 | M-open minus CE_fit | -0.0168 | -0.0168 [-0.0185, -0.0150], 0/12 | CONFIRMED |
| 2 | documented - shuffled | +0.0005 [-0.0000, +0.0010] | author's shuffle file: +0.00052 [-0.00003, +0.00102]; my own shuffle (own RNG, SEEDS[0:5]): +0.00053 [-0.00003, +0.00106] | CONFIRMED |
Folds positive: 11 of 12 for both (Cervix has 4 attribution lines, so it counts as not positive; max = 11).
Shuffled-minus-M0 is about +0.0153 per seed, so a random medium keeps ~96% of the gain. The CI lower bound is
marginally below 0 in both versions: W2-SHUF (negative control not beaten) stands.

## Checks
3. TG fold rule: DESIGN §9 TG (3) says >= 9 of 12 LOLO folds, a fold with < 5 estimable lines is "not positive",
   denominator stays 12. w2_score.py has TG_MIN_FOLDS = 9, FOLD_MIN_LINES = 5, delta >= 0.005, CI lower > 0. MATCHES.
4. Validation parquet w2_val_scores.parquet: columns ModelID, entrez, H1, W2:M2(f=1)[H1], H2, W2:M2(f=1)[H2]. Scores
   only, no label columns, no NaN, values in [0,1]. 198 lines = exactly the validation-tier lines with an expression
   matrix (200 validation lines, 198 with expression, none rare); 2866 genes/line = U; no duplicates. H1 and H2 equal
   w1.generic_table e_H1/e_H2 exactly (max abs diff 0). 3.2% of entries differ between W2 and H1. PASS.
5. Eligible set fixed before scoring: NOT VERIFIABLE from history; partial support only.
   - w2_score.py and w2_score_train.py were untracked until commit 80c992d (after scoring); no earlier commit.
   - PROVENANCE shows two scoring runs: 00:59 (HEAD 97fa948) and 01:39 (HEAD 8d2eb26). w2_score.py mtime is 01:24 and
     w2_score_train.py mtime is 01:37, both after the first run, and PROVENANCE records no script hashes, so the code
     of run 1 is not preserved and I cannot show ELIGIBLE was unchanged between runs. Run-1 outputs were overwritten.
   - Supporting: the set (M1, M2(f), gate-free, M3 excluded) follows the DESIGN W2 ladder/"gate-free" text and A9i
     (commit 63a51df, 10-27, before scoring).
   - Impact is small: all five eligible steps pass TG and M2(f=1) wins by the one-SE rule. The README sentence
     "fixed in w2_score.py before scoring" should be softened to "consistent with DESIGN; not hash-verified".

## Other problems
- The shuffle file is permuted over all 835 attribution lines (valid within lineage x serum-free, checked on all 135
  groups), so within the 466 training lines the medium multiset is not preserved. Only medium text is mixed, not
  labels; my own within-training shuffle gives the same answer. Minor.
- M-open passes TG (+0.0195, larger than M2(f=1)) but is -0.0168 minus CE_fit: the gain is mostly removal of
  common-essential false positives, which a random medium also gets. Supports the W2-MOPEN flag.
- README sigma*sqrt(2 ln m) uses m=5 (0.0047) while the spec uses m=1 (A11g convention); a labelling difference only.

## Training aggregates viewed
results/w2/train/step_scores.tsv (attribution H1 rows M2(f=1), M-open only), results/w2/train/shuffled_control.tsv
(pooled M2(f=1) H1 row), plus my own aggregates in results/w2/review/rv_w2_score.json. No per-line or per-gene
label tables opened.
