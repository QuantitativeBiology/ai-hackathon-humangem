# M5 A12(b) panel pre-check: non-author review (2026-10-28)

Reviewer script: results/m5/review/rv_m5_score.py (output results/m5/review/rv_m5_result.json). It does not import
m5_score_train or w4_stage_a_*; it uses nextgem.labels (train(), load_split, gene_universe), nextgem.metrics
(tie-aware per_line tAP, bootstrap_mean) and the Stage-A panel file. Training tier only, L_prob, universe U, 30 panel lines.

Label probe (head -c 50 on the raw gene-effect CSV): BLOCKED by the label-guard hook, as required.

## Recomputed numbers (lineage-equal mean, stratified paired bootstrap, B=2000, seed SEEDS[0])
| Contrast | Reviewer | Own bootstrap CI (cross-check) | Author |
|---|---|---|---|
| F2 R1+B-loo vs H1+B-loo | +0.00306 [+0.00166, +0.00451], 8/10 | [+0.00171, +0.00454] | +0.0031 [0.0017, 0.0045] 8/10 |
| F4 R1+B-loo vs H1+B-loo | +0.00440 [+0.00355, +0.00520], 10/10 | [+0.00361, +0.00519] | +0.0044 [0.0035, 0.0052] 10/10 |
| sanity: H1+B-loo vs generic H1 | +0.03319 [+0.03096, +0.03543], 10/10 | [+0.03106, +0.03540] | +0.0332 |

Digits match train_decision.json to 6+ places (0.003064, 0.004397, 0.033188). Fold positive = lineage mean > 0 with >= 2 lines.

## Checks
- Score: e = round(1 - ratio_r, 6); NaN ratio -> 1 (e = 0); U genes without a Human1 gene -> 0. Human1 Entrez map = sort by
  gene, drop duplicate Entrez (w4_analysis.entrez_gene_map / W1 convention). Reimplemented independently; same results.
- Failed builds: builds.tsv has 0 failed of 61 and all 60 line KO tables exist. The author's rule (failed -> all-tied
  score 0 in both arms of a contrast, §11.4) is therefore never triggered; my script checks for none rather than
  re-implementing the rule.
- Provenance: results/m5/MANIFEST.tsv rows for m5_build.py and m5_score_train.py both show commit 4729547db7a5..., dirty = False
  (pinned worktree /data/tmpA/nextgem_runs/4729547).
- The author's nested-LOLO step is trivial (one fixed config per family; outer delta asserted equal to the plain delta),
  so a plain paired delta is equivalent.
- F4 fails the build gate (8/30 pass, GR-units artefact), so it is descriptive regardless; this does not affect the decision.

## Decision
Holds: neither F2 (+0.0031) nor F4 (+0.0044) reaches +0.005, so the context lever is not selected (Stage B/C skipped).
Caveat: F4's CI upper bound (+0.0052) straddles the bar; its point estimate misses it by 0.0006. Under the pre-registered
point-estimate rule the verdict is unambiguous.

## Training aggregates viewed
results/m5/train/train_decision.json (full). My script printed per-arm lineage-equal tAP (generic H1 0.2982, H1+B-loo 0.3314,
F2 0.3345, F4 0.3358) and the contrasts above. I did not open train_per_line.tsv, train_folds.tsv, any label file or data/locked.
(README.md and the code in scripts/m5_score_train.py were read; the README's contrast table repeats the headline numbers.)
