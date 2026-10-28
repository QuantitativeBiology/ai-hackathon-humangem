# Label-exposure register

Every viewing of DepMap label *values* (or label-derived summaries) by any agent or the PI. Required by DESIGN.md §11.5. Agents listed here must not author or approve a priori lists, candidate sets or decision thresholds.

> **Note (2026-10-26):** the "must not author or approve" rule below is superseded in part by the PI decision of 2026-10-26 at the end of this file (disqualifying exposure = gene × line values or post-freeze validation/test summaries).

## Entries (as of DESIGN v2, 2026-10-26)

- the baseline review read the CRISPR header only;
- the literature agent saw the all-line DepMap common-essential count;
- the **round-1 statistics review**:
  - evaluated the retired seed 20261026 under the v0 rules;
  - named the old test tier (Pancreas, Soft Tissue, PNS, Prostate, Eye) and validation tier (CNS/Brain, Ovary/Fallopian Tube, Bladder/Urinary Tract, Cervix);
  - computed label summaries on the old tier: positives per line, eligible-gene counts and a prevalence/SD remark;
  - also computed all-line prevalences and the common-essential share of positives;
- the **round-1 data reviewer** computed ASNS prevalence by medium and per-library positive rates on all lines;
- the **design author**, while inspecting a file header with `head`, printed about 10 L-prob values for one line (ACH-002730);
- the **all-line label copy** `ge_U.parquet` (CRISPRGeneEffect restricted to U, all 1,361 lines) and other derived scratch files (`testtier.py`, `splits.pkl`). These were written by a round-1 reviewer in the shared scratchpad and **deleted by the PI on 2026-10-26**;
- the **PI's header-only read** of a label file while testing the hook.

**Consequences.** The retired seed's tiers are known, which is why the seed is replaced by a user-supplied one (§4.1). All other exposures are aggregate (no gene × line information on the new test tier) and method-symmetric. Agents in the register **must not author or approve** a priori lists, candidate sets or decision thresholds; the gene panels and the drug list are produced by unexposed agents at M1 (open question F-4 on the lists already in this document). Whoever writes an a priori list attests in the commit message that they are not in the register.

## Later entries

- 2026-10-26, PI: two hook-liveness probes after restart were blocked (no values seen). Earlier pre-activation probes printed the gene-name header row only.
- 2026-10-26, PI: read reviewer summaries containing all-line aggregate label statistics (common-essential share of positives, ASNS essentiality by base medium, per-library positive rates; values redacted here, see the reviews). No gene × line values. The PI therefore does not author a priori lists (F-4); an unexposed agent re-approves them.
- 2026-10-26, round-3 verifier: read this register and the round-2 statistics review, which contained all-line aggregate label statistics (no gene × line values, nothing on the new test tier). Probe of raw labels blocked.
- 2026-10-26, PI decision (resolves DESIGN F-4): **aggregate all-line label statistics do not disqualify** an agent from authoring or approving a priori lists. Disqualifying exposure = gene × line label values, or any label summary on the new (post-freeze) validation/test tiers. The design author (~10 values for one line) therefore does not author lists; the round-3 verifier (aggregates only) may, and authored `notes/apriori_lists_approved.md`.
- 2026-10-26, PI: deleted remaining label-reading scratch scripts (testtier.py, stats1.py, prev.py, splits_sim.py, apnull.py, genes2.py) left by round-1 agents.
- 2026-10-26, **user decision on F-4** (asked by the PI after the M1 process review flagged that the PI relaxed the rule post hoc): re-check. A fully unexposed agent re-derives the a priori lists blind (notes/apriori_lists_blind.md, without reading DESIGN, LITERATURE, notes/ or the existing lists); differences with notes/apriori_lists_approved.md are resolved on biology and logged. Future changes to leakage/exposure rules require user ratification.
- 2026-10-26, W1 analyst (M2): computed and viewed training-tier label aggregates (per-line and per-lineage tAP, per-gene training prevalence, binary confusion summaries). No gene × line values; no validation or test labels; raw probe blocked.
- 2026-10-26, PI: read the W1 summary (training-tier aggregates only).
- 2026-10-26, M2 W1 reviewer (notes/reviews/M2_W1_review.md): recomputed and viewed training-tier aggregates via labels.train()/fold() only (lineage-equal tAP for H1/H2/B2, per-lineage sign of H2−H1, P57×L-ge06 MCC, CE_fit size, positive-pair count, fraction outside the addressable set); read W1 training outputs. No gene × line values; no validation/test labels or summaries; raw probe blocked. Label-free dry run of evaluate_val.py on validation line IDs printed no metric values.
- 2026-10-26, PI: read the M2 W1 review (training-tier aggregates only).
- 2026-10-27, PI: ran the single W1 validation look (scalars only: G* = Human1, Δ̄_lin(H2−H1) and CI, validation tAP and MCC summaries). The validation-look outputs are committed (09b69a2); any agent reading results/w1/val_look_* sees validation-tier aggregates.
- 2026-10-27, M3 W4 reviewer: label-free review; incidentally opened results/w1/val_look_W1_M2_baseline_gstar.json (committed validation aggregates; not used in any check).
- 2026-10-27, W2 engineer: read FROZEN.tsv (label.primary κ = 0.934, a label-agreement aggregate from M1). Otherwise label-free.
- 2026-10-27, M2 process reviewer: process documents only, including the committed training and validation aggregates in RESEARCH_TREE.md.
- 2026-10-27, W4 engineer (Stage-A analysis): computed and viewed training-tier aggregates (nested-LOLO tAP_line, Δ̄ vs H1/H2, MCC per family and threshold on the 30-line panel). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-27, PI: read the W4 Stage-A analysis summary (training aggregates only).
- 2026-10-27, W4 Stage-A reviewer: recomputed training-tier aggregates (per-line tAP compared by max-abs difference only, nested Δ̄ per family, fold choices, variance partition) via labels.train()/fold(). No gene × line values written; no validation or test; raw probe blocked.
- 2026-10-27, W5 engineer: training-tier aggregates only (LOLO Δ̄, CIs, fold counts, MCC, per-lineage Δ for the W5 arms; training positive fraction of ~100 changed genes; all-training co-essentiality r for 13 candidate pairs). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-27, W5 reviewer: training-tier aggregates only (LOLO Δ̄, CIs, fold counts for 8 contrasts under L_prob and L_ge05). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-27, W3 engineer: training-tier aggregates only (per-mode and per-line tAP, LOLO Δ̄ with CIs and fold means, context p/q values, per-task call precision, MCC/F1/precision/recall, A8 gene-pFBA deltas). Saw the symmetric matched-call-rate result before choosing the literal DESIGN control (deviation D-1). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-27, PI: read the W3 report (the same training aggregates as listed by the engineer, summary level).
- 2026-10-28, W3 reviewer: training-tier aggregates only (mode_scores.tsv rows for T0, both shortlisted B-loo and L-T10 on H1/H2; own recomputed tAP/Δ/CI/folds for 6 cells; post-hoc random-removal control). No per-line or gene × line values; no validation or test; raw probe blocked.
- 2026-10-28, PI attribution of unnamed tripwire entries flagged by the W3 reviewer: 2026-10-27 19:30:19 pid 1173952 = pytest toy fixtures (deletions of copies outside the data roots); 19:30–19:37 pids 1173994, 1174464, 1175510 and 20:22:14 pid 1187213 = W5 engineer/reviewer scoring runs via labels.train() (already registered above). All training-tier; no validation or test opens.
- 2026-10-28, harness reviewer (evaluate_val shortlist kind): synthetic dry-run scalars, the G* identity, and two W3 training-tier aggregates printed inside the winners_curse_note of results/w3/train/shortlist.json. No validation values (dry-run training positive rates hidden and the output deleted); raw probe blocked.
- 2026-10-28, PI: ran validation look W3_M4_shortlist (scalar-only outputs by design: lineage-weighted tAP/Δ̄ with CIs per cell, decision file). No per-line or gene × line validation values.
- 2026-10-28, M5 engineer: one accidental view of a single training-tier per-line aggregate row in results/w3/train/per_line.tsv (one line, mode T1: n, n_pos, tAP, e_random). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-28, M5 engineer (panel pre-check run): training-tier panel aggregates only (Δ̄_lin, CIs, fold means for 5 contrasts; per-arm lineage-equal tAP, mΔAP_gene, MCC on L_prob/L_ge06; minus-CE and expressed-pairs deltas). Did not open train_per_line.tsv. No validation or test.
- 2026-10-28, PI: read the M5 pre-check report (same aggregates, summary level).
- 2026-10-28, M5 reviewer: training-tier panel aggregates only (train_decision.json, own recomputed Δ̄/CI/folds and per-arm lineage-equal tAP). No per-line tables; no validation or test; raw probe blocked.
- 2026-10-28, W2 M4 engineer: training-tier aggregates only (step_scores for 28 steps × 2 parents and strata; shuffled-medium control pooled and per-lineage for 5 lineages; attribution.json; shortlist; ASNS by medium; binary metrics; sensitivities). Synthetic dry run only; no validation or test; raw probe blocked.
- 2026-10-28, W2 reviewer: training-tier aggregates only (attribution H1 rows for M2(f=1) and M-open; pooled shuffled-control row; own recomputed Δ̄/CI/folds). No per-line or per-gene tables; no validation or test; raw probe blocked.
- 2026-10-28, M6 engineer: training-tier aggregates only (per-cell tAP/Δ̄/CI/folds incl. minus-CE and expressed pairs, pairwise contrasts, interaction, §6.8 binary tables); accidentally n/n_pos/prevalence for 3 training lines in a header check of results/w1/train_per_line.tsv. Synthetic dry run only; no validation or test.
- 2026-10-28, PI: read the M6 training report (same aggregates, summary level).
- 2026-10-28, PI: ran validation look M6_factorial (scalar-only outputs: lineage-weighted tAP/Δ̄ with CIs per cell, decision file).
- 2026-10-28, T1-prep engineer: training-tier aggregates only (skew of P1 − G* per-line ΔtAP; lineage-equal mean ΔtAP for C1/C1b on training, matching M6). Synthetic dry run only; no validation or test values.
- 2026-10-28, W7 engineer: training-tier aggregates only (tAP_lin/CI per idea and seed, paired Δ vs current best/P1/B2/G*, folds, nested β, gene-grouped CV tAP). No gene × line values; no validation or test; raw probe blocked.
- 2026-10-28, W7 reviewer: training-tier aggregates only (autoresearch/results.tsv, screen_r2b/genecv, own LOLO aggregates, per-line training tAP in memory for comparison). Synthetic dry run only; no validation or test.
- 2026-10-28, PI: ran validation look W7_shortlist_1 (scalar-only outputs).
- 2026-10-28, PI: ran validation look W7_P2_select (scalar-only outputs).
- 2026-10-28, REPORT drafter: committed training and validation aggregates only (as in README/results tables); no label files; raw probe blocked.
- 2026-10-28, PI: ran the rehearsal (look 6/40, REHEARSAL_P1_M6): scalar outputs on validation revealed (results/m6/rehearsal/); per-unit outputs sealed by the harness (hashes only in the ledger).
