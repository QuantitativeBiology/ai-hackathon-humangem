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
