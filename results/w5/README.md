# W5 — GPR curation (training tier; no validation look)

Author: W5 engineer (Claude Opus), 2026-10-27. Design: DESIGN.md v3 §9 W5 (+ §5–6, §10.3, §11, §13; amendments A1–A10;
a priori lists per A2 = `notes/apriori_lists_final.md`). Nothing here is committed (the PI commits).

**Label hygiene.** Session start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label hook.
Labels were read only through `nextgem.labels.train()` and `labels.fold(k)` (via `nextgem.cv.lolo`) inside
`scripts/w5_score_train.py` and `scripts/w5_val_dryrun.py`; `labels.validation()`/`labels.test()` were never called and
`data/locked/` was never touched. De Kegel 2021 / Ito 2021 (`data/paralogs/*.xlsx`) were never opened (the validator
refuses them). Label-copy scanner over `results/w5`, the W5 code and `notes/USER_DOWNLOADS.md`: **0 hits**.

## Headline (training tier, 12 LOLO lineages, 686 expression lines, L_prob, all U genes; Δ̄_lin vs the generic parent)

| Arm | H1 (= G\*) Δ̄_lin [95% CI] | folds > 0 | H2 Δ̄_lin [95% CI] | folds > 0 | TG |
|---|---|---|---|---|---|
| TYMS only (MAR04644 irreversible; H1-only, mandated by §3.1) | **+0.0007** [+0.0004, +0.0010] | 9/12 | — | — | fail (< 0.005) |
| **W5a** (literal DESIGN rules; H1 incl. TYMS) | **−0.0286** [−0.0289, −0.0282] | 0/12 | **−0.0293** [−0.0295, −0.0291] | 0/12 | fail |
| W5a GPR part only (H1 released) | −0.0297 [−0.0300, −0.0295] | 0/12 | — | — | fail |
| W5a_pid30 (sensitivity: paralog = ≥ 30 % identity) | −0.0021 [−0.0024, −0.0018] | 1/12 | −0.0031 [−0.0032, −0.0030] | 0/12 | fail |
| **W5b** (W5a + in-fold co-essentiality) | identical to W5a (0 fold edits) | | identical to W5a | | fail |
| Negative controls (5 seeds, matched random edits) | −0.0057 … −0.0145 | 0/12 | −0.0031 … −0.0151 | 0/12 | — |

Without CE_fit (LOLO): TYMS +0.0021 [+0.0018, +0.0024] (12/12); W5a −0.0008 (H1), −0.0032 (H2); W5a_pid30 +0.0013
[+0.0010, +0.0017] (10/12, H1), −0.0011 (H2). Fixed-cutoff MCC (ratio_r < 0.01 × L_prob): H1 0.3218 → W5a 0.2815
(Δ −0.0403), TYMS +0.0015, W5a_pid30 −0.0014; H2 0.3184 → W5a 0.2757 (Δ −0.0427). L_ge05 / L_ge10 / U∩ agree in sign.
All numbers: `train_deltas.tsv`, `train_grid.tsv`, `train_binary.tsv`, `train_summary.json`.

**Shortlist (≤ 3): empty.** No curated arm passes TG on its generic parent (DESIGN W5 decision rule). Δ_curation is
reported regardless (above). The H1 TYMS item is not a candidate: §3.1 puts it in every optimised H1 pipeline; its
training effect is small and positive. W5 therefore proposes **no validation look** (a contingent config is provided,
§8 below).

**Why curation hurts.** The literal AND → OR rule ("paralogs and not co-members of any complex") merges distant Ensembl
`other_paralog` subunits of RNA polymerases I/II/III (e.g. POLR1A/POLR2A/POLR3A, 21–25 % identity) inside the
RNA-synthesis GPR (MAR07161/2), making 15 POLR genes non-lethal; the 13 in U include 12 with training positive fraction
0.85–1.00 (POLR1D 0.07; POLR2A/POLR2J have no CRISPR column). Restricting paralogs to ≥ 30 % identity removes
≈ 0.026 of the −0.029 (H1). W5a is worse than every random-edit control. With the
≥ 30 %-identity paralog definition the remaining edits roughly cancel: RRM2B (FP fixed, as §1.3 predicts) is paid for by
RRM2 becoming non-lethal (positive in 99.9 % of training lines), COX7B (80 %) is lost to its COX7B2 buffer, and the
pan-cancer-absence removals make IDI1, SC5D and FLAD1 lethal although they are positive in only 2–6 % of lines.
No OR → AND edit qualified in either parent (only 0–1 OR nodes in the targets have co-complex evidence for every
cross-branch pair), and no W5b candidate reached Chronos r ≥ 0.4 in any fold (H1 candidates on all training lines: |r| < 0.1).

## 1. What was done

| Step | Script (library) | Outputs (results/w5/) |
|---|---|---|
| Label-free targets: every single edit per GPR group, local re-evaluation, rank, top 200 | `scripts/w5_enumerate.py` (`src/nextgem/w5.py`) | `targets_<parent>.tsv`, `enum_<parent>.json`, `tmp/enum_<parent>.parquet`, `tmp/base_<parent>.parquet` |
| W5a evidence edits + evidence-ID validator (+ W5a_pid30 sensitivity; W5b candidates) | `scripts/w5_curate.py` (`src/nextgem/w5_evidence.py`) | `edits_W5a[_pid30]_<parent>.tsv`, `curated_W5a[_pid30]_<parent>.tsv`, `or_nodes_*`, `and_paralog_pairs_*`, `blocked_*`, `evidence_validation.tsv`, `apriori_targets_check.tsv`, `w5b_candidates.json`, `curation_summary.json` |
| KO sweeps of all arms + controls (M0, GLPK) | `scripts/w5_apply.py` | `ko_<parent>.parquet`, `gene_changes.tsv`, `neg_control_edits.tsv`, `pos_control_v1v2.tsv`, `apply_summary.json` |
| 2 × 2 factorial setup (constraint arm stubbed) | `scripts/w5_factorial.py` | `factorial_2x2_setup.json` |
| Training scoring, LOLO, W5b refit in every fold (reuses `nextgem.w1`, `nextgem.metrics`) | `scripts/w5_score_train.py` | `train_grid.tsv`, `train_deltas.tsv`, `train_per_lineage_delta.tsv`, `train_gene_changes.tsv`, `train_binary.tsv`, `w5b_fold_edits.tsv`, `train_summary.json` |
| Proposed validation-look config + label-free dry run (A6) | `scripts/w5_val_dryrun.py` | `proposed_val_look_W5.json`, `val_dryrun.json` |
| Ordered clean re-run + determinism check + runtime | `scripts/w5_run_all.py` | `rerun_sha_check.tsv`, `runtime.tsv`, `MANIFEST.tsv` |
| Tests (stable names) | `tests/test_w5.py` (16 tests) | — |

**KO engine.** Same semantics as `nextgem.knockouts` (closing every reaction whose GPR is False with the gene off),
GLPK only, M0, biomass MAR13082, ratio_r rounded to 6 decimals; exact zero-flux shortcut; identical closed sets share one
LP; per-LP GLPK limit 60 s with one from-scratch retry (presolve, 300 s), then an error (0 retries were needed in the
final runs). Checks: the full sweep equals `results/baseline/<parent>_ko.parquet` for every gene (max |Δ| = 0, both
parents); a full sweep of the curated H1 model equals the local re-evaluation exactly (`factorial_2x2_setup.json`);
toy known-answer test (OR → AND makes an isozyme lethal; AND → OR rescues a subunit; equals cobra's own knock-out).

## 2. Label-free targets (`targets_<parent>.tsv`)

| | Human2 v2.1.1 | Human1 v1.19.0 (+ TYMS base) |
|---|---|---|
| reactions with a GPR / distinct canonical GPR groups | 7,746 / 1,885 | 8,043 / 1,726 |
| single edits enumerated (OR→AND / AND→OR / remove branch) | 9,617 (972 / 3,663 / 4,982) | 8,307 (834 / 3,527 / 3,946) |
| (edit, gene) pairs whose closed set changes / distinct LPs solved | 11,801 / 694 | 10,448 / 523 |
| groups where some edit flips a binary call / changes any ratio | 62 / 109 | 51 / 92 |
| **targets** (top 200 ∩ any ratio change) — groups / reactions | **109 / 928** | **92 / 1,064** |

Top targets (both parents): RNA synthesis MAR07161/2 (37 genes), MAR13086 (32), the large acyl-CoA/CYP/UGT OR blocks,
MAR07160, the RNR block (rank 14/16), SEC11A/SEC11C (15/17). Binary call = ratio_r < 0.01 on M0.

## 3. W5a evidence-driven edits (`edits_W5a_<parent>.tsv`)

Rules exactly as the DESIGN W5 table, applied bottom-up to the target GPR groups (an edit applies to every reaction of
the group). Evidence: Ensembl 116 paralogs (8,971 pairs among model genes; 856 within_species, 8,115 other_paralog),
CORUM 5.1 + Complex Portal (930 model genes in ≥ 1 complex), ProCan 20250211 (948 lines; 1,653 model proteins) then
Nusinow (375 lines; 1,966), DepMap 26Q3 expression (1,292 screened lines with a default profile).

| | H2 | H1 |
|---|---|---|
| groups / reactions edited | 29 / 341 | 20 / 405 |
| remove OR branch (pan-cancer absent) | 49 | 26 |
| AND → OR (paralogs, no shared complex) | 11 | 12 |
| OR → AND | 0 | 0 |
| protected OR nodes (§3 list) skipped | 12 of 93 | 10 of 79 |
| candidate blocked by protection | PRPS1L1 (MAR04052, PRPS1/PRPS2 OR) | same |
| evidence IDs validated (both parents, both variants) | 538 / 538 exist locally | |

Removed (pan-cancer absent) genes: PAPOLB, CYP3A43, CYP2C8, CYP1A2, CYP2C19, CYP4Z1, CYP2F1, CYP2A13, CYP2A7, CYP2A6,
CYP26C1 (H2), CYP7A1, UGT2B10, UGT2B28, UGT3A1, UGT2A1, RBP3, ALAS2, IDI2, LENEP, TECRL, DAO, ACSBG2, ACER1, FTMT (H2),
ACOT6, AGXT2 (H2). AND → OR merges: POLR (5 cliques), RRM2/RRM2B (2 nodes), SEC11A/SEC11C, COX7B/COX7B2, HADHA/HADH
(other_paralog, 13 % id), HADHB/ACAA2 (H1).

A priori items (`apriori_targets_check.tsv`, §1.3/§3 "curation targets"): recovered — RNR AND → `RRM1 and (RRM2 or
RRM2B)`, SEC11A/SEC11C, TECR or TECRL (TECRL removed). Not recovered — ATP5MC1/2/3 (co-members of the ATP-synthase
complex, so the AND → OR rule cannot fire; this is the §3 "modified" note), MAR06435/6 (targeted, no evidence),
ALDH block and MAR01074 (not targets), HMGCS1/2 (no shared GPR).

The UniProt "unrelated function" removal path is **not implemented**: no UniProt file is on disk; listed as item 12
in `notes/USER_DOWNLOADS.md` (not downloaded).

## 4. W5b (`w5b_candidates.json`, `w5b_fold_edits.tsv`)

Candidate OR nodes after W5a: H1 11 (0 with complex anchoring), H2 17 (1 anchored). Rule (primary): all cross-branch
pairs non-paralog **and** co-complex **and** (co-abundance ≥ 0.5 or Chronos r ≥ 0.4 on the fold's fitting lines, ≥ 30
lines). Sensitivity "W5b_loose": non-paralog + co-essentiality only. Recomputed inside every LOLO fold from
`labels.fold(k)`: **0 edits in all 12 folds for both variants and both parents**; on all training lines the 13 H1
candidate pairs have |r| < 0.1 (`val_dryrun.json`). W5b is therefore identical to W5a.

## 5. Controls

- *Negative* (`neg_control_edits.tsv`): per seed, the same number and type of single edits as W5a (49/11 H2, 26/12 H1),
  one per GPR group, among target groups with no W5a edit and non-protected OR nodes (fallback to non-target groups
  when the target pool is exhausted, logged as `pool`). Δ̄_lin −0.003 to −0.015. W5a (−0.029) is below all ten.
- *Positive* (`pos_control_v1v2.tsv`; possibly label-informed, PI decision 16): of 1,064 H1 target reactions, 169
  changed GPR between v1.19.0 and v2.1.1, almost all by gene-set change (153; 1 logic-only). W5a(H1) reproduces the H2
  rule exactly for 1/169 (0.6 %) and edits 33 % of them in some way.

## 6. Gene-level changes (`gene_changes.tsv` label-free; `train_gene_changes.tsv` with training prevalence)

Binary flips vs generic: H2 W5a 22 (pid30 7), H1 W5a 24 (pid30 9; TYMS alone 2). Per gene (training L_prob positive
fraction in brackets): POLR1A/1B/1C/1H/2B/2C/2G/2I/3A/3B/3H/3K lethal → viable [0.85–1.00] (wrong), POLR1D [0.07]
(right), POLR2A/POLR2J (no CRISPR column); RRM2 lethal → viable [0.999] (wrong), RRM2B [0.001] (right); SEC11A [0.006], SEC11C [0.00] (right); COX7B
partial → viable [0.80] (wrong), COX7B2 [0.001] (right); IDI1 [0.055], SC5D [0.036], FLAD1 [0.022] viable → lethal
(wrong for most lines); H1 TYMS: UMPS [0.58] and DHODH [0.42] viable → lethal.

## 7. CAD rescued by orotate but not uridine (W2 review item)

Not a GPR or reaction-direction error. `(S)-dihydroorotate` (MAM00180c) is a **biomass constituent** in both parents:
MAR10064 ("small metabolite pool for biomass reaction") consumes 0.0043 per unit. Its only producers are CAD's
dihydroorotase (MAR04608, reversible) and DHODH (MAR04575, reversible). With uridine, UMP is made by UCK and growth
rises (H2 WT 74.8 → 92.1), but no dihydroorotate can be made without CAD, so the CAD KO stays lethal; closing any one of
CAD's steps MAR04034/MAR04050/MAR04608 is lethal under uridine (MAR09802, glutaminase "CAD or GLS2", is not).
Orotate rescues because DHODH runs in reverse (orotate → dihydroorotate). DHODH/UMPS KOs are rescued by uridine (W2)
because CAD still supplies the pool. GPR curation cannot fix this; the candidate fix is a biomass-composition edit
(drop dihydroorotate from MAR10064, or a W3 objective variant), outside W5's remit. Same in H1 (H1 has no
dihydroorotate exchange; H2 has MAR20202, closed in all media).

## 8. Validation look (proposed, not run) and dry run (A6)

`proposed_val_look_W5.json` (look id `W5_M3_gpr_curation`, status **NOT TO BE RUN under the frozen rule**): methods
H1_TYMS, H1_W5a, H1_W5a_pid30, H1_W5b vs G\* = H1; comparators H1, H2, B2, B2m, B8±/B8′ (H1); cells L_prob/L_ge05/
L_ge10 × {all, noCE} × {U, U∩} × {none, expressed pairs}; scalars only; §8 gate (paired Δ̄_lin vs G\* CI lower bound
> 0) then the paired one-SE rule; W5b refitted on all training lines (§5.3). It needs a new `--look` branch in the
frozen `evaluate_val.py`, i.e. a harness change with a full independent code review (A6 item 2). Dry run
(`val_dryrun.json`): 198 validation expression lines × 2,866 U genes for every method, pairing asserts pass, no NaN,
score SHA-256 recorded, **no validation label read**.

## 9. 2 × 2 factorial (`factorial_2x2_setup.json`)

Arms {generic, curated GPR} × {generic, best constraint pipeline}; curated GPRs are substituted into the reactions
present in an *already-extracted* model (`nextgem.w5.curated_for_extracted`), no re-extraction. The generic-constraint row
is computed (the arms above); the constraint arm is a stub (`best_constraint_pipeline` raises NotImplementedError) until
W2/W3/W4 have shortlists. Breakdown by edit class: RMBRANCH, AND2OR, OR2AND, TYMS.

## 10. Candidate FROZEN.tsv rows (PI action; SHA-256 of the final re-run)

`candidate_frozen_rows.tsv` (FROZEN.tsv columns; `frozen_on`/`amendment` left for the PI): the edit lists
`edits_W5a_<parent>.tsv` (H2 `d9e32079…`, H1 `86d489d0…`), the curated rules `curated_W5a_<parent>.tsv` (H2 `d617fa26…`,
H1 `48deec88…`), the targets, the W5a_pid30 sensitivity lists, `w5b_candidates.json`, the arm KO tables, the two W5
library files, and a value row for the H1 TYMS item (MAR04644 bounds (0, 1000)). Freezing matters only if a W5 arm
enters the §10.1 factorial or a look; none is shortlisted.

## 11. Deviations and ambiguities for the PI

1. **Edit unit = canonical GPR group, not single reaction** (*ours*). An edit applies to every reaction with the same
   canonical GPR. The literal per-reaction reading scores 0 flips for enzymes with several individually essential
   reactions (the 14 RNR reactions, the a priori §1.3 item), so they could never be targets.
2. **"Top 200" = top 200 ∩ any KO effect.** Only 109 (H2) / 92 (H1) groups have any single edit that changes any gene's
   ratio; below that every group ties at (0, 0, 0), so filling 200 would be arbitrary.
3. **Ranking key** (*ours*): number of distinct genes flipped by any single edit of the group (union), then the
   largest single-edit flip count, then genes with any ratio change. The union favours large AND groups (rank 1 =
   RNA synthesis, 37 genes, max single edit 4).
4. **AND → OR granularity** (*ours*): enumerated as pairwise merges of AND siblings; applied as cliques of qualifying
   pairs (non-clique components logged, none occurred).
5. **"Ensembl paralogs" is literal (any Compara pair, incl. `other_paralog` at 13–25 % identity).** This is what
   produces the harmful POLR and HADHA/HADH merges. W5a_pid30 (≥ 30 % identity both ways, W1's convention) was declared
   in code before scoring as a sensitivity row and is **not** shortlist-eligible without a PI decision (it fails TG
   anyway).
6. **"Soluble enzyme" proxy** (no UniProt file): Human-GEM `genes.tsv` location without "Cell membrane" and no transport
   reaction in either parent. "Never detected" = no non-missing value in ProCan averaged matrix and Nusinow. "Lines" for
   TPM < 1 = the 1,292 screened lines with a default expression profile (all tiers; label-free).
7. **OR → AND evidence** is required for every cross-branch pair of distinguishing genes (genes not in all branches).
   Co-abundance: ProCan if ≥ 30 shared lines, else Nusinow (fallback on availability, not on r).
8. **W5b reading** (*ours*): co-essentiality replaces the co-abundance leg (complex + non-paralog still required);
   "W5b_loose" (co-essentiality + non-paralog) as sensitivity. Moot here: both produce 0 edits.
9. **Protected list**: §3 encoded verbatim (`PROTECTED_S3`, every fragment tested against the file). Listed reactions:
   all their OR nodes are protected; gene groups: any OR with ≥ 2 branches carrying genes of one group. "ATP6V0A1-A4"
   resolves A3 through the alias of TCIRG1. Pan-cancer-absence removal inside protected ORs (§3 "tissue-restricted
   isoforms" note allows it) was **not** applied, per the brief ("never edit the isozyme-exception list"); one case
   (PRPS1L1) is logged in `blocked_*.tsv`.
10. **H1 enumeration on the TYMS base** (*ours*), since every optimised H1 carries it; the engine check is on the
    released H1.
11. **Provenance with uncommitted code**: `allow_dirty=True` (the brief forbids commits; W2 precedent); every MANIFEST row
    records the SHA-256 of every W5 code file. Rows before the `w5_run_all.py` run are superseded development runs.
12. **Negative-control pool**: target groups without W5a edits first, then non-target groups; random AND → OR edits in
    the pool are not paralog-restricted (matched on type and number only).
13. Tests were run with `NEXTGEM_NO_TEST_RECORD=1` (no rows added to `results/test_runs.tsv`); the PI's full run records them.

## 12. Runtime and provenance

Final ordered re-run (`scripts/w5_run_all.py`, `runtime.tsv`):

| Step | Wall (s) | CPU (s) | Workers |
|---|---|---|---|
| enumerate H2 | 118 | 234 | 6 |
| enumerate H1 | 149 | 230 | 6 |
| curate (W5a + pid30, validator) | 52 | 51 | 1 |
| apply (all arms + controls) | 220 | 226 | 6 |
| factorial check (H1 full sweep) | 69 | 69 | 1 |
| score_train (LOLO, W5b refits) | 119 | 203 | 6 |
| val dry run | 25 | 25 | 1 |
| **Total** | **752 (12.5 min)** | **1,038 (0.29 CPU-h)** | ≤ 6 |

All W5 work including development runs: 0.69 CPU-h recorded in `MANIFEST.tsv` (+ < 0.1 CPU-h unrecorded exploration),
against the §13 W5 budget of 60 CPU-h. Peak RSS 2.3 GB. `MANIFEST.tsv` (argv, commit, dirty flag, input/output
SHA-256, versions, seeds, solver, workers, wall/CPU, RSS, W5 code SHAs). All LPs on GLPK; **no Gurobi session was
used**; ≤ 6 worker processes. `rerun_sha_check.tsv`: all 48 outputs of the final re-run are bitwise identical to the previous run except the 7 summary JSONs, which embed wall/CPU times.

## 13. Label-exposure statement

The W5 engineer (Claude Opus) saw: (a) training-tier aggregates written by W1 (its README headline numbers); (b) the
training-tier outputs of this workstream: LOLO Δ̄_lin/CI/fold counts, MCC, per-lineage Δ on training lineages, the
training positive fraction of the ≈ 100 genes whose score changes, and all-training co-essentiality r for the 13 H1 W5b
candidate gene pairs. No gene × line label value was printed or written; no validation- or test-tier label, label summary or
sealed output was read; the raw label files were never opened (hook probe blocked). Validation/test ModelIDs were read
from the label-free split file only for the dry-run score matrices. Exposure class: training-tier aggregates only.

## Independent review and PI notes (2026-10-27)
- **Review:** notes/reviews/M3_W5_review.md (GO). The reviewer's own code gives the same curated closure sets, identical KO ratios and the same Δ̄_lin (H1 W5a −0.02855, H2 W5a −0.02928, H1 TYMS-only +0.00069, 9/12 folds). There is no leakage: W5b co-essentiality is computed inside `cv.lolo` folds, and the scanner found 0 hits.
- **Edit review (DESIGN-required; 20% sample plus all POLR merges, 28 edits):** 6 accept, 14 reject, 8 questionable.
  - All 12 POLR AND→OR merges are rejected on biology. Pol I/II/III subunits sit in distinct, non-redundant complexes, and MAR07161 is a lumped pseudo-reaction that needs all three. These merges account for ≈ 0.026 of the loss: without them W5a gives −0.0019 (H1) and −0.0030 (H2).
- **Documented deviations found by the review:**
  - The "soluble enzyme" proxy (no Cell membrane localisation, no transport reaction) lets about 17 ER-membrane CYP/UGT genes through. This does not meet DESIGN's "soluble" qualifier. Under a strict reading only the CYP7A1 edit changes a score (SC5D becomes lethal), which is far below the TG margin.
  - Proteomics detection is tested by exact UniProt accession only. CYP1A2 is detected in Nusinow under the isoform accession P05177-2, and UGT2A1 under a superseded accession (Q9Y4X1-4, 1 peptide), so both removals are wrong. Neither changes a score.
  - The validator does not re-check the ≥ 95% TPM < 1 fraction. The reviewer re-checked it, and it holds for all edits.
- **Wording of the result:** "W5a as specified (without the UniProt path) fails TG on both parents". This is not a claim that GPR curation cannot help. The ≥ 30% identity sensitivity also fails (−0.0021 / −0.0031), and the H1 TYMS fix helps but stays below the 0.005 SESOI.
- **Provenance:** the MANIFEST rows are dirty because the code was uncommitted at run time. The reviewer confirmed that the code hashes recorded in them equal the committed files.
