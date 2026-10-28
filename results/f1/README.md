# F1: GPR curation with the stricter paralog rule, on top of P1

**Post hoc / exploratory (A15).** This follow-up comes after the T2 unseal (DESIGN §17 A15; RESEARCH_TREE D17). It never replaces T1/T2 or the headline, and it gets no test evaluation. The F1 agent did not open any T1/T2 output, `results/m6/test/` or `data/locked/`.

Author: F1 engineer agent (Claude Opus), 2026-10-28. Status: **complete. F1 fails TG on training (curated P1 − P1 = −0.0006, 0/12 folds), so no validation look is prepared.** The code was committed at bc8089e before any label read. The clean re-runs and the scoring ran from the clean main tree (HEAD f57bffe for the re-runs, 6d8b109 for the scoring).

## Label hygiene

- The session-start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label hook.
- No label has been read so far. The only scoring run used `--synthetic-labels` (Bernoulli(0.12); no label file opened; output in `tmp/score_dry/`, git-ignored).
- Exposure carried in: the W5 README and review (training aggregates, including W5's per-gene training positive fractions for RRM2/RRM2B/SEC11A/SEC11C/COX7B). D17 fixed the rule before this agent started. The interpretation choices below were made on biology and annotation grounds, and are listed here so that the PI can ratify or change them before the first label read.

## 1. The D17 rule as implemented (`scripts/f1_common.py`)

Only one edit type is used: **AND → OR merges** of AND-ed *leaf* sibling genes. A connected component of qualifying pairs is merged only if it is a clique; otherwise it is logged as blocked. This is the W5 AND → OR edit.

A pair qualifies iff all three hold:
1. It is an Ensembl 116 paralog pair with min(%id a→b, %id b→a) ≥ 30.
2. **Same compartment**: the Human1 `genes.tsv` `compartments` sets are equal and non-empty. Order is ignored, and a missing annotation counts as "not same".
3. The two genes share no CORUM 5.1 / Complex Portal complex. This W5 leg is kept: D17 tightens W5's paralog rule and does not drop its complex check.

Reaction-level exclusions apply per reaction. The edit unit is the single reaction, not W5's GPR group, and every GPR reaction is examined, not only W5's M0 top-200 targets:

| Code | Exclusion | Definition |
|---|---|---|
| E1 | exchange | `reaction.boundary`, or subsystem "Exchange/demand reactions" |
| E2 | transport | subsystem "Transport reactions", or metabolites in more than one compartment |
| E3 | **pseudo/lumped** (label-free detection rule, from Human-GEM annotations) | **(a)** subsystem "Artificial reactions" or "Pool reactions", or **(b)** composition-weighted lumping: one side has ≥ 2 metabolites with non-integer coefficients that sum to 1. Example: MAR07161, 0.18 ATP + 0.3 CTP + 0.34 GTP + 0.18 UTP → RNA + PPi. Rule (b) does not flag CYOOm3i (MAR13081, 1.96 H2O + 0.02 O2⁻) or the 0.5 O2 P450 reactions. |
| E4 | large complex | some AND node of the canonical GPR spans ≥ 6 distinct genes |

**Plus the TYMS fix:** MAR04644 bounds go from (−1000, 1000) to (0, 1000), as in Human2 (DESIGN §3.1).

**Not applied (outside D17):**
- W5's OR → AND edit;
- W5's pan-cancer-absent OR-branch removals (the ER-membrane "soluble"-proxy issue raised by the W5 review).

**Interpretations for the PI to confirm before scoring:**
- AND→OR only.
- W5's complex leg is kept.
- "Same compartment" means equal `genes.tsv` sets.
- "≥ 6 genes" means distinct genes under any AND node.
- E3 is defined by rules (a) and (b).
- E2 includes multi-compartment reactions. COX7B/COX7B2 lands in E2 because complex IV moves protons across the membrane; it also fails the compartment test.
- The edit unit is the single reaction.

### Lumped list

`lumped_reactions_human1.tsv` lists the 76 E3 reactions. Nine of them carry a GPR:
- MAR07160–MAR07164 (rule b: DNA/RNA synthesis and degradation);
- MAR00010, MAR00013, MAR00014 and MAR00031 (Pool reactions).

## 2. Edits (`edits_human1.tsv`, `curated_rules_human1.tsv`)

| | count |
|---|---|
| reactions with a GPR (Human1 v1.19.0) | 8,043 |
| GPR reactions excluded (E1 / E2 / E3 / E4; one reaction can have several) | 2,670 (0 / 2,659 / 9 / 13) |
| AND-sibling Ensembl paralog pairs (reaction × pair; distinct gene pairs) | 467; 73 |
| pairs passing all three pair criteria (distinct gene pairs) | 99; 5 |
| **applied edits** (reaction × merged gene set) | **20 AND→OR in 20 reactions** (34 AND nodes) **+ 1 TYMS** |

The 20 AND→OR edits merge three gene sets:

| Merged pair | Reactions |
|---|---|
| RRM2 / RRM2B | 14 RNR reactions |
| SEC11A / SEC11C | 2 reactions (MAR06978, MAR08029) |
| CYP4F2 / CYP4F3 | 4 ω-oxidation reactions |

Evidence re-check (`evidence_check_human1.tsv`): 136/136 OK. It uses the W5 validator for paralog rows and no shared complex, plus a raw-csv re-read of min %id and of the `genes.tsv` compartment sets.

### What the rule excludes

Taken from `edits_summary.json` and `paralog_pairs_human1.tsv`.

**POLR (lumped) merges.**
- Under **W5's literal rule** (any Ensembl paralog + no shared complex), all **12 POLR gene-pair merges** (24 reaction × pair rows on MAR07161/MAR07162) are excluded by E3 (lumped). E4 (AND ≥ 6) also applies to them.
- Under the **≥ 30 % rule alone**, 2 POLR pairs would merge: POLR3B/POLR2B and POLR3B/POLR1B. Both are excluded by E3/E4. POLR3B/POLR2B passes all three pair criteria, so the reaction-level exclusion is what blocks it. POLR3B/POLR1B also fails the compartment test.
- No other lumped reaction holds a qualifying pair.

**Other exclusions.**

| Gene pair | Reactions | Excluded by |
|---|---|---|
| SLCO1A2/SLCO1B3 | 63 | E2 transport (passes all pair criteria) |
| ABCD1/ABCD2 | 5 | E2 transport (also fails the compartment test) |
| SLC7A9/SLC7A6 | 81 | E2 transport (also fails the compartment test) |
| COX7B/COX7B2 | 2 | E2 transport + E4 + compartment |

**Compared with W5a_pid30.** The compartment criterion removes HADHB/ACAA2 (46 reactions), GLB1L/GLB1 (7) and COX7B/COX7B2. HADHA/HADH (13 %) and HADHA/EHHADH (28 %) already fail the %id test.

### Label-free consequence

Six genes change their KO closure: RRM2, RRM2B, SEC11A, SEC11C, CYP4F2 and CYP4F3. `closure_changes_human1.tsv` lists them.

## 3. KO sweeps (`sweeps/`; `scripts/f1_sweeps.py`)

**Setup.**
- Base: Human1 + W3 B-loo:cofactor_pool_biomass.
- Media: all **20** Human1 M2(f=1) bounds signatures used by any line. Bounds were rebuilt with the M6 code and matched to the W2 sweep ids.
- Solver: GLPK, A10 limits (300 s per LP + 1 retry; 1,800 s CPU per job).
- 65 jobs, **0 failures**.

**Method per arm.**
- TYMS: full sweeps.
- GPR and F1: only the 6 changed-closure genes are re-solved, on the P1 base (GPR) or the TYMS base (F1). The wild-type LP does not change.

**Spot checks, all bitwise equal** (`spot_checks.tsv`):
- the P1 full sweep equals the pinned M6/W3 sweeps (2 signatures);
- the F1 full sweep equals the composed F1 (2 signatures);
- the GPR full sweep equals the composed GPR (1 signature).

**Storage.** `ratio_r_delta_human1.parquet` holds only the (arm, signature, gene) ratios that differ from P1. The arms are rebuilt exactly by `f1_sweeps.arm_ratios`, and a round-trip check is asserted.

**Label-free effects vs P1.**

| Arm | Effect |
|---|---|
| GPR | RRM2 and RRM2B go from lethal to viable on all 20 signatures. SEC11A/C and CYP4F2/3 change no ratio on any M2(f=1) medium. |
| TYMS | UMPS and DHODH go from viable to lethal on the 6 thymidine-containing signatures (DMEM:F12 family, and RPMI + HT): 65 of 1,361 lines. 88 genes change some ratio, mostly sterol-pathway shifts of ≤ 0.003. |
| F1 | GPR + TYMS: 90 genes, 52 binary flips across signatures |

## 4. Training scoring: code ready, not run on labels

Script: `scripts/f1_score_train.py`. Contrasts were pre-declared in the code before any label read.
- **Primary (TG):** F1 − P1 on L_prob, all U genes. TG passes if Δ̄_lin ≥ 0.005, the CI lower bound is > 0, and ≥ 9/12 folds are positive.
- **Components:** TYMS − P1, GPR − P1, F1 − TYMS.
- **Descriptive:** every arm − G*.
- **Label cells:** all genes, minus CE_fit (fold-scoped LOLO via `cv.lolo`), and expressed pairs.
- **Binary:** fixed-cutoff MCC / F1-score.

The dry run (`--synthetic-labels`) passes:
- P1's score matrix equals the M6 cell H1+Bloo/M2 bitwise on the 686 lines;
- no arm has a failed line.

The real run also asserts that P1's per-line tAP equals `results/m6/train/per_line.tsv`.

`scripts/f1_val_spec.py` writes `results/m10/val/look_spec_M10_F1_gpr_curation.json` and its scores. It refuses to run unless the real training summary shows that F1 passes TG. Its label-free score builder already reproduces the M6 validation P1/H1/H2 columns. `tests/test_f1.py` (16 tests, stable names) checks that the spec passes `evaluate_val.check_spec`.


## 4b. Training result (real labels; `train/`, run from the clean tree 6d8b109, MANIFEST dirty=False, 42 CPU-s)

Setup: 686 lines, 12 lineages, L_prob, all U genes; lineage-stratified bootstrap, B = 2,000. Verification passed:
- P1's score matrix equals M6 bitwise;
- P1's tAP_line equals M6 `train/per_line.tsv` (max diff 5e-7, rounding).

| Contrast | Δ̄_lin [95% CI] | folds > 0 | TG | without CE_fit: Δ̄_lin [CI], folds |
|---|---|---|---|---|
| **curated P1 (F1) − P1** | **−0.00063 [−0.00071, −0.00055]** | **0/12** | **fail** | +0.00206 [+0.00195, +0.00217], 12/12 |
| TYMS only − P1 | −0.00002 [−0.00009, +0.00005] | 3/12 | fail | +0.00009 [+0.000002, +0.00019], 5/12 |
| GPR rule only − P1 | −0.00061 [−0.00064, −0.00058] | 0/12 | fail | +0.00197 [+0.00192, +0.00202], 12/12 |
| F1 − TYMS (GPR increment) | −0.00061 [−0.00064, −0.00058] | 0/12 | fail | +0.00197 [+0.00192, +0.00202], 12/12 |

Other aggregates:
- **Expressed pairs:** F1 − P1 = −0.00046.
- **tAP_lin:** G* 0.2996, P1 0.3529, TYMS 0.3529, GPR 0.3523, F1 0.3523.
- **tAP_lin without CE_fit:** P1 0.1861, F1 0.1882.
- **Descriptive vs G\*:** P1 +0.0534, F1 +0.0527.
- **Fixed-cutoff binary vs P1, L_ge06:** MCC −0.0016 [−0.0017, −0.0016]; F1-score −0.0032.
- **TYMS binary vs P1:** MCC +0.00002; F1-score +0.0001.

Reading:
- The GPR rule changes only RRM2/RRM2B. RRM2 goes from lethal to viable on every medium, which costs about 0.0006 on all genes. The gain without CE_fit is consistent (12/12 folds) but small (+0.002), below the 0.005 SESOI.
- The TYMS fix acts on only 65 of 1,361 lines, the thymidine media, and is null on P1.
- The stricter rule therefore removes the W5 POLR damage (W5a −0.029). With an all-genes Δ of −0.0006, it neither passes TG nor helps on all genes.

**Decision:** no validation spec is prepared. `scripts/f1_val_spec.py` refuses to run on this summary, and nothing was written under `results/m10/`.

## 5. Compute

Development runs are listed in `MANIFEST.tsv`; all are dirty because the code is uncommitted.

| Run | CPU-h |
|---|---|
| f1_edits ×2 | 0.013 |
| f1_sweeps (365 s wall, 8 workers) | 0.45 (1,609 CPU-s) |
| synthetic score dry run | 0.005 |
| **Total (development)** | **≈ 0.47** |
| clean re-runs (edits + sweeps, from f57bffe) | 0.46 |
| training scoring (6d8b109) | 0.012 |
| **Total F1** | **≈ 0.95** |

No single run exceeded 1 CPU-h. No Gurobi LP was solved.

## Files

| File | Content |
|---|---|
| `edits_human1.tsv` | applied edits with evidence IDs (20 AND→OR + TYMS) |
| `curated_rules_human1.tsv` | old/new GPR of the 20 edited reactions |
| `paralog_pairs_human1.tsv` | every AND-sibling paralog pair in every reaction: pair criteria, reaction exclusions, applied |
| `excluded_reactions_human1.tsv` | GPR reactions excluded by E1–E4 (all reasons) |
| `lumped_reactions_human1.tsv` | the E3 list (76 reactions; 9 with a GPR) |
| `evidence_check_human1.tsv` | raw-file evidence re-check |
| `closure_changes_human1.tsv` | the 6 genes whose KO closure changes |
| `edits_summary.json` | counts, including exclusions by reason and POLR vs other |
| `sweeps/` | ratio deltas vs P1, job table, spot checks, summary |
| `MANIFEST.tsv` | provenance rows |
