# Design review round 1: algorithms and compute

Reviewer: CBM algorithms / HPC reviewer (not the author). Reviewed DESIGN.md v0.1, PROMPT.md, CLAUDE.md, results/baseline/README.md, src/nextgem/, the installed `.venv`, and data/Human-GEM-main/code.
Scope: whether the algorithms are correctly specified and feasible, and what they cost. I did not compute any essentiality performance against DepMap. The only DepMap file read was expression, for 12 lines (first rows of the file), as a label-free feasibility probe.

## 0. Verdict

The design is sound in structure, but **as written it cannot run in a multi-day budget**, and three algorithm specifications fail on contact with Human2:
1. The medium ladder M1/M2 gives WT no-growth.
2. F0/F5 pruning without repair gives WT no-growth in 12/12 lines tested.
3. ftINIT is specified in a different configuration from its positive control, and its cost is under-estimated about 3–10×.

All of these can be fixed inside the PROMPT scope. With the fixes, my estimate is **≈ 2,200–2,700 CPU-h, about 70–85 h of compute wall-clock**. Of that, about 1,000 CPU-h is Gurobi MILP time bound by the 16 licence sessions, and the rest is LP work that can run on GLPK on the other cores. That fits an 8-day plan only if per-line ftINIT is built once, for at most one configuration.

## 1. What is installed vs MATLAB-only (checked)

| Item | Status |
|---|---|
| `raven-toolbox 3.0.0b1` (.venv) | **Has ftINIT**: `raven_toolbox.init.{prep_init_model, ftinit, run_ftinit, fill_tasks, gene_scores_from_expression, score_reactions_from_genes, remove_low_score_genes, get_init_steps}`. Series available: `1+0`, `1+1`, `2+0`, `2+1`, `full`. Extras not in RAVEN: `resolve_ties`, `prove_abs_gap`. No loopless option (documented). Task tools: `raven_toolbox.tasks.check`, `tasks.tasklist`. |
| Human-GEM Python ports (`code/test/`) | `estimateEssentialGenes.py` ports `prepHumanModelForftINIT` (drug and AA-triplet reaction removal, MAR13081 removal, spontaneous/custom-ignore lists) and runs ftINIT `1+1`, `removeGenes`, `useScoresForTasks`, with `resolve_ties=True`, `prove_abs_gap=1.0`, `time_limit=900`. Also ported: `taskEssentialGenes.py` (task-essential genes with distinct-KO dedup and process pool) and `gradedEssentiality.py` (sharding, `--prep-cache`, `--model-dir`, checkpointing, Gurobi-session handling). **Reuse these rather than re-porting.** |
| MATLAB only | `code/tINIT/*.m` (already ported in Python, above). `code/DepMapGeneEss/generate_DepMap_models_ftINIT.m`: the DepMap protocol is **`1+0` by default**, `removeGenes=false`, `useScoresForTasks=true`. `evaluateDepMapEssentialityPredictions.m`, `getTaskEssentialGenes.m`, `checkTasksGenes.m`. None is needed if the Python ports are used. |
| Repository-measured ftINIT cost | `estimateEssentialGenes.py` header: prep ≈ 30 min once per template. Five Hart lines: summed **92 min** without ties, **124 min** with `resolve_ties` (HCT116 7→22 min, RPE1 38→59 min). That is about 18–25 CPU-min per line, with a tail of about 1 h. |
| cobra 0.32.1 | `moma`, `add_moma`, `room`, `add_room`, `loopless_solution`, `add_loopless`, `fastcc`, `gapfill`. No GIMME, iMAT or FASTCORE. |
| Not installed (PyPI reachable) | `corda 0.5.1`: pure Python, needs only `cobra>=0.25`, should install cleanly. `troppo 0.0.7`: pins `cobamp==0.2.1`, an old stack, risky next to numpy 2.4. `memote 0.17.0`: pins `pandera~=0.18`, `PyGithub~=1.55`, `jsonschema~=4.20`, so **use a separate venv**. FROG (`fbc_curation`) is not installed. HiGHS (`highspy`) is not installed; optlang's `hybrid_interface` could use it. |

## 2. Timing probes (Human2 v2.1.1, Ham's M0, Gurobi 1 thread, this host)

Scripts and JSON are in `/tmp/claude-10105/-data-tmpA-nextgem/c5c0f45b-d54c-4e57-bdef-3d0657510b12/scratchpad/review_r1_compute/` (`probe1.py`, `probe2.py`, `probe1.json`, `probe2_12.json`).

| Quantity | Measured |
|---|---|
| SBML load | 18.9 s |
| `model.copy()` | 15.1 s |
| pickle load | 14.9 s (15 MB) |
| WT LP | 1.1 s (Gurobi), 5.4 s (GLPK, cold) |
| pFBA | 0.5 s |
| Genes whose KO closes ≥ 1 reaction | 1,128 / 2,870 |
| Genes whose closed reactions carry pFBA flux (must be solved under the zero-flux shortcut) | **320** |
| Per-KO time on those 320 genes, Gurobi | 0.01–3.3 s, mean ≈ 0.48 s (n = 15) |
| Same KOs, GLPK warm | 0.04–1.1 s, mean ≈ 0.33 s (n = 10) |
| Baseline full sweep (README), summed | 209 s; most of it is in the ~320 flux-carrying genes |
| **`cobra.flux_analysis.moma(linear=True)` per KO** | **≈ 62 s** (it rebuilds the problem on every call) |
| SBML write | 8.2 s, 50.7 MB raw, **2.14 MB gzipped** (full Human2; context models are smaller) |
| RSS of one process holding Human2 plus one copy | ≈ 1.45 GB |
| Ham's drop-one-input (biomass ratio = 0 when dropped) | linoleate, linolenate, α-tocopherol, γ-tocopherol, retinoate, lipoic acid, aquacob(III)alamin, thiamin, riboflavin, pyridoxine, pantothenate, folate, biotin, Fe2+, Pi, O2, and the 9 essential amino acids. Thymidine drop gives 0.63 (it is a bulk carbon source at bound 1000). |
| F0 probe: close every reaction whose GPR is False with TPM < 1, 12 DepMap lines | **WT growth = 0 in 12/12.** Each line switches off 5–9 generic-lethal genes: CYP4F2, CYP26A1/B1, CYP24A1, CYP1A2, BBOX1, TMLHE, ACADL. |

Implications:
- **The zero-flux shortcut is exact for FBA but saves almost nothing**, about 5 s of 209 s, because the cheap genes were already cheap. Do not budget savings from it.
- The real savings come from three things:
  - (a) never copying, unpickling or re-reading a model per line;
  - (b) running LPs on GLPK on non-licensed cores;
  - (c) caching and post-processing, described below.

## 3. Must-fix

**MF1. The medium ladder is infeasible as specified (W2 M1, M2a, M2b).** Human2 biomass needs linoleate, linolenate, α/γ-tocopherol, retinoate, lipoic acid and cobalamin.
- RPMI and DMEM lack the lipids, tocopherols, retinoate and lipoate.
- DMEM also lacks B12.
- The serum lists S1 and S2 contain no tocopherol, retinoate or lipoate.

So M1, M2a and M2b give WT no-growth for essentially every line. By §6.1 every line then falls back to the parent and the method "fails". Changes:
- (i) Add to `serum_model.tsv` a fixed **"biomass micronutrient floor"** (α/γ-tocopherol, retinoate or retinol, lipoate, linoleate, linolenate, cobalamin), sourced as serum or albumin carriers. It enters S1 and S2.
- (ii) M1 (no serum) is then expected to be non-growing for most media. Report it as a pre-registered negative step, or define M1 = base medium plus the floor at the serum bound. Decide this now.
- (iii) Add a mandatory **WT-growth precheck per (medium × model)** before any sweep. Its output is a table, and the precheck is a gate.
- (iv) Run the same drop-one-input analysis on Human1, whose biomass requirements may differ.

**MF2. F0 and F5 (and the genomic constraint applied "to every family") need an a priori repair rule (W4).** Pure GPR pruning killed WT in 12/12 lines. The culprits are biomass-artefact lethals: vitamin/xenobiotic CYPs and carnitine synthesis. Generic FBA calls these genes lethal, yet they are unexpressed in most lines. Specify one of the following, fixed now:
- (a) **Minimal re-opening.** After pruning, solve one LP or MILP that re-opens the cheapest set of closed reactions (cost = expression deficit) needed for WT growth ≥ 1% of the generic value, and log the reactions per line. This is GIMME-like and costs about 1–2 LPs.
- (b) **A frozen protection list.** Never close reactions required for WT growth on the selected medium in the parent.

Whichever is chosen, the same rule applies to the genomic constraint, and the "failed model" rate is reported per family. Also flag these genes in the W1 failure account ("biomass-artefact lethal"), and tie them to W3's B-loo_j: retinoate and tocopherol components plausibly explain several generic false positives.

**MF3. ftINIT (F1) specification, cost and acceptance test (W4 F1, §10, §12, M3).**
- (a) **Configuration.** F1 says "Human-GEM protocol 1+0". That matches the MATLAB DepMap script (`1+0`, removeGenes = false). But the positive control "reproduces the repository AUPRC 0.25–0.33". Those numbers (`data/testResults/gene-essential_summary.md`) come from the **Python `1+1`, removeGenes = true, resolve_ties = true, prove_abs_gap = 1, time_limit = 900** pipeline. Pin one F1 configuration, with every parameter listed: series, removeGenes, useScoresForTasks, resolve_ties, prove_abs_gap, time_limit, mip_gap_abs, big_m, threshold. Run the Hart control **with that same configuration**, plus the repository configuration as a reference.
- (b) **Acceptance.** With `resolve_ties` the build is deterministic, so ±0.05 AUPRC is far too loose. Require that running the repository's own `gradedEssentiality.py` and our wrapper on the same prep gives per-line essential-gene Jaccard ≥ 0.95 and growth-AUPRC within 0.01. Our F1 config vs the repository config is then a reported difference, not a pass/fail.
- (c) **Cost.** Use the repository-measured 18–25 CPU-min per line (tail about 60 min) plus a 30 min prep per (template × task list), not 6 min. As written, the M3 gate ("a family is kept only if ≤ 15 CPU-min per model") would **exclude ftINIT outright**, contradicting the overrun rule in §12. Replace the per-model cap with a per-family **total budget**: F1 ≤ 700 CPU-h overall. F1 runs on the variance panel, on group models, and per line for at most one configuration.
- (d) **Medium.** ftINIT protects the 57 tasks, which use their own inputs, so the extracted model is not guaranteed to grow on the line's documented medium. Add a "growth on medium m" task to `fill_tasks`: one prep per distinct medium, about 15 × 30 min, or gap-fill only. Otherwise the 5% failure rule will trip.
- (e) **Stability study.** With `resolve_ties`, seeds barely matter. Measure stability as (i) `resolve_ties` off × SEEDS[0:3], and (ii) a perturbation of gene scores by expression noise. The raven-docs note that equal optima differ by about 10% of essential calls.

**MF4. The W4 run matrix contradicts §12 and is about 6–9× over budget.** "Each family × threshold × modality runs at the per-line level on all training lines (LOLO for F5 and R3)" multiplies out to:
- 6 families × 3 thresholds × (RNA ± genomic) × ~950 lines ≈ **34,000 builds**, plus the protein arms (ProCan, 949 lines, has now arrived);
- the F5 grid: 48 × 950 ≈ **45,600 sweeps**.

That is about 8,000–12,000 CPU-h, about 1,900 of it ftINIT MILP alone (119 h on 16 sessions). §12 budgets 5,400 builds. Replace it with a staged plan:
- **Stage A (screen), on the 100-line variance panel, which is already built:** choose ≤ 2 thresholds per family and the modality, using panel-only training mAP_line. Builds for the panel's proteomic lines are added for the protein arm.
  - R3 needs **no LOLO**: it is a label-free expression statistic, so compute it once from training lines.
  - F5 grid: run it on the panel (48 × 100), or as a coarse 12-point grid plus local refinement.
- **Stage B:** at most 4 configurations go per line on all screened lines, with at most 1 of them ftINIT.
- **Stage C:** group models (lineage, background, combination; about 150 models).
- **Go/no-go for per-line F1:** run it only if, on the panel, per-line F1 beats per-lineage F1 on training. Otherwise per-lineage F1 is the F1 deliverable. **TG for MILP families** is then computed on panel lines (state this).

**MF5. MOMA is specified in a way that costs about 5.5 CPU-h per model (W3 O-MOMA, §6.1).** `cobra.moma(linear=True)` costs about 62 s per KO. Over 320 flux-carrying genes that is about 5.5 CPU-h per model. Required changes:
- Use `add_moma(model, solution=ref, linear=True)` **once per (model, reference)**, then toggle KO bounds. Estimate: ≤ 2 s per KO, about 10 CPU-min per model.
- Linear MOMA is degenerate, so its biomass flux is not unique and gives a noisy ratio. Either add a lexicographic second step (fix the MOMA distance, then maximise biomass) or use QP MOMA (Gurobi QP, unique optimum).
- Define the reference as pFBA WT flux on the same medium.
- Verify the zero-flux shortcut separately for MOMA. It is valid only if the check uses the MOMA reference vector, in which case the distance is 0 and the ratio is 1.

**MF6. The task-as-constraint and NGAM modes are ill-posed as written (W3).**
- **C-T modes.** RAVEN tasks are closed-system tests with their own inputs. SU/IC tasks are conversions with fixed input fluxes. They cannot be added as "output ≥ 1% of WT max" constraints to a growth LP on a culture medium. Restrict C-modes to tasks expressible as a **product demand**: BS tasks, plus ER as an ATP/NAD(P)H demand. Implement each as a demand reaction with lb = 1% of its max under the same medium with biomass free. Add all such demands jointly, after checking that WT is feasible with all of them (drop and log infeasible ones). Note also that a 1% floor makes the C-modes nearly equivalent to the criterion modes for lethal calls; say so, and keep the C-mode only as the continuous score.
- **O-NGAM.** "Fraction of WT ATP-hydrolysis flux" is ill-defined: in an FBA optimum, flux through MAR03964 is often 0 or arbitrary (alternate optima). Define NGAM_k = k × the **maximum** MAR03964 flux on the medium (computed once), and verify the reaction ID in both models.

## 4. Should-fix

**SF1. Engine efficiency (§10 KO engine, §12).**
- Keep one parent model per worker for the whole run. Build each line's model by bound changes inside `with model:`, which costs about 1 s instead of 15–19 s for copy, unpickle or SBML read.
- Dedupe genes by their closed-reaction set (as `taskEssentialGenes.py` does).
- Store context models as reaction-ID sets plus GPR diffs, **never pickles**: 15 MB × thousands of models would be tens of GB on a 98%-full NFS. The repository's `--model-dir` pickling is fine only for ≤ 50 models.
- Keep the zero-flux shortcut, since it is exact, but drop its claimed budget effect.

**SF2. Licence and memory plan (§12).**
- Reserve Gurobi's 16 WLS sessions for the MILPs (ftINIT, F5 if MILP, CORDA if it needs one).
- Run LP KO sweeps on **GLPK**: comparable per-KO time in the probe, and already verified equal to Gurobi on 200 KOs. HiGHS is a possible alternative. Re-verify ratio_r equality on 1 full sweep per family. Fix the solver per method, so that ratio_r rounding edge cases cannot differ between runs.
- RAM: 68 GB was available at review time (57 GB used by others). Budget about 1.5 GB per LP worker. Measure ftINIT RSS in the M3 pilot; deepcopy-heavy, expect 3–4 GB. Cap the total at about 50 GB, e.g. 12 MILP + 16 LP workers, adjusted after the pilot.
- WLS needs internet at solve time, so use per-line checkpointing (the `gradedEssentiality.py --checkpoint-dir` pattern) and session-acquire retries. Do not probe more than 16 sessions without PI approval, because the licence is shared.

**SF3. One task capacity matrix per model (W3).**
- Compute a gene × task matrix of max-output capacity once per (parent, medium) for the 57 Essential and 195 CellFie tasks. Use dedup and a per-task zero-flux shortcut against each task's own solution. Estimate: about 2–5 CPU-h for both models.
- T0–T5, T4_k, O-cap and Negative-2 then become **post-processing with zero extra LPs**. T5 on generic models is a per-line lookup over that matrix.
- Random pseudo-tasks need their own matrices: 5 seeds × 11 groups, small.

**SF4. Local re-evaluation in W5.** A GPR edit on reaction r changes only the KOs of genes in r's GPR. Re-solve only those genes, with dedup. Enumeration then costs about 2–5 CPU-h, not 70.
- In the 2 × 2 factorial, **apply curated GPRs to the already-extracted context models** (fixed reaction sets).
- Re-extracting ftINIT with curated scores confounds Δ_curation with MILP alternate-optimum noise; the raven-docs report about 10% of calls differ. Treat template-level re-extraction as a sensitivity analysis, run only for P1.

**SF5. The M6 combination factorial (§9.1) costs 24 × 950 sweeps, ≈ 760 CPU-h.** Run the ≤ 24 combinations on the 100-line panel and take the top 3 to all training lines, for about 175 CPU-h. Medium, tasks and curation are post-extraction, so cache the extractions and re-sweep only.

**SF6. W6 "knockout flux profiles".** X_KO = the e-score vector is nearly binary and mostly identical across lines (about 150 lethal genes, ratio_r bimodal). PROMPT asks about "knockout flux profiles". Add X_KOflux: for an a priori list of about 50 genes, namely the metabolic-drug targets fixed in W6.2, store the pFBA flux of the KO model (or MOMA), normalised by biomass flux, then PCA. Cost: about 50 × 1,361 LPs, about 20 CPU-h. Normalise X_flux per line by biomass flux, since units are arbitrary at bound 1000. pFBA is loop-free by construction, so no loopless step is needed.

**SF7. W7 H-c (flux-cone learning).** The pilot threshold of "≤ 2 CPU-h per line" implies up to 2,700 CPU-h, 9× the W7 cap. Restrict H-c to the generic and group models (≤ 150), or cap it at 100 CPU-h.

**SF8. W8 resource.**
- **Reuse the W4/T1 builds**: cache by (config hash, input hash). Build only the about 770 lines that were not screened. That is about 320 CPU-h if P1 is ftINIT, versus about 875 for rebuilding all 2,100. W8 can start right after T1, overlapping W7, because P1 is frozen.
- Disk: about 2.1 MB gzipped per full Human2 model was measured, so the resource is ≤ 4.5 GB, not 10 GB.
- Gates:
  - Run "libSBML 0 errors" on the **parent** first. If Human-GEM itself emits errors, gate on "no errors beyond the parent's".
  - memote: separate venv; 10 models.
  - FROG: run `fbc_curation` on 10 models only. FROG runs FVA over all reactions, about 30 CPU-min per model.
- Validation (write, read and 20 KOs) costs about 1.5 CPU-min per model, about 55 CPU-h.

**SF9. Extraction family specs (W4).**
- **F4 is E-Flux (Colijn 2009).** Name it, and specify:
  - GPR-less reactions keep their default bounds;
  - exchanges stay at the medium bounds, not the expression bounds;
  - with internal bounds at least 10 and uptakes at 1000, growth is set by the scaled internal bounds — state this as intended.
- **F2 GIMME:** fix its required-function fraction (e.g. 0.9 × max biomass) and its expression cutoff for each threshold rule.
- **F3 FASTCORE:** state where the code comes from. troppo is a pinning risk, so prefer your own toy-tested implementation using cobra `fastcc`. **CORDA (F3b)** installs cleanly, but its runtime on Human-GEM is unmeasured (likely 30–120 min). Put it in the M3 pilot with a cap.

## 5. Optional

- **Loopless.** Not needed for KO ratios: loops cannot change maximal biomass. Document this so nobody adds `loopless_solution` to sweeps. The ftINIT loop caveat is documented upstream.
- **ftINIT scores.** These are 5·ln(TPM/threshold). R2 and R3 plug in as the reference level. Write down the exact call.
- **H1 replication of ftINIT.** It needs its own prep, about 30 min. Run H1 ftINIT only on the panel unless P1 is F1.
- **W2 shuffled-medium control.** It needs **no new sweeps**: generic results are cached per medium, so shuffling is a relabelling. Count it as 0 CPU-h.

## 6. Compute estimate

Unit costs:

| Unit | CPU time | Basis |
|---|---|---|
| LP sweep, full Human2 | 3.5 CPU-min | baseline |
| LP sweep, context model | ≈ 2–4 CPU-min | estimate |
| ftINIT | 20–25 CPU-min (+ 30 min prep per template/task set) | repository |
| GIMME / E-Flux build | < 1 min | estimate |
| FASTCORE | ≈ 5 min | estimate |
| CORDA | 30–120 min | unmeasured |

Totals:

| Block | DESIGN §12 | Literal §W4 text | **Proposed (with fixes)** | of which Gurobi MILP |
|---|---|---|---|---|
| W1 + W2 (per medium, incl. MF1 precheck) | 16 | 16 | 20 | 0 |
| W3 (capacity matrix, add_moma, NGAM, B-loo, pseudo-tasks) | 50 | 50 | 20 | 0 |
| W5 (local re-eval + factorial post-extraction) | 70 | 70 | 55 | 0 |
| W4 prep + Hart control + pilots | 5 | 5 | 20 | 10 |
| W4 variance panel (2,000 builds; 400 ftINIT) | 130 | 130 | 230 | 135 |
| W4 Stage A screen (protein arm, F5 on panel) | — | ~4,000 | 130 | 20 |
| W4 per-line (≤ 4 configs × 1,361, 1 ftINIT) | 540 | ~6,000 | 730 | 455 |
| W4 group models (≈ 150 × 4) + stability (60 ftINIT) | 30 | 30 | 100 | 70 |
| M6 factorial | — (missing) | 760 | 175 | 0 |
| W6 | 20 | 20 | 40 | 0 |
| W7 | 300 | 300 (+ H-c ≤ 2,700) | 300 | ≤ 50 |
| W8 (≈ 770 new builds + validation of 2,100) | 245 | 245 | 380 | 320 |
| Re-runs and bug fixes (20%) | — | — | 450 | 210 |
| **Total CPU-h** | **≈ 1,400** | **≈ 12,000** | **≈ 2,650** | **≈ 1,270** |

**Wall-clock.**
- **Gurobi MILP:** 1,270 CPU-h / 14 effective sessions ≈ **90 h**.
- **LP:** about 1,400 CPU-h / 24 GLPK cores ≈ **58 h**, running in parallel with the MILP work.

So compute is about 3.5–4 days of pure wall time. The critical path is ftINIT:

M3 prep and Hart control (≈ 1.5 h) → panel F1 (≈ 10 h) → H4 verdict and Stage A → per-line F1 (≈ 33 h) → M6 and T1 → W8 F1 builds (≈ 23 h), overlapping W7 (≤ 20 h).

This fits D+8 only if per-line F1 runs once.

If the panel shows that per-line F1 is not better than per-lineage F1 (MF4 go/no-go), the per-line and W8 F1 blocks drop out. That removes about 775 MILP CPU-h, leaving a total of about 1,700 CPU-h and about 50 h wall.

**Disk:**
- KO parquets: ≤ 2 GB.
- Reaction-set and GPR-diff storage: < 0.5 GB.
- SBML resource: ≤ 4.5 GB.
- Checkpoints: < 1 GB.

Total ≤ 10 GB, well under the 30 GB ceiling, as long as no per-line pickles are written (SF1).

**Licence:** 16 sessions verified. Plan for 14 MILP sessions plus 2 for interactive and verification use, and run LPs on GLPK. Do not exceed 16 without approval.

## 7. DESIGN.md edits requested (by section)

- §2.3 / §6.1: a failed context model's fallback must follow the MF2 repair rule. Report the failure rate per family and threshold.
- §8 W2: MF1 (micronutrient floor, WT-growth precheck, M1 redefinition); the shuffled control costs 0 CPU-h.
- §8 W3: MF5 (add_moma, uniqueness), MF6 (C-modes as product demands; NGAM from max capacity), SF3 (capacity matrix).
- §8 W4: MF2, MF3 (pinned F1 config, Hart control on the same config, medium task, stability definition), MF4 (staged A/B/C; R3 without LOLO; F5 on the panel; TG on the panel for MILP families), SF9.
- §8 W5: SF4. §8 W6: SF6. §8 W7: SF7. §8 W8: SF8.
- §9.1: SF5 (factorial on the panel).
- §10: add the exactness test for MOMA add_moma vs moma (5 genes), the GLPK = Gurobi ratio_r test per family, and the medium WT-growth precheck table.
- §12: replace the whole table with §6 above. Replace the per-model 15-min gate with per-family totals; add memory and session allocation.
- §13 M3: add the pilots for ftINIT (`1+0` vs `1+1`, RSS, time), CORDA, the GLPK sweep, add_moma, and the F0 repair rule.
