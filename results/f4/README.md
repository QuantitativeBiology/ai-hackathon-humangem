# F4: line-specific constraints beyond the medium, on top of P1

**Post hoc / exploratory (A15).** This follow-up comes after the T2 unseal (DESIGN §17 A15; RESEARCH_TREE D19). It never replaces T1/T2 or the headline, and it gets no test evaluation. The F4 agent did not open any T1/T2 output, `results/m6/test/`, `results/f3/` or `data/locked/`.

Author: F4 engineer agent (Claude Opus), 2026-10-28. **Status: complete. Both arms fail.**
- **F4a (genomic loss)** fails TG on training: Δ̄ = −0.00013 [−0.00026, −0.00002], 3/12 folds.
- **F4b (E-Flux on P1)** fails its panel pre-check: Δ̄ = −0.0007 [−0.0017, +0.0002], 4/10 folds, below +0.005. It stops there, as D19 requires.

No validation spec is prepared, and nothing was written under `results/m10/`.

## Provenance

- **Code:** committed before any label read (the chain before b62bd80): `scripts/f4_common.py`, `f4_losses.py`, `f4_sweeps.py`, `f4_score_train.py`, `f4_eflux_build.py`, `f4_eflux_score.py`; `tests/test_f4.py` (9 tests pass).
- **Production runs:** all ran from the pinned clean worktree `/data/tmpA/nextgem_runs/b62bd80`. Their MANIFEST rows have dirty = False.
- **Development row:** the first `f4_losses.py` row in MANIFEST.tsv (dirty = True) is a development run of the same code. Its output is identical to the clean run except for the per-line `cpu_s` column, and the clean outputs replaced it.
- **Solver and limits:** GLPK only, at most 12 workers, BLAS at 1 thread. A10 limits: 300 s per LP + 1 retry; 1,800 s CPU per job.
- **Failed builds and timeouts:** 0.

## Label hygiene

- **Probe:** the session-start probe `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label hook.
- **Label access:** labels were read only by `f4_eflux_score.py` and `f4_score_train.py`.
  - Both read through `labels.train()` and `cv.lolo` / `labels.fold()` (fold-scoped CE_fit).
  - Training tier only; no validation, test or `data/locked/` access.
  - `f4_score_train.py` was first dry-run on synthetic Bernoulli labels (`tmp/score_dry`); its checks passed.
- **Label-copy scanner:** 0 hits over `results/f4` and `scripts`.
- **Training aggregates viewed by the author:**
  - the contrast tables below (Δ̄_lin, CIs, folds positive, the minus-CE_fit and expressed-pairs cells, the changed-lines subset);
  - per-arm lineage-equal tAP_lin, mΔAP_gene and MCC;
  - the binary paired Δ;
  - the F4b per-fold (panel lineage) mean Δ;
  - the verification max |Δ|.
- **Not opened by the author:** `train/per_line.tsv` and `eflux_panel/train/per_line.tsv`.
- **Exposure carried in:** the M5, F1 and RESEARCH_TREE training aggregates, including the one-line F3 summary in RESEARCH_TREE.

## F4a: genomic loss on P1

**Rule (D19 = DESIGN §W4 loss rule).**
- **Lost gene (rule r2-11):** damaging code 2; or code 1 and (relative CN < 0.5 or log2(TPM+1) < 1); or relative CN < 0.3. This is `context_inputs.lost_genes`.
- **Acceptance:** the order-independent drop `context.accept_losses` runs on the P1 parent. That parent is Human1 + B-loo on the line's M2(f=1) bounds; the signature comes from `results/m6/scores/line_assignment.tsv`, which equals W2's `line_sweeps.tsv` (asserted).
- **Closure:** the accepted genes are knocked out, so their GPR closure (`context.loss_closure`) is closed.
- **KO sweep:** re-run with the P1 engine (`w3.ko_sweep`) under the A4 convention:
  - a lost gene is non-functional everywhere;
  - the KO of its isozyme partner closes the shared reaction;
  - the lost gene's own KO changes nothing.
- **Repair:** none. The generic parent has nothing closed by extraction, so the repair has nothing to repair.

### Label-free loss counts (`losses/`)
686 training lines in 12 lineages, on 15 M2(f=1) signatures.

| | lines |
|---|---|
| ≥ 1 lost model gene | **665** |
| ≥ 1 accepted loss | 665 |
| the accepted losses close ≥ 1 reaction | 544 |
| ≥ 1 model-contradicted loss (lost, but lethal alone in the model) | 15 |
| needing the expression-ordered joint drop | 0 |
| no CN (mutations only, flagged) | 32 |

Further counts:
- **Lost genes per line:** median 4, mean 7.2, max 95.
- **Closed reactions per line:** median 2, max 258.
- **Distinct lost genes:** 1,302.

Most often lost (lines): USP9Y 473, ACOT1 305, UGT2B17 236, MTAP 90, PTEN 53, STS 32, APOB 24, HACD4 23, MPC1L 21, USP26 20 (`losses/gene_counts.tsv`). Several of these close no reaction alone; USP9Y and UGT2B17 are among them.

**Model-contradicted losses** (an error-analysis output, as in DESIGN §W4):
- FASN ×2 and NUP43 ×2;
- DTYMK, TYMS, PFAS, PSAT1 and TECR;
- NUP160, NUP188, NUP205 and NUP210;
- POLR1E, POLR1H, POLR2E, POLR2J, POLR2L, POLR3K and POLRMT.

In every one of the 15 signatures, the P1 parent's WT growth equals the pinned P1 sweep's.

### Sweeps (`sweeps/`)
- **Jobs:** 610, one per distinct (signature, accepted set): 514 full sweeps and 96 partial ones. Partial jobs close no reaction, so only the genes whose KO closure changes are re-solved.
- **Status:** 0 failed. 23,123 CPU-s, 3,018 s wall on 8 workers.
- **Checks:** all pass, bitwise (`spot_checks.tsv`).
  - An empty-loss job equals the pinned P1 sweep on 2 signatures.
  - A partial job swept in full equals the composed partial job.
  - On 2 jobs, a cross-check against the W4 engine (`context.finish` with the accepted losses + `context.ko_sweep`) gives max |Δratio_r| = 0.
- **Label-free effect:**
  - Only **95/610 jobs (97 lines)** change any KO ratio vs P1, across 178 genes. In most lines the lost reactions carry no flux that matters for any KO.
  - **22 (job, gene) flips to lethal:** CTPS1 ×5, SGPL1 ×4, POLR3GL ×2, UAP1 ×2, GFPT1, TXNRD1, TXNRD3, CMPK1, PTDSS1, PISD, UBAC2, LCAT, PCYT2.
  - Most of these are isozyme partners of a lost paralog.
  - **0 flips to viable.**

### Training result (`train/`; 686 lines, 12 lineages, L_prob, all U genes, lineage-stratified bootstrap B = 2,000)

| Contrast | Δ̄_lin [95% CI] | folds > 0 | TG |
|---|---|---|---|
| **F4a − P1** | **−0.00013 [−0.00026, −0.00002]** | **3/12** | **fail** |
| F4a − P1, minus CE_fit | −0.00015 [−0.00030, −0.00003] | 3/12 | — |
| F4a − P1, expressed pairs | −0.00012 [−0.00024, −0.00001] | 3/12 | — |
| F4a − P1, 97 changed lines only (descriptive) | −0.00085 [−0.00174, −0.00005] | 2/12 | — |
| P1 − G* (descriptive) | +0.0534 [+0.0523, +0.0545] | 12/12 | — |
| F4a − G* (descriptive) | +0.0532 [+0.0522, +0.0544] | 12/12 | — |

Other aggregates:
- **tAP_lin:** G* 0.2996, P1 0.3529, F4a 0.3528.
- **Fixed-cutoff binary F4a − P1:** MCC +0.00002 (L_ge06 and L_prob, CIs include 0); F1-score +0.00004.

Verification (`train/verification.tsv`), all pass:
- the P1 score matrix equals the M6 cell H1+Bloo/M2 bitwise;
- F4a equals P1 on every line without an accepted loss;
- P1's tAP_line equals the M6 train `per_line` (max diff 5e-7).

**Reading:** the pre-registered loss rule barely changes the KO landscape on P1, and the small change it makes is slightly negative. The new lethal calls, such as CTPS1 when CTPS2 is lost, are paralog-loss predictions. On these training labels they do not add precision on average. F4a stays out of the method.

## F4b: E-Flux on P1 (30-line panel pre-check)

**Builds (`eflux_panel/builds/`, label-free).**
- **Configuration:** Stage A's carried F4 configuration, F4|R1|rna|real|none, built on the P1 parent with each panel line's M2(f=1) medium. The W4 builders are unchanged (`context.extract` / `finish` / `ko_sweep` / `task_profile`), one fresh process per build.
- **Run:** 30 F4b builds plus 4 generic P1 builds (one per panel signature). 0 failed, 0 timeouts. Mean cost 126 CPU-s per build.
- **Generic check:** the W4-engine P1 builds equal the pinned P1 sweeps bitwise (4/4).

**Gate statistics.** The known E-Flux units artefact (results/w4, results/m5) recurs:
- WT growth is 0.011–0.026 of the parent, i.e. 0.50–1.23 in scaled units; 25/30 builds grow below 1.
- The GR task (biomass ≥ 1) fails in **22/30** builds. No build fails an ER task.
- **gate_all passes in 8/30 builds**, so the configuration fails the gate.
- Tasks passed: 51–53 of 57.
- The lethal set equals the line's P1 lethal set in **30/30** builds (103–113 genes). F4b's signal is in its partial ratios only.

**Pre-check (`eflux_panel/train/`).** This uses the M5 / Stage-A nested panel LOLO with one fixed configuration, so the outer Δ equals the plain Δ (asserted).

| Contrast | Δ̄_lin [95% CI] | folds > 0 | minus CE_fit | expressed pairs | ΔMCC |
|---|---|---|---|---|---|
| **F4b − P1 (decision)** | **−0.0007 [−0.0017, +0.0002]** | **4/10** | −0.0069 | −0.0004 | 0.0000 |
| P1 − G* (sanity) | +0.0523 [+0.0479, +0.0568] | 10/10 | +0.0394 | +0.0480 | +0.0553 |
| F4b − G* | +0.0516 [+0.0469, +0.0564] | 10/10 | +0.0325 | +0.0476 | +0.0553 |

Other aggregates:
- **Panel tAP_lin:** G* 0.2982, P1 0.3505, F4b 0.3498.
- **mΔAP_gene:** P1 −0.0013, F4b +0.0007.
- **Fold means of F4b − P1:** −0.0041 (Esophagus/Stomach) to +0.0037 (Head and Neck).

Verification: P1's panel tAP equals the M6 per_line, and generic H1 equals W1 (max diff 5e-7).

**Decision (D19):** Δ̄ = −0.0007 < +0.005, so **F4b does not proceed to all training lines**.

Note the contrast with M5. On H1 + B-loo / M0, E-Flux was +0.0044 vs its parent. On P1, which uses the line media, that small gain disappears.

## Compute

| Run | CPU-h |
|---|---|
| Development (loss run; loss, E-Flux and sweep pilots; dry runs; under `tmp/` except the dirty loss row) | ≈ 0.45 |
| f4_losses (clean) | 0.14 |
| f4_eflux_build (clean, 4 workers) | 1.40 |
| f4_sweeps (clean, 8 workers) | 6.42 |
| f4_eflux_score + f4_score_train | 0.02 |
| **Total** | **≈ 8.4** |

No Gurobi session was opened.

## Files

| Path | Content |
|---|---|
| `losses/line_losses.tsv` | per training line: signature, flags, lost / contradicted / dropped / accepted genes, closed reactions, sweep job |
| `losses/gene_counts.tsv` | per lost gene: lines lost, accepted, contradicted; reactions closed alone |
| `losses/parent_check.tsv`, `losses/summary.json` | P1-parent WT check; counts |
| `sweeps/ratio_r_delta_human1.parquet` | (job, gene) ratios that differ from P1; rebuilt by `f4_sweeps.job_ratios` |
| `sweeps/jobs.tsv`, `spot_checks.tsv`, `sweeps.json` | job status and cost; checks; summary |
| `train/` | F4a training cells, contrasts, binary, verification, summary (`per_line.tsv`: tAP per line, no label counts) |
| `eflux_panel/builds/` | F4b builds, gates, KO table, generic check, build summary |
| `eflux_panel/train/` | F4b pre-check contrasts, folds, config scores, verification, decision |
| `MANIFEST.tsv` | provenance rows (CPU-s in `cpu_s`) |
