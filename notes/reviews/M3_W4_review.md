# M3 / W4 independent review: context-specific model code before Stage A

**Reviewer:** independent review agent (not the author). **Date:** 2026-10-27.
**Code reviewed:** commit 7010b07 (`src/nextgem/context.py`, `context_inputs.py`, `context_run.py`, `gurobi_pool.py`;
`scripts/w4_prepare.py`, `w4_pilots.py`, `w4_stage_a.py`, `w4_summarise.py`; `tests/test_context.py`,
`tests/test_gurobi_pool.py`). HEAD during the review was 6bf19d2, which changes no W4 file.
**Spec:** DESIGN.md §9 W4, §13 (Stage-A table A1–A7), §11.4, §17 amendments A4 and A5;
`results/w4/README_ratified_A4.md`.
**Review runtime:** ≈ 1.9 h wall, ≈ 5 CPU-h, ≤ 8 processes, ≤ 2 Gurobi sessions held at once (pool).
**Recompute scripts and outputs:** `results/w4/review/` (scripts `r0`–`r7`, `own.py`, `with_slot.py`; JSON/TSV
outputs; regenerable intermediates in `results/w4/review/tmp/`, git-ignored).

## 0. Session probe, scope and label exposure
- **Probe.** `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label-guard hook
  ("Blocked by label guard: shell command references raw label files").
- **Label-free.** No `labels.*` call, no raw label file, nothing under `data/locked/`. Inputs used: DepMap
  expression, WGS CN, damaging-mutation matrix, ProCan (through the author's reader), `Model.csv` (through the
  author's ProCan reader), the canonical split file, and Human-GEM.
- **Incidental exposure (reported, not used).** While checking whether Stage A should use Human1 or Human2, I opened
  `results/w1/val_look_W1_M2_baseline_gstar.json`. It is the committed, already-revealed W1 validation-look record:
  G\* = Human1, with the validation Δ̄_lin(H2 − H1) and its CI. This is aggregate, validation-tier and already
  public in the commit log. It played no part in any check below. DESIGN §3.1 makes Human2 the default parent for
  optimisation whatever G\* is, so Stage A on Human2 is correct.
- **Gurobi.** I took sessions only through `nextgem.gurobi_pool` (n = 4), and held at most 2 at once. Every GLPK
  process set `cobra.Configuration().solver = "glpk"` before it loaded a model, so it opened no WLS session.
- **Side effect.** The requested `pytest` run (item 9) appended one row to the tracked `results/test_runs.tsv` and
  created `results/test_runs/cbf132dd74ef24e5.tsv`. This is the repo's automatic test-run record (a non-full run,
  23 passed). I left both in place and committed nothing; the PI may keep or discard them.

## Verdict summary

| # | Item | Verdict |
|---|---|---|
| 1 | Algorithm correctness: GIMME, E-Flux and FASTCORE rebuilt and compared | **PASS** (with NOTEs 1a–1d) |
| 2 | GPR evaluation and the A4 (a) KO convention | **PASS** |
| 3 | Repair rule and loss rule | **PASS** for the rules; their determinism falls under B2 |
| 4 | A5 grid | **PASS** |
| 5 | Leakage and tiers | **PASS** (NOTE 5a) |
| 6 | Solver discipline and the Gurobi pool | **FIX-REQUIRED** (B1), NOTEs 6a–6c |
| 7 | Reproducibility | **FIX-REQUIRED** (B2, B3), NOTEs 7a–7c |
| 8 | Compute projection and the A6 drop | **PASS** (NOTE 8a) |
| 9 | Tests | **PASS**: 23 passed (NOTE 9a) |

**GO / NO-GO: NO-GO until B1–B3 are fixed.** All three are small, mechanical changes. The algorithms are correct.

---

## 1. Algorithm correctness: PASS

**Method.** I wrote my own minimal GIMME and E-Flux (`results/w4/review/own.py`). They use cobra's GPR AST and
`GPR.eval`, plain cobra LPs, and an own brute-force KO sweep under the A4 (a) convention, with no zero-flux
shortcut. I rebuilt 3 panel lines (ACH-000957, ACH-002024, ACH-003383) × {F2, F3, F4} × R1 × real RNA × no genomic
constraint, each on its declared solver: F2 and F4 on GLPK, F3 on Gurobi. My own signed scores, from my own reader
(log2(TPM+1) − 1), equal the author's for all genes. Script: `r1_build.py`. Outputs: `r1_*.tsv`.

| Line | Family | Reaction set, own vs author | KO lethal (ratio_r < 0.01), author / own / shared | max \|Δ ratio_r\| over 2,870 genes |
|---|---|---|---|---|
| ACH-000957 | F2 GIMME | identical (11,269 open); same closed set; objective 1633.66476 in both | 168 / 168 / 168 | 0 |
| ACH-002024 | F2 | identical (11,028); objective 2408.19898 | 169 / 169 / 169 | 0 |
| ACH-003383 | F2 | identical (11,126); objective 2066.11128 | 167 / 167 / 167 | 0 |
| ACH-000957 | F4 E-Flux | bounds identical on all 7,629 bounded reactions (0 differ); cap_max identical | 154 / 154 / 154 | 0 |
| ACH-002024 | F4 | bounds identical (7,629) | 154 / 154 / 154 | 0 |
| ACH-003383 | F4 | bounds identical (7,629) | 154 / 154 / 154 | 0 |
| ACH-000957 | F3 FASTCORE | (author's support; own KO) 6,721 open | 173 / 173 / 173 | 0 |
| ACH-002024 | F3 | 6,423 (second build in process) / 6,429 (fresh process); see B2 | 182 / 182 / 182 | 0 |
| ACH-003383 | F3 | 6,587 | 176 / 176 / 176 | 0 |

- **GIMME** (`context.py:303-326`):
  - The required function is 0.9 × the generic WT biomass on the line's medium (`parent.wt_generic`, i.e. the
    maximum biomass), set as the biomass lower bound. This matches DESIGN and Becker & Palsson 2008.
  - The penalty is c_r = max(0, −s_r), the threshold minus the presence score, on fwd + rev variables.
  - After optimisation, a reaction is kept if s_r ≥ 0, if it has no evidence, or if |v| > 1e-6; the rest are closed.
  - My independent LP gives the same closed set and the same objective, to 1e-10.
  - Repair re-opened 0 reactions for F2 on all 3 lines, as expected: the GIMME solution already grows at 90%.
- **E-Flux** (`context.py:467-486`, `finish` `:677-685`, `closure_map` `:735-742`):
  - Bounds are max(10, 1000·cap/cap_max) on internal GPR reactions; boundary and GPR-less reactions keep the parent
    or medium bounds.
  - Reversible reactions get the symmetric bound (−b, b). Irreversible reactions keep their lb (Human2 has no
    internal GPR reaction with lb > 0 or ub < 0, so lb > ub cannot occur).
  - GPR aggregation is AND = min and OR = sum. An unknown gene is ignored under AND and makes an OR unknown, which
    keeps the default bound (README §5 item 5). Off genes have linear value 0.
  - A KO sets the gene's linear value to 0 and re-derives the bound. A reaction whose GPR becomes False is closed.
  - All of this is reproduced exactly by my implementation.
- **FASTCORE** (`context.py:329-422`):
  - The core is {s_r ≥ 0} ∩ fastcc-consistent; ε = 1 (cobra's default); the support tolerance is 1e-6. LP7/LP9 and
    the flip/singleton logic follow Vlassis 2014, Algorithm 1.
  - cobra `fastcc` (GLPK, flux_threshold 1) finds the author's output fully consistent:
    - ACH-000957 support: 6,630 / 6,630 consistent;
    - ACH-000957 final (after repair): 6,721 / 6,721;
    - further lines are listed in §10.
  - The support never leaves the cached consistent set (0 reactions outside it).
  - My brute-force KO on the author's F3 models equals the author's sweep exactly. This also validates the zero-flux
    shortcut and the pFBA support on real models.

**NOTE 1a. FASTCORE drops consistent core reactions, and the ratified README says it does not.**
- The ratified README §5 item 3 (`README_ratified_A4.md:198`) states that unreachable core reactions are "logged,
  not fatal. The pilots had 0."
- This is wrong. `pilot_builds_gurobi_lp.tsv` records `unreachable_core` = 13–40 per build, and my rebuilds give 28,
  25 and 23. Of these, 24, 21 and 19 consistent core reactions are absent from the support.
- FVA (`r6_core_fva.json`) shows that **all 24** missing core reactions of ACH-000957 have a maximum |flux| < 0.99
  on M0 (0.0015–0.86). Their subsystems are vitamin B12, vitamin E and retinol metabolism, plus transport and
  glycerophospholipid reactions.
- They are unreachable by construction at ε = 1, the ratified choice. This is a systematic exclusion of
  micronutrient-limited pathways (about 0.5% of the core), not a bug.
- The README statement should be corrected, and the effect stated in the W4 report. If the PI wants full core
  coverage, an ε below the smallest core flux capacity would be needed; that is an amendment, not a code fix.

**NOTE 1b. E-Flux binary KO profiles equal the generic model's** (`r7_generic_vs_f4.json`).
- On all 3 lines, F4's lethal set (154 genes) equals the generic Human2/M0 lethal set exactly. Only 62–64 genes
  have a partial ratio, against 62 in the generic model.
- With WT growth at about 1–1.6% of generic, the floor of 10 almost never binds after a KO. The binary F4 KO
  profile is therefore line-invariant, which follows directly from the ratified scaling and floor.
- Consequences:
  - in H4.3 (binary KO profiles), F4 will be "method-dominated" by construction;
  - F4's line-specific signal sits only in the partial ratios.
- Under the literal reading of A4 (a) ("a KO zeroes the gene capacity", with the floor still applied and no GPR
  closure), none of the 154 lethal genes would be lethal. So the author's reading (close a GPR-False reaction) is
  the only one that gives E-Flux any lethal calls. I agree with that reading, but it should be stated explicitly.
- This is not a code defect. The PI should know it before Stage A is interpreted.

**NOTE 1c. FASTCORE support tolerance.** The support is |v| > 1e-6 (README §5 item 3, ratified). Vlassis' code
uses |v| ≥ 0.99·ε, so this support is more inclusive. Not a defect.

**NOTE 1d. FASTCORE closes every reaction outside the support, including exchanges.** F3 lines therefore import
only the medium components that FASTCORE's support (plus repair) used. This is standard FASTCORE behaviour, but it
makes F3's effective medium line-specific. State it in the W4 report.

## 2. GPR evaluation and the A4 (a) KO convention: PASS
Script: `r2_gpr_toys.py`; output: `r2_gpr_toys.json`.
- **Random nested GPRs.** 5,000 random nested AND/OR GPRs (depth ≤ 3; 25% of genes unknown; scores N(0,1)):
  - the author's `presence_score` equals my independent evaluator in every case;
  - `s_r < 0` ⇔ cobra `GPR.eval(knockouts = off)` is False, with 0 violations.
- **Hand-made toy** (isozymes, nested (A∧B)∨(C∧D), a complex with a bypass, an AND with a nested OR holding an
  unknown gene), through the author's `closure_map` and `ko_sweep` against my brute force. All match:
  - **Off isozyme does not buffer:**
    - gI1 ∨ gI2 with gI2 off: KO gI1 → ratio 0;
    - the generic parent is buffered: ratio 1.
  - **Nested:** with gC off, the pruned GPR is gA ∧ gB. KO gA or gB closes the reaction; KO gD has no effect.
  - **Kept or re-opened complex with a False pruned GPR keeps its full GPR:** KO of either subunit, including the
    off one, closes it. This is also checked with a re-opened reaction under a genomic loss.
  - **Unknown genes:** gU ∧ (gV ∨ gW) with gV off and gW unknown stays open, and KO gW closes it.
  - **E-Flux:**
    - toy bounds as expected (1000, 500, 10, 500, ±500; cap_max 4);
    - KO of the expressed isozyme drops the bound to the floor (10), because the off isozyme keeps the GPR True;
    - KO of a single-gene reaction closes it;
    - the author's sweep equals my brute force.
- The code path is `context.py:719-743` (`closure_map`); `presence_score` is at `:108-123`.

## 3. Repair rule and loss rule: PASS (rules); the determinism caveat goes to B2
Script: `r3_loss_repair.py` (GLPK); output: `r3_loss_repair.tsv`.

**Loss rule** (`context_inputs.py:216-242`, `context.py:510-542`):
- My own r2-11 rule, on my own reads of the CN and damaging matrices, gives the same lost-gene lists as the author's.
- The author's drop gives identical results under 6 input orders.
- My independent drop (cobra GPR closures, GLPK growth, 1% floor, descending-expression / ascending-Entrez order)
  gives the same accepted and contradicted sets.
- Both equal `results/w4/tmp/accepted_losses_panel.json` and the smoke runs' `accepted_losses.json`.
- Per-line numbers are in §10.
- Order independence holds by construction: the input is canonicalised with `sorted(set(...))` (`context.py:518`).

**Repair** (`context.py:567-615`, `finish` `:686-696`), F0, R1, on the none and loss variants:
- WT before repair ≈ 0; after repair 65.45 on ACH-000957 (0.876 of generic). Restored on every line.
- 35 reactions re-opened (0.27% of reactions); the size is reported in `n_reopened` / `repair_frac_reactions`.
- **Repeating `finish` on the same extraction:**
  - identical re-opened set on ACH-000957 (none and loss variants);
  - **not identical on ACH-002024 and ACH-003383**: same count (35 and 39), different reactions. These repeats ran
    after the loss-rule LPs in the same GLPK process;
  - in a clean fresh process (`r8_glpk_history_ACH-002024.json`), two consecutive repeats were identical (symdiff 0;
    KO tables identical).
  So the GLPK repair LP's alternate optimum depends on the process history (see B2). The rule itself is
  implemented as specified: union of the growth L1 LP and the 10 task L1 LPs, weights 1 + deficit, 1% floor.
- Re-opened ∩ accepted-loss closures = ∅. 14 loss closures on ACH-000957 stay closed.
- My own growth L1 LP support (33 reactions) is a subset of the author's re-opened union (growth + 10 task LPs).
- All 10 viability tasks (GR + 9 ER) pass on the repaired line model, as on the generic parent.

## 4. A5 grid (A4 (d)): PASS
Script: `r0_inputs.py`; outputs: `r0_inputs.json`, `a5_dedupe_check.tsv`.
- **Grid.** `LCP_GRID` (`context_inputs.py:309`) is exactly the 24 points w_r = 1, w_p ∈ {0, 0.5, 1},
  w_c ∈ {0, 1}, τ ∈ {−1.5, −1, −0.5, 0}.
- **z statistics are training-only.** I recomputed gene mean and SD (ddof 1) over the 714 training lines with
  expression, and the ProCan floor-filled mean and SD over the 363 training ProCan lines. Both equal
  `threshold_gene_means.parquet` exactly (max |Δ| = 0). Including validation and test lines would shift gene means
  by up to 0.62, so the check is sensitive.
- **Scores.** loss = −1 for accepted losses; z_prot = 0 for lines without ProCan and for non-soluble or
  non-quantified genes. My independent z-scores reproduce the author's `lcp_scores` on all 720 (line, point)
  pairs, to 1e-9.
- **Deduplication** (`w4_stage_a.py:296-311`):
  - the 720 pairs map to **400** distinct score vectors (10 ProCan lines × 24 + 20 lines × 8), the README figure;
  - every line has ≥ 1 accepted loss, so w_c is never redundant;
  - the mapping is by identical score vector, which is conservative: it refines the 352 distinct off-sets, and
    builds with the same off-set but different scores can differ through the repair weights;
  - `a5_grid_map.tsv` gets a row for every (line, point).
- **Solver.** F5 is declared GLPK (`context.py:72`).
- **NOTE.** For F5 the genomic loss is a soft term. With w_c = 1 and τ = −1.5, a lost gene with z_RNA ≥ −0.5
  stays on, and with w_c = 0 accepted losses are ignored. This follows A4 (d) exactly; I note it only so the A5
  results are not read as "with genomic constraint".

## 5. Leakage and tiers: PASS
- **Panel.** My own implementation of README §5 item 8 reproduces `stage_a_panel.tsv` exactly:
  - eligible lineages = the 12 training lineages with ≥ 15 lines, sorted;
  - `default_rng(11)` permutes them and the first 10 are taken;
  - 3 lines are drawn per lineage from the sorted expression-profiled training lines.
  All 30 lines are training-tier with expression; n_p = 10 and all 30 have CN. Membership is also asserted at run
  time (`w4_stage_a.py:237-238`, `context_inputs.py:282`).
- **Thresholds.** R1 is a constant. R2 (P25 = 0.36952) and R3 (gene means, clipped to [P25, P75]) are recomputed
  from the 714 training lines, and the protein statistics from the 363 training ProCan lines. They match exactly.
- **Inputs.** `w4_prepare.py:41-55` uses `tier == "train"` only. Stage A reads omics for panel lines only
  (`w4_stage_a.py:165-170`).
- **Labels.** No W4 module or script imports `labels`, and no label is used in any Stage-A step.

**NOTE 5a. The repair pilot used validation and test omics.** `repair12` used the expression of 1 test and 2
validation lines (the compute probe's first 12 rows). It is a label-free feasibility pilot, not a fitting step, and
the IDs are masked in `repair_probe12.tsv`. The masks (`masked_test_1`, …) keep the row position, though, so the IDs
are re-derivable. The content is label-free model properties, so this is harmless; avoid it in future pilots.

## 6. Solver discipline (A4 (c)): FIX-REQUIRED

**What passes.**
- `DECLARED_SOLVER` (`context.py:72`) is GLPK for F0, F2, F4 and F5, and Gurobi for F3 and F3b. The loss rule
  is GLPK (`w4_stage_a.py:68`).
- Each worker re-checks per job (`w4_stage_a.py:127-130`). `--solver-override` requires `--pilot`, and `--pilot`
  requires `--out`.

**B1 (blocking). GLPK workers and the parent process open unpooled Gurobi WLS sessions.**
- `context_run.init_worker` calls `load_model()` with cobra's default solver, Gurobi, and only then switches to
  GLPK (`context_run.py:47-49`). Building the cobra model calls `gurobipy.Model()`, which starts the
  process-default WLS environment.
- That environment stays alive for the life of the process, outside the slot pool. The same happens in the parent
  (`Inputs.__init__`, `w4_stage_a.py:156`).
- Evidence:
  - `results/w4/logs/glpk.log` shows 6 "Academic license … WLS" banners for a run with 5 GLPK workers plus the
    parent;
  - `logs/prepare_a5.log` shows one for a script that needs no solver.
- With the production command (20 GLPK workers), Stage A would hold **≈ 21 unpooled WLS sessions plus 4 pooled
  ones**, while W2 also uses Gurobi. This breaks §11.4 / A4 (c) and risks licence-limit errors, which per §11.4
  abort jobs.
- **Fix:**
  - set `cobra.Configuration().solver = solver` before `load_model` in `init_worker`, and `"glpk"` in
    `Inputs.__init__` and `w4_prepare.py`, as `w1_inputs.py` and `w2_runtime.py` already do;
  - add a test that a GLPK worker init never creates a gurobipy model.
- I verified the fix pattern: none of my GLPK review processes printed a WLS banner.

**NOTE 6a.** `--pilot` alone disables the per-job solver assertion (`w4_stage_a.py:129`), even without
`--solver-override`. Better to assert `solver_of[family]` always, with overrides applied, and to compare
`solver_name(parent.model)`, not the `_W["solver"]` string.

**NOTE 6b. Pool semantics** (`r5_pool_crash.py`, `r5_pool.json`):
- A SIGKILLed holder frees its slot. The lock files left in /dev/shm are not stale locks, because flock state dies
  with the process.
- **However,** the pool size is per caller: `hold_slot(n)` with n = 6 takes a 5th slot while 4 are held.
  `--gurobi-slots` (`w4_stage_a.py:205`) lets a run exceed the licence pool.
- Fix the size globally (`NEXTGEM_GUROBI_SLOTS`, default 4) and refuse a CLI value above it.
- `gurobi_pool` is used only by W4 so far; grep finds no W1 or W2 use. The pool limits nothing unless every
  workstream adopts it.

**NOTE 6c. The pool test checks semantics only.** `tests/test_gurobi_pool.py` covers limits, release,
cross-process blocking and the env override in a temp dir. It does not cover crash release, nor that W4 workers
actually hold a slot, nor B1.

## 7. Reproducibility: FIX-REQUIRED

**What passes.**
- **Seeds.** The panel and the permutation use SEEDS[0] = 11. A7 uses SEEDS[0:3]: GIMME 1e-3·U(0,1) cost noise,
  and FASTCORE N(0, 0.25) log2 noise (README §5 item 10). Gurobi runs with Threads = 1 and Seed = 11. No unseeded
  RNG.
- **Dirty-tree refusal** (`w4_stage_a.py:144-149, 252-254`). On a throwaway repo, untracked or modified files under
  `scripts/` give dirty = True; changes only under `results/` do not. The script exits before any work.
- **Partly history-stable on GLPK.** F2 built fresh and F2 built after F0 in the same process have identical masks
  and KO hashes. The F0 repair, however, is not always stable (§3, B2).

**B2 (blocking). Builds depend on the worker's history (F3 on Gurobi; the F0 repair on GLPK), and GLPK LPs can stall.**
- ACH-002024 F3 with identical inputs, code and solver gave three different results:
  - fresh process: 6,429 open, mask 694001…, reproduced bitwise by a second fresh process;
  - after another F3 build in the same process: 6,423 open, mask 991b13…;
  - the pilot worker: 6,432 open.
  Reaction-set Jaccard between these is ≈ 0.985, and the KO tables differ.
- The cause is Gurobi warm-starting from the previous basis and model state, combined with FASTCORE's
  alternate optima.
- **GLPK repair is affected too** (§3): two consecutive F0 `finish()` calls after other LPs in the same GLPK process
  re-opened different reaction sets of the same size, on 2 of 3 lines. In a clean process the repeats were
  identical. The repair weights tie on many reactions (w = 1 + deficit), so alternate L1 optima are common.
- **GLPK liveness.** In two attempts of `r8_glpk_history.py`, the third consecutive F0 finish + sweep in one GLPK
  process ran for more than 8 and more than 30 CPU-min, against 56–59 s for the first two, and I killed it.
  Nothing bounds a GLPK solve:
  - no GLPK time limit;
  - `guard()` (`w4_stage_a.py:270-275, 348`) runs only after a job completes;
  - so a stalled LP would hang its worker and the run.
- ProcessPoolExecutor assigns jobs to workers dynamically (`w4_stage_a.py:265`), so Stage-A F3 results, including
  the A7 FASTCORE stability row, would depend on scheduling. A rerun or resume would not reproduce them, and
  scheduling noise would enter the method-vs-tissue variance partition.
- **Fix, either of:**
  - (i) `max_tasks_per_child=1` for the Gurobi pool, i.e. one build per fresh process. It costs about 20 s of
    model load per F3 job, ≈ 1.5 CPU-h over the 270 F3 jobs (A1 120, A2 30, A3 60, A7 60; the A4 re-sweeps ride
    on the same jobs);
  - (ii) reset the solver state at the start of every build, e.g. `model.solver.problem.reset(1)` for Gurobi, and
    verify.
- Apply the chosen fix to both pools. For the GLPK pool, option (i) costs ≈ 20 s × ≈ 1,090 GLPK jobs ≈ 6 CPU-h.
  Option (ii) needs a basis reset for GLPK, e.g. `glp_std_basis`.
- Add a GLPK time limit (optlang `configuration.timeout`) and treat a timeout like a numerical non-optimal status
  (§11.4: retry, then failed solve). Add a per-job wall-clock watchdog to the parent loop.
- Either way, add a test: the same F3 build, and the same F0 build, run fresh and after a different build give a
  bitwise-identical open mask and KO table.

**B3 (blocking). Provenance (§11.4).**
- `w4_stage_a.py` does not run under `nextgem.run.provenance()`. The function exists: it was committed in 67f5f4f,
  an ancestor of 7010b07. The README's "does not exist yet" (`README.md:302`) is stale.
- The hand-rolled MANIFEST row (`w4_stage_a.py:367-377`) lacks input SHAs, output SHAs, solver and package
  versions, and the author. The run inputs that are results files would not be pinned: `stage_a_panel.tsv`,
  `threshold_stats.json`, `threshold_gene_means.parquet` and `fastcc_human2_M0.json`. The dirty check ignores
  `results/`.
- §11.4 requires every result-producing script to run under provenance(). Stage A (≈ 60 CPU-h) is also a
  pinned-worktree job. In a worktree, `.venv/bin/python` imports the editable install from the main checkout unless
  guarded, and provenance's `nextgem.__file__` assertion is what catches that.
- **Fix:** wrap `main()` in `provenance(...)` with these inputs plus the omics files, and register the outputs.

**NOTE 7a. Checkpoints.** Checkpoints are keyed by build_id only (`w4_stage_a.py:187-188, 314`), and
`accepted_losses.json` is reused if present (`:279-280`). A resume after a code change would silently mix
versions. Record the commit in TMP and refuse a resume across commits.

**NOTE 7b. ratio_r is stored as float32** (`context_run.py:105`). Pilot-vs-rerun ratio differences of ≤ 3e-8 are
float32 rounding, not model differences. Prefer float64, or document that ratio_r is float32 (6-decimal rounding
survives float32 near 1, so this is cosmetic).

**NOTE 7c. A7 perturbations.**
- A7 uses `lines[:20]` (`w4_stage_a.py:108`), the first 20 panel rows. That is 7 lineages; Head and Neck, Uterus,
  Cervix and part of Esophagus/Stomach are excluded. A stratified 2-per-lineage choice would be more defensible.
- Each seed applies the identical noise vector to every line (`:180`, `gimme_closed` `:310`), so perturbations
  are correlated across lines. Acceptable, but state it.

## 8. Compute projection and A6: PASS
I recomputed the arithmetic from `pilot_summary.json`.
- **Per-family Stage-A CPU** (A1–A4, plus A7 for F2 and F3): F0 10.7, F2 14.6, F3 14.2, F4 19.7 CPU-h.
- **GLPK lane** at the measured GLPK/Gurobi ratios (0.87, 0.98, 0.57) plus A5 (400 builds × 1.72 min = 11.5 CPU-h):
  46.3 CPU-h. **Gurobi lane** (F3): 14.2 CPU-h. **Total ≈ 60.5 CPU-h** ✓.
- **Wall time:** max(14.2/4, 46.3/20) = 3.55 h ✓.
- **Memory:** 20 × 1.56 + 4 × 1.43 = 36.9 GB ✓.
- **My own builds are consistent.** Author build CPU without tasks: F2 82–93 s, F4 75–86 s, F0 83 s on GLPK;
  F3 102–139 s on Gurobi. The pilot's like-for-like figures are 85, 93, 91 and 160 s.
- **A6 drop.** The CORDA pilot cost 1.74 and 1.52 CPU-h per model. Both are above 1 CPU-h, so A6 is dropped under
  its frozen rule, as A4 (e) records ✓.

**NOTE 8a. The wall projection is optimistic.**
- It assumes all 4 Gurobi sessions go to W4. W2 is running concurrently, and with 2 sessions the Gurobi lane takes
  ≈ 7 h.
- The fix for B2 (i) adds ≈ 1.5 CPU-h to the Gurobi lane.
- The parent process (≈ 3.5 GB) is not in the 37 GB.
- Everything remains well inside the 100 CPU-h row budget.

## 9. Tests: PASS
`.venv/bin/python -m pytest -q tests/test_context.py tests/test_gurobi_pool.py`, run while holding a pool slot:
**23 passed in 4.3 s.**

**NOTE 9a. Coverage gaps.** The suite has no test for B1 (no Gurobi session in GLPK workers), none for B2 (history
independence), and none that `--pilot` is the only override path (the existing test checks only the refusal).

## Other documentation NOTEs
- `README_ratified_A4.md:151` and `README.md:213` contain a garbled paragraph: "WALL_" followed by a duplicate of
  §3.5. A4 pins this file by hash, so leave it, but correct the working README.
- `README_ratified_A4.md:198`: "The pilots had 0" unreachable core reactions is false (NOTE 1a).
- A4 pin check: `sha256(results/w4/README_ratified_A4.md)` = 16f6aac0…34a4, which matches the pin "16f6aac0…a4" ✓.

## 10. Per-line detail (loss rule, repair, fastcc)

| Line | Lost model genes (own = author) | Contradicted / joint-dropped | Accepted (= pilot file) | 6 orders identical | F0 WT after repair | Re-opened | Loss closures (none re-opened) | Viability tasks (line / generic) | fastcc on F3 output |
|---|---|---|---|---|---|---|---|---|---|
| ACH-000957 | 23 | 0 / 0 | 23 ✓ | ✓ | 65.45 | 35 | 14 | 10 / 10 | support 6,630/6,630; final 6,721/6,721 consistent |
| ACH-002024 | 4 | 0 / 0 | 4 ✓ | ✓ | 46.74 | 35 | 2 | 10 / 10 | support 6,331/6,331 consistent |
| ACH-003383 | 1 | 0 / 0 | 1 ✓ | ✓ | 46.76 | 39 | 0 | 10 / 10 | final 6,587/6,587 consistent |

Joint drop (step 4) is not exercised on these lines, nor on any panel line (loss-rule pilot); the toy test covers it.

## Blocking fixes (in order)
1. **B1:** set the cobra solver before loading a model in every GLPK process (`context_run.init_worker`,
   `w4_stage_a.Inputs`, `w4_prepare`), so that no unpooled WLS session is opened, and add a test.
2. **B2:** make builds independent of worker history in both pools (F3 on Gurobi; the F0/F2/F5 repair on GLPK):
   one build per process (`max_tasks_per_child=1`), or a basis reset per build. Add a GLPK time limit and a per-job
   watchdog. Add a test that the same build run fresh and after another build is bitwise-identical.
3. **B3:** run `w4_stage_a.py` under `nextgem.run.provenance()`, with the results-file inputs (panel, threshold
   statistics, FASTCC cache) and the omics files hashed, outputs registered, and the worktree assertion in force.
   Update README §6 ("provenance does not exist yet").

**Recommended before launch (non-blocking):**
- NOTE 6a: always assert the declared solver.
- NOTE 6b: make the pool size global.
- NOTE 7a: commit-stamped checkpoints.
- NOTE 1a: correct the README text on `unreachable_core`.
- NOTE 1b: tell the PI that F4's binary KO profile equals the generic model's.

After B1–B3 are fixed and re-tested (pytest plus one smoke run per solver lane checking: no WLS banner in GLPK
workers, an identical F3 mask across job orders, a MANIFEST row with input SHAs), Stage A can be launched:
**GO conditional on B1–B3.**
