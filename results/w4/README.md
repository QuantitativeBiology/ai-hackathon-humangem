# W4: context-specific model extraction, M3 pilots and Stage-A preparation

**Author:** W4 engineer agent (label-free; never opened a label file, nor `labels.py`). **Date:** 2026-10-26.
**Status:** pilots done; PI decisions recorded as **amendment A4** (DESIGN.md §17) and implemented (§0).
**The Stage-A sweep has not been run.** The code is uncommitted, so the production script refuses to start
until it is committed (§11.4). The ratified version of this README (sha256 16f6aac0…a4, the hash A4 pins)
is kept byte-identical as `README_ratified_A4.md`; this file has since changed (§0 added, projections updated).

## 0. Amendment A4 (PI decisions, 2026-10-27) and their implementation

| A4 | Decision | Implementation |
|---|---|---|
| (a) | KO gene-state convention ratified | unchanged (`context.closure_map`) |
| (b) | §5 items 2–12 ratified as written | unchanged; ratified text in `README_ratified_A4.md` |
| (c) | one declared solver per family: GLPK for F0, F2, F4; Gurobi for F3 (and F3b), via a flock slot pool | `context.DECLARED_SOLVER`. **F5 is not named in A4**: it is declared GLPK because it uses the F0 machinery (pruning + repair), which passed GLPK = Gurobi; PI to confirm. The loss rule (phase 0) is `context.LOSS_RULE_SOLVER` = GLPK. Every job asserts that its parent model runs the declared solver, always (also under `--pilot`); only `--solver-override` (which needs `--pilot` and `--out`) changes the declaration. Pool: the shared module `src/nextgem/gurobi_pool.py`, with one global size (`POOL_SLOTS` = 8, or `NEXTGEM_GUROBI_SLOTS`), no per-caller size; `--gurobi-workers` above it is refused |
| (d) | A5: fixed 24-point grid: w_r = 1; w_p ∈ {0, 0.5, 1}; w_c ∈ {0, 1}; τ ∈ {−1.5, −1, −0.5, 0}; z_prot = 0 without proteomics; identical builds deduplicated | `context_inputs.LCP_GRID`, `lcp_scores`; row A5 in `w4_stage_a.py`; `a5_grid_map.tsv` maps all 720 (line, point) pairs to their build |
| (e) | A6 dropped and reported | superseded by **A7** (user decision D8): CORDA runs as row A6′, see §0.2 |
| (f) | Stage A on M0 | unchanged |

**A5 score definition (implemented; for the reviewer).**
- A gene is off iff `w_r·z_RNA + w_p·z_prot + w_c·loss < τ`. Internally the signed score is
  `s_g = w_r·z_RNA + w_p·z_prot + w_c·loss − τ`, and the F0 machinery (pruning, repair, KO convention)
  runs on it.
- **z_RNA** = (log2(TPM+1) − gene mean) / gene SD, both over the 714 training lines with expression
  (the same fitting set as R3; ddof 1). A gene with no column or zero SD gets z = 0.
- **z_prot** = (log2 ProCan − gene mean) / gene SD over the 363 training lines with ProCan. "Not detected"
  is set to the protein floor value (training minimum − 1, README §5 item 7) both in the line and in the
  training distribution.
  - It is used only for soluble enzymes quantified in ProCan; other genes get z_prot = 0.
  - z_prot = 0 for every gene of a line without proteomics.
- **loss** = −1 for genes in the line's accepted-loss set (loss rule after the order-independent drop),
  else 0. The sign makes the term push lost genes towards "off"; the frozen formula leaves it implicit.
- The statistics are in `threshold_gene_means.parquet` (`rna_mean`, `rna_sd`, `protein_mean_floor`,
  `protein_sd_floor`), written by `w4_prepare.py`. Regenerating them left `stage_a_panel.tsv` and
  `threshold_stats.json` byte-identical.
- A5 builds have no A4 genomic re-sweep: loss enters through w_c.
- **Deduplication:** grid points whose score vectors are identical on a line (to 1e-9) are built once.
  This happens for w_p on lines without proteomics, and for w_c on lines with no accepted loss.

### 0.1 Review M3_W4 (NO-GO) fixes, 2026-10-27
- **B1 (no unpooled WLS sessions).**
  - `context_run.init_worker` sets `cobra.Configuration().solver` before loading the model.
  - The Stage-A parent (`Inputs`, `a6p_subset`), `w4_prepare.py` and `w4_pilots.py` read Human-GEM under GLPK.
  - FASTCC in `w4_prepare.py` and every Gurobi pilot worker now hold a pool slot.
  - Test `test_glpk_worker_never_touches_gurobi`: a subprocess replaces `gurobipy.Model` and `gurobipy.Env`
    with functions that abort, then runs a full GLPK worker init and a parent-style SBML read. Both succeed,
    and no WLS banner is printed.
- **B2 (history independence, liveness).**
  - **One build per fresh spawned process** (`context_run.Runner`; the loss-rule jobs too), so no build sees
    another's solver state.
  - Every LP has a 300 s limit (`context.LP_TIME_LIMIT_S`: Gurobi `TimeLimit`, GLPK optlang timeout). A
    non-optimal solve is retried once (Gurobi NumericFocus 2; GLPK from the standard basis), then counted as
    a failed solve (§11.4).
  - **Per-job watchdog** in the parent. It starts once the child's solver is ready, so slot waiting and model
    load are excluded. It kills a job above its CPU limit (`--job-timeout`, default 1,800 s; CORDA 5 h), or
    above twice that in wall time. The job is recorded as `failed-timeout` and the run continues.
  - Any other worker exception aborts the run with no consolidated output.
  - Test `test_build_independent_of_job_order`: the same build run first, and run after another build in the
    same lane, gives a bitwise-identical open mask and KO table, for F0 on GLPK (repair) and F3 on Gurobi
    (FASTCORE).
- **B3 (provenance).**
  - `w4_stage_a.main()` runs under `nextgem.run.provenance()`. The manifest row
    (`stage_a/MANIFEST.tsv`) hashes every input and output.
  - Inputs hashed: the omics files, Human-GEM, the task list, the split, the panel, the threshold statistics,
    the FASTCC cache and `a6prime_lines.tsv`.
  - Outputs registered: `builds.tsv`, `ko_*.parquet`, `accepted_losses.json`, `a5_grid_map.tsv`,
    `a6prime_lines.tsv`.
  - It refuses a dirty tree unless `--pilot`, and refuses to run if the imported `nextgem` is not the package
    next to the script (pinned-worktree guard). The hand-rolled MANIFEST was removed.
- **Non-blocking items, all done:**
  - checkpoints live in `tmp/jobs_<sha of the W4 code files + run parameters>/`, so a resume after a code
    change starts afresh;
  - KO `ratio_r` is stored as float64;
  - `--pilot` alone keeps the solver assertion;
  - the pool has one global size;
  - A7 uses the first 2 panel lines of each of the 10 panel lineages, with noise drawn independently per line
    (seed sequence `[seed, panel position]`).

### 0.2 Amendment A7 (user decision D8): CORDA row A6′
- **Builds:** CORDA × {R1, R3} × real RNA on 10 representative panel lines, 20 builds, on Gurobi through the pool.
- **Line choice (label-free), `a6prime_lines.tsv`:**
  - one line per panel lineage: the panel line with the smallest Euclidean distance (over the 2,854 model
    genes with an expression column) to the lineage's median log2(TPM+1) profile;
  - the median is taken over **all training lines of that lineage with expression** (18–89 lines), not only
    its 3 panel lines;
  - ties are broken by ModelID; there is a test (`test_a6prime_line_choice`).
  - Chosen lines: Bone ACH-001283, Bowel ACH-003383, Cervix ACH-002020, Esophagus/Stomach ACH-000736,
    Head and Neck ACH-001692, Kidney ACH-002528, Myeloid ACH-000113, PNS ACH-000345, Skin ACH-000899,
    Uterus ACH-000407.
- **CPU limit:** 5 CPU-h per build (3× the pilot cost); above it the build is `failed-timeout`.
- **Tagging:** CORDA is tagged "descriptive, not eligible for P1 (amendment A7)" in `builds.tsv`
  (`status_tag`, `p1_eligible = False`), `ko_A6p.parquet`, `a6prime_lines.tsv` and the comparison.
- **Label-free comparison** (`w4_summarise.py --stage-a <dir>` writes `<dir>/a6prime_comparison.tsv`): each
  CORDA build against F0/F2/F3/F4 from A1 on the same line and threshold, reporting:
  - reaction-set Jaccard;
  - essential-set (ratio_r < 0.01) Jaccard and set sizes;
  - Spearman correlation of ratio_r;
  - WT growth ratio, re-opened reactions and build CPU.

  The label-based tAP comparison is not done here; it comes later through `labels.train()/fold()`.

### 0.3 Smoke test of every row (2 panel lines, 1 CORDA build; `stage_a_smoke/`, `--pilot`, dirty tree)
- **Setup:** first 2 panel lines (ACH-000957, ACH-002024; both Bowel), every row (A1, A2, A3, A4, A5, A6′, A7),
  1 CORDA build (`--max-corda 1`), 2 Gurobi + 5 GLPK workers, one build per process.
- **Result:** **109 builds + 48 A4 re-sweeps; 0 failed, 0 timeouts.**
  - Wall 5,840 s, dominated by the single CORDA build (95.6 CPU-min); 7.7 CPU-h in all (provenance, children
    included). Without CORDA the run took about 3,400 s.
  - A5: 48 grid points gave 32 distinct builds.
  - A6′ chose ACH-000957 (the Bowel line closest to the lineage median of its 2-line subset).
- **Solvers:**
  - Every KO row carries its declared solver (F3 and F3b on Gurobi; F0, F2, F4, F5 and the loss rule on GLPK).
  - **The log has exactly 23 WLS banners = 22 F3 processes + 1 CORDA process**, so no GLPK process (85 builds
    + 2 loss-rule jobs) and not the parent opened a Gurobi session.
- **Outputs:**
  - KO ratios are float64.
  - `status_tag` = "descriptive, not eligible for P1 (amendment A7)" on CORDA; "label-calibrated (F5)" on A5.
  - `MANIFEST.tsv` (provenance) row: commit, dirty = True (pilot), 13 input SHA-256s (omics, Human-GEM, task
    list, split, panel, threshold statistics, FASTCC cache), 11 output SHA-256s, solver map, versions, wall,
    CPU and RSS.
- **Repair:**
  - max re-opened fraction 1.6% (F0), 0.9% (F3), 2.4% (F5); all below the 5% gate.
  - FASTCORE `unreachable_core` 28–34 per build (see §5 item 3).
- **A6′ comparison** (label-free, `stage_a_smoke/a6prime_comparison.tsv`), ACH-000957, R1, CORDA against
  F0 / F2 / F3 / F4:

| Measure | F0 | F2 | F3 | F4 |
|---|---|---|---|---|
| reaction-set Jaccard | 0.59 | 0.59 | 0.83 | 0.56 |
| essential-set Jaccard | 0.97 | 0.97 | 0.98 | 0.90 |
| ratio_r Spearman | 0.96 | 0.90 | 0.92 | 0.89 |
| essential genes | 172 vs 168 | 172 vs 168 | 172 vs 173 | 172 vs 154 |
| WT / generic | 0.870 vs 0.876 | 0.870 vs 0.999 | 0.870 vs 0.968 | 0.870 vs 0.013 |
| build CPU-min | 95.6 vs 2.1 | 95.6 vs 1.8 | 95.6 vs 2.6 | 95.6 vs 2.0 |

  CORDA keeps 7,126 reactions open and re-opens 0.

### 0.4 Projection with A5 and A6′ (declared solvers, Gurobi pool = 8)
From `pilot_summary.json → stage_a_final` (pilot unit costs; per-job process overhead 35.6 CPU-s, measured
in the smoke test):

| Lane | Jobs | CPU-h |
|---|---|---|
| Gurobi: F3 (A1–A4, A7) + A6′ CORDA (20 × 1.63 CPU-h = 32.6) + overhead | 290 | 49.7 |
| GLPK: F0, F2, F4 (A1–A4) + F2 A7 + A5 (400 builds, 11.1 CPU-h) + loss rule + overhead | 1,120 | 57.0 |
| **Total** | 1,410 | **≈ 107 CPU-h** |

- **Wall ≈ 6.2 h** with 8 Gurobi + 20 GLPK workers. The Gurobi lane gives 49.7 / 8 h; the 20 CORDA builds
  alone need ≥ 3 rounds × 1.63 h ≈ 4.9 h; the GLPK lane gives 57 / 20 ≈ 2.9 h.
- **Memory** ≈ 47.5 GB: 8 × 1.6 + 20 × 1.56 + 3.5 GB for the parent. Launch only with ≥ 68 GB available,
  to keep 20 GB free.
- **Range:** in the smoke test the F0 and F4 GLPK builds cost up to ≈ 30% more per build than in the pilot,
  so an upper estimate is ≈ 120 CPU-h.
- **Row budget:** the total is slightly above the ≈ 100 CPU-h Stage-A row because A6′ adds 33 CPU-h
  (A7 amendment) and the per-process isolation adds ≈ 14 CPU-h (review B2).
- **Budget flag:** `--budget-cpu-h 110` aborts at 220 h of build CPU. The guard counts in-job build CPU only,
  not the process overhead.

**Production command** (after the commit and re-review, from the pinned worktree):
```
.venv/bin/python scripts/w4_stage_a.py --gurobi-workers 8 --glpk-workers 20 --budget-cpu-h 110
```

## 1. What was built

| File | Content |
|---|---|
| `src/nextgem/context.py` | GPR machinery (presence min/max and E-Flux min/sum scores), families F0, F2 GIMME, F3 FASTCORE, F3b CORDA, F4 E-Flux and the F5 builder; genomic constraint with the order-independent drop; the fixed repair rule; line-model gene-state convention; KO sweep with exact zero-flux shortcut; 57-task capability profile |
| `src/nextgem/context_inputs.py` | label-free inputs: RNA, WGS CN, damaging mutations, ProCan (via UniProt); R1/R2/R3 thresholds and signed scores for the RNA, protein and RNA+protein modalities; loss rule (r2-11, v1 and CN sensitivities); Stage-A panel draw; permuted-expression control |
| `src/nextgem/context_run.py` | worker runner: extract → (genomic, A4 post-extraction) → repair → sweep → tasks → open-reaction bitmask |
| `scripts/w4_prepare.py` | panel, threshold statistics, FASTCC cache |
| `scripts/w4_pilots.py` | `repair12`, `lossrule`, `cost`, `glpk` pilots |
| `scripts/w4_stage_a.py` | **ready-to-run Stage-A sweep** (§4; A5 and declared solvers per §0) |
| `src/nextgem/gurobi_pool.py` | shared Gurobi flock slot pool (A4 c), for every workstream |
| `scripts/w4_summarise.py` | writes `pilot_summary.json` |
| `tests/test_context.py` | 18 toy and known-answer tests (§12 row W4); each solver test runs on GLPK and on Gurobi |

**Reuse (§2).** cobra `fastcc` and its LP7 helpers (FASTCORE), raven `apply_task_constraints` (repair
task LPs, task profile), `corda` 0.5.1 (F3b, unmodified apart from a NumPy shim for the removed
`np.in1d`). **troppo does not install:** troppo 0.0.7 needs `cobamp==0.2.1`, which pins `pandas==1.3.5`;
that version does not build on Python 3.11 and would downgrade the pinned pandas. So F2 GIMME is our
own implementation, as DESIGN allows ("otherwise own"), and it passes its toy test. `corda==0.5.1` was
installed and `requirements.txt` regenerated (one added line).

## 2. Toy tests (§12 W4): all pass

`NEXTGEM_NO_TEST_RECORD=1 .venv/bin/python -m pytest tests/test_context.py tests/test_gurobi_pool.py` gives
26 passed (23 + 3): the A5 score known-answer test; the refusal of an undeclared solver (also under
`--pilot`); B1, no gurobipy in a GLPK worker; B2, job-order independence for F0/GLPK and F3/Gurobi
(≈ 5.6 min); the A6′ line rule; the global pool size and crash release. **Full suite (2026-10-27, not
recorded to `test_runs.tsv`): 439 passed, 0 failed.** Known answers:

- **F0:** prunes the two dead routes. The repair rule re-opens the cheaper one (weight 1.5 against 3) and
  growth is restored. In KO, a re-opened reaction keeps its full GPR.
- **Gene-state convention:** an off isozyme does not buffer, so KO of the only expressed isozyme is lethal.
  The generic model is buffered.
- **GIMME:** a low-expression reaction that is required is kept; unused low-expression reactions are closed;
  objective = 4.5.
- **FASTCORE:** exact support, including a reversible core reaction that is used only after a flip, and the
  blocked reaction removed by fastcc.
- **E-Flux:** bounds 1000 and 750 (OR = sum); isozyme KOs give ratios of 1/3 and 2/3; the floor of 10 is hit.
- **CORDA:** the low-confidence path is chosen over the negative one.
- **Loss rule:** a single lethal loss is model-contradicted. An isozyme pair that is lethal only jointly is
  resolved by dropping the tied member with the lower Entrez ID. The result is identical over 20 shuffled
  input orders. Accepted-loss closures are never re-opened.
- **Other checks:** the presence score matches boolean pruning (300 random cases); the sweep shortcut equals
  brute force; known answers for the loss-rule definition, including the v1 rule, CN 0.5 and mutation-only
  lines; panel determinism and training-only membership; the permutation preserves each gene's values.

## 3. M3 pilot results (Human2, medium M0, Gurobi 1 thread, Seed 11)

### 3.1 Repair rule on the 12 lines the compute probe killed under F0 (`repair_probe12.tsv`)
These are the first 12 default rows of the expression file. 8 are unscreened; the 1 test and 2 validation
lines are masked in the output.
- WT growth before repair was 0 in **12/12**, and repair restored growth in **12/12**.
- WT after repair was 0.621–0.876 of the generic WT.
- **25–43** reactions were re-opened (1.2–2.8% of the closed reactions).
- The repair takes a median of 10.5 CPU-s per line: 11 LPs (growth + GR + 9 ER).

### 3.2 Loss rule (`loss_rule_panel.tsv`)
On the 30 panel lines, the drop gave identical results under **6 input orders per line, for all 30 lines**.
- Lost model genes per line: median 4, max 23.
- 2 lines had a model-contradicted loss (ACADL, BBOX1, both generic-lethal biomass artefacts).
- No line needed a joint drop; step 4 is exercised by the toy test only.
- Cost: ≤ 0.73 CPU-s per line.

### 3.3 Per-line cost, memory and the permuted control (`pilot_builds_gurobi_lp.tsv`)
3 panel lines × {real, permuted} RNA, R1, plus the genomic re-sweep on real builds: 36 builds, **0 failed**,
including all 12 permuted-expression control builds.

| Family | build CPU-min (mean / max) | A4 re-sweep CPU-min | extract / repair / sweep / tasks CPU-s | re-opened (max % of rxns) | open rxns | WT / generic | KO lethal | tasks passed (of 57) | peak RSS GB |
|---|---|---|---|---|---|---|---|---|---|
| F0 | 1.80 / 2.35 | 1.75 | 0 / 11.8 / 65.2 / 29.9 | 25–40 (0.31%) | 10,769–11,262 | 0.625–0.876 | 165–172 | 44–56 | 1.42 |
| F2 | 2.13 / 2.66 | 2.00 | 1.3 / 11.5 / 83.4 / 29.7 | 0 | 10,793–11,268 | 0.982–0.999 | 165–170 | 44–56 | 1.42 |
| F3 | 2.46 / 3.00 | 1.25 | 69.3 / 16.2 / 45.0 / 15.9 | 88–93 (0.73%) | 6,378–6,721 | 0.968–0.984 | 173–182 | 40–52 | 1.43 |
| F4 | 3.32 / 3.60 | 3.23 | 0 / 1.9 / 160.6 / 34.9 | — (not repaired) | 12,734–12,748 | 0.011–0.017 | 154–155 | 51–53 | 1.42 |
| F3b | 104.6 / 91.1 (1.74 / 1.52 CPU-h) | — | 6,214 and 5,402 / 18 / 29 / 18 | 0 | 6,819–7,126 | 0.870–0.875 | 172–178 | 45–52 | 1.59 |

- The DESIGN unit cost of ≈ 1.5 CPU-min per line does not hold. The measured cost is 1.8–3.3 CPU-min per
  build, including the 57-task profile (0.3–0.6 min), which the method-vs-tissue partition needs.
- **One-off cost:** the FASTCC consistent set takes 24.6 CPU-min per parent × medium (cached; 8,997 of
  12,748 reactions consistent on M0). Every W2 medium that F3 uses needs its own FASTCC.
- **Repair-size gate (≤ 5% of reactions):** the maximum is 0.73% (F3).
- **F4 binary KO profiles equal the generic model's** under the ratified scaling. The reviewer found the
  F4 lethal set (154 genes) identical to the generic Human2/M0 set on 3 lines: with WT growth this low, the
  floor of 10 almost never binds after a KO. F4's line-specific signal is therefore in its partial ratios
  only, and in H4.3 (binary KO profiles) F4 is "method-dominated" by construction. F4's lethal calls come only
  from closing GPR-False reactions (A4 a).
- **E-Flux WT growth is 1–2% of generic.** With 1000·cap/max-cap scaling and a floor of 10, most bounds
  sit at 10. F4 growth is in its own units, so the F4 failure rule is "no growth (≤ 1e-6)", not the
  1%-of-generic repair floor (§6).

### 3.4 GLPK = Gurobi per family (`glpk_vs_gurobi.tsv`, line ACH-000957)
| Family | full pipeline: reaction-set Jaccard / genes with different ratio_r | same model swept on GLPK: genes differing | GLPK/Gurobi build CPU (like-for-like) | Verdict |
|---|---|---|---|---|
| F0 | 0.9998 / 0 | 0 | 0.87 | **equal; GLPK allowed** |
| F2 | 1.0 / 0 | 0 | 0.98 | **equal; GLPK allowed** |
| F3 | 0.9873 / 93 (max \|Δ\| 0.052) | 0 | 1.53 | **not equal** (FASTCORE alternate optima): Gurobi only |
| F4 | 1.0 / 1 (\|Δ\| = 1e-6, one rounding unit) | (n/a) | 0.57 | **not equal under strict equality**: Gurobi only, unless the PI accepts \|Δratio_r\| ≤ 1e-6 |

- Sweeps on an identical model are equal on GLPK for every repaired family; the differences come from
  alternate optima in extraction and repair.
- `w4_stage_a.py` therefore fixes **one solver per family**, so optima are never mixed within a family
  (now `context.DECLARED_SOLVER`, A4 c).
- Every KO table carries a `solver` column.
- **Superseded by A4 (c):** the declared solvers are GLPK for F0, F2, F4 (and F5) and Gurobi for F3; see §0.

### 3.5 CORDA (F3b) cost band (`pilot_builds_gurobi_corda.tsv`)
Two panel lines (ACH-000957, ACH-002024), R1 real RNA, 2 Gurobi processes run in parallel.
- **CORDA extraction alone took 1.73 and 1.50 CPU-h; per model, 1.74 and 1.52 CPU-h**, so the band is
  **> 1 CPU-h per model → A6 is DROPPED** under the frozen row-A6 rule, and the drop must be reported.
- Both builds are valid: no failure, WT 0.870–0.875 of generic, 0 reactions re-opened, 6,819–7,126 open
  reactions, peak RSS 1.59 GB.
- Full A6 would cost ≈ 210 CPU-h (120 models); the 15-line fallback ≈ 52 CPU-h.
- The method-vs-tissue test therefore rests on F2 GIMME, F3 FASTCORE and F4 E-Flux, plus F0 as the naive
  reference. That still meets "≥ 3 established non-ftINIT families" (D2).
- If the PI wants CORDA anyway, the r3 verifier's recommendation (a subset rather than a drop) would need a
  §17 amendment. *Superseded by A7: CORDA runs as row A6′ (§0.2).*

### 3.6 RSS → worker counts
- **Per-worker peak RSS:** 1.43 GB on Gurobi and 1.56 GB on GLPK, measured as VmHWM. `ru_maxrss` cannot be
  used here: on Linux it survives fork + exec, so spawned workers report their parent's peak (3.4 GB).
- **Host:** 134.8 GB total, but only **63–77 GB available** across this session (shared host). Keeping
  20 GB free (§11.4) leaves 43–57 GB, which is **27–35 workers maximum**.
- **Recommendation, frozen as config:** **4 Gurobi workers (F3, F4, F3b) + 20 GLPK workers (F0, F2)**.
  This needs ≈ 37 GB and stays at 24 LP workers until availability is re-checked at launch. Do not go to
  28 unless ≥ 64 GB is available (24 × 1.56 + 4 × 1.43 = 43 GB, plus 20 GB reserve).
- The script aborts below 10 GB available.

### 3.7 Stage-A projection (frozen table, n_panel = 30, n_p = 10 with ProCan)
Builds:
- A1: 480; A2: 120; A3: 240; A4: 600 re-sweeps; A7: 120.
- A5: not run (§6).
- A6: per the CORDA band below.

CPU, using the measured per-family costs (protein builds assumed to cost the same as RNA builds):

| Row | CPU-h |
|---|---|
| A1 | 19.4 |
| A2 | 4.9 |
| A3 | 9.7 |
| A4 | 20.6 |
| A7 | 4.6 |
| **LP total** | **59.2** (+ A6) |

This is within the 100 CPU-h Stage-A budget. Wall clock is set by the Gurobi lane (lane CPU-h / workers):

| Solver assignment | Gurobi lane CPU-h | GLPK lane CPU-h | Wall (4 Gurobi + 20 GLPK) |
|---|---|---|---|
| strict: F0, F2 on GLPK | 33.9 | 23.6 | **8.5 h** (+ A6) |
| F4 also on GLPK (if the 1e-6 tolerance is accepted) | 14.2 | 34.8 | **3.5 h** (+ A6) |
| all Gurobi, 4 sessions | 59.1 | — | 14.8 h |

A6 was dropped by A4 (e) and is superseded by A6′ (A7, §0.2); the table above predates A5 and A6′ (current projection: §0.4).

## 4. Ready-to-run Stage-A sweep
```
# 1. commit src/nextgem/context*.py, scripts/w4_*.py, tests/test_context.py; create the pinned worktree (§11.4)
# 2. (once) .venv/bin/python scripts/w4_prepare.py      # panel, thresholds, FASTCC cache (already in results/w4/)
.venv/bin/python scripts/w4_stage_a.py --dry-run
.venv/bin/python scripts/w4_stage_a.py --gurobi-workers 8 --glpk-workers 20 --budget-cpu-h 110   # see §0.4
```
- **Outputs:** `results/w4/stage_a/{builds.tsv, ko_A1|A2|A3|A4|A5|A6p|A7.parquet, accepted_losses.json,
  a5_grid_map.tsv, a6prime_lines.tsv, MANIFEST.tsv}`.
- **Checkpoints:** per build in `stage_a/tmp/jobs_<code hash>/` (git-ignored). A rerun with the same code
  resumes.
- **Guards:**
  - refuses a dirty tree (provenance), and refuses pilot options without `--pilot` and `--out`;
  - aborts at 2× the declared CPU budget, or below 10 GB available memory;
  - Gurobi processes hold a slot of the shared pool (global size 8);
  - one build per process, with a per-job watchdog (§0.1);
  - panel lines are asserted to be training-only.
- **Smoke test:** §0.3. The pre-review smoke directories (`stage_a_smoke_a5/` and the first `stage_a_smoke/`) were replaced by the post-fix run.

## 5. Implementation choices made where DESIGN is silent (for PI ratification)

1. **KO gene-state convention** (presence families). Genes lost under the genomic constraint are
   non-functional. An *open* reaction whose GPR is still True without the off genes uses that pruned GPR,
   so an off isozyme does not buffer; this mirrors ftINIT `removeGenes = true` in Human-GEM's DepMap
   pipeline. An open reaction kept or re-opened despite a False pruned GPR keeps its full GPR. Without
   this convention, the a priori isozyme-switch expectations (HK2, LDHA, CDS2, NAMPT) cannot appear.
   E-Flux keeps full GPRs; a KO zeroes the gene's capacity.
2. **"No evidence" in GPR scores:**
   - AND = min over known genes;
   - OR = max over known genes, but unknown if any child is unknown and every known child is < 0.
   This makes `s_r < 0` ⇔ GPR False under the off genes (tested).
3. **FASTCORE:**
   - The core is exactly "reactions scoring ≥ threshold". Biomass is not forced into the core; the repair
     rule restores growth (88–93 reactions, 0.7%).
   - ε = 1 (cobra fastcc default); support tolerance 1e-6.
   - A core reaction that cannot be activated even as a flipped singleton is dropped from J and counted
     (`unreachable_core`), not fatal. **Correction (review NOTE 1a): the pilots had 21–40 per build**
     (real RNA: 28, 25 and 23 on ACH-000957, ACH-002024 and ACH-003383; permuted: 32, 40 and 21), out of
     4,555–4,983 core reactions (≈ 0.5%).
     - FVA shows that these consistent core reactions cannot reach |v| ≥ 0.99 on M0. They are
       micronutrient-limited: vitamin B12, vitamin E, retinol, some transport and glycerophospholipid
       reactions.
     - They are therefore excluded by construction at the ratified ε = 1. Full coverage would need a smaller
       ε, which would be an amendment.
     - Also, FASTCORE closes every reaction outside its support, exchanges included, so each F3 line's
       effective medium is line-specific (review NOTE 1d).
4. **CORDA:**
   - Gene confidence bins: s ≥ 1 → 3; [0, 1) → 2; [−1, 0) → 1; < −1 → −1; no evidence → 0.
   - Package defaults n = 3, penalty 100, support 5.
   - Biomass is not forced (conf 0).
5. **E-Flux:**
   - Linear value 2^s (expression relative to threshold); 0 if off.
   - Bound = max(10, 1000·cap/max cap), max taken within the line.
   - An OR with an unknown gene keeps the default bound.
6. **Failed-build rule for F4:** WT growth ≤ 1e-6 (no repair rule; scaled units).
7. **Thresholds:**
   - R2 = P25 of all (model gene × training line) log2(TPM+1) values = 0.370, which is more permissive
     than R1 (30% of values fall below R1).
   - R3 = training mean clipped to [0.370, 5.227].
   - Protein: R2/R3 on log2 ProCan intensities from 363 training lines; protein R1 = "detected".
   - "Not detected" counts as absent only for soluble enzymes. The proxy for soluble is: no transport
     reaction and no Cell membrane / Inner mitochondria localisation (2,148 of 2,870 genes).
8. **Panel draw:**
   - `default_rng(SEEDS[0])` permutes the 12 eligible training lineages (sorted by name) and takes the
     first 10.
   - 3 lines are then drawn from each lineage's expression-profiled training lines (sorted).
   - No fill was needed.
   - Lineages: Bone, Bowel, Cervix, Esophagus/Stomach, Head and Neck, Kidney, Myeloid, PNS, Skin, Uterus.
     n_p = 10 lines with ProCan, and all 30 have CN.
9. **Permuted control:** each gene's log2 values are permuted across the 30 panel lines (SEEDS[0]); the
   real-data thresholds are used.
10. **A7:**
    - GIMME objective perturbation adds 1e-3·U(0,1) to each reaction cost.
    - FASTCORE noise is N(0, 0.25) on log2 values (`--noise-sd`).
    - DESIGN gives no magnitudes for either.
11. **Loss rule ("kills WT growth"):** growth < 1% of generic WT (the repair floor). The drop is evaluated on
    the parent with the line's medium, so it is family-independent.
12. **Medium:** everything is on M0 (Ham's) until W2 fixes per-line media. The repair, loss-rule, GIMME and
    FASTCC inputs are all per parent × medium.

## 6. Open items and infeasibilities
- **Resolved by A4 (d), see §0.** *Original note:* **A5 (F5 LCP) cannot run as frozen.** DESIGN fixes "24 points (coarse 12 + local refinement 12)" but no
  grid values, and the local refinement depends on labels. The builder exists (`extract(..., "F5",
  off_override=...)`); the grid needs a PI specification. A5's 18 CPU-h are not in the projection.
- **Resolved by A4 (c): F4 is declared GLPK.** *Original note:* **The F4 GLPK tolerance** decides between ≈ 8.5 h and ≈ 3.5 h of Stage-A wall time (plus A6).
- The 57-task profile per build uses closed boundaries (Human-GEM task convention).
- **P57 per-line cost** (an M3 row item) was not measured here; it belongs to the W3/§6.8 owner.
- **Provenance:** `w4_stage_a.py` runs under `nextgem.run.provenance()` (§0.1, B3). The M3 pilots
  (`w4_pilots.py`) ran from a dirty tree; their numbers are pilot measurements, not production outputs.
- `/data/tmpA` is shared. Stage-A outputs are projected at about 15 MB of KO parquet plus a few MB of
  summaries.

## 7. Runtimes of this session's pilots (wall; CPU in the logs under `logs/`)

| Step | Wall | Notes |
|---|---|---|
| `w4_prepare` | 25 min | mostly FASTCC: 24.6 CPU-min |
| `repair12` | 226 s | 1 worker |
| `lossrule` | 97 s | 1 worker |
| `cost` (4 LP families) | 1,799 s | 4 workers, 36 builds |
| `glpk` | 372 s | 5 GLPK workers |
| CORDA | 6,291 s | 2 workers, 1 model each |
| first smoke test (pre-review) | 441 s | 2 workers |
| A5 smoke (pre-review) | 989 s | 4 GLPK workers |
| full smoke, all rows (post-fix) | 5,840 s (CORDA build 95.6 CPU-min) | 2 Gurobi + 5 GLPK workers |
