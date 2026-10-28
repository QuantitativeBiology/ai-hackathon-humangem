# W8: the P1 resource (DESIGN §9 W8; RESEARCH_TREE D11)

**Author:** W8 resource engineer agent (Claude Opus). **Date:** 2026-10-28.
**Status:** code and tests pass, two pilots pass every gate (4 models), and MEMOTE ran on a 4-model sample. The full build is about 3.3 CPU-h, above the 1 CPU-h limit, so it must run from a pinned clean worktree (§11.4). It has **not** been run, and the author has committed nothing. DESIGN.md, FROZEN.tsv, the ledgers and `evaluate_*.py` were not edited.

## Status update (2026-10-28, after the pinned build; supersedes the projections in §§4–6 where they differ)
- **Full build** (pinned b200324, committed 6a71f7f; `results/w8/p1/summary.json`).
  - 40 models, 1,361 lines × 2 parents.
  - 2.76 CPU-h (9,922 CPU-s), 1,610 s wall at 8 GLPK workers; 202–252 CPU-s per model; 77.9 MB of models.
  - **39/40 models pass all six gates.** The one failure is `medium-20cc933a4b852b6a_human2-2.1.1_P1` (RPMI + testosterone, 1 line). It fails `random_ko_match` on gene ENSG00000099817: the in-memory model gives 0.0 and the reloaded model gives −8.04e-9, above the 1e-9 absolute floor (RESEARCH_TREE 9.1).
- **Gate-miss investigation** (`scripts/w8_ko_noise_check.py` → `results/w8/ko_noise/`, 212 CPU-s). The gate's own 20-gene sample was re-solved on both models.

  | Solver | ENSG00000099817 KO, in-memory model | ENSG00000099817 KO, reloaded model | Max \|Δ\| over all 20 genes |
  |---|---|---|---|
  | GLPK, default tolerances | −1.5e-11 | −8.04e-9 (reproduced) | 8.0e-9 |
  | GLPK, feasibility tolerance 1e-10 | 3.3e-13 | −2.5e-13 | 1.0e-12 |
  | GLPK exact (rational arithmetic) | **0** | **0** | 0 |
  | Gurobi, tolerances 1e-9 | **0** | **0** | 0 |

  - The −8.04e-9 value is a violation of the biomass lower bound (0). That violation is below GLPK's default primal feasibility tolerance of 1e-7, so it is LP noise. The KO is exactly lethal in both models, and the two models are equivalent.
  - **The gate was not changed.** The miss stays recorded as a documented gate failure with its cause. Changing the rule needs a PI decision; options include tying the floor to the solver feasibility tolerance, clipping negative growth to 0 as `ratio_r` already does, or re-solving mismatches exactly.
- **Exit code.** The committed script returns 3 when any model fails a gate: `run()` returns 3 when `n_models_pass < n_models`, `main()` returns that value and `SystemExit(main())` exits with it. So the observed "exit 0" most likely came from the invoking wrapper, for example a pipe to `tee` or a later command in the same shell line. I could not see that command.
  - The exit logic is now explicit: `resource.gate_exit_code(models)` lists the failing models on stderr and returns 3.
  - Both build scripts use it and print `exit code N (k/n …)`.
  - Tests: `test_w8_gate_exit_code` and `test_w8_build_main_propagates_exit_code`.
- **Extension to DESIGN's scope** (`scripts/w8_extend_unscreened.py`, `results/w8/unscreened/`; dirty tree; 3 runs totalling 2,712 CPU-s ≈ 0.75 CPU-h).
  - Scope: **753** unscreened DepMap 26Q3 models with a default expression profile and a cancer lineage, exactly DESIGN's `resource_models_not_benchmarkable`.
  - Media come from W2's frozen rule, `media.assign_line_media`, called on these lines; `media.py` was not edited. That gives 666 onboarded and 87 imputed lines. By stratum: onboarded_sensitivity 409, unmappable 194, imputed 87, fallback_sheet_missing 56, mixture_on_disk 7. Every non-documented line gets the flagged classic RPMI + 10% FBS.
  - 17 medium keys. Every key shared with the screened cohort reproduces that cohort's W2 signature (asserted).
  - 4 new signatures per parent, covering 6 lines each:
    - DMEM + pyruvate + 2-ME;
    - RPMI + cortisol + estradiol;
    - 2 serum-free DMEM:F12 variants.
  - Each new signature got a B-loo reference sweep with the M6 engine (8/8 ok, 565 CPU-s) and a model; **8/8 pass all six gates** (2,144 CPU-s). The 2 serum-free models per parent keep the mevalonate core lethal, so they carry no W3-MEV.
  - The resource now holds **48 models (93.3 MB), 2,114 lines and 4,228 manifest rows** (`screened` = true for 1,361 lines, false for 753).
  - Screened rows and the 40 existing files are unchanged, which was verified. `resource/README.md` was regenerated.
  - The 8 new files were built from a dirty tree (`git_dirty: true` in their notes). For a clean release, re-run `--step sweeps` and `--step build` from a pinned worktree (≈ 0.75 CPU-h).
- **Determinism finding.** cobra keeps group members in a Python set, so the SBML bytes and sha256 of an identical model differ between runs. A pilot file and the released file differ only in member order, apart from notes.
  - `write_sbml_gz` now sorts group members (`canonical_group_members`; test `test_w8_canonical_group_members`).
  - The released files are content-identical to what the fixed writer produces, but not byte-reproducible.
- **MEMOTE on 10 models** (both parents × the 5 most-used signatures: classic RPMI, DMEM, IMDM, DMEM:F12, DMEM:F12 + cortisol + estradiol). Run in `.venv_memote` with the same skips as before (7 heavy tests plus the unscored `test_matrix`).
  - 6 models run on the release files: 1,670 CPU-s, 28 min. The other 4 are the pilot runs, whose model content is identical to the release files once group order is canonicalised and notes are ignored (checked).
  - Scores: **Human2 0.638** for all 5 media; **Human1 0.630** for all 5. Consistency 0.425; gene annotation 0.467; SBO 0.818.
  - Summaries: `results/w8/memote/memote_summary.json` (4 pilot models) and `memote_summary_release6.json`.
- **FROG** (`fbc_curation` 0.3.2 in `.venv_memote`, with pymetadata pinned to 0.5.12 because 0.8.2 breaks its import; cobrapy curator only).
  - Classic-RPMI Human2 model: **complete**, 2,383 CPU-s (0.66 CPU-h), 39.8 min wall, 1 process.
    - Timings: objective 16 s, gene deletions 197 s, reaction deletions 826 s, FVA 1,337 s.
    - Every LP returned optimal: 2,870 gene deletions, 12,748 FVA rows, and all reaction deletions.
    - Report: `results/w8/frog/medium-4cd587facad2c6fb_human2-2.1.1_P1/` (FROG TSVs + metadata, 2.8 MB).
  - Known bug in `scripts/w8_frog.py`: `frog_summary.json` reads the wrong keys from `FrogReport.to_dfs()` (`objective`, `gene_deletion`, `reaction_deletion`), so it shows `objective` and two deletion counts as null. The TSVs are complete. The fix is deferred because code edits are paused.
  - The Human1 model is still to run: ≈ 0.7 CPU-h.

**Label hygiene.**
- The session-start probe (`head -c 50` of the raw CRISPR gene-effect file) was **blocked** by the hook.
- W8 is label-free. It never opens `labels.py`, a label file, a label-derived list or `data/locked/`.
- From `data/depmap/` it reads only `Model.csv`, and only the columns `ModelID`, `OncotreeLineage` and `OncotreePrimaryDisease`.
- The label-copy scanner over `results/w8`, the new code and the session scratchpad found **0 hits** in 1,113 files.
- The scanner's `detect()` finds no ModelID × gene layout in `manifest.tsv` or `models.tsv`, so the pre-commit test-keyed rule does not apply.

## 1. What P1 is and what is exported
**P1** (D11) = generic parent + W3 `B-loo:cofactor_pool_biomass` on the W2 documented medium M2(f=1), gate-free.
- The variant removes `MAM10012c` from `MAR13082`; the spec is in `results/w3/structure/variants_<parent>.json`.
- There is no context family, so nothing is pruned or repaired.
- P1(H1) is the selected pipeline. P1(H2) is its Human2 counterpart.

**Scope.** Every line with a medium assignment: the **1,361** lines in `results/m6/scores/line_assignment.tsv`, for both parents (2,722 line × parent rows).
- The lines map to **20 distinct M2(f=1) medium signatures per parent, so 40 models**: 19 serum media and 1 serum-free.
- **Release order:** the build runs Human2 first, then Human1 (DESIGN W8).

**Model construction** (`resource.build_medium_variant_model`):
- `parent.copy()`, then `model.medium = bounds` (every other uptake is closed; secretion is unchanged);
- then `Variant.apply` (permanent), then objective `MAR13082`.
- These are the same steps, in the same order, as `w3.ko_sweep` uses inside its context.

**Where the bounds come from.**
- They are rebuilt with `w2_engine.build_sweeps` through `m6_sweeps.bounds_by_id`, exactly as M6 did.
- Each job asserts that the bounds hash to the W2 signature.
- Closed reactions are kept with bounds (0, 0); 8 per model come from the parent.
- `version` is the Human-GEM release: the SBML `version` note and the `version` column.

## 2. Storage scheme (choice and justification)
**Chosen: one gzipped SBML per (parent, medium signature), plus a line manifest.**

| Path | Content |
|---|---|
| `resource/v1.0.0/medium-<sig>_<parent>-<ver>_P1.xml.gz` | 40 models: SBML L3V1 FBC v2, deterministic gzip |
| `resource/models.tsv` | one row per model: sha256, size, counts, flags, the 5 DESIGN gates + `gate_reload_m6` (pass/FAIL), job status, CPU-s |
| `resource/manifest.tsv` | one row per (line, parent), 2,722 rows. Each row gives the model file and its sha256, plus the line's own metadata: medium key, documented text, source, **stratum = tier** (attribution / onboarded / mixture / imputed / sheet missing / unmappable), fallback flag ("classic RPMI + 10% FBS"), serum-free flag, documented serum %, lineage, `cancer` (DESIGN scope rule) and `screened` |
| `resource/README.md` | layout, line access, licence (generated by the build) |

**Per-line access** (DESIGN naming `{ModelID}_{parent}-{version}_P1.xml.gz`):
- `scripts/w8_line_sbml.py --lines ACH-… --parent human2`, or `nextgem.resource.load_line("resource", ModelID, parent)`.
- Either one reads the shared file (sha256-checked against the manifest) and adds the line's notes (`model_id`, `line_medium_*`, `shared_model_file/sha256`). Nothing else changes.

**Why not one file per line.**
- **Disk.** Per-line files would take 2,722 × ≈ 2 MB ≈ 5.3 GB. That fits the 10 GB cap, but on a 98%-full shared volume (89 GB free), and it would be 34× duplication of 40 distinct models.
- **Compute.** Materialising measured 34 s wall per line, so ≈ 26 CPU-h in total, plus gating 2,722 files instead of 40.
- **Honesty about specificity.** PROMPT.md asks to "deliver a resource of optimised, curated GEMs for DepMap 26Q3 cancer cell lines" and to "deliver generated context-specific models in SBML". P1 was selected over the context families (M5: not selected), so its only line-specific input is the culture medium. Shipping 2,722 near-identical files would suggest more line specificity than there is (cf. flag W2-SHUF).
- **Every line still gets an SBML model.** The manifest gives each line its own row, its own medium metadata and a sha-pinned file, and a per-line file is one command away.
- Total resource size: **≈ 78 MB** (the 4 pilot models average 1.95 MB) plus ≈ 1 MB of manifest.

## 3. SBML notes (every model)
- `pipeline` (P1 per D11, full description), `pipeline_id` = P1;
- `human_gem_release` / `parent_version` / `version`;
- `medium_signature`, `w2_sweep_id`, `medium_step` (M2(f=1), gate-free), `medium_keys`, `medium`;
- `medium_source_tier` / `line_strata`: tier counts of the lines sharing the file, e.g. `attribution:539;fallback_sheet_missing:210;…`;
- `serum`, `n_uptake_open`, `n_lines`;
- `w3_variant` (JSON spec);
- `git_commit` + `git_dirty`;
- `input_shas`: SHA-256 of the parent SBML, task file and variant spec; the media library (compositions, exchange map, line media); the floor tables; line_sweeps / line_assignment; the M6 catalog, the MEV table and the reference sweeps; `Model.csv`;
- `flags`, `mevalonate_core_lethal`;
- `repaired_reactions` = "none (generic parent; no context family, nothing pruned or re-opened)";
- `licence`, `release_status`, `build_time_s`.

**Flags.** `W3-D1;W2-SHUF;W2-MOPEN` on every model. `W3-MEV` is set from the model's own M6 sweep (`results/m6/scores/w3_mev.tsv`, mevalonate core non-lethal). It therefore lands on the 19 serum media and not on the serum-free DMEM:F12 model.

## 4. Validation gates (`resource.validate_model`; every model)
1. **libSBML.** No errors beyond the raw parent; both parents have 0 errors.
2. **cobra round-trip.** Identical reaction, metabolite and gene sets, bounds and GPRs.
3. **WT growth.** Within 1e-6 relative.
4. **Random KOs.** 20 KOs (SEEDS[0]) within 1e-6 relative or 1e-9 absolute.
5. **Viability tasks.** GR + 9 ER on the reloaded model.
6. **`reload_vs_m6`** (requested here). The reloaded model's WT growth and KO ratio vector on a 40-gene sample are compared with the in-memory model and with the M6 sweep of the same (parent, B-loo, signature).
   - The sample is 25 genes with M6 ratio < 1 and 15 at ratio 1 (SEEDS[0]).
   - Ratios must match at 1e-6 (the sweeps' 6-decimal rounding).
   - WT must match M6 at 1e-5 relative, because M6 stores WT at 6 significant digits.

**Stub fix.** `viability_tasks` keyed its results by the raven task id, but raven ids are categories: the 9 ER tasks all have id "ER". So the stub (and the demo) checked only the last ER task. Keys are now ER1..ER9 plus GR (`task_keys`), and a test asserts all 10 are checked and pass on Human2.

**A10 limits.** GLPK. Each LP has a 300 s limit with one retry (`_growth` now goes through `w3.solve_status`), and each model job has 1,800 s CPU (`ITIMER_PROF`). A timeout or any exception is recorded as a failed build in `models.tsv` and is never re-run with relaxed limits.

### Pilot results (dirty tree, `results/w8/tmp/pilot/`, tagged "pilot, not a release")
| Pilot | Models | Gates | CPU (MANIFEST) | Wall |
|---|---|---|---|---|
| `--pilot 1 --workers 2` | H2 + H1 on classic RPMI (1,050 lines each) | 2/2 pass all 6 | 874 s | 716 s |
| `--pilot 2 --workers 2` | + H2/H1 on DMEM (169 lines) | **4/4 pass all 6** | 1,548 s | 798 s |

**Reload gate detail.**
- Max |Δratio| is 0 both between the reloaded and in-memory models and between the reloaded model and M6.
- WT on classic RPMI is 46.909155 (H2) and 46.863208 (H1); M6 has 46.9092 and 46.8632.
- All 10 viability tasks pass.

**Per-model cost.** About 240–305 CPU-s, of which libSBML `checkConsistency` is ≈ 125 s. The parent libSBML checks add ≈ 200 s wall, run in parallel for the two parents.

## 5. Full build (to run from a pinned clean worktree; ≈ 3.3 CPU-h)
```
# after the PI commits the files below; data/ symlinked as for earlier pinned runs
cd /data/tmpA/nextgem_runs/<commit>
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=$PWD/src \
  /data/tmpA/nextgem/.venv/bin/python -u scripts/w8_build_p1.py --workers 8 --release 1.0.0
```
- **Projection.** 40 jobs × ≈ 270 CPU-s, plus the parent checks and 16 worker loads: **≈ 3.3 CPU-h, ≈ 35–45 min wall** at 8 workers, with peak RSS ≈ 2.6 GB per worker.
- **Outputs.** In the worktree: `resource/` (≈ 79 MB), `results/w8/p1/{gates.jsonl,summary.json}` and a MANIFEST row in `results/w8/MANIFEST.tsv`.
- **Guards.**
  - `--allow-dirty` is refused without `--pilot`.
  - A non-pilot run refuses an imported `nextgem` that is not next to the script.
  - `provenance()` refuses a dirty `src/`, `scripts/` or `tests/`.
- **Exit code.** 0 only if every model passes every gate; otherwise 3.

## 6. MEMOTE and FROG
**MEMOTE** is installed in a separate venv, `/data/tmpA/nextgem/.venv_memote`, with memote 0.17.0 and its own cobra 0.32.1, optlang and swiglpk; it takes 520 MB.
- The main `.venv` is untouched. `.venv_memote/` was added to `.gitignore`.
- `scripts/w8_memote.py` refuses to run under any other interpreter and imports no project code.
- **Skipped tests** (`SKIP`, documented in the script):
  - the FVA- and per-metabolite LP sweeps and the MILP (several CPU-h per model on Human-GEM);
  - the unscored `test_matrix` module (dense rank/nullspace, ≈ 25 CPU-min of the 33-min timing run).
- One skipped test is scored (`test_find_reactions_unbounded_flux_default_condition`, consistency). MEMOTE counts it as 0, so the totals are a **lower bound** under this skip list.

**Results.**
- Timing run, H2 classic-RPMI P1 model, matrix tests included: 2,008 CPU-s. **Total score 0.638**:
  - consistency 0.425;
  - annotation: metabolites 0.777, reactions 0.764, genes 0.467, SBO 0.818.
- 28 of 77 tests failed, mostly Human-GEM annotation/SBO conventions and biomass checks (GAM, NGAM, precursors). The failures are listed in `memote_summary.json`.
- **Sample of 4** (the pilot P1 models, H2/H1 × classic RPMI/DMEM; final skip list): `results/w8/memote/` (`memote_summary.json` + one full report `*.memote.json.gz` per model, 2.6 MB).
  1,295 CPU-s (0.36 CPU-h), 21.8 min wall, peak RSS 2.8 GB.

  | Model | Total | Consistency | Ann. met | Ann. rxn | Ann. gene | SBO | Tests passed / failed / skipped |
  |---|---|---|---|---|---|---|---|
  | H2 P1, classic RPMI | 0.638 | 0.425 | 0.777 | 0.764 | 0.467 | 0.818 | 29 / 28 / 15 |
  | H2 P1, DMEM | 0.638 | 0.425 | 0.777 | 0.764 | 0.467 | 0.818 | 29 / 28 / 15 |
  | H1 P1, classic RPMI | 0.630 | 0.425 | 0.705 | 0.719 | 0.467 | 0.818 | 28 / 26 / 15 |
  | H1 P1, DMEM | 0.630 | 0.425 | 0.705 | 0.719 | 0.467 | 0.818 | 28 / 26 / 15 |

  - Scores do not depend on the medium.
  - The score equals the timing run's, which confirms that the matrix module is unscored.
  - The failures are parent-model properties, not the P1 edits. They include annotation/SBO coverage, charge/mass balance, stoichiometric consistency, GAM/NGAM and biomass precursors, and dead-ends/orphans.

**FROG.** Not run; this is the plan.
- The FROG/fbc_curation report needs the objective, FVA on every reaction (≈ 26k LPs), all gene deletions (2.9k) and all reaction deletions (13k). That is ≈ 2–4 CPU-h per Human-GEM model on GLPK.
- Plan: install `fbc_curation` into `.venv_memote` (a separate venv, ≈ +50 MB), then run FROG on the 2 classic-RPMI models (H1, H2) from the pinned worktree, ≈ 6 CPU-h.
- DESIGN's "10 models" would be ≈ 20–40 CPU-h, which needs a PI decision against the remaining budget.
- MEMOTE on 10 models with this skip list is ≈ 1.4 CPU-h, which also goes to the worktree.

## 7. Deviations and open points for the PI
1. **Scope: 1,361 screened lines, not DESIGN's 2,030.** DESIGN's 2,030 = default expression ∩ cancer, including 753 unscreened lines.
   - W2's medium assignment exists only for the 1,361 CRISPR lines, and the brief said "lines that have a medium assignment".
   - The unscreened lines would need a W2 assignment from `Model.csv` onboarded media. Most would map to existing signatures, but any new signature needs a new model and its own M6 reference sweep.
   - `manifest.tsv` carries `cancer`, so the 1 Fibroblast line (and any non-cancer line) is visible rather than dropped; `screened` is "true" for all rows.
2. **Per-signature files instead of per-line files** (§2). The per-line DESIGN naming is available on demand.
3. **Built before T1.** DESIGN says the build starts right after T1. P1 has no fitted parameters (the medium assignment is label-free and B-loo is fixed), so these models are exactly the ones T1 scores.
   - Default `release_status`: "release candidate (P1 per D11; built before T1)". Version 1.0.0 is the author's choice and can be changed with `--release`.
4. **Gate 6 is added**, and the stub's ER-task bug is fixed (§4). With the fix, the demo's viability gate was checking only the last ER task. The demo exports were not re-run.
5. **Model-level `benchmarked` = true** (every line of every signature is screened). The line-level flags are in `manifest.tsv`.
6. **The licence NOTICE for HPLM/Plasmax is not needed.** P1 uses documented media only.
7. **Two pilot runs from the dirty tree** (0.67 CPU-h together) plus the MEMOTE timing run (0.56 CPU-h) and the MEMOTE sample, each under 1 CPU-h. Their outputs are git-ignored (`results/w8/tmp/`), except the MEMOTE sample in `results/w8/memote/` (small JSON).

## Files
- Code: `src/nextgem/resource.py`, extended with the P1 build, reload gate, line manifest, `load_line`/`materialise_line`, `task_keys`, notes `extra`/`flags` and A10 `_growth`. Scripts: `scripts/w8_build_p1.py`, `scripts/w8_line_sbml.py`, `scripts/w8_memote.py`.
- Tests: `tests/test_w8_p1.py`, 11 tests. Together with `tests/test_w8_resource.py`, 25 pass in 44 s.
- Provenance: `results/w8/MANIFEST.tsv`, one row per run with CPU-s. The stub demo's `PROVENANCE.tsv` is unchanged.

## Runtime
| Run | Wall | CPU |
|---|---|---|
| `w8_build_p1.py --pilot 1 --workers 2` | 716 s | 874 s |
| `w8_build_p1.py --pilot 2 --workers 2` | 798 s | 1,548 s |
| MEMOTE timing run (1 model, matrix tests included) | 2,011 s | 2,008 s |
| MEMOTE sample (4 models, final skip list) | 1,309 s | 1,295 s |
| tests (25) | 44 s | — |
| per-line materialisation check (1 line) | 34 s | 40 s |
| MEMOTE venv install | ≈ 5 min | — |
