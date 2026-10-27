---
artifact: scripts/evaluate_val.py (+ tests/test_evaluate_val.py) before freezing; W1/M2 training results (results/w1/, src/nextgem/w1.py, scripts/w1_*.py)
design_ref: DESIGN.md v3 (frozen) + §17 amendments A1-A3
reviewed_head_commit: 3f11827 (all W1 code and evaluate_val.py untracked at review time)
author: W1 analyst (Claude Opus)
reviewer: independent reviewer (Claude Opus, M2 review agent); not the author
independent_code_path: results/w1/review/recompute_w1.py (own tAP, KO join, MCC, LOLO B2, CE_fit/E[AP_random], GPR addressable set); results/w1/review/dryrun_evaluate_val.py (real-data path of evaluate_val.py on label-free inputs + a pseudo-look on 4 training lineages; no ledger, no validation labels)
verdict: evaluate_val.py ACCEPT WITH FIXES (blocking list at the end); W1 training results ACCEPT (numbers reproduced), with provenance notes
date: 2026-10-26
---

# M2 review: evaluate_val.py and the W1 training results

## 0. Leakage hygiene of this review

- **Probe.** `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **blocked** by the PreToolUse hook ("Blocked by label guard ..."), before any other action.
- I read labels only through `nextgem.labels.train()` and `labels.fold(k)`, and only in the two review scripts under `results/w1/review/`.
- I did not run `scripts/evaluate_val.py`. I did not call `labels.validation()` or `labels.test()`, and I did not touch `data/locked/`.
- The dry run imports `evaluate_val.py`'s functions. Its validation-line part is label-free: the score matrices and the B6 predictions on the real validation lines. Its "look" part runs on 4 **training** lineages standing in for validation, writes no ledger row and no output table, and prints only shapes, pass/fail and runtimes.
- I wrote nothing to the scratchpad. Review outputs are in `results/w1/review/` (two scripts and two small JSON files, all under 10 kB).
- Tests: `.venv/bin/python -m pytest -q tests/` gives **409 passed** in 134 s.

## A. `scripts/evaluate_val.py` (pre-freeze code review)

### A1. Scalars only (PI decision 21): **PASS, with a NOTE**

- **What it writes.** It writes exactly two files: `val_look_<id>.tsv` (grid and §6.8 rows) and `..._gstar.json` (`evaluate_val.py:196-200`). It prints `res`, which holds the paths, G\*, the SHAs and the wall time (`:259`).
- **Grid rows.** They come from `w1.method_summary` (`w1.py:322-350`). Every column is an aggregate over lines, genes or pairs: `tap_lin` with CI, `e_random_lin`, lift, `delta_null_lin`, line-weighted means, `prevalence_lin`, `n_lines`, `n_lines_estimable`, `n_lineages`, `n_genes_estimable`, `m_delta_ap_gene` and the pooled values. Deltas come from `w1.delta_summary` (`w1.py:353-360`): `delta_lin`, its CI and `delta_line`.
- **Binary rows.** They come from `metrics.summarise_binary`, which gives lineage-equal, line and median summaries plus counts (`metrics.py:1172-1194`). The per-line tables `lt`, `gt` and `pl` stay in memory and are never written.
- **Never written:** per-lineage values, bootstrap replicates (`BootstrapResult.rep`) and per-gene values.
- **Exception messages.** The `raise` messages in `metrics.py` (`PairingError`, `ValueError`) carry only counts or column names, for example `metrics.py:358-363` and `:381`. `assert_scalar_only` names only a column (`evaluate_val.py:178-182`). The ledger rows carry only file names and SHAs (`:191`, `:202`). The tripwire log (`logs/label_access.log`) records access events, not values.
- **NOTE (weak guard).** `assert_scalar_only` (`:174-182`) catches ModelID strings and four forbidden column names only. A future numeric per-unit column, or a per-lineage row whose key is a lineage *name* (for example in a `method` or `stratum` field), would pass it. Suggested hardening: assert that the output row count equals labels × gene sets × universes × strata × methods + the expected binary rows, and that no string cell equals a lineage name.
- **NOTE (volume).** The look reveals about 36 cells × about 20 methods × about 25 columns, plus deltas against H1 and H2. That is several thousand validation scalars, including B2, B6 and the B8 variants. DESIGN W1 ("the full grid is computed on training and validation") allows this, but every revealed scalar adds to adaptive-reuse risk for later validation decisions made against G\*.
  - B6 is provisional: it must be refit on W2's medium annotation (README caveat). Revealing it now means the later B6 needs another look.
  - Consider dropping B6, or flagging it "provisional" in the look's detail.

### A2. Ledger ordering and look cap: **PASS**

- **Ordering.** Fitting on training happens first (`run_look`, `:189`). Next, `ledger.look()` appends and **commits** `started` (`ledger.py:269-272`; commit `--only` at `:245-256`, with rollback and `LedgerError` if the commit fails). Only then is `load_val()` called (`:194`). After that come `output` (`:202`) and `completed` (`ledger.py:275`), and `completed` is written only if the body did not raise.
  - The toy test checks the ordering at the moment labels are read, by reading `HEAD:results/validation_ledger.tsv` inside `load_val` (`tests/test_evaluate_val.py:94-101`).
- **Cap.** The cap is enforced at `started` (`ledger.py:190-197`). It counts every distinct look id that has a `started` row, so crashed looks count. The W7 sub-cap is ≤ 20. The cap is raised only by a *committed* `val_cap_raise_<n>.md` (`ledger.py:159-168`).
  - Tests: `test_cap_refuses_before_labels` (cap hit before any label read), `test_repeat_look_refused_without_reading_labels`, and in `tests/test_ledger.py:80-131` the cap, the W7 sub-cap and crash/resume.
- **NOTE.** A resumed crashed look (`ledger.py:182-189`) returns the old `started` row and **appends nothing**. Re-reads of validation labels after a crash are therefore invisible in the committed ledger; only the local, git-ignored `logs/label_access.log` shows them.
  - Deterministic re-runs reveal nothing new, but the ledger should still record each re-entry.
  - Suggested fix: an `output`/`detail` row "resumed", or a `resume` event.

### A3. Label access path: **PASS**

- `labels.validation()` opens scope `val` (`labels.py:494-497`).
- The tripwire allows `val` only for the entry point `evaluate_val.py`, read-only (`tripwire.py:79-89`). The entry point must sit in `<checkout>/scripts/` of the same checkout (`:264-274`), and both `labels.py` and `evaluate_val.py` must match `HEAD:FROZEN.tsv` (`trusted`, `:367-383`).
- `labels.test()` from `evaluate_val.py` is denied: policy `{"read": {"val"}}`. This is tested in `tests/test_tripwire.py:154-157` and `:188-192`.
- Training labels come via `labels.train()` (public fold-scoped API).
- `real_context` also reads the DepMap common-essential list with `purpose="reporting"` (`evaluate_val.py:224`). Its use for the `noDepMapCE` gene set is a reporting stratum, which §5.3 and §6.6 allow.
- **Consequence for the freeze.** Until `evaluate_val.py` has a row in the **committed** FROZEN.tsv, the label read is denied. That denial would happen *after* `started` is committed (see A5).

### A4. Are these the right scalars, computed with the right code? **PASS for G\* and CE_fit; FIX for freezing the code**

- **G\* rule.** DESIGN §6.3 (line 280) says: "G\* = the better of generic H1 and H2 on validation Δ̄_lin, fixed at M2".
  - The code compares the lineage-equal mean `tap_lin` of H1 and H2 in the primary cell: L_prob, all genes, U, no stratum, validation expression lines (`evaluate_val.py:162-169`).
  - Methods are paired on identical pairs (`metrics.paired_line_delta` asserts identical `n`, `n_pos` and `estimable`; `metrics.py:468-478`), so tap_lin(H2) − tap_lin(H1) equals Δ̄_lin(H2 − H1) up to floating-point rounding. The rule is equivalent.
  - "Generic" means plain FBA e = 1 − ratio_r (B4/B5), not the task-augmented score. Correct.
  - NOTE 1: a tie goes to H2 (the default parent, §3.1). This is not pre-registered but is harmless. To match the DESIGN wording exactly, decide on the sign of `delta_lin_vs_H1` rather than comparing two separately rounded means.
  - NOTE 2: there is no one-SE or CI condition. DESIGN asks only for "the better", so a point estimate is correct.
- **CE_fit for validation.** It comes from `prepare_w1` → `w1.fit_components(train_long)` on **all training** lines (`evaluate_val.py:90-95`, `w1.py:116-130`). That matches §5.3 ("training for validation looks"). B2, B2m, B2c, B6 priors and the B6 model also use training only (`b6_full.pkl` from `labels.train()`, `w1_b6.py:132-135`).
- **Same metric code.** Validation scoring uses the same `w1.method_summary`, `w1.delta_summary`, `w1.binary_endpoints` and `metrics.*` functions as `w1_score_train.py` (`w1_score_train.py:120,129`). There is no reimplementation. The B and seed are the same: `metrics.B_OTHER` = 2,000 and SEEDS[0] = 11 (§6.5).
  - BCa is not triggered for H2 − H1 (training skew 0.149, `train_summary.json`).
- **FIX-REQUIRED: the computing code is not frozen.**
  - `FROZEN.tsv` has **no row for `src/nextgem/metrics.py`**, nor for the known-answer tests `tests/test_metrics.py` and `tests/test_metrics_binary.py`. DESIGN §11.3 says FROZEN.tsv covers `metrics.py` and the known-answer tests and tolerances. `results/w1/README.md` ("metrics.py (frozen)") is wrong on this point.
  - `w1.py` and `evaluate_val.py` are also unfrozen and untracked.
  - The tripwire checks only `labels.py` and the entry script. `config_sha` (`evaluate_val.py:73`) records the SHAs of `w1.py`, `metrics.py` and `evaluate_val.py`, but nothing *refuses* to run on a modified `w1.py` or `metrics.py`.
- **NOTE.** `data/controls/gene_panels.tsv` is an a priori list that §11.3 says belongs in FROZEN.tsv. It is an input of the look (its SHA is in `config_sha`) but is not frozen. SESOI_C1b is likewise "frozen" only in a JSON file.
- **NOTE.** The look scores validation **expression** lines only. That matches the training primary set and §4.3, but there is no "all lines" row for the generic models, while the training grid has one. This is consistent and needs no fix.

### A5. Determinism, seeds and failure modes: **PASS for determinism; FIX-REQUIRED for failure modes**

- **Deterministic.**
  - The bootstrap uses SEEDS[0] through `metrics.stratified_bootstrap_counts`.
  - B1 is seeded per (seed, line) (`w1.py:163-165`).
  - The B6 model is a fixed pickle, `predict_proba` is deterministic, and its output is rounded to 1e-8.
  - The author's full re-run reproduced every training output bitwise (`rerun_sha_check.tsv`), and I checked that the MANIFEST output SHAs equal the current files.
- **Crash mid-look burns the look.** A crash after `started` leaves the look open. It can resume only with the same `config_sha` (`ledger.py:186-189`). That hash covers `evaluate_val.py`, `w1.py`, `metrics.py` and every input SHA. The look id is hard-wired (`argparse choices=["W1_M2_baseline"]`, `:242`).
  - So any crash that needs a **code or input change** to fix makes `W1_M2_baseline` unrecoverable. Recovery would need an amended `evaluate_val.py`, a new look id and a second look against the cap.
  - Several label-free failure points currently sit **inside** the look, after `started`:
    - `ctx["b6"](lines, fitted)` (B6 feature matrix and `predict_proba`), called in `w1_scalars` (`:106-107`);
    - `w1.score_matrices`, the `stratum_specs` application and the output writes;
    - the tripwire trust check: if `evaluate_val.py` is not yet in **committed** FROZEN.tsv, `check_frozen.py` still passes (it checks only listed files), `started` is committed, and then `labels.validation()` is denied.
  - The toy tests never exercise the real-data path: `ctx["b6"] = None` and toy strata (`tests/test_evaluate_val.py:84-87`).
- **NOTE.** `check_frozen.py` is run with `--no-test-runs` (`:250`). That skips the `results/test_runs.tsv` check, which DESIGN §11.3 uses to catch disappeared or downgraded known-answer tests. The look should run the full check.
- **NOTE.** There is no provenance wrapper. DESIGN §11.4 requires `nextgem.run.provenance()`, which does not exist. Without it, nothing refuses a dirty `src/`, `scripts/` or `tests/`, and no MANIFEST row is written for the look.
  - The ledger records `head_commit` and `config_sha`, but the config JSON (input and code SHAs) is never written, so the hash cannot be audited.
  - Suggested fix: write `val_look_<id>_config.json` (label-free) before `started`, and refuse to run if `git status --porcelain src scripts tests` is non-empty.
- **Reviewer dry run** (`results/w1/review/dryrun_evaluate_val.py`; results in `dryrun_evaluate_val.json`):
  - `config_sha` = `0325ec02…a7acc104`, matching the README.
  - **Label-free part, real validation lines.** 200 validation lines, 198 of them with expression. All 25 score matrices, including the B6 full model, are finite on the 198 lines. No line has a medium class missing from the B6 class codes. Runtime: 7 s.
  - **Pseudo-look on 4 training lineages** (184 lines; fit on the remaining training lines). `w1_scalars` completes, with 856 rows × 47 columns in the sections `tAP_grid` and `binary_6.8`. The scalar guard passes, and G\* has the expected keys. Runtime: 196 s.
  - The real path therefore runs end to end on real-shaped inputs. Moving the label-free steps before `started` is still required, so that an input change (for example a regenerated `tmp/`) cannot strand the look.

## B. W1/M2 training results

### B1. Independent recomputation: **PASS (reproduced exactly)**

- **AP implementation.** My tie-aware AP is derived from linearity: precision at a positive's in-block rank is linear in the number of block positives above it.
  - On 300 random small tied cases, I checked it against **brute-force enumeration** of every distinct block permutation: max |diff| 2.2e-16.
  - `metrics.tap` gives the same 2.2e-16 against brute force, and 2.2e-16 against mine on 200 cases of size 50 to 3,000.
  - A constant scorer equals my closed-form E[AP_random] exactly.
- **Own pipeline.** I used my own Entrez join of `results/baseline/*_ko.parquet`, with e = 1 − clip(round(ratio_r, 6), 0, 1) and absent genes set to 0. Labels came via `labels.train()`. The 686 lines are the eligible training lines (not `rare`) with `has_expression` from the split file. That set is identical to the author's `has_expression_matrix` set.

| Quantity (training, L_prob, U, 686 lines, ≥ 5 positives) | Reviewer | Author (`train_grid.tsv` / `train_binary.tsv`) |
|---|---|---|
| H1 tAP_line, lineage-equal | 0.2995756 | 0.299576 |
| H2 tAP_line, lineage-equal | 0.2916272 | 0.291627 |
| Δ̄_lin(H2 − H1) | −0.0079483 (negative in 12/12 lineages) | −0.0079 [−0.0083, −0.0075]; 12/12 |
| MCC P57 × L-ge06, parent genes, H1 (lineage-equal / line mean) | 0.328646 / 0.329957 | 0.328646 / 0.329957 |
| MCC P57 × L-ge06, parent genes, H2 | 0.331475 / 0.332747 | 0.331475 / 0.332747 |
| MCC difference H2 − H1 (lineage-equal) | +0.00283 | +0.0028 |
| Lines with NaN MCC | 0 / 0 | 0 / 0 |
| B2 (L_prob), LOLO via `labels.fold(k)` | 0.90013 | 0.900 |

- The parent gene sets are U ∩ H1 = 2,825 and U ∩ H2 = 2,807 genes.
- **Conclusion:** H1 > H2 on training tAP is real, with a small effect. H1.2 is not supported on training. The paper protocol (P57) gives a negligible reversal. G\* is still decided on validation only, which is correct.

### B2. Baselines B2, B2m, B2c and B6 (tAP ≈ 0.90) are fold-scoped: **PASS**

- **B2, B2m, B2c and CE_fit.**
  - They are fitted in `w1_score_train.py:60-80` through `cv.lolo`. `cv.lolo` hands the fit function only `labels.fold(k)`, the training labels minus lineage k, and asserts that no held-out ID is in the fitting set (`cv.py:36-44`).
  - Scores for lineage k use the components fitted without k (`w1.score_matrices`, `w1.py:153-191`).
  - B2m's class membership (`base_medium_class`) is label-free.
  - The `noCE` and selective masks are also per-fold (`:73-78`, `:100-103`).
  - My own LOLO B2, built directly from `labels.fold(k)`, gives 0.90013. The high value is therefore not a leak: gene identity (common essentials) dominates, as the literature summary in RESEARCH_TREE anticipated.
- **B6 (`w1_b6.py`).**
  - It is LOLO through `cv.lolo` (`:127`).
  - Inside a fold, the B2 and B2m features of fitting rows are *stacked* from the fitting set minus the row's own lineage (`stacked_priors`, `:37-75`).
  - Held-out rows get priors from fitting lines only (`w1.heldout_priors`).
  - The other features (expression, CN, mutations, neighbours, paralogs, isozymes, medium class, serum) are label-free.
  - I found no use of a scored line's own labels, or of its lineage's labels.
- The author's fold-scoped memoriser control returns exactly E[AP_random] on every fold (`verification.tsv` C5), which corroborates this.

### B3. SESOI_C1b = 0.005723 and the M2 gate: **PASS, with a NOTE**

- **Recomputation.** CE_fit is the set of genes with mean L_prob ≥ 0.9 over all 767 training lines: 167 genes. With E[AP_random] in closed form on lines with ≥ 5 positives, I get:
  - mean E[AP_random] on U = 0.1235308;
  - mean E[AP_random] on U minus CE_fit = 0.0706943;
  - SESOI = **0.0057228**, identical to `sesoi_c1b.json`.
- **NOTE (ambiguous averaging).** "Mean" is not specified in §6.4. The author used the line mean over all 767 training lines, including the rare and no-expression lines. The alternatives give nearly the same value:
  - line mean over the 737 eligible lines: 0.0057222;
  - lineage-equal mean over the 686 expression lines: 0.0057170.
  The choice is immaterial (< 0.1%). It should still be recorded as the implementer's choice and frozen in FROZEN.tsv as a `value` row before the look.
- **NOTE (timing).** "Before any model is scored": SESOI is computed in the same run as the scoring (`w1_score_train.py:209-226`), but it depends only on line sizes and positive counts, not on any model. The value cannot have been influenced by model scores.
- **Gate.** The H2 lift on training is 2.359 (lineage-equal) and 2.374 (pooled), both > 1, so the gate passes. My recomputed H2 tAP agrees, and its lift is 0.29163 / 0.12361 = 2.36.

### B4. P57 vs Human-GEM `taskEssentialGenes`, and the addressable set: **PASS**

- **P57.** `w1_verify.py --task-equality` (`:250-273`) runs Human-GEM's own `find_task_essential_categories` (`code/test/taskEssentialGenes.py`) on the same loaded models with Gurobi. The result is equal gene sets: 226/226 for H1 and 229/229 for H2 (`verification_tasks.tsv`).
  - This validates the reimplementation against the repository's Python code. It does not validate it against the papers' MATLAB/RAVEN `checkTasks` (NOTE, acceptable).
  - The task files match Robinson S3 on all 57 (id, description) pairs. IN/OUT sets match for 55/57; the README attributes the other two to metabolite renamings, which I did not check independently.
- **Addressable set.** I recomputed it independently: load each parent with cobra, and call a gene addressable if knocking it out makes ≥ 1 of its reactions' GPR false.
  - Within U: 1,064 genes in H1 and 1,098 in H2, identical to `gene_features.tsv`. DESIGN's 1,120 for H2 counts all model genes, not only U.
  - There are 235,187 positive L_prob pairs on the 686 lines. The fraction outside the addressable set is **0.4232 (H2)** and **0.4513 (H1)**, identical to `train_addressable.tsv`.
  - The definition matches §5.4: the share of positives outside the set is the recall ceiling's complement.
  - This is plausible: about 62% of U genes (1,768 of 2,866 for H2) are not addressable, either because they are absent from the parent or because an isozyme covers their KO.

### B5. Label copies in committed outputs: **PASS**

- `scripts/scan_label_copies.py`, run as a milestone scan with fingerprints loaded, reports **0 hits** over 668 files.
- I checked ModelIDs in the non-`tmp/` files of `results/w1`:
  - `train_per_line.tsv` (2.6 MB) and `verification.tsv` hold **training** IDs only: 686 and 5.
  - `line_features.tsv` and `background_members.tsv` hold all 1,361 IDs, but their columns are label-free: split covariates, medium, serum, expression flag, and genomic group membership as 0/1.
  - The test-keyed pre-commit rule (`scan.test_keyed`) does not trigger on any of them.
- **NOTE.**
  - `train_per_gene.tsv` (training per-gene prevalence, n_pos and `frac_pos_train`) and `train_per_line.tsv` (training per-line n_pos, tAP) are training-tier label aggregates. DESIGN W1 (the per-gene table, ext-PI SF1) and §6.6 require them, so they are fine to commit. Readers are exposed to training aggregates, which is non-disqualifying under the register rule.
  - `train_per_line.tsv` is the largest committed file at 2.6 MB. If compactness matters, it could be Parquet under `tmp/` plus a summary.

### B6. MANIFEST.tsv and runtime.tsv: **NOTE (FIX before M2 sign-off, not blocking for the look)**

- **What MANIFEST.tsv records:** timestamp, script, argv, commit, a dirty flag, output SHAs, seeds, solver name, workers, wall time, CPU time, peak RSS and author.
- I checked the output SHAs against the current files (`train_grid`, `per_line`, `per_gene`, `sesoi`, `binary`, both `p57` files, `verification`, `b6_full.pkl`, `gene_panels.tsv`): all match.
- **Gaps against §11.4 and §6.7:**
  - (a) **no input SHAs** (the KO parquets, `labels_train` hash, expression, CN and mutation sources);
  - (b) output SHAs are truncated to 16 hex characters;
  - (c) the `commit` field is 3f11827 while all W1 code is untracked (`dirty=True`), so the code that produced the outputs is not recoverable from git. Re-run, or record per-script SHAs, after the code commit;
  - (d) no solver version (Gurobi 13.x, cobra/optlang);
  - (e) several outputs have no row: `train_per_lineage.tsv`, `train_addressable.tsv`, `train_summary.json`, `b6_run.json`, `p57_*.json`, `inputs_summary.json`, `runtime.tsv`, `rerun_sha_check.tsv`;
  - (f) there is no per-line runtime (§6.7: "per line and in total").
- `runtime.tsv` has wall time and a Gurobi flag only; CPU time and RSS are in MANIFEST.
- The README's runtime table is adequate for reporting.

## Other observations (non-blocking)

- B6 < B2 on training (0.888 vs 0.900): the §12 W7 positive control fails at the M2 defaults. The author flagged it correctly and did not tune B6 (tuning is W7's).
- B2's mΔAP_gene of −0.033 is explained as a LOLO artefact. The explanation is plausible: the out-of-fold prior is anti-correlated with fold prevalence. I did not verify it.
- The failure-mode traces (TYMS, β-oxidation, CYP4F2/BBOX1/TMLHE) and the v1.3.0 sensitivity are not done. They are W1 items still open at M2.

## Blocking issues before freezing evaluate_val.py

1. **Freeze the code that computes the look, and commit it before the look.** Add FROZEN.tsv rows for:
   - `scripts/evaluate_val.py` (the tripwire requires it in **committed** FROZEN.tsv);
   - `src/nextgem/w1.py`;
   - **`src/nextgem/metrics.py`**, which §11.3 requires and which is currently missing although the README says it is frozen;
   - the known-answer tests `tests/test_metrics.py`, `tests/test_metrics_binary.py` and `tests/test_evaluate_val.py`;
   - `data/controls/gene_panels.tsv` and the SESOI_C1b value.

   Commit the W1 code, so that `head_commit` in the ledger identifies the code that ran.
2. **Move every label-free step before the `started` row**, and add a pre-flight check. A crash after `started` that needs a code or input fix permanently burns the hard-wired look id (A5).
   - Compute `real_context`, the B6 predictions and `w1.score_matrices` for the validation lines in `run_look` *before* `ledger.look(...)`, and pass `S` into `w1_scalars`.
   - Pre-flight, before `started`, that `tripwire.trusted()` would pass: `evaluate_val.py` and `labels.py` match HEAD:FROZEN.tsv, and so do `w1.py` and `metrics.py`.
   - Add a toy test with `ctx["b6"]` set (a stub callable), so the B6 branch of `w1_scalars` is covered.
3. **Run `check_frozen.py` without `--no-test-runs`** at the start of the look (§11.3: the evaluate scripts check that known-answer tests have not disappeared or been downgraded). Otherwise, record in DESIGN §17 why the test-runs check is skipped.

Should-fix (not blocking):
- Write the config JSON (label-free) next to the output before `started`.
- Refuse to run on a dirty `src/`, `scripts/` or `tests/` (the §11.4 provenance wrapper is missing).
- Log resumed looks in the ledger.
- Harden `assert_scalar_only` (expected row count; no lineage-name cells).
- Decide G\* on the sign of `delta_lin_vs_H1`.
- Consider leaving the provisional B6 out of the W1 look.
- Complete MANIFEST provenance (B6).

## Reviewer label exposure (for notes/label_exposure.md)

2026-10-26, M2/W1 reviewer:
- **What I computed and viewed:** training-tier aggregates only, via `labels.train()` and `labels.fold(k)`:
  - lineage-equal tAP for H1, H2 and LOLO B2;
  - Δ̄_lin(H2 − H1) and its sign per lineage (12/12 negative; per-lineage values were not printed);
  - P57 × L-ge06 MCC means;
  - CE_fit size (167) and mean E[AP_random] on U and on U minus CE_fit;
  - the count of positive pairs (235,187) and the fraction outside the addressable set.
- **What I read:** the author's training outputs: the README tables, including per-lineage ranges and per-stratum training values, `train_summary.json`, the primary rows of `train_grid` and `train_binary`, and the headers of `train_per_line`, `train_per_gene` and `train_per_lineage`.
- **What I did not see:**
  - no gene × line label values;
  - no validation or test labels or summaries (the dry run's validation part is label-free; its pseudo-look is on training lineages and printed no metric values);
  - no raw files (the probe was blocked).

---

## Re-check of the author's fixes (2026-10-27, same reviewer)

- **Probe.** The raw-label probe (`head -c 50` on the raw CRISPR gene-effect file) was **blocked** by the hook.
- **What I did not do:** run the look, call `labels.validation()` or `labels.test()`, touch `data/locked/`, edit FROZEN.tsv, DESIGN or the ledgers, or commit.
- **What I ran:**
  - `--print-config` and `--print-frozen-candidates`, both label-free;
  - the full test suite;
  - `results/w1/review/recompute_w1.py`, re-run;
  - a new label-free dry run of the revised script, `results/w1/review/dryrun_evaluate_val_v2.py` → `dryrun_evaluate_val_v2.json`.
- **Tests:** `pytest -q tests/` gives **435 passed** (219 s). The author reported 430; the extra tests come from other new test files (`test_gurobi_pool.py`, `test_run.py`).

### Verdict on each claimed fix

| Claim | Verdict | Evidence |
|---|---|---|
| Pre-flight before anything is written: committed HEAD:FROZEN.tsv rows, a clean git state, and check_frozen with test runs | **PASS** | See below this table. |
| Label-free work and the config JSON run before `started` | **PASS** | See below this table. |
| B6 excluded, both branches toy-tested | **PASS** | See below this table. |
| A resumed look appends an `output` row with `{"event": "resumed"}` before labels are re-read, and is counted once | **PASS** | See below this table. |
| `assert_scalar_only` hardened | **PASS** | See below this table. |
| G\* decided on the sign of Δ̄_lin(H2 − H1), tie to H2 | **PASS, with a NOTE** | See below this table. |
| `run.py` `provenance()`, MANIFEST regenerated | **PASS** | See (e) below. |
| Training outputs bitwise unchanged | **PASS** | See (d) below. |

Evidence for each row:

- **Pre-flight.** `main` (`evaluate_val.py:418-422`) runs `preflight(required_files(cfg))` and `check_frozen_full()` before anything is written.
  - `check_frozen_full` (`:188-191`) runs `check_frozen.py` with no `--no-test-runs`.
  - `preflight` (`:157-185`) requires, for every required file, a `file` row in `git show HEAD:FROZEN.tsv` whose SHA equals the current file. Files that git does not ignore must also be tracked and unmodified; ignored inputs are checked by SHA only.
  - `provenance(allow_dirty=False)` (`run.py:131-132`) then refuses a dirty `src/`, `scripts/` or `tests/`, and untracked files count as dirty. All of this happens before `started`.
  - The toy test is `test_preflight_frozen_and_clean`.
- **Label-free work first.** Inside `provenance`, the order is `real_context` → `prepare_look` (fit on `labels.train()`, then every validation-line score matrix with a finiteness check) → `write_config` → `run_look` (`:433-436`).
  - `run_look` refuses to run if the config JSON is missing or stale (`:338-340`).
  - Tests: `test_config_json_required_before_started` and `test_label_free_failure_happens_before_started`. In both, the ledger stays empty and there are 0 label reads.
- **B6.** `INCLUDE_B6 = False` (`:56`), so the B6 input files drop out of the config.
  - Tests: `test_b6_excluded_by_default` (stub never called) and `test_b6_branch_when_enabled` (12 B6 rows), plus the label-free failure test above.
- **Resume marker.** The marker is appended inside `ledger.look` after `started` and before `load_val()` (`:346-351`). It is committed through `ledger.append`, so a failed commit aborts before any read. Test: `test_resumed_look_is_logged_and_counted_once`.
- **Scalar guard** (`:281-321`). It now checks:
  - row caps per section;
  - one row per scalar family;
  - no column name or string cell that is a line ID, a lineage name, an Entrez ID or a gene symbol;
  - a whitelist for numeric columns.

  Six toy tests cover it. On the real context (pseudo-look) it passes with 720 + 124 rows. The binary cap 2 × 5 × 6 × 2 + 4 = 124 is exact. String cells are fixed by the configuration, not by the data, so a pass here predicts a pass on the real look.
- **G\*.** `decide_gstar` (`:270-278`) applies the sign rule, and `test_gstar_sign_rule` covers −, 0 and +. NOTE: a NaN Δ̄_lin would silently give Human2 (`d < 0` is False). Raising on a non-finite `d` would be better. In practice it cannot be NaN here: the dry run gives 198 validation expression lines, and every training lineage has estimable lines.

### (a) Can any validation label value, or anything derived from one, be computed or printed before `started`? **No path found.**

Everything that runs before `ledger.look` touches only label-free files or training labels:
- `w1_config`: SHAs of label-free inputs and code;
- `frozen_candidates` and `preflight`: SHAs and git only;
- `check_frozen_full`: hashes protected files through `labels.sha256_protected`, and prints at most 2,000 characters of problem lines, which hold SHAs and paths and never values;
- `real_context`: split, line and gene features, omics, the gene universe from the baseline mappings, gene symbols, and the DepMap common-essential list (`purpose="reporting"`), which is an existing DepMap all-line product and not a statistic we compute on validation;
- `prepare_look`: `labels.train()` only;
- `provenance` entry: a git status check.

`labels.validation()` is called in exactly one place (`:351`), after `started`, and after the resume marker if there is one. The tripwire would refuse any other caller anyway (val read only for the frozen `evaluate_val.py`; test denied).

After `started`, exception messages carry column names or counts but no values. The guard's messages name the column and the *kind* of offending value ("holds a lineage name value"), not the value itself.

### (b) Does the candidate list cover every file whose change could alter the look's numbers? **Not quite: 3 gaps.**

These are already covered:

| File | How it is covered |
|---|---|
| `labels.py` | already frozen; in `CODE_FILES` |
| `config.py`, `paths.py`, `ledger.py`, `tripwire.py`, `frozen.py`, `cv.py`, `cohort.py`, `splits.py` | frozen; enforced by `check_frozen_full`. `cv`, `cohort` and `splits` are not imported by the look. |
| `labels_train.parquet`, `labels_val.parquet` | `label_file` rows, hashed by `check_frozen` |
| `results/splits/lineage_split.tsv` | tracked and clean; `labels.load_split()` refuses it unless it matches the committed `SPLIT_SHA256`. Adding it to the required list would be belt-and-braces. |
| `mut_dam`, `mut_hot`, `gene_neighbours`, `b6_full` | read or loaded, but cannot affect any number while B6 is off |

These are missing from `required_files`, the config and the candidate list:

1. **`results/baseline/human1_gene_mapping.tsv` and `results/baseline/human2_gene_mapping.tsv`.** They define U and U∩ through `labels.gene_universe()` (`labels.py:127-141`), so they enter **every** number, including G\*. They are tracked but not frozen, and a working-tree edit would pass the pre-flight. **Blocking.**
2. **`data/labels/label_derived_lists.parquet`** (`depmap_ce`, used by the `noDepMapCE` cells). It is unfrozen and not in the config. It is a protected file, so the pre-flight's plain `_sha` cannot hash it. It needs a `label_file` FROZEN row, which `check_frozen` hashes through `labels.sha256_protected`. It does not affect G\*. **Blocking** for the reported grid, but a one-row fix.
3. **`src/nextgem/__init__.py`.** It is executed on every `nextgem` import and is not frozen. The risk is low (a lazy loader). Add a row.

Notes:
- `KNOWN_ANSWER_TESTS` appear in the candidate list but not in `required_files`, so the pre-flight does not insist on them. Once the PI adds them, `check_frozen_full` enforces them.
- Library versions (numpy 2.4.6, pandas 2.3.3, scipy) are not pinned: `requirements.txt` is modified in the working tree. They are recorded after the fact in the MANIFEST `versions` column. Acceptable, given the bitwise-reproducibility check.

### (c) Can the "resumed" marker be used to evade the look cap? **No.**

- `ledger._validate` (`ledger.py:198-204`) accepts an `output` row only for a look id that already has a `started` row and is not completed.
- The cap counts **distinct look ids** with a `started` row (`ledger.py:143-156, 190-197`). A marker therefore cannot create, remove or rename a look, and cannot reopen a completed one.
- Resuming requires the identical `config_sha`, which covers the code and input SHAs (`ledger.py:186-189`). A resume can only re-run the identical, deterministic computation, and every re-entry is now a committed marker row.
- The toy test checks this end to end: 39 other looks plus a crashed look give 40; the resume is allowed and not double-counted; a new look id is then refused with "cap reached".

### (d) Training outputs vs my `recompute_w1.json`: **unchanged.**

The SHA-256 values of `train_grid`, `train_per_line`, `train_per_gene`, `train_per_lineage`, `sesoi_c1b`, `train_binary`, `train_binary_paired`, `train_addressable`, both `p57` parquets and `tmp/b6_oof` equal the values I checked at the first review. The author's `rerun_sha_check.tsv` agrees; only `verification_tasks.tsv`, which holds elapsed seconds, differs.

My re-run of `recompute_w1.py` reproduces every figure:

| Quantity | Value |
|---|---|
| H1 tAP_lin | 0.2995756 |
| H2 tAP_lin | 0.2916272 |
| Δ̄_lin(H2 − H1) | −0.0079483 |
| MCC, H1 | 0.328646 |
| MCC, H2 | 0.331475 |
| B2 LOLO | 0.900133 |
| SESOI_C1b | 0.0057228 |
| Fraction of positives outside the H2 addressable set | 0.42325 |

### (e) Regenerated MANIFEST: **adequate, with one remaining NOTE.**

- All 11 rows have the 17 `run.py` columns.
- Every row has full SHA-256 values for inputs and outputs.
- Every row records seeds, workers, wall time, CPU time (including children), peak RSS, the solver and its version (Gurobi 13.0.3, GLPK 5.0), and package versions (python 3.11.17, cobra 0.32.1, optlang 1.9.1, numpy 2.4.6, pandas 2.3.3, sklearn 1.9.1).
- Label inputs are recorded through `extra.label_files_frozen_sha256`.
- Every file in `results/w1/` except `README.md` and `MANIFEST.tsv` itself is an output of some row.
- **NOTE:** every row still says `commit = 3f11827, dirty = True`, because the W1 code is untracked. Once the code is committed, re-run `scripts/w1_run_all.py` (about 41 min, bitwise-deterministic) so that the training outputs carry a clean commit. This is not needed for the look, whose inputs are pinned by SHA in FROZEN.tsv and in its config.

### Dry run of the revised script (label-free; `dryrun_evaluate_val_v2.json`)

- `config_sha` = `331a2c8d…a60e5`, with `include_b6 = false`.
- **Real validation lines:** 200 lines, 198 with expression. `prepare_look` builds all 24 score matrices (B6 excluded) with finite scores in 2.4 s. This is exactly the pre-`started` work of the real look.
- **Pseudo-look on 4 training lineages** (184 lines): `w1_scalars` completes in 189 s. It produces 720 grid rows and 124 binary rows × 47 columns. The hardened guard passes, and G\* is consistent with the sign rule.

### Re-check decision: **GO, conditional on PI freeze actions only (no further code change required)**

Remaining blockers. All of them are FROZEN.tsv or commit actions by the PI, and all are enforced automatically by the script's pre-flight, `check_frozen_full` or `provenance` **before** `started`:

1. Add FROZEN.tsv rows for everything `--print-frozen-candidates` lists: `evaluate_val.py`, `w1.py`, `metrics.py`, `run.py`, the 10 inputs, the 4 known-answer tests, and `sesoi.C1b`. **Also add rows the list misses:**
   - `results/baseline/human1_gene_mapping.tsv` and `results/baseline/human2_gene_mapping.tsv` (they define U);
   - `labels/label_derived_lists.parquet` as a **`label_file`** row (the DepMap CE gene set);
   - `src/nextgem/__init__.py`.

   Optionally, also add the two mapping files to `input_files()` so that `config_sha` captures them. That is a one-line change and must be made before `evaluate_val.py` is frozen.
2. Commit FROZEN.tsv, all W1 code and the tracked look inputs (`results/w1/*` and `results/baseline/*`). Also commit, or remove, every other untracked file under `src/`, `scripts/` and `tests/` (W2, W4, context, gurobi_pool); otherwise `provenance(allow_dirty=False)` refuses to run. Commit `results/test_runs.tsv` so that `check_frozen` passes with test runs.
3. Run `scripts/evaluate_val.py --look W1_M2_baseline --preflight-only`. It must print "pre-flight OK" before the real look is run, once.

Should-fix (not blocking): make `decide_gstar` raise on a non-finite Δ̄_lin; after the code commit, re-run `w1_run_all.py` for clean provenance.

Reviewer exposure in this re-check: the same training-tier aggregates as before (re-computed values listed in (d)). The dry run printed no metric values. No validation or test labels or summaries; no raw files.
