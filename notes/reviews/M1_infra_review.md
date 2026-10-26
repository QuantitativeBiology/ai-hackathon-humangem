---
artifact: M1 infrastructure (labels, splits, leakage guards). Uncommitted working tree on top of HEAD 0c06695: src/nextgem/{labels,splits,cohort,cv,ledger,frozen,scan,tripwire,paths,config}.py, scripts/{build_labels,make_splits,scan_label_copies,check_frozen,verify_ledgers,pre_commit_guard,install_tripwire,install_git_hooks,git_locked}.py, tests/ (infra files), FROZEN.tsv, .venv sitecustomize.py, .git/hooks/pre-commit and pre-push, results/m1_infra/README.md, plus the diffs to scripts/design_numbers.py, scripts/checksum_data.py, src/nextgem/__init__.py, .gitignore, requirements.txt
reviewed_state: HEAD 0c0669529f42; FROZEN.tsv sha256 49ff7f8e22ef83c7...; check_frozen.py reports "FROZEN.tsv OK (30 rows)", so every reviewed src/scripts/tests file matches its FROZEN.tsv SHA-256; sitecustomize 77c0246008855249...; pre-commit 2b57258c73e4c5ae...; pre-push 5aa44604ae70c682...; README c464bc9378bbc9da...
design_ref: DESIGN.md FROZEN v3 (61b0beb): sections 3, 4/4.1, 5, 11; section 12 rows "Splits", "Agent guards", "Labels"
author: M1 infra builder
reviewer: independent reviewer (Claude agent, Opus; did not write the infrastructure)
date: 2026-10-26
verdict: ACCEPT WITH FIXES
independent_code_path: yes
label_access: none. The liveness probes (head of the raw gene-effect file; Read of data/locked/HOOK_PROBE) were both denied. Every bypass probe ran against a toy label store built by tests/toy.py in the scratchpad. No draw was evaluated for any rng seed in 38..137.
files_edited: none (no commit). Scratch scripts are in the session scratchpad rev/ folder (indep_split.py, bypass.py, scanprobe.py, hookcases.py).
---

# M1 infrastructure review

## Verdict

**ACCEPT WITH FIXES.**
- The label definitions, the NaN accounting, fold scoping, the tier files and the split algorithm match DESIGN v3.
- An independent implementation of the draw, written from the DESIGN text, gives identical per-line tiers for 20 simulated seeds.
- The full suite passes.

Six must-fixes (M1–M6) remain. All are cheap. They must land before the production draw or before `build_labels.py` runs, whichever applies. They concern:
- production overrides in `make_splits.py`;
- sandbox access to the production seed;
- split inputs that are not yet pinned;
- trust in tripwire entry points;
- the amendment-id parser;
- a scanner oracle.

## 1. Fidelity to DESIGN

| Item | Finding |
|---|---|
| Labels (§5.1, §6.8) | Correct, including the strict boundaries: L_prob = dep > 0.5; L_ge05 < −0.5; L_ge06 < −0.6; L_ge10 < −1. A NaN source gives an NA label (nullable Int8). |
| κ rule | κ is computed once, on training lines pooled over U pairs. It is logged as a `non_look` row in the validation ledger (committed), and `label.primary` goes to FROZEN.tsv in the same run. κ without CE_fit is informational only. |
| Join | Entrez only. Build refuses on unequal ID sets, duplicate Entrez IDs, U genes without a CRISPR column, and label rows that are not in the split. |
| NaN accounting | `dropcounts.tsv` has per-gene, per library-class × tier and per-tier counts, plus the lines missing from the labels. |
| Tier files | `labels_train.parquet` in data/labels; validation and test files in data/locked. The hashes go to FROZEN.tsv. A rebuild must reproduce them and overwrites nothing. |
| Fold scoping | `fold(k)` and nested `fold(k, inner)` exclude the held-out lineage. Rare lines are never a fold. `cv.lolo` asserts that there is no overlap. Scoping is a convention, which is acceptable: training labels are public through `train()`, and `score()` sees the held-out labels. |
| Sealed test | The tripwire allows test reads only to `evaluate_test.py`. That script also has validation reads, which the refit on training ∪ validation (§4.2) needs. |
| Split (§4.1 steps 1–6) | Eligibility ≥ 15; ranking by size descending, then alphabetical; strata 8/8/7; permutation order S1, S2, S3; test 2/2/3 and validation 1/1/2. Label-blind checks (a)–(d). `default_rng(S + attempt)`, attempts 0–99. Retired-range refusal for S in 20260927…20261125, repeated per attempt. Refuses if `SPLIT_SHA256` exists, and the ledger refuses a second `split_started`. The attempt log keeps every rejected draw. The covariates match step 6. |
| Cross-checks | 20 reference draws; acceptance statistics for 1–1000; equality with the verbatim AST-extracted `draw`/`checks` of design_numbers.py at 61b0beb, for 1–1000 and for S…S+99. |
| design_numbers.py refactor | It imports `splits.draw`/`checks` and asserts that the constants are equal. Behaviour is unchanged (the README reports identical JSON). Not done: the §11.1 note "refactored to use `labels.ids()`". It is moot, because the script opens no label file. |
| checksum_data.py | Skips data/locked and data/labels. This is fine: their hashes live in FROZEN.tsv. |

**Deviations I accept:**
- The fingerprints are per cell, not per row (README item 2). But see M6.
- Per-tier NaN counts include the validation and test tiers. This is missingness only, and DESIGN asks for the counts.

## 2. Independent split implementation (20 simulated seeds)

My script reads CRISPRScreenMap, Model.csv and the default expression rows itself. It takes per-lineage proteomics coverage from the frozen JSON (`coverage_by_lineage`). From the DESIGN text it implements:
- the strata;
- the draw;
- checks (a)–(d);
- the attempt loop;
- the per-line tier assignment.

| Check | Result |
|---|---|
| Lineage counts, per-lineage expression and strata equal the frozen JSON | true |
| My reference draws for seeds 1–20 equal the JSON | true |
| The builder's per-lineage expression and proteomics equal mine | true / true |
| Per-line tiers and attempt counts equal the builder's (`cohort.load` + `splits.run_attempts` + the make_splits assignment, in memory), seeds 900000…900962 step 37 | **20 / 20 identical** (12 accepted at attempt 0, 8 needing one redraw) |

Nothing was written to disk, and no tiers were printed. Seed 38 and the range 38–137 were never evaluated.

## 3. The proposed §4.1 amendment (`--code-commit C`)

**It is sound and necessary for the split itself.** The literal V-1 check cannot pass, because 61b0beb contains no split code.

**What the amendment gets right.** The production tiers become a function of pre-reveal artefacts only. These are:
- the strata from the 61b0beb JSON;
- the verbatim `draw`/`checks` of 61b0beb, for every attempt seed;
- the code and FROZEN state at C, with C..HEAD clean;
- DESIGN.md, which may only gain lines.

**Four surfaces are still not pinned, so it is not yet airtight.** Each could move the tiers after S became known:
1. `expr_by_lin`/`prot_by_lin` come from the new `cohort.py`. They are checked only through aggregate acceptance statistics. `harmonized_MS_CCLE_Gygi.csv` is read "if it exists", and it is not in the freeze `data/checksums.tsv`.
2. `run_attempts`, the loop that decides which attempt is accepted, is new code and is not cross-checked.
3. The line → lineage map is checked only by its counts.
4. The production run accepts overrides that bypass every pin:
   - `--reference-commit` (any commit becomes "the frozen reference");
   - `--repo` (git checks run on another clone while the code runs from this one);
   - `NEXTGEM_DATA_ROOT`/`NEXTGEM_STATE_ROOT`.

These are closed by M1 and M3.

**What it cannot do.** It cannot restore commit-reveal protection for anything written after 5a8c218. Since the reveal, anyone with the repository can compute the tiers from 61b0beb's design_numbers.py and S = 38. Everything authored in that window could in principle have been tuned to known tier names:
- this infrastructure;
- the metrics library;
- the D5 blind re-derivation of the a priori lists, which is in progress;
- any amendment.

The label guards limit this to tier *names*, not label values. Lineage-specific a priori content (background groups, D5) is the exposed surface.

**Proposed tightening:**
- (a) Close M1 and M3.
- (b) Restrict post-freeze DESIGN.md additions to hunks inside the §17 amendment table. "Only gained lines" allows any added sentence.
- (c) Draw immediately after C, to shrink the window.
- (d) Every post-reveal author, including the D5 list author, attests in the commit or review that they never evaluated `draw()` on seeds 38–137 and were not told the tiers. I attest this for myself.
- (e) The §17 row logs both the `--code-commit` rule and the fact that S = 38 was stated by the user *before* the freeze (RESEARCH_TREE D4). §4.1 step 2 requires "after", so this is a deviation to record. It does not invalidate the design, because the design authors did not know S.

## 4. Guard bypass probes (toy store; "works" = toy test-tier rows obtained)

**Tripwire**

| Result | Probes |
|---|---|
| Blocked | `open`/pandas/`pq.read_table`; `python -I`; `np.fromfile`; `np.memmap`; symlink made by shell `ln -s`; multiprocessing spawn child; `np.load` of fingerprints.npz |
| **Works** | `python -S` with site-packages appended to `sys.path`; system `python3` (documented residual) |
| **Works** | A `sitecustomize.py` on `PYTHONPATH` shadows the tripwire |
| **Works** | `pyarrow.fs.LocalFileSystem().open_input_file`; `pyarrow.OSFile`; ctypes libc `open` |
| **Works** | `subprocess` shell glob (`cat …/lock*/labels_te*`): the argument regex misses globs |
| **Works** | Hard link made by shell `ln`: inode protection covers raw files only, not tier or locked files |
| **Works** | `cp` with a glob in the shell, then a normal read of the copy |
| **Works** | Spoofing `sys.modules['__main__'].__file__ = …/scripts/evaluate_test.py` before `labels.test()` |
| **Works** | **Copying `src/` to any directory and adding `scripts/evaluate_test.py` next to it.** The labels.py check uses only the last three path components, and the entry point is matched by name. Writing `scripts/evaluate_test.py` in the project works the same way, because the file does not exist and is not frozen yet. This is the most plausible "make it work" bypass. |
| **Works** | Setting `nextgem_tripwire._state.busy = True`, or monkeypatching `nextgem_tripwire.check` (deliberate) |

**Content scanner** (with toy fingerprints)

| Result | Probes |
|---|---|
| Flagged | A wide CSV copy of test labels |
| **Not flagged** | The same data as `.npy`, pickle, JSON or `.tsv.bz2` |
| **Not flagged** | ModelIDs without the `ACH-` prefix |
| **Not flagged** | 3,000 padding rows on top (only the first 3,000 rows are read) |
| **Not flagged** | Negated values with bare-Entrez columns |
| **Works as an oracle** | A long table of 30 known training cells plus one test cell, with 60 candidate values in separate columns. The per-column "n/m cells match" counts revealed the toy test cell's value exactly (−0.36). |

**Git hooks**
- `git commit --no-verify` and `core.hooksPath` bypass both hooks.
- Pre-push scans working-tree versions of tracked files, not the pushed commits. A label copy committed and later deleted would still be pushed in history.

**PreToolUse hook** (PI-owned; decisions checked by feeding it events, nothing executed)
- Allowed:
  - `cp data/lock*/labels_te*`;
  - `data/dep[m]ap/CRISPRGene?ffect.csv`;
  - shell-variable concatenation;
  - a `python3 -c` that builds the file name;
  - `labels_train`, which is not in the pattern.
- Together with the system `python3`, which has no tripwire, this reads raw labels. The §11.1 "requested hardening" is still open.

## 5. Test suite

- Run with 4 cores and `NEXTGEM_NO_TEST_RECORD=1` (so `results/test_runs.tsv` was not modified): **143 passed, 0 failed**.
  - 82 are infrastructure tests (frozen 6, hook 7, labels 12, ledger 9, precommit 7, scanner 4, splits 22, tripwire 15), matching the README.
  - The metrics tests are now 61 (13 binary + 48 tAP). They were 56 when the README was written, because the metrics fixes were committed in 0c06695.
- An earlier run, before 0c06695, gave 138 passed, matching the README.
- `check_frozen.py`: OK (30 rows).

## 6. Hooks versus legitimate PI work

Hooks I saw pass in practice:
- the PI's real commits fe0fc46 and 0c06695, made with the hooks installed.

Simulated in a scratch clone with the hooks installed, all passing:
- the M1 code commit (src, scripts, tests, FROZEN.tsv, results), about 1.5 s;
- a DESIGN §17 amendment commit;
- a SEED_SOURCE commit (simulated seed);
- a ledger `non_look` append (`git commit --only`);
- pre-push (0.7 s).

In the clone, `git_checks` passes with `--code-commit` set to the amendment commit. It refuses with `--code-commit` set to the code commit when the amendment came later, which is correct: the README says the amendment must be in or before C.

**Over-blocking found (S1).** One failed test in the latest full run blocks every commit. In the clone I recorded a fake full run with one failing test. After that:
- a notes-only commit was rejected;
- `ledger.append('started')` was rolled back, so the look aborted.

**Future friction.** The shape rule (b) and the test-keyed rule will flag legitimate per-line binary prediction tables for validation and test lines (P57 calls). These tables should go under `sealed/` or `results/**/tmp/`.

## Must-fix

1. **M1, production overrides in make_splits.py.** Without `--sandbox`, refuse:
   - `--reference-commit` ≠ the freeze commit;
   - `--repo` ≠ PROJECT_ROOT;
   - `paths.is_redirected()`.
2. **M2, sandbox access to the production seed.** `--sandbox --seed 38` is accepted today: SEED_SOURCE.md does not exist yet, and the seed is extracted with a loose regex.
   - Require sandbox seeds ≥ 500000, as the README states.
   - Or refuse any sandbox seed whose range [s, s+99] meets [38, 137].
   - Exclude 38–137 from tests that assert per-seed draws. On failure, `test_acceptance_stats_and_frozen_function_equality` would print both tiers for that seed.
3. **M3, pin the remaining draw inputs.**
   - Require exact equality of `expr_by_lin`/`prot_by_lin` with `coverage_by_lineage` in the 61b0beb JSON. They are equal today.
   - Compute the accepted attempt independently with the frozen `fdraw`/`fchecks`, and assert that it equals `run_attempts`.
   - Check the SHA-256 of Model.csv, CRISPRScreenMap.csv, the expression file and the proteomics inputs against `data/checksums.tsv` at 61b0beb. Make the Gygi file mandatory, or record its hash now.
4. **M4, tripwire entry-point trust.** Grant raw, validation, test, sealed and fingerprint access only if both hold:
   - the labels.py frame is the real project file, or a pinned worktree whose labels.py SHA equals HEAD:FROZEN.tsv;
   - the entry script's SHA is in the committed FROZEN.tsv.

   `evaluate_val.py` and `evaluate_test.py` must be frozen before first use anyway (§11.3).
5. **M5, amendment-id parser.** `frozen.amendment_ids` accepts the table's header words as amendment IDs: Date, Section, Change, Reason, Before, after, the, affected, analysis. So `check_frozen.py --add X --amendment Change` changes a frozen row with an empty log. Accept only IDs from data rows: skip the header and separator rows, and use an explicit ID token such as `A1`.
6. **M6, scanner oracle** (before `build_labels.py`).
   - Fingerprint the validation and test cells separately from training, so that training padding cannot reach the thresholds. Or fingerprint held-out cells only.
   - Never print match counts.
   - Treat fingerprints.npz as sealed-equivalent: with the salt, a 32-bit hash of 2-dp values is brute-forceable into a near-copy of all tiers.

## Should-fix

- **S1, test-run check.**
  - Scope the downgrade check to commits that touch src/, scripts/, tests/ or FROZEN.tsv.
  - Compare against a frozen test-ID baseline, not just the previous full run. Today a downgrade is forgotten after one more full run, and test_runs.tsv is an untracked, editable file.
- **S2, split marker.** Create `SPLIT_SHA256` only after `split_started` commits, or remove it if that commit fails. At present a failed ledger commit permanently blocks the production draw.
- **S3, tripwire hardening.**
  - Load the tripwire from a `.pth` file as well as sitecustomize.
  - Wrap `pyarrow.fs.LocalFileSystem.open_input_file`/`open_input_stream`, `pyarrow.OSFile` and `pyarrow.ipc`.
  - Extend inode protection to the tier files and to data/locked.
  - Make the process-argument check glob-aware.
- **S4, scanner and pre-push.**
  - Staged or pushed `.npy`/`.npz`/`.pkl`/`.h5`, JSON and `.bz2`/`.xz`/`.zip` files under results/ count as hits unless allow-listed.
  - Read every row of staged tables up to 50 MB.
  - Match ModelIDs against the known ID set, with or without the prefix.
  - Pre-push scans the blobs of the pushed commit range.
- **S5, build_labels.py overrides.** `--data-root`/`--state-root` are accepted while the raw files are the real ones (symlinks or inodes).
  - With them, the tier files and a fingerprints.npz that no longer carries a protected name can be written outside the real data root.
  - With `--state-root` alone, the real tier files are rewritten with a new salt, because the FROZEN lookup goes to the other root.
  - Refuse unless the raw files are not the real ones, or the real FROZEN.tsv has no tier hashes.
- **S6, DESIGN additions.** Restrict post-freeze DESIGN.md additions to the §17 table (see §3 b).

---

## Fix verification (round 2, 2026-10-26)

**Reviewed state**
- HEAD 57670ff, working tree uncommitted.
- File hashes: FROZEN.tsv `b5ff0c4a…`; make_splits.py `6d2a2c77…`; splits.py `9fc97d6a…` (unchanged); tripwire.py `742734d0…`; scan.py `b98a1e5e…`; amendment_A1.md `ee52ce9b…`.

**Method**
- The same rules applied: no seed in 38–137 was evaluated in any mode, nothing was committed, and no reviewed file was edited.
- The probes reran against a fresh toy store.
- Commit flow was simulated in a scratch clone with the hooks installed.

**Verdict: GO for the real draw**, once the A1 order of operations is followed. No blocker remains. One cheap pre-draw fix is recommended (R1 below).

### Must-fixes

| # | Re-test | Result |
|---|---|---|
| M1 | `make_splits.py` without `--sandbox`, seed 900000: tried `--repo /tmp`, `--reference-commit 0c06695`, and `NEXTGEM_STATE_ROOT=/tmp` | All three refused before any work. A production run on the dirty tree is refused at the git checks, before any ledger row or marker is written. |
| M2 | `--sandbox` with seeds 499999 and 1000 | Both refused (≥ 500000 is required) before any draw. The tests exclude 38–137 from per-seed loops. The 1–1000 acceptance statistics run only in production. |
| M3 | Exact `coverage_by_lineage` check; `split.*` pins in FROZEN.tsv; inputs against the freeze `data/checksums.tsv`; accepted attempt compared with a loop on the frozen functions | **Fixed.** In memory against the real data: coverage check passes; pins equal the committed rows (clone after C). My own line → lineage map, from CRISPRScreenMap and Model.csv, hashes to the same value as `split.line_lineage_sha256`. Gygi and CTRP, which the freeze checksums lack, are covered by `split.inputs_sha256`. |
| M4 | Probes P10 (spoofed `__main__` naming the unfrozen `evaluate_test.py`) and P11 (copied `src/` plus own `scripts/evaluate_test.py`) | **Both now denied** ("entry point differs from its FROZEN.tsv SHA-256"). Spoofing `__main__` to a *frozen* script (`build_labels.py`) still works; this is a deliberate act, listed as residual. |
| M5 | `amendment_ids` on the current DESIGN.md, then with the A1 row pasted | Current file: ∅. With A1: {A1}. Header words no longer count. **Fixed.** |
| M6 | Oracle probe S9: 30 training cells plus 1 test cell with 60 candidates | The output is now one file-level line, `fingerprint: matches label values`, with no counts or columns. The hit comes from the training cells alone. Held-out hits need ≥ 5 matches and ≥ 50 %, so the oracle no longer reveals test values. **Fixed.** |

**Other re-tests**
- P9, a `subprocess` shell glob, is now denied: the process-argument check is glob-aware.
- A shell *script file* run with `sh file` still copies and hard-links (P14, P16). This is a residual: links outside the real data root are listed.
- The S2–S8 copies are now flagged: pickle, JSON, `.bz2`, unprefixed IDs, negated values, padding.
- Exception: a 5 × 70 unkeyed `.npy` is still not flagged. It falls under the ≥ 10-row minimum and is toy-sized; a real copy would be flagged.
- Unchanged and accepted as residual (README, "Accepted residual risk"): PYTHONPATH sitecustomize, `pyarrow.fs` and `OSFile`, ctypes, `-S`, tripwire-state tampering.

### Amendment A1 text

**Sound.** The row is one line with 5 cells, and its only amendment id is A1. It states:
- the `--code-commit` rule;
- the additions-only-to-§17 rule;
- the frozen-function and accepted-attempt cross-check;
- the input pins;
- the production and sandbox refusals;
- the pre-freeze seed deviation;
- the attestation;
- the test retirements.

In the clone, with A1 pasted and committed as C (DESIGN.md, FROZEN.tsv, src, scripts, tests, results), plus a SEED_SOURCE with a *simulated* seed:
- `git_checks(…, code=C, 900000)` returns `[]`;
- `design_additions_failures` returns `[]`;
- `check_frozen` reports OK (36 rows). The two retired test IDs are accepted through `results/test_retirements.tsv`, which carries A1.

**Wording.**
- "S was kept out of all files until 5a8c218" omits the amended-away slip that RESEARCH_TREE D4 records. Add "(one unpushed slip amended away; LESSONS.md)".
- The attestation is written on behalf of other authors. The PI should confirm it with each of them, including the D5 list author.

**Expected state before C.** In the project working tree, `check_frozen.py` currently reports the two retired tests as missing. That is correct: A1 is not in DESIGN.md yet. The check passes once A1 is in C.

### Pre-commit, pre-push and test suite

| Check | Result |
|---|---|
| Commit C (whole M1 code + A1) | passes in 1.7 s |
| Notes-only record commit | passes in 1.1 s |
| Ledger `non_look` append (`git commit --only`) | passes in 1.1 s |
| SEED_SOURCE commit | passes |
| Test-run check on record commits | skipped (S1 fixed) |
| `SPLIT_SHA256` marker | now created after `split_started` is committed (S2 fixed) |
| Pre-push | scans the pushed commit range and runs the full suite (~85 s). That is acceptable, because pushes are rare. |
| Full suite, 4 cores, not recorded | **155 passed, 0 failed** (84 s) |

### Remaining recommendation (not a blocker)

**R1. Seed check in `make_splits.py`.**
- The SEED_SOURCE seed check only looks for the seed as a standalone number anywhere in the file.
- The template contains the full freeze SHA, so a mistyped `--seed` such as 6, 8, 55, 61, 62, 96, 780 or 2881 would pass it. The one-shot draw would then run with the wrong seed.
- Fix: parse exactly `Seed S = (\d+)` and require `--seed` to equal it.
- Without that fix, the PI must use the literal command in the A1 instructions (`--seed 38`).
