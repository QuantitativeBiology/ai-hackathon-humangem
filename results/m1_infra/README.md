# M1 infrastructure: labels, splits, guards

Built on 2026-10-26 by the M1 infrastructure agent, following DESIGN.md v3 (frozen at 61b0beb): §3, §4.1, §5, §11 and §12. The review fixes from `notes/reviews/M1_infra_review.md` (ACCEPT WITH FIXES) have been applied. The draft of amendment A1 is in `amendment_A1.md`.

Nothing here has been run on real labels or with the production seed. Simulated seeds 38–137 are never evaluated outside the production draw.

## What was built

| Component | Files | DESIGN |
|---|---|---|
| Paths | `src/nextgem/paths.py`, which honours `NEXTGEM_DATA_ROOT` and adds `NEXTGEM_STATE_ROOT` for toy/sandbox isolation; `src/nextgem/config.py` | §7, §11.1 |
| Label door | `src/nextgem/labels.py`, `scripts/build_labels.py`, `src/nextgem/cv.py` (LOLO, nested) | §5, §11.1 |
| Split | `src/nextgem/splits.py` (the draw, imported by design_numbers.py), `src/nextgem/cohort.py` (label-free covariates; every input is mandatory), `scripts/make_splits.py` | §4.1 + A1 |
| Tripwire | `src/nextgem/tripwire.py`, loaded by `.venv/.../sitecustomize.py` (**installed**; `scripts/install_tripwire.py`) | §11.1 |
| Scanner | `src/nextgem/scan.py`, `scripts/scan_label_copies.py` | §11.1 |
| Ledgers | `src/nextgem/ledger.py`, `scripts/verify_ledgers.py`, `scripts/git_locked.py` | §11.3 |
| Frozen manifest | `FROZEN.tsv` (36 rows: M1 code and tests, `notes/apriori_lists_approved.md`, five `split.*` input pins), `src/nextgem/frozen.py`, `scripts/check_frozen.py` | §11.3 |
| Git hooks (**installed**) | `.git/hooks/pre-commit` runs `scripts/pre_commit_guard.py`: fast, scanner + check_frozen only, no tests. `.git/hooks/pre-push` runs `scripts/pre_push_guard.py`: scanner on the blobs of the pushed commits, check_frozen, and the full test suite. Both are installed by `scripts/install_git_hooks.py`. | §11.1, §11.3 |

## Review fixes

**M1 — production overrides.** Without `--sandbox`, `make_splits.py` refuses `--repo`, any `--reference-commit` other than the freeze commit, and `NEXTGEM_DATA_ROOT`/`NEXTGEM_STATE_ROOT`.

**M2 — sandbox seeds.** Sandbox seeds must be ≥ 500000. The 1–1000 checks skip seeds 38–137 outside production.

**M3 — pinned draw inputs.**
- Input files are checked against `data/checksums.tsv` at 61b0beb. The Gygi and CTRP files, which that file lacks, are covered by the `split.inputs_sha256` pin.
- Per-lineage expression and proteomics counts must equal `coverage_by_lineage` in the frozen JSON.
- FROZEN.tsv holds `split.*` pins for the input manifest, the line→lineage map, the per-lineage counts and the draw parameters. `make_splits.py --compute-pins` writes them; it draws nothing.
- The accepted attempt must equal a loop run on the frozen 61b0beb functions.

**M4 — tripwire trust.**
- Access to raw, validation, test, sealed and fingerprint files, and hashing, needs `labels.py` and the entry script to carry their SHA-256 from FROZEN.tsv.
- Real data uses the *committed* `HEAD:FROZEN.tsv` of the project. Toy data may also use the working tree.
- A copied `src/` with a fake `scripts/evaluate_test.py` is denied. So is a spoofed entry name of a script that is not frozen.

**M5 — amendment ids.** Only explicit `A<n>` ids in data rows of the §17 table count. Header words and the placeholder do not.

**M6 — scanner oracle.**
- The output is `HIT <path> <reason class>` only: no counts, values or columns.
- Training and held-out cells are fingerprinted separately. A hit needs ≥ 5 held-out matches and ≥ 50% of held-out-keyed cells, or ≥ 20 training matches and ≥ 50%, so training padding cannot trigger a hit on held-out cells.
- `fingerprints.npz` is sealed-equivalent: only a frozen `scan_label_copies.py` may read it.

**Item 7 — more formats.**
- Read: npy/npz (`allow_pickle=False`), pickles (restricted unpickler: pandas/numpy containers only), json, and `.gz/.bz2/.xz/.zip`.
- ModelIDs are recognised without the `ACH-` prefix: zero-padded strings, or unique integer columns matching the known ID set.
- Values are matched as given and negated.
- Files up to 200 MB are read completely, and all held-out rows are checked.
- Unkeyed Chronos-like arrays are flagged.
- In staged/pushed mode, binary containers under `results/` are hits unless listed in `results/scanner_allowlist.tsv`.

**Item 8 — fast pre-commit; suite in pre-push.**
- The pre-commit test-run check applies only to commits that touch `src/`, `scripts/`, `tests/` or `FROZEN.tsv`. Notes, results and ledger rows are never blocked.
- The downgrade baseline is the union of all earlier full runs (S1). Retirement needs `results/test_retirements.tsv` plus an amendment id.

**Item 9 — amendment A1 as tightened.** After 61b0beb, DESIGN.md may only gain rows of the §17 table; the placeholder row may be replaced. This was checked against the real DESIGN.md with the drafted A1 row.

**Also fixed.**
- S2: the `SPLIT_SHA256` marker is created only after `split_started` is committed.
- S5: `build_labels` refuses redirected roots that point to real data.
- S3 (part): the process-argument check is glob-aware, and inode protection covers the real `data/locked` and `data/labels`.

## Accepted residual risk (review §4; not closed by design)

These bypasses remain open. Mitigation is the PreToolUse hook, the content scanner (milestones, pre-commit, pre-push on pushed history), the label-exposure register and non-author review.

**Tripwire.** It covers the project venv's normal launches only. Ways around it:
- `python -S` or `-I` with `site-packages` appended to `sys.path`;
- system `python3`, `Rscript` and `sqlite3`;
- a `sitecustomize.py` earlier on `PYTHONPATH` (note: `PYTHONPATH` includes the current directory);
- native readers that open files themselves: `pyarrow.fs.LocalFileSystem.open_input_file`, `pyarrow.OSFile`, ctypes `open`;
- hard links to files created after the process started, and hard links outside the real data root;
- shell `cp`, `cat` or `ln` run outside Python;
- deliberate in-process tampering: `nextgem_tripwire._state`, monkeypatching `check`, faking `__main__.__file__` to a *frozen* entry script, or faking a code filename.

**Git.**
- `git commit --no-verify` and `core.hooksPath` skip the hooks. Pre-push still scans the pushed history, but `git push --no-verify` skips that too.
- Untracked files are only found by milestone scans.

**Scanner.** Copies it cannot detect:
- encrypted, transformed (scaled, ranked, permuted) or base64-encoded copies;
- copies without ModelID/Entrez keys whose shape is not Chronos-like;
- pickles that need classes outside pandas/numpy;
- HDF5/joblib files outside `results/`.

**PreToolUse hook** (owned by the PI). These are still allowed:
- globs (`data/lock*/…`, `dep[m]ap`);
- shell variables;
- built file names;
- the `labels_train` name.

The §11.1 hardening request is still open.

**Commit-reveal.** S was revealed in 5a8c218, before the code existed. Tier *names* for S = 38 are computable by anyone holding 61b0beb. A1 pins the draw to pre-reveal artefacts. Each post-reveal author attests that they did not evaluate `draw()` on seeds 38–137.

## How to run

```
.venv/bin/python -m pytest tests -q --basetemp=<scratch dir>     # full suite; appends results/test_runs.tsv
.venv/bin/python scripts/scan_label_copies.py                     # milestone scan (project + /tmp/claude-*)
.venv/bin/python scripts/check_frozen.py [--staged] [--no-test-runs]
.venv/bin/python scripts/check_frozen.py --add <path...>          # refresh rows not yet committed
.venv/bin/python scripts/verify_ledgers.py
.venv/bin/python scripts/make_splits.py --seed <>=500000> --freeze-commit 61b0beb --sandbox <tmp dir>   # dry run
```

The production order (commit C with A1, then SEED_SOURCE.md, then the draw, then the labels) is in `amendment_A1.md` §3. The draw command is:

```
.venv/bin/python scripts/make_splits.py --seed 38 --freeze-commit 61b0beb --code-commit <C>
```

## Test results and runtime (after the review fixes, 2026-10-26)

Full suite: **155 passed, 0 failed**, in 85 s wall-clock and 7.7 CPU-min. The run is recorded in `results/test_runs.tsv`.

The 155 tests are 94 M1-infrastructure tests plus 61 metrics tests. The infrastructure tests by file:

| File | Tests |
|---|---|
| splits | 27 |
| tripwire | 18 |
| labels | 12 |
| ledger | 9 |
| precommit | 9 |
| frozen | 7 |
| hook | 7 |
| scanner | 5 |

Two pre-review tests were replaced. They are listed in `results/test_retirements.tsv` under A1, and their retirement takes effect once the A1 row is pasted.

New coverage:
- the production-override and sandbox-seed refusals;
- the pins and coverage equality, including tampered-pin and tampered-input negatives;
- the attempt loop against the frozen loop (200 simulated seeds ≥ 600000, including redraws);
- §17-only DESIGN additions, checked against the real DESIGN.md with the A1 row: it passes, and an added sentence fails;
- copied-checkout and spoofed entry points (denied);
- glob arguments (denied);
- amendment ids taken only from data rows;
- the union test-run baseline and retirement;
- record-only commits not blocked by a failing test;
- pre-push catching a label copy that was committed and then deleted;
- the scanner, flagging copies in pkl, json, tsv.bz2, npy, unprefixed+negated and 3,000-row-padded forms, with coarse output: the oracle probe gives identical output whichever column holds the true value.

Other checks:
- **Milestone scan** (project and `/tmp/claude-*`): 0 hits, 20 s. Fingerprints are not built yet.
- **`check_frozen.py`**: OK (36 rows).
- **`make_splits.py --compute-pins`**: 7 s; draws nothing.
- **Sandbox draw**: ≈ 7 s.
- **Live hook probe**: `Read data/locked/HOOK_PROBE` was denied.
