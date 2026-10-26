# Design review, round 2: AI-agent failure modes

Reviewer: the same independent agent as round 1 (not an author of DESIGN.md). Date: 2026-10-26.

**Reviewed:**
- DESIGN.md v1 (d84a593);
- notes/design_review_r1/PI_decisions.md;
- `.claude/settings.json` and `.claude/hooks/guard_labels.py` (8fa72d5);
- `scripts/design_numbers.py`, the `draw()` function and its docstring only;
- the file list of the shared session scratchpad.

**Label safety.** I read no label values. Hook tests used synthetic JSON. The audit-hook probe used a 5-row toy file in my scratchpad, which I deleted afterwards. I did **not** evaluate the split draw for any seed.

## Verdict: REVISE. Four new must-fix items, all cheap

Three of the four are about things that already happened or are already on disk, not hypothetical agent behaviour:
- the test tier is effectively known;
- a full all-line label copy sits under an unprotected name;
- the hook is not active in the running session.

Everything else from round 1 is adequately addressed on paper.

---

## 1. Round-1 must-fixes: status in v1 + hook

| R1 | Status | Remaining gap |
|---|---|---|
| **MF1** seal raw labels | **Partly addressed.** §11.1 has a single door, an audit-hook tripwire, three label files and a raw-level canary. Encryption was rejected (PD15); the hook was added instead. | The hook is **not live** in the running session and fails open (N3). It does not protect `labels_val`, `data/locked/` or cached copies, and plain `grep -r` over `data/` gets through (§3). The tripwire misses non-venv interpreters and direct pyarrow reads, and its frame-based allow list lets any caller of `labels.py` internals through (§4). An all-line label copy already exists outside the door (N2). |
| **MF2** label-derived denylist | **Addressed** (§11.2, W5, W7 whitelist, Human2 caveat). | none |
| **MF3** fold-scoped access, nested CV | **Addressed** (`labels.fold`, `cv.lolo`, nested choices; F5 nested; CE_fit per fold). | none |
| **MF4** ledgers, caps, idempotent looks, rehearsal, split run-once | **Addressed on paper** (§11.3, §10.2, §4.1). | "Refuses to overwrite" does not stop **pre-computing** the draw, which has already happened (N1). There is a commit-rights conflict (S3). |
| **MF5** frozen manifest, test ledger | **Addressed** (§11.3). | none |
| **MF6** provenance, clean tree, pinned worktree, atomic writes | **Addressed** (§11.4). | Scope the dirty-tree check to code (S4). |
| **MF7** solver failures abort | **Addressed** (§11.4: slot pool, SolverInfraError, solver assertion, MILP gap rule, no-information score). | none (see R2-2 below) |

---

## 2. NEW must-fix

### N1. The test tier is already known: re-seed the split with an input that no agent can preview (§4.1, §7, §11.5)
This concurs with the round-2 statistics and external-PI reviews. I confirmed the mechanism from the code alone.
- **The draws are nested.** `design_numbers.py::draw()` calls `rng.permutation` once per stratum, in the order S1, S2, S3, and takes test then validation as prefixes. The permutations for seed 20261026 + attempt 0 are therefore identical under v0 and v1. v1's test tier is v0's test tier plus perm_S1[1] and perm_S3[2], both of which were in v0's validation tier.
- **Those tiers are in a committed file.** `notes/design_review_r1/statistics.md`, committed in fe15541 and read by the author and the PI, names both v0 tiers. Its scratch script `testtier.py` computed label aggregates on them: positives per line, eligible-gene counts, and a lineage prevalence/SD remark.
- **The design text is false on this point.**
  - §4.1 says "The real draw is unseen … no one knows the test lineages before M1".
  - §11.5 says these exposures are "not line-level test information under the 7-lineage split, which has not been drawn".
- **The a priori lists in v1 were written with this knowledge in reach.** That includes the background groups, gene panels, drug classes and the ProCan-primary choice.

**Fix.**
1. Retire 20261026. Add a §12 check that `make_splits.py` refuses it.
2. Take the new base seed from an input that **cannot be ground by the agent who freezes the design**. Options, in order of preference:
   - (a) a number the user sends at M1;
   - (b) the drand beacon round fixed as "the first round after the M1 freeze commit's timestamp", recorded with its signature.

   A seed derived from the freeze-commit SHA alone, as the external PI proposed, is reproducible but **grindable**. The committing agent can amend the message or timestamp and preview draws locally. If it is used anyway, `make_splits.py` must create the freeze commit and draw in one process, and log every local commit (reflog) as part of the attempt log.
3. Correct §4.1 and §11.5. Log the predicted-tier exposure in the register (S5) as "v0 tiers named and label aggregates computed on them; superseded by the re-seed".

The exposures themselves were aggregate, so a new seed fully neutralises them. No line-level re-masking is needed.

### N2. Purge the all-line label caches that live outside the label door, and detect new ones by content, not by name (§11.1, the pre-commit rule)
- **What is on disk.** The shared session scratchpad, which every sub-agent is told to use, holds `ge_U.parquet` (38 MB). Its producer `prev.py` shows that it is `CRISPRGeneEffect.csv` restricted to U, for **all 1,361 lines**. The same directory has `testtier.py` and a `splits.pkl`.
  - The scratchpad path is `/tmp/claude-10105/-data-tmpA-nextgem/c5c0f45b-d54c-4e57-bdef-3d0657510b12/scratchpad/`.
  - Neither the hook (which matches filenames) nor the tripwire (which matches paths) protects this file.
  - It is the easiest label source an agent could find: a 38 MB parquet with a helpful name, in the directory the agent is told to use.
- **The R1 SF6 pre-commit rule would miss it.** That rule matches test ModelIDs together with *label-like column names*. A wide label matrix has gene columns (`ACACA (31)`), not a column called `label`.

**Fix.**
1. Before M1, the PI deletes every label-derived scratch artefact: `ge_U.parquet`, and any pkl or csv built from CRISPR files by the round-1 reviewers or the design author. The deletion is logged in the register.
2. Add a content scanner, `scripts/scan_label_copies.py`, run by `status.py`, by the pre-commit hook and at each milestone. It looks under `results/`, `notes/`, the scratchpad roots and `/tmp/claude-*`. It flags any table that meets all of these:
   - its index or a column is mostly `ACH-` IDs;
   - it has ≥ 50 columns matching `^\S+ \(\d+\)$`;
   - its numeric values include negatives below −0.3.

   The negative values distinguish Chronos from expression and CN. The same rule flags 0/1 matrices whose rows are val/test ModelIDs and whose columns are gene IDs. Hits block the commit, and in the scratchpad they must be deleted.
3. Add one line to CLAUDE.md: "Never write label values outside `data/locked/` or `labels_train.parquet`; cache derived *statistics*, not label matrices."

### N3. The PreToolUse hook is not active, and it fails open (`.claude/settings.json`, `guard_labels.py`; §11.1)
- **Not active.** A live probe from this sub-agent was allowed through. I asked the Read tool for the non-existent `notes/design_review_r2/labels_test_hook_probe_nonexistent.txt`, and it returned "file does not exist" rather than the guard's denial. Claude Code snapshots hooks at session start, so the PI session and every sub-agent it spawns run **without** the guard until the session is restarted and the hook is approved in `/hooks`.
- **Fails open.** On malformed or unexpected input, the script raises and exits with code 1 (tested). Claude Code treats exit codes other than 0 and 2 as non-blocking, so the tool call proceeds. As far as I know, `"onFailure": "block"` is not a documented hook key; verify it, and do not rely on it.

**Fix.**
1. Restart the PI session, or approve the hook in `/hooks`, before any further sub-agent is launched.
2. Add a **liveness probe** to the session-start routine and to each milestone checklist: a Read of the sentinel `data/locked/HOOK_PROBE` must be *denied*. Record the result in RESEARCH_TREE.md.
3. Wrap `main()` in `try/except` and `deny` on any exception (fail-closed).
4. Add `data/locked`, `labels_val` and `labels_train` to `PROTECTED`. Direct reads of training labels bypass fold scoping, and validation labels must be reached only through `evaluate_val.py`.

### N4. T1 "sealed" outputs need a sealed location and a pre-commit window that lasts until T2 (§10.2, §11.1)
- **Why it matters.** §10.2 seals per-lineage, per-gene, attribution and failure-account outputs until P2 is frozen. This is the stats MF6 / PD21 protection of C4 and C5.
- **The gap.** v1 does not say **where** the sealed outputs are written. The pre-commit rule stops guarding test labels "before the T1 ledger row exists", so it switches off exactly when the sealed tables are created.
  - A sealed per-gene test table written under `results/` and committed in a routine PI commit is readable by every W7 agent.
  - W7 agents could then tailor P2 to the test errors.

**Fix.**
1. `evaluate_test.py` writes sealed outputs only to `data/locked/sealed/<look>/`, which is git-ignored and protected by the hook and tripwire. Only the SHA-256 goes into the ledger.
2. Extend the pre-commit window to "before the **T2 unseal** ledger row".
3. Have `scan_label_copies.py` (N2) also flag tables keyed by test ModelIDs × genes during that window.

---

## 3. Hook bypasses (tested with synthetic tool-call JSON against `guard_labels.py`)

| Case | Result | Matters in practice? |
|---|---|---|
| `Read data/depmap/CRISPRGeneEffect.csv` | denied | — |
| `duckdb …'data/depmap/CRISPRGeneEffect.csv'`; `awk … data/depmap/Screen*.csv`; `for f in data/depmap/*` | denied | — |
| `grep -r ACH-000001 data/` (Bash), and the Grep tool with `path=/data/tmpA/nextgem/data` | **allowed**; prints the label rows of that line | **Yes.** Searching all data for a ModelID or symbol is a routine debugging move. |
| `cd data/depmap && head -2 CRISPR*` | **allowed** (the wildcard regex needs a `depmap/` prefix) | **Yes, plausible.** The design author's accidental print of L-prob values was this kind of action. |
| Grep tool `path="data/depmap/."` | allowed | low (needs intent) |
| Read or any command on `data/locked/labels_val.parquet` | **allowed** (only `labels_test` is protected) | **Yes.** It breaks "validation reveals scalars only". |
| Read or pandas on the cached copy `ge_U.parquet` | **allowed** | **Yes. It exists now (N2).** |
| `cp -r data/depmap /tmp/dd`, then a script reads the copy | allowed | medium (backups, "work on a local copy") |
| `find data -name '*Effect*' -exec head -2 {} \;` | allowed | low to medium |
| `git diff --no-index /dev/null 'data/depmap/CRISPRGeneEffect.csv'` (quoted path removed by the git exemption) | allowed | low (needs intent) |
| `.venv/bin/python scripts/x.py`, where x.py reads labels | allowed (by design) | covered only by the tripwire (§4) |
| Malformed or unexpected hook input | hook crashes with exit 1, so the call is **allowed** | **Yes** (N3) |

**Cheap hardening** (should-fix, beyond N3):
- Deny Bash and Grep calls that search recursively from `data`, `data/` or the repo root without an explicit `--exclude-dir=depmap --exclude-dir=locked`.
- Deny `cd` into `data/depmap` or `data/locked`.
- Deny `cp -r`, `rsync` and `tar` of `data/depmap` or `data/locked`.

I would not try to defeat deliberate obfuscation such as string concatenation or variables. The threat model is careless agents, not adversarial ones.

## 4. Tripwire (`labels.py` + audit hook) design bypasses

1. **The allow rule is by frame, not by entry point.** "Any open from a frame outside `labels.py`" lets any script that calls a private raw reader in `labels.py`, for example `labels._load_raw()` "for a quick look", read every line.
   - **Fix:** raw and val/test reads in `labels.py` also check `Path(__main__.__file__).name` ∈ {`build_labels.py`, `evaluate_val.py`, `evaluate_test.py`}.
   - The raw build is further allowed only while FROZEN.tsv has no label hashes.
   - **Matters: yes.** This is the most likely way an honest agent gets around the door.
2. **Interpreters outside the venv are not covered.**
   - System `python3` (miniconda base, pandas 1.5.3) does not load the venv `sitecustomize.py`. The hook script itself is run with `python3`, so agents will copy that habit.
   - `Rscript` and `sqlite3` are also on the PATH.
   - **Fix:** the hook denies Bash commands that run `python3`, `python`, `Rscript` or `sqlite3` without the `.venv/bin/python` path, when the command or heredoc mentions `data/` or `.parquet`. CLAUDE.md already mandates the venv.
   - **Matters: yes.**
3. **Some readers raise no audit event.**
   - Tested in the venv:
     - `pyarrow.parquet.read_table(path)` raised **0** audit `open` events;
     - `pd.read_parquet`, `pd.read_csv` (C and pyarrow engines) and `np.loadtxt` each raised 1.
   - **Fix:** `sitecustomize.py` also wraps `pyarrow.parquet.read_table`, `pyarrow.parquet.ParquetFile`, `pyarrow.dataset.dataset` and `pyarrow.csv.read_csv` with the same path check. It also handles the `subprocess.Popen` audit event, checking its args for protected paths.
   - **Matters: medium** (an agent "optimising" parquet I/O).
4. **Derived copies end the protection.** Once a label matrix is written elsewhere, nothing guards it. N2's content scanner is the only general answer.
5. **Stronger, cheap option: obfuscate at rest with an agent-readable key** (should-fix S1). This answers PD15's objection, because no user is involved.
   - After `build_labels.py`, store `labels_val`, `labels_test` and the raw CRISPR, ScreenGeneEffect and Dependency files as Fernet-encrypted blobs.
   - The key sits in `data/locked/.key`. It is opened only by `labels.py`, the opener is checked by entry point, and the key is hook-protected.
   - Delete the plaintext raw files; they can be re-downloaded and their checksums are in DATA.md.
   - Every accidental path in §3 and §4 (grep, system python, R, pyarrow, `cp -r`) then yields ciphertext.
   - This is not security against a determined agent, and it does not need to be.

---

## 5. Open questions in my area

- **R2-11: do reviewers' all-line label views compromise the test tier?**
  - **Aggregate all-line statistics on their own do not.** These are pooled prevalences, the common-essential share, ASNS prevalence by medium and per-library positive rates. They reveal nothing about any line, and they apply equally to every method, so they cannot favour a pipeline. They can inform design choices such as the SESOI or the NA rule, but those choices are method-symmetric.
  - **The real compromise is the combination** of (i) the real-seed draw being computed and published with (ii) label aggregates computed on that predicted test tier. The test lineages were then known while the a priori lists were written. N1 (re-seed) neutralises this completely, because the exposures were aggregate.
  - **The design author's ~10 L-prob values for ACH-002730** are negligible after the re-seed. Record the lineage in the register.
  - **Rule from now on:** no reviewer or agent computes any label statistic on non-training lines. After the split, label-informed review work uses `labels.train()` only.
- **R2-2: no-information score for failed lines.**
  - Keep it. Excluding the line rewards failure; the generic fallback hides failure; a per-gene e = 0 is a pseudo-prediction.
  - Add two things:
    - (a) every table shows the failed-line count per method;
    - (b) a sensitivity row that excludes lines failed by **any** compared method, identically for all methods.
  - The flag at > 5% failed lines stays.
- **R2-4: κ rule.** Acceptable. It is deterministic and pre-specified, it uses training lines only, it is symmetric across methods, and it is executed by `labels.py` at M1 before any model is scored.
  - The outcome is written to FROZEN.tsv in the same run.
  - No agent inspects κ and then decides; if one does, it is an amendment.
- **R2-9: `openpyxl` and guide maps.**
  - Fine from my angle, provided `requirements.txt` is regenerated and committed in the same commit as the install.
  - The guide maps are label-free.
- **R2-3: compute, from the cost-blow-up angle only.** Whichever option is chosen, enforce the per-job 2× kill and the `budget.tsv` check before Stage B launches. B7 on all lines is the single largest MILP block, so pilot its per-line cost on 10 lines first (§11.5 rule).

---

## 6. Should-fix (new or refined)

- **S1.** Fernet obfuscation at rest with an agent-readable, hook-protected key (§4.5).
- **S2.** Tripwire: an entry-point check, pyarrow wrappers, `subprocess.Popen` audit and a non-venv interpreter rule in the hook (§4.1–4.3). The hook hardening in §3 also belongs here.
- **S3.** The ledger commit conflicts with the single-committer rule (§11.3 vs §11.5).
  - `ledger.append()` is the **only** non-PI commit path. It uses `git commit --only <ledger>` under a lock file shared with the PI's commits.
  - The `started` row must be committed **before** any number is computed. If that commit fails, abort, so no unlogged reveal can happen.
- **S4.** Scope the provenance "clean tree" check to `src/`, `scripts/` and `tests/`. Otherwise sub-agent outputs under `results/` keep the tree dirty and the check gets whitelisted into meaninglessness.
  - In pinned-worktree runs, assert that `nextgem.__file__` lies inside the worktree. The editable install could otherwise shadow it.
- **S5.** `notes/label_exposure.md` does not exist yet, although §11.5 cites its entries. Create it now, with these entries:
  - the round-1 statistics predicted tiers and aggregates;
  - the data reviewer's aggregates;
  - `ge_U.parquet`;
  - ACH-002730.
- **S6.** The rehearsal look "on validation" (§10.2) counts as a validation look and gets a ledger row.

## Items checked, no change needed
- **Leakage controls:**
  - MF2 denylist and W7 whitelist;
  - fold-scoped API and nested selection;
  - CE_fit per fold;
  - B2 and B2m from the same label;
  - R3 as a label-free statistic;
  - the W5 reviewer working without label access.
- **Ledgers and freezing:**
  - hash-chained ledgers;
  - caps in code;
  - idempotent looks;
  - first T1 number as headline;
  - FROZEN.tsv and the test-run ledger.
- **Solver failures:** the SolverInfraError path, the GLPK declaration with an equality check, and the MILP gap rule.
- **Reporting and process:**
  - claims register and `verify_report`/`verify_citations`;
  - review front matter;
  - job registry and `status.py`;
  - 1.5× schedule buffer.
