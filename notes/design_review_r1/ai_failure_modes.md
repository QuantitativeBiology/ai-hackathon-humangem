# Design review, round 1: AI-agent failure modes

Reviewer: independent agent (did not write DESIGN.md). Date: 2026-10-26. Reviewed: DESIGN.md v0.1, PROMPT.md, RESEARCH_TREE.md, LESSONS.md, CLAUDE.md, DATA.md, notes/review_baseline_ko.md, `.gitignore`, `.claude/settings.json`, `src/nextgem/knockouts.py` (status handling only), and a listing of `data/depmap/`.
No label values were read, and no essentiality performance was computed.

## Verdict: REVISE (7 must-fix, 9 should-fix, 5 optional)

§11 names the right failure modes. Most of its guards, however, are **policies that a careless agent can bypass without noticing**, and it does not cover the main leak.

- **The locked test file is undermined by the raw label files.** `data/locked/labels_test.parquet` (mode 0400) is locked, but the raw label sources stay world-readable in `data/depmap/` and contain every line's labels:
  - `CRISPRGeneEffect.csv`;
  - `ScreenGeneEffect.csv`;
  - `CRISPRGeneDependency.csv`, once it arrives.

  The agents run as the same Unix user, so mode 0400 stops nobody. Gene priors, co-essentiality, W5b, the κ check, the "selective" stratum and every W7 feature will be written by agents who naturally reach for the raw CSV.
- **The other guards are greps and spot-checks.**
- **The ledgers are editable TSVs.**

The fixes below replace these policies with mechanisms. They are cheap to build now, before any label is read. Today nothing has been scored and the split has not been drawn, which is the last cheap moment to install them.

Ranking:
- **must-fix**: a plausible path to a wrong, leaked or overstated headline result, with no mechanical guard in the current design;
- **should-fix**: a real risk with a partial guard;
- **optional**: hardening.

---

## MUST-FIX

### MF1. Seal the raw label sources, not just the derived test file (§4.2 "Physical separation", §10 Splits row, §11 Leakage row)

**Failure.** An agent computes a gene prior, a co-essentiality matrix or W7 features with `pd.read_csv("data/depmap/CRISPRGeneEffect.csv")`. That read includes the test and validation lines.
- The §10 leakage canary does not catch this, because it perturbs only `labels_test.parquet`.
- The §11 grep for `data/locked/` and `CRISPRInferredCommonEssentials` does not catch it either. Paths are built from variables, and the leaking file is a different one.
- `ScreenGeneEffect.csv` (per-screen Chronos, 585 MB) is a second, unnamed copy of the labels.

**Change.**
1. **One label door.** `src/nextgem/labels.py` is the only module that may open label-bearing files. It exposes tier-scoped accessors:
   - `labels.train()` returns training lines only;
   - `labels.fold(k)`: see MF3;
   - validation is reachable only via `scripts/evaluate_val.py`;
   - test is reachable only via `scripts/evaluate_test.py`.

   Every data path comes from `src/nextgem/paths.py`, which honours an env var `NEXTGEM_DATA_ROOT`, so the canary can redirect it. A test fails if the string `data/depmap` or a label filename appears outside `paths.py`.
2. **Split the parquet files.** `build_labels.py` writes three files, not two:
   - `labels_train.parquet`;
   - `labels_val.parquet`, in `data/locked/`;
   - `labels_test.parquet`, in `data/locked/`.

   The validation cap (§9.2) is unenforceable while validation labels sit in a `trainval` file that every agent loads.
3. **Seal the raw files after `build_labels.py`.** Move them into an encrypted archive with `age` or `gpg -c`, using a passphrase that **the user holds** and that is never written to disk or to the conversation, and delete the plaintext:
   - `CRISPRGeneEffect.csv`;
   - `ScreenGeneEffect.csv`;
   - `CRISPRGeneDependency.csv`;
   - `labels_test.parquet`.

   Record the SHA-256 of each in DATA.md. Each file can be re-downloaded from DepMap, so nothing irreplaceable is lost. At T1 and T2 the user decrypts. That makes the user's action the authorisation for each test look, and it removes any reliance on an agent-written "PI authorisation line" (§9.3).
   - *Fallback, if the user declines to hold a key:* move the files to `data/locked/raw/` and apply the next item. That is weaker, but it still catches the realistic failure, which is careless rather than adversarial access.
4. **Access tripwire, in two layers.**
   - (a) A `sitecustomize.py` in `.venv/lib/python3.11/site-packages/` registers `sys.addaudithook` for `open` events. Any open of a path under `data/locked/` or of a label filename is appended to `logs/label_access.log` (timestamp, argv, cwd, git commit, PID). It raises `PermissionError` unless the calling script is `build_labels.py`, `evaluate_val.py` or `evaluate_test.py`. CLAUDE.md already mandates the venv interpreter, so this covers all project Python.
   - (b) A Claude Code `PreToolUse` hook in `.claude/settings.json` denies Read, Grep and Bash calls that name those paths. This needs the user's approval, because it is a configuration change.
5. **Canary at the raw level (§10).** Point `NEXTGEM_DATA_ROOT` at a copy in which the **test- and validation-line rows of the raw CRISPR files** are replaced by noise. Every training-tier output must then be bitwise identical. Run this at M4 and M6 on every shortlisted pipeline. It is cheap if it is restricted to the label-consuming steps; FBA sweeps are label-free.

### MF2. A denylist of label-derived external resources (§2.2, W5 evidence table, W7, §11)

**Failure.** Several datasets that the design plans to use, or that agents are likely to pick up, were **derived from CRISPR screens of the same lines, test lineages included**. Using them in fitting or curation leaks labels with no file in `data/locked/` ever being touched.

| Source | Where it enters | Why it leaks |
|---|---|---|
| **De Kegel 2021** paralog dependency predictions (download #11) | W5 GPR evidence (AND→OR) | trained on DepMap CRISPR essentiality across all lines |
| Ito 2021 paralog SL screens | W5 GPR evidence | CRISPR screens in cancer lines that overlap DepMap (check before use) |
| Project Score (Behan 2019), CoRe, Hart CEG/NEG, `AchillesCommonEssentialControls.csv`, older DepMap releases, DEMETER2 RNAi | W7 features, gene priors, "known essential" lists an agent adds "to help" | essentiality on the same or overlapping lines |
| `CRISPRInferredCommonEssentials.csv` | reporting (allowed); a tempting feature | inferred from all 26Q3 lines |
| PRISM, GDSC and CTRP responses of test lines | W7 features | phenotypes of the evaluated lines, closely tied to dependency; fine as W6 *targets*, not as predictor inputs |
| Human-GEM v2 GPR fixes (Luo 2026) | the generic parent H2, and the W5 positive control | LITERATURE §346 says Luo benchmarked on DepMap and Project Score; if GPRs were tuned to DepMap essentiality, generic H2 is already label-informed on every lineage |

**Change.**
1. Add a `label_derived` column to DATA.md with values `no`, `yes (lines: …)` or `unknown`. Add `LABEL_DERIVED_SOURCES` to `src/nextgem/config.py`. The `paths.py` accessor for such a source refuses unless the caller passes `purpose="reporting"` or `purpose="external_validation"`.
2. W5: move the De Kegel tables, and Ito unless it is shown to use no DepMap lines, out of W5a ("label-free") and into W5b. There they must be **recomputed within folds**, which is impossible for published tables, so in practice they are excluded and only cited.
3. W7: add a closed feature whitelist. The only allowed sources are expression, CN, mutation, proteomics, GEM outputs and the within-tier gene prior. Any idea in `ideas.tsv` that adds a feature source not on the list needs an amendment.
4. §2.1/§16: record as a stated limitation whether Human2's curation used DepMap labels (Luo SI, paper request already open). If it did, the H2-vs-H1 comparison and the W5 "recovery of v1→v2 fixes" control are partly label-informed. The design should say so up front, not discover it in REPORT.md.

### MF3. Fold-scoped label access for every LOLO fit (§4.2 Training row, §8 TG, W2 M3, W3 T5/Neg-2, W4 F5/R3, W5b, W6 PCA, B2)

**Failure.** The design asks for many statistics to be "recomputed inside each LOLO fold". The convenient shortcut is to compute them once on all training lines and then run LOLO. That leaks the held-out training lineage into its own prediction and inflates every training gate. Affected statistics:
- the gene prior;
- R3 local thresholds;
- CellFie percentiles;
- PCA loadings;
- W5b co-essentiality;
- W3 Negative-2's gene-prior ranking;
- F5's (w, τ).

In addition, choices are made by best training-LOLO mAP_line, and the TG is then computed on the same folds. That is selection on the evaluation folds, and the TG lower bound becomes optimistic. The affected choices are:
- F5's 48-point grid;
- M2a vs M2b;
- "best ladder variant";
- the W7 screen.

**Change.** Provide one harness, `nextgem.cv.lolo(method, fold_ids)`. Each method implements `fit(train_view)` and `score(line_ids)`. The harness builds `train_view` from `labels.fold(k)`, which physically excludes fold k, and a method never receives a raw label handle.
- Any hyperparameter choice inside a method runs an **inner** LOLO on `train_view` (nested).
- The TG for a candidate that was chosen from a set uses the nested outer-fold estimate.
- Add a toy test in which a method that memorises its training labels gets AP = π on the held-out fold.

### MF4. Make the ledgers and the look caps mechanical and tamper-evident (§9.2, §9.3, §4.2)

**Failure.** The ledgers are TSVs that any agent can edit or truncate.
- After a context reset, a PI can re-run T1 not knowing it already happened.
- A crash halfway through `evaluate_test.py` leaves an ambiguous state.
- The "PI authorisation line" is free text that an agent can write itself.
- Exceeding 60 validation looks needs only "a PI note".

**Change.**
1. **Hash-chained, append-only ledgers.** Every append goes through `nextgem.ledger.append()`, which takes a local lock, writes a row containing the SHA-256 of the previous row, and **commits the ledger immediately** with a fixed message. `scripts/verify_ledgers.py` checks:
   - the chain;
   - that `git log -p` of each ledger contains only additions;
   - that every validation and test output file has a ledger row.

   It runs at every milestone and inside every `evaluate_*` script.
2. **Validation cap in code.** `evaluate_val.py` refuses look 61 and later, and refuses W7 look 31 and later, unless `results/authorisations/val_cap_raise_<n>.md` exists and is committed. For W7 looks it prints only the scalar and its CI (W7 step 3); per-gene and per-line outputs are not written.
3. **Idempotent test looks.** Look IDs (T1, T2, T1-postfix) are unique. `evaluate_test.py` does the following:
   - writes a `started` row;
   - computes;
   - writes the output with its SHA;
   - writes a `completed` row.

   If any row for that look ID exists, it refuses and prints the path of the existing output. A crashed look can be resumed only with the identical config hash, and it still counts as the look.
4. **Dress rehearsal.** Before T1, run `evaluate_test.py --rehearsal` against a synthetic `labels_test` with the same shape but random labels, and then against validation. This catches bugs in the evaluation script before the real look, which is the main reason a "third look for a bug fix" (§9.3) would otherwise be needed.
5. **Freeze the split mechanically.** `make_splits.py` and its acceptance constraints are committed **before** the first run, and the output records the script's SHA. The script refuses to run if `results/splits/SPLIT_SHA256` exists. This guards against re-rolling once an agent sees which lineages landed in test.

### MF5. Frozen-artefact manifest, and a test ledger against "fixing tests to pass" (§11 Metric drift and Unverified fixes rows, §10, §16)

**Failure.** §11 freezes only `metrics.py`. Other artefacts determine results just as directly, and an agent can edit them while "fixing a bug":
- the label builder;
- the U definition and join;
- the evaluation scripts;
- known-answer test tolerances;
- `exchange_map.tsv`;
- `serum_model.tsv`;
- the W3 task lists;
- the W4 background groups;
- the W6 compound and subtype lists.

Agents also resolve failing tests by loosening tolerances, marking tests `xfail` or `skip`, or deleting them. The baseline's "36/37 checks pass, one known" shows how easily an expected failure becomes routine.

**Change.**
1. **Frozen manifest.** Keep `FROZEN.tsv` at the repo root with path, SHA-256, frozen-at milestone and amendment ID. It covers:
   - `metrics.py`, `labels.py`, `splits.py`, `paths.py`, `cv.py`;
   - `evaluate_val.py` and `evaluate_test.py`;
   - every known-answer test file and the tolerance constants module;
   - the split and label parquet SHAs;
   - each a priori candidate list, when it is fixed.
2. **Pre-commit check.** `scripts/check_frozen.py` runs in a git pre-commit hook and at the start of every `evaluate_*` script. A frozen file may change only if §16 has an amendment row naming the new SHA, and that row's commit must precede the first ledger row or result that uses it. The check reads both from git history.
3. **Test ledger.** Each pytest run of the verification suite appends to `results/test_runs.tsv`: commit, n_pass, n_fail, n_xfail, n_skip, and the hash of the sorted test-ID list. `check_frozen.py` fails if a known-answer test ID disappears, or changes to xfail or skip, without an amendment.
4. **One rule in CLAUDE.md.** "When a test fails, never edit the test or its tolerance in the same commit as the code fix." This makes test edits visible in review.

### MF6. Provenance wrapper, clean-tree execution and pinned code for long runs (§6.6, §11, §14)

**Failure.**
- **Mixed code versions.** `nextgem` is installed in editable mode, so several agents edit `src/nextgem/` while multi-hour sweeps run. W4 alone needs 34 h of wall time. Workers spawned later import the new code, and a single parquet then mixes code versions, with nothing to show it.
- **Untraceable numbers.** Agents also produce numbers with ad-hoc `python -c` snippets or notebooks, and those numbers flow into the tree and REPORT.md.
- **Truncated files.** The shared NFS volume was 98% full, so a disk-full mid-write can leave truncated parquets that later read as data.

**Change.**
1. **Wrapper for result-producing scripts.** Every script that produces a result runs through `nextgem.run.provenance()`, which extends the planned `nextgem.timing` manager. It refuses to start on a dirty tree, apart from a whitelist of untracked output dirs. It writes a `MANIFEST.tsv` row for each output:
   - script and argv;
   - git commit;
   - input SHA-256s;
   - seeds;
   - solver name and version;
   - number of workers;
   - wall and CPU time;
   - peak RSS;
   - output SHA-256;
   - author agent ID.
2. **Pinned worktrees.** Jobs longer than 1 CPU-h run from `git worktree add /data/tmpA/nextgem_runs/<commit>` with a non-editable `PYTHONPATH` pointing at that worktree's `src/`. Editing the main checkout then cannot change a running job.
3. **Atomic writes.** Write to `*.tmp`, fsync and rename. Before any job with an estimated output above 1 GB, check that at least 20 GB is free.
4. **Notebooks are exploratory.** A number cited in RESEARCH_TREE.md or REPORT.md must come from a `MANIFEST.tsv` output (see MF7).

### MF7. Solver and licence failures must fail loudly, never become predictions (§6.1 NaN and failed context models, §12)

**Failure.** The design maps two kinds of failure to predictions:
- NaN solves become e = 0 (non-essential);
- failed context builds fall back to the generic parent.

Several sub-agents launching 15 to 16 Gurobi workers each will exceed the WLS session limit. A brief network outage does the same, because the WLS licence needs internet access at solve time. Either produces licence or environment errors. Depending on where the error surfaces, it becomes one of three outcomes:
- (a) a crash;
- (b) a NaN that is silently scored 0, or a "build error" that falls back to generic;
- (c) cobra silently selecting GLPK as the default solver for that process.

Both (b) and (c) produce plausible-looking but wrong results. The failure rate is driven by load, not by the method, so it also differs between methods. In addition, ftINIT and other MILPs that stop at a time limit with a non-zero gap return a usable-looking but suboptimal model.

**Change.**
1. **Global Gurobi slot pool.** Use `nextgem.solver.slot()` with N = 16 lock files on **local** disk (`/dev/shm/nextgem_gurobi/slot_i`, `fcntl.flock`; not on the NFS volume, where flock is unreliable). Every worker acquires a slot before it creates a Gurobi env. A job that cannot get a slot waits; it never falls back.
2. **Separate failure classes.** Licence, env or network exceptions are a distinct class, `SolverInfraError`. The job aborts, writes no output and records the failure in the job registry. Only *mathematical* non-optimality (numeric or time-limit status after a successful env) may become NaN under §6.1. The NaN-to-0 rule then applies only to that class.
3. **Assert the solver.** Every KO parquet must have `solver == "gurobi"` on every row (the column already exists in `knockouts.py`), unless the run was explicitly declared a GLPK cross-check. Every loader checks this.
4. **MILP status.** Builds record MILP status, gap and time-limit flag. A time-limited build with gap > 1e-4 counts as a **failed build** in the §6.1 5% rule, and is not silently accepted.

---

## SHOULD-FIX

### SF1. Label-exposure register (§4, §16)
Create `notes/label_exposure.md` and append to it every time label *values* (not headers) are seen, by whom and on which lines. As of today the entries are:
- the header only (baseline review);
- the all-line common-essential count 419/2,798 (literature agent);
- nothing else.

Whoever writes the a priori candidate lists (W3 tasks, W4 background groups, W6 compound list, serum list) states in the commit message that they have no entries in the register. This turns "fixed a priori" into something that can be checked.

### SF2. Two label-informed decision points invite "fix until it passes" (§13 M2, §9.3)
- **M2 STOP rule.** "H2 lift ≤ 1 → label/join audit" lets an agent audit until lift > 1.
  - Specify the audit as **label-free** checks only: Entrez join on known genes, NaN counts, orientation of the sign of Chronos.
  - The audit is run by an agent that has not seen the metric.
  - Report the pre-audit and post-audit numbers together.
- **Bug-fix third look.** The post hoc third look after a T1 "bug fix" is an outcome-dependent option.
  - The **first** T1 number stays the headline in REPORT.md's summary.
  - Any fix must be applied identically to P1 and to every baseline.
  - The bug must be shown by a label-free failing regression test.
  - Bugs found by inspecting T1 errors are reported, not re-evaluated, unless the PI writes down why the bug would have been found anyway.

### SF3. Hand-transcribed inputs are the most likely hallucination path (W2 inputs, W5 evidence table, §11 Hallucination row)
The planned 20% spot-check is too weak for the inputs that drive the W2 headline. For these, require **double independent extraction**:
- media compositions transcribed from vendor PDFs;
- FBS reference concentrations;
- `exchange_map.tsv`;
- UniProt "function unrelated" judgements.

Two agents, without access to each other's output, each produce the table from the source file in `data/media/` or `papers/`. A script diffs the two tables, and every discrepancy is resolved against the source with a page or row reference. Every row carries `source_file` and `page/table`.
- W5: generate GPR edits **by script** from the local CORUM, Complex Portal, Ensembl and Nusinow files. A validator checks that every evidence ID exists locally and lists the genes involved.
- UniProt judgements must cite an accession and an entry version from a downloaded UniProt file. Add that file to DATA.md, since it is not there now.
- Freeze `serum_model.tsv` and `exchange_map.tsv` in `FROZEN.tsv` (MF5) **before** the first W2 scoring. An edit after scoring is an amendment.

### SF4. Context loss and concurrent writers (§13, §14; RESEARCH_TREE "the tree is your memory")
Concurrent writers to the record files need rules:
- Only the PI writes RESEARCH_TREE.md, LESSONS.md, DESIGN.md, DATA.md and `data/checksums.tsv`.
- Sub-agents write only under their own `results/<ws>/` and return a summary. Today RESEARCH_TREE.md already has uncommitted edits, which shows the hazard.
- Only the PI commits, and never with `git add -A`. A sweeping add pulls another agent's half-written files into a commit.
- Shared TSV appends go through the locked `ledger.append` (MF4).

To support recovery after a context reset:
- **Job registry.** Add `results/jobs.tsv` with job ID, owner agent, PID, host, command, worktree commit, expected outputs, CPU budget and status.
- **Status script.** Add `scripts/status.py`. It reads the registry, the ledgers and the frozen manifest, and checks whether PIDs are alive, outputs exist and SHAs match.
- After a reset, the PI runs `status.py` before launching anything. Without this, a 34 h W4 build can be launched twice, or orphaned.
- Sub-agent summaries are saved as `results/<ws>/agent_summaries/<date>_<agent>.md`, so they survive the PI's context.

### SF5. Cost and time guards (§12)
- **CPU budget.** Each registry job declares a CPU-h budget, and the runner kills it at 2× the declared budget.
- **Pilot first.** Any step estimated above 10 CPU-h must first run a ≤ 10-line pilot whose measured cost is recorded.
- **Budget ledger.** `results/budget.tsv` holds cumulative CPU-h against the §12 plan, checked at each milestone. The W7 300 CPU-h cap is read from this ledger, not estimated.
- **Agent concurrency.** Cap the number of sub-agents running at once (suggested 4) and give each task card a wall-clock limit.
- **Packing.** The §12 wall-clock figures assume perfect 16-way Gurobi packing on a shared host. Add a 1.5× buffer to the milestone dates, or the W7 stop rule ("24 h before deadline") will bite early.

### SF6. A pre-commit check that test labels never enter git (§14, `.gitignore`)
`results/` is committed (D1). If any agent writes a per-pair table with labels for test lines, for example a "full W1 benchmark" run before T1, and commits it, the labels stay readable in git history forever. Once committed, the lock cannot be restored.

The pre-commit hook therefore rejects:
- any committed table or parquet that contains test-tier ModelIDs together with a column named like `label|y|effect|dependency|chronos`, before the T1 ledger row exists;
- any file under `data/` or `papers/`;
- any file above 50 MB (D1).

### SF7. The ftINIT acceptance gate must be one-shot (W4, §10 W4 row, §13 M3)
The Hart 2015 lines (HCT116, DLD1, HeLa, a GBM line) belong to DepMap lineages, so one may be a test lineage. Tuning the RAVEN-python port until its AUPRC lands within ±0.05 would be label-tuning, and it would also game the gate.
- Port parameters are fixed to the repository defaults and committed before the run.
- The reference values (0.25–0.33) are read by script from `data/Human-GEM-main/data/testResults/gene-essential_summary.md`, not typed.
- One attempt is logged. A failure demotes F1, as the design already says.

### SF8. Overstatement controls in REPORT.md (§6.4, §9.4, §13 M9, §11 Selective reporting)
1. **Claims register.** Add `results/claims.tsv`: claim ID, sentence, evidence output SHA, contrast (C1–C5 or `exploratory`), the CI used, and the reviewer.
2. **Number check.** `scripts/verify_report.py` extracts every number in REPORT.md and checks it against:
   - `results/report/numbers.json`, produced by the report generator from the manifest outputs; or
   - a literature number tagged with a citation key.

   An unmatched number fails the M9 gate.
3. **Wording rules.**
   - "Transfers to unseen tissue" may be stated only with the two-stage lineage bootstrap CI (§6.4), not the line bootstrap.
   - Per-lineage results on a single test lineage are labelled descriptive.
   - "Gate not passed" (§9.1.4) and "method-dominated" (H4.2) flags are carried into REPORT.md automatically by the generator.
4. **Citations.** `scripts/verify_citations.py` resolves **every** DOI in LITERATURE.md and REPORT.md via Crossref or PubMed and fuzzy-matches the title and first author. This replaces the 20% citation spot-check. LITERATURE.md has 88 references; the check is cheap and catches the commonest LLM error.

### SF9. Review sign-off as files the gate scripts can check (§10, PROMPT "an agent that didn't produce a result reviews it")
Every `REVIEW.md` starts with front-matter:
- `reviewed_output_sha256`;
- `reviewed_commit`;
- `author_agent`;
- `reviewer_agent`;
- `verdict`;
- `independent_code_path: yes/no`.

`scripts/check_reviews.py` fails the milestone if any of the following holds:
- an output cited by a [x] or [!] tree node or by `claims.tsv` has no ACCEPT review whose SHA matches the current output (a re-run invalidates the review);
- the author and reviewer are the same agent.

The reviewer of W5 evidence works without label access, so it checks evidence and does not judge edits by their effect.

The failure tiers also need stating. W1 failure-mode test 3 (errors stratified by CN, expression and so on) and every error inspection before T1 are explicitly **training-tier only**, because §8 W1 does not say so.

---

## OPTIONAL

- **O1.** W7 ideas that name a specific lineage or its biology must name a training lineage. The PI flags any idea that targets a test-tier lineage's known biology, since tier identities are visible in `lineage_split.tsv`.
- **O2.** Q3: recompute a common-essential list from training lines with DepMap's rule and use it as a reporting sensitivity. This removes the last all-line label artefact from the tables.
- **O3.** Reviewers receive code and outputs, not the author's conclusion or a hypothesis to confirm. Where possible, use a different model from the author (Opus reviews Sonnet and vice versa).
- **O4.** At M4 and M6, give a fresh agent the task "find any path by which test or validation labels influence a training-tier output". It gets the code and the access log (MF1) as a red-team audit.
- **O5.** At each milestone, `scripts/reproduce.py --sample 10%` re-runs randomly chosen cheap manifest entries from a clean worktree and compares output SHAs. Solver outputs are compared on `ratio_r`.

## Items checked where no change is needed
- **Seeds (§5).** The fixed list, `SeedSequence.spawn()` and the unseeded-RNG grep are adequate. Extend the grep to `random_state=None` and to sklearn estimators without `random_state`.
- **Ties and rounding (§6.1).** Consistent with the baseline review M1, and covered by toy tests.
- **Missing-gene and NaN *label* policy (§2.3, §3).** It is explicit and counted. The *prediction* NaN policy is covered by MF7.
- **Silent drops.** `dropcounts.tsv` and the hard count assertions in §11 are adequate if the expected counts are put in `FROZEN.tsv` (MF5), not written as "≈".
- **DATA.md completeness.** One small gap: `AchillesNonessentialControls.csv` (arrived 14:49) is not in `data/checksums.tsv` yet. Re-run `checksum_data.py`.
