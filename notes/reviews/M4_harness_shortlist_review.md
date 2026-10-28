# M4 review: evaluate_val.py "shortlist" look kind (A11f, under A6(2))

Reviewer: independent non-author agent (Claude Opus 5.5), 2026-10-28. Change under review: worktree
`agent-a2cf78e901e1ecc6f`, uncommitted. evaluate_val.py 1fccbe0c…a045; test_evaluate_val_shortlist.py 033dc58d…875a;
look_spec_W3_M4_shortlist.json b1b9d9af…1a5c; w3_val_scores.parquet aeff6564…ad70.

**Verdict: APPROVE WITH FIXES.**
- No leakage path found. Nothing writes or prints non-scalars. Validation labels cannot be read before the spec, scores-SHA, G* and frozen checks pass.
- W1 is unchanged apart from the harness SHA.
- **The spec as written cannot pass pre-flight (F1).** F2–F4 are cheap harness hardenings.

**Exposure.**
- The probe `head -c 50 …/CRISPRGeneEffect.csv` was **BLOCKED** by the hook.
- I read no label file, did not touch the locked store and ran no real look. I used only `--print-config`, `--print-frozen-candidates`, `--preflight-only` and `--dry-run`, all in a scratch clone.
- I saw synthetic dry-run scalars and the G* identity (Human1, already in the spec).
- Beyond those, I saw **two training-tier aggregates**: they sit in the `winners_curse_note` string of `results/w3/train/shortlist.json`, which I printed to check candidate names.
- I saw no validation value. The dry run's training positive rates were filtered from my view and the file deleted.

## Runs
- **Tests.**
  - New file: 40/40 pass. Existing suites (evaluate_val, metrics, metrics_binary, labels, frozen, run, tripwire): 160 pass.
  - 3 tripwire tests fail (sitecustomize_installed_and_current, links_to_real_raw_files_denied, access_log_records_denials). They fail identically with the **HEAD** evaluate_val.py, because the tripwire is anchored to /data/tmpA/nextgem (its PROJECT_ROOT, logs and data/depmap). The author's artefact claim is **confirmed**.
- **W1 equivalence** (my scratch tests).
  - HEAD and the new module run the full toy W1 look, with and without B6. Outputs are byte-identical: TSV, G* JSON, config JSON, config_sha and ledger rows.
  - `w1_config()` and `required_files()` are equal.
  - Real `--print-config` differs from results/w1/val_look_W1_M2_baseline_config.json only in code.evaluate_val.py (and so config_sha).
  - The `tap_grid` body is line-for-line the old loop.
- **Statistics by hand.** Full shortlist look on a toy with lineages of 2, 3 and 9 lines and distinct scores.
  - **Δ̄_lin** equals the lineage-mean of per-line sklearn-AP differences (1e-9), and is not the line-weighted mean.
  - **Bootstrap SD** matches the analytic √(Σ_k s_k²/n_k)/K within 8%. So the resampling unit is lines within lineage.
  - **CI and gate.** delta_ci equals the bootstrap_mean(B=2000, seed=11) percentiles, and eligible ⇔ lo > 0.
  - **One-SE rule.** diff_vs_best equals the hand paired difference. The simpler candidate was picked at 0.047 ≤ SE 0.054.
  - **Winner's curse.** term = σ√(2 ln 2), with σ the top candidate's analytic SD.
- **W3 look.**
  - `--print-config`: config_sha fd8a5f20…eddf.
  - `--dry-run`: 24 rows (6 methods × 2 genesets × 2 strata) and the scalar guard passes. H1/H2 equal the generic scores; m = 1 so term = 0. No ledger row and no file written. `--preflight-only` demands FROZEN and tracked status for the spec, scores, G* file and provenance files.

## Findings
**F1 (must fix; blocks the look). Stale provenance pins.**
- The spec pins `results/w3/train/shortlist.json` = 65a369… and `scripts/w3_score_train.py` = 73241….
- On main (efc4138) they are c54c653… and 20c9975…, so `pinned_problems()` (evaluate_val.py:274) fails pre-flight. I reproduced this.
- The parquet's MANIFEST row (results/w3/MANIFEST.tsv, 2026-10-27T23:20:08) has dirty=True and inputs={}. It predates the last w3_score_train.py edit.
- `scripts/w3_common.py`, which both W3 scripts import, is unpinned.
- Fix: (1) regenerate the parquet from committed code on a clean tree and confirm the SHA is still aeff6564…; (2) re-pin both files and add w3_common.py; (3) have the new spec SHA re-reviewed. Its candidates already match shortlist.json: one H1 gate candidate.

**F2 (should fix). The W7 sub-cap can be dodged.** check_spec (174-197) accepts any `workstream`. ledger.py:196 counts
the ≤ 20 W7 looks by that string. A spec in results/w7/val/ that declares "W3" escapes the cap. Fix: require
`W\d+` and that it equals the spec's `results/<ws>/val` directory.

**F3 (should fix). The primary label is not enforced.** decide_shortlist gates on `labels[0]` (502), and check_spec (188)
accepts any order. A spec with `["L_ge05","L_prob"]` would gate on L_ge05, but §8 and the G* rule use L_prob. Fix:
require `labels[0] == "L_prob"`.

**F4 (should fix). Some label-free errors surface only after `started` and burn the look id.**
- A method name equal to a gene symbol or lineage passes check_spec (182). It then fails assert_scalar_only in shortlist_scalars (490), after `started`.
- `B: 2000.0` or `seed: 11.0` pass `==` (191) and then break the bootstrap.
- A retry needs a new config_sha, which the ledger refuses. A6's dry run would catch both, but the harness does not enforce a dry run.
- Fix: run the guard's `bad()` on method names in prepare_shortlist, and require `type(...) is int` for B and seed.

**F5 (low).**
- `--print-frozen-candidates` crashes (FileNotFoundError, line 294) when a provenance file is missing. Reproduced.
- A missing or invalid G* file raises a bare error, not PreflightError (248).
- Duplicate look_ids across results/*/val are silently shadowed (162). Pre-flight fails safe, but they should be refused.
- `scores.path` is unconstrained (244). A path to the locked labels is denied by the tripwire (no labels.py scope) and by the column checks, so there is no leak. Restricting it to tracked files under results/ is cheap.
- The line_gene "binary" check (224) is cosmetic: label + 1e-9 passes. Tripwire plus provenance are the real guard.
- `--dry-run` prints training positive rates (659). Drop or round them so dry-run reviewers stay at synthetic-only exposure.

**F6 (PI to ratify; the code is right under each reading).**
- (a) Order of gate and selection: the harness applies the §8 gate first, then the §9 one-SE rule among eligible candidates; §10.1 (P1) applies the one-SE rule first, then the gate. They differ only if a simpler candidate within 1 SE fails the gate. This W3 look has one candidate, so it is unaffected.
- (b) "Simplest" means the spec's candidate order. The spec reviewer must check it against generic < medium < tasks < ….
- (c) m and σ: m is the number of gate candidates; σ is the bootstrap SD of the top candidate's Δ̄_lin vs G*.
  - m = 1 gives a term of 0. This is correct, because validation is independent of the training-stage selection.

**F7 (process).**
- The branch is based on 7b22d8b; main is at efc4138. All harness files and tests are identical, so the patch applies cleanly.
- The parquet is already committed on main with the same SHA, so drop the worktree copy.

## Leakage checklist (all pass)
- **Outputs.**
  - The TSV goes through assert_scalar_only (numeric whitelist; no ID, gene or lineage strings; row cap; n_bin=0). The decision JSON holds fixed keys with scalars, spec method names and the spec's rule text. The per-line `cells` never leave memory. Error messages carry names, not values.
- **No validation read in the non-look modes.** None of print-config, frozen-candidates, preflight-only or dry-run calls labels.validation(). The look reads it only after the committed `started` row, and the tripwire also denies it until evaluate_val.py is frozen.
- **Scores input.** Extra, missing or label-like columns, non-float columns, IDs, digits and duplicates are all refused before any label read, and tested. H1/H2 must equal the W1 generic scores, and the genes must equal U.
- **Ledger.**
  - The started row precedes any label read. A completed look id is refused, and the 40-look cap applies.
  - A copied spec is a new look: it needs its own FROZEN rows and counts against the cap.
- **Freeze and pins.**
  - The spec SHA and its fields are both in config_sha.
  - required_files covers the code, W1 context inputs, spec, scores, G* file and provenance files.
  - Post-freeze edits fail on the HEAD:FROZEN SHA, git cleanliness and the pins.

## Test quality
- **Coverage.** The tests target the right properties: refusals happen before label access, in every mode; known-answer Δ values (worst-case AP and constant AP); the dry run makes no validation read and writes no ledger row; the gate, one-SE rule and winner's curse, including m = 1; W1 config equality against the committed record.
- **Vacuous tests:** none.
- **Gaps:**
  - no test of pinned_problems or of the composition of required_files; no test of `main --dry-run`; no test of F2–F4 once fixed.
- **Test IDs.** Parametrised IDs are positional (change0…, over0…). Add `ids=` before freezing, because the test-run check treats a reordered ID as a disappearance.

## FROZEN rows to add (amendment A11; file rows, item = path)
Add these after the F1 regeneration and re-pin. If F2–F5 or the tests change, take each new SHA from
`--print-frozen-candidates` on main.

| path | sha256 | note |
|---|---|---|
| scripts/evaluate_val.py (replaces 8f85a666… row) | 1fccbe0c37d5a6ca3c8c95935ee431f03d81795b0c0c1f23c31f317f01bda045 | W3 validation look: code (A11f) |
| tests/test_evaluate_val_shortlist.py | 033dc58d29592903fc4e19824cdec2d2ae15569101bd9fbbf8a300995c02875a | known-answer test (DESIGN 11.3) |
| results/w3/val/look_spec_W3_M4_shortlist.json | new SHA after re-pin (b1b9d9af… is stale) | W3 validation look: input |
| results/w3/val/w3_val_scores.parquet | aeff6564aeede5859a7408924d63765e2367a47418aafd88b668d4bd065aad70 | W3 validation look: input |
| results/w1/val_look_W1_M2_baseline_gstar.json | 43663450ec5af961982f84971540928d2ca84ad7d0efa51eb40d77952106afa9 | W3 validation look: input (G*) |
| results/w3/train/shortlist.json | c54c653cf84b2d3aea3a64b75cd54ae862624e3aaa9dd7e686e0388a390087f9 | provenance |
| results/w3/train/mode_scores.tsv | b6f2b52922abbe2b497080959960e2a16047bc8b75b3f587b2ff7629e0e2739d | provenance |
| src/nextgem/w3.py | 77cb2b75723f3f7edd0fecb4feffabea3ba14d406ac4378c64fc415c860224ae | provenance |
| scripts/w3_score_train.py | 20c9975d1ea6cc1112e9aa3a7ab1ecef7fece91ffff47c20e65a07f9c52932ff | provenance |
| scripts/w3_val_dryrun.py | b237a86681c0bad7e34e6014310f2d6adab819b4118e5333b95d9eb89fd503b1 | provenance |
| scripts/w3_common.py (once pinned in the spec) | ae45d3dbecbb05092f1b347b3f919fb1afef5d00baab057f1886bc3c2f71be88 | provenance |

The W1 context inputs and labels/w1/metrics/run.py are already frozen at matching SHAs. Run `--preflight-only` and
`--dry-run` on main before the look.

## Re-check (2026-10-28; worktree rebased on main c0032b8; exposure: synthetic dry-run scalars only)
**Verdict: APPROVE.** I confirmed every reported SHA: evaluate_val.py deb82172…, test file fee116ee…, spec 2200dea8…, config_sha bf52cc8b…. All runs used a fresh scratch clone of c0032b8 plus the 3 worktree files.
- **F1 fixed.** All 6 provenance pins (now including w3_common.py) equal the committed files. The scores file 3b537c02… comes from a clean-tree MANIFEST row (7eb1f2a, dirty=False), and its kept columns equal the old file exactly. Pre-flight now reports only the missing FROZEN rows, the uncommitted files and my symlink artefacts; there are no pin problems.
- **F2–F5 fixed as recommended, each with a named test.** The spec checks now require: workstream = `results/<ws>/val`; labels[0] = L_prob; integer B and seed. Pre-start checks: unit-like method names are refused in prepare_shortlist (the old guard `bad()` extracted verbatim), and a missing or invalid G* file is a PreflightError. Path checks: duplicate look_ids are refused; scores.path must be `results/<ws>/val/*.parquet`; provenance must sit under results/, src/ or scripts/, with `..` and absolute paths refused. Output fixes: frozen-candidates prints `# MISSING` instead of crashing, and the dry run hides training rates. F6 is ratified in A11(g), and F7 is resolved: the worktree no longer carries a parquet copy.
- **Tests and runs.**
  - New file: 57/57 pass. The test IDs are now named; no recorded run contains the old IDs, so the rename is safe.
  - Existing suites: 177 pass, plus the same 3 test_tripwire environment failures.
  - My independent checks still pass against the new file: W1 toy look byte-identical to HEAD (with and without B6), and the hand statistics.
  - `--dry-run`: 16 rows, scalar guard passes, 0 ledger rows, no file written, m = 1, term = 0.
- **Regressions and leakage:** none found.
  - A stray duplicate spec now refuses every `--look`, W1 included. This fails closed and is acceptable.
  - The scores MANIFEST row still has `inputs={}`, so provenance rests on the SHA pins, which is adequate.

**FROZEN rows to add** (tab-separated; replace `<date>`; the first row replaces the 8f85a666… row):
```
scripts/evaluate_val.py	file	scripts/evaluate_val.py	deb821720cd6ddf51f6858fa2060b4b13253ba941e7e4cee9131cb5b8f65ef4c	<date>	A11	W3 validation look: code (A11f; M4 harness review) [replaces the 8f85a666… row]
tests/test_evaluate_val_shortlist.py	file	tests/test_evaluate_val_shortlist.py	fee116eebebb63c449935ac32e04a84be1b8e229fcacae1a7451e85e1d11784f	<date>	A11	known-answer test (DESIGN 11.3)
results/w3/val/look_spec_W3_M4_shortlist.json	file	results/w3/val/look_spec_W3_M4_shortlist.json	2200dea8d5cdc0d6abd1e56e4c39fe10be5a5419351a0d08259a18e072f6dbaf	<date>	A11	W3 validation look: input
results/w3/val/w3_val_scores.parquet	file	results/w3/val/w3_val_scores.parquet	3b537c027a0553c06622ac8443d0085a95c5015639341bcc8e6ee458dbb60020	<date>	A11	W3 validation look: input
results/w1/val_look_W1_M2_baseline_gstar.json	file	results/w1/val_look_W1_M2_baseline_gstar.json	43663450ec5af961982f84971540928d2ca84ad7d0efa51eb40d77952106afa9	<date>	A11	W3 validation look: input (G*)
results/w3/train/shortlist.json	file	results/w3/train/shortlist.json	c54c653cf84b2d3aea3a64b75cd54ae862624e3aaa9dd7e686e0388a390087f9	<date>	A11	W3 validation look: provenance
results/w3/train/mode_scores.tsv	file	results/w3/train/mode_scores.tsv	b6f2b52922abbe2b497080959960e2a16047bc8b75b3f587b2ff7629e0e2739d	<date>	A11	W3 validation look: provenance
src/nextgem/w3.py	file	src/nextgem/w3.py	77cb2b75723f3f7edd0fecb4feffabea3ba14d406ac4378c64fc415c860224ae	<date>	A11	W3 validation look: provenance
scripts/w3_score_train.py	file	scripts/w3_score_train.py	20c9975d1ea6cc1112e9aa3a7ab1ecef7fece91ffff47c20e65a07f9c52932ff	<date>	A11	W3 validation look: provenance
scripts/w3_val_dryrun.py	file	scripts/w3_val_dryrun.py	b237a86681c0bad7e34e6014310f2d6adab819b4118e5333b95d9eb89fd503b1	<date>	A11	W3 validation look: provenance
scripts/w3_common.py	file	scripts/w3_common.py	ae45d3dbecbb05092f1b347b3f919fb1afef5d00baab057f1886bc3c2f71be88	<date>	A11	W3 validation look: provenance
```
Before the look: commit the harness, test and spec, add the rows, then run `--preflight-only` and `--dry-run` on main.
