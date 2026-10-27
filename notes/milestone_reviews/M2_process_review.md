# M2 process review (generic baseline fixed, 2026-10-27)

Fresh-process reviewer. Read: CLAUDE.md, PROMPT.md, RESEARCH_TREE.md, LESSONS.md, the M1 review, M2_W1_review.md, git log, DESIGN.md §17 (A1-A5), skims of results/w1 and results/w4 READMEs. No label files opened; nothing re-run.

## 1. Stop doing

1. **Freeze ceremony for decisions the data barely move.** The W1 look chose G* on a gap of -0.0071 tAP (09b69a2); training already had H1 ahead in all 12 lineages (RESEARCH_TREE 3.3). It cost two review passes (M2_W1_review.md, ~5.5k words), amendment A5, a freeze (67f5f4f) and a snapshot commit of unreviewed W2/W4 code "so the tree is clean" (7010b07), which coupled unrelated workstreams. Do instead: reuse the reviewed look harness so a new look needs a config row, not a review or amendment. Pre-register "when arms are within SESOI on training, decide on training and spend the look on the confirmatory contrast only."
2. **Amendments for guard plumbing.** A5 amends test names (DESIGN.md:1138); A3 exists because an integer covariate looked like a gene column (:1135). Do instead: log guard-internal changes in results/test_retirements.tsv or a changelog when the SHAs of labels, metrics, splits and lists in FROZEN.tsv are unchanged. Reserve amendments for metrics, labels, splits, endpoints, a priori lists.
3. **Duplicated and superseded documents.** results/w4/README.md (318 lines) plus README_ratified_A4.md (256 lines, kept only to pin a hash); three a priori lists (20+40+52 KB); design_review r1-r3 (~400 KB); DESIGN.md at 23k words. Do instead: pin by hash line, delete the byte copy, move superseded rounds to notes/archive/ with an index, and point every agent brief at named DESIGN sections.
4. **Second written review rounds on mechanical fixes.** The M2 blocking list ended as three PI actions, all refused automatically by the script's pre-flight. Review once; verify fixes by dry run.
5. **Tree bloat.** RESEARCH_TREE.md (1,935 words) repeats seed history (D3-D5) and carries numeric detail (3.3). PROMPT wants a lean resume point; keep pointers.
6. **Further guard hardening** beyond new-look needs. Residual bypasses are accepted (results/m1_infra/README.md); the assert_scalar_only NOTE can wait until a look is added.

## 2. On track for PROMPT.md?

- **Compute budget is safe**: W1 about 3.5 CPU-h; Stage A projected 60 CPU-h, 3.5 h wall (results/w4/README.md §0). **Calendar and agent time is the risk**: 11 h wall (75fe8e3 13:30 to f679744 00:55), and about 25 of 40 commits are design, guard, ledger or review infrastructure. Delivered science: one baseline table.
- Open: medium (4), tasks (5), GPR (7), flux information (8), SBML resource (9), REPORT.md (10). No SBML export script or REPORT stub exists; the auto-research loop (W7, look sub-cap 20) has no code yet. These are the critical path.
- Scientific risk: gene-prior B2 (0.90) is 3x generic GEM (0.30), and 42% of positives lie outside the addressable set. A context GEM beating B2 is unlikely, so the auto-research loop must be allowed to combine GEM features with priors, and REPORT must say this plainly.
- Recommend: stub REPORT.md and SBML export now using the 30-line pilot models; run W2/W3/W5/W6 as parallel label-free work; allocate the remaining 39 looks per workstream now; cap each milestone at one design check, one code review, one process review.

## 3. Continue

1. **Non-author recomputation.** The reviewer reproduced every headline number (H1 0.2996, H2 0.2916, B2 0.9001, 42.3% outside) with its own code, and the label-free dry run found steps after the ledger `started` row and missing frozen inputs that toy tests could not (M2_W1_review.md A5).
2. **Mechanical enforcement**: ledger-before-labels, look cap, tripwire, frozen manifest. The look ran in 191 s, scalars only, with started/output/completed rows (a4c7a57, a9fd478, da1510b). No leakage incident since M1.
3. **Label-free parallel work** while gated steps wait: W4 pilots found repair restores 12/12 lines and priced Stage A before spending it.
4. **Pre-declared decision rules** (G* sign rule, SESOI frozen from training, CORDA dropped by its cost rule) and the per-run MANIFEST (SHAs, seeds, solver, wall/CPU).

## 4. LESSONS.md pruning (21 to 14 lines)

- Merge lines 8, 10, 11, 12, 16 into two: "exposure policy first, changes ratified by non-PI" and "secrets: commit-reveal, reveal after the using code is reviewed; no label aggregates in committed docs or scratch".
- Merge 5 and 9 (mechanisms, live probe, fail closed). Merge 18, 19, 20 (pre-freeze: label-free dry run on real IDs, final test names, commit all workstreams' code). Fold 6 into 2.
- Move 17 (hook and multi-line commit messages), 21 (SendMessage resume), 13 (blocked downloads) to CLAUDE.md; check 14 (swapped Model.csv columns) is in DATA.md, then keep one generic "check contents" line or drop it.
- Keep 1, 7, 15.
- **Add three rules:** (a) amend DESIGN only when a reportable number, label path or frozen SHA changes (A3, A5); (b) pin documents by hash, never keep byte copies (results/w4); (c) stub every PROMPT deliverable (SBML export, REPORT.md) at the first milestone with a model.
