# Review round 2 — PI decisions (2026-10-26)

Reviews: `external_pi.md`, `statistics.md`, `ai_failure_modes.md`, `algorithms_compute.md`, `data.md`, `prior_work.md`, `domain_metabolism.md` (this folder). All seven verdicts: revise minimally, then freeze. Every round-2 must-fix is accepted except where modified below. Should-fix items are at the author's discretion; record each one in the change log.

## Splits (ext-PI N1, stats NM1, AI N1)
1. **Retire seed 20261026.** The round-1 statistics review published its draw, so v1's test tier is known.
   - The new seed is supplied by the user **after** the freeze commit (commit-reveal). It is not a commit SHA, which an agent can game.
   - `make_splits.py` refuses 20261026. It records the seed source, the freeze commit and every redraw.
   - Correct §4.1 (no "unseen" claim for the v1 tier).
   - Update the §11.5 exposure register with:
     - the round-1 statistics review's label summaries on the old tier;
     - the design author's `head` print;
     - the all-line label copy `ge_U.parquet` and other derived scratch files, written by a round-1 reviewer and deleted by the PI on 2026-10-26;
     - the PI's header-only read while testing the hook.

## Gatekeeping and baselines
2. **C6 goes last:** C1 → C1b → C3 → C4 → C5 → C6.
   - This follows the prior-work placement over the statistics suggestion of "after C1b". C6 compares against another method, so a C6 failure must not block the core scientific claims.
   - B7g is the only confirmatory comparator for C6; B7t is an estimate.
3. **Rename B7** to "Human-GEM repository ftINIT protocol (Python, 1+1)". It is not the published MATLAB protocol.
   - Run it on all validation and test lines. Training lines are deferred until P1 is known, and the F1 training-gate gap must be stated.
4. **B8:** report both tie-break directions.
5. **SESOI:** restate the basis for 0.01, and add C1b's SESOI with its own basis.

## Solver rule (ext-PI N3, compute NMF2)
6. **Time limits.** For every MILP build (B7, F1, any other ftINIT/MILP arm), a feasible time-limited incumbent is **accepted** and logged with its gap, matching the repository pipeline.
   - Only no-incumbent, infeasible, licence or solver errors count as a failed build.
   - A failed build gets the no-information score in every arm, B7 included.
   - Licence and solver errors still abort the run (round-1 decision 20).

## Medium (domain MF1, compute NMF1, ext-PI N2, data N1)
7. **Requirement definition.**
   - Micronutrients/vitamins required by biomass: a fixed, non-limiting floor, not scaled by f.
   - Serum components the model can synthesise (cholesterol, fatty acids, carnitine, …): req = pFBA uptake at ≥99% of WT growth with the component open; then apply f.
   - True auxotrophs: FVA-min uptake.
   - Define the requirement on a stated reference medium.
8. **Transporter gate:** a separate on/off option. "Medium alone" is measured on gate-free variants. In attribution the gate counts as **context**, not medium.
9. **Media sheets.** The medium-attribution set is the lines whose base-medium sheet is on disk at M3 (currently RPMI, DMEM, DMEM/F12, IMDM; ~870 lines). The others are flagged as fallback.
   - Recompute the two-media lineage counts on that basis via `design_numbers.py`.
   - The user has been asked for MEM/EMEM, F12, McCoy's 5A and L-15.
10. **F1 repair.** F1 models scored on any medium other than Ham's go through the same repair step.

## Genomic constraint (domain MF2, data should-fix)
11. **Loss** = homozygous damaging (matrix code 2), or heterozygous damaging with LOH, low CN, or TPM < 1.
    - Make the drop-and-log rule order-independent: test each loss alone on the parent; then drop every loss that kills growth.

## Compute (compute NMF4–5)
12. **Freeze an explicit Stage-A build table** (rows: config × line set × seeds; columns: builds, expected CPU-h). The total must match §12.
13. **Workers:** start at 12 MILP + 16 LP. Final counts are set from ftINIT peak memory measured at M3, keeping ≥20 GB free on the shared host. The MILP lane starts at M2.

## Agent guards (AI N2–N4)
14. **Content-based label-copy scanner** (`scripts/scan_label_copies.py`). It finds files anywhere in the project, results and scratch whose numeric content matches label matrices (e.g. hashed row fingerprints).
    - Run it at every milestone and before every push.
15. **Tripwire.** Install the tripwire for every interpreter launch of `.venv` via `sitecustomize`. It covers pandas, pyarrow and plain `open`.
    - The residual risk of using a different interpreter is documented. The PreToolUse hook (now fail-closed; also guards `labels_val`, `data/locked/` and recursive data searches) blocks direct agent reads.
    - Every new agent session first runs a probe proving the hook is active.
16. **Sealed outputs.** T1 sealed outputs go to `data/locked/` (hook-protected, git-ignored). Only their hashes are committed.
    - A pre-commit guard stays on until the T2 unseal.

## Approvals
17. Install `openpyxl` in `.venv` and pin it in `requirements.txt`. Guide maps are approved at low priority; no claim may depend on them.
