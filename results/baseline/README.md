# Baseline single-gene knockouts, generic Human1 and Human2

Predictions only. Nothing here is scored against DepMap CRISPR (labels and held-out lineages are fixed later in DESIGN.md);
the only DepMap file read is the header of `CRISPRGeneEffect.csv`, to report gene-ID coverage.

## Scoring rule: use `ratio_r`
**All downstream scoring (AUPRC, lethal calls, rankings) must use `ratio_r`, never `ratio`.**
`ratio` is the raw KO/WT value. About 2650 genes per model sit at ratio ≈ 1, spread over about 60 distinct floats in [1 − 1e-13, 1 + 3e-13] that change between runs (only 356/2887 and 441/2870 raw values were bitwise equal across two runs). That noise makes spurious ranks among true ties and moved a simulated AP by up to 0.018 (review M1).
`ratio_r` = `ratio` rounded to 6 decimals, then clipped to [0, 1] (no negative zero). Two full runs gave bitwise-identical `ratio_r` for every gene in both models (`scripts/compare_ko_runs.py`). Ties in `ratio_r` are true ties and must be scored as such; sklearn's `average_precision_score` treats equal scores as a single threshold.
Choosing the continuous essentiality score and the operating point is left to DESIGN.md.

## Status policy
- An infeasible KO LP is recorded as `ko_growth` 0 and `ratio`/`ratio_r` 0, and its status stays `infeasible`.
- Any other non-optimal KO or wild-type status (time limit, numeric, unbounded) raises `NonOptimalError`, so the run fails.
- If the wild type does not grow, ratios are NaN with status `wt_no_growth`. `essential_genes()` raises on NaN unless `nan="exclude"` is passed explicitly.
- `run_baseline_ko.py` also refuses to write outputs if any status is not optimal/infeasible or any ratio is NaN.
- In this run, all 5757 solves were optimal: there were no infeasible KOs and no NaN ratios.

## How produced
One command: `.venv/bin/python scripts/run_baseline_ko.py` (about 2.5 min, 15 workers). Prerequisite: `scripts/download_human1.py`.
Verification: `.venv/bin/python scripts/verify_knockouts.py` (about 3.5 min; writes `verification.tsv`).
Determinism check: `run_baseline_ko.py --out <dir>`, then `scripts/compare_ko_runs.py results/baseline <dir>`.
`scripts/download_human1.py` checks every file against pinned upstream SHA-256 values and fails on mismatch. All 8 files matched on a fresh download on 2026-10-26.

- Models: Human1 = Human-GEM **v1.19.0** (commit 022ed5c, last 1.x tag; v2.0.0 follows), `data/Human1/`; Human2 = Human-GEM **v2.1.1**, `data/Human-GEM-main/`.
- Medium: Ham's = inputs of the `GR` task in each release's own `metabolicTasks_Essential.txt` (derived from the model: `load_model` tags the release, and untagged models must pass `tasks_file`) (44 exchanges, uptake bound 1000, all other uptakes closed). The GR task is identical in both releases and all 44 inputs resolve to the same exchange IDs in both models.
- Objective: biomass `MAR13082` (same ID, same 13 species and coefficients in both models).
- KO: `cobra.Gene.knock_out()` inside a context (closes every reaction whose GPR evaluates False; isozymes buffer, complex subunits do not), then FBA. `ratio = ko_growth / wt_growth`. Infeasible KO LP = growth 0 (none occurred); non-optimal other status = NaN (none occurred).
- Solver: Gurobi 13 (WLS), 1 thread per worker; 16 concurrent WLS sessions worked with no licence error (tested; higher not probed on the shared licence). GLPK agrees with Gurobi (see verification).
- Code: `src/nextgem/knockouts.py` (medium, objective, solver, processes are parameters), `src/nextgem/genes.py`, `src/nextgem/model.py`.

## Files
- `{human1,human2}_ko.parquet`: gene (Ensembl), entrez, symbol, wt_growth, ko_growth, ratio (raw), **ratio_r** (score to use), status, solver, runtime (s per KO solve), map_status.
- `{human1,human2}_gene_mapping.tsv`: one row per model gene (no drops): symbol, entrez (primary), entrez_all, n_entrez, map_status, n_reactions, depmap_column (matched on Entrez), depmap_symbol_differs, depmap_symbol_candidate.
- `summary.tsv`, `model_differences.txt`, `verification.tsv`.

## Results (Ham's medium, biomass)
| | Human1 v1.19.0 | Human2 v2.1.1 |
|---|---|---|
| genes (all KO'd) | 2887 | 2870 |
| WT growth (a.u.) | 65.573 | 74.757 |
| lethal, ratio_r < 0.01 | 149 | 154 |
| ratio_r < 0.5 | 149 | 154 |
| partial (0.01 <= ratio_r < 1) | 74 | 91 |
| of which ratio_r >= 0.99 | 21 | 29 |
| ratio_r = 1 | 2664 | 2625 |
| distinct ratio_r values | 24 | 32 |
| genes whose KO closes no reaction (ratio 1 by construction) | 1801 | 1750 |
| solver status | all optimal | all optimal |
| KO wall time (15 workers) | 53 s | 56 s |
| summed solve time | 166 s | 209 s |

Lethal sets: 141 shared, 8 Human1-only, 13 Human2-only (by symbol). Loading each SBML takes about 15 s; each worker also unpickles the model once (about 14 s).
The ratio distribution is essentially bimodal: no gene lies in 0.01–0.5, and the smallest non-lethal ratio is 0.505. The partial genes (OXPHOS, ATP synthase, TYMP, SUCLG1) are mostly an artefact of unbounded uptake, where thymidine is used as a bulk carbon source. They should be treated as medium-dependent (review m1, for DESIGN.md).

## Gene mapping (model's own `genes.tsv`; Ensembl -> Entrez/symbol)
| | Human1 | Human2 |
|---|---|---|
| one_to_one | 2874 | 2857 |
| one_to_many (several Entrez; first listed used as primary, all in `entrez_all`) | 11 | 11 |
| unmapped (no Entrez: GK3P, SLC37A4 [ENSG00000281500]) | 2 | 2 |
| many_to_one | 0 | 0 |
| primary Entrez found in DepMap CRISPR columns | 2825 | 2807 |
| not in DepMap CRISPR (incl. the 2 unmapped) | 62 | 63 |
| DepMap symbol differs (renamed symbols, same Entrez) | 10 | 7 |

The 11 one-to-many genes (CYP2D6, PDE10A, USP4, USP9Y, GNMT, DHRS4, PDXK, AQP7, KLHL23, USP17L2, FNTB) carry extra LOC/pseudogene Entrez IDs; the first is the protein-coding gene. KLHL23 and USP17L2 are not in DepMap CRISPR. The two unmapped genes (GK3P pseudogene; SLC37A4 as the alternative-locus Ensembl ID ENSG00000281500, which has no Entrez in the table and no SLC37A4 column in DepMap CRISPR) get no Entrez and no symbol-based rescue; whether to curate them is a DESIGN.md decision. Gene IDs are matched on Entrez, never on symbol.

## Sanity notes
- Verification: 45/46 checks pass. The 1 known model issue is Human1 UMPS (below), and there are 0 unexpected failures.
  - The toy network gives exact results, including the status policy and `ratio_r`.
  - GPR closure equals an independent evaluator for all genes in both models.
  - FVA consistency holds, and CAD/UMPS are lethal on Human2.
  - In the genome-scale infeasible-KO test (biomass lb = 1), CAD KO is `infeasible` with ratio 0 in both models.
  - The task file is derived from each release, and unknown IDs are rejected.
  - Gurobi and GLPK agree on 200 random KOs (max |ratio diff| 4e-13).
  - Serial and parallel runs agree, bitwise on `ratio_r`.
- Independent review (`notes/review_baseline_ko.md`): an independent recomputation matched all 2887 Human1 and 2870 Human2 genes.
- **Human1 UMPS is not lethal (and DHODH is not either)**: Human1 v1.19.0 has thymidylate synthase MAR04644 reversible; Ham's contains thymidine, so dTMP -> dUMP -> deoxyuridine -> uracil -> UMP (UPRT, MAR04343) bypasses UMPS. Setting the MAR04644 lower bound to 0 (as in Human2) makes exactly UMPS and DHODH lethal and changes no other lethal call. It does change other numbers: WT growth drops by 0.25% (65.573 → 65.412), and 65 other ratios shift by ≤ 0.002 (OXPHOS/ATP synthase, TYMP; TYMS 0.9975 → 1). This was a diagnostic only; the baseline keeps the model as released. This is a Human1/Human2 difference to remember when interpreting pyrimidine genes.
- 1801 of 2887 Human1 and 1750 of 2870 Human2 genes close no reaction when knocked out (fully buffered by isozymes or all reactions already blocked/GPR-redundant), so most "predictions" are ratio = 1 by construction.
- Growth is in arbitrary units (uptake bound 1000); only ratios are meaningful.
