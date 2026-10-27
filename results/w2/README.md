# W2 — medium and exchange bounds: library, floor, requirements, M3 pilots (label-free)

Author: W2 engineer agent, 2026-10-26; **updated 2026-10-27 for amendment A9** (DESIGN.md §17, PI decisions after
the blind two-agent extraction and the independent review `notes/reviews/M3_W2_review.md`, verdict GO pending
B1–B3, all applied here). **Label-free**: the label probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was blocked
at the start of both sessions; no label file and no `labels.py` was opened. Nothing here is committed by the author;
DESIGN.md, FROZEN.tsv and the ledgers were not edited.

## How to reproduce (provenance rows in `PROVENANCE.tsv`; pre-A9 runs in `MANIFEST.tsv`)
Every W2 script runs under `nextgem.run.provenance()` (DESIGN §11.4), via `w2_runtime.Manifest`, which writes one row
per run to `PROVENANCE.tsv`: argv, commit, dirty flag, full SHA-256 of inputs and outputs, solver and versions,
workers, wall/CPU, peak RSS. `allow_dirty=True` is used because the code is run uncommitted (W2 and W4 edits);
the rows record `dirty=True`. Gurobi is used only through `nextgem.gurobi_pool` (one slot per Gurobi process, ≤ 4).

| Step | Command | A9 run: wall / CPU | Solver, workers |
|---|---|---|---|
| 1 library, line assignment | `.venv/bin/python scripts/w2_build_library.py` | 33 s / 38 s | GLPK (no slot), 1 |
| 2 floor + req_c | `.venv/bin/python scripts/w2_medium_floor.py` | 149 s / 149 s (reproducibility rerun, see §4) | Gurobi, 1 slot |
| 3 prechecks, floor non-binding, calibration grid | `.venv/bin/python scripts/w2_prechecks.py --workers 4 [--solver glpk]` | 233 s / 855 s (Gurobi); final 338 s / 1,210 s (GLPK) | 4 |
| 4 KO sweeps (only sweeps not already cached) | `.venv/bin/python -u scripts/w2_ko_sweeps.py --workers 4 --verify 10 [--solver glpk]` | 196 new: 5,574 s / 21,841 s (Gurobi); 20 new: 513 s / 1,757 s (GLPK) | 4 |
| 4b GLPK = Gurobi cross-check | `... w2_ko_sweeps.py --solver glpk --recompute <ids>` | 143 s / 259 s | GLPK, 2 |
| 5 controls, per-line lookup | `.venv/bin/python scripts/w2_controls.py` | 28 s | none |
| tests | `NEXTGEM_NO_TEST_RECORD=1 .venv/bin/python -m pytest tests/test_w2_media.py -q` | 4 s (254 tests) | — |

The last increment (20 sweeps) and the final prechecks ran on GLPK because W4's Stage A held all 8 pool slots. GLPK
reproduces cached Gurobi sweeps **bitwise** on ratio_r (2 full sweeps, 5,757 genes; `ko/solver_crosscheck.tsv`),
consistent with the reviewer's 128-pair GLPK check. `ko/sweeps.tsv` records the solver of every sweep (452 Gurobi,
20 GLPK).

Code: `src/nextgem/media.py` (sheets, component/exchange maps, string parser, line assignment, bound scheme, gate,
shuffle), `src/nextgem/w2_engine.py` (sweep definitions, exact KO engine), `src/nextgem/w2_runtime.py` (pool-aware
parent loading, provenance), `scripts/w2_*.py`, `tests/test_w2_media.py`. `pypdf==5.1.0` is in `requirements.txt`
(the Plasmax table is a PDF).

## A9 changes applied (vs the reviewed version)
| A9 | Change | Effect |
|---|---|---|
| (a) | Thiamine HCl → 2 Cl⁻; Plasmax "alpha-Aminobutyrate" (symbol-font α kept); IMDM cystine·2HCl at MW 313 → **292 µM** (`compositions.tsv` column `correction`); mapping `w2-media-map-v2` | no bound changes (availability bounds) |
| (b) | `CYSTEINE_RULE = "design_literal"`: cysteine closed unless 2-ME, in every medium incl. HPLM/Plasmax. Documented free cysteine = sensitivity (`SENS:DMEM:F12-cysteine-documented`, `M-HPLM|cys-documented`, `M-Plasmax|cys-documented`) | DMEM:F12-based keys and M-HPLM/M-Plasmax repointed |
| (c) | Standard glutamine for plain DMEM (4 mM) / DMEM:F12 (2.5 mM) kept; sensitivities `SENS:DMEM-no-standard-Gln`, `SENS:DMEM:F12-no-standard-Gln` (new) | — |
| (d) | RPMI base = **classic RPMI 1640** (ATCC sheet without pyruvate and HEPES; glucose 11.1 mM and NaHCO3 23.8 mM from Cantor 2017 Table S1 "Classic synthetic media", Gibco 11875). The on-disk ATCC sheet = `RPMI-ATCC`, sensitivity `SENS:RPMI-ATCC`. Pyruvate added by the text becomes explicit (`RPMI+pyruvate`, `RPMI+alanine+pyruvate`) | RPMI-derived keys repointed. **req_c reference medium stays the on-disk ATCC sheet** (DESIGN: "RPMI-1640 … sheet on disk"), so the floor and req_c are unchanged |
| (e) | H2 floor at 100 × req_c accepted (25 binding cases, ≤ 1.87%; retinoate) | reported |
| (f) | `uptake_bounds(..., serum_pct_scaling=True)`: f_eff = f × serum% / 10 (defined and unit-tested; to be run at the TG-chosen f only) | not run |
| (g) | Uridine expectation is Plasmax-only (HPLM has none) | controls |
| (h) | ATP supplement → adenosine + Pi (MAR09254, MAR09072); ATP removed from the exchange map | 1 line, new key `RPMI+adenine+adenosine+hypoxanthine+thymidine` |
| (i) | M3 uncalibrated, excluded | unchanged |
| (j) | Controls: equal values → `non_informative`; sterol and CAD/uridine → `model_structure_failure (W3/W5)` | controls |
| fix | **2-ME bug**: with 2-ME the cysteine exchange now opens whenever cystine is present (before, it opened only if cysteine itself was listed, so "RPMI + 2-ME" equalled RPMI). Found while repointing; test added | `RPMI+cys(2-ME)` (2 lines): 20 new sweeps |
| F3 | Golden-table parser test: 16 representative strings with full expected parses | tests |

## 1. Media library (`media_library/`)
- `compositions.tsv` (**SHA-256 `1ae56834254ba1991d25461e1d950962e00f6cd23a5f2ce0c72832d2b2cd8326`**): RPMI (classic,
  derived) 40, RPMI-ATCC 42, DMEM 30, DMEM/F12 50, IMDM 41, HPLM 70, Plasmax 80 rows.
- `exchange_map.tsv` (**SHA-256 `b38111abb246e55955b4ca6b65869a3d8e0d84136679ab7c07f7549bf87a4c02`**,
  `w2-media-map-v2`): raw component → canonical metabolite(s) with stoichiometry → exchange id per parent, each id
  verified in Human1 v1.19.0 and Human2 v2.1.1. Two-agent record: `review/extraction_diff.tsv` (A9a).
- Unmapped nutrients (no exchange in either parent): RPMI and RPMI-ATCC 4-aminobenzoate, nitrate; DMEM nitrate;
  DMEM/F12 nitrate; IMDM nitrate, selenite; HPLM 4-aminobenzoate, nitrate; Plasmax Mn2+, methyl acetoacetate,
  nitrate, selenite, vanadate. Non-nutrients: HEPES, phenol red (`unmapped_by_medium.tsv`).
- `formulation_strings.tsv` (157 strings used by screened lines), `distinct_media.tsv` (**22** medium keys),
  `serum_model.tsv` (FBS reference concentrations still pending; they enter no bound).

## 2. Per-line medium assignment and attribution set (`line_media.tsv`)
Unchanged by A9 (independently re-derived by the reviewer with 0 differences): attribution **835** (RPMI 557, DMEM 180,
DMEM/F12 60, IMDM 38), onboarded sensitivity 35, mixtures on disk 6, fallback (sheet missing) 210, imputed 74,
unmappable 201. Imputed, fallback and unmappable lines get classic RPMI + 10% FBS (flagged). Gate: open 1,133,
closed 159, no expression 69 (open by default, flagged).

## 3. Uptake-bound scheme (`media.uptake_bounds`)
M-open: all exchanges 1000. M0: Ham's GR inputs at 1000. M1: mapped components present → 1000, absent → 0; O2, H2O
→ 1000; **cysteine closed unless 2-ME (then opened alongside cystine)**; floor at `floor_mult × req_c` (1000 if the
medium supplies it). M2(f): M1 + f × req_c for the six serum components (+ Fe/Zn/Cu at 1000), serum media only.
Optional f × pct/10. M3: M2 + glucose/glutamine/O2 at calibration rates. +g: Asp/Glu uptake 0 when the gate is
closed. M-HPLM/M-Plasmax: the composition under M2(f) + floor. Sensitivities: as in the A9 table.

## 4. Floor and requirement table (`floor/`)
Floor (both parents): Fe2+, α-/γ-tocopherol, biotin, cobalamin, linoleate, linolenate, lipoate, retinoate. req_c on
the on-disk RPMI (ATCC) sheet + floor at 1000. **req_c > 0 for every serum component** (none flagged).
The canonical tables are the reviewed ones. An A9 rerun of `w2_medium_floor.py` (`floor/rerun_check/`) reproduces
the floor set and every non-numeric field exactly; req_c agrees to 8.3e-12 relative (LP path noise). Because sweep
signatures hash bounds to 12 significant digits, the reviewed tables are kept, not overwritten.
`floor_mult.tsv`: H1 10×, H2 100× (25 cases still bind, ≤ 1.87%, accepted by A9e).

## 5. M3 pilots: prechecks (`prechecks/`)
- **WT precheck: all 472 distinct sweeps (558 ladder labels) pass**, WT growth 0.61–1.90 × M0.
- **M3 uncalibrated** (A9i): growth ≈ 12.3 × the O2 bound regardless of glucose/glutamine; window 0.01–0.06 h⁻¹ only at
  O2 ≈ 0.001–0.005. Excluded from selection.

## 6. Generic KO sweeps (`ko/`)
472 distinct sweeps (236 per parent). `ratio_r_{parent}.parquet`, `sweeps.tsv` (incl. solver), `sweep_labels.tsv`,
`m0_equality.tsv`, `solver_crosscheck.tsv`, **`label_changes_A9.tsv`** (every ladder label: old → new sweep).
- **M0 equals `results/baseline/{parent}_ko.parquet` bitwise** (2,887 / 2,870 genes), re-verified after every A9 run.
- A9 relabelling (558 labels): 134 unchanged, 272 repointed, 152 added, 60 removed. 220 labels point to sweeps computed
  in the A9 runs (216 new sweeps: 196 + 20); the others reuse sweeps that already existed. Sweeps no longer labelled
  (124 pre-A9 signatures) were dropped from the parquets.
- Cumulative KO runtime: pre-A9 380 sweeps 2.77 h / 11.0 CPU-h; A9 increments 1.69 h / 6.6 CPU-h (216 sweeps).
- Lethal counts (ratio_r < 0.01) over the 22 documented media: H1 M1 146–159, M2(f = 1) 136–149; H2 M1 151–162,
  M2(f = 1) 141–154; M0 149 / 154; M-open 91 / 97.
- NOTE (review N4): at M2(f = 0.01), ten serum-rescued genes sit at ratio_r 0.0099–0.0100, so binary calls at the
  0.01 threshold are knife-edge for that f.

## 7. Controls (`controls/`; column `status`)
| Control (label-free) | H1 | H2 |
|---|---|---|
| ASNS: DMEM (no Asn) lethal, RPMI ratio 1 (M1) | pass | pass |
| Mevalonate core lethal in every documented/physiological label | pass | pass |
| Purines rescued in HPLM, and in Plasmax, vs RPMI | pass, pass | pass, pass |
| DHODH/UMPS rescued by uridine, **Plasmax only** (HPLM has no uridine) | pass | pass |
| CAD rescued by uridine (Plasmax) | model_structure_failure (W3/W5) | same |
| Sterol branch rescued by serum cholesterol (4 f levels) | model_structure_failure (W3/W5) | same |
| GOT1, GOT2, SLC1A3 with the gate; SLC7A11 cystine vs cysteine; pyruvate (RPMI-ATCC vs classic RPMI) | non_informative (equal values) | same |
| M-open collapse | pass (91 < 149) | pass (97 < 154) |

Sterol/CAD diagnosis (review §2.6–2.7): the biomass cofactor pool requires vitamin D3 catabolites reachable only from
lathosterol via 7-dehydrocholesterol (cholesterol rescues only DHCR7), and the CAD KO is rescued by dihydroorotate/
orotate but not by uridine. Both are GEM-structure issues for W3/W5, not fixable through the medium.
Shuffled medium (SEEDS[0:5]): 44–46% of attribution lines change medium; mean per-line |Δratio_r| ≈ 1.1e-3.

## 8. Remaining open items
- FBS reference concentrations + citations in `serum_model.tsv` (documentation only).
- FROZEN.tsv rows for `compositions.tsv` / `exchange_map.tsv` (SHAs above): PI action.
- The serum-% sensitivity runs only once TG has chosen f.
- The DESIGN script name `medium_floor.py` is implemented as `w2_medium_floor.py` (accepted in review).
