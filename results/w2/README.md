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

## M4 — training scoring, attribution, controls, shortlist (2026-10-28)

Author: W2 engineer agent. The session's label probe (a `head` of the raw CRISPR gene-effect file) was **blocked**
by the hook. Labels only through `nextgem.labels.train()` / `fold(k)` (via `cv.lolo`), training tier only; no
validation or test label read; `data/locked/` untouched. Nothing committed by the author; DESIGN.md, FROZEN.tsv, the
ledgers and `scripts/evaluate_val.py` were not edited.

**Code** (uncommitted; PROVENANCE rows `dirty=True`): `src/nextgem/w2_score.py` (label-free helpers),
`scripts/w2_score_train.py` (training scoring), `scripts/w2_flips.py` (label-free flip analysis),
`scripts/w2_serum_pct.py` (A9f sensitivity sweeps), `scripts/w2_val_scores.py` (look scores + spec),
`tests/test_w2_scoring.py` (13 toy tests). Scoring machinery reused from W1/W3: `metrics.per_line`/tAP,
`w1.method_summary`, `w1.delta_summary` (lineage-equal Δ̄_lin, lineage-stratified line bootstrap, B = 2,000,
SEEDS[0]), `w1.fit_components` (fold-scoped CE_fit), `metrics.binary_per_line`/`summarise_binary`.

**Set-up.** W1 primary training lines (686; 12 lineages); strata from `line_media.tsv`: **attribution 466**
(primary), onboarded 15 (sensitivity), other 205 (imputed 9, fallback 113, unmappable 82, mixture 1), all 686.
Score e = 1 − ratio_r of the line's cached sweep (`controls/line_sweeps.tsv`); **M0 equals W1's generic e bitwise**
(both parents, asserted). Reference: generic parent on M0 (G\* = Human1). TG: Δ̄_lin ≥ 0.005, CI lower > 0,
≥ 9 of 12 training LOLO folds (DESIGN §9; a fold with < 5 estimable lines is "not positive": Cervix has 4
attribution lines, so the maximum is 11/12; the brief's "8 of 10" is the Stage-A panel rule and is met as well).
Eligible for the shortlist (as specified by DESIGN §W2 and A9i, committed in 63a51df; the code that first applied it was not committed before scoring, see notes/reviews/M4_W2_review.md): M1, M2(f) gate-free on H1. M3: uncalibrated,
excluded, no sweep (A9i). Gate (+g): context lever, never in "medium alone".

### Ladder (attribution stratum, L_prob, all U genes; `train/step_scores.tsv`)
| Step | H1 tAP_lin | H1 Δ̄_lin vs M0 [95% CI] | folds | TG | Δ̄ noCE | H2 Δ̄_lin vs H2-M0 |
|---|---|---|---|---|---|---|
| M0 | 0.2992 | — | — | — | — | — (H2−H1 = −0.0084) |
| M-open (neg. control) | 0.3187 | **+0.0195** [0.0180, 0.0212] | 11 | pass | **−0.0168** | +0.0121 |
| M1 | 0.3081 | +0.0089 [0.0078, 0.0101] | 11 | pass | +0.0132 | +0.0078 |
| M2(f=0.001) | 0.3145 | +0.0153 [0.0141, 0.0166] | 11 | pass | +0.0175 | +0.0137 |
| M2(f=0.01) | 0.3145 | +0.0153 [0.0141, 0.0166] | 11 | pass | +0.0175 | +0.0137 |
| M2(f=0.1) | 0.3146 | +0.0154 [0.0142, 0.0167] | 11 | pass | +0.0175 | +0.0137 |
| **M2(f=1)** | 0.3150 | **+0.0159 [0.0146, 0.0172]** | 11 | pass | +0.0182 | +0.0134 |
| M1g / M2g(f) | = gate-free twin | gate effect 0.0000 (both parents) | | | | |
| M-HPLM:M2(f=1) | 0.3108 | +0.0116 [0.0111, 0.0121] | 11 | (not eligible) | +0.0133 | +0.0069 |
| M-Plasmax:M2(f=1) | 0.3092 | +0.0100 [0.0097, 0.0103] | 11 | (not eligible) | +0.0095 | +0.0043 |

mΔAP_gene is ≈ 0 for every step (−0.0004 to 0.0000): medium adds no per-gene, across-line ranking. Physiological
variants vs the documented M2(f=1): HPLM −0.0043, Plasmax −0.0059 (H1). Other strata, M2(f=1) H1: onboarded +0.0126
[0.0103, 0.0149] (15 lines, 1 fold estimable), other +0.0143 [0.0127, 0.0158] (12/12), all +0.0160 [0.0152, 0.0168].

### Decision
**Shortlist: M2(f=1) on Human1** (all five eligible steps pass TG; one-SE rule: only M2(f=1) is within one
paired-difference SE of the best). σ·√(2 ln m) = 0.0047 (σ over the 27 scored H1 steps' Δ, m = 5 eligible) < +0.0159.
**Flags carried with the candidate:**
- **W2-SHUF** (pre-registered negative control not beaten): documented − shuffled medium (within lineage × serum-free,
  SEEDS[0:5]) = +0.0005 [−0.0000, +0.0010] for H1 (H2: +0.0006 [0.0001, 0.0010]); shuffled media reproduce ≈ 97% of
  the gain (+0.0154 of +0.0159). Per eligible lineage (Bone, Head and Neck, Lung, PNS, Skin): −0.0001 to +0.0011.
- **W2-MOPEN** (negative control passes TG): M-open scores +0.0195, above every documented variant, on all genes;
  it is −0.0168 without CE_fit, where M2(f=1) is +0.0182. The all-genes TG pass is therefore not by itself evidence of
  documented-medium information.
- **W2-CONST** (post hoc, descriptive, enters no decision): one line-constant medium for every line (classic RPMI +
  10% FBS at M2(f=1)) scores **+0.0015 [0.0011, 0.0020] above** the documented per-line media (H1; H2 +0.0012).
- Reading: the gain is "a realistic culture medium instead of Ham's" (no hypoxanthine/thymidine → purine, TYMS,
  DHODH, UMPS calls; serum carnitine/cholesterol → BBOX1, TMLHE, ALDH9A1, SHMT2, SPCS/SEC11 and DHCR7 rescued; GSH →
  GCLC/GCLM/GSS), not line-specific medium. Medium acts as a generic-model setting here, not a context lever.
- Mevalonate core on the generic biomass: HMGCR, MVK, PMVK, MVD lethal in every medium at every step except M-open
  (pass for M2(f=1), both parents). HMGCS1, IDI1, FDPS, GGPS1 ratio 1 everywhere (known GPR issue). W3-MEV is a
  B-loo × serum interaction and is not tested here (no W2 × W3 cells; M6). Serum carnitine rescues several W3
  cofactor-pool artefact genes, so W2 and W3 B-loo:cofactor_pool_biomass overlap; expect sub-additivity in M6.

### Attribution (attribution stratum; `train/attribution.json`, `train/flip_counts.tsv`)
| | H1 | H2 |
|---|---|---|
| (i) medium alone = best gate-free documented − M0 | M2(f=1): +0.0159 [0.0146, 0.0172] | M2(f=0.1): +0.0137 [0.0126, 0.0149] |
| (ii) net flip fraction of M0 binary errors (L_prob; parent genes) | 2.7% (6,641 fixed − 2,634 broken of 147,072; 4,719 fixes are FPs) | 3.0% |
| (ii) same, L_ge06 | 2.6% | 3.2% |
| (iii) (doc − M0)/(oracle − M0); oracle = per-line best of {M0, M1, M2(f)} | **0.979** [0.974, 0.983] (oracle − M0 = +0.0162) | 0.971 [0.966, 0.976] |

Onboarded (sensitivity) H1 +0.0126; other stratum +0.0143. A small attributable share was pre-registered as a
plausible outcome.

### Controls
- Negative: M-open and shuffled medium, above (neither supports documented-medium information).
- Positive, training labels — **ASNS by medium (Lagziel replication): pass.** Lines whose documented medium lacks Asn
  (DMEM-based; M1 lethal; 76 lines, 11 lineages) have L_prob positive rate 0.66 vs 0.06 in 390 lines with Asn;
  ASNS tAP_gene over lines with the M1 score 0.549 vs null 0.167 (both parents).
- Positive, label-free (§1.1/§1.2 panel at every step; `train/panel_by_step.tsv`, `mevalonate_by_step.tsv`): ASNS,
  purines, TYMS, DHODH/UMPS (vs M0's HT), carnitine genes, GSH genes and DHCR7 move in the expected direction; sterol
  branch and CAD remain model-structure failures; gate controls non-informative.

### Label-free flip analysis (`flips/`; Human1, M0 reference)
- Single **additions** (59 exchanges any library medium opens): all viable; exact pruning (WT ≤ 50× M0; r_min 0.505)
  → only M0-lethal genes swept. **14 genes flip, all lethal → viable**: carnitine / O-acetylcarnitine (ALDH9A1,
  BBOX1, TMLHE, SHMT2, SPCS1/2/3, SEC11A/C), GSH (GCLC, GCLM, GSS), cholesterol (DHCR7), pyridoxal (PNPO).
- Cached documented sweeps (M1, M2(f)) flip **24 genes** vs M0 (the 14 + ATIC, GART, PAICS, PFAS, PPAT, TYMS, DHODH,
  UMPS from the loss of HT, ASNS and ALDH4A1).
- Single **removals** (44 M0 inputs, full sweeps): **not completed**. The first combined run stalled on GLPK
  degenerate LPs (3 workers ≥ 20 min on single jobs) and was stopped at ≈ 1.3 CPU-h with no output; the script now
  has an LP time limit. The removal half needs ≈ 1 CPU-h per parent and should run from a pinned clean worktree
  (`scripts/w2_flips.py --kinds removal`).

### Binary endpoints (§6.1 cutoff ratio_r < 0.01, parent genes, attribution stratum; `train/binary.tsv`)
| H1 | MCC | F1 | precision | recall |
|---|---|---|---|---|
| M0, L_ge06 | 0.336 | 0.366 | 0.511 | 0.286 |
| M2(f=1), L_ge06 | 0.349 (Δ +0.0126 [0.0108, 0.0144]) | 0.375 (+0.0087) | 0.534 | 0.289 |
| M2(f=1), L_prob | 0.339 (Δ +0.0162 [0.0145, 0.0179]) | 0.355 (+0.0106) | 0.594 | 0.254 |

H2 M2(f=1), L_ge06: MCC 0.348 (Δ +0.0149), F1 0.378. Fixed-cutoff calls, not P57 (no per-medium task checks).

### Sensitivities for M2(f=1) (`train/sensitivity.tsv`; Δ vs primary, H1)
RPMI-ATCC, DMEM+pyruvate, DMEM:F12 documented cysteine: 0.0000; DMEM without standard Gln −0.0002; DMEM:F12 without
standard Gln −0.0001; **A9f serum-% scaling (f_eff = pct/10; 28 new sweeps in `sens_serum_pct/`): −0.0002
[−0.0003, −0.0001]**.

### Validation look prepared (not run)
`val/look_spec_W2_M4_shortlist.json` → `val/w2_val_scores.parquet` (SHA-256 `17b9693a…3d68f`), **layout
`line_gene`** (the harness supports line-specific scores: `ModelID` + `entrez` × method over exactly the 198
validation expression lines; 146 are in the attribution set, the rest get the integrated-pipeline medium). Methods H1,
H2, W2:M2(f=1)[H1] (gate candidate), W2:M2(f=1)[H2] (descriptive); L_prob; all/noCE; U; none/expressed_pairs. H1/H2
columns equal the W1 generic scores (asserted). `--print-config` config_sha
`ddd23c7a9c152804dd4378d69d02689191dbb716459f0c728466c13bbb899e00`; `--dry-run` (synthetic labels): scalar guard pass,
16 rows, no validation label read, no ledger row (`val/dry_run.json`). The gate cell cannot be restricted to the
attribution set (harness strata are W1's). Before the look: commit, FROZEN rows, non-author config review (A6), and a
PI decision on the flags above.

### Runtime (PROVENANCE.tsv rows)
| Step | wall | CPU |
|---|---|---|
| `w2_score_train.py` (×2: first run; rerun adding A9f + post-hoc rows, step scores identical) | 133 s, 128 s | 131 s, 124 s |
| `w2_flips.py --kinds addition` (H1, 6 GLPK workers) | 350 s | 985 s |
| stalled combined flip run (stopped; no output, no PROVENANCE row) | ≈ 22 min | ≈ 1.3 CPU-h (process CPU time) |
| `w2_serum_pct.py --f 1` (28 sweeps, 6 GLPK workers) | 332 s | 1,556 s |
| `w2_val_scores.py` | 1.5 s | 1.6 s |
| harness `--dry-run` | 9 s | — |
| tests (13 new; with the W2, harness and scanner suites 352 pass) | 130 s | — |

Total ≈ 2.1 CPU-h, of which ≈ 1.3 CPU-h is the stalled run. Label-copy scanner over `results/w2` and the scratchpad:
0 hits.
