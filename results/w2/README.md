# W2 — medium and exchange bounds: library, floor, requirements, M3 pilots (label-free)

Author: W2 engineer agent, 2026-10-26. **Label-free**: no label-bearing file and no `labels.py` was opened; the
liveness probe `head -c 20 data/depmap/CRISPRGeneEffect.csv` was blocked by the label guard before any work.
Design: DESIGN.md v3 (FROZEN) W2, §10, §12 row W2, §13–14 (M3 pilots). Nothing here is committed; nothing is frozen
(FROZEN.tsv, DESIGN.md and the record files were not edited).

## How to reproduce (in order; runtimes from `MANIFEST.tsv`)
| Step | Command | Wall | CPU | Workers |
|---|---|---|---|---|
| 1 library, line assignment | `.venv/bin/python scripts/w2_build_library.py` | 41 s | 44 s | 1 (GLPK, no Gurobi) |
| 2 floor + req_c | `.venv/bin/python scripts/w2_medium_floor.py` | 199 s | 199 s | 1 Gurobi |
| 3 prechecks, floor non-binding, calibration | `.venv/bin/python scripts/w2_prechecks.py --workers 4` | 231 s | 826 s | 4 Gurobi |
| 4 KO sweeps (380) | `.venv/bin/python -u scripts/w2_ko_sweeps.py --workers 4 --verify 10` | 9,974 s (2.77 h) | 39,658 s (11.0 CPU-h) | 4 Gurobi |
| 5 controls, per-line lookup | `.venv/bin/python scripts/w2_controls.py` | 20 s | 21 s | 0 |
| tests | `NEXTGEM_NO_TEST_RECORD=1 .venv/bin/python -m pytest tests/test_w2_media.py -q` | 5 s | — | — |

Code: `src/nextgem/media.py` (sheets, component/exchange maps, string parser, line assignment, bound scheme, gate,
shuffle), `src/nextgem/w2_engine.py` (sweep definitions, KO engine), `src/nextgem/w2_runtime.py` (parent loading,
runtime manifest), `scripts/w2_*.py`, `tests/test_w2_media.py` (232 tests: every one of the 212 distinct DepMap
medium strings, toy nutrient removal, bound scheme, assignment counts = design numbers, req_c ≠ 0).
`pypdf==5.1.0` was installed into `.venv` (Plasmax table is a PDF) and added to `requirements.txt`.
`nextgem.run.provenance()` (§11.4) does not exist; each script appends wall/CPU/RSS/solver/workers/output SHA to
`MANIFEST.tsv` instead (`code_dirty=True`: the W2 code is uncommitted).

## 1. Media library (`media_library/`)
- `compositions.tsv`: every row of the six formulations on disk (RPMI ATCC 42 rows, DMEM high-glucose 30,
  DMEM/F12 50, IMDM 41, HPLM 70 = Cantor 2017 Table S1 sheet "HPLM formulation", Plasmax 80 = Vande Voorde 2019
  Table S1, Plasmax column, parsed from the PDF text).
- `exchange_map.tsv` (`w2-media-map-v1`): raw component → canonical metabolite(s) with stoichiometry (salts split
  into ions) → exchange id in **each parent**, with the parent's metabolite name; every mapped id was verified to be
  an exchange with the expected metabolite in Human1 v1.19.0 and Human2 v2.1.1 (same ids in both). Non-literal
  mappings carry a note (cyanocobalamin → aquacob(III)alamin, NH4+ → NH3, pyroglutamate → 5-oxoproline, ...).
- **Unmapped components per medium** (`unmapped_by_medium.tsv`; no exchange in either parent, or non-nutrient):

  | Medium | Nutrients with no exchange | Non-nutrients |
  |---|---|---|
  | RPMI | 4-aminobenzoate, nitrate | HEPES, phenol red |
  | DMEM | nitrate | — |
  | DMEM/F12 | nitrate | phenol red |
  | IMDM | nitrate, **selenite** | HEPES |
  | HPLM | 4-aminobenzoate, nitrate | phenol red |
  | Plasmax | Mn2+, methyl acetoacetate, nitrate, **selenite**, vanadate | phenol red |

  Supplements: selenite/selenium (ITS, N2, B27), BrdU and 2-mercaptoethanol have no exchange (2-ME only switches the
  cysteine rule); insulin, transferrin, BSA, growth factors, antibiotics, kits are non-nutrients.
- `formulation_strings.tsv`: all 157 strings used by screened lines (+ the imputed text) with base, mixture
  fractions, serum %, parsed supplements and flags. The parser covers all 212 distinct strings in
  ModelCondition/Model.csv (unit-tested one by one).
- `distinct_media.tsv`: 20 bound-relevant medium keys (base sheet + mapped additions, serum flag).
- `serum_model.tsv`: f-scaled palmitate, stearate, oleate, cholesterol, carnitine, ethanolamine; trace metals
  Fe2+/Fe3+/Zn2+/Cu2+ at 1000; selenium unrepresentable; no nucleosides/hypoxanthine. **FBS reference
  concentrations and citations are PENDING** (domain agent, DESIGN) — they do not enter any bound.

## 2. Per-line medium assignment and attribution set (`line_media.tsv`, `library_summary.json`)
Tier logic identical to `design_numbers.py` (asserted in tests): screen 1,225, screen_mixed 6, onboarded 56,
imputed 74. Content assertions check the Model.csv / ModelCondition column swap.

| Stratum | Lines | Medium used |
|---|---|---|
| attribution (screen tier, sheet on disk) | **835** (RPMI 557, DMEM 180, DMEM/F12 60, IMDM 38) = [dn] | documented |
| onboarded_sensitivity (onboarded tier, sheet on disk) | 35 | documented |
| mixture_on_disk (IMDM:RPMI, DMEM:RPMI) | 6 | linear blend, flagged |
| fallback_sheet_missing (MEM/EMEM, F12, McCoy's, L-15; + 12 mixtures with a missing sheet) | 210 | RPMI + 10% FBS, flagged |
| imputed | 74 | RPMI + 10% FBS, flagged |
| unmappable (organoid/stem/other + 6 screen_mixed) | 201 | RPMI + 10% FBS, flagged (own stratum) |

Gate (+g, SLC1A1/2/3 log2(TPM+1) ≥ 1, default expression rows): open 1,133, closed 159, no expression 69 (gate left
open, flagged). MF-ID cross-check (`mf_crosscheck.tsv`): one non-organoid disagreement (MF-024: L-15 vs NSA, 1 line).
`shuffled_media.tsv`: medium keys permuted within lineage × serum-free among attribution lines, SEEDS[0:5].

## 3. Uptake-bound scheme (`media.uptake_bounds`)
M-open: all exchanges 1000. M0: Ham's GR inputs 1000. M1: every mapped component present → 1000, absent → 0;
O2, H2O → 1000; cystine (cysteine only if 2-ME or the sheet lists free cysteine); floor at `floor_mult × req_c`
(1000 if the medium supplies it). M2(f): M1 + f × req_c for the six serum components (1000 if the medium has them)
+ trace metals 1000, serum-containing media only. M3: M2 + glucose/glutamine/O2 at calibration rates. +g: Asp/Glu
uptake 0 when the gate is closed (a context lever; never in "medium alone"). M-HPLM/M-Plasmax: the composition as a
serum-containing medium under M2(f) + floor. Sensitivities: RPMI without pyruvate, DMEM + pyruvate, DMEM without the
standard glutamine default, DMEM/F12 with cysteine closed (DESIGN-literal).

## 4. Floor and requirement table (`floor/`)
Drop-one on M0 (both parents): 25 essential inputs. Floor = union over the four sheets of the essential inputs they
lack (greedy-pruned): **Fe2+, α-tocopherol, γ-tocopherol, biotin, cobalamin, linoleate, linolenate, lipoate,
retinoate** (same for H1 and H2; cholesterol never). Beyond the compute probe's list: Fe2+ (RPMI and IMDM carry no
iron; Fe3+ in DMEM substitutes), biotin and cobalamin (absent from DMEM high-glucose). Pyridoxine is not needed in
IMDM (pyridoxal).
`requirements.tsv` on the reference medium (RPMI sheet + floor at 1000): floor req_c = FVA-min uptake (cofactor-pool
items 0.00192, linoleate 5.80, linolenate 11.97, retinoate 0.0268 at WT 46.4); serum req_c = pFBA uptake at ≥ 99%
WT with c open (H2: palmitate 5.30, stearate 11.18, oleate 8.01, cholesterol 1.30, carnitine 0.0019, ethanolamine
1.76; H1 within 2%). **req_c ≠ 0 for every serum component in both parents (none flagged).**

## 5. M3 pilots: prechecks (`prechecks/`)
- **WT precheck: all 380 distinct sweeps (466 ladder labels) pass**, WT growth 0.61–1.90 × M0 (M0 = 65.57 H1,
  74.76 H2). Every documented medium, HPLM, Plasmax and the sensitivity media grow under M1 and all M2(f), with and
  without the gate. Nothing is blocked.
- **Floor non-binding**: H1 passes at 10 × req_c in all 128 medium × variant cases. H2 binds (rel. 0.26–1.9%) in
  DMEM/F12 (± serum), RPMI + HT/adenine/ATP and the DMEM/F12 sensitivity → raised to **100 × req_c for H2** and
  logged (`floor/floor_mult.tsv`); the same 20 cases still bind at 100× (≤ 1.9%): retinoate and Fe2+ are consumed
  beyond the biomass need when growth exceeds the reference (a model artefact; floor bounds are capped at 1000).
- **Calibration: M3 UNCALIBRATED.** `data/media/calibration.tsv` does not exist. The grid (RPMI M2, f ∈ {0.01, 1})
  shows WT growth ≈ 12.3 × O2 bound (h⁻¹ per mmol O2 gDW⁻¹ h⁻¹), independent of glucose and glutamine, because all
  other medium components stay at 1000 and fuel growth: the 0.01–0.06 h⁻¹ window is reached only at O2 ≈
  0.001–0.005, at any glucose/glutamine. As specified (only glucose, glutamine, O2 bounded), M3 cannot be calibrated
  at physiological O2 rates; per DESIGN it is reported uncalibrated and excluded from selection.
- M-open growth 124.9 (H1) / 138.5 (H2), ≈ 1.9 × M0.

## 6. Generic KO sweeps (`ko/`)
380 distinct sweeps (190 per parent: M-open, M0, 20 documented media × {M1, M2(f) × 4} × gate where the medium has
Asp/Glu, M-HPLM/M-Plasmax × M2(f) × gate, 4 sensitivity media). `ratio_r_{parent}.parquet` (sweeps × all model
genes, float64), `sweeps.tsv`, `sweep_labels.tsv` (466 labels → sweep), `m0_equality.tsv`.
- Engine = `nextgem.knockouts` semantics with three exact shortcuts (no-closure, zero-flux, identical closure set).
  **M0 sweep equals `results/baseline/{parent}_ko.parquet` bitwise for all 2,887 / 2,870 genes**; 10 shortcut genes
  re-solved per sweep: max |ratio − 1| = 1.8e-11.
- Per sweep: 191–310 LPs (≈ 104 s CPU); 0 non-optimal statuses; 15 infeasible KO LPs in total (recorded as 0).
- **Runtime: 2.77 h wall, 11.0 CPU-h on 4 Gurobi workers** (≈ 1.7 CPU-min per sweep, i.e. per distinct medium ×
  variant × gate × parent; the §13 estimate was 3.5 CPU-min). A first run was killed by the agent's own 12,000 s
  timeout before writing output; the script now checkpoints each sweep (`results/w2/tmp/ko_parts/`, deleted after
  assembly) and resumes.
- Lethal counts (ratio_r < 0.01) over the 20 documented media: H1 M1 146–159 (median 154), M2(f ≥ 0.1) 136–149;
  H2 M1 151–162 (158), M2(f ≥ 0.1) 141–154. M0: 149 / 154. M-open: 91 / 97.

## 7. Controls (`controls/`)
`positive_controls.tsv`, `negative_controls.tsv`, `panel_ratio_r.tsv` (ladder label × mevalonate core and a priori
§1.2 genes), `line_sweeps.tsv` (ModelID × parent × ladder step → sweep_id; the +g steps point to the gated sweep
only for gate-closed lines), `ladder_summary.tsv`.

| Control (label-free) | H1 | H2 |
|---|---|---|
| ASNS: DMEM (no Asn) lethal, RPMI ratio 1 (M1) | pass | pass |
| Mevalonate core (HMGCR, MVK, PMVK, MVD) lethal in every documented/physiological label (201) | pass | pass |
| Purines (PPAT, GART, PFAS, PAICS, ATIC) rescued in HPLM and Plasmax vs RPMI | pass (1 vs 0) | pass |
| DHODH, UMPS rescued by uridine in Plasmax (CAD stays lethal) | pass | pass |
| **Sterol branch: lethal serum-free, rescued by serum cholesterol** | **FAIL** | **FAIL** |
| GOT1/GOT2/SLC1A3 with the gate; SLC7A11 under cystine; pyruvate relaxes GOT/MDH/OXPHOS | no effect (ratio 1 either way; non-informative) | same |
| M-open collapse (lethal M-open < M0) | 91 vs 149 pass | 97 vs 154 pass |

- **Sterol failure**: serum-free, 7/12 (H1) and 8/12 (H2) sterol genes are lethal. Serum cholesterol (even at f = 1)
  rescues **only DHCR7** (ratio ≈ 0.6–0.7 at f = 1, ≈ f × 0.6–0.7 below), and FDFT1, SQLE, LSS, NSDHL, EBP, DHCR24
  (+ HSD17B7 in H2) stay lethal. The parents therefore need a sterol-branch product other than cholesterol, so the
  §1.2 "serum LDL rescues the sterol branch" expectation is not representable. This is a model issue for W3/W5,
  not a medium issue. CYP51A1, TM7SF2, MSMO1, SC5D are ratio 1 everywhere (GPR/bypass).
- The gate changes no GOT1/GOT2/SLC1A3 ratio (Asp is synthesised either way), and SLC7A11 is 1 under M1 cystine (an
  alternative cystine route in the model). These positive controls are non-informative as built.
- Shuffled medium (SEEDS[0:5], within lineage × serum-free): 44–46% of attribution lines change medium; mean per-line
  |Δratio_r| ≈ 1.1e-3 (max 6e-3). Medium choice moves few genes, so the label-based shuffled-medium null will be
  tight. The label-based evaluation is not run here.

## 8. DESIGN ambiguities and deviations (for the PI)
1. **Cysteine rule**: DESIGN M1 says cysteine closed unless 2-ME; the DMEM/F12 sheet documents 0.1 mM free cysteine.
   Default keeps documented cysteine; the DESIGN-literal variant is the sensitivity `SENS:DMEM:F12-cysteine-closed`.
2. **Glutamine-free sheets on disk** (DMEM, DMEM/F12 "NoGlutamine"): plain "DMEM"/"DMEM:F12" texts get the standard
   formulation's glutamine (flag `gln_standard_default`; sensitivity without it).
3. **The RPMI sheet is ATCC-modified** (4.5 g/L glucose, HEPES, 1 mM pyruvate; Gibco 11875 has no pyruvate):
   used as on disk; sensitivity `SENS:RPMI-no-pyruvate`.
4. **Floor contains biotin, cobalamin and Fe2+** (required by biomass, absent from DMEM or RPMI/IMDM). This removes
   the biotin/B12 contrasts of the a priori panel (§1.2: MMAA/MMAB/MMUT, ACACA/PC/PCCA/PCCB in DMEM).
5. **H2 floor still binds at 100 × req_c** in 4 media (≤ 1.9%); DESIGN gives no step beyond 100×.
6. **M3 infeasible as specified** (above); `calibration.tsv` missing.
7. **Serum % is parsed but unused**: M2 bounds are f × req_c regardless of 5/10/20% serum (DESIGN gives no % scaling).
8. **Serum hemin** (a priori §1.2 heme genes) and serum nucleosides are not in the DESIGN serum model; heme-pathway
   genes cannot be rescued by any variant.
9. **HPLM has no uridine** in Cantor 2017 Table S1 (a priori lists say it does); Plasmax has 3 µM uridine.
10. **Two-agent extraction** (DESIGN: compositions and exchange map by two independent agents, script-diffed, frozen
   in FROZEN.tsv, stored in `data/media/`) not done: single extraction here, written to `results/w2/media_library/`
   because `data/` is read-only. FBS reference concentrations pending.
11. Unmappable and screen_mixed lines get the imputed medium (flagged); DESIGN only says "own stratum".
12. Two-sheet mixtures: on-disk mixtures (6 lines) are built but not in the attribution set (same as
    design_numbers); mixtures with a missing sheet (12) are `fallback_sheet_missing` (design_numbers' 198 counts
    single-sheet bases only).
13. Horse serum / Nu-Serum / FCS count as serum (design_numbers counts "FBS" only: 160 vs 159 serum-free screened).
14. Script name: DESIGN says `scripts/medium_floor.py`; implemented as `scripts/w2_medium_floor.py`.
