# M3 / W2 independent review: media library, bounds, floor, KO sweeps, controls

**Reviewer:** independent review agent (not the author). **Date:** 2026-10-27.
**Code reviewed (uncommitted, as on disk):** `src/nextgem/media.py`, `src/nextgem/w2_engine.py`,
`src/nextgem/w2_runtime.py`, `scripts/w2_build_library.py`, `scripts/w2_medium_floor.py`, `scripts/w2_prechecks.py`,
`scripts/w2_ko_sweeps.py`, `scripts/w2_controls.py`, `tests/test_w2_media.py`; outputs under `results/w2/`
(starting from `results/w2/README.md`). HEAD 6bf19d2.
**Spec:** DESIGN.md §9 W2 (lines 452–528), §10, §12 row W2 (line 995).
**Review outputs:** `results/w2/blind/` (Part 1, blind second extraction) and `results/w2/review/` (diff, KO
verification, line-assignment re-derivation, sterol trace, diagnostics). All are label-free.

## 0. Session probe, order of work and label exposure
- **Probe.** `head -c 50 data/depmap/CRISPRGeneEffect.csv` was **BLOCKED** by the label guard ("Blocked by label
  guard: shell command references raw label files ...").
- **Blind order.** Part 1 (`results/w2/blind/`) was written before I opened any W2 output or W2 code. The
  attestation is in `results/w2/blind/README.md`. The only project files I read before that were CLAUDE.md,
  DESIGN.md, DATA.md, data/README.md, `src/nextgem/model.py` and `src/nextgem/gurobi_pool.py`.
- **Label exposure: none.** No `labels.*` call, no raw label file, nothing under `data/locked/`. DepMap files read:
  `Model.csv`, `ModelCondition.csv`, `ScreenSequenceMap.csv`, `CRISPRScreenMap.csv` (ScreenID column only; a
  screen map, not a label) for the line-assignment re-derivation. The guard did not block any of them.
- **Gurobi.** Sessions were taken only through `nextgem.gurobi_pool`, one at a time, for:
  - the two full no-shortcut sweeps;
  - the final sterol rescue scan;
  - one aborted demand scan.

  Every other run used GLPK, with `cobra.Configuration().solver` set before any model was loaded.
- **Side effect.** The requested pytest run used `NEXTGEM_NO_TEST_RECORD=1`, so no test-run row was written.
- **Review runtime:** ≈ 1.5 h wall, ≈ 1.4 CPU-h, ≤ 1 Gurobi session at a time.
  - KO verification: 7.8 min.
  - Sterol scans: ≈ 70 min. Two over-broad GLPK and Gurobi scans were stopped and replaced by a targeted 5.5-min
    scan plus two 1-min follow-ups.
  - Diagnostics: ≈ 3 min. Extraction and diff: ≈ 2 min. Tests: 5 s.

## Verdict summary

| # | Item | Verdict |
|---|---|---|
| 1 | Two-agent extraction: compositions and exchange map (script diff) | **PASS**: 312/313 components identical, µM identical; 5 judgement calls adjudicated (§1) |
| 2a | Formulation-string parser | **PASS** with a NOTE: the per-string unit test asserts too little (§2.1) |
| 2b | Line → medium assignment (835 attribution lines) | **PASS**: independent per-line re-derivation, 0 differences over 1,361 lines |
| 2c | Bound scheme M1 / M2(f) / gate | **FIX-REQUIRED**: the cysteine default deviates from the frozen text (F1). ATP uptake opened (F5, non-blocking) |
| 2d | Floor and req_c | **PASS** with NOTE (retinoate is the binding component at 100×; floor masks only biologically empty contrasts) |
| 2e | KO sweep engine and its three shortcuts | **PASS**: 128 random pairs (GLPK, no shortcuts) and 2 full sweeps (Gurobi, no shortcuts) all agree exactly; M0 equals `results/baseline` bitwise |
| 2f | Controls table | **FIX-REQUIRED** (F4): 9 of 28 "passed" rows compare equal values and are non-informative; they must not be scored as passes |
| 2g | Sterol-branch FAIL | **Confirmed**: a GEM-structure issue, not a medium issue (§2.7) |
| 3 | DESIGN ambiguities (README §8) | NOTE: recommendations in §3 |
| 4 | M3 calibration | **Confirmed uncalibrated.** The spec itself produces this (all other uptakes at 1000). Report as uncalibrated and exclude M3 (§4) |
| 5 | `pytest tests/test_w2_media.py` | **PASS**: 232 passed in 3.2 s |

**Final: GO for using the W2 library and the M1/M2(f) sweeps downstream once blocking fixes B1–B3 are done** (§6).
None of them requires re-running the KO sweeps, except B2 if the PI picks the DESIGN-literal cysteine rule for
Human2 DMEM/F12 (2 media × 10 sweeps ≈ 20 sweeps, ≈ 35 CPU-min; the sensitivity sweep already exists, so this is
only a relabelling).

---

## 1. Extraction diff (Part 1 vs W2) — `results/w2/review/diff_extractions.py` → `extraction_diff.tsv`

Levels compared: (1) component lists per medium, (2) concentrations after unit normalisation (µM, rel. tol 1e-6),
plus the vendor-sheet mg/L and MW values, (3) exchange ids and stoichiometry per raw component and parent,
(4) per medium × exchange summed µM. The script exits non-zero if any difference is unadjudicated; it exits 0.

| Level | Result |
|---|---|
| Components | 313 rows on both sides (RPMI 42, DMEM 30, DMEM/F12 50, IMDM 41, HPLM 70, Plasmax 80). 312 match by name. 1 name-only difference (below) |
| Concentrations | **312/312 identical µM**; vendor mg/L and MW identical for every row |
| Exchange ids, both parents | identical for every mapped metabolite except the 5 judgement calls below; both parents use the same ids |
| Per-medium exchange sets | differ only through those 5 components |

| Difference | Blind (B) | W2 (A) | Correct (against source / model) |
|---|---|---|---|
| Plasmax "α-Aminobutyrate" name | `alpha-Aminobutyrate` | `Aminobutyrate` | **B** (the PDF carries a symbol-font α, U+F061 U+F02D, which W2 strips). Name only: same 41 µM, same MAR10206 |
| Acetylglycine (HPLM 90 µM, Plasmax 70 µM) | unmapped | MAR10196 "Acetyl-Glycine" | **A**. The exchange exists in both parents; my name search missed the hyphen. Blind miss |
| NH4Cl / ammonium chloride / NH4VO3 | NH4+ (MAR11420) | NH3 (MAR09073) | **A acceptable.** Both NH4+[e] and NH3[e] have transport reactions in both parents. NH3 is the conventional Human-GEM route (the GR task) and equivalent at availability bounds |
| Thiamine hydrochloride Cl⁻ stoichiometry | 2 Cl⁻ | 1 Cl⁻ (`media.py:165`) | **B**: thiamine chloride hydrochloride, MW 337. **No bound effect** (chloride is open at 1000). Trivial fix F6 |
| Plasmax methyl acetoacetate (41 µM) | acetoacetate (non-literal) | unmapped | **A** as the default (literal; no ester exchange). Acetoacetate is at most a sensitivity |

**Shared caveats (same on both sides, from the sheets):**
- IMDM cystine·2HCl uses MW 240 (free cystine), giving 381 µM. With the salt MW (313) it is 292 µM, which matches
  Gibco's IMDM. IMDM Tyr·2Na uses MW 225 (anhydrous).
- Availability-only bounds make these irrelevant now. Any quantitative use (M3 or later) should use 292 µM.

The script also records the model checks (`level = model_check`): NH4+[e] and NH3[e] connectivity, and Acetyl-Glycine[e] →
MAR09995.
**Unmapped nutrients agree** on both sides: 4-aminobenzoate, nitrate, selenite, Mn²⁺, vanadate, methyl
acetoacetate (+ HEPES and phenol red as non-nutrients). None is a biomass precursor.

**Extraction verdict: PASS.** After the trivial thiamine fix, the W2 tables can be frozen as the reconciled
extraction (B1).

## 2. Code and results review

### 2.1 Formulation-string parser — PASS with NOTE
- `media.parse_formulation` (`media.py:544–573`) and `_parse_supplement` (`media.py:522–538`) parse all 212
  distinct strings and raise on any unknown token (`tests/test_w2_media.py:89–91`).
- I inspected all 77 buildable strings, their medium keys and flags (`formulation_strings.tsv`). Every key is
  correct, for example:
  - `RPMI + 10% FBS + 100uM NEAA` → `RPMI+alanine` (RPMI lacks only Ala);
  - `DMEM:F12 + … NEAA` → `DMEM:F12`;
  - mixtures → linear blends;
  - 2-ME → `+cys(2-ME)`.
- Serum-free classification: the 20 serum-free strings (159 lines) are the organoid, stem and defined media, plain
  `DMEM:F12`, and one F12 + ITS string. No text says "no FBS". Horse serum is counted as serum (flagged).
- **NOTE N1 (test strength, `tests/test_w2_media.py:44–58`).** The per-string test checks only that the parse
  succeeds and that the serum flag is consistent. Only 8 strings have expected values (`:61–86`). DESIGN line 476
  asks for "every distinct string has a unit test". A golden table (the frozen `formulation_strings.tsv` columns
  base / parts / serum_pct / supplements → canonical) compared by the test would meet that cheaply. Non-blocking;
  recommended before the freeze (F3).
- NOTE: "1mM NEAA" is read as 1 mM of each amino acid. This is irrelevant under availability bounds.

### 2.2 Line → medium assignment — PASS
- `media.assign_line_media` (`media.py:718–791`) was re-derived independently in
  `results/w2/review/verify_assignment.py`. That script has its own tier logic over ScreenSequenceMap (QC-passing,
  non-excluded, non-pDNA, in CRISPRScreenMap), ModelCondition and Model.csv.
- **0 per-line differences over 1,361 screened lines.** Tiers are screen 1,225 / onboarded 56 / imputed 74 /
  screen_mixed 6. The attribution set is **835** (RPMI 557, DMEM 180, DMEM/F12 60, IMDM 38), equal to the design
  numbers.
- The modal-text rule counts replicate rows rather than screens (`media.py:739`). This is immaterial today
  (0 models have several texts with one base).

### 2.3 Bound scheme M1 / M2(f) / gate — FIX-REQUIRED (F1), otherwise PASS
`media.uptake_bounds` (`media.py:658–701`) implements the ladder as written:
- M1 availability at 1000, absent = 0;
- O2 and H2O open;
- floor at `floor_mult × req_c` where the medium lacks the component (`:682–684`);
- M2 adds f × req_c for the six synthesisable serum components unless the medium supplies them (`:688–690`), plus
  Fe²⁺/Fe³⁺/Zn²⁺/Cu²⁺ at 1000 (`:691–692`), for serum media only;
- the gate removes Asp/Glu uptake (`:698–700`);
- secretion is never constrained.

The toy tests (`tests/test_w2_media.py:178–225`) cover each rule.

- **F1 (blocking decision): cysteine rule.** `CYSTEINE_RULE = "documented"` (`media.py:72`, applied at
  `:629–630`) opens cysteine wherever a sheet lists free cysteine (DMEM/F12, HPLM, Plasmax). The frozen DESIGN
  text (line 485) says "cysteine closed unless 2-ME is present". This is an undeclared deviation from a frozen
  rule, and the biology favours the DESIGN-literal rule (§3.1).
  - Impact (cached sweeps): Human2 DMEM/F12 changes 83–92 gene ratios, with 0 lethal flips. Human1 is unchanged.
  - Fix: make `design_literal` the default and keep `documented` as the sensitivity. This is a relabelling: the
    sweeps exist as `SENS:DMEM:F12-cysteine-closed`, plus the 2-ME keys. Alternatively, the PI logs an amendment
    keeping "documented".
- **F5 (non-blocking): ATP supplement.** `"ATP": ("MAR00569", …)` (`media.py:109`, parsed at `:507`) opens ATP
  uptake at 1000. This raises WT growth by 43% from ATP import alone (RPMI M1: 46.4 → 66.3, both parents;
  `diag_m3_cys.tsv`), i.e. free energy.
  - Extracellular ATP is hydrolysed by ecto-nucleotidases and is not taken up intact, so it should map to
    adenosine/adenine or be dropped.
  - It affects 1 medium key with 1 line, which is outside the attribution set. Fix before the freeze.
- **NOTE N2: trace components at 1000.** Availability bounds turn µM-level components into bulk carbon sources:
  - GSH 3 µM, hydroxyproline and nucleobases in RPMI;
  - under M3 at O2 = 1, Asp, Gly and hydroxyproline supply ≈ 95% of the carbon in H1 (`diag_m3_cys.tsv`).

  This is inherent to the M1 spec (DESIGN line 485). It shifts WT growth, rarely ratio_r. It is not a fix; it is
  the reason M3 fails (§4).
- The gate: `gate_matters` (`w2_engine.py:46–47`) builds +g only for media with Asp/Glu. Lines with no expression
  keep the gate open, flagged (69 lines). PASS.

### 2.4 Floor and req_c — PASS with NOTE
- **Floor derivation** (`scripts/w2_medium_floor.py:57–86`): drop-one inputs on M0, then for each sheet the missing
  essentials, greedily pruned. The floor is Fe²⁺, α/γ-tocopherol, biotin, cobalamin, linoleate, linolenate,
  lipoate, retinoate in both parents; cholesterol is excluded (`:86`). It matches DESIGN plus the three documented
  additions.
  - Fe²⁺: RPMI and IMDM carry no iron.
  - Biotin and B12: absent from DMEM high-glucose. FBS supplies both in practice, so the floor reproduces what
    serum does.
- **req_c** (`:88–118`) follows DESIGN lines 502–506 exactly:
  - reference medium = the RPMI sheet + floor open at 1000;
  - FVA-min for floor components;
  - pFBA at 99% with c open for serum components.

  All 6 serum req_c values are > 0 in both parents (`results/w2/prechecks/req_nonzero.tsv`; test
  `tests/test_w2_media.py:240–249`).
- **Floor non-binding** (`scripts/w2_prechecks.py:71–87`): Human1 is non-binding at 10×. Human2 binds, is raised
  to 100× per DESIGN, and 20 cases still bind (≤ 1.9%). I traced the binding component (`diag_floor_cad.tsv`).
  - Human2 DMEM/F12 M1 at 100× req = 74.567 vs 74.757 with the floor at 1000.
  - Raising **only retinoate** to 1000 recovers the full difference; the tocopherols and linolenate change nothing.
  - So the model uses retinoate beyond the biomass need (a catabolic sink). This is a 1000-bound artefact, not an
    undersupplied floor. W2's reading is correct.
  - NOTE N3: DESIGN has no rule beyond 100×. Recommend logging "accepted, ≤ 1.9%, retinoate catabolism" rather
    than raising further.

### 2.5 KO sweep engine and the three exact shortcuts — PASS
Engine: `w2_engine.ko_task` (`w2_engine.py:159–205`). Closures are computed once per parent (`:118–126`); they are
medium-independent because exchanges carry no GPR. The three shortcuts:
- (a) no closure → ratio 1;
- (b) every closed reaction carries |v| ≤ 1e-9 in the WT optimum → WT stays feasible → ratio 1 (`:175–177`);
- (c) identical closure sets share one LP (`:178–190`).

All three are exact in exact arithmetic. Verification, `results/w2/review/verify_ko.py` → `ko_verify_summary.tsv`:

| Check | Result |
|---|---|
| M0 sweep vs `results/baseline/{parent}_ko.parquet` ratio_r | **bitwise equal** (`np.array_equal` and byte equality), 2,887 (H1) / 2,870 (H2) genes |
| Sweep definitions rebuilt with W2's `sweep_context` + `build_sweeps` | all 190 + 190 parquet sweep ids reproduced, none extra |
| 128 random gene × sweep pairs, 16 random sweeps, **no shortcut, GLPK** (a different solver), stratified: 48 ratio-1 genes with a non-empty closure (the zero-flux shortcut class), 48 ratio < 1, 16 shared-closure, 16 uniform | **all 128 identical** to W2 (max \|Δ\| = 0) |
| Full no-shortcut sweep, RPMI M2 f = 0.1, Gurobi (1 pool session) | H1 2,887 genes, **0 differences**; H2 2,870 genes, **0 differences** (≈ 160 s each) |

The status policy (infeasible → 0; non-optimal → raise) matches `nextgem.knockouts`. No issues found.

### 2.6 Controls table — FIX-REQUIRED (F4)
`scripts/w2_controls.py` builds `results/w2/controls/positive_controls.tsv`.
- **ASNS** (DMEM 0 vs RPMI 1): PASS in both parents.
- **Mevalonate core** (0/201 labels lose HMGCR/MVK/PMVK/MVD): PASS.
- **Purines rescued in HPLM/Plasmax** (1 vs 0): PASS.
- **M-open collapse**: PASS.

- **F4.** Several rows are scored `passed = True` because a `>=` / `<=` comparison holds with equal values. Those
  rows are non-informative, not passes:
  - GOT1, GOT2, SLC1A3 with the gate (1.0 vs 1.0; `:91–95`);
  - SLC7A11 cysteine vs cystine (1.0 vs 1.0; `:96–100`);
  - "pyruvate relaxes GOT/MDH/OXPHOS" (identical values; `:111–116`);
  - CAD within "pyrimidine core rescued by uridine" (0 vs 0; `:101–110`).

  The README table already calls them non-informative, but the TSV (which downstream code reads) says pass.
  Fix: add a `non_informative` state when the two sides are equal and the expectation is a strict change.
- Why they are non-informative (model facts, for W5):
  - **SLC7A11:** cystine import in Human2 is `SLC7A11 and SLC3A2` OR `SLC7A9 and SLC3A1` (rBAT), and others
    (MAR04931, MAR06891/2, MAR08687). The KO is GPR-buffered whether or not GSH is open (`diag_m3_cys.tsv`).
  - **GOT1/GOT2:** Asp is synthesised either way.
  - **Pyruvate:** RPMI vs RPMI without pyruvate changes **0 genes** in either parent (§3.3).
  - **CAD in Plasmax:** CAD stays lethal although uridine is open at 1000. The KO is rescued by supplying
    dihydroorotate or orotate, but **not** by UMP or uridine (`diag_floor_cad.tsv`). So the model needs
    DHODH/orotate flux for something other than UMP. This is a GEM-structure issue, like ADSL, and goes to W5 /
    the failure account.
- **HPLM uridine.** The a priori list (`notes/apriori_lists_final.md` §1.2, CAD/DHODH/UMPS row) says "HPLM and
  Plasmax contain uridine". Cantor 2017 Table S1 (on disk) has **no uridine**. W2's test (`test_w2_media.py:109`)
  and both extractions agree. The expected-direction check must name Plasmax only.

### 2.7 Sterol-branch FAIL — confirmed: GEM structure, not medium
W2's observation reproduces on Ham's M0 + cholesterol uptake at 1000, which is cleaner than the full medium
(`results/w2/review/sterol_trace.py` → `sterol_trace.tsv`).
- These KOs stay lethal with cholesterol: FDFT1, SQLE, LSS, NSDHL, EBP, DHCR24 (both parents), plus HSD17B7 (H2
  only; ratio 1 in H1 through its GPR).
- DHCR7 is rescued (0 → 1).
- CYP51A1, TM7SF2, MSMO1 and SC5D are ratio 1 throughout (GPR or bypass).

**Missing sterol product: 7-dehydrocholesterol (provitamin D3) → the vitamin D3 catabolites in the biomass
cofactor pool.**
1. **Rescue scan.** For each lethal KO, I sink-supplied each named post-squalene intermediate (c and r
   compartments).
   - **Lathosterol rescues every one** in both parents.
   - Zymosterol and lanosterol rescue only the KOs upstream of them.
   - Nothing downstream of 7-DHC, cholesterol included, rescues.
2. **Flux trace** (`sterol_flux.py` → `sterol_flux.tsv`, FDFT1 KO + cholesterol + lathosterol sink). The supplied
   lathosterol is used through:
   - SC5D lathosterol oxidase → **provitamin D3** (MAR01557);
   - → vitamin D3 by **MAR02114**, the non-enzymatic UV-B photolysis step, which has no GPR;
   - → calcidiol (CYP2R1, MAR02115) → calcitriol (MAR02118) → 24-/23-hydroxylated catabolites.
3. **Where the requirement enters biomass.** The cofactor pool (MAR10065, both parents) demands
   0.0345 mmol/gDW each of **24-oxo-1α,23,25-trihydroxyvitamin D3** and **25-hydroxyvitamin D3-26,23-lactone**.
   The pool also holds about 20 retinoid and tocopherol catabolites, which is also why retinoate and the tocopherols
   are floor components. Cholesterol cannot be converted back to 7-DHC (DHCR7 is irreversible), so:
   - any KO upstream of 7-DHC is lethal whatever the medium;
   - DHCR7 KO is rescued, because 7-DHC is still made and cholesterol comes from serum.
4. **Confirmation** (`sterol_vitd.py` → `sterol_vitd.tsv`).
   - Opening **calcidiol** uptake (MAR09215) with cholesterol rescues FDFT1, SQLE, LSS, NSDHL, HSD17B7, EBP and
     DHCR24 to ratio 1.0 in both parents. HMGCR stays lethal, as expected (prenyl/ubiquinol/dolichol).
   - Opening vitamin D3 uptake (MAR09442) does not rescue: the imported D3 does not reach the CYP2R1 route.

**Classification: GEM-structure issue, not a medium-model issue.** The biomass treats vitamin D3 catabolites as
obligate proliferative cofactors, synthesised de novo through a UV-dependent photolysis step. Proliferating
cultured cells need neither, and cholesterol uptake by LDL rescues loss of the post-squalene pathway in culture:
SQLE-deficient ALCL lines are cholesterol auxotrophs grown on serum/LDL (Garcia-Bermudez et al. 2019, Nature,
[10.1038/s41586-019-0945-5](https://doi.org/10.1038/s41586-019-0945-5); PubMed).

- The medium model is correct. Adding calcidiol to the serum model would mask the artefact, and the frozen serum list
  does not contain it, so it must **not** be added there.
- The correct lever is W3's biomass B-loo/curation of the vitamin D (and retinoid) cofactor-pool terms, plus W5. It
  should be logged in the failure account as "biomass-artefact lethal".
- For W2, the sterol positive control is **non-evaluable by construction** in both parents, so it is not a W2 failure
  (F4/B3).
- The mevalonate-core diagnostic (HMGCR, MVK, PMVK, MVD lethal in all 201 labels) remains valid.

## 3. DESIGN ambiguities (W2 README §8): biology, recommended default and sensitivity

PubMed was the source for the citations below.

| # | Ambiguity | Recommended default | Sensitivity | Basis |
|---|---|---|---|---|
| 3.1 | **Cysteine rule** (DMEM/F12, HPLM and Plasmax list free cysteine) | **DESIGN-literal: cysteine closed unless 2-ME** (reverse W2's default, F1) | documented free cysteine open (sweeps exist) | Cysteine in culture medium is rapidly (Cu-catalysed) autoxidised to cystine, which is why media "usually contain cystine but not cysteine" (Ishii & Bannai 1985, [10.1002/jcp.1041250119](https://doi.org/10.1002/jcp.1041250119); Ishii, Sugita & Bannai 1987, [10.1002/jcp.1041330217](https://doi.org/10.1002/jcp.1041330217)). 2-ME supports growth by reducing cystine extracellularly (same papers). Impact: H2 DMEM/F12 83–92 ratios shift, 0 lethal flips; H1 none |
| 3.2 | **Glutamine default** for the "NoGlutamine" DMEM and DMEM/F12 sheets | **add standard Gln** (DMEM 4 mM, DMEM/F12 2.5 mM; W2's choice, `media.py:263`) | no added Gln (exists) | Standard DMEM (Gibco 11965) and DMEM/F12 (11320) contain L-Gln, and lines are never propagated Gln-free. Only GLUL flips (lethal without Gln), which is biologically true only in Gln-free culture. Needs an amendment-log line, because "documented" here means the standard catalogue formulation, not the sheet on disk |
| 3.3 | **ATCC RPMI sheet** (4.5 g/L glucose, 10 mM HEPES, 1 mM pyruvate) | **classic RPMI-1640, no pyruvate** (Moore 1967 formulation, PMID 4960081; Gibco 11875) | ATCC-modified with pyruvate | DepMap texts add pyruvate explicitly to RPMI ("+ 1mM / 2mM Sodium pyruvate", "+ 25mM Glucose", "+ 25mM HEPES"; 11 strings, 53 lines; pyruvate alone in 5 strings, 12 lines). That implies the base RPMI that labs record is the classic one. **Impact: 0 genes change in either parent** (all M1/M2 levels), so the choice is cosmetic for the sweeps. It matters only for medium keys and M3. A Gibco 11875 sheet should be added to `data/media/` for provenance |
| 3.4 | **Floor masks biotin/B12** (DMEM lacks both) | **keep the floor** | DESIGN's B-loo (W3); no extra variant needed | DMEM lines are cultured with FBS, which supplies biotin and B12. Biologically there is no gene-specific essentiality to recover: without B12 or biotin, every biotin-/B12-dependent route fails at once (an auxotrophic death, not an MMUT- or PCCA-specific dependency), and MMUT/MMAA/MMAB are not proliferation-essential in any medium. The lost "contrasts" were not expected to exist. Label floor items `floor_source = serum` in serum media, `modelling` otherwise |
| 3.5 | **HPLM uridine** | **HPLM has no uridine** (Cantor 2017 Table S1, on disk, [10.1016/j.cell.2017.03.023](https://doi.org/10.1016/j.cell.2017.03.023)); Plasmax 3 µM uridine + 2 µM uracil (Vande Voorde 2019, Sci Adv, doi 10.1126/sciadv.aau7314, Table S1) | none | Correct the a priori expectation (§2.6). HPLM's pyrimidine-relevant feature is uric acid (350 µM), which inhibits UMPS; the model cannot represent that inhibition |
| 3.6 | **Serum % parsed but unused** | **unused** (as DESIGN; f spans 1,000× while serum % spans 4×, 5–20%) | M2(f\*) with f_eff = f\* × pct/10 at the TG-selected f only (no new grid) | FBS total cholesterol is ≈ 34 mg/dL (≈ 0.87 mM; Cantor 2017 Table S1, "Media profiling_non-polar"), so 10% FBS gives ≈ 87 µM in the medium. Lipid supply is roughly proportional to %, but it is already absorbed by f |

**Other README §8 items.**
- 5 (floor at 100×): accept and log (§2.4).
- 6 (M3): §4.
- 8 (no hemin in the serum model): accept. Heme-pathway genes stay lethal under every variant, which is the
  pre-registered failure class.
- 10 (two-agent extraction): done by this review (§1).
- 11–13 (unmappable → imputed medium, mixture handling, horse serum as serum): consistent with DESIGN's strata.
  Accept.
- 14 (script name `w2_medium_floor.py` vs DESIGN's `medium_floor.py`): accept and log.

**NOTE N4 (important for TG and §6.8 binary calls): f = 0.01 collides with the 0.01 lethality threshold.**
- At M2(f = 0.01), ten serum-rescued genes in each parent land at ratio_r 0.00994 in some media and 0.01000 in
  others, so their binary call (< 0.01) flips on 6th-decimal WT-growth differences.
  - H1: ALDH9A1, BBOX1, TMLHE, SHMT2 (carnitine), DHCR7 (cholesterol), SPCS1/2/3, SEC11A/C.
  - H2: the same, with SHMT1 in place of SHMT2.
- No gene is within (0.005, 0.02) at f = 0.001, 0.1 or 1.
- Continuous ratio_r and AP are unaffected. Binary metrics at f = 0.01 are knife-edge and should be reported with
  this caveat. No rule change is proposed.

## 4. M3 calibration — confirmed uncalibrated; it is the spec, not a solver artefact

- **Confirmed** on W2's grid (`results/w2/prechecks/calibration_grid.tsv`):
  - growth / O2 = **12.33 (H1) and 12.30 (H2) exactly**, for every O2 ≤ 0.03, at every glucose (0.05–2) ×
    glutamine (0.01–0.1) × f (0.01, 1);
  - at O2 = 1–10 growth is still 3–17 (model units ≈ h⁻¹);
  - the 0.01–0.06 h⁻¹ window is reached only at O2 ≈ 0.001–0.005 mmol gDW⁻¹ h⁻¹, about 100–1,000× below measured
    cancer-cell O2 uptake.
- **Mechanism** (`results/w2/review/diag_m3_cys.tsv`, RPMI M3 with glucose 0.5 and glutamine 0.1):
  - At low O2, growth is capped by the biosynthetic O2 demand of biomass: fatty-acyl desaturases (MAR00181,
    MAR02293, MAR00145) and the sterol oxygenases (CYP51A1 MAR01497, SQLE, MSMO1).
  - That demand is a fixed ≈ 0.08 mmol O2 per gDW, hence the exact linearity.
  - At O2 = 1, growth (≈ 3 h⁻¹) is fuelled by uptakes the spec leaves at 1000: Asp 47%, Gly 31%, hydroxyproline
    20% of carbon in H1; Glu 37%, pyruvate 30%, Asp 24% in H2. Glucose (≤ 2%) and glutamine (< 1%) are irrelevant.
- **Classification.** This is a direct consequence of the frozen spec: "glucose, glutamine and O2 bounded …, the
  rest as in the best M2" (DESIGN line 487), where the rest are availability bounds of 1000. It is not a numerical
  or coding artefact; the code implements the spec (`media.py:693–697`). It is a modelling artefact *of the spec*:
  with ≥ 30 carbon sources unbounded, no choice of three rates can reach the window at physiological O2.
- **Per the frozen rule** (DESIGN line 508): report M3 as uncalibrated and exclude it from selection. No rescue is
  proposed. For the record, W2 should still write `data/media/calibration.tsv` (or its results-side copy) with the
  literature rates and citations, and evaluate the pre-registered point. The grid shows the outcome, but DESIGN asks
  for the verdict on fixed values (F7, non-blocking).

## 5. Tests
`NEXTGEM_NO_TEST_RECORD=1 .venv/bin/python -m pytest -q tests/test_w2_media.py` → **232 passed in 3.20 s**
(4.8 s wall).

## 6. GO / NO-GO

**GO** for downstream use of the W2 media library (compositions, exchange map, line assignment, attribution set
of 835) and of the M0/M1/M2(f)/gate KO sweeps. The sweep engine is exact (§2.5) and the extraction is reproduced
independently (§1). It is conditional on the blocking fixes below, all of which are cheap.

**Blocking fixes (before the W2 freeze and the first W2 scoring):**
1. **B1 — Freeze the reconciled extraction.**
   - Apply F6 (thiamine HCl → 2 Cl⁻; `media.py:165`) and the Plasmax α-aminobutyrate name.
   - Write the reconciled `compositions.tsv` and `exchange_map.tsv` with their SHA-256 to FROZEN.tsv. The
     DESIGN's `data/media/` location is read-only, so record the results path in an amendment-log row.
   - Attach `results/w2/review/extraction_diff.tsv` as the two-agent diff record (DESIGN line 477).
2. **B2 — Cysteine rule (F1).** Either:
   - switch the default to the DESIGN-literal rule (`CYSTEINE_RULE = "design_literal"`; ladder labels for
     DMEM/F12-based keys point to the cysteine-closed sweeps; documented-cysteine becomes the sensitivity), or
   - log a PI amendment keeping "documented".

   The review recommends the former (§3.1).
3. **B3 — Controls TSV (F4).** Score the equal-value rows as `non_informative`, not `passed`. Mark the sterol-branch
   and CAD-in-Plasmax rows as "model-structure failure (W5)", not medium failures. Correct the HPLM-uridine
   expectation to Plasmax-only.

**Recommended, non-blocking:**
- F2/§3: amendment-log lines for the glutamine default, the RPMI base (adopt classic, no pyruvate; 0 genes change)
  and the 100× floor acceptance.
- F3: a golden-table per-string parser test.
- F5: drop or remap the ATP uptake.
- F7: write `calibration.tsv` and record the M3 verdict on it.
- Fill the FBS reference concentrations in `serum_model.tsv` (documentation only; no bound uses them).
- N4: flag f = 0.01 binary calls as knife-edge.
