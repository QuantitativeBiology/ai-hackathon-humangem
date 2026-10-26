# Design review r2: data correctness

DESIGN.md v1 (d84a593), `scripts/design_numbers.py`, `results/design/design_numbers.json`. No label values were read in this round.

## 1. Round-1 must-fixes
| R1 item | Status | Evidence |
|---|---|---|
| 1 Medium tiers | Addressed | §8 W2 has `screen` 1,225 / `screen_mixed` 6 / `onboarded` 56 / `imputed` 74. It states 115 of 118 KY-only lines lack a screen medium and that 57 of 1,220 disagree (my R1 count was 66 of 1,231, because I counted mixed models; the script's definition is cleaner). |
| 2 Library and organoid strata | Addressed | §4.1 step 6 writes `library_class` (1,048/118/194/1), `n_screens`, `screen_type`, `model_type` (108 organoids). Reporting strata at §6 line 207. |
| 3 Evaluation set and fallbacks | Addressed | §4.3: the 69 no-expression lines are scored for generic models and baselines only; the 140 no-CN lines omit the CN constraint (flagged). |
| 4 Background definitions | Mostly addressed | Rule is now "relCN<0.3 or damaging >=1", with 0.1/0.5 sensitivities and a ploidy sensitivity. IDH dropped; groups are listed at M1 from counts. Remaining: see N1. |
| 5 Numbers | Addressed | U = 2,866, U∩ = 2,766, Nusinow 375/309, metabolomics 927/652 are all in the JSON and match my R1 values. |
| 6 Medium-aware baseline | Addressed | B2m added; medium class plus serum flag in B6. |
Licence item (R1 should-fix 11) is in §W8 (NOTICE, separate HPLM/Plasmax release, re-verification before M8).

## 2. `design_numbers.py` check
- Label policy holds. `CRISPRGeneEffect.csv` and `CRISPRGeneDependency.csv` are touched only through `usecols=[0]` and `nrows=0`. The gene-effect matrix is never loaded as values, and the JSON `_meta` says so. No prevalence or performance is computed. Passed.
- The split simulation uses seeds 1..1000, which cannot equal 20261026 + attempt. Passed.
- My independent numbers agree: coverage (1,292 expr, 1,221 CN, 1,214 both, 309 Nusinow, 652 metabolomics, 735 PRISM), 23 eligible lineages with 1,331 lines, rare lines 30, background counts (MTAP 233, STK11 59, ASS1 161, IDH 8, KRAS 209).
- Small points, not blocking:
  - The `screen` tier compares text of `ModelCondition.FormulationID` across screens and calls a model "mixed" if base media differ. A model with several formulations that share a base (e.g. different supplements) takes `texts[0]` arbitrarily; use the modal text or flag `screen_multi_supplement`.
  - `PassesQC == True` is applied at replicate-row level; a screen with one failed replicate still contributes via its passing rows. Fine, but say so.
  - GDSC overlap is null. My count via `Model.COSMICID` (xlsx read by zip parsing) was 600 screened lines, 538 with expression. Keep that as a cross-check once openpyxl is installed.
  - ProCan (947 mapped, 598 screened) joins on `SangerModelID`; confirm the count of unmapped ProCan IDs is printed.
  - The `ne` sets (`AchillesNonessentialControls`) are IDs only; fine.

## 3. Open question R2-9
- **Guide maps: approve, low priority.** The Avana and KY maps (DepMap `AvanaGuideMap.csv`; Behan 2019 / Project Score supplementary guide table) allow a multi-mapping-guide stratum, which is the main non-biological source of false-positive dependencies for CN-amplified or multicopy genes. The CN>=4 and multicopy strata already planned cover most of that. Do not make any claim depend on them. Note that the combined Chronos matrix has already dropped some multi-targeting guides, so the stratum may be thin.
- **openpyxl: approve.** Environment change is trivial and needed for GDSC (`.xlsx`) and for the Gibco/HPLM xlsx sheets now in `data/media/`. Pin it in `requirements.txt` and record it in LESSONS. Without it, GDSC overlap cannot be computed (I used zipfile parsing as a stopgap; not recommended for the pipeline).

## 4. New issues
### Must-fix
**N1. Media sheet coverage is overstated.** §W2 defines "sheet-mappable" as RPMI, DMEM, DMEM:F12, MEM/EMEM, IMDM, F12, McCoy's 5A, L-15, and counts 1,018 (1,034 with mixtures) screen-tier lines. On disk (`data/media/`) there are only RPMI (ATCC), DMEM high glucose, DMEM/F12 and IMDM sheets plus HPLM (mmc1) and the Plasmax PDF. MEM/EMEM (115 lines), F12 (46), McCoy's 5A (21) and L-15 (16) = 198 lines have no composition source yet. State in §W2 which bases are currently usable (RPMI 582, DMEM 185, DMEM:F12 64, IMDM 39 = 870 lines), and make the attribution set and the "documented medium" claim defined on the sheets actually present at M3 (the dn JSON should count by `sheets_on_disk`). Otherwise the stated 1,018 cannot be met and the W2 fallback trigger is ambiguous. If the missing four are not downloaded, the remaining ~198 lines go to `unmappable` and the 13 lineages with two media (n_lineages_with_2_media_ge8) must be recomputed on that basis.

### Should-fix
- N2. The mutation matrix coding (0/1/2) is "verified from documentation at M1". It is DepMap's per-gene count of damaging/hotspot alleles, 0 none, 1 heterozygous, 2 homozygous/two hits. Say this now, since the genomic-constraint rule uses `>=1` (any damaging allele). Loss of one allele (heterozygous) is not a loss of function for an enzyme; consider `>=2` (or `>=1` with LOH from CN) as the primary and `>=1` as sensitivity. Counts at `>=2`: STK11 40 (vs 59 at `>=1` or deletion), KEAP1 18.
- N3. KY-only lines (118) are the lines most often lacking expression (62 of the 69). §4.3 handles the exclusion, but the report should say that the evaluated context-specific set under-represents Sanger-only lines, and the library stratum cannot be evaluated for context models on them.
- N4. `ME2_loss` (27 lines, 29% feasibility), `NAPRT_silenced` (30 lines, 49%), `FH_SDH_loss` (48%) and `PHGDH_amp` (47%) are small; keep them as exploratory and do not count them in any confirmatory family.

## Verdict
Revise minimally, then freeze: only N1 is a must-fix, and it is a text/number correction, not a design change.
