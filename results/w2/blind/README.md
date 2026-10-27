# W2 blind second extraction (media compositions and exchange map)

The second of the two independent extractions that DESIGN §9 W2 ("Compositions") requires. It was produced on
2026-10-27 by an independent review agent (M3 W2 review), label-free.

## Attestation
I attest that before writing every file in this folder (`compositions.tsv`, `exchange_map.tsv`,
`medium_exchanges.tsv`, `unmapped.tsv`, both scripts and this README) I read **none** of these:
- any file under `results/w2/` other than this folder (I created it; I did not list or open its siblings);
- `src/nextgem/media.py`, `src/nextgem/w2_engine.py`, `src/nextgem/w2_runtime.py`;
- `scripts/w2_*.py`;
- `tests/test_w2_media.py`;
- `notes/apriori_lists_*.md`.

I read only these: CLAUDE.md, DESIGN.md (W2 section, §10, §12), DATA.md, data/README.md, `src/nextgem/model.py`,
`src/nextgem/gurobi_pool.py`, and the six source documents in `data/media/`. No label file, no `labels.*` call and
no `data/locked/` access. The label-hook probe was blocked as expected.

## Sources (data/media/, as downloaded 2026-10-26)
| Medium | File | Units in source |
|---|---|---|
| RPMI | `rpmi1640_ATCC.xlsx` (ATCC 30-2001-type modified RPMI: 4.5 g/L glucose, 2.05 mM Gln, 10 mM HEPES, 1 mM pyruvate) | mg/L + MW + mM |
| DMEM | `DMEM_HighGlucose_NoGlutamine_NoPhenolRed.xlsx` (no pyruvate) | mg/L + MW + mM |
| DMEM/F12 | `DMEM_F12_NoGlutamine_NoHEPES.xlsx` | mg/L + MW + mM |
| IMDM | `IMDM_NoPhenolRed.xlsx` (contains 4 mM Gln, 25 mM HEPES, 1 mM pyruvate) | mg/L + MW + mM |
| HPLM | `mmc1.xlsx`, sheet "HPLM formulation" (Cantor et al. 2017 Table S1), rows 7-81; the derived "Small ion" block is excluded | µM |
| Plasmax | `aau7314_sm.pdf` Table S1 pp. 2-3 (Vande Voorde et al. 2019), Plasmax column; "NA" entries excluded | µM (the µ glyph is lost in the PDF text; "α-" is in a symbol font) |

## Files
- `compositions.tsv`: one row per medium × component, in source order (313 rows: RPMI 42, DMEM 30, DMEM/F12 50,
  IMDM 41, HPLM 70, Plasmax 80). `conc_uM` is the normalised concentration: for vendor sheets the sheet's own
  mM × 1000 (every sheet row has mM = mg/L ÷ MW within 1%, `sheet_mM_consistent`), for HPLM/Plasmax the stated µM.
- `exchange_map.tsv`: one row per medium × component × canonical metabolite, with stoichiometry per formula unit
  (`stoich`), `met_uM` = conc_uM × stoich, and the exchange id in Human1 v1.19.0 and Human2 v2.1.1 (exact
  metabolite-name match in compartment `e` among `model.exchanges`; both parents give identical ids for every mapped
  metabolite). `mapping_type`: `literal`, `salt_dissociation`, `nonliteral`, `unmapped`, `unmapped_ion`.
- `medium_exchanges.tsv`: per medium × parent, the set of exchange ids with summed µM over components.
- `unmapped.tsv`: unmapped nutrients per medium.
- `build_blind.py` (the extraction), `dump_exchanges.py` (the exchange dump to `tmp/`, git-ignored, regenerable).

## Mapping decisions (non-literal or judgement calls)
- **Vitamin B12** (cyanocobalamin) → `aquacob(III)alamin` (MAR09269), the only cobalamin exchange. Non-literal.
- **Methyl acetoacetate** (Plasmax) → acetoacetate (MAR09132). Non-literal; the alternative is unmapped.
- **Salts** dissociate with stoichiometry: Ca-pantothenate → 2 pantothenate + Ca²⁺; hemicalcium pantothenate
  (Plasmax) → 1 pantothenate + 0.5 Ca²⁺; thiamine HCl (MW 337, thiamine chloride hydrochloride) → thiamin + 2 Cl⁻;
  cystine·2HCl → cystine + 2 Cl⁻; Tyr·2Na → tyrosine + 2 Na⁺; Na₂HPO₄ → 2 Na⁺ + Pi; Fe(NO₃)₃ → Fe³⁺ + 3 NO₃⁻;
  "Ferric sulfate" (Plasmax) assumed Fe₂(SO₄)₃ → 2 Fe³⁺ + 3 SO₄²⁻ (trace; µM given per formula unit).
- **Isomer assumptions:** lactate → L-lactate; 3-hydroxybutyrate → (R)-3-hydroxybutanoate; α-aminobutyrate →
  (S)-2-aminobutanoate; hydroxyproline → trans-4-hydroxy-L-proline; pyroglutamate → 5-oxoproline; inositol = myo.
- **IMDM pyridoxal** (not pyridoxine) is kept literal: IMDM uses pyridoxal HCl.
- **Sheet anomalies, kept as the sheet states** (flagged in `note`): IMDM cystine·2HCl uses MW 240 (free cystine),
  giving 381 µM; with the salt MW (313) it is 292 µM, which matches the Gibco IMDM sheet (0.29 mM). IMDM Tyr·2Na
  uses MW 225 (anhydrous); with the dihydrate it is 398 µM. IMDM "Asparagine (freebase)" uses MW 150 (monohydrate).
  For availability-only bounds (M1) none of these matters; for any quantitative use, the IMDM cystine value should
  be 292 µM.
- **RPMI sheet** is the ATCC modified formulation (4.5 g/L glucose, pyruvate and HEPES), not the classic Gibco
  11875 (2 g/L glucose, no pyruvate, no HEPES). See the review for the consequence.

## Unmapped nutrients per medium
| Medium | Unmapped components (no exchange in either parent) | Unmapped counter-ions in mapped salts |
|---|---|---|
| RPMI | p-aminobenzoic acid; HEPES (buffer); phenol red (indicator) | nitrate (Ca(NO₃)₂) |
| DMEM | — | nitrate (Fe(NO₃)₃) |
| DMEM/F12 | phenol red | nitrate (Fe(NO₃)₃) |
| IMDM | HEPES | nitrate (KNO₃); selenite (Na₂SeO₃; a selenate exchange exists but is not substituted) |
| HPLM | p-aminobenzoate; N-acetylglycine; phenol red | nitrate (Ca(NO₃)₂) |
| Plasmax | N-acetylglycine; phenol red | vanadate (NH₄VO₃); nitrate (Fe(NO₃)₃); Mn²⁺ (MnCl₂); selenite (Na₂SeO₃) |

None of the unmapped items is a biomass precursor in either parent. HPLM contains **no uridine** in Table S1
(Plasmax has 3 µM uridine and 2 µM uracil).

Water, O₂, CO₂ and H⁺ are not medium components; they are model-level boundary assumptions outside these tables.
