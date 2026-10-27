"""Blind second extraction of W2 media compositions and the component -> exchange map.

Written without reading any W2 output or W2 code (see README.md attestation).
Inputs: data/media/* (sheets as downloaded), results/w2/blind/tmp/exchanges.tsv (dump_exchanges.py).
Outputs: compositions.tsv, exchange_map.tsv, medium_exchanges.tsv, unmapped.tsv (this folder).
"""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pypdf

ROOT = Path(__file__).resolve().parents[3]
MEDIA = ROOT / "data" / "media"
OUT = Path(__file__).resolve().parent

VENDOR = {
    "RPMI": ("rpmi1640_ATCC.xlsx", "ATCC RPMI-1640 (modified: 4.5 g/L glucose, 2 mM Gln, 10 mM HEPES, 1 mM pyruvate)"),
    "DMEM": ("DMEM_HighGlucose_NoGlutamine_NoPhenolRed.xlsx", "DMEM high glucose, no glutamine, no phenol red"),
    "DMEM/F12": ("DMEM_F12_NoGlutamine_NoHEPES.xlsx", "DMEM/F12 (1:1), no glutamine, no HEPES"),
    "IMDM": ("IMDM_NoPhenolRed.xlsx", "IMDM, no phenol red"),
}
SECTIONS = {"Amino Acids", "Vitamins", "Inorganic Salts", "Other Components"}


# ----------------------------------------------------------------------------- compositions
def vendor_rows(medium: str, fname: str) -> list[dict]:
    df = pd.read_excel(MEDIA / fname, header=None)
    assert list(df.iloc[0]) == ["Components", "Molecular Weight", "Concentration (mg/L)", "mM"], fname
    rows, section = [], None
    for i in range(1, len(df)):
        name, mw, mgl, mm = df.iloc[i]
        name = str(name).strip()
        if name in SECTIONS:
            section = name
            continue
        mw, mgl, mm = float(mw), float(mgl), float(mm)
        rows.append(dict(medium=medium, section=section, component=name, conc_source=mgl, unit_source="mg/L",
                         mw_source=mw, mM_source=mm, conc_uM=mm * 1000.0,
                         source=f"data/media/{fname} row {i + 1}"))
    return rows


def hplm_rows() -> list[dict]:
    df = pd.read_excel(MEDIA / "mmc1.xlsx", sheet_name="HPLM formulation", header=None)
    assert str(df.iloc[4, 0]).strip() == "All concentrations in uM"
    rows, section = [], "Glucose"
    for i in range(6, len(df)):
        name, conc = df.iloc[i, 0], df.iloc[i, 1]
        if pd.isna(name):
            if i > 79:
                break  # end of formulation block; rows below are the derived 'Small ion' summary
            continue
        name = str(name).strip()
        if pd.isna(conc):
            section = name
            if name == "Small ion":
                break
            continue
        rows.append(dict(medium="HPLM", section=section, component=name, conc_source=float(conc),
                         unit_source="uM", mw_source=None, mM_source=None, conc_uM=float(conc),
                         source=f"data/media/mmc1.xlsx sheet 'HPLM formulation' row {i + 1}"))
    return rows


def plasmax_rows() -> list[dict]:
    r = pypdf.PdfReader(str(MEDIA / "aau7314_sm.pdf"))
    text = r.pages[1].extract_text() + "\n" + r.pages[2].extract_text()
    assert "Table S1. Comparison between the formulations of Plasmax and HPLM" in text
    rows, section = [], None
    sec_names = ["Proteinogenic Amino Acids", "Non-proteinogenic Amino Acids", "Amino Acids Derivatives",
                 "Other Components", "Inorganic Salts", "Trace Elements", "Vitamins"]
    num = r"(NA|[0-9]*\.?[0-9]+\*?)"
    for raw in text.splitlines():
        line = raw.strip()
        for s in sec_names:  # section header can share a line with the first component
            if line.startswith(s):
                section, line = s, line[len(s):].strip()
        m = re.match(rf"^(.*?)\s+{num}\s+{num}\s*$", line)
        if not m or section is None:
            continue
        name, plasmax, hplm = m.group(1).strip(), m.group(2), m.group(3)
        name = name.replace("\uf061\uf02d", "alpha-")  # PDF symbol-font glyphs for 'α-'
        assert all(ord(ch) < 0xE000 for ch in name), repr(name)
        footnote = ""
        fm = re.match(r"^(.*?\D)\s?(\d)$", name)  # trailing stock-solution footnote digit
        if fm and not name.endswith("B12"):
            name, footnote = fm.group(1).strip(), fm.group(2)
        elif name.startswith("Vitamin B12") and name != "Vitamin B12":
            name, footnote = "Vitamin B12", name[len("Vitamin B12"):].strip()
        if name.startswith("N-Trimethylglycine"):
            name, footnote = "N-Trimethylglycine (betaine)", "2"
        if plasmax == "NA":
            continue
        rows.append(dict(medium="Plasmax", section=section, component=name, conc_source=float(plasmax),
                         unit_source="uM", mw_source=None, mM_source=None, conc_uM=float(plasmax),
                         source=f"data/media/aau7314_sm.pdf Table S1 (p2-3); stock footnote {footnote or '-'}"))
    return rows


# ----------------------------------------------------------------------------- component -> metabolites
# (model metabolite name [e], stoichiometry per formula unit). mapping_type per component.
L, S, N, U = "literal", "salt_dissociation", "nonliteral", "unmapped"
CMAP: dict[str, tuple[str, list[tuple[str, float]], str]] = {
    # amino acids
    "glycine": (L, [("glycine", 1)], ""),
    "l-alanine": (L, [("alanine", 1)], ""), "alanine": (L, [("alanine", 1)], ""),
    "l-arginine": (L, [("arginine", 1)], ""), "arginine": (L, [("arginine", 1)], ""),
    "l-arginine hydrochloride": (S, [("arginine", 1), ("chloride", 1)], ""),
    "l-asparagine": (L, [("asparagine", 1)], ""), "asparagine": (L, [("asparagine", 1)], ""),
    "l-asparagine-h2o": (L, [("asparagine", 1)], "monohydrate"),
    "l-asparagine (freebase)": (L, [("asparagine", 1)], "sheet MW 150 is the monohydrate MW (anhydrous 132)"),
    "l-aspartic acid": (L, [("aspartate", 1)], ""), "aspartate": (L, [("aspartate", 1)], ""),
    "l-cystine 2hcl": (S, [("cystine", 1), ("chloride", 2)], ""),
    "l-cystine": (L, [("cystine", 1)], ""), "cystine": (L, [("cystine", 1)], ""),
    "l-cysteine hydrochloride-h2o": (S, [("cysteine", 1), ("chloride", 1)], ""),
    "l-cysteine": (L, [("cysteine", 1)], ""), "cysteine": (L, [("cysteine", 1)], ""),
    "l-glutamic acid": (L, [("glutamate", 1)], ""), "l-glutamate": (L, [("glutamate", 1)], ""),
    "glutamate": (L, [("glutamate", 1)], ""),
    "l-glutamine": (L, [("glutamine", 1)], ""), "glutamine": (L, [("glutamine", 1)], ""),
    "l-histidine": (L, [("histidine", 1)], ""), "histidine": (L, [("histidine", 1)], ""),
    "l-histidine hydrochloride-h2o": (S, [("histidine", 1), ("chloride", 1)], ""),
    "l-hydroxyproline": (L, [("trans-4-hydroxy-L-proline", 1)], "trans-4-hydroxy-L-proline"),
    "4-hydroxyproline": (L, [("trans-4-hydroxy-L-proline", 1)], "trans-4-hydroxy-L-proline"),
    "4-hydroxy-l-proline": (L, [("trans-4-hydroxy-L-proline", 1)], "trans-4-hydroxy-L-proline"),
    "l-isoleucine": (L, [("isoleucine", 1)], ""), "isoleucine": (L, [("isoleucine", 1)], ""),
    "l-leucine": (L, [("leucine", 1)], ""), "leucine": (L, [("leucine", 1)], ""),
    "l-lysine hydrochloride": (S, [("lysine", 1), ("chloride", 1)], ""),
    "l-lysine": (L, [("lysine", 1)], ""), "lysine": (L, [("lysine", 1)], ""),
    "l-methionine": (L, [("methionine", 1)], ""), "methionine": (L, [("methionine", 1)], ""),
    "l-phenylalanine": (L, [("phenylalanine", 1)], ""), "phenylalanine": (L, [("phenylalanine", 1)], ""),
    "l-proline": (L, [("proline", 1)], ""), "proline": (L, [("proline", 1)], ""),
    "l-serine": (L, [("serine", 1)], ""), "serine": (L, [("serine", 1)], ""),
    "l-threonine": (L, [("threonine", 1)], ""), "threonine": (L, [("threonine", 1)], ""),
    "l-tryptophan": (L, [("tryptophan", 1)], ""), "tryptophan": (L, [("tryptophan", 1)], ""),
    "l-tyrosine disodium salt dihydrate": (S, [("tyrosine", 1), ("Na+", 2)], ""),
    "l-tyrosine disodium salt": (S, [("tyrosine", 1), ("Na+", 2)],
                                 "sheet MW 225 = anhydrous Na2-Tyr; Gibco IMDM lists the dihydrate (MW 261)"),
    "l-tyrosine": (L, [("tyrosine", 1)], ""), "tyrosine": (L, [("tyrosine", 1)], ""),
    "l-valine": (L, [("valine", 1)], ""), "valine": (L, [("valine", 1)], ""),
    # non-proteinogenic amino acids and derivatives
    "alpha-aminobutyrate": (L, [("(S)-2-aminobutanoate", 1)], "L-2-aminobutyrate (Plasmax PDF: symbol-font alpha)"),
    "l-citrulline": (L, [("citrulline", 1)], ""), "citrulline": (L, [("citrulline", 1)], ""),
    "l-homocysteine": (L, [("homocysteine", 1)], ""),
    "l-ornitine": (L, [("ornithine", 1)], "source typo 'Ornitine'"), "ornithine": (L, [("ornithine", 1)], ""),
    "l-pyroglutamate": (L, [("5-oxoproline", 1)], "pyroglutamate = 5-oxoproline"),
    "l-acetyl glycine": (U, [], "N-acetylglycine: no exchange in either parent"),
    "acetylglycine": (U, [], "N-acetylglycine: no exchange in either parent"),
    "l-carnosine": (L, [("carnosine", 1)], ""),
    "glutathione (reduced)": (L, [("GSH", 1)], ""), "glutathione": (L, [("GSH", 1)], "reduced (GSH)"),
    "taurine": (L, [("taurine", 1)], ""),
    "n-trimethylglycine (betaine)": (L, [("betaine", 1)], ""), "betaine": (L, [("betaine", 1)], ""),
    # vitamins
    "biotin": (L, [("biotin", 1)], ""), "d-biotin": (L, [("biotin", 1)], ""),
    "choline chloride": (S, [("choline", 1), ("chloride", 1)], ""), "choline": (L, [("choline", 1)], ""),
    "d-calcium pantothenate": (S, [("pantothenate", 2), ("Ca2+", 1)], "Ca(pantothenate)2, MW 477"),
    "d-pantothenic acid hemicalcium": (S, [("pantothenate", 1), ("Ca2+", 0.5)], "hemicalcium salt: per pantothenate"),
    "pantothenate": (L, [("pantothenate", 1)], "HPLM lists free pantothenate (= 2 x RPMI Ca-pantothenate molarity)"),
    "folic acid": (L, [("folate", 1)], ""), "folate": (L, [("folate", 1)], ""),
    "niacinamide": (L, [("nicotinamide", 1)], ""),
    "para-aminobenzoic acid": (U, [], "PABA: no exchange in either parent (not a human vitamin)"),
    "p-aminobenzoate": (U, [], "PABA: no exchange in either parent (not a human vitamin)"),
    "pyridoxine hydrochloride": (S, [("pyridoxine", 1), ("chloride", 1)], ""),
    "pyridoxine": (L, [("pyridoxine", 1)], ""),
    "pyridoxal hydrochloride": (S, [("pyridoxal", 1), ("chloride", 1)], "IMDM uses pyridoxal, not pyridoxine"),
    "riboflavin": (L, [("riboflavin", 1)], ""),
    "thiamine hydrochloride": (S, [("thiamin", 1), ("chloride", 2)], "thiamine chloride HCl (MW 337): 2 Cl-"),
    "thiamine": (L, [("thiamin", 1)], ""),
    "vitamin b12": (N, [("aquacob(III)alamin", 1)], "cyanocobalamin -> aquacob(III)alamin (only cobalamin exchange)"),
    "vitamin b-12": (N, [("aquacob(III)alamin", 1)], "cyanocobalamin -> aquacob(III)alamin (only cobalamin exchange)"),
    "i-inositol": (L, [("inositol", 1)], "myo-inositol"), "myo-inositol": (L, [("inositol", 1)], ""),
    "ascorbate": (L, [("ascorbate", 1)], ""),
    # salts and trace elements
    "calcium nitrate (ca(no3)2 4h2o)": (S, [("Ca2+", 1), ("nitrate", 2)], "nitrate has no exchange"),
    "ca(no3)2.4h2o": (S, [("Ca2+", 1), ("nitrate", 2)], "nitrate has no exchange"),
    "calcium nitrate": (S, [("Ca2+", 1), ("nitrate", 2)], "nitrate has no exchange"),
    "magnesium sulfate (mgso4) (anhyd.)": (S, [("Mg2+", 1), ("sulfate", 1)], ""),
    "mgso4": (S, [("Mg2+", 1), ("sulfate", 1)], ""), "magnesium sulfate": (S, [("Mg2+", 1), ("sulfate", 1)], ""),
    "potassium chloride (kcl)": (S, [("K+", 1), ("chloride", 1)], ""),
    "kcl": (S, [("K+", 1), ("chloride", 1)], ""), "potassium chloride": (S, [("K+", 1), ("chloride", 1)], ""),
    "sodium bicarbonate (nahco3)": (S, [("Na+", 1), ("HCO3-", 1)], ""),
    "nahco3": (S, [("Na+", 1), ("HCO3-", 1)], ""), "sodium bicarbonate": (S, [("Na+", 1), ("HCO3-", 1)], ""),
    "sodium chloride (nacl)": (S, [("Na+", 1), ("chloride", 1)], ""),
    "nacl": (S, [("Na+", 1), ("chloride", 1)], ""), "sodium chloride": (S, [("Na+", 1), ("chloride", 1)], ""),
    "sodium phosphate dibasic (na2hpo4) anhydrous": (S, [("Na+", 2), ("Pi", 1)], ""),
    "na2hpo4": (S, [("Na+", 2), ("Pi", 1)], ""),
    "sodium phosphate monobasic (nah2po4-h2o)": (S, [("Na+", 1), ("Pi", 1)], ""),
    "sodium phosphate monobasic": (S, [("Na+", 1), ("Pi", 1)], "Plasmax; HPLM column marks dibasic"),
    "calcium chloride (cacl2) (anhyd.)": (S, [("Ca2+", 1), ("chloride", 2)], ""),
    "cacl2": (S, [("Ca2+", 1), ("chloride", 2)], ""), "calcium chloride": (S, [("Ca2+", 1), ("chloride", 2)], ""),
    'ferric nitrate (fe(no3)3"9h2o)': (S, [("Fe3+", 1), ("nitrate", 3)], "nitrate has no exchange"),
    "ferric nitrate": (S, [("Fe3+", 1), ("nitrate", 3)], "nitrate has no exchange"),
    "ferric sulfate": (S, [("Fe3+", 2), ("sulfate", 3)], "Fe2(SO4)3 assumed (source says 'Ferric Sulfate'; "
                                                            "Advanced DMEM/F12 lists ferrous sulfate)"),
    "ferrous sulfate (feso4-7h2o)": (S, [("Fe2+", 1), ("sulfate", 1)], ""),
    "cupric sulfate (cuso4-5h2o)": (S, [("Cu2+", 1), ("sulfate", 1)], ""),
    "cupric sulfate": (S, [("Cu2+", 1), ("sulfate", 1)], ""),
    "magnesium chloride (anhydrous)": (S, [("Mg2+", 1), ("chloride", 2)], ""),
    "mgcl2": (S, [("Mg2+", 1), ("chloride", 2)], ""), "magnesium chloride": (S, [("Mg2+", 1), ("chloride", 2)], ""),
    "zinc sulfate (znso4-7h2o)": (S, [("zinc", 1), ("sulfate", 1)], ""),
    "zinc sulfate": (S, [("zinc", 1), ("sulfate", 1)], ""),
    "potassium nitrate (kno3)": (S, [("K+", 1), ("nitrate", 1)], "nitrate has no exchange"),
    "sodium selenite (na2seo3-5h20)": (S, [("Na+", 2), ("selenite", 1)],
                                       "selenite has no exchange (selenate MAR09169 exists; not substituted)"),
    "sodium selenite": (S, [("Na+", 2), ("selenite", 1)],
                        "selenite has no exchange (selenate MAR09169 exists; not substituted)"),
    "nh4cl": (S, [("NH4+", 1), ("chloride", 1)], ""), "ammonium chloride": (S, [("NH4+", 1), ("chloride", 1)], ""),
    "ammonium metavanadate": (S, [("NH4+", 1), ("vanadate", 1)], "vanadate has no exchange"),
    "manganous chloride": (S, [("Mn2+", 1), ("chloride", 2)], "Mn2+ has no exchange"),
    # other components
    "d-glucose (dextrose)": (L, [("glucose", 1)], ""), "glucose": (L, [("glucose", 1)], ""),
    "d-glucose": (L, [("glucose", 1)], ""),
    "hepes": (U, [], "buffer; not a model metabolite"),
    "phenol red": (U, [], "pH indicator; not a model metabolite"),
    "sodium pyruvate": (S, [("pyruvate", 1), ("Na+", 1)], ""), "pyruvate": (L, [("pyruvate", 1)], ""),
    "hypoxanthine na": (S, [("hypoxanthine", 1), ("Na+", 1)], ""), "hypoxanthine": (L, [("hypoxanthine", 1)], ""),
    "linoleic acid": (L, [("linoleate", 1)], ""),
    "lipoic acid": (L, [("lipoic acid", 1)], ""),
    "putrescine 2hcl": (S, [("putrescine", 1), ("chloride", 2)], ""),
    "thymidine": (L, [("thymidine", 1)], ""),
    "2-hydroxybutyrate": (L, [("2-hydroxybutyrate", 1)], ""),
    "3-hydroxybutyrate": (L, [("(R)-3-hydroxybutanoate", 1)], "D-(R)-beta-hydroxybutyrate assumed"),
    "3-hydroxyisobutyrate": (L, [("3-hydroxyisobutyrate", 1)], ""),
    "acetate": (L, [("acetate", 1)], ""), "acetone": (L, [("acetone", 1)], ""),
    "acetyl carnitine": (L, [("O-acetylcarnitine", 1)], ""),
    "carnitine": (L, [("L-carnitine", 1)], ""),
    "citrate": (L, [("citrate", 1)], ""), "creatine": (L, [("creatine", 1)], ""),
    "creatinine": (L, [("creatinine", 1)], ""), "formate": (L, [("formate", 1)], ""),
    "fructose": (L, [("fructose", 1)], ""), "galactose": (L, [("galactose", 1)], ""),
    "glycerol": (L, [("glycerol", 1)], ""),
    "lactate": (L, [("L-lactate", 1)], "L-lactate assumed"),
    "malonate": (L, [("malonate", 1)], ""),
    "methyl acetoacetate": (N, [("acetoacetate", 1)], "methyl ester: no exchange; acetoacetate after hydrolysis "
                                                          "(non-literal; alternative = unmapped)"),
    "succinate": (L, [("succinate", 1)], ""),
    "uracil": (L, [("uracil", 1)], ""),
    "urate": (L, [("urate", 1)], ""), "uric acid": (L, [("urate", 1)], ""),
    "urea": (L, [("urea", 1)], ""), "uridine": (L, [("uridine", 1)], ""),
}

# ions that are counter-ions only: absence of their exchange is not a nutrient gap
COUNTER_ONLY_UNMAPPED = {"nitrate", "selenite", "vanadate", "Mn2+"}


def main() -> None:
    comp = []
    for med, (fname, _desc) in VENDOR.items():
        comp += vendor_rows(med, fname)
    comp += hplm_rows()
    comp += plasmax_rows()
    comp = pd.DataFrame(comp)
    comp.insert(1, "order", comp.groupby("medium").cumcount() + 1)
    # internal consistency check of vendor sheets: mM vs mg/L / MW
    v = comp.unit_source == "mg/L"
    rel = (comp.loc[v, "conc_source"] / comp.loc[v, "mw_source"] - comp.loc[v, "mM_source"]).abs() / comp.loc[v, "mM_source"]
    comp["sheet_mM_consistent"] = ""
    comp.loc[v, "sheet_mM_consistent"] = (rel < 0.01).map({True: "yes", False: "NO"})
    comp["note"] = ""

    def note(medium, component, text):
        m = (comp.medium == medium) & (comp.component == component)
        assert m.sum() == 1, (medium, component)
        comp.loc[m, "note"] = text

    note("IMDM", "L-Cystine 2HCl", "sheet MW 240 is free cystine; with the 2HCl MW (313) 91.4 mg/L = 292 uM "
         "(Gibco IMDM: 91.24 mg/L, 0.29 mM). conc_uM keeps the sheet value 381 uM")
    note("IMDM", "L-Tyrosine disodium salt", "sheet MW 225 (anhydrous); with dihydrate MW 261, 398 uM")
    note("IMDM", "L-Asparagine (freebase)", "MW 150 is the monohydrate")
    note("RPMI", "D-Glucose (Dextrose)", "ATCC modified RPMI: 4500 mg/L (standard Gibco 11875: 2000 mg/L)")
    note("HPLM", "Pantothenate", "free pantothenate (RPMI vitamin mix Ca-pantothenate x2)")
    note("HPLM", "Na2HPO4", "Plasmax table footnote: HPLM phosphate present as dibasic salt")
    note("Plasmax", "D-Pantothenic acid hemicalcium", "hemicalcium salt; concentration per pantothenate")
    comp.to_csv(OUT / "compositions.tsv", sep="\t", index=False)

    # exchange map
    ex = pd.read_csv(OUT / "tmp" / "exchanges.tsv", sep="\t")
    rows, missing_keys = [], []
    for r in comp.itertuples():
        key = r.component.strip().lower()
        if key not in CMAP:
            missing_keys.append((r.medium, r.component))
            continue
        mtype, mets, cnote = CMAP[key]
        if not mets:
            rows.append(dict(medium=r.medium, component=r.component, conc_uM=r.conc_uM, mapping_type=mtype,
                             metabolite="", stoich=None, met_uM=None, human1_rxn="", human2_rxn="", note=cnote))
            continue
        for met, st in mets:
            ids = {}
            for parent in ("human1", "human2"):
                hit = ex[(ex.parent == parent) & (ex.name == met) & (ex.comp == "e")].rxn.unique()
                assert len(hit) <= 1, (parent, met, hit)
                ids[parent] = hit[0] if len(hit) else ""
            mt = mtype
            if not ids["human1"] and not ids["human2"]:
                mt = "unmapped_ion" if met in COUNTER_ONLY_UNMAPPED else "unmapped"
            rows.append(dict(medium=r.medium, component=r.component, conc_uM=r.conc_uM, mapping_type=mt,
                             metabolite=met, stoich=st, met_uM=r.conc_uM * st,
                             human1_rxn=ids["human1"], human2_rxn=ids["human2"], note=cnote))
    assert not missing_keys, missing_keys
    em = pd.DataFrame(rows)
    em.to_csv(OUT / "exchange_map.tsv", sep="\t", index=False)

    # per medium x parent exchange set with summed concentration
    me = []
    for parent in ("human1", "human2"):
        col = f"{parent}_rxn"
        x = em[em[col] != ""].groupby(["medium", col, "metabolite"], as_index=False).agg(
            uM=("met_uM", "sum"), components=("component", lambda s: "; ".join(sorted(set(s)))))
        x = x.rename(columns={col: "rxn"})
        x.insert(1, "parent", parent)
        me.append(x)
    pd.concat(me).sort_values(["medium", "parent", "rxn"]).to_csv(OUT / "medium_exchanges.tsv", sep="\t", index=False)

    # unmapped nutrients per medium (component-level, plus unmappable ions within mapped salts)
    un = em[em.mapping_type.isin(["unmapped", "unmapped_ion"])][["medium", "component", "metabolite", "conc_uM",
                                                                 "mapping_type", "note"]]
    un.to_csv(OUT / "unmapped.tsv", sep="\t", index=False)
    print(comp.groupby("medium").size())
    print(em.mapping_type.value_counts())
    print(un.to_string())
    print(comp[comp.sheet_mM_consistent == "NO"][["medium", "component", "conc_source", "mw_source", "mM_source"]])


if __name__ == "__main__":
    main()
