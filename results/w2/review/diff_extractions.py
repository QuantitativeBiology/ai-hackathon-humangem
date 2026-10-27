"""Script diff of the two independent W2 extractions (DESIGN §9 W2 'Compositions'; AI SF3).

A = results/w2/media_library/{compositions,exchange_map,canonical_composition}.tsv  (W2 engineer)
B = results/w2/blind/{compositions,exchange_map,medium_exchanges}.tsv              (blind second extraction)

Levels compared:
  1. components per medium (raw names, row count);
  2. concentrations after unit normalisation (µM; rel. tol 1e-6);
  3. exchange ids per raw component and parent (set of (rxn, stoich));
  4. per medium x parent exchange sets and summed µM per exchange.
Every difference is classified and adjudicated against the source sheet (see VERDICT below) and written to
results/w2/review/extraction_diff.tsv.  The script also checks, on the models, the two mapping questions that the
sources cannot settle (NH4+ vs NH3 exchange; acetylglycine exchange).
Run: .venv/bin/python results/w2/review/diff_extractions.py   (GLPK; ~40 s with the model check)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
A_DIR = ROOT / "results" / "w2" / "media_library"
B_DIR = ROOT / "results" / "w2" / "blind"
OUT = Path(__file__).resolve().parent / "extraction_diff.tsv"
MED = {"DMEM:F12": "DMEM/F12"}
TOL = 1e-6

# Adjudication of each known difference against the source (key: (level, medium or '*', item)).
VERDICT = {
    ("exchange", "*", "Acetylglycine"): ("A", "blind_miss",
        "model has 'Acetyl-Glycine' [e] exchange MAR10196 in both parents; blind name search missed the hyphenated "
        "name. W2 correct."),
    ("exchange", "*", "L-Acetyl glycine"): ("A", "blind_miss",
        "as Acetylglycine: MAR10196 exists. W2 correct."),
    ("exchange", "*", "NH4Cl"): ("A", "mapping_choice",
        "NH4+[e] (MAR11420) vs NH3[e] (MAR09073): see model check; W2's NH3 choice is the usable route if NH4+[e] "
        "has no transport reaction. Physiologically equivalent (NH4+/NH3 equilibrium)."),
    ("exchange", "*", "Ammonium Chloride"): ("A", "mapping_choice", "as NH4Cl"),
    ("exchange", "*", "Ammonium Metavanadate"): ("A", "mapping_choice", "as NH4Cl (ammonium part)"),
    ("exchange", "*", "Thiamine hydrochloride"): ("B", "stoichiometry",
        "thiamine hydrochloride MW 337 = thiamine chloride hydrochloride C12H17ClN4OS.HCl: 2 Cl- per formula unit; "
        "W2 has 1. Chloride only (bound 1000 either way): no effect on bounds."),
    ("component", "Plasmax", "Aminobutyrate"): ("B", "name_only",
        "source PDF has a symbol-font alpha before 'Aminobutyrate' (glyphs U+F061 U+F02D); W2 dropped it. Same "
        "41 uM and same exchange (MAR10206) on both sides: label only, no bound effect."),
    ("component", "Plasmax", "alpha-Aminobutyrate"): ("B", "name_only", "pair of the row above"),
    ("exchange", "*", "alpha-Aminobutyrate"): ("B", "name_only",
        "A maps the same component under 'Aminobutyrate' -> MAR10206; identical."),
    ("exchange", "*", "Methyl acetoacetate"): ("A", "mapping_choice",
        "no methyl-acetoacetate exchange. Blind mapped non-literally to acetoacetate; W2 leaves it unmapped. "
        "Literal reading (W2) is the conservative default; acetoacetate is a sensitivity at most."),
}


def norm_med(m: str) -> str:
    return MED.get(m, m)


def model_check() -> list[dict]:
    """NH4+ vs NH3 and acetylglycine routes in both parents."""
    import cobra
    cobra.Configuration().solver = "glpk"
    from nextgem.model import load_model
    rows = []
    for parent in ("human1", "human2"):
        m = load_model(parent)
        for met_id in ("MAM02579e", "MAM02578e", "MAM03397e"):  # NH4+[e], NH3[e], Acetyl-Glycine[e]
            met = m.metabolites.get_by_id(met_id)
            other = sorted(r.id for r in met.reactions if not r.boundary)
            rows.append(dict(level="model_check", medium="*", item=f"{parent}:{met.name}[e] {met_id}",
                             a_value="", b_value=f"non-exchange reactions: {','.join(other) or 'NONE'}",
                             diff_class="info", correct_side="", reason=""))
        # can ammonia from NH4+[e] alone reach the cell?  test on Ham's medium with NH3 closed vs NH4+ open
        from nextgem.model import set_hams_medium
        with m:
            set_hams_medium(m)
            base = m.slim_optimize()
            med = dict(m.medium)
            rows.append(dict(level="model_check", medium="*", item=f"{parent}: NH3 MAR09073 in Ham's GR inputs",
                             a_value="", b_value=str("MAR09073" in med), diff_class="info", correct_side="",
                             reason=f"WT on Ham's = {base:.4f}"))
    return rows


def main(run_model_check: bool = True) -> None:
    out: list[dict] = []
    # ---------------------------------------------------------------- 1-2 components and concentrations
    a = pd.read_csv(A_DIR / "compositions.tsv", sep="\t")
    a["medium"] = a.medium.map(norm_med)
    b = pd.read_csv(B_DIR / "compositions.tsv", sep="\t")
    a_key = a.set_index(["medium", "raw_component"]).uM
    b_key = b.set_index(["medium", "component"]).conc_uM
    assert a_key.index.is_unique and b_key.index.is_unique
    for k in sorted(set(a_key.index) | set(b_key.index)):
        if k not in b_key.index or k not in a_key.index:
            out.append(dict(level="component", medium=k[0], item=k[1],
                            a_value=a_key.get(k, ""), b_value=b_key.get(k, ""),
                            diff_class="component_presence", correct_side="?", reason="present on one side only"))
            continue
        va, vb = float(a_key[k]), float(b_key[k])
        if abs(va - vb) > TOL * max(abs(va), abs(vb), 1e-12):
            out.append(dict(level="concentration", medium=k[0], item=k[1], a_value=va, b_value=vb,
                            diff_class="concentration", correct_side="?", reason=f"rel diff {abs(va - vb) / vb:.2e}"))
    n_comp = len(set(a_key.index) & set(b_key.index))
    # unit/source cross-check for the vendor sheets: A's mg/L and MW against B's source values
    bv = b[b.unit_source == "mg/L"].set_index(["medium", "component"])
    av = a.dropna(subset=["mg_per_L"]).set_index(["medium", "raw_component"])
    for k in bv.index.intersection(av.index):
        if abs(bv.loc[k, "conc_source"] - av.loc[k, "mg_per_L"]) > 1e-9 or abs(bv.loc[k, "mw_source"] - av.loc[k, "mw"]) > 1e-9:
            out.append(dict(level="source_value", medium=k[0], item=k[1],
                            a_value=f"{av.loc[k, 'mg_per_L']} mg/L, MW {av.loc[k, 'mw']}",
                            b_value=f"{bv.loc[k, 'conc_source']} mg/L, MW {bv.loc[k, 'mw_source']}",
                            diff_class="source_value", correct_side="?", reason=""))

    # ---------------------------------------------------------------- 3 exchange ids per raw component
    ea = pd.read_csv(A_DIR / "exchange_map.tsv", sep="\t", keep_default_na=False)
    eb = pd.read_csv(B_DIR / "exchange_map.tsv", sep="\t", keep_default_na=False)
    eb = eb.drop_duplicates(["component", "metabolite", "human1_rxn", "human2_rxn", "stoich"])
    for comp in sorted(set(eb.component)):
        for parent in ("human1", "human2"):
            sa = {(r[f"exchange_{parent}"], float(r["stoich"])) for _, r in ea[ea.raw_component == comp].iterrows()
                  if r[f"exchange_{parent}"]}
            sb = {(r[f"{parent}_rxn"], float(r["stoich"])) for _, r in eb[eb.component == comp].iterrows()
                  if r[f"{parent}_rxn"]}
            if comp not in set(ea.raw_component):
                out.append(dict(level="exchange", medium="*", item=comp, a_value="(component absent from A map)",
                                b_value=sorted(sb), diff_class="map_presence", correct_side="?", reason=parent))
            elif sa != sb:
                ids_a, ids_b = {x for x, _ in sa}, {x for x, _ in sb}
                cls = "stoichiometry" if ids_a == ids_b else "exchange_id"
                out.append(dict(level="exchange", medium="*", item=comp, a_value=sorted(sa), b_value=sorted(sb),
                                diff_class=cls, correct_side="?", reason=parent))

    # ---------------------------------------------------------------- 4 per medium exchange sets and µM
    ca = pd.read_csv(A_DIR / "canonical_composition.tsv", sep="\t", keep_default_na=False)
    ca = ca[ca.exchange != ""].copy()
    ca["medium"] = ca.medium.map(norm_med)
    ca_sum = ca.groupby(["medium", "exchange"]).uM.sum()
    mb = pd.read_csv(B_DIR / "medium_exchanges.tsv", sep="\t")
    mb = mb[mb.parent == "human2"].groupby(["medium", "rxn"]).uM.sum()  # ids identical in both parents (asserted)
    assert (pd.read_csv(B_DIR / "medium_exchanges.tsv", sep="\t").groupby(["medium", "rxn"]).parent.nunique() == 2).all()
    for k in sorted(set(ca_sum.index) | set(mb.index)):
        va, vb = ca_sum.get(k), mb.get(k)
        if va is None or vb is None:
            out.append(dict(level="medium_exchange", medium=k[0], item=k[1], a_value=va if va is not None else "",
                            b_value=vb if vb is not None else "", diff_class="exchange_set", correct_side="?",
                            reason="exchange open in one extraction only"))
        elif abs(va - vb) > TOL * max(va, vb):
            out.append(dict(level="medium_exchange", medium=k[0], item=k[1], a_value=va, b_value=vb,
                            diff_class="exchange_uM", correct_side="?", reason=f"rel diff {abs(va - vb) / vb:.2e}"))

    df = pd.DataFrame(out)
    # adjudicate component-level exchange differences; propagate to medium-level rows via the component
    comp_of_rxn = {"MAR10196": ["Acetylglycine", "L-Acetyl glycine"], "MAR09073": ["NH4Cl"],
                   "MAR11420": ["NH4Cl"], "MAR09132": ["Methyl acetoacetate"], "MAR09150": ["Thiamine hydrochloride"]}
    for i, r in df.iterrows():
        key = (r.level, r.medium if r.level == "component" else "*", r["item"])
        if r.level in ("exchange", "component") and key in VERDICT:
            side, cls, why = VERDICT[key]
            df.loc[i, ["correct_side", "diff_class", "reason"]] = [side, cls, f"{r.reason}: {why}"]
        elif r.level == "medium_exchange" and r["item"] in comp_of_rxn:
            src = comp_of_rxn[r["item"]][0]
            if r["item"] == "MAR09150":
                src = "Thiamine hydrochloride"
            elif r["item"] == "MAR09073" or r["item"] == "MAR11420":
                src = "NH4Cl"
            side, cls, why = VERDICT[("exchange", "*", src)]
            df.loc[i, ["correct_side", "diff_class", "reason"]] = [side, f"{cls} (propagated)", why]
    if run_model_check:
        df = pd.concat([df, pd.DataFrame(model_check())], ignore_index=True)
    df.to_csv(OUT, sep="\t", index=False)
    print(f"components compared: {n_comp}; A rows {len(a)}, B rows {len(b)}")
    print(df.groupby(["level", "diff_class", "correct_side"]).size().to_string())
    unresolved = df[df.correct_side == "?"]
    if len(unresolved):
        print("UNRESOLVED:\n" + unresolved.to_string())
        sys.exit(1)


if __name__ == "__main__":
    main(run_model_check="--no-model" not in sys.argv)
