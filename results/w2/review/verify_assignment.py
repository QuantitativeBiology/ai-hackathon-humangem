"""Independent per-line re-derivation of the W2 screen-tier medium assignment (review M3-W2; label-free).

Own implementation of DESIGN W2 "Medium source tiers" from ScreenSequenceMap / ModelCondition / Model.csv metadata
(no label file), compared line by line with results/w2/media_library/line_media.tsv.
Writes results/w2/review/assignment_diff.tsv (differences only) and prints a summary.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DM = ROOT / "data" / "depmap"
ON_DISK = {"RPMI", "DMEM", "DMEM:F12", "IMDM"}


def main() -> None:
    lm = pd.read_csv(ROOT / "results" / "w2" / "media_library" / "line_media.tsv", sep="\t")
    ssm = pd.read_csv(DM / "ScreenSequenceMap.csv", usecols=["ScreenID", "ModelConditionID", "ModelID", "ScreenType",
                                                              "ExcludeFromCRISPRCombined", "PassesQC"])
    try:
        screens = set(pd.read_csv(DM / "CRISPRScreenMap.csv", usecols=["ScreenID"]).ScreenID)
    except Exception as e:  # noqa: BLE001 - map file may be guarded; fall back to all screens
        print("CRISPRScreenMap not read:", e)
        screens = None
    mc = pd.read_csv(DM / "ModelCondition.csv", usecols=["ModelConditionID", "FormulationID"])
    mo = pd.read_csv(DM / "Model.csv", usecols=["ModelID", "FormulationID"]).set_index("ModelID").FormulationID
    s = ssm[(ssm.PassesQC.astype(str) == "True") & (ssm.ExcludeFromCRISPRCombined.astype(str) == "False")
            & (ssm.ScreenType != "pDNA")]
    if screens is not None:
        s = s[s.ScreenID.isin(screens)]
    s = s.merge(mc, on="ModelConditionID", how="left")
    rows = []
    for mid in lm.ModelID:
        texts = s[s.ModelID == mid].FormulationID.dropna()
        if len(texts):
            bases = {t.split(" + ")[0].strip() for t in texts}
            if len(bases) > 1:
                tier, text = "screen_mixed", None
            else:
                tier, text = "screen", texts.value_counts().index[0]
        elif pd.notna(mo.get(mid)):
            tier, text = "onboarded", mo[mid]
        else:
            tier, text = "imputed", None
        base = text.split(" + ")[0].strip() if text else None
        attrib = tier == "screen" and base in ON_DISK
        rows.append(dict(ModelID=mid, tier_rev=tier, text_rev=text, base_rev=base, attribution_rev=attrib))
    rv = pd.DataFrame(rows).merge(lm[["ModelID", "tier", "text", "attribution_set", "base_class"]], on="ModelID")
    diff = rv[(rv.tier_rev != rv.tier) | (rv.text_rev.fillna("") != rv.text.fillna(""))
              | (rv.attribution_rev != rv.attribution_set)]
    diff.to_csv(Path(__file__).resolve().parent / "assignment_diff.tsv", sep="\t", index=False)
    att = rv[rv.attribution_rev]
    print("lines", len(rv), "| tiers", rv.tier_rev.value_counts().to_dict())
    print("attribution (review)", len(att), att.base_rev.value_counts().to_dict(), "| W2", int(lm.attribution_set.sum()))
    print("per-line differences:", len(diff))
    print(diff.head(20).to_string())


if __name__ == "__main__":
    main()
