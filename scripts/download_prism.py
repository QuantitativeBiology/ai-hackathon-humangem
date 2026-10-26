"""PRISM Repurposing (Corsello et al. 2020 Nat Cancer), Broad/DepMap, via figshare (CC BY 4.0).
 - Repurposing Public 24Q2 (figshare article 25917643): latest PRISM Repurposing public release found on figshare
   (single-dose primary screen, collapsed log-fold-change; treatment/cell-line metadata). Skips the 0.5-0.7 GB replicate-level matrices.
 - PRISM Repurposing 20Q2 secondary screen (article 20564034): multi-dose dose-response curve parameters (AUC/IC50) for ~1,400 compounds.
Files are resolved from the figshare API (no hardcoded download ids).
"""
import requests
from _download_utils import DATA, UA, fetch

SETS = {
    "repurposing_24Q2": (25917643, ["README.txt", "Repurposing_Public_24Q2_Cell_Line_Meta_Data.csv",
        "Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv", "Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv",
        "Repurposing_Public_24Q2_LFC_COLLAPSED.csv", "Repurposing_Public_24Q2_Treatment_Meta_Data.csv",
        "Repurposing_Public_24Q2_Readme.txt"]),
    "repurposing_20Q2_secondary": (20564034, ["prism-repurposing-20q2-secondary-screen-cell-line-info.csv",
        "prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv",
        "prism-repurposing-20q2-secondary-screen-replicate-collapsed-treatment-info.csv",
        "prism-repurposing-20q2-secondary-screen-readme"]),
}
for sub, (aid, wanted) in SETS.items():
    meta = requests.get(f"https://api.figshare.com/v2/articles/{aid}", headers=UA, timeout=60).json()
    print(sub, meta["title"], meta["published_date"], meta["license"]["name"])
    (DATA / "prism" / sub).mkdir(parents=True, exist_ok=True)
    (DATA / "prism" / sub / "FIGSHARE_ARTICLE.txt").write_text(f"{meta['title']}\n{meta['published_date']}\n{meta['license']['name']}\n{meta['url_public_html']}\n")
    files = {f["name"]: f["download_url"] for f in meta["files"]}
    for w in wanted:
        fetch(files[w], DATA / "prism" / sub / w)
