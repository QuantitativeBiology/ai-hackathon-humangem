"""Nusinow et al. 2020 Cell (Quantitative proteomics of the Cancer Cell Line Encyclopedia), TMT proteomics of 375 lines.
Source: Gygi lab, https://gygi.hms.harvard.edu/publications/ccle.html (raw: MassIVE MSV000085836).
Saves to data/proteomics_nusinow/.  Usage: .venv/bin/python scripts/download_proteomics_nusinow.py
"""
from _download_utils import DATA, fetch

BASE = "https://gygi.hms.harvard.edu/data/ccle/"
FILES = [
    "Table_S1_Sample_Information.xlsx",
    "protein_quant_current_normalized.csv.gz",  # normalized protein quant (equivalent to Table S2)
    "Table_S4_Protein_RNA_Correlation_and_Enrichments.xlsx",
]
for f in FILES:
    fetch(BASE + f, DATA / "proteomics_nusinow" / f)
