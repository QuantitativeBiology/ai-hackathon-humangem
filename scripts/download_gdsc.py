"""GDSC1 & GDSC2 drug sensitivity, Genomics of Drug Sensitivity in Cancer release 8.5 (Oct 2023), Sanger.
Source: https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/  (https://www.cancerrxgene.org/downloads)
"""
from _download_utils import DATA, fetch

BASE = "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/"
for f in ["GDSC1_fitted_dose_response_27Oct23.xlsx", "GDSC2_fitted_dose_response_27Oct23.xlsx",
          "screened_compounds_rel_8.5.csv", "Cell_Lines_Details.xlsx"]:
    fetch(BASE + f, DATA / "gdsc" / f)
