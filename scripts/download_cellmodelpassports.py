"""Sanger Cell Model Passports model list (latest; model annotations incl. Sanger/COSMIC/DepMap IDs). Public download area."""
from _download_utils import DATA, fetch

fetch("https://cog.sanger.ac.uk/cmp/download/model_list_latest.csv.gz", DATA / "cellmodelpassports" / "model_list_latest.csv.gz", overwrite=True)
