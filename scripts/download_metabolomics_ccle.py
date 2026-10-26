"""CCLE metabolomics (Li et al. 2019 Nat Med; 225 metabolites x 928 lines), Broad CCLE legacy download area.
Saves to data/metabolomics_ccle/.  Cell line annotations give the CCLE_ID<->lineage mapping; the matrix also carries DepMap_ID.
"""
from _download_utils import DATA, fetch

BASE = "https://data.broadinstitute.org/ccle/"
for f in ["CCLE_metabolomics_20190502.csv", "Cell_lines_annotations_20181226.txt"]:
    fetch(BASE + f, DATA / "metabolomics_ccle" / f)
