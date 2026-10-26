"""Protein complexes: EMBL-EBI Complex Portal (human, 9606; CC0) and CORUM 5.1 (CC BY 4.0, Zenodo mirror 17419058 of
https://mips.helmholtz-muenchen.de/corum/download, which is unreachable from this host)."""
from _download_utils import DATA, fetch

fetch("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv", DATA / "complexes" / "complexportal_9606.tsv", overwrite=True)
Z = "https://zenodo.org/api/records/17419058/files/{}/content"
for f in ["corum_humanComplexes.txt", "corum_uniprotCorumMapping.txt"]:
    fetch(Z.format(f), DATA / "complexes" / f)
