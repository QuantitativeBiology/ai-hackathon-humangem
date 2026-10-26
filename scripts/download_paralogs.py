"""Human paralog pairs from Ensembl Compara via BioMart (Ensembl current release; recorded in RELEASE.txt),
plus gene id/symbol/Entrez map. Queried per chromosome-free batch of attributes."""
import time
import requests
from _download_utils import DATA, UA

URL = "https://www.ensembl.org/biomart/martservice"
OUT = DATA / "paralogs"; OUT.mkdir(parents=True, exist_ok=True)

def q(dataset_attrs, filters=""):
    attrs = "".join(f'<Attribute name="{a}"/>' for a in dataset_attrs)
    xml = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" '
           f'uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default">{filters}{attrs}</Dataset></Query>')
    for _ in range(8):
        r = requests.get(URL, params={"query": xml}, headers=UA, timeout=900)
        if r.ok and not r.text.startswith("Query ERROR") and "\t" in r.text.split("\n")[0]:
            return r.text
        print("retry (server busy, sleeping 90 s):", r.status_code, r.text[:80].replace("\n", " "))
        time.sleep(90)
    raise SystemExit("BioMart failed")

(OUT / "ensembl_gene_map.tsv").write_text(q(["ensembl_gene_id", "external_gene_name", "entrezgene_id", "gene_biotype", "chromosome_name"]))
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
PARA = ["ensembl_gene_id", "hsapiens_paralog_ensembl_gene", "hsapiens_paralog_associated_gene_name", "hsapiens_paralog_perc_id",
        "hsapiens_paralog_perc_id_r1", "hsapiens_paralog_subtype", "hsapiens_paralog_orthology_type"]
parts = []
for i, c in enumerate(CHROMS):  # chunk by chromosome: the whole-genome query times out server-side
    cache = OUT / f"_chr{c}.tsv"  # per-chromosome cache so an interrupted run resumes
    if not cache.exists():
        cache.write_text(q(PARA, f'<Filter name="chromosome_name" value="{c}"/><Filter name="with_hsapiens_paralog" excluded="0"/>'))
        time.sleep(10)
    t = cache.read_text()
    lines = t.splitlines()
    parts.extend(lines if i == 0 else lines[1:])
    print(c, len(lines) - 1, flush=True)
(OUT / "ensembl_human_paralogs.tsv").write_text("\n".join(parts) + "\n")
for c in CHROMS:
    (OUT / f"_chr{c}.tsv").unlink()
ver = requests.get(URL, params={"type": "datasets", "mart": "ENSEMBL_MART_ENSEMBL"}, headers=UA, timeout=60).text
line = [l for l in ver.splitlines() if "hsapiens_gene_ensembl" in l][0]
(OUT / "RELEASE.txt").write_text("BioMart hsapiens_gene_ensembl dataset line (contains assembly/release):\n" + line + "\n")
print(line)
