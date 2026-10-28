"""W6 inputs (label-free): compound table, screen matching, responses, assay media, omics and comparator blocks,
subtype task list (generated before any fitting), coverage.

    .venv/bin/python scripts/w6_inputs.py [--parent human1] [--allow-dirty]

Committed outputs (results/w6/inputs/): compounds.tsv, positive_control_pairs.tsv, compound_matches.tsv,
assay_media.tsv, subtype_tasks.tsv, coverage.json. Regenerable caches (results/w6/tmp/inputs/): lines, RNA,
proteomics, reaction-expression blocks, responses, metabolomics, subtype labels (per-line tables stay in tmp).

No CRISPR label file, no label-derived list and no data/locked/ path is opened.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per process (<= 6 worker processes in total; DESIGN §11.4)

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w6_common as C  # noqa: E402

from nextgem import paths, w6  # noqa: E402

DATA = paths.data_root()
IN_DIR = C.OUT_DIR / "inputs"
TMP = C.tmp_dir("inputs")
PRISM_SEC = DATA / "prism/repurposing_20Q2_secondary/prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv"
PRISM_PRI = DATA / "prism/repurposing_24Q2/Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv"
PRISM_PRI_LIST = DATA / "prism/repurposing_24Q2/Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv"
GDSC2 = DATA / "gdsc/GDSC2_fitted_dose_response_27Oct23.xlsx"
GDSC_CMP = DATA / "gdsc/screened_compounds_rel_8.5.csv"
GDSC_LINES = DATA / "gdsc/Cell_Lines_Details.xlsx"
CTRP_AUC = DATA / "ctrp/CTRPAUCMatrix.csv"
CTRP_COND = DATA / "ctrp/CTRPLog2ViabilityCollapsedConditions.csv"
METAB = DATA / "metabolomics_ccle/CCLE_metabolomics_20190502.csv"
MODEL_CSV = DATA / "depmap/Model.csv"
MOLSUB = DATA / "depmap/OmicsInferredMolecularSubtypes.csv"
EXPR = DATA / "depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
PROCAN = DATA / "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"
LINE_MEDIA = C.PROJECT / "results/w2/media_library/line_media.tsv"
SPLIT = C.PROJECT / "results/splits/lineage_split.tsv"
PROT_MIN_DETECTED = 0.8  # proteins detected in >= 80% of cohort ProCan lines; the rest imputed per line (own floor)


def gene_tsv(parent: str) -> Path:
    return DATA / ("Human1/model/genes.tsv" if parent == "human1" else "Human-GEM-main/model/genes.tsv")


def read_lines() -> pd.DataFrame:
    lm = pd.read_csv(LINE_MEDIA, sep="\t", usecols=["ModelID", "medium_key", "medium_source", "serum_free"])
    sp = pd.read_csv(SPLIT, sep="\t", usecols=["ModelID", "lineage", "tier"])
    df = lm.merge(sp, on="ModelID", how="inner")
    assert len(df) == len(lm) == len(sp), "line_media and split must cover the same cohort"
    gd = pd.read_excel(GDSC_LINES, sheet_name=0)
    gd = gd[["COSMIC identifier", "Screen Medium"]].dropna()
    mdl = pd.read_csv(MODEL_CSV, usecols=["ModelID", "COSMICID", "SangerModelID"])
    mdl = mdl.dropna(subset=["COSMICID"])
    mdl["COSMICID"] = mdl.COSMICID.astype(int)
    gd["COSMIC identifier"] = gd["COSMIC identifier"].astype(int)
    g = mdl.merge(gd, left_on="COSMICID", right_on="COSMIC identifier")[["ModelID", "Screen Medium"]]
    g = g.drop_duplicates("ModelID")
    df = df.merge(g.rename(columns={"Screen Medium": "gdsc_screen_medium"}), on="ModelID", how="left")
    return df.set_index("ModelID").sort_index()


def read_rna(lines) -> pd.DataFrame:
    import pyarrow.csv as pacsv

    tab = pacsv.read_csv(EXPR)
    df = tab.to_pandas()
    df = df[df.IsDefaultEntryForModel == "Yes"]
    meta = [c for c in df.columns if not c.endswith(")")]
    df = df.set_index("ModelID")
    assert df.index.is_unique
    df = df.drop(columns=[c for c in meta if c != "ModelID" and c in df.columns])
    df = df.reindex([l for l in lines if l in df.index])
    return df.astype("float32")


def read_procan(lines) -> pd.DataFrame:
    head = pd.read_csv(PROCAN, sep="\t", header=None, nrows=1, dtype=str)
    uniprot = head.iloc[0, 2:].tolist()
    body = pd.read_csv(PROCAN, sep="\t", header=None, skiprows=3, dtype={0: str, 1: str}).set_index(1).drop(columns=0)
    body.columns = uniprot
    mdl = pd.read_csv(MODEL_CSV, usecols=["ModelID", "SangerModelID"]).dropna()
    mdl = mdl[mdl.ModelID.isin(set(lines)) & mdl.SangerModelID.isin(set(body.index))]
    X = body.loc[mdl.SangerModelID].astype("float64")
    X.index = pd.Index(mdl.ModelID.to_numpy(), name="ModelID")
    X = X[~X.index.duplicated()].sort_index()
    det = X.notna().mean()
    X = X.loc[:, det >= PROT_MIN_DETECTED]
    X = X.T.groupby(level=0).mean().T  # duplicated accessions
    return w6.impute_line_floor(X).astype("float32")


def model_blocks(parent: str, rna: pd.DataFrame, prot: pd.DataFrame):
    import cobra

    cobra.Configuration().solver = "glpk"
    from nextgem.context import compile_gprs
    from nextgem.genes import parse_depmap_columns
    from nextgem.model import load_model

    m = load_model(parent)
    gprs = compile_gprs(m)
    gt = pd.read_csv(gene_tsv(parent), sep="\t", dtype=str).set_index("genes")
    model_genes = sorted(g.id for g in m.genes)
    entrez = gt.geneEntrezID.reindex(model_genes).str.split(";").str[0]
    dm = parse_depmap_columns(rna.columns)
    col_of = dict(zip(dm.entrez, dm.column))
    lin = {g: 2.0 ** rna[col_of[e]].astype("float64") - 1.0 for g, e in entrez.items() if isinstance(e, str) and e in col_of}
    rna_lin = pd.DataFrame(lin, index=rna.index)
    rxn_rna = w6.log_capacity(w6.reaction_capacity(gprs, rna_lin))
    up = gt.geneUniProtID.reindex(model_genes)
    plin = {}
    for g, u in up.items():
        if isinstance(u, str):
            us = [x for x in u.split(";") if x in prot.columns]
            if us:
                plin[g] = (2.0 ** prot[us].astype("float64")).max(axis=1)
    prot_lin = pd.DataFrame(plin, index=prot.index)
    rxn_prot = w6.log_capacity(w6.reaction_capacity(gprs, prot_lin))
    return rxn_rna.astype("float32"), rxn_prot.astype("float32"), rna_lin.shape[1], prot_lin.shape[1]


# ----------------------------------------------------------------------------------------------- responses
def prism_tables(lines):
    sec = pd.read_csv(PRISM_SEC, usecols=["depmap_id", "name", "auc", "screen_id"])
    sec = sec[sec.depmap_id.isin(set(lines))].dropna(subset=["auc", "name"])
    sec_long = sec.groupby(["name", "depmap_id"], as_index=False).auc.median()
    sec_idx = w6.name_index((n, n) for n in sec_long.name.unique())
    sec_n = sec_long.groupby("name").depmap_id.nunique().to_dict()
    lst = pd.read_csv(PRISM_PRI_LIST).drop_duplicates("IDs")
    pri = pd.read_csv(PRISM_PRI, index_col=0)
    pri = pri[[c for c in pri.columns if c in set(lines)]]
    pri_idx = w6.name_index(list(zip(lst.IDs, lst["Drug.Name"])) + list(zip(lst.IDs, lst.Synonyms)))
    pri_n = pri.notna().sum(axis=1).to_dict()
    return sec_long, sec_idx, sec_n, pri, pri_idx, pri_n


def gdsc2_table(lines):
    g = pd.read_excel(GDSC2, usecols=["DATASET", "SANGER_MODEL_ID", "DRUG_ID", "DRUG_NAME", "AUC"])
    g = g[g.DATASET == "GDSC2"]
    mdl = pd.read_csv(MODEL_CSV, usecols=["ModelID", "SangerModelID"]).dropna()
    g = g.merge(mdl, left_on="SANGER_MODEL_ID", right_on="SangerModelID")
    g = g[g.ModelID.isin(set(lines))]
    g = g.groupby(["DRUG_ID", "ModelID"], as_index=False).agg(AUC=("AUC", "median"), DRUG_NAME=("DRUG_NAME", "first"))
    syn = pd.read_csv(GDSC_CMP, usecols=["DRUG_ID", "DRUG_NAME", "SYNONYMS"])
    ids = set(g.DRUG_ID)
    pairs = [(i, n) for i, n in zip(g.DRUG_ID, g.DRUG_NAME)] + [(i, s) for i, s in zip(syn.DRUG_ID, syn.SYNONYMS) if i in ids]
    return g, w6.name_index(pairs), g.groupby("DRUG_ID").ModelID.nunique().to_dict()


def ctrp_table(lines):
    auc = pd.read_csv(CTRP_AUC, index_col=0)
    auc = auc.loc[[l for l in auc.index if l in set(lines)]]
    cond = pd.read_csv(CTRP_COND, usecols=["CompoundName", "CompoundID"]).drop_duplicates()
    cond = cond[cond.CompoundID.isin(auc.columns)]
    return auc, w6.name_index(zip(cond.CompoundID, cond.CompoundName)), auc.notna().sum().to_dict()


def main(args) -> int:
    text = (C.PROJECT / w6.APRIORI_FILE).read_text()
    table = w6.compound_table(text)
    C.atomic_write_tsv(table, IN_DIR / "compounds.tsv", header_lines=[w6.PROVISIONAL_TAG + " run inputs; list = "
                                                                       f"{w6.APRIORI_FILE} §4 (A2), verbatim"])
    pcp = w6.positive_control_pairs_long()
    C.atomic_write_tsv(pcp, IN_DIR / "positive_control_pairs.tsv")

    lines_df = read_lines()
    lines = list(lines_df.index)
    print(f"cohort lines {len(lines)}", flush=True)
    rna = read_rna(lines)
    print(f"RNA {rna.shape}", flush=True)
    prot = read_procan(lines)
    print(f"ProCan {prot.shape}", flush=True)
    rxn_rna, rxn_prot, n_g_rna, n_g_prot = model_blocks(args.parent, rna, prot)
    print(f"rxnRNA {rxn_rna.shape}, rxnProt {rxn_prot.shape}", flush=True)
    feat_lines = set(rna.index)  # "lines with features": cohort lines with an expression profile (and a medium)

    culture = lines_df.medium_key
    gdsc_med = lines_df.gdsc_screen_medium
    maps = {k: w6.assay_medium_map(k, lines, culture, gdsc_med) for k in w6.ASSAY_MAPS}
    am = pd.concat({k: v.value_counts() for k, v in maps.items()}).rename("n_lines").reset_index()
    am.columns = ["assay_map", "medium_key", "n_lines"]
    am["source"] = am.assay_map.map({"prism": w6.PRISM_MEDIUM_SOURCE,
                                     "gdsc": "GDSC Cell_Lines_Details.xlsx 'Screen Medium' (R -> RPMI, D/F12 -> DMEM:F12)",
                                     "culture": "frozen W2 library line_media.tsv medium_key (A9)"})
    C.atomic_write_tsv(am, IN_DIR / "assay_media.tsv")

    # ---- responses
    sec_long, sec_idx, sec_n, pri, pri_idx, pri_n = prism_tables(lines)
    g2, g2_idx, g2_n = gdsc2_table(lines)
    ct, ct_idx, ct_n = ctrp_table(lines)
    toks = list(table.compound)
    m_sec = w6.match_compounds(toks, sec_idx, sec_n).assign(dataset="PRISM", source="PRISM 20Q2 secondary AUC")
    m_pri = w6.match_compounds(toks, pri_idx, pri_n).assign(dataset="PRISM", source="PRISM 24Q2 primary LFC")
    m_g2 = w6.match_compounds(toks, g2_idx, g2_n).assign(dataset="GDSC2", source="GDSC2 AUC")
    m_ct = w6.match_compounds(toks, ct_idx, ct_n).assign(dataset="CTRP", source="CTRP AUC")
    resp, match_rows = [], []

    def add(mrow, values: pd.Series):
        v = values.dropna()
        v = v[v.index.isin(lines)]
        dmap = maps[w6.DATASET_MAP[mrow["dataset"]]]
        n_feat = int(sum(1 for l in v.index if l in feat_lines and l in dmap.index))
        r = dict(mrow)
        r.update(n_lines_response=len(v), n_lines_features=n_feat, eligible=n_feat >= w6.MIN_LINES)
        match_rows.append(r)
        if n_feat >= w6.MIN_LINES:
            resp.append(pd.DataFrame({"dataset": mrow["dataset"], "compound": mrow["compound"], "source": mrow["source"],
                                      "ModelID": v.index, "value": v.to_numpy(float)}))

    for (_, a), (_, b) in zip(m_sec.iterrows(), m_pri.iterrows()):
        done = False
        if a.screen_id is not None:
            vals = sec_long[sec_long.name == a.screen_id].set_index("depmap_id").auc
            n_feat = sum(1 for l in vals.dropna().index if l in feat_lines)
            if n_feat >= w6.MIN_LINES or b.screen_id is None:
                add(a, vals)
                done = True
        if not done and b.screen_id is not None:
            add(b, pri.loc[b.screen_id])
        if a.screen_id is None and b.screen_id is None:
            match_rows.append(dict(a, n_lines_response=0, n_lines_features=0, eligible=False))
    for _, r in m_g2.iterrows():
        if r.screen_id is None:
            match_rows.append(dict(r, n_lines_response=0, n_lines_features=0, eligible=False))
        else:
            add(r, g2[g2.DRUG_ID == r.screen_id].set_index("ModelID").AUC)
    for _, r in m_ct.iterrows():
        if r.screen_id is None:
            match_rows.append(dict(r, n_lines_response=0, n_lines_features=0, eligible=False))
        else:
            add(r, ct[r.screen_id])
    matches = pd.DataFrame(match_rows).merge(table[["compound", "class", "off_model_control"]], on="compound")
    C.atomic_write_tsv(matches, IN_DIR / "compound_matches.tsv")
    responses = pd.concat(resp, ignore_index=True)

    # ---- metabolomics
    met = pd.read_csv(METAB).drop(columns=["CCLE_ID"]).groupby("DepMap_ID").mean()
    met = met.loc[[l for l in met.index if l in set(lines)]]

    # ---- subtypes (task list generated here, before any fitting)
    mdl = pd.read_csv(MODEL_CSV, usecols=["ModelID", "OncotreeSubtype"]).set_index("ModelID")
    ms = pd.read_csv(MOLSUB).set_index("ModelID")
    lab_rows, task_rows = [], []
    lin = lines_df.lineage
    rna_lines = [l for l in lines if l in feat_lines]
    schemes = {"OncotreeSubtype": mdl.OncotreeSubtype.reindex(rna_lines)}
    pam = [c for c in ms.columns if c.startswith("PAM50")]
    P = ms[pam].reindex(rna_lines)
    one = P.fillna(False).astype(bool)
    pam_lab = pd.Series(np.where(one.sum(1) == 1, one.idxmax(1), None), index=rna_lines)
    pam_lab[P.isna().all(1)] = None
    schemes["Molecular:PAM50"] = pam_lab
    for c in ms.columns:
        if c in pam:
            continue
        s = ms[c].reindex(rna_lines)
        schemes[f"Molecular:{c}"] = s.map({True: "True", False: "False", "True": "True", "False": "False"})
    for scheme, lab in schemes.items():
        for L, idx in lab.groupby(lin.reindex(rna_lines)).groups.items():
            s = lab.loc[idx].dropna()
            vc = s.value_counts()
            keep = vc[vc >= w6.SUBTYPE_MIN_CLASS]
            if len(keep) >= 2:
                tid = f"{L}|{scheme}"
                s = s[s.isin(keep.index)]
                task_rows.append({"task_id": tid, "lineage": L, "scheme": scheme, "n_lines": len(s), "n_classes": len(keep),
                                  "classes": ";".join(f"{k}:{v}" for k, v in keep.sort_index().items())})
                lab_rows.append(pd.DataFrame({"task_id": tid, "ModelID": s.index, "label": s.to_numpy()}))
    tasks = pd.DataFrame(task_rows).sort_values("task_id")
    C.atomic_write_tsv(tasks, IN_DIR / "subtype_tasks.tsv",
                       header_lines=["generated by scripts/w6_inputs.py before any fitting (DESIGN W6 analysis 1)"])
    sub_labels = pd.concat(lab_rows, ignore_index=True)

    # ---- caches (per-line tables: tmp only)
    TMP.mkdir(parents=True, exist_ok=True)
    C.atomic_write_parquet(lines_df, TMP / "lines.parquet", register=False)
    C.atomic_write_parquet(rna, TMP / "rna.parquet", register=False)
    C.atomic_write_parquet(prot, TMP / "prot.parquet", register=False)
    C.atomic_write_parquet(rxn_rna, TMP / f"rxn_rna_{args.parent}.parquet", register=False)
    C.atomic_write_parquet(rxn_prot, TMP / f"rxn_prot_{args.parent}.parquet", register=False)
    C.atomic_write_parquet(responses, TMP / "responses.parquet", register=False)
    C.atomic_write_parquet(met.astype("float64"), TMP / "metabolomics.parquet", register=False)
    C.atomic_write_parquet(sub_labels, TMP / "subtype_labels.parquet", register=False)

    elig = matches[matches.eligible]
    cov = {
        "label": w6.PROVISIONAL_TAG,
        "cohort_lines": len(lines), "lines_with_rna": len(rna), "lines_with_procan": len(prot),
        "procan_proteins_used": prot.shape[1], "rxn_rna_reactions": rxn_rna.shape[1], "rxn_prot_reactions": rxn_prot.shape[1],
        "model_genes_with_rna": n_g_rna, "model_genes_with_protein": n_g_prot,
        "compounds_listed": len(table), "off_model_controls_listed": int(table.off_model_control.sum()),
        "units_eligible": {d: int((elig.dataset == d).sum()) for d in ("PRISM", "GDSC2", "CTRP")},
        "units_eligible_excl_off_model": {d: int(((elig.dataset == d) & ~elig.off_model_control).sum()) for d in ("PRISM", "GDSC2", "CTRP")},
        "compounds_with_any_eligible_unit": int(elig.compound.nunique()),
        "prism_sources": elig[elig.dataset == "PRISM"].source.value_counts().to_dict(),
        "lines_per_dataset_features_median": elig.groupby("dataset").n_lines_features.median().to_dict(),
        "metabolomics_lines": int(len(met)), "metabolomics_lines_with_rna": int(met.index.isin(feat_lines).sum()),
        "metabolites": int(met.shape[1]),
        "subtype_tasks": int(len(tasks)), "subtype_tasks_oncotree": int((tasks.scheme == "OncotreeSubtype").sum()),
        "subtype_lineages": int(tasks.lineage.nunique()),
        "target_genes": w6.target_genes(table),
    }
    C.atomic_write_json(cov, IN_DIR / "coverage.json")
    print(cov, flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", default="human1")
    ap.add_argument("--allow-dirty", action="store_true", help="label-free pilot only; recorded in the manifest")
    a = ap.parse_args()
    ins = [LINE_MEDIA, SPLIT, MODEL_CSV, MOLSUB, METAB, GDSC2, GDSC_CMP, GDSC_LINES, CTRP_AUC, CTRP_COND, PRISM_PRI,
           PRISM_PRI_LIST, PRISM_SEC, PROCAN, EXPR, C.PROJECT / w6.APRIORI_FILE]
    sys.exit(C.run_main(lambda: main(a), inputs=ins, solver="glpk", workers=1, allow_dirty=a.allow_dirty,
                        budget_cpu_h=0.25, note={"step": "inputs", "pilot": a.allow_dirty}))
