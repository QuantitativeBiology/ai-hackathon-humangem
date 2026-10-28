"""W6 analyses, controls and decision statistics (DESIGN §9 W6, §6.5).

    .venv/bin/python scripts/w6_analyze.py --config provisional [--pilot] [--allow-dirty]

Analyses (same folds, matched capacity: every block reduced to PCA-20 (PCA-50 sensitivity) fitted on training lines
only, ridge with lineage one-hot covariate for drugs/metabolites, kNN for subtypes):
  1. subtypes: within-lineage tasks from results/w6/inputs/subtype_tasks.tsv (written before fitting); stratified
     5-fold, SEEDS[0:3]; kNN macro-F1 (decision metric) and silhouette;
  2. metabolic drugs: units = (compound of notes/apriori_lists_final.md §4, dataset) with >= 100 lines with features;
     LOLO held-out Spearman rho (pooled; within-lineage mean as secondary); off-model controls reported separately;
     fit-free within-lineage rho of the target's e and expression vs sensitivity; positive-control pairs;
  3. metabolomics: the 225 CCLE metabolites, LOLO held-out rho.
Assay media: drugs per dataset (PRISM SOP medium, GDSC Screen Medium, CTRP culture medium); metabolomics and subtypes
use the line's culture medium.

Controls: flux permuted within lineage (SEEDS[0:5]); synthetic positive (tanh-injected signal carried only by the flux
block); synthetic negative (a block that is a random linear map of expression).

Decision (per arm): one-sided Wilcoxon of the per-unit Δ per analysis, Holm across the 3 analyses (alpha 0.025
one-sided), and the real median Δ above every permuted-seed median, for
  H_flux: Δ = metric(RNA + X_flux) - metric(RNA);
  H_sim:  Δ = metric(X_flux) - metric(X_rxnRNA)  (and X_KO - X_rxnRNA).
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # one BLAS thread per process (<= 6 worker processes in total; DESIGN §11.4)

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w6_common as C  # noqa: E402

from nextgem import paths, w6  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

PILOT_UNITS = {"drugs": 6, "metabolomics": 10, "subtypes": 10, "synthetic": 4}


# ------------------------------------------------------------------------------------------------ blocks
def omics_blocks(tmp_in: Path, parent: str):
    rna = pd.read_parquet(tmp_in / "rna.parquet").astype("float64")
    prot = pd.read_parquet(tmp_in / "prot.parquet").astype("float64")
    rxn_rna = pd.read_parquet(tmp_in / f"rxn_rna_{parent}.parquet").astype("float64")
    rxn_prot = pd.read_parquet(tmp_in / f"rxn_prot_{parent}.parquet").astype("float64")
    from nextgem.genes import parse_depmap_columns

    gt = pd.read_csv(paths.data_root() / ("Human1/model/genes.tsv" if parent == "human1" else "Human-GEM-main/model/genes.tsv"),
                     sep="\t", dtype=str)
    model_entrez = set(gt.geneEntrezID.dropna().str.split(";").str[0])
    dm = parse_depmap_columns(rna.columns)
    model_cols = list(dm.column[dm.entrez.isin(model_entrez)])
    neg = w6.linear_of_expression(rna[model_cols], SEEDS[0])
    blocks = {"RNA": w6.Block.from_matrix("RNA", rna), "PROT": w6.Block.from_matrix("PROT", prot),
              "RXN_RNA": w6.Block.from_matrix("RXN_RNA", rxn_rna), "RXN_PROT": w6.Block.from_matrix("RXN_PROT", rxn_prot),
              "NEG": w6.Block.from_matrix("NEG", neg)}
    sym_col = {c.rsplit(" (", 1)[0]: c for c in rna.columns}
    return blocks, rna, sym_col


def arm_blocks(feat: Path, tmp: Path, arm: str):
    st = pd.read_csv(feat / "states.tsv", sep="\t", comment="#")
    st = st[st.arm == arm].set_index("state")
    mu = st.wt_growth
    fl = pd.read_parquet(feat / "flux_long.parquet")
    fl = fl[fl.state.isin(st.index)]
    F = fl.pivot_table(index="state", columns="rxn", values="v", fill_value=0.0).reindex(st.index).fillna(0.0)
    F = np.arcsinh(F.div(mu, axis=0))
    ko = 1.0 - pd.read_parquet(feat / "ko_ratio.parquet").T.reindex(st.index)  # file: genes x states
    out = {"FLUX": (F, w6.Block.from_matrix(f"{arm}:FLUX", F)), "KO": (ko, w6.Block.from_matrix(f"{arm}:KO", ko))}
    kp = tmp / "koflux_long.parquet"
    if kp.exists():
        kl = pd.read_parquet(kp)
        kl = kl[kl.state.isin(st.index)]
        kl["col"] = kl.target + ":" + kl.rxn
        K = kl.pivot_table(index="state", columns="col", values="v", fill_value=0.0).reindex(st.index).fillna(0.0)
        K = np.arcsinh(K.div(mu, axis=0))
        out["KOFLUX"] = (K, w6.Block.from_matrix(f"{arm}:KOFLUX", K))
    R = pd.DataFrame(1.0, index=st.index, columns=["all_reactions_present"])  # generic: every reaction present
    out["RXN"] = (R, w6.Block.from_matrix(f"{arm}:RXN", R))
    return out, ko


class Views:
    """All views (block x assay map x permutation) with their PCA universes and a shared PCA cache."""

    def __init__(self, omics: dict, arms: dict, maps: dict, medium_state: pd.DataFrame, lineage: pd.Series,
                 rna_lines: list):
        self.omics, self.arms, self.maps, self.lineage = omics, arms, maps, lineage
        self.ms = medium_state
        self.rna_lines = set(rna_lines)
        self.v: dict[str, w6.View] = {}
        self.universe: dict[str, list] = {}
        self.cache = w6.PCACache()

    def omic(self, name: str) -> w6.View:
        tag = f"{name}|-|-"
        if tag not in self.v:
            b = self.omics[name]
            s = pd.Series(b.states, index=b.states)
            self.v[tag] = w6.View(b, s, tag)
            self.universe[tag] = list(b.states)
        return self.v[tag]

    def arm(self, arm: str, block: str, amap: str, perm: int | None = None) -> w6.View:
        tag = f"{arm}:{block}|{amap}|{perm if perm is not None else '-'}"
        if tag not in self.v:
            b = self.arms[arm][block][1]
            m = self.maps[amap]
            st = self.ms[self.ms.arm == arm].set_index("medium_key").state
            s = m.map(st).dropna()
            s = s[s.index.isin(self.rna_lines)]
            if perm is not None:
                s = w6.permute_within_lineage(s, self.lineage, perm)
            self.v[tag] = w6.View(b, s, tag)
            self.universe[tag] = list(s.index)
        return self.v[tag]

    def lines(self, views) -> set:
        out = None
        for v in views:
            s = set(v.state_of.index)
            out = s if out is None else out & s
        return out or set()


# ------------------------------------------------------------------------------------------------ units
def feature_sets(V: Views, arm: str, amap: str, analysis: str) -> dict[str, list]:
    """name -> list of views. Primary: RNA, RNA+FLUX, FLUX, RXN_RNA; the rest secondary or controls."""
    fs = {
        "RNA": [V.omic("RNA")],
        "RNA+FLUX": [V.omic("RNA"), V.arm(arm, "FLUX", amap)],
        "RNA+KO": [V.omic("RNA"), V.arm(arm, "KO", amap)],
        "RNA+RXN_RNA": [V.omic("RNA"), V.omic("RXN_RNA")],
        "RNA+RXN": [V.omic("RNA"), V.arm(arm, "RXN", amap)],
        "FLUX": [V.arm(arm, "FLUX", amap)],
        "KO": [V.arm(arm, "KO", amap)],
        "RXN_RNA": [V.omic("RXN_RNA")],
        "RNA+NEG": [V.omic("RNA"), V.omic("NEG")],
    }
    if analysis == "drugs" and "KOFLUX" in V.arms[arm]:
        fs["RNA+KOFLUX"] = [V.omic("RNA"), V.arm(arm, "KOFLUX", amap)]
        fs["KOFLUX"] = [V.arm(arm, "KOFLUX", amap)]
    for s in w6.PERM_SEEDS:
        fs[f"RNA+FLUX@perm{s}"] = [V.omic("RNA"), V.arm(arm, "FLUX", amap, s)]
        fs[f"FLUX@perm{s}"] = [V.arm(arm, "FLUX", amap, s)]
        fs[f"KO@perm{s}"] = [V.arm(arm, "KO", amap, s)]
    return fs


PROT_SETS = ("PROT", "PROT+FLUX", "RXN_PROT", "PROT+RXN_PROT")


def prot_sets(V: Views, arm: str, amap: str) -> dict[str, list]:
    return {"PROT": [V.omic("PROT")], "PROT+FLUX": [V.omic("PROT"), V.arm(arm, "FLUX", amap)],
            "RXN_PROT": [V.omic("RXN_PROT")], "PROT+RXN_PROT": [V.omic("PROT"), V.omic("RXN_PROT")]}


def run_regression_unit(y: pd.Series, V: Views, arm: str, amap: str, analysis: str, k_sens: bool = True) -> list[dict]:
    rows = []
    fs = feature_sets(V, arm, amap, analysis)
    base_lines = V.lines(fs["RNA+FLUX"])
    yy = y[y.index.isin(base_lines)].dropna()
    for name, views in fs.items():
        pred = w6.lolo_ridge(yy, views, V.universe, V.lineage, V.cache, k=w6.PCA_K)
        rows.append({"set": name, "k": w6.PCA_K, "n": len(yy), **w6.heldout_scores(yy, pred, V.lineage)})
    if k_sens:
        for name in ("RNA", "RNA+FLUX", "FLUX", "RXN_RNA"):
            pred = w6.lolo_ridge(yy, fs[name], V.universe, V.lineage, V.cache, k=w6.PCA_K_SENS)
            rows.append({"set": name, "k": w6.PCA_K_SENS, "n": len(yy), **w6.heldout_scores(yy, pred, V.lineage)})
    ps = prot_sets(V, arm, amap)
    yp = y[y.index.isin(V.lines(ps["PROT+FLUX"]))].dropna()
    if len(yp) >= w6.MIN_LINES:
        for name, views in ps.items():
            pred = w6.lolo_ridge(yp, views, V.universe, V.lineage, V.cache, k=w6.PCA_K)
            rows.append({"set": name, "k": w6.PCA_K, "n": len(yp), **w6.heldout_scores(yp, pred, V.lineage)})
    return rows


def run_subtype_unit(labels: pd.Series, V: Views, arm: str) -> list[dict]:
    rows = []
    fs = feature_sets(V, arm, "culture", "subtypes")
    lab = labels[labels.index.isin(V.lines(fs["RNA+FLUX"]))]
    vc = lab.value_counts()
    lab = lab[lab.isin(vc[vc >= w6.SUBTYPE_MIN_CLASS].index)]
    if lab.nunique() < 2:
        return rows
    cache: dict = {}
    for name, views in fs.items():
        rows.append({"set": name, "k": w6.PCA_K, "n": len(lab), **w6.subtype_scores(lab, views, cache=cache)})
    ps = prot_sets(V, arm, "culture")
    lp = lab[lab.index.isin(V.lines(ps["PROT+FLUX"]))]
    vc = lp.value_counts()
    lp = lp[lp.isin(vc[vc >= w6.SUBTYPE_MIN_CLASS].index)]
    if lp.nunique() >= 2:
        for name, views in ps.items():
            rows.append({"set": name, "k": w6.PCA_K, "n": len(lp), **w6.subtype_scores(lp, views, cache=cache)})
    return rows


# ------------------------------------------------------------------------------------------------ decisions
def deltas(units: pd.DataFrame, metric: str, a: str, b: str, k: int = w6.PCA_K) -> pd.Series:
    u = units[units.k == k].pivot_table(index="unit", columns="set", values=metric)
    if a not in u or b not in u:
        return pd.Series(dtype=float)
    return (u[a] - u[b]).dropna()


def decision_table(per_analysis: dict[str, tuple[pd.DataFrame, str]], arm: str) -> pd.DataFrame:
    out = []
    specs = [
        ("H_flux", "X_flux", "RNA+FLUX", "RNA", lambda s: f"RNA+FLUX@perm{s}", "RNA"),
        ("H_sim", "X_flux", "FLUX", "RXN_RNA", lambda s: f"FLUX@perm{s}", "RXN_RNA"),
        ("H_sim", "X_KO", "KO", "RXN_RNA", lambda s: f"KO@perm{s}", "RXN_RNA"),
    ]
    for hyp, block, a, b, permname, pb in specs:
        d, perm = {}, {}
        for an, (units, metric) in per_analysis.items():
            d[an] = deltas(units, metric, a, b)
            perm[an] = {s: deltas(units, metric, permname(s), pb) for s in w6.PERM_SEEDS}
        t = w6.decision_rows(d, perm)
        t.insert(0, "contrast", f"{a} - {b}")
        t.insert(0, "block", block)
        t.insert(0, "hypothesis", hyp)
        t.insert(0, "arm", arm)
        t["metric"] = [per_analysis[an][1] for an in t.analysis]
        out.append(t)
    return pd.concat(out, ignore_index=True)


def secondary_table(per_analysis: dict[str, tuple[pd.DataFrame, str]], arm: str) -> pd.DataFrame:
    rows = []
    pairs = [("RNA+KO", "RNA"), ("RNA+KOFLUX", "RNA"), ("RNA+RXN_RNA", "RNA"), ("RNA+RXN", "RNA"), ("RNA+NEG", "RNA"),
             ("PROT+FLUX", "PROT"), ("PROT+RXN_PROT", "PROT"), ("KOFLUX", "RXN_RNA")]
    for an, (units, metric) in per_analysis.items():
        for a, b in pairs:
            d = deltas(units, metric, a, b)
            if len(d):
                rows.append({"arm": arm, "analysis": an, "contrast": f"{a} - {b}", "k": w6.PCA_K, "n_units": len(d),
                             "median_delta": float(d.median()), "p_wilcoxon_one_sided": w6.wilcoxon_greater(d)})
        for a, b in [("RNA+FLUX", "RNA"), ("FLUX", "RXN_RNA")]:
            d = deltas(units, metric, a, b, k=w6.PCA_K_SENS)
            if len(d):
                rows.append({"arm": arm, "analysis": an, "contrast": f"{a} - {b}", "k": w6.PCA_K_SENS, "n_units": len(d),
                             "median_delta": float(d.median()), "p_wilcoxon_one_sided": w6.wilcoxon_greater(d)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ main
def main(a) -> int:
    t0 = time.time()
    cfg, cfg_path = C.load_config(a.config)
    parent = cfg.arms[0].parent
    run_tag = cfg.tag + ("_pilot" if a.pilot else "")
    feat = C.run_dir(run_tag) / "features"
    out = C.run_dir(run_tag) / "analysis"
    tmp_in, tmp = C.tmp_dir("inputs"), C.tmp_dir(run_tag)
    label = cfg.label + (" — PILOT (cost/plumbing check, not a result)" if a.pilot else "")
    hdr = [f"{label}; config {cfg_path.name}"]

    lines_df = pd.read_parquet(tmp_in / "lines.parquet")
    lineage = lines_df.lineage
    omics, rna, sym_col = omics_blocks(tmp_in, parent)
    ms = pd.read_csv(feat / "medium_state.tsv", sep="\t", comment="#")
    arms, ko_tables = {}, {}
    for arm in cfg.arms:
        arms[arm.name], ko_tables[arm.name] = arm_blocks(feat, tmp, arm.name)
    lines = list(lines_df.index)
    maps = {k: w6.assay_medium_map(k, lines, lines_df.medium_key, lines_df.gdsc_screen_medium) for k in w6.ASSAY_MAPS}
    V = Views(omics, arms, maps, ms, lineage, list(rna.index))

    responses = pd.read_parquet(tmp_in / "responses.parquet")
    met = pd.read_parquet(tmp_in / "metabolomics.parquet")
    sub = pd.read_parquet(tmp_in / "subtype_labels.parquet")
    tasks = pd.read_csv(C.OUT_DIR / "inputs" / "subtype_tasks.tsv", sep="\t", comment="#")
    table = pd.read_csv(C.OUT_DIR / "inputs" / "compounds.tsv", sep="\t", comment="#").set_index("compound")
    units_drug = responses[["dataset", "compound"]].drop_duplicates().sort_values(["dataset", "compound"])
    metabolites = list(met.columns)
    task_ids = list(tasks.task_id)
    if a.pilot:
        units_drug = units_drug.groupby("dataset").head(2).head(PILOT_UNITS["drugs"])
        metabolites = metabolites[: PILOT_UNITS["metabolomics"]]
        task_ids = task_ids[: PILOT_UNITS["subtypes"]]

    all_dec, all_sec, summ = [], [], {"label": label, "arms": {}}
    timing = {}
    for arm in cfg.arms:
        an = arm.name
        ta = time.time()
        # ---- 2. drugs
        rows = []
        for _, u in units_drug.iterrows():
            y = responses[(responses.dataset == u.dataset) & (responses.compound == u.compound)].set_index("ModelID").value
            amap = w6.DATASET_MAP[u.dataset]
            for r in run_regression_unit(y, V, an, amap, "drugs"):
                rows.append({"unit": f"{u.dataset}|{u.compound}", "dataset": u.dataset, "compound": u.compound,
                             "off_model_control": bool(table.loc[u.compound, "off_model_control"]), **r})
        drugs = pd.DataFrame(rows)
        timing[f"{an}:drugs"] = time.time() - ta
        # ---- 3. metabolomics
        tb = time.time()
        rows = []
        for mname in metabolites:
            for r in run_regression_unit(met[mname], V, an, "culture", "metabolomics"):
                rows.append({"unit": mname, **r})
        metab = pd.DataFrame(rows)
        timing[f"{an}:metabolomics"] = time.time() - tb
        # ---- 1. subtypes
        tc = time.time()
        rows = []
        for tid in task_ids:
            lab = sub[sub.task_id == tid].set_index("ModelID").label
            for r in run_subtype_unit(lab, V, an):
                rows.append({"unit": tid, **r})
        subt = pd.DataFrame(rows)
        timing[f"{an}:subtypes"] = time.time() - tc
        # ---- synthetic positive control (culture map, metabolomics lines, flux-only nonlinear signal)
        td = time.time()
        F = arms[an]["FLUX"][0]
        sv = V.arm(an, "FLUX", "culture")
        syn_lines = [l for l in met.index if l in sv.state_of.index]
        ss = np.random.SeedSequence(SEEDS[0]).spawn(w6.N_SYNTH)
        rows = []
        for i, sq in enumerate(ss[: (PILOT_UNITS["synthetic"] if a.pilot else w6.N_SYNTH)]):
            seed = int(sq.generate_state(1)[0])
            z = w6.random_state_direction(F, seed)
            y = w6.synthetic_positive_response(z, sv.state_of.loc[syn_lines], seed)
            for nm in ("RNA", "RNA+FLUX"):
                views = [V.omic("RNA")] if nm == "RNA" else [V.omic("RNA"), sv]
                pred = w6.lolo_ridge(y, views, V.universe, V.lineage, V.cache)
                rows.append({"unit": f"synthetic{i}", "set": nm, "k": w6.PCA_K, "n": len(y), **w6.heldout_scores(y, pred, V.lineage)})
        synth = pd.DataFrame(rows)
        timing[f"{an}:synthetic"] = time.time() - td

        drugs_primary = drugs[~drugs.off_model_control]
        per_analysis = {"subtypes": (subt, "knn_macro_f1"), "drugs": (drugs_primary, "rho"), "metabolomics": (metab, "rho")}
        dec = decision_table(per_analysis, an)
        all_dec.append(dec)
        sec = secondary_table(per_analysis, an)
        sec_off = secondary_table({"drugs_off_model": (drugs[drugs.off_model_control], "rho")}, an)
        d_off = deltas(drugs[drugs.off_model_control], "rho", "RNA+FLUX", "RNA")
        extra = pd.DataFrame([{"arm": an, "analysis": "drugs_off_model", "contrast": "RNA+FLUX - RNA", "k": w6.PCA_K,
                               "n_units": len(d_off), "median_delta": float(d_off.median()) if len(d_off) else np.nan,
                               "p_wilcoxon_one_sided": w6.wilcoxon_greater(d_off)}])
        d_sil = deltas(subt, "silhouette", "RNA+FLUX", "RNA")
        sil = pd.DataFrame([{"arm": an, "analysis": "subtypes(silhouette)", "contrast": "RNA+FLUX - RNA", "k": w6.PCA_K,
                             "n_units": len(d_sil), "median_delta": float(d_sil.median()) if len(d_sil) else np.nan,
                             "p_wilcoxon_one_sided": w6.wilcoxon_greater(d_sil)}])
        d_syn = deltas(synth, "rho", "RNA+FLUX", "RNA")
        syn_row = pd.DataFrame([{"arm": an, "analysis": "synthetic_positive", "contrast": "RNA+FLUX - RNA", "k": w6.PCA_K,
                                 "n_units": len(d_syn), "median_delta": float(d_syn.median()) if len(d_syn) else np.nan,
                                 "p_wilcoxon_one_sided": w6.wilcoxon_greater(d_syn)}])
        all_sec.append(pd.concat([sec, sec_off, extra, sil, syn_row], ignore_index=True))

        # ---- fit-free test and positive-control pairs
        rows = []
        ko = ko_tables[an]
        s2g = w6.symbol_map(paths.data_root() / ("Human1/model/genes.tsv" if parent == "human1" else "Human-GEM-main/model/genes.tsv"), list(ko.columns))
        st_map = ms[ms.arm == an].set_index("medium_key").state
        for _, u in units_drug.iterrows():
            y = responses[(responses.dataset == u.dataset) & (responses.compound == u.compound)].set_index("ModelID").value
            y = y[y.index.isin(rna.index)]
            sens = -y
            lmap = maps[w6.DATASET_MAP[u.dataset]].map(st_map).dropna()
            yl = sens[sens.index.isin(lmap.index)]
            tg = [t for t in str(table.loc[u.compound, "targets"]).split(";") if t and t != "nan"]
            for t in tg:
                g = s2g.get(t)
                e = pd.Series(ko[g].reindex(lmap.loc[yl.index]).to_numpy(), index=yl.index) if g in ko.columns else None
                rho_e, nl_e = w6.fit_free_rho(e, yl, lineage) if e is not None else (np.nan, 0)
                ex = rna[sym_col[t]].reindex(yl.index) if t in sym_col else None
                rho_x, nl_x = w6.fit_free_rho(ex, yl, lineage) if ex is not None else (np.nan, 0)
                e_const = bool(e is not None and e.groupby(lineage.reindex(e.index)).nunique().max() <= 1)
                rows.append({"dataset": u.dataset, "compound": u.compound, "target": t, "n_lines": len(yl),
                             "rho_e_within": rho_e, "n_lineages_e": nl_e, "e_constant_within_every_lineage": e_const,
                             "rho_expr_within": rho_x, "n_lineages_expr": nl_x})
        ff = pd.DataFrame(rows)
        pcp = w6.positive_control_pairs_long()
        held = drugs[(drugs.k == w6.PCA_K) & drugs.set.isin(["RNA", "RNA+FLUX"])].pivot_table(
            index=["dataset", "compound"], columns="set", values="rho").reset_index()
        pc = pcp.merge(ff.rename(columns={"target": "gene"}), on=["gene", "compound"], how="left").merge(
            held, on=["dataset", "compound"], how="left")
        pc["status"] = np.where(pc.dataset.isna(), "no eligible screen unit",
                                np.where(pc.e_constant_within_every_lineage.fillna(True).astype(bool),
                                         "non_informative (target e constant within lineages)", "evaluated"))
        sub_dir = out / an.replace("+", "_")
        for name, df in (("units_drugs", drugs), ("units_metabolomics", metab), ("units_subtypes", subt),
                         ("units_synthetic_positive", synth), ("fit_free", ff), ("positive_control_pairs", pc)):
            C.atomic_write_tsv(df, sub_dir / f"{name}.tsv", header_lines=hdr + [f"arm {an}"])
        summ["arms"][an] = {"n_units": {"drugs": int(drugs.unit.nunique()), "drugs_primary": int(drugs_primary.unit.nunique()),
                                        "metabolomics": int(metab.unit.nunique()) if len(metab) else 0,
                                        "subtypes": int(subt.unit.nunique()) if len(subt) else 0},
                            "fit_free_e_defined": int(ff.rho_e_within.notna().sum()) if len(ff) else 0,
                            "fit_free_rows": len(ff)}
        print(f"[{an}] done in {time.time() - ta:.0f}s", flush=True)

    dec = pd.concat(all_dec, ignore_index=True)
    sec = pd.concat(all_sec, ignore_index=True)
    C.atomic_write_tsv(dec, out / "decisions.tsv", header_lines=hdr + [
        "one-sided Wilcoxon of per-unit delta; Holm across the 3 analyses within each (arm, hypothesis, block); "
        f"alpha {w6.ALPHA_ONE_SIDED}; 'above_permuted' = real median > every permuted-seed median (SEEDS[0:5])"])
    C.atomic_write_tsv(sec, out / "secondary.tsv", header_lines=hdr + ["exploratory / controls; no multiplicity correction"])
    summ["timing_s"] = {k: round(v, 1) for k, v in timing.items()}
    summ["wall_s"] = round(time.time() - t0, 1)
    summ["pca_cache_entries"] = len(V.cache._c)
    C.atomic_write_json(summ, out / "summary.json")
    print(dec.to_string(), flush=True)
    print(sec.to_string(), flush=True)
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="provisional")
    ap.add_argument("--pilot", action="store_true", help="a few units per analysis (cost/plumbing check)")
    ap.add_argument("--budget-cpu-h", type=float, default=2.0)
    ap.add_argument("--allow-dirty", action="store_true", help="label-free pilot only; recorded in the manifest")
    a = ap.parse_args()
    cfg, cfg_path = C.load_config(a.config)
    feat = C.run_dir(cfg.tag + ("_pilot" if a.pilot else "")) / "features"
    ins = [cfg_path] + [C.OUT_DIR / "inputs" / f for f in ("compounds.tsv", "subtype_tasks.tsv", "compound_matches.tsv")]
    ins += [feat / f for f in ("states.tsv", "flux_long.parquet", "ko_ratio.parquet", "medium_state.tsv")]
    sys.exit(C.run_main(lambda: main(a), inputs=ins, seeds=list(SEEDS[0:5]), solver="none", workers=1,
                        allow_dirty=a.allow_dirty, budget_cpu_h=a.budget_cpu_h,
                        note={"step": "analyze", "config": cfg.tag, "label": cfg.label, "pilot": a.pilot}))
