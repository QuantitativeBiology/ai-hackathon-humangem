"""Shared helpers for the W7 scripts (label-free). Inputs, idea registry checks, budget log, atomic writes."""

from __future__ import annotations

import datetime as dt
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "results" / "w7"
TMP = OUT / "tmp"  # git-ignored
VAL = OUT / "val"
MANIFEST = OUT / "MANIFEST.tsv"
IDEAS_TSV = PROJECT / "autoresearch" / "ideas.tsv"
BUDGET = PROJECT / "results" / "budget.tsv"
AUTHOR = "W7 auto-research engineer (Claude Opus)"
ART_VARIANT = "W7:artefact_bloo"
W1TMP = PROJECT / "results" / "w1" / "tmp"


def atomic_write_parquet(df: pd.DataFrame, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_parquet(tmp, index=False, compression="zstd")
    with open(tmp, "rb") as fh:
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return path


def atomic_write_tsv(df: pd.DataFrame, path: Path, float_format="%.6g") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False, float_format=float_format)
    os.replace(tmp, path)
    return path


def atomic_write_text(text: str, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
    return path


def frozen_label_hashes() -> dict:
    out = {}
    for ln in (PROJECT / "FROZEN.tsv").read_text().splitlines()[1:]:
        f = ln.split("\t")
        if len(f) >= 4 and f[1] == "label_file":
            out[f[0]] = f[3]
    return out


def logged_ideas() -> pd.DataFrame:
    return pd.read_csv(IDEAS_TSV, sep="\t")


def assert_logged(idea_id: str) -> None:
    """Protocol 1: an idea runs only after its row exists in autoresearch/ideas.tsv."""
    if idea_id not in set(logged_ideas().id):
        raise SystemExit(f"{idea_id} is not logged in {IDEAS_TSV}; log it before running")


def register_ideas() -> None:
    """Ideas beyond the seed set are registered in nextgem.w7.IDEAS (code) and logged in ideas.tsv (text)."""
    from nextgem import w7  # noqa: F401  (IDEAS lives in the module)


def log_budget(job: str, cpu_h: float, wall_h: float, note: str = "") -> None:
    row = "\t".join([dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%MZ"), "W7", AUTHOR, job,
                     f"{cpu_h:.3f}", f"{wall_h:.3f}", note]) + "\n"
    if not BUDGET.exists():
        BUDGET.write_text("date\tworkstream\tagent\tjob\tcpu_h\twall_h\tnote\n")
    with open(BUDGET, "a") as fh:
        fh.write(row)


def input_files() -> list[Path]:
    ins = [W1TMP / f"{k}.parquet" for k in ("expr", "cn", "mut_dam", "mut_hot")]
    ins += [PROJECT / "results" / "w1" / "line_features.tsv", PROJECT / "results" / "w1" / "gene_features.tsv",
            W1TMP / "gene_neighbours.json", PROJECT / "results" / "splits" / "lineage_split.tsv",
            PROJECT / "results" / "m6" / "scores" / "SHA256SUMS.tsv",
            PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv",
            PROJECT / "results" / "w2" / "media_library" / "line_media.tsv"]
    art = OUT / "sweeps" / "ratio_r_human1.parquet"
    if art.exists():
        ins.append(art)
    return ins


def art_matrix(lines, U) -> pd.DataFrame | None:
    """lines x U e-matrix of H1 + artefact-groups B-loo on the line's M2(f=1) medium (W7-03; label-free)."""
    sys.path.insert(0, str(PROJECT / "scripts"))
    import m6_common as mc

    from nextgem import paths, w2_score as ws

    f = OUT / "sweeps" / "ratio_r_human1.parquet"
    if not f.exists():
        return None
    R = pd.read_parquet(f).set_index("sweep")
    ko = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "entrez", "ratio_r"])
    E = mc.e_from_ratio(R, ws.gene_to_entrez(ko), U)
    ls = pd.read_csv(PROJECT / "results" / "w2" / "controls" / "line_sweeps.tsv", sep="\t")
    sid = mc.line_sweep_ids(ls, "human1", mc.M2).reindex(list(lines))
    keys = sid.map(lambda s: mc.sweep_key(ART_VARIANT, s))
    M = E.loc[keys.to_numpy()]
    M.index = list(lines)
    if M.isna().any().any():
        raise RuntimeError("artefact variant: NaN scores (failed sweep)")
    return M.astype("float32")


def deviation_features(omics: dict, aggs, genes, split: pd.DataFrame) -> dict:
    """Line-relative omics features (label-free): expression minus the gene's mean over training lines (gene norms
    from training lines only, so no validation/test line shapes a fitted feature), CN minus the line's median CN,
    and the same gene-centring for neighbour and paralog expression."""
    train = set(split.loc[split.tier == "train", "ModelID"])
    out = {}
    for name, M in (("expr_dev", omics["expr"]), ("neigh_dev", aggs[0]), ("paralog_dev", aggs[1])):
        M = M.reindex(columns=genes)
        mu = M.loc[[l for l in M.index if l in train]].mean(axis=0)
        out[name] = (M - mu).astype("float32")
    C = omics["cn"].reindex(columns=genes)
    out["cn_dev"] = C.sub(C.median(axis=1), axis=0).astype("float32")
    return out


def w2_medium_features(lines, genes) -> dict:
    """Line-constant W2 tiered medium annotation (label-free; results/w2/media_library/line_media.tsv): medium key
    (categorical code), serum %, serum-free flag, documented-medium flag."""
    lm = pd.read_csv(PROJECT / "results" / "w2" / "media_library" / "line_media.tsv", sep="\t").set_index("ModelID")
    lm = lm.reindex(list(lines))
    codes = {k: i for i, k in enumerate(sorted(lm.medium_key.dropna().unique()))}
    cols = {"w2_medium": lm.medium_key.map(codes).fillna(-1).to_numpy(float),
            "w2_serum_pct": lm.serum_pct.to_numpy(float),
            "w2_serum_free": lm.serum_free.map({True: 1.0, False: 0.0}).to_numpy(float),
            "w2_documented": lm.medium_source.astype(str).str.startswith("documented").astype(float).to_numpy()}
    return {k: pd.DataFrame(np.repeat(v[:, None], len(genes), 1).astype("float32"), index=list(lines), columns=genes)
            for k, v in cols.items()}


def build_inputs(tiers=("train",), art: bool = True):
    """Label-free Inputs for nextgem.w7 (all lines of ``tiers`` for the GEM features)."""
    from nextgem import labels, w1, w7

    U = labels.gene_universe()
    genes = list(U.entrez.astype(int))
    omics = w1.load_omics()
    lf = w1.line_features()
    gf = w1.gene_features()
    nb = json.loads((W1TMP / "gene_neighbours.json").read_text())
    aggs = w1.expression_aggregates(omics["expr"], genes, nb["neighbours_h2"], nb["paralogs"])
    class_codes = {c: i for i, c in enumerate(sorted(lf.base_medium_class.dropna().unique()))}
    gem_line = w7.load_gem_line(tiers)
    if art:
        A = art_matrix(list(gem_line["e_P1"].index), genes)
        if A is not None:
            gem_line["e_ART"] = A
            gem_line["art_lethal"] = (A >= w7.CALL).astype("float32")
    split = labels.load_split()
    gem_line.update(deviation_features(omics, aggs, genes, split))
    gem_line.update(w2_medium_features(sorted(set(gem_line["e_P1"].index) | set(omics["expr"].index)), genes))
    rx = TMP / "h1_reaction_counts.parquet"
    if not rx.exists():
        atomic_write_parquet(w7.human1_reaction_counts().reset_index(), rx)
    rxn = pd.read_parquet(rx).set_index("entrez")
    gen = w1.generic_table(U.entrez)
    return w7.Inputs(genes=genes, omics=omics, aggs=aggs, gfeat=gf, lf=lf, class_codes=class_codes, gem_line=gem_line,
                     gem_gene=w7.gene_gem_features(genes, gen, rxn), tier_of=split.set_index("ModelID").tier,
                     lineage_of=split.set_index("ModelID").lineage)
