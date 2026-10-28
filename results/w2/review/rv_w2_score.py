"""Independent M4 W2 recomputation (reviewer). Does not import w2_score / w2_score_train.
Labels only via nextgem.labels (train tier) and cv.lolo (fold-scoped CE_fit). Usage: OMP_NUM_THREADS=1 python rv_w2_score.py"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from nextgem import labels, metrics, paths, w1, cv
from nextgem.config import SEEDS

W2 = Path("/data/tmpA/nextgem/results/w2"); OUT = W2 / "review"
B = 2000; S0 = SEEDS[0]

split = labels.load_split().set_index("ModelID")
lf = w1.line_features()
U = list(labels.gene_universe().entrez.astype(int))
elig = split[(split.tier == "train") & (split.stratum != "rare")]
E = sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])
lineage = split.lineage
lm = pd.read_csv(W2 / "media_library/line_media.tsv", sep="\t").set_index("ModelID")
A = [l for l in E if lm.stratum.get(l) == "attribution"]
print("n E", len(E), "n A", len(A), flush=True)
ls = pd.read_csv(W2 / "controls/line_sweeps.tsv", sep="\t")
sl = pd.read_csv(W2 / "ko/sweep_labels.tsv", sep="\t")

# H1 score matrix per sweep (own implementation)
ko = pd.read_parquet(paths.baseline_dir() / "human1_ko.parquet", columns=["gene", "entrez"]).dropna(subset=["entrez"])
ko["entrez"] = ko.entrez.astype(int); ko = ko.drop_duplicates("entrez").set_index("gene")
ratio = pd.read_parquet(W2 / "ko/ratio_r_human1.parquet")
r = ratio.reindex(columns=ko.index)
Esw = (1.0 - r).fillna(0.0).round(6); Esw.columns = ko.entrez.to_numpy()
Esw = Esw.reindex(columns=U, fill_value=0.0)

def sweeps(step):
    return ls[(ls.parent == "human1") & (ls.ladder_step == step)].set_index("ModelID").sweep_id
def S_of(sw, lines):
    M = Esw.loc[sw.reindex(lines).to_numpy()]; M.index = list(lines); return M

# labels
def fit_predict(fit_long, ids):
    return w1.fit_components(fit_long, lf.base_medium_class, labels_=("L_prob",)).ce["L_prob"], ids
folds = cv.lolo(fit_predict, lambda pr, held: pr, columns=["gene_effect", "L_prob"])
fold_names = list(folds); line_ce = {}
for k, (ce, ids) in folds.items():
    for l in ids: line_ce[l] = ce
tl = labels.train(columns=["gene_effect", "L_prob", "L_ge06"])
Y = w1.wide(tl, "L_prob", A, U); del tl
Ynoce = Y.where(w1.ce_mask(A, U, line_ce))
lin = lineage.reindex(A)

def lt(S, Yc): return metrics.per_line(S, Yc)
def delta(a, b):
    d = metrics.paired_line_delta(a, b); bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=B, seed=S0)
    return d, bt
def fpos(d):
    g = d.groupby(lineage.reindex(d.index)).agg(["mean", "size"])
    return int(sum((g.loc[f, "size"] >= 5) and (g.loc[f, "mean"] > 0) for f in fold_names if f in g.index)), {f: round(float(g.loc[f, "mean"]), 4) for f in g.index}

res = {}; T = {}
for step in ("M0", "M2(f=1)", "M-open"):
    S = S_of(sweeps(step), A)
    T[step] = {"ce": lt(S, Y), "noce": lt(S, Ynoce)}
for step in ("M2(f=1)", "M-open"):
    out = {}
    for c in ("ce", "noce"):
        d, bt = delta(T[step][c], T["M0"][c]); n, fm = fpos(d)
        out[c] = dict(delta=bt.estimate, lo=bt.ci_low, hi=bt.ci_high, folds_pos=n, n_lines=len(d), fold_means=fm)
    res[step] = out
    print(step, {c: (round(v["delta"], 4), round(v["lo"], 4), round(v["hi"], 4), v["folds_pos"]) for c, v in out.items()}, flush=True)

# shuffled control
doc_sw = sweeps("M2(f=1)")
lt_ = sl[(sl.parent == "human1") & (sl.variant == "M2") & np.isclose(sl.f.astype(float), 1.0) & (sl.gate_open.astype(str) == "True")]
lookup = lt_.drop_duplicates("medium_key").set_index("medium_key").sweep_id
assert (lm.medium_key.reindex(A).map(lookup).to_numpy() == doc_sw.reindex(A).to_numpy()).all(), "documented sweep != medium_key lookup"
sf = lm.serum_free.reindex(A).astype(str)
sh_author = pd.read_csv(W2 / "media_library/shuffled_media.tsv", sep="\t").set_index("ModelID")
def run_shuffle(keys_by_seed, tag):
    lts = []
    for ks in keys_by_seed:
        sw = ks.map(lookup); assert not sw.isna().any()
        lts.append(lt(S_of(sw, A), Y))
    m = lts[0].copy(); m["tap"] = np.mean([x.tap.to_numpy() for x in lts], axis=0)
    d, bt = delta(T["M2(f=1)"]["ce"], m)
    per_seed = [float(metrics.lineage_equal_mean(metrics.paired_line_delta(T["M2(f=1)"]["ce"], x), lineage)) for x in lts]
    sm0 = [float(metrics.lineage_equal_mean(metrics.paired_line_delta(x, T["M0"]["ce"]), lineage)) for x in lts]
    return dict(tag=tag, delta=bt.estimate, lo=bt.ci_low, hi=bt.ci_high, per_seed=per_seed, shuffled_minus_m0=sm0)
# (a) author's shuffled file, validated as a within-group permutation
ok = True
for sd in SEEDS[0:5]:
    col = sh_author[f"shuffled_seed{sd}"].reindex(A); grp = lm.lineage.reindex(A).astype(str) + "|" + sf
    for g, idx in grp.groupby(grp).groups.items():
        ok &= sorted(col.loc[idx]) == sorted(lm.medium_key.reindex(idx))
print("author shuffle is within-group permutation:", ok)
res["author_file_is_permutation"] = bool(ok)
res["shuffle_author_file"] = run_shuffle([sh_author[f"shuffled_seed{sd}"].reindex(A) for sd in SEEDS[0:5]], "author_file")
print(res["shuffle_author_file"], flush=True)
# (b) own independent shuffle
own = []
for sd in SEEDS[0:5]:
    rng = np.random.default_rng(sd); k = lm.medium_key.reindex(A).copy()
    grp = lm.lineage.reindex(A).astype(str) + "|" + sf
    for g, idx in grp.groupby(grp).groups.items():
        idx = list(idx); k.loc[idx] = rng.permutation(lm.medium_key.reindex(idx).to_numpy())
    own.append(k)
res["shuffle_own"] = run_shuffle(own, "own_rng")
print(res["shuffle_own"], flush=True)
# fraction of lines whose medium actually changed under shuffle
res["frac_changed_own"] = float(np.mean([(o != lm.medium_key.reindex(A)).mean() for o in own]))
(OUT / "rv_w2_score.json").write_text(json.dumps(res, indent=1, default=float))
