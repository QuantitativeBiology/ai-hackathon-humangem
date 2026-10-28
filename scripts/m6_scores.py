"""M6 label-free score vectors for every factorial cell and every line (training, validation AND test tiers).

Test-line scores are model outputs, not labels. They are written only under results/m6/scores/ and never next to a
label file. No label is read here (the split file, line media and the W1 generic table are label-free).

Per parent, the sweep e-matrix (sweep key x U; e = 1 − ratio_r, W1 Entrez convention, absent gene -> 0, 6 decimals)
collects:
  T0 | M0 and T0 | every M2(f=1) signature  : cached W2 sweeps (results/w2/ko)
  B-loo | M0 and B-loo | classic-RPMI M2(f=1): cached W3 sweeps (results/w3/sweeps)
  B-loo | other M2(f=1) signatures           : results/m6/sweeps (scripts/m6_sweeps.py)
A line's score under a cell is the row of the sweep (parent, variant, the line's medium at the cell's step); the
line -> medium mapping is W2's (results/w2/controls/line_sweeps.tsv, gate-free M2(f=1); imputed / fallback /
unmappable lines carry W2's integrated-pipeline medium). A failed sweep (none expected) would give its lines a NaN
row, which the scoring scripts replace by the §11.4 no-information score.

Outputs:
  results/m6/scores/sweep_e.parquet          long (key, parent, variant, w2_sweep_id, entrez, e): every cell sweep
  results/m6/scores/line_assignment.tsv      ModelID, parent, step, w2_sweep_id, medium_key, stratum (no tier column)
  results/m6/scores/scores_<tier>.large.parquet  line_gene layout (ModelID, entrez, one column per method) for all
                                             lines of the tier (git-ignored *.large.*; regenerable by this script)
  results/m6/scores/SHA256SUMS.tsv           SHA-256 + shape of every scores file (committed)
  results/m6/val/m6_val_scores.parquet       the validation look's scores (the 198 validation expression lines)
  results/m6/scores/w3_mev.tsv               W3-MEV diagnostic (A11b): HMGCR/MVK/PMVK/MVD KO ratio_r per cell sweep

Checks (abort otherwise): T0|M0 equals results/baseline bitwise; H1 / H2 columns equal the W1 generic scores; every
line of the split has a sweep for every cell; no NaN score unless its sweep failed.

Usage: .venv/bin/python scripts/m6_scores.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m6_common as mc  # noqa: E402
from m6_sweeps import W2, W3, reused_rows  # noqa: E402

TIERS = ("train", "validation", "test")


def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sweep_ratios() -> tuple[dict, pd.DataFrame]:
    """parent -> ratio_r frame (sweep key x model genes) over every key any cell uses; and the M6 job table."""
    reused = reused_rows()
    meta = pd.read_csv(mc.SWDIR / "sweeps.tsv", sep="\t")
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t")
    out = {}
    for p in mc.PARENTS:
        new = pd.read_parquet(mc.SWDIR / f"ratio_r_{p}.parquet").set_index("sweep")
        R = pd.concat([reused[p], new.reindex(columns=reused[p].columns)])
        keys = cat[cat.parent == p].key
        miss = sorted(set(keys) - set(R.index))
        if miss:
            raise RuntimeError(f"{p}: sweeps missing: {miss[:3]}")
        out[p] = R.loc[sorted(set(keys))]
    return out, meta


def e_matrices(U, ls: pd.DataFrame) -> tuple[dict, dict, pd.DataFrame]:
    """parent -> e-matrix (sweep key x U); parent -> Entrez map; the M6 job table. Asserts T0|M0 == baseline."""
    from nextgem import paths, w2_score as ws

    R, meta = sweep_ratios()
    E, emap = {}, {}
    for p in mc.PARENTS:
        ko = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez", "ratio_r"])
        emap[p] = ws.gene_to_entrez(ko)
        base = ko.set_index("gene").ratio_r
        m0 = sorted(set(mc.line_sweep_ids(ls, p, mc.M0)))
        k0 = mc.sweep_key(mc.T0, m0[0])
        if len(m0) != 1 or not np.array_equal(R[p].loc[k0, base.index].to_numpy(), base.to_numpy(), equal_nan=True):
            raise RuntimeError(f"{p}: T0|M0 differs from results/baseline")
        E[p] = mc.e_from_ratio(R[p], emap[p], U)
    return E, emap, meta


def mev_table(R: dict, split: pd.DataFrame, ls: pd.DataFrame) -> pd.DataFrame:
    """W3-MEV diagnostic (A11b; label-free): mevalonate-core KO ratio_r under every sweep of every cell, with the
    sweep's serum flag and its line counts per tier."""
    from nextgem import paths

    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t").set_index("key")
    rows = []
    for c in mc.CELLS:
        sym = pd.read_parquet(paths.baseline_dir() / f"{c.parent}_ko.parquet", columns=["gene", "symbol"])
        gid = {s: sym.gene[sym.symbol == s].iloc[0] for s in mc.MEVALONATE_CORE}
        keys = mc.cell_sweep_keys(c, ls)
        tier = split.tier.reindex(keys.index)
        for k in sorted(set(keys)):
            n = tier[keys == k].value_counts()
            r = {"parent": c.parent, "cell": c.key, "method": c.method, "key": k, "serum": bool(cat.loc[k, "serum"]),
                 "medium_keys": cat.loc[k, "medium_keys"], "n_train": int(n.get("train", 0)),
                 "n_validation": int(n.get("validation", 0)), "n_test": int(n.get("test", 0))}
            for s, g in gid.items():
                r[s] = float(R[c.parent].loc[k, g])
            r["core_lethal"] = bool(all(r[s] < 0.01 for s in mc.MEVALONATE_CORE))
            rows.append(r)
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    inputs = [W2 / "controls" / "line_sweeps.tsv", mc.SWDIR / "sweeps.tsv", mc.SWDIR / "catalog.tsv"]
    inputs += [W2 / "ko" / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    inputs += [W3 / "sweeps" / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    inputs += [mc.SWDIR / f"ratio_r_{p}.parquet" for p in mc.PARENTS]
    with mc.provenance(inputs=inputs, solver="none", workers=1, allow_dirty=a.allow_dirty) as prov:
        outs, note = run()
        prov.output(*outs)
        prov.note(**note)
    return 0


def run():
    from nextgem import labels, w1

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    gen = w1.generic_table(Udf.entrez)
    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    E, emap, meta = e_matrices(U, ls)
    mev = mev_table(sweep_ratios()[0], split, ls)
    outs = [mc.atomic_write_tsv(mev, mc.SCDIR / "w3_mev.tsv")]
    failed_keys = set(meta.loc[(meta.purpose == "new") & (meta.status != "ok"), "key"])

    # sweep-level e (long; keyed by sweep, not by line)
    rows = []
    for p in mc.PARENTS:
        X = E[p]
        long = X.stack(future_stack=True).rename("e").reset_index()
        long.columns = ["key", "entrez", "e"]
        long.insert(1, "parent", p)
        long.insert(2, "variant", long.key.str.split("|").str[0])
        long.insert(3, "w2_sweep_id", long.key.str.split("|").str[1])
        rows.append(long)
    se = pd.concat(rows, ignore_index=True)
    se["entrez"] = se.entrez.astype(np.int64)
    outs.append(mc.atomic_write_parquet(se, mc.SCDIR / "sweep_e.parquet"))

    # line -> medium assignment (label-free; no tier column)
    lm = pd.read_csv(W2 / "media_library" / "line_media.tsv", sep="\t").set_index("ModelID")
    la = []
    for p in mc.PARENTS:
        for step in mc.STEPS:
            s = ls[(ls.parent == p) & (ls.ladder_step == step)][["ModelID", "sweep_id", "medium_key", "stratum"]]
            la.append(s.rename(columns={"sweep_id": "w2_sweep_id"}).assign(parent=p, step=step))
    la = pd.concat(la, ignore_index=True)[["ModelID", "parent", "step", "w2_sweep_id", "medium_key", "stratum"]]
    if set(la.ModelID) != set(split.index):
        raise RuntimeError("line assignment does not cover exactly the split's lines")
    outs.append(mc.atomic_write_tsv(la.sort_values(["parent", "step", "ModelID"]), mc.SCDIR / "line_assignment.tsv"))

    # per-tier line_gene score files (all lines of the tier) + the validation look file (expression lines)
    methods = [c.method for c in mc.CELLS]
    c3 = mev[(mev.cell == "H1+Bloo/M2") & mev.serum]
    sums, note = [], {"failed_sweeps": sorted(failed_keys),
                      "w3_mev_cell3_serum_sweeps_core_lethal": f"{int(c3.core_lethal.sum())}/{len(c3)}"}

    def tier_frame(lines):
        df = pd.DataFrame({"ModelID": np.repeat(lines, len(U)), "entrez": np.tile(np.asarray(U, dtype=np.int64), len(lines))})
        for c in mc.CELLS:
            M = mc.cell_matrix(c, ls, E[c.parent], lines)
            v = M.to_numpy(dtype=np.float64)
            bad = np.isnan(v).any(axis=1)
            keys = mc.cell_sweep_keys(c, ls).reindex(lines)
            if bad.any() and not set(keys[bad]) <= failed_keys:
                raise RuntimeError(f"{c.key}: NaN scores outside failed sweeps")
            df[c.method] = v.ravel()
        return df

    for tier in TIERS:
        lines = sorted(split.index[split.tier == tier])
        df = tier_frame(lines)
        for s in ("H1", "H2"):  # generic columns = W1 generic scores (line-constant)
            g = np.tile(gen[f"e_{s}"].reindex(U).to_numpy(), len(lines))
            if not np.array_equal(df[s].to_numpy(), g):
                raise RuntimeError(f"{tier}: {s} differs from the W1 generic e_{s}")
        f = mc.atomic_write_parquet(df, mc.SCDIR / f"scores_{tier}.large.parquet")
        sums.append({"file": str(f.relative_to(mc.PROJECT)), "sha256": _sha(f), "tier": tier, "n_lines": len(lines),
                     "n_genes": len(U), "n_rows": len(df), "methods": ";".join(methods), "layout": "line_gene",
                     "tracked": False})
        if tier == "validation":
            V = [l for l in lines if lf.loc[l, "has_expression_matrix"]]
            dv = df[df.ModelID.isin(set(V))].reset_index(drop=True)
            fv = mc.atomic_write_parquet(dv, mc.VALDIR / "m6_val_scores.parquet")
            sums.append({"file": str(fv.relative_to(mc.PROJECT)), "sha256": _sha(fv), "tier": "validation (expression lines)",
                         "n_lines": len(V), "n_genes": len(U), "n_rows": len(dv), "methods": ";".join(methods),
                         "layout": "line_gene", "tracked": True})
            outs.append(fv)
            note["n_val_expr_lines"] = len(V)
        outs.append(f)
        note[f"n_{tier}_lines"] = len(lines)
    for f in ("sweep_e.parquet", "line_assignment.tsv"):
        p_ = mc.SCDIR / f
        sums.append({"file": str(p_.relative_to(mc.PROJECT)), "sha256": _sha(p_), "tier": "all (sweep-keyed)" if f.endswith(
            "parquet") else "all (assignment)", "n_lines": np.nan, "n_genes": len(U), "n_rows": len(pd.read_parquet(p_))
                     if f.endswith("parquet") else len(la), "methods": "", "layout": "sweep" if f.endswith("parquet")
                     else "assignment", "tracked": True})
    outs.append(mc.atomic_write_tsv(pd.DataFrame(sums), mc.SCDIR / "SHA256SUMS.tsv"))
    print(pd.DataFrame(sums)[["file", "tier", "n_lines", "n_rows", "sha256"]].to_string(index=False))
    return outs, note


if __name__ == "__main__":
    raise SystemExit(main())
