"""W5 label-free targets (DESIGN §9 W5 "Targets"): enumerate every single GPR edit per GPR group, re-evaluate only the
genes of the edited GPR (deduplicated closed sets, compute SF4), rank groups by the number of genes whose binary
call (ratio_r < 0.01, M0) flips, and keep the top 200.

Parents: Human2 v2.1.1 as released; Human1 v1.19.0 on its curation base (MAR04644 irreversible, DESIGN §3.1,
because every optimised H1 carries it). The engine is checked first against results/baseline (generic parents).
GLPK only (declared solver); no labels are read.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w5_enumerate.py --parent human2 [--workers 6]
Outputs: results/w5/targets_<parent>.tsv, results/w5/tmp/enum_<parent>.parquet, results/w5/enum_<parent>.json,
         results/w5/tmp/base_<parent>.parquet (base ratio_r per gene).
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys

import cobra

cobra.Configuration().solver = "glpk"  # before any model is loaded (gurobi_pool.py docstring)

import pandas as pd  # noqa: E402

from nextgem import paths, w5  # noqa: E402

G: dict = {}


def _solve(cs):
    return cs, G["eng"].ratio(cs)


def run_lps(eng, closed_sets, workers):
    todo = [cs for cs in closed_sets if cs not in eng.cache
            and not all(abs(eng.flux[r]) <= w5.ZERO_FLUX_TOL for r in cs) and cs]
    todo = sorted(set(todo), key=lambda s: (len(s), sorted(s)))
    G["eng"] = eng
    if workers > 1 and len(todo) > 50:
        with mp.get_context("fork").Pool(workers) as pool:
            for cs, v in pool.imap_unordered(_solve, todo, chunksize=8):
                eng.cache[cs] = v
    else:
        for cs in todo:
            eng.ratio(cs)
    return len(todo)


def main(prov, a) -> int:
    t = w5.Timer()
    parent = a.parent
    tyms = parent == "human1"
    out = w5.w5_dir()
    # --- engine check on the released parent
    m0 = w5.load_parent(parent)
    eng0 = w5.KOEngine(m0)
    gp0 = w5.model_gprs(m0)
    genes = sorted(g.id for g in m0.genes)
    base0 = w5.full_sweep(eng0, gp0, genes)
    bl = pd.read_parquet(paths.baseline_dir() / f"{parent}_ko.parquet").set_index("gene").ratio_r
    prov.input(paths.baseline_dir() / f"{parent}_ko.parquet")
    diff = (base0 - bl.reindex(base0.index)).abs()
    check = {"max_abs_diff_vs_baseline": float(diff.max()), "n_differ": int((diff > 0).sum()), "n_genes": len(genes)}
    if check["n_differ"]:
        raise AssertionError(f"KO engine differs from results/baseline: {check}")
    # --- curation base
    if tyms:
        m = w5.load_parent(parent, tyms_irreversible=True)
        eng = w5.KOEngine(m)
    else:
        m, eng = m0, eng0
    gp = w5.model_gprs(m)
    base = w5.full_sweep(eng, gp, genes)
    cl = w5.closures(gp)
    groups = w5.groups_of(gp)
    gt = w5.gene_table(parent)
    prot = w5.resolve_protected(parent)
    # --- enumerate
    rows = []
    needed = []
    for gi, (rule, rxns) in enumerate(groups.items()):
        tree = gp[rxns[0]]
        R = frozenset(rxns)
        gg = w5.genes_of(tree)
        for e in w5.enumerate_edits(tree):
            new = w5.apply_edit(tree, e)
            node = w5.get_node(tree, e.path)
            prot_reason = prot.or_node_protected(rxns, node) if node[0] == "or" else None
            for g in sorted(gg):
                c0 = cl.get(g, frozenset())
                c1 = w5.edited_closure(c0, R, new, g)
                if c1 == c0:
                    continue
                rows.append((gi, rule, e.key(), e.kind, w5.to_rule(new), g, c1, prot_reason))
                needed.append(c1)
    n_lp = run_lps(eng, needed, a.workers)
    df = pd.DataFrame(rows, columns=["group", "rule", "edit", "kind", "new_rule", "gene", "closed", "protected"])
    df["n_closed"] = df.closed.map(len)
    df["ratio_new"] = df.closed.map(eng.ratio)
    df["ratio_base"] = df.gene.map(base)
    df["call_base"] = df.ratio_base < w5.LETHAL
    df["call_new"] = df.ratio_new < w5.LETHAL
    df["flip"] = df.call_base != df.call_new
    df["ratio_change"] = df.ratio_new != df.ratio_base
    df["entrez"] = df.gene.map(gt.entrez)
    df["symbol"] = df.gene.map(gt.symbol)
    df = df.drop(columns="closed")
    # --- rank groups
    edits_per_group = {gi: len(w5.enumerate_edits(gp[rx[0]])) for gi, (r, rx) in enumerate(groups.items())}
    recs = []
    for gi, (rule, rxns) in enumerate(groups.items()):
        sub = df[df.group == gi]
        fl = sub[sub.flip]
        per_edit = fl.groupby("edit").gene.nunique()
        recs.append({
            "group": gi, "rule": rule, "reactions": ";".join(rxns), "n_reactions": len(rxns),
            "n_genes": len(w5.genes_of(gp[rxns[0]])), "n_edits": edits_per_group[gi],
            "n_flip_genes": fl.gene.nunique(), "n_flip_genes_U": fl[fl.entrez.notna()].gene.nunique(),
            "max_flip_single_edit": int(per_edit.max()) if len(per_edit) else 0,
            "n_ratio_change_genes": sub[sub.ratio_change].gene.nunique(),
            "flip_genes": ";".join(sorted(set(fl.symbol.fillna(fl.gene)))),
            "edit_kinds_flipping": ";".join(sorted(set(fl.kind))),
            "protected_or_nodes": int(sub.protected.notna().any()),
        })
    tg = pd.DataFrame(recs)
    tg = tg.sort_values(["n_flip_genes", "max_flip_single_edit", "n_ratio_change_genes", "reactions"],
                        ascending=[False, False, False, True], kind="mergesort").reset_index(drop=True)
    tg["rank"] = range(1, len(tg) + 1)
    # top 200, restricted to groups where some single edit changes some gene's KO ratio: below that point every
    # group is tied at (0, 0, 0) and the order would be arbitrary (reaction id); see README "Deviations"
    tg["impact"] = (tg.n_flip_genes > 0) | (tg.n_ratio_change_genes > 0)
    tg["target"] = (tg["rank"] <= w5.N_TARGETS) & tg.impact
    cut = [int(x) for x in tg.loc[w5.N_TARGETS - 1, ["n_flip_genes", "max_flip_single_edit", "n_ratio_change_genes"]]] \
        if len(tg) >= w5.N_TARGETS else None
    out.mkdir(parents=True, exist_ok=True)
    f_t = w5.atomic_tsv(tg, out / f"targets_{parent}.tsv")
    f_e = w5.atomic_parquet(df, w5.tmp_dir() / f"enum_{parent}.parquet")
    bdf = pd.DataFrame({"gene": base.index, "ratio_r": base.values, "entrez": base.index.map(gt.entrez),
                        "symbol": base.index.map(gt.symbol)})
    f_b = w5.atomic_parquet(bdf, w5.tmp_dir() / f"base_{parent}.parquet")
    summ = {
        "parent": parent, "curation_base": "MAR04644 irreversible (0,1000)" if tyms else "as released",
        "engine_check_vs_baseline": check, "wt_growth_base": eng.wt, "n_reactions_with_gpr": len(gp),
        "n_gpr_groups": len(groups), "n_edits": int(sum(edits_per_group.values())),
        "n_edits_by_kind": {k: int(v) for k, v in pd.Series([e.kind for r, rx in groups.items() for e in w5.enumerate_edits(gp[rx[0]])]).value_counts().items()},
        "n_edit_gene_rows_changed_closure": len(df), "n_unique_lps_solved": n_lp, "n_flip_rows": int(df.flip.sum()),
        "n_groups_with_flip": int((tg.n_flip_genes > 0).sum()),
        "n_targets": int(tg.target.sum()), "target_reactions": int(tg[tg.target].n_reactions.sum()),
        "cut_key_at_rank_200": cut,
        "n_groups_tied_with_cut": int((tg[["n_flip_genes", "max_flip_single_edit", "n_ratio_change_genes"]].values.tolist().count(cut))) if cut else 0,
        "n_groups_with_any_ratio_change": int(tg.impact.sum()),
        "n_base_lethal": int((base < w5.LETHAL).sum()), "protected_unresolved_symbols": prot.unresolved,
        "workers": a.workers, **t.lap(),
    }
    (out / f"enum_{parent}.json").write_text(json.dumps(summ, indent=1, default=str))
    prov.output(f_t, f_e, f_b, out / f"enum_{parent}.json")
    print(json.dumps(summ, indent=1, default=str))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parent", choices=w5.PARENTS, required=True)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    sys.exit(w5.run_main(lambda prov: main(prov, args), solver="glpk", workers=args.workers))
