"""W5 KO sweeps of the curated and control models on M0 (label-free, GLPK).

Arms per parent (ratio_r for every model gene; local re-evaluation of the genes of edited GPRs on the curation base,
every other gene keeps its base value; the base sweep is a full sweep, checked against results/baseline):
  generic           released parent (= results/baseline)
  base              curation base: H2 released; H1 with MAR04644 irreversible ("H1_TYMS", DESIGN §3.1)
  W5a               W5a curated GPRs (results/w5/curated_W5a_<parent>.tsv) on the curation base
  W5a_pid30         sensitivity: W5a with paralogs restricted to >= 30% identity (declared before scoring)
  W5a_gpr_only      (H1 only) W5a GPRs on the released H1 (decomposition: GPR part of Δ_curation)
  neg_s<seed>       negative control (DESIGN W5): random single edits of matched number and type, one per GPR group,
                    among target groups with no W5a edit (non-protected OR nodes), SEEDS[0:5]; falls back to
                    non-target groups only if the target pool is too small (logged)
Positive control (no KO needed): recall of the Human1 v1.19.0 -> Human2 v2.1.1 GPR changes on the H1 target
reactions by W5a(H1) (possibly label-informed: Luo et al. may have tuned GPRs on DepMap; DESIGN W5).

Outputs: results/w5/ko_<parent>.parquet (one column per arm, baseline-file row order), results/w5/gene_changes_<parent>.tsv,
results/w5/neg_control_edits.tsv, results/w5/pos_control_v1v2.tsv, results/w5/apply_summary.json
Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w5_apply.py [--workers 6]
"""

from __future__ import annotations

import argparse
import json
import sys

import cobra

cobra.Configuration().solver = "glpk"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem import paths, w5  # noqa: E402
from nextgem.config import SEEDS  # noqa: E402

def local_arm(eng, gp, base, new_trees, workers):
    return w5.local_sweep(eng, gp, base, new_trees, workers)


def curated_trees(parent: str, gp: dict, variant: str = "W5a") -> dict:
    c = pd.read_csv(w5.w5_dir() / f"curated_{variant}_{parent}.tsv", sep="\t")
    new = {}
    for r in c[c.edited].itertuples():
        t = w5.parse_rule(r.curated_rule)
        for rx in r.reactions.split(";"):
            assert w5.to_rule(gp[rx]) == r.rule, (rx, r.rule)
            new[rx] = t
    return new


def random_edits(parent, gp, seed, n_by_kind: dict, prot) -> tuple[dict, list]:
    """Matched random single edits (one per group) among target groups without a W5a edit."""
    tg = pd.read_csv(w5.w5_dir() / f"targets_{parent}.tsv", sep="\t")
    cur = pd.read_csv(w5.w5_dir() / f"curated_W5a_{parent}.tsv", sep="\t")
    edited = set(cur[cur.edited].group)
    rng = np.random.default_rng([seed, 5, 0 if parent == "human1" else 1])
    pools = {"target": tg[tg.target & ~tg.group.isin(edited)], "nontarget": tg[~tg.target]}
    used, new, log = set(), {}, []
    for kind in sorted(n_by_kind):
        need = n_by_kind[kind]
        for pool_name in ("target", "nontarget"):
            if need == 0:
                break
            cands = []
            for r in pools[pool_name].itertuples():
                if r.group in used:
                    continue
                rx = r.reactions.split(";")
                tree = gp[rx[0]]
                for e in w5.enumerate_edits(tree):
                    if e.kind != kind:
                        continue
                    node = w5.get_node(tree, e.path)
                    if node[0] == "or" and prot.or_node_protected(rx, node):
                        continue
                    cands.append((r.group, r.reactions, e))
            order = rng.permutation(len(cands))
            for i in order:
                if need == 0:
                    break
                grp, rxs, e = cands[i]
                if grp in used:
                    continue
                used.add(grp)
                rx = rxs.split(";")
                t2 = w5.apply_edit(gp[rx[0]], e)
                for x in rx:
                    new[x] = t2
                log.append({"parent": parent, "seed": seed, "kind": kind, "group": grp, "pool": pool_name, "edit": e.key(),
                            "old_rule": w5.to_rule(gp[rx[0]]), "new_rule": w5.to_rule(t2), "reactions": rxs})
                need -= 1
        if need:
            log.append({"parent": parent, "seed": seed, "kind": kind, "group": None, "pool": "SHORTFALL", "edit": str(need)})
    return new, log


def main(prov, a) -> int:
    t = w5.Timer()
    out = w5.w5_dir()
    summary, neg_log, changes = {}, [], []
    for parent in w5.PARENTS:
        bl = pd.read_parquet(paths.baseline_dir() / f"{parent}_ko.parquet")
        prov.input(paths.baseline_dir() / f"{parent}_ko.parquet", out / f"curated_W5a_{parent}.tsv",
                   out / f"targets_{parent}.tsv", out / f"edits_W5a_{parent}.tsv")
        m = w5.load_parent(parent, tyms_irreversible=(parent == "human1"))
        eng = w5.KOEngine(m)
        gp = w5.model_gprs(m)
        genes = list(bl.gene)
        base = w5.full_sweep(eng, gp, genes)
        b0 = pd.read_parquet(w5.tmp_dir() / f"base_{parent}.parquet").set_index("gene").ratio_r
        assert (base - b0.reindex(base.index)).abs().max() == 0, "base sweep differs from w5_enumerate"
        arms = {"generic": bl.set_index("gene").ratio_r.reindex(genes), "base": base}
        new = curated_trees(parent, gp)
        arms["W5a"] = local_arm(eng, gp, base, new, a.workers)
        arms["W5a_pid30"] = local_arm(eng, gp, base, curated_trees(parent, gp, "W5a_pid30"), a.workers)
        prov.input(out / f"curated_W5a_pid30_{parent}.tsv")
        if parent == "human1":
            m0 = w5.load_parent(parent)
            eng0 = w5.KOEngine(m0)
            gp0 = w5.model_gprs(m0)
            arms["W5a_gpr_only"] = local_arm(eng0, gp0, arms["generic"], new, a.workers)
        edits = pd.read_csv(out / f"edits_W5a_{parent}.tsv", sep="\t")
        n_by_kind = edits.kind.value_counts().to_dict() if len(edits) else {}
        prot = w5.resolve_protected(parent)
        for sd in SEEDS[0:5]:
            rnd, log = random_edits(parent, gp, sd, n_by_kind, prot)
            neg_log += log
            arms[f"neg_s{sd}"] = local_arm(eng, gp, base, rnd, a.workers)
        df = pd.DataFrame({"gene": genes, "entrez": bl.entrez.to_numpy(), "symbol": bl.symbol.to_numpy(),
                           **{f"ratio_{k}": v.reindex(genes).to_numpy() for k, v in arms.items()}})
        f = w5.atomic_parquet(df, out / f"ko_{parent}.parquet")
        prov.output(f)
        # gene-level changes (label-free)
        for k in [c for c in df.columns if c.startswith("ratio_") and c not in ("ratio_generic",)]:
            d = df[df[k] != df.ratio_generic]
            for r in d.itertuples():
                changes.append({"parent": parent, "arm": k[6:], "gene": r.gene, "entrez": r.entrez, "symbol": r.symbol,
                                "ratio_generic": r.ratio_generic, "ratio_base": r.ratio_base, "ratio_arm": getattr(r, k),
                                "call_generic": r.ratio_generic < w5.LETHAL, "call_arm": getattr(r, k) < w5.LETHAL})
        summary[parent] = {
            "wt_growth_base": eng.wt, "n_genes": len(genes), "n_lp": eng.n_lp, "n_retry": eng.n_retry,
            "W5a_reactions_edited": len(new), "W5a_edits_by_kind": n_by_kind,
            **{f"n_changed_vs_generic_{k}": int((df[f"ratio_{k}"] != df.ratio_generic).sum()) for k in arms if k != "generic"},
            **{f"n_flip_vs_generic_{k}": int(((df[f"ratio_{k}"] < w5.LETHAL) != (df.ratio_generic < w5.LETHAL)).sum())
               for k in arms if k != "generic"},
            **{f"n_lethal_{k}": int((df[f"ratio_{k}"] < w5.LETHAL).sum()) for k in arms},
        }
    gc = pd.DataFrame(changes)
    nl = pd.DataFrame(neg_log)
    prov.output(w5.atomic_tsv(gc, out / "gene_changes.tsv"), w5.atomic_tsv(nl, out / "neg_control_edits.tsv"))
    # ---- positive control: H1 -> H2 GPR changes on H1 target reactions
    pc, pcs = positive_control()
    prov.output(w5.atomic_tsv(pc, out / "pos_control_v1v2.tsv"))
    summary["positive_control_v1v2"] = pcs
    summary["runtime"] = t.lap()
    (out / "apply_summary.json").write_text(json.dumps(summary, indent=1, default=str))
    prov.output(out / "apply_summary.json")
    print(json.dumps(summary, indent=1, default=str))
    return 0


def positive_control():
    from nextgem.model import load_model

    h1 = load_model("human1")
    h2 = load_model("human2")
    g1, g2 = w5.model_gprs(h1), w5.model_gprs(h2)
    tg = pd.read_csv(w5.w5_dir() / "targets_human1.tsv", sep="\t")
    cur = pd.read_csv(w5.w5_dir() / "curated_W5a_human1.tsv", sep="\t").set_index("group")
    rows = []
    for r in tg[tg.target].itertuples():
        for rx in r.reactions.split(";"):
            t1 = g1[rx]
            t2 = g2.get(rx)
            cr = w5.parse_rule(cur.loc[r.group, "curated_rule"])
            changed = t2 is None or w5.to_rule(t2) != w5.to_rule(t1)
            rows.append({"reaction": rx, "group": r.group, "h1_rule": w5.to_rule(t1), "h2_rule": w5.to_rule(t2) if t2 else "",
                         "changed_v1_v2": changed,
                         "change_class": ("removed_or_no_gpr" if t2 is None else
                                          "logic_only" if changed and w5.genes_of(t1) == w5.genes_of(t2) else
                                          "gene_set" if changed else "unchanged"),
                         "W5a_edited": w5.to_rule(cr) != w5.to_rule(t1),
                         "W5a_equals_h2": t2 is not None and w5.to_rule(cr) == w5.to_rule(t2)})
    df = pd.DataFrame(rows)
    ch = df[df.changed_v1_v2]
    lo = ch[ch.change_class == "logic_only"]
    s = {"caveat": "possibly label-informed (Human2 GPR curation may have used DepMap; DESIGN W5, §11.2)",
         "h1_target_reactions": len(df), "changed_v1_v2": len(ch), "changed_logic_only": len(lo),
         "changed_gene_set": int((ch.change_class == "gene_set").sum()),
         "recall_exact_all_changed": float(ch.W5a_equals_h2.mean()) if len(ch) else None,
         "recall_exact_logic_only": float(lo.W5a_equals_h2.mean()) if len(lo) else None,
         "recall_any_edit_changed": float(ch.W5a_edited.mean()) if len(ch) else None,
         "W5a_edited_reactions": int(df.W5a_edited.sum()),
         "W5a_edited_that_equal_h2": int((df.W5a_edited & df.W5a_equals_h2).sum())}
    return df, s


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args()
    sys.exit(w5.run_main(lambda prov: main(prov, args), solver="glpk", workers=args.workers, seeds=list(SEEDS[0:5])))
