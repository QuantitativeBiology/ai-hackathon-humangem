"""W5a evidence-driven GPR edits (DESIGN §9 W5 table), label-free, on the top-200 target GPR groups of each parent.

Reads results/w5/targets_<parent>.tsv (scripts/w5_enumerate.py) and the local evidence files (Ensembl 116 paralogs,
CORUM 5.1, Complex Portal, ProCan 20250211, Nusinow 2020, DepMap 26Q3 expression via results/w1/tmp/expr.parquet).
Never edits the isozyme exceptions of notes/apriori_lists_final.md §3. De Kegel 2021 / Ito 2021 are not read.

Outputs (results/w5/):
  edits_W5a_<parent>.tsv        applied edits, one row each, with evidence ids and details
  curated_W5a_<parent>.tsv      target groups: old rule, curated rule, reactions
  or_nodes_<parent>.tsv         every OR node examined (pair evidence flags; W5b input)
  and_paralog_pairs_<parent>.tsv every AND-ed paralog sibling pair examined
  blocked_<parent>.tsv          candidate edits not applied (protected relation, non-clique, all branches absent)
  evidence_validation.tsv       one row per (edit, evidence id): exists locally?
  apriori_targets_check.tsv     the §3 "curation targets (not exceptions)" and §1.3 RNR item: targeted? edited?
  curation_summary.json
  w5b_candidates.json       W5b candidate OR nodes (paths in the curated trees, gene pairs)

Usage: .venv/bin/python scripts/w5_curate.py
"""

from __future__ import annotations

import json
import sys

import cobra

cobra.Configuration().solver = "glpk"

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem import w5, w5_evidence as ev  # noqa: E402

# notes/apriori_lists_final.md §3 "Curation targets (not exceptions)" and §1.3 "Model GPR error to fix" (checked only)
APRIORI_TARGETS = [
    ("TECR or TECRL (MAR02194)", ["MAR02194"], ["TECR", "TECRL"]),
    ("PTGR1 or HSD17B12 or TECRL (MAR01074)", ["MAR01074"], ["PTGR1", "HSD17B12", "TECRL"]),
    ("ALDH3A2/ALDH2/ALDH1B1/ALDH9A1/ALDH7A1 block", [], ["ALDH3A2", "ALDH2", "ALDH1B1", "ALDH9A1", "ALDH7A1"]),
    ("SEC11A and SEC11C", [], ["SEC11A", "SEC11C"]),
    ("MAR06435/MAR06436", ["MAR06435", "MAR06436"], []),
    ("RNR AND (RRM1 and (RRM2 or RRM2B))", ["MAR02352", "MAR02357", "MAR04611", "MAR04612", "MAR04614"], ["RRM2", "RRM2B"]),
    ("ATP5MC1/2/3 triplet (AND -> OR)", [], ["ATP5MC1", "ATP5MC2", "ATP5MC3"]),
    ("HMGCS1/HMGCS2 relation", [], ["HMGCS1", "HMGCS2"]),
]


# W5a = literal DESIGN rule (any Ensembl 116 paralog pair). W5a_pid30 = sensitivity declared before any scoring: a
# "paralog" must have >= 30% identity in both directions (W1's paralog convention), because the literal rule merges
# distant 'other_paralog' pairs (e.g. POLR1A/POLR2A/POLR3A, 21-25% identity). W5b builds on W5a.
VARIANTS = (("W5a", 0.0), ("W5a_pid30", ev.PARALOG_PID_SENSITIVITY))


def main(prov) -> int:
    t = w5.Timer()
    out = w5.w5_dir()
    models = {p: w5.load_parent(p) for p in w5.PARENTS}
    db = ev.build_db(models)
    for f in ev.evidence_inputs():
        prov.input(f)
    summary = {"evidence": {
        "model_genes_union": len(db.genes), "paralog_pairs_among_model_genes": len(db.paralog),
        "paralog_pairs_by_type": pd.Series([v[0] for v in db.paralog.values()]).value_counts().to_dict(),
        "genes_in_any_complex": sum(1 for v in db.complexes.values() if v),
        "complex_members_unmapped": db.unmapped_complex_members,
        "procan_lines": int(db.procan.shape[0]), "procan_proteins_model": int(db.procan.notna().any().sum()),
        "nusinow_lines": int(db.nusinow.shape[0]), "nusinow_proteins_model": int(db.nusinow.notna().any().sum()),
        "expr_lines": int(db.expr_n.max()), "transport_genes": len(db.transport_genes)}}
    all_edits, w5b, ap_rows, files = [], {}, [], []
    for variant, min_pid in VARIANTS:
        db.min_pid = min_pid
        summary[variant] = {}
        for parent in w5.PARENTS:
            tg = pd.read_csv(out / f"targets_{parent}.tsv", sep="\t")
            prov.input(out / f"targets_{parent}.tsv") if variant == "W5a" else None
            prot = w5.resolve_protected(parent)
            cur = ev.Curation()
            rows = []
            if variant == "W5a":
                w5b[parent] = {}
            for r in tg[tg.target].itertuples():
                rx = r.reactions.split(";")
                tree = w5.parse_rule(r.rule)
                ctx = {"parent": parent, "group": int(r.group), "rank": int(r.rank), "reactions": r.reactions}
                new = ev.curate_tree(tree, rx, db, prot, cur, ctx)
                rows.append({"parent": parent, "group": int(r.group), "rank": int(r.rank), "reactions": r.reactions,
                             "n_reactions": len(rx), "rule": r.rule, "curated_rule": w5.to_rule(new),
                             "edited": w5.to_rule(new) != r.rule})
                cands = ev.w5b_candidates(new, rx, db, prot) if variant == "W5a" else []
                if cands:
                    w5b[parent][str(int(r.group))] = {"curated_rule": w5.to_rule(new), "reactions": rx,
                                                      "candidates": [{**c, "path": list(c["path"]), "pairs": [list(p) for p in c["pairs"]]}
                                                                     for c in cands]}
            cdf = pd.DataFrame(rows)
            edf = pd.DataFrame(cur.edits)
            if len(edf):
                edf.insert(0, "edit_id", [f"{variant}_{w5.SHORT[parent]}_{i + 1:03d}" for i in range(len(edf))])
            gt = w5.gene_table(parent)
            sym = lambda s: ";".join(str(gt.symbol.get(g, g)) for g in str(s).split(";"))  # noqa: E731
            if len(edf):
                edf["symbols"] = edf.genes.map(sym)
            all_edits.append(edf)
            files += [w5.atomic_tsv(edf, out / f"edits_{variant}_{parent}.tsv"),
                      w5.atomic_tsv(cdf, out / f"curated_{variant}_{parent}.tsv"),
                      w5.atomic_tsv(pd.DataFrame(cur.or_nodes), out / f"or_nodes_{variant}_{parent}.tsv"),
                      w5.atomic_tsv(pd.DataFrame(cur.and_pairs), out / f"and_paralog_pairs_{variant}_{parent}.tsv"),
                      w5.atomic_tsv(pd.DataFrame(cur.removals_blocked), out / f"blocked_{variant}_{parent}.tsv")]
            orn = pd.DataFrame(cur.or_nodes)
            summary[variant][parent] = {
                "target_groups": int(tg.target.sum()), "target_reactions": int(tg[tg.target].n_reactions.sum()),
                "groups_edited": int(cdf.edited.sum()), "reactions_edited": int(cdf[cdf.edited].n_reactions.sum()),
                "edits_by_kind": edf.kind.value_counts().to_dict() if len(edf) else {},
                "or_nodes_examined": len(orn), "or_nodes_protected": int((orn.protected != "").sum()) if len(orn) else 0,
                "or_nodes_all_nonparalog": int(orn.all_nonparalog.sum()) if len(orn) else 0,
                "or_nodes_all_cocomplex": int(orn.all_cocomplex.sum()) if len(orn) else 0,
                "or_nodes_all_coab": int(orn.all_coab.sum()) if len(orn) else 0,
                "and_paralog_pairs": len(cur.and_pairs), "blocked": len(cur.removals_blocked),
                "w5b_candidate_or_nodes": sum(len(v["candidates"]) for v in w5b.get(parent, {}).values()) if variant == "W5a" else None,
                "w5b_candidate_or_nodes_anchored": sum(sum(c["anchored"] for c in v["candidates"]) for v in w5b.get(parent, {}).values()) if variant == "W5a" else None,
                "protected_unresolved": prot.unresolved,
            }
            # a priori curation targets: are they among the targets, and did W5a edit them?
            sym2g = {}
            for g, s in gt.symbol.items():
                sym2g.setdefault(s, set()).add(g)
            allg = pd.read_csv(out / f"targets_{parent}.tsv", sep="\t")
            for name, rxs, syms in APRIORI_TARGETS:
                gs = set().union(*(sym2g.get(s, set()) for s in syms)) if syms else set()
                hit = allg[allg.reactions.str.split(";").map(lambda v: bool(set(v) & set(rxs))) |
                           allg.rule.map(lambda r: len(gs) > 1 and len(set(w5.genes_of(w5.parse_rule(r))) & gs) >= 2)]
                ed = edf[edf.group.isin(hit.group)] if len(edf) else edf
                ap_rows.append({"variant": variant, "parent": parent, "item": name, "n_groups": len(hit), "best_rank": int(hit["rank"].min()) if len(hit) else None,
                                "n_target_groups": int(hit.target.sum()), "n_groups_with_flip": int((hit.n_flip_genes > 0).sum()),
                                "edits": ";".join(f"{k}:{s}" for k, s in zip(ed.kind, ed.symbols)) if len(ed) else ""})
    edits = pd.concat([e for e in all_edits if len(e)], ignore_index=True) if any(len(e) for e in all_edits) else pd.DataFrame()
    val = ev.validate(edits, db.genes) if len(edits) else pd.DataFrame(columns=["edit_id", "evidence_id", "ok", "check"])
    files.append(w5.atomic_tsv(val, out / "evidence_validation.tsv"))
    files.append(w5.atomic_tsv(pd.DataFrame(ap_rows), out / "apriori_targets_check.tsv"))
    (w5.w5_dir() / "w5b_candidates.json").write_text(json.dumps(w5b, indent=0, default=bool))
    files.append(w5.w5_dir() / "w5b_candidates.json")
    summary["validation"] = {"n_ids": len(val), "n_ok": int(val.ok.sum()) if len(val) else 0,
                             "n_failed": int((~val.ok).sum()) if len(val) else 0}
    summary["runtime"] = t.lap()
    (out / "curation_summary.json").write_text(json.dumps(summary, indent=1, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)))
    files.append(out / "curation_summary.json")
    prov.output(*files)
    print(json.dumps(summary, indent=1, default=str))
    if summary["validation"]["n_failed"]:
        print("EVIDENCE VALIDATION FAILED", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(w5.run_main(main, solver="glpk", workers=1))
