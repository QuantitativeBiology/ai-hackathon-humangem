"""F1 step 1 (label-free): enumerate the D17 GPR edits on Human1 v1.19.0 with their evidence (A15, post hoc).

Rule: scripts/f1_common.py (module docstring). Evidence: Ensembl 116 paralogs (%id both directions), Human-GEM
genes.tsv compartments (Human1), CORUM 5.1 + Complex Portal (no shared complex), read with the W5 readers
(``nextgem.w5_evidence``; complex membership mapped through the H1 ∪ H2 gene tables exactly as W5). Every reaction
with a GPR is examined (not only W5's top-200 targets: on P1's media the M0 target ranking does not apply).

Outputs (results/f1/):
  edits_human1.tsv            one row per applied edit (AND -> OR clique merges, plus the TYMS bounds edit)
  curated_rules_human1.tsv    rxn, old_rule, new_rule for every reaction whose GPR changes
  paralog_pairs_human1.tsv    every AND-ed leaf-sibling Ensembl paralog pair in every reaction, with the pair
                              criteria, the reaction exclusions and whether a merge was applied
  excluded_reactions_human1.tsv  every GPR reaction excluded by E1-E4, with all reasons
  lumped_reactions_human1.tsv every reaction (with or without GPR) flagged pseudo/lumped by E3 (the label-free list)
  evidence_check_human1.tsv   evidence re-check from the raw files (independent code paths)
  closure_changes_human1.tsv  genes whose single-gene KO closure set changes under the curated GPRs
  edits_summary.json          counts (incl. the merges excluded by E3/E4/E2, POLR vs other)

Usage: .venv/bin/python scripts/f1_edits.py [--allow-dirty]
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import f1_common as fc  # noqa: E402


def gene_tables():
    from nextgem import w5

    g1 = w5.gene_table("human1")
    union = pd.concat([g1, w5.gene_table("human2")])
    union = union[~union.index.duplicated(keep="last")]  # as w5_evidence.build_db (H2 rows win)
    return g1, union


def enumerate_edits(model, g1, paralog, complexes):
    from nextgem import w5

    def facts(a, b):
        v = paralog.get(frozenset((a, b)))
        sc = tuple(sorted(complexes.get(a, set()) & complexes.get(b, set())))
        return fc.PairFacts(paralog=v is not None, ptype=v[0] if v else None, pid=v[1] if v else float("nan"),
                            shared_complexes=sc, comp_a=fc.comp_set(g1.compartments.get(a)),
                            comp_b=fc.comp_set(g1.compartments.get(b)))

    recs, excl_rows, lumped_rows, curated = [], [], [], []
    for r in model.reactions:
        stoich = {m.id: c for m, c in r.metabolites.items()}
        comps = {m.compartment for m in r.metabolites}
        lr = fc.lumped_reason(r.subsystem, stoich)
        if lr:
            lumped_rows.append({"rxn": r.id, "name": r.name, "subsystem": r.subsystem, "reason": lr,
                                "n_genes": len(r.genes), "equation": r.reaction})
        if not r.gene_reaction_rule.strip():
            continue
        tree = w5.parse_rule(r.gene_reaction_rule)
        ex = fc.reaction_exclusions(boundary=r.boundary, subsystem=r.subsystem, stoich=stoich, met_comps=comps,
                                    tree=tree)
        if ex:
            excl_rows.append({"rxn": r.id, "name": r.name, "subsystem": r.subsystem, "n_genes": len(r.genes),
                              "max_and_genes": fc.max_and_genes(tree), "exclusions": "; ".join(ex)})
        ctx = {"rxn": r.id, "rxn_name": r.name, "subsystem": r.subsystem, "rxn_excluded": "; ".join(ex)}
        new = fc.curate_tree(tree, facts, recs, ctx, apply=not ex)
        if new != tree:
            curated.append({"rxn": r.id, "old_rule": r.gene_reaction_rule, "new_rule": w5.to_rule(new),
                            "old_canonical": w5.to_rule(tree)})
    return recs, pd.DataFrame(excl_rows), pd.DataFrame(lumped_rows), pd.DataFrame(curated)


def build_edit_table(recs, g1, model) -> pd.DataFrame:
    from nextgem.w5_evidence import pair_key

    sym = g1.symbol.to_dict()
    rows = []
    merges, seen = [], {}
    for r in recs:  # one edit per (reaction, merged gene set); a GPR with several AND nodes holding the pair counts once
        if r["kind"] == "merge" and r["applied"]:
            k = (r["rxn"], r["genes"])
            if k in seen:
                seen[k]["n_and_nodes"] += 1
                continue
            seen[k] = {**r, "n_and_nodes": 1}
            merges.append(seen[k])
    for i, m in enumerate(merges, 1):
        genes = m["genes"].split(";")
        ids = []
        for a, b in combinations(genes, 2):
            ids += [f"ENSEMBL116_PARALOG:{pair_key(a, b)}", f"NO_SHARED_COMPLEX:{pair_key(a, b)}"]
        rows.append({"edit_id": f"F1_H1_{i:03d}", "kind": "AND2OR", "rxn": m["rxn"], "rxn_name": m["rxn_name"],
                     "subsystem": m["subsystem"], "genes": m["genes"],
                     "symbols": ";".join(str(sym.get(g, "")) for g in genes),
                     "compartments": " | ".join(";".join(sorted(fc.comp_set(g1.compartments.get(g)))) for g in genes),
                     "pairs": m["pairs"], "n_and_nodes": m["n_and_nodes"], "old_node": m["old_node"],
                     "evidence_ids": ";".join(ids),
                     "rule": "D17: Ensembl paralog min %id >= 30 + same genes.tsv compartment set + no shared complex; "
                             "reaction not exchange/transport/pseudo-lumped/AND>=6"})
    t = model.reactions.get_by_id(fc.TYMS_REACTION)
    rows.append({"edit_id": "F1_H1_TYMS", "kind": "TYMS_IRREVERSIBLE", "rxn": t.id, "rxn_name": t.name,
                 "subsystem": t.subsystem, "genes": ";".join(sorted(g.id for g in t.genes)),
                 "symbols": ";".join(sorted(str(sym.get(g.id, "")) for g in t.genes)), "compartments": "",
                 "pairs": "", "n_and_nodes": 0, "old_node": f"bounds {t.bounds}", "evidence_ids": "DESIGN §3.1; Human2 v2.1.1 MAR04644 bounds (0, 1000)",
                 "rule": f"bounds {t.bounds} -> {fc.TYMS_BOUNDS}"})
    return pd.DataFrame(rows)


def raw_check(edits: pd.DataFrame, genes_union: pd.DataFrame) -> pd.DataFrame:
    """Independent re-check from the raw files: W5 validator (paralog row, no shared complex) plus min %id >= 30 from
    the raw Ensembl rows and the compartment sets from the raw Human1 genes.tsv (csv module, not pandas/w5)."""
    from nextgem import paths
    from nextgem.w5_evidence import PARALOGS, validate

    ae = edits[edits.kind == "AND2OR"]
    v = validate(ae[["edit_id", "evidence_ids"]], genes_union)
    v["check_source"] = "w5_evidence.validate"
    need = {}
    for _, e in ae.iterrows():
        gs = e.genes.split(";")
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                need.setdefault(frozenset((gs[i], gs[j])), []).append(e.edit_id)
    best = {}
    with open(paths.data_root() / PARALOGS) as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        ia, ib = hdr.index("Gene stable ID"), hdr.index("Human paralogue gene stable ID")
        p1 = hdr.index("Paralogue %id. target Human gene identical to query gene")
        p2 = hdr.index("Paralogue %id. query gene identical to target Human gene")
        for row in rd:
            k = frozenset((row[ia], row[ib]))
            if k in need:
                best[k] = max(best.get(k, -1.0), min(float(row[p1]), float(row[p2])))
    comp = {}
    with open(fc.PROJECT / "data" / "Human1" / "model" / "genes.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            comp[row["genes"]] = frozenset(t.strip() for t in row["compartments"].split(";") if t.strip())
    rows = []
    for k, eids in need.items():
        a, b = sorted(k)
        for eid in eids:
            rows.append({"edit_id": eid, "evidence_id": f"PID>={fc.MIN_PID:g}:{a}|{b}", "ok": best.get(k, -1) >= fc.MIN_PID,
                         "check": f"raw Ensembl min %id = {best.get(k, float('nan')):.2f}", "check_source": "raw csv"})
            rows.append({"edit_id": eid, "evidence_id": f"SAME_COMPARTMENTS:{a}|{b}",
                         "ok": bool(comp.get(a)) and comp.get(a) == comp.get(b),
                         "check": f"{';'.join(sorted(comp.get(a, [])))} vs {';'.join(sorted(comp.get(b, [])))}",
                         "check_source": "raw csv"})
    return pd.concat([v, pd.DataFrame(rows)], ignore_index=True)


def closure_changes(model, rules: dict) -> pd.DataFrame:
    from nextgem.w2_engine import gene_closures

    base = gene_closures(model)
    saved = {rid: model.reactions.get_by_id(rid).gene_reaction_rule for rid in rules}
    try:
        fc.apply_arm(model, fc.ARM["GPR"], rules)
        new = gene_closures(model)
    finally:
        for rid, rule in saved.items():
            model.reactions.get_by_id(rid).gene_reaction_rule = rule
    rows = []
    for g in sorted(base):
        a, b = set(base[g]), set(new[g])
        if a != b:
            rows.append({"gene": g, "n_closed_before": len(a), "n_closed_after": len(b),
                         "no_longer_closed": ";".join(sorted(a - b)), "newly_closed": ";".join(sorted(b - a))})
    return pd.DataFrame(rows, columns=["gene", "n_closed_before", "n_closed_after", "no_longer_closed", "newly_closed"])


def summarise(recs, excl, lumped, curated, edits, g1, closures) -> dict:
    sym = g1.symbol.to_dict()
    P = pd.DataFrame([r for r in recs if r["kind"] == "pair"])
    P["polr"] = [str(sym.get(a, "")).startswith("POLR") and str(sym.get(b, "")).startswith("POLR")
                 for a, b in zip(P.a, P.b)]
    P["excluded"] = P.rxn_excluded.fillna("") != ""
    P["w5_literal"] = P.complex_ok  # W5a AND -> OR: any Ensembl paralog + no shared complex
    P["w5_pid30"] = P.complex_ok & P.pid_ok

    def by_reason(mask):
        sub = P[mask & P.excluded]
        out = {}
        for tag in ("E1", "E2", "E3", "E4"):
            m = sub.rxn_excluded.str.contains(tag)
            out[tag] = {"pair_rows": int(m.sum()), "POLR_pair_rows": int((m & sub.polr).sum()),
                        "distinct_gene_pairs": int(sub[m][["a", "b"]].drop_duplicates().shape[0]),
                        "reactions": int(sub[m].rxn.nunique())}
        return out

    def count(mask):
        s = P[mask]
        return {"pair_rows": int(len(s)), "distinct_gene_pairs": int(s[["a", "b"]].drop_duplicates().shape[0]),
                "reactions": int(s.rxn.nunique()), "POLR_pair_rows": int(s.polr.sum()),
                "POLR_distinct_gene_pairs": int(s[s.polr][["a", "b"]].drop_duplicates().shape[0])}

    merges = [r for r in recs if r["kind"] == "merge"]
    blocked = [r for r in recs if r["kind"] == "blocked"]
    ae = edits[edits.kind == "AND2OR"]
    return {
        "parent": "human1 (v1.19.0)",
        "rule": {"min_pid": fc.MIN_PID, "same_compartment": "equal non-empty genes.tsv compartment sets",
                 "no_shared_complex": True, "max_and_genes_excluded_at": fc.MAX_AND_GENES,
                 "lumped": "subsystem in Artificial/Pool reactions, or >= 2 non-integer coefficients on one side summing to 1",
                 "edit_types": "AND2OR only (+ TYMS bounds)", "edit_unit": "single reaction"},
        "reactions_with_gpr": None,  # filled by run()
        "lumped_reactions_total": int(len(lumped)), "lumped_reactions_with_gpr": int((lumped.n_genes > 0).sum()),
        "lumped_reactions_with_gpr_ids": sorted(lumped[lumped.n_genes > 0].rxn),
        "excluded_gpr_reactions": {t: int(excl.exclusions.str.contains(t).sum()) for t in ("E1", "E2", "E3", "E4")},
        "excluded_gpr_reactions_total": int(len(excl)),
        "and_sibling_paralog_pairs_all": count(P.paralog_type.notna()),
        "w5_literal_rule_pairs": count(P.w5_literal),
        "w5_literal_rule_pairs_excluded_by_reaction": by_reason(P.w5_literal),
        "w5_pid30_rule_pairs": count(P.w5_pid30),
        "w5_pid30_rule_pairs_excluded_by_reaction": by_reason(P.w5_pid30),
        "pair_criteria_qualified": count(P.pair_qualifies),
        "pair_qualified_but_excluded_by_reaction": by_reason(P.pair_qualifies),
        "pair_fail_counts_among_w5_pid30": {
            "compartment": int((P.w5_pid30 & ~P.comp_ok).sum())},
        "merges_applied": {"edits_rxn_x_geneset": int(len(ae)), "and_node_merges": int(ae.n_and_nodes.sum()), "reactions": int(ae.rxn.nunique()),
                           "distinct_gene_sets": int(ae.genes.nunique()),
                           "gene_sets": sorted(set(ae.symbols))},
        "merges_qualified_not_applied_reaction_excluded": int(sum(1 for m in merges if not m["applied"])),
        "blocked_non_clique_components": int(len(blocked)),
        "reactions_with_changed_gpr": int(len(curated)),
        "genes_with_changed_ko_closure": int(len(closures)),
        "genes_with_changed_ko_closure_symbols": sorted(str(sym.get(g, g)) for g in closures.gene),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-dirty", action="store_true")
    a = ap.parse_args(argv)
    from nextgem import w5_evidence as we
    from nextgem.model import GEMS

    inputs = [we.data(we.PARALOGS), we.data(we.CORUM), we.data(we.CPX), GEMS["human1"].genes_tsv,
              GEMS["human2"].genes_tsv, GEMS["human1"].sbml]
    with fc.provenance(inputs=inputs, solver="none", workers=1, allow_dirty=a.allow_dirty) as prov:
        outs, summ = run()
        prov.output(*outs)
        prov.note(**{k: summ[k] for k in ("lumped_reactions_with_gpr", "excluded_gpr_reactions_total",
                                           "reactions_with_changed_gpr", "genes_with_changed_ko_closure")},
                  merges_applied=summ["merges_applied"]["edits_rxn_x_geneset"])
    return 0


def run():
    import cobra

    from nextgem import w5_evidence as we
    from nextgem.model import load_model

    t0 = time.time()
    cobra.Configuration().solver = "glpk"
    model = load_model("human1")
    if tuple(model.reactions.get_by_id(fc.TYMS_REACTION).bounds) == fc.TYMS_BOUNDS:
        raise RuntimeError("released H1 MAR04644 is already irreversible")
    g1, gu = gene_tables()
    paralog = we._read_paralogs({g.id for g in model.genes})
    complexes, _ = we._read_complexes(gu)
    recs, excl, lumped, curated = enumerate_edits(model, g1, paralog, complexes)
    edits = build_edit_table(recs, g1, model)
    chk = raw_check(edits, gu)
    rules = dict(zip(curated.rxn, curated.new_rule)) if len(curated) else {}
    clos = closure_changes(model, rules)
    pairs = pd.DataFrame([r for r in recs if r["kind"] in ("pair", "blocked")])
    sym = g1.symbol.to_dict()
    if len(pairs):
        pairs.insert(pairs.columns.get_loc("a") + 2, "symbols", [f"{sym.get(x, '')}|{sym.get(y, '')}" if isinstance(x, str) else ""
                                                             for x, y in zip(pairs.a, pairs.b)])
    clos.insert(1, "symbol", clos.gene.map(sym))
    summ = summarise(recs, excl, lumped, curated, edits, g1, clos)
    summ["reactions_with_gpr"] = int(sum(1 for r in model.reactions if r.gene_reaction_rule.strip()))
    summ["evidence_check_all_ok"] = bool(chk.ok.all())
    summ["evidence_check_n"] = int(len(chk))
    summ["wall_s"] = round(time.time() - t0, 1)
    D = fc.F1
    outs = [fc.atomic_write_tsv(edits, D / "edits_human1.tsv"),
            fc.atomic_write_tsv(curated[["rxn", "old_rule", "new_rule", "old_canonical"]] if len(curated) else
                                pd.DataFrame(columns=["rxn", "old_rule", "new_rule", "old_canonical"]), fc.curated_rules_path()),
            fc.atomic_write_tsv(pairs, D / "paralog_pairs_human1.tsv"),
            fc.atomic_write_tsv(excl, D / "excluded_reactions_human1.tsv"),
            fc.atomic_write_tsv(lumped, D / "lumped_reactions_human1.tsv"),
            fc.atomic_write_tsv(chk, D / "evidence_check_human1.tsv"),
            fc.atomic_write_tsv(clos, D / "closure_changes_human1.tsv"),
            fc.atomic_write_text(json.dumps(summ, indent=1, default=float) + "\n", D / "edits_summary.json")]
    print(json.dumps(summ, indent=1, default=float))
    if not summ["evidence_check_all_ok"]:
        raise RuntimeError("evidence check failed")
    return outs, summ


if __name__ == "__main__":
    raise SystemExit(main())
