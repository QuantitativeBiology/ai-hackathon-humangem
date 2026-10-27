"""W5 independent review, part 1 (LABEL-FREE). Reviewer's own code; does not import nextgem.w5 / w5_evidence.

(a) Independent evidence check of every W5a edit against the raw local files (Ensembl 116 paralogs, CORUM, Complex
    Portal, ProCan, Nusinow, DepMap 26Q3 expression read from the raw CSV with IsDefaultEntryForModel == "Yes").
(b) Own application of the W5a edit list to the released GPR strings (cobra GPR AST), compared, gene by gene, with the
    author's curated rules (single-gene KO closure sets).
(c) KO ratios (cobra Gene.knock_out semantics, GLPK, M0 = Ham's, biomass MAR13082) for every gene of an edited GPR:
    generic, W5a (H1 on the TYMS base), W5a without the RNA-polymerase merges (attribution diagnostic), each single
    edit alone (flip counts), plus a full single-gene sweep of H1 + TYMS (MAR04644 irreversible) for TYMS-only.
(d) CAD / dihydroorotate check (biomass pool MAR10064) under Ham's +/- uridine and +/- orotate.

Usage: OMP_NUM_THREADS=1 .venv/bin/python results/w5/review/r1_labelfree.py
Writes results/w5/review/r1_*.{tsv,json,parquet}.
"""
from __future__ import annotations

import ast
import gzip
import json
import math
import re
import sys
import time
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from nextgem.knockouts import configure  # noqa: E402  (shared infra: Ham's medium + biomass objective)
from nextgem.model import load_model, GEMS  # noqa: E402

ROOT = Path("/data/tmpA/nextgem")
W5 = ROOT / "results" / "w5"
OUT = W5 / "review"
D = ROOT / "data"
T0, C0 = time.time(), time.process_time()


# ----------------------------------------------------------------------------- GPR trees from cobra's AST
def tree_of(rule: str):
    """cobra GPR -> nested ('and'|'or', [children]) / gene string (no canonicalisation beyond flattening)."""
    if not rule.strip():
        return None
    g = cobra.core.gene.GPR.from_string(rule)

    def rec(n):
        if isinstance(n, ast.Expression):
            return rec(n.body)
        if isinstance(n, ast.Name):
            return n.id
        if isinstance(n, ast.BoolOp):
            op = "and" if isinstance(n.op, ast.And) else "or"
            ch = []
            for v in n.values:
                c = rec(v)
                if not isinstance(c, str) and c[0] == op:
                    ch.extend(c[1])
                else:
                    ch.append(c)
            return (op, ch)
        raise ValueError(type(n))

    return rec(g.body)


def ev(t, off):
    if isinstance(t, str):
        return t not in off
    if t[0] == "and":
        return all(ev(c, off) for c in t[1])
    return any(ev(c, off) for c in t[1])


def genes(t):
    if t is None:
        return set()
    if isinstance(t, str):
        return {t}
    return set().union(*(genes(c) for c in t[1]))


def leaves(node):
    return {c for c in node[1] if isinstance(c, str)}


def apply_and2or(t, gs: set):
    """Every AND node whose direct leaf children include all of ``gs``: replace those leaves by (OR of gs)."""
    n_hit = 0

    def rec(x):
        nonlocal n_hit
        if isinstance(x, str):
            return x
        ch = [rec(c) for c in x[1]]
        if x[0] == "and" and gs <= {c for c in ch if isinstance(c, str)}:
            n_hit += 1
            ch = [c for c in ch if not (isinstance(c, str) and c in gs)] + [("or", sorted(gs))]
        return (x[0], ch)

    return rec(t), n_hit


def apply_rmbranch(t, g: str):
    """Every OR node with leaf child ``g`` (and >= 2 children): drop that leaf."""
    n_hit = 0

    def rec(x):
        nonlocal n_hit
        if isinstance(x, str):
            return x
        ch = [rec(c) for c in x[1]]
        if x[0] == "or" and g in ch and len(ch) >= 2:
            n_hit += 1
            ch = [c for c in ch if c != g]
            if len(ch) == 1:
                return ch[0]
        return (x[0], ch)

    return rec(t), n_hit


def to_str(t):
    if t is None:
        return ""
    if isinstance(t, str):
        return t
    return "(" + f" {t[0]} ".join(to_str(c) for c in t[1]) + ")"


def closes(t, g):
    return t is not None and g in genes(t) and not ev(t, {g})


# ----------------------------------------------------------------------------- KO engine (cobra semantics, GLPK)
class KO:
    def __init__(self, model):
        self.m = model
        sol = model.optimize()
        self.wt = float(sol.objective_value)
        self.flux = sol.fluxes.to_dict()
        assert self.wt > 1e-6
        self.cache = {}
        self.n_lp = 0

    def ratio(self, closed):
        closed = frozenset(closed)
        if not closed:
            return 1.0
        if closed in self.cache:
            return self.cache[closed]
        if all(abs(self.flux[r]) <= 1e-9 for r in closed):  # WT optimum stays feasible: exact shortcut
            return 1.0
        with self.m:
            for r in closed:
                self.m.reactions.get_by_id(r).bounds = (0.0, 0.0)
            v = self.m.slim_optimize(error_value=float("nan"))
            st = self.m.solver.status
        self.n_lp += 1
        if st == "infeasible":
            v = 0.0
        elif st != "optimal":
            raise RuntimeError(f"status {st}")
        r = float(min(1.0, max(0.0, round(v / self.wt, 6))) + 0.0)
        self.cache[closed] = r
        return r


def closure_map(rules: dict):
    """reaction -> tree  ==>  gene -> frozenset(closed reactions)."""
    out = {}
    for rid, t in rules.items():
        for g in genes(t):
            if closes(t, g):
                out.setdefault(g, set()).add(rid)
    return {g: frozenset(v) for g, v in out.items()}


def load(parent, tyms=False):
    m = load_model(parent)
    configure(m)
    if tyms:
        m.reactions.get_by_id("MAR04644").bounds = (0.0, 1000.0)
    return m


# ----------------------------------------------------------------------------- evidence (raw files)
def read_paralog_pairs(need: set):
    cols = ["Gene stable ID", "Human paralogue gene stable ID", "Paralogue %id. target Human gene identical to query gene",
            "Paralogue %id. query gene identical to target Human gene", "Human paralogue homology type",
            "Paralogue last common ancestor with Human"]
    rows = []
    for ch in pd.read_csv(D / "paralogs/ensembl_human_paralogs.tsv", sep="\t", usecols=cols, dtype=str, chunksize=1_000_000):
        k = ch[cols[0]] + "|" + ch[cols[1]]
        rows.append(ch[k.isin(need)])
    df = pd.concat(rows)
    df["key"] = [ "|".join(sorted((a, b))) for a, b in zip(df[cols[0]], df[cols[1]])]
    df["pid_min"] = np.minimum(df[cols[2]].astype(float), df[cols[3]].astype(float))
    return df.groupby("key").agg(type=(cols[4], "first"), pid_min=("pid_min", "max"), lca=(cols[5], "first"))


def complexes_by_gene(sym: dict, acc: dict):
    """gene -> set of complex ids from CORUM (all organisms in the human file) and Complex Portal (raw parse)."""
    out = {}
    c = pd.read_csv(D / "complexes/corum_humanComplexes.txt", sep="\t", dtype=str).fillna("")
    name_col = "complex_name" if "complex_name" in c.columns else c.columns[1]
    cname = {}
    for _, r in c.iterrows():
        us = {x.strip() for x in r["subunits_uniprot_id"].split(";")}
        gs = {x.strip() for x in r["subunits_gene_name"].split(";")}
        for g in sym:
            if sym[g] in gs or (acc.get(g, set()) & us):
                out.setdefault(g, set()).add("CORUM:" + r["complex_id"])
        cname["CORUM:" + r["complex_id"]] = r[name_col]
    p = pd.read_csv(D / "complexes/complexportal_9606.tsv", sep="\t", dtype=str).fillna("")
    for _, r in p.iterrows():
        toks = {re.split(r"[\(\-]", t)[0] for t in r["Expanded participant list"].split("|")}
        for g in sym:
            if acc.get(g, set()) & toks:
                out.setdefault(g, set()).add("CPX:" + r["#Complex ac"])
        cname["CPX:" + r["#Complex ac"]] = r["Recommended name"]
    return out, cname


def proteomics_detection(gs: dict, acc: dict):
    """gene -> dict(procan_acc, procan_sym, nusinow_acc, nusinow_sym): detected (any non-missing value)."""
    f = D / "proteomics_goncalves/Protein_matrix_averaged_20250211.tsv"
    two = pd.read_csv(f, sep="\t", nrows=2, header=None)
    accs, syms = list(two.iloc[0]), list(two.iloc[1])
    want_idx = [i for i, (a, s) in enumerate(zip(accs, syms)) if i >= 2 and (
        any(str(a).split("-")[0] in acc.get(g, set()) for g in gs) or s in set(gs.values()))]
    pc = pd.read_csv(f, sep="\t", header=None, skiprows=3, usecols=want_idx)
    det_pc_acc, det_pc_sym = {}, {}
    for i in want_idx:
        nonmiss = bool(pc[i].notna().any())
        a, s = str(accs[i]).split("-")[0], syms[i]
        for g in gs:
            if a in acc.get(g, set()):
                det_pc_acc[g] = det_pc_acc.get(g, False) or nonmiss
            if s == gs[g]:
                det_pc_sym[g] = det_pc_sym.get(g, False) or nonmiss
    nu = pd.read_csv(D / "proteomics_nusinow/protein_quant_current_normalized.csv.gz")
    ch = [c for c in nu.columns if "_TenPx" in c and not c.endswith("_Peptides")]
    nu["acc0"] = nu.Uniprot_Acc.astype(str).str.split("-").str[0]
    nu["nonmiss"] = nu[ch].notna().any(axis=1)
    det_nu_acc = {g: bool(nu.loc[nu.acc0.isin(acc.get(g, set())), "nonmiss"].any()) for g in gs}
    det_nu_sym = {g: bool(nu.loc[nu.Gene_Symbol == gs[g], "nonmiss"].any()) for g in gs}
    return {g: {"procan_acc": det_pc_acc.get(g, False), "procan_sym": det_pc_sym.get(g, False),
                "nusinow_acc": det_nu_acc[g], "nusinow_sym": det_nu_sym[g]} for g in gs}


def expression_fraction(entrez: dict, screened: set):
    f = D / "depmap/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    hdr = pd.read_csv(f, nrows=0).columns
    col = {}
    for c in hdr:
        m = re.match(r"^(.+) \((\d+)\)$", c)
        if m and int(m.group(2)) in set(entrez.values()):
            col[int(m.group(2))] = c
    df = pd.read_csv(f, usecols=["ModelID", "IsDefaultEntryForModel"] + list(col.values()))
    df = df[df.IsDefaultEntryForModel == "Yes"].drop_duplicates("ModelID").set_index("ModelID")
    out = {}
    for g, e in entrez.items():
        c = col.get(e)
        if c is None:
            out[g] = (np.nan, 0, np.nan, 0)
            continue
        v = df[c].dropna()
        vs = v[v.index.isin(screened)]
        out[g] = (float((v < 1.0).mean()), int(len(v)), float((vs < 1.0).mean()), int(len(vs)))
    return out


# ----------------------------------------------------------------------------- main
def main():
    from nextgem import labels  # ids() only: ModelIDs and gene header, no values

    screened = set(labels.ids()["model_ids"])
    report = {}
    gene_rows, edit_rows, cmp_rows, alone_rows = [], [], [], []
    ko_tables = {}
    for parent in ("human2", "human1"):
        S = "H2" if parent == "human2" else "H1"
        gt = pd.read_csv(GEMS[parent].genes_tsv, sep="\t", dtype=str)
        sym = dict(zip(gt.genes, gt.geneSymbols))
        ent = {g: int(e.split(";")[0]) for g, e in zip(gt.genes, gt.geneEntrezID) if isinstance(e, str) and e.split(";")[0].isdigit()}
        acc = {g: set(u.split(";")) if isinstance(u, str) else set() for g, u in zip(gt.genes, gt.geneUniProtID)}
        edits = pd.read_csv(W5 / f"edits_W5a_{parent}.tsv", sep="\t")
        cur = pd.read_csv(W5 / f"curated_W5a_{parent}.tsv", sep="\t")
        m = load(parent, tyms=(parent == "human1"))
        rules = {r.id: tree_of(r.gene_reaction_rule) for r in m.reactions if r.gene_reaction_rule.strip()}
        # ---- (b) apply edits myself, group by group (an edit lists all reactions of its GPR group)
        new = {}
        seen = set()
        for e in edits.itertuples():
            gs = set(e.genes.split(";"))
            key = (e.kind, e.genes, e.reactions)
            if key in seen:  # same edit at a second node of the same tree (RRM): applied to all matching nodes above
                continue
            seen.add(key)
            for rx in e.reactions.split(";"):
                t = new.get(rx, rules[rx])
                if e.kind == "AND2OR":
                    t2, hit = apply_and2or(t, gs)
                elif e.kind == "RMBRANCH":
                    (g,) = gs
                    t2, hit = apply_rmbranch(t, g)
                else:
                    raise ValueError(e.kind)
                if hit == 0:
                    raise AssertionError(f"{e.edit_id} {rx}: no node matched")
                new[rx] = t2
        # RRM-style duplicate edits: an AND2OR edit listed twice for two nodes is applied to all matching nodes at
        # the first listing; the second is then a no-op on the already-merged tree -> handled since hit counts >0 above
        author = {rx: tree_of(r.curated_rule) for r in cur[cur.edited].itertuples() for rx in r.reactions.split(";")}
        assert set(author) == set(new), (set(author) ^ set(new))
        n_mis = 0
        for rx in new:
            for g in genes(rules[rx]) | genes(new[rx]) | genes(author[rx]):
                a, b = closes(new[rx], g), closes(author[rx], g)
                if a != b:
                    n_mis += 1
                    cmp_rows.append({"parent": parent, "reaction": rx, "gene": g, "mine": a, "author": b})
        report[f"{S}_edited_reactions"] = len(new)
        report[f"{S}_closure_mismatches_vs_author"] = n_mis
        # ---- (c) KO ratios for genes of edited GPRs
        eng = KO(m)
        base_cl = closure_map(rules)
        cur_rules = {**rules, **new}
        cur_cl = closure_map(cur_rules)
        polr = edits[(edits.kind == "AND2OR") & edits.symbols.str.contains("POLR")].edit_id.tolist()
        # W5a without the RNA-polymerase merges (attribution diagnostic)
        new_np = {}
        seen = set()
        for e in edits[~edits.edit_id.isin(polr)].itertuples():
            gs = set(e.genes.split(";"))
            if (e.kind, e.genes, e.reactions) in seen:
                continue
            seen.add((e.kind, e.genes, e.reactions))
            for rx in e.reactions.split(";"):
                t = new_np.get(rx, rules[rx])
                t2, hit = apply_and2or(t, gs) if e.kind == "AND2OR" else apply_rmbranch(t, next(iter(gs)))
                new_np[rx] = t2
        np_cl = closure_map({**rules, **new_np})
        touched = sorted(set().union(*(genes(rules[rx]) for rx in new)))
        for g in touched:
            gene_rows.append({"parent": parent, "gene": g, "symbol": sym.get(g), "entrez": ent.get(g),
                              "ratio_generic_glpk": eng.ratio(base_cl.get(g, ())),
                              "ratio_W5a": eng.ratio(cur_cl.get(g, ())),
                              "ratio_W5a_noPOLR": eng.ratio(np_cl.get(g, ()))})
        # each edit alone on the parent (flip = change of the ratio_r < 0.01 call)
        for e in edits.itertuples():
            gs = set(e.genes.split(";"))
            one = {}
            for rx in e.reactions.split(";"):
                one[rx] = (apply_and2or(rules[rx], gs) if e.kind == "AND2OR" else apply_rmbranch(rules[rx], next(iter(gs))))[0]
            cl1 = closure_map({**rules, **one})
            tg = sorted(set().union(*(genes(rules[rx]) for rx in one)))
            flips, changed = [], []
            for g in tg:
                r0, r1 = eng.ratio(base_cl.get(g, ())), eng.ratio(cl1.get(g, ()))
                if r0 != r1:
                    changed.append(sym.get(g, g))
                if (r0 < 0.01) != (r1 < 0.01):
                    flips.append(sym.get(g, g))
            alone_rows.append({"parent": parent, "edit_id": e.edit_id, "kind": e.kind, "symbols": e.symbols,
                               "n_flips_alone": len(flips), "flipped_alone": ";".join(flips),
                               "n_ratio_changed_alone": len(changed)})
        # group-level flips (all W5a edits of the group together)
        for grp, sub in edits.groupby("group"):
            rx = set(sub.reactions.iloc[0].split(";"))
            tg = sorted(set().union(*(genes(rules[r]) for r in rx)))
            fl = [sym.get(g, g) for g in tg if (eng.ratio(base_cl.get(g, ())) < 0.01) != (eng.ratio(cur_cl.get(g, ())) < 0.01)]
            for eid in sub.edit_id:
                alone_rows.append({"parent": parent, "edit_id": eid, "kind": "GROUP", "symbols": "",
                                   "n_flips_alone": len(fl), "flipped_alone": ";".join(fl), "n_ratio_changed_alone": np.nan})
        # ---- (a) evidence
        need_pairs, gset = set(), {}
        for e in edits.itertuples():
            gs = e.genes.split(";")
            for g in gs:
                gset[g] = sym.get(g)
            if e.kind == "AND2OR":
                for i in range(len(gs)):
                    for j in range(i + 1, len(gs)):
                        need_pairs |= {f"{gs[i]}|{gs[j]}", f"{gs[j]}|{gs[i]}"}
        report.setdefault("_need", {})[parent] = (need_pairs, gset, ent, acc, edits)
        if parent == "human1":
            # full sweep of H1 + TYMS (TYMS-only arm) with the same engine
            full = {}
            for g in sorted(set(base_cl) | {x.id for x in m.genes}):
                full[g] = eng.ratio(base_cl.get(g, ()))
            ko_tables["H1_TYMS_full"] = pd.DataFrame({"gene": list(full), "ratio_tyms": list(full.values())})
        report[f"{S}_n_lp"] = eng.n_lp
        report[f"{S}_wt"] = eng.wt
        print(parent, "done", round(time.time() - T0), "s", flush=True)

    # ---- evidence checks (union over parents)
    allpairs, allg, allent, allacc = set(), {}, {}, {}
    for parent, (np_, gs, ent, acc, _) in report["_need"].items():
        allpairs |= np_
        allg.update(gs)
        allent.update({g: ent[g] for g in gs if g in ent})
        allacc.update({g: acc.get(g, set()) for g in gs})
    par = read_paralog_pairs(allpairs)
    cx, cname = complexes_by_gene(allg, allacc)
    det = proteomics_detection(allg, allacc)
    ex = expression_fraction(allent, screened)
    for parent, (np_, gs, ent, acc, edits) in report["_need"].items():
        for e in edits.itertuples():
            gl = e.genes.split(";")
            row = {"parent": parent, "edit_id": e.edit_id, "kind": e.kind, "symbols": e.symbols}
            if e.kind == "AND2OR":
                prs = []
                ok = True
                shared_any = []
                for i in range(len(gl)):
                    for j in range(i + 1, len(gl)):
                        k = "|".join(sorted((gl[i], gl[j])))
                        p = par.loc[k] if k in par.index else None
                        sh = sorted(cx.get(gl[i], set()) & cx.get(gl[j], set()))
                        prs.append(f"{allg[gl[i]]}|{allg[gl[j]]}:{'none' if p is None else p['type']}"
                                   f"{'' if p is None else ',pid=%.1f,lca=%s' % (p['pid_min'], p['lca'])}"
                                   f",shared={len(sh)}")
                        ok &= (p is not None) and not sh
                        shared_any += sh
                cxs = sorted(set().union(*(cx.get(g, set()) for g in gl)))
                row.update({"evidence_ok": ok, "pairs": ";".join(prs),
                            "complexes_of_genes": ";".join(f"{c}={cname.get(c, '')[:40]}" for c in cxs[:12])})
            else:
                g = gl[0]
                d = det[g]
                fr, n, frs, ns = ex.get(g, (np.nan, 0, np.nan, 0))
                ok = (not any(d.values())) and frs >= 0.95
                row.update({"evidence_ok": ok, **d, "frac_tpm_lt1_all_default": fr, "n_all_default": n,
                            "frac_tpm_lt1_screened": frs, "n_screened_with_expr": ns})
            edit_rows.append(row)
    report.pop("_need")
    ev_df = pd.DataFrame(edit_rows)
    al = pd.DataFrame(alone_rows)
    gk = pd.DataFrame(gene_rows)
    ev_df.to_csv(OUT / "r1_evidence_check.tsv", sep="\t", index=False)
    al.to_csv(OUT / "r1_edit_flips.tsv", sep="\t", index=False)
    gk.to_csv(OUT / "r1_ko_affected_genes.tsv", sep="\t", index=False)
    pd.DataFrame(cmp_rows).to_csv(OUT / "r1_closure_mismatches.tsv", sep="\t", index=False)
    ko_tables["H1_TYMS_full"].to_parquet(OUT / "r1_h1_tyms_full.parquet", index=False)
    report["n_edits_evidence_ok"] = int(ev_df.evidence_ok.sum())
    report["n_edits"] = int(len(ev_df))
    report["wall_s"] = round(time.time() - T0, 1)
    report["cpu_s"] = round(time.process_time() - C0, 1)
    (OUT / "r1_summary.json").write_text(json.dumps(report, indent=1, default=float))
    print(json.dumps(report, indent=1, default=float))


if __name__ == "__main__":
    sys.exit(main())
