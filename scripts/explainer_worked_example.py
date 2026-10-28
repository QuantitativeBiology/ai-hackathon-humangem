"""Explainer worked example: trace one real training cell line through every step of P1 (label-free model outputs).

P1 = Human1 (Human-GEM v1.19.0) + removal of the cofactor/vitamin pool from the biomass reaction
(``B-loo:cofactor_pool_biomass``) + the line's documented culture medium with serum (M2, f = 1).

This script re-runs the knockout (KO) sweeps for one cell line with the project's own engine
(``nextgem.w3.ko_sweep``, GLPK) at each step and checks them against the committed P1 scores
(results/m6/scores/sweep_e.parquet). It reads NO label file:
  * per-line tAP of the chosen line comes from the committed results/m6/train/per_line.tsv (and G* from
    results/w1/train_per_line.tsv), which are aggregates over genes;
  * per-gene essentiality context comes from the committed, unsealed test summary
    results/m6/test/unsealed/T1/per_gene.tsv (fraction of the 380 test lines in which the gene is essential),
    which is an aggregate over lines of *other* (test) lineages.
No per-line x per-gene label value is read or written.

Steps traced (cells of the M6 factorial):
  0  generic Human1 on Ham's medium (the reference G*)          key T0|<Ham's signature>
  1  + cofactor-pool removal (still Ham's)                       key B-loo|<Ham's signature>
  (side) medium only: generic Human1 on the line's medium       key T0|<line M2 signature>
  2  + documented medium and serum = the P1 model                key B-loo|<line M2 signature>
  3  P1 scores read as a ranking with ties (tAP)

Usage:
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/explainer_worked_example.py \
      --line ACH-000463 --out <dir>
"""

from __future__ import annotations

import os

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

R = PROJECT / "results"
CALL = 0.99  # e >= 0.99 <=> KO/WT ratio < 0.01: "called lethal" (DESIGN §6.1 fixed cutoff)
BLOO = "B-loo:cofactor_pool_biomass"
POOL_MET = "MAM10012c"  # cofactor_pool_biomass, reactant of the biomass reaction MAR13082
POOL_RXN = "MAR10065"  # assembles the pool
BIOMASS = "MAR13082"

# Illustrative panel: symbol -> (pathway, why it is in the panel)
PANEL = {
    "CYP26A1": ("retinoic-acid catabolism", "biomass artefact: the pool demands retinoid catabolites"),
    "CYP24A1": ("vitamin-D catabolism", "biomass artefact: the pool demands vitamin-D catabolites"),
    "BBOX1": ("carnitine synthesis", "artefact twice over: pool demands carnitine; serum supplies it"),
    "COQ2": ("coenzyme Q synthesis", "real cofactor, removed together with the pool"),
    "HMBS": ("heme synthesis", "real cofactor, removed together with the pool"),
    "PPAT": ("de novo purine synthesis", "masked on Ham's by hypoxanthine"),
    "DHODH": ("de novo pyrimidine synthesis", "masked on Ham's by thymidine (Human1 TYMS is reversible)"),
    "UMPS": ("de novo pyrimidine synthesis", "masked on Ham's by thymidine (Human1 TYMS is reversible)"),
    "TYMS": ("thymidylate synthesis", "masked on Ham's by thymidine"),
    "DHCR7": ("cholesterol synthesis (last step)", "rescued by serum cholesterol"),
    "GCLC": ("glutathione synthesis", "rescued by glutathione in the documented medium/serum model"),
    "HMGCR": ("mevalonate pathway", "known cost: non-lethal once the pool is gone and serum is present"),
    "CAD": ("de novo pyrimidine synthesis (first step)", "lethal at every step (control)"),
}

STEPS = [
    ("step0", "Generic Human1 on Ham's medium (G*)", "T0", "M0"),
    ("step1", "+ cofactor pool removed from biomass (Ham's)", BLOO, "M0"),
    ("medium_only", "Side cell: generic Human1 on the line's medium + serum", "T0", "M2(f=1)"),
    ("step2", "+ documented medium and serum = P1 model", BLOO, "M2(f=1)"),
]
PER_LINE_CELL = {"step0": "H1", "step1": "H1+Bloo/M0", "medium_only": "H1/M2", "step2": "H1+Bloo/M2"}


def gene_symbols() -> pd.DataFrame:
    mp = pd.read_csv(R / "baseline" / "human1_gene_mapping.tsv", sep="\t")
    return mp


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="ACH-000463", help="DepMap ModelID of a training line")
    ap.add_argument("--out", required=True, help="output directory (JSON + TSV + Markdown table)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    # ---------------------------------------------------------------- 1. the line (metadata only)
    model_meta = pd.read_csv(PROJECT / "data" / "depmap" / "Model.csv",
                             usecols=["ModelID", "CellLineName", "OncotreeLineage", "OncotreeSubtype"])
    split = pd.read_csv(R / "splits" / "lineage_split.tsv", sep="\t", usecols=["ModelID", "tier", "lineage"])
    media_tab = pd.read_csv(R / "w2" / "media_library" / "line_media.tsv", sep="\t",
                            usecols=["ModelID", "text", "stratum", "medium_source", "serum_pct", "medium_key"])
    assign = pd.read_csv(R / "m6" / "scores" / "line_assignment.tsv", sep="\t")
    meta = model_meta[model_meta.ModelID == a.line].iloc[0]
    sp = split[split.ModelID == a.line].iloc[0]
    if sp.tier != "train":
        raise SystemExit(f"{a.line} is in tier {sp.tier}; the worked example uses a training line only")
    med = media_tab[media_tab.ModelID == a.line].iloc[0]
    la = assign[(assign.ModelID == a.line) & (assign.parent == "human1")].set_index("step")
    sweep_id = {"M0": la.loc["M0", "w2_sweep_id"], "M2(f=1)": la.loc["M2(f=1)", "w2_sweep_id"]}
    n_same_m2 = int(((assign.parent == "human1") & (assign.step == "M2(f=1)")
                     & (assign.w2_sweep_id == sweep_id["M2(f=1)"])).sum())
    line_info = {
        "ModelID": a.line, "name": meta.CellLineName, "lineage": sp.lineage, "subtype": meta.OncotreeSubtype,
        "tier": sp.tier, "documented_medium_text": med.text, "medium_source": med.medium_source,
        "assignment_stratum": med.stratum, "serum_pct": float(med.serum_pct), "medium_key": med.medium_key,
        "w2_sweep_id": sweep_id, "n_screened_lines_sharing_this_P1_model": n_same_m2,
    }

    # ---------------------------------------------------------------- 2. model, bounds, variant
    import cobra

    cobra.Configuration().solver = "glpk"
    import m6_sweeps  # noqa: E402  (scripts/; label-free)
    from nextgem import w3
    from nextgem.model import load_model

    t0 = time.time()
    bounds, binfo = m6_sweeps.bounds_by_id({"human1": {"M0": [sweep_id["M0"]], "M2(f=1)": [sweep_id["M2(f=1)"]]}})
    model = load_model("human1")
    genes = sorted(g.id for g in model.genes)
    closures = w3.gene_closures(model)
    t_load = time.time() - t0
    variant_spec = m6_sweeps.variant_json("human1", BLOO)
    variants = {"T0": w3.Variant("T0"), BLOO: w3.Variant.from_json(variant_spec)}

    # biomass and pool structure (label-free)
    bm = model.reactions.get_by_id(BIOMASS)
    pool = model.reactions.get_by_id(POOL_RXN)
    pool_coef = bm.get_coefficient(POOL_MET)
    biomass_reactants = sorted(((m.id, m.name, round(-c, 6)) for m, c in bm.metabolites.items() if c < 0),
                               key=lambda x: -x[2])
    pool_species = sorted(m.name for m, c in pool.metabolites.items() if c < 0)

    # medium contrast (label-free)
    def name_of(ex_id):
        r = model.reactions.get_by_id(ex_id)
        return list(r.metabolites)[0].name

    b0, b2 = bounds[sweep_id["M0"]], bounds[sweep_id["M2(f=1)"]]
    only_hams = sorted(name_of(x) for x in set(b0) - set(b2))
    only_line = sorted(name_of(x) for x in set(b2) - set(b0))
    limited_line = sorted((name_of(x), round(v, 4)) for x, v in b2.items() if v < 1000)
    medium_info = {
        "hams_open_exchanges": len(b0), "line_open_exchanges": len(b2),
        "open_in_hams_not_in_line": only_hams, "open_in_line_not_in_hams": only_line,
        "line_bounded_below_1000 (serum model and micronutrient floor)": limited_line,
        "line_serum": binfo[sweep_id["M2(f=1)"]]["serum"],
        "line_medium_keys_sharing_signature": binfo[sweep_id["M2(f=1)"]]["medium_keys"],
    }

    # ---------------------------------------------------------------- 3. KO sweeps per step
    mp = gene_symbols()
    one = mp[mp.entrez.notna()].copy()  # W1 convention: the primary Entrez of every mapped gene
    one["entrez"] = one.entrez.astype(int)
    ens2ent = dict(zip(one.gene, one.entrez))
    sweep_e = pd.read_parquet(R / "m6" / "scores" / "sweep_e.parquet")
    universe = sorted(sweep_e.entrez.unique())  # U: 2,866 model genes with a CRISPR column
    ent2sym = (pd.Series(dict(zip(one.entrez, one.symbol))))
    catalog = pd.read_csv(R / "m6" / "sweeps" / "catalog.tsv", sep="\t")

    steps_out, e_by_step = {}, {}
    for key, label, var, step in STEPS:
        sid = sweep_id[step]
        t1, c1 = time.time(), time.process_time()
        r = w3.ko_sweep(model, genes, closures, bounds[sid], variants[var])
        cpu = time.process_time() - c1
        e_gene = pd.Series(1.0 - r["ratio_r"], index=genes)
        e_u = pd.Series(0.0, index=universe)  # W1 convention: a U gene absent from the model scores 0
        for g, v in e_gene.items():
            ent = ens2ent.get(g)
            if ent is not None and ent in e_u.index:
                e_u[ent] = round(float(v), 6) + 0.0
        sk = f"{var}|{sid}"
        committed = sweep_e[sweep_e.key == sk].set_index("entrez").e.reindex(universe)
        max_dev = float(np.abs(committed.to_numpy() - e_u.to_numpy()).max())
        cat = catalog[catalog.key == sk]
        e_by_step[key] = e_u
        lethal = e_u >= CALL
        partial = (e_u > 0) & (e_u < CALL)
        steps_out[key] = {
            "label": label, "variant": var, "medium_step": step, "sweep_key": sk,
            "wt_growth_au": round(float(r["wt_growth"]), 4),
            "wt_growth_committed_au": float(cat.wt_growth.iloc[0]) if len(cat) else None,
            "n_genes_scored_U": len(universe), "n_lethal": int(lethal.sum()), "n_partial": int(partial.sum()),
            "n_no_effect": int((e_u == 0).sum()), "n_distinct_scores": int(e_u.nunique()),
            "n_lp_solved": int(r["n_solved"]), "n_genes_needing_lp": int(r["n_genes_lp"]),
            "cpu_s": round(cpu, 1), "wall_s": round(time.time() - t1, 1),
            "max_abs_dev_vs_committed_sweep_e": max_dev,
        }
        print(f"{key}: WT {r['wt_growth']:.4f}, lethal {int(lethal.sum())}, partial {int(partial.sum())}, "
              f"LPs {r['n_solved']}, {cpu:.0f} CPU-s, max dev vs committed {max_dev:.2e}", flush=True)

    # ---------------------------------------------------------------- 4. what changed between steps
    def sym(ents):
        return sorted(str(ent2sym.get(x, x)) for x in ents)

    def flips(a_key, b_key):
        ea, eb = e_by_step[a_key] >= CALL, e_by_step[b_key] >= CALL
        return {"uncalled (lethal -> not)": sym(ea.index[ea & ~eb]), "newly called (not -> lethal)": sym(ea.index[~ea & eb])}

    transitions = {
        "step0->step1 (pool removed)": flips("step0", "step1"),
        "step1->step2 (medium + serum added)": flips("step1", "step2"),
        "step0->medium_only (medium alone)": flips("step0", "medium_only"),
        "step0->step2 (all of P1)": flips("step0", "step2"),
    }

    # ---------------------------------------------------------------- 5. panel genes
    pg = pd.read_csv(R / "m6" / "test" / "unsealed" / "T1" / "per_gene.tsv", sep="\t")
    pg = pg[pg.method == "P1[H1]"].set_index("entrez")
    sym2ens = {s: g for g, s in zip(one.gene, one.symbol)}
    sym2ent = {s: e for e, s in zip(one.entrez, one.symbol)}
    p1 = e_by_step["step2"]
    p1_levels = sorted(p1.unique(), reverse=True)
    rows = []
    for s, (pathway, why) in PANEL.items():
        if s not in sym2ent:
            raise SystemExit(f"{s} is not a mapped Human1 gene")
        ent = sym2ent[s]
        in_model = sym2ens[s] in set(genes)
        row = {"symbol": s, "ensembl": sym2ens[s], "entrez": int(ent), "in_human1": bool(in_model),
               "pathway": pathway, "role_in_example": why}
        for key, *_ in STEPS:
            row[f"e_{key}"] = float(e_by_step[key][ent])
        v = p1[ent]
        higher = int((p1 > v).sum())
        tied = int((p1 == v).sum())
        row["P1_rank_block"] = f"ranks {higher + 1}-{higher + tied} of {len(p1)} (tied with {tied - 1} genes)"
        if ent in pg.index:
            row["test_frac_lines_essential"] = round(float(pg.loc[ent, "prevalence"]), 4)
            row["test_n_lines"] = int(pg.loc[ent, "n"])
        rows.append(row)
    panel = pd.DataFrame(rows)

    # ---------------------------------------------------------------- 6. per-line tAP (committed aggregates)
    pl = pd.read_csv(R / "m6" / "train" / "per_line.tsv", sep="\t")
    pl = pl[(pl.ModelID == a.line)]
    w1 = pd.read_csv(R / "w1" / "train_per_line.tsv", sep="\t")
    w1 = w1[(w1.ModelID == a.line) & (w1.method == "H1") & (w1.geneset == "all")]
    tap = {}
    for key, cell in PER_LINE_CELL.items():
        r_all = pl[(pl.cell == cell) & (pl.label_cell == "all")]
        r_no = pl[(pl.cell == cell) & (pl.label_cell != "all")]
        tap[key] = {"cell": cell, "tAP_all_genes": round(float(r_all.tap.iloc[0]), 6)}
        if len(r_no):
            tap[key]["tAP_without_common_essentials"] = round(float(r_no.tap.iloc[0]), 6)
            tap[key]["label_cell_noCE_name"] = str(r_no.label_cell.iloc[0])
    tap["random_ranking_expectation"] = round(float(pl[(pl.cell == "H1") & (pl.label_cell == "all")].e_random.iloc[0]), 6)
    tap["G_star_from_w1"] = round(float(w1.tap.iloc[0]), 6)
    tap["G_star_check_equal"] = bool(abs(float(w1.tap.iloc[0]) - tap["step0"]["tAP_all_genes"]) < 1e-6)
    tap["P1_minus_Gstar"] = round(tap["step2"]["tAP_all_genes"] - tap["step0"]["tAP_all_genes"], 6)
    # where does this line sit among the 686 training lines?
    allp = pl_all = pd.read_csv(R / "m6" / "train" / "per_line.tsv", sep="\t",
                                usecols=["cell", "label_cell", "ModelID", "tap"])
    allp = allp[allp.label_cell == "all"].pivot(index="ModelID", columns="cell", values="tap")
    d = allp["H1+Bloo/M2"] - allp["H1"]
    tap["training_lines_P1_minus_Gstar_median"] = round(float(d.median()), 6)
    tap["training_lines_P1_minus_Gstar_IQR"] = [round(float(d.quantile(.25)), 6), round(float(d.quantile(.75)), 6)]
    tap["training_lines_with_gain"] = f"{int((d > 0).sum())}/{len(d)}"
    del pl_all

    # ---------------------------------------------------------------- 7. write
    res = {
        "what": "P1 worked example: one training line traced through each step (label-free model outputs; "
                "per-line tAP and per-gene test fractions from committed aggregate tables)",
        "script": "scripts/explainer_worked_example.py",
        "line": line_info,
        "biomass": {"reaction": BIOMASS, "pool_metabolite": POOL_MET, "pool_coefficient": round(-pool_coef, 6),
                    "n_biomass_reactants": len(biomass_reactants), "biomass_reactants": biomass_reactants,
                    "pool_reaction": POOL_RXN, "n_pool_species": len(pool_species), "pool_species": pool_species,
                    "edit": variant_spec},
        "medium": medium_info,
        "steps": steps_out,
        "transitions": transitions,
        "panel": rows,
        "tAP": tap,
        "runtime_s": {"load_and_bounds": round(t_load, 1), "total": round(time.time() - t_start, 1)},
        "call_rule": "lethal = e >= 0.99 (KO/WT growth ratio < 0.01); partial = 0 < e < 0.99",
        "label_note": "No label file read. tAP from results/m6/train/per_line.tsv; G* from results/w1/train_per_line.tsv; "
                      "test fractions from results/m6/test/unsealed/T1/per_gene.tsv (method P1[H1], 380 test lines).",
    }
    (out / "worked_example.json").write_text(json.dumps(res, indent=1, default=str))
    panel.to_csv(out / "worked_example_panel.tsv", sep="\t", index=False)
    steps_df = pd.DataFrame([{"step": k, **{kk: vv for kk, vv in v.items() if kk not in ("sweep_key",)},
                              "tAP": tap[k]["tAP_all_genes"]} for k, v in steps_out.items()])
    steps_df.to_csv(out / "worked_example_steps.tsv", sep="\t", index=False)

    # small Markdown table
    hdr = "| Gene | Pathway | Step 0: Ham's | Step 1: − pool | Medium only | Step 2: P1 | Test lines essential |\n"
    hdr += "|---|---|---|---|---|---|---|\n"
    body = "".join(
        f"| {r['symbol']} | {r['pathway']} | {r['e_step0']:.3f} | {r['e_step1']:.3f} | {r['e_medium_only']:.3f} | "
        f"{r['e_step2']:.3f} | {100 * r.get('test_frac_lines_essential', float('nan')):.1f}% |\n" for r in rows)
    foot = "| **Genes called lethal** | | " + " | ".join(
        str(steps_out[k]["n_lethal"]) for k in ("step0", "step1", "medium_only", "step2")) + " | |\n"
    foot += "| **Line's tAP** | | " + " | ".join(
        f"{tap[k]['tAP_all_genes']:.3f}" for k in ("step0", "step1", "medium_only", "step2")) + " | |\n"
    (out / "worked_example_table.md").write_text(hdr + body + foot)
    print(json.dumps({"line": line_info["name"], "tAP": {k: tap[k]["tAP_all_genes"] for k in PER_LINE_CELL},
                      "lethal": {k: v["n_lethal"] for k, v in steps_out.items()}}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
