"""W2 M4 on the training tier: the medium ladder, attribution, controls and the shortlist decision
(DESIGN §9 W2 "Metrics", "Controls", "Decision"; §6; §9 TG). Labels only through ``nextgem.labels`` (``train()``,
and ``fold(k)`` via ``nextgem.cv.lolo`` for fold-scoped CE_fit). No validation or test label is read.

Reuses ``nextgem.w1`` (wide, fit_components, method_summary, delta_summary, ce_mask) and ``nextgem.metrics``
(tie-aware AP, lineage-stratified bootstrap, binary metrics), as W1/W3, so the numbers are comparable.

Lines: the W1 primary set (training lines of the 12 eligible lineages with an expression profile; 686), split into
the reporting strata of DESIGN W2: **attribution** (screen tier, sheet on disk; primary), **onboarded**
(sensitivity), **other** (imputed, fallback_sheet_missing, unmappable, mixture_on_disk; separate stratum), and
**all** (the integrated pipeline view). Primary cell: L_prob, all U genes, no stratum; secondary: minus CE_fit (LOLO).

Scores (``nextgem.w2_score``): e = 1 − ratio_r of the line's cached sweep for the ladder step
(``results/w2/controls/line_sweeps.tsv``); M0 equals W1's generic e (asserted). Reference: the generic parent on M0
(G* = Human1); every Human2 step is also compared with Human1 M0.

TG (§9): Δ̄_lin vs M0 >= 0.005, paired line-bootstrap CI lower bound > 0, Δ > 0 in >= 9 of 12 training LOLO folds
(fold with < 5 estimable lines = not positive). Eligible for the W2 shortlist (pre-registered in
``nextgem.w2_score``): M1, M2(f) gate-free on G* = Human1 (M3 excluded, A9i). Decision (≤ 1): one-SE rule.

Attribution (all on identical sets, attribution stratum): (i) medium alone = Δ̄_lin(best gate-free documented
variant − M0); (ii) net flip fraction of generic binary calls (L_prob primary, L_ge06 secondary; parent genes);
(iii) training-only medium oracle: per line the best of {M0, M1, M2(f)} by its own tAP; (doc − M0)/(oracle − M0).
Controls: M-open (negative); shuffled medium within lineage × serum-free (SEEDS[0:5], cached sweeps; pooled and per
eligible lineage); positive label-free panel (approved lists §1.1/§1.2) and the mevalonate core at every step;
ASNS by medium on training labels (Lagziel replication).

Outputs (results/w2/train/): step_scores.tsv, per_line.parquet, attribution.json, flip_counts.tsv,
shuffled_control.tsv, panel_by_step.tsv, mevalonate_by_step.tsv, asns_by_medium.tsv, binary.tsv,
sensitivity.tsv, shortlist.json, train_summary.json.

Usage: OMP_NUM_THREADS=1 .venv/bin/python scripts/w2_score_train.py
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / "scripts"))

W2 = PROJECT / "results" / "w2"
TRDIR = W2 / "train"
LABEL_COLS = ["gene_effect", "L_prob", "L_ge06"]
PARENTS = ("human1", "human2")
SHORT = {"human1": "H1", "human2": "H2"}
MEV_CORE = ("HMGCR", "MVK", "PMVK", "MVD")
MEV_ALL = ("HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "GGPS1")
SENS_BASE = {"SENS:RPMI-ATCC": "RPMI", "SENS:DMEM+pyruvate": "DMEM", "SENS:DMEM-no-standard-Gln": "DMEM",
             "SENS:DMEM:F12-no-standard-Gln": "DMEM:F12", "SENS:DMEM:F12-cysteine-documented": "DMEM:F12"}


def load_inputs():
    from nextgem import labels, paths, w1, w2_score as ws

    split = labels.load_split().set_index("ModelID")
    lf = w1.line_features()
    Udf = labels.gene_universe()
    U = list(Udf.entrez.astype(int))
    elig = split[(split.tier == "train") & (split.stratum != "rare")]
    E = sorted(l for l in elig.index if lf.loc[l, "has_expression_matrix"])
    lm = pd.read_csv(W2 / "media_library" / "line_media.tsv", sep="\t").set_index("ModelID")
    ls = pd.read_csv(W2 / "controls" / "line_sweeps.tsv", sep="\t")
    sl = pd.read_csv(W2 / "ko" / "sweep_labels.tsv", sep="\t")
    emap, Esw, sym = {}, {}, {}
    for p in PARENTS:
        ko = pd.read_parquet(paths.baseline_dir() / f"{p}_ko.parquet", columns=["gene", "entrez", "symbol"])
        emap[p] = ws.gene_to_entrez(ko)
        sym[p] = ko.set_index("gene").symbol
        Esw[p] = ws.sweep_e_matrix(pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet"), emap[p], U)
    gen = w1.generic_table(Udf.entrez)
    return dict(split=split, lf=lf, U=U, E=E, lm=lm, ls=ls, sl=sl, emap=emap, Esw=Esw, sym=sym, gen=gen,
                lineage=split.lineage, Udf=Udf)


def main() -> int:
    from nextgem import cv, labels, metrics, w1, w2_score as ws
    from nextgem.config import SEEDS
    from nextgem.w2_runtime import Manifest, atomic_write_tsv
    from w3_common import frozen_label_hashes

    t0 = time.time()
    with Manifest("w2_score_train.py", solver="none", workers=1,
                  inputs=[W2 / "ko" / f"ratio_r_{p}.parquet" for p in PARENTS]
                  + [W2 / "controls" / "line_sweeps.tsv", W2 / "ko" / "sweep_labels.tsv",
                     W2 / "media_library" / "line_media.tsv", W2 / "media_library" / "shuffled_media.tsv",
                     W2 / "controls" / "panel_ratio_r.tsv"]) as man:
        man.extra.update(label_files_frozen_sha256=frozen_label_hashes(),
                         label_api="nextgem.labels.train()/fold(k) only (training tier)", seeds=list(SEEDS[0:5]))
        G = load_inputs()
        U, E, lineage, ls = G["U"], G["E"], G["lineage"], G["ls"]
        steps = ws.step_order(ls.ladder_step.unique())
        sweep_of = {(p, s): g.set_index("ModelID").sweep_id for (p, s), g in ls.groupby(["parent", "ladder_step"])}

        def S_of(p, s, lines=E):
            return ws.line_matrix(sweep_of[(p, s)], G["Esw"][p], lines)

        # ---- equality: M0 == W1 generic e (both parents)
        checks = {}
        for p in PARENTS:
            v = G["gen"][f"e_{SHORT[p]}"].reindex(U).to_numpy()
            checks[f"M0_equals_w1_e_{SHORT[p]}"] = bool(np.array_equal(S_of(p, "M0").to_numpy(), np.tile(v, (len(E), 1))))
        if not all(checks.values()):
            raise RuntimeError(f"M0 differs from W1 generic e: {checks}")

        # ---- labels (training tier only) and fold-scoped CE_fit
        line_class = G["lf"].base_medium_class

        def fit_predict(fit_long, ids):
            return w1.fit_components(fit_long, line_class, labels_=("L_prob",)).ce["L_prob"], ids

        folds = cv.lolo(fit_predict, lambda pr, held: pr, columns=["gene_effect", "L_prob"])
        fold_names = list(folds)
        line_ce = {}
        for k, (ce, ids) in folds.items():
            for l in ids:
                line_ce[l] = ce
        train_long = labels.train(columns=LABEL_COLS)
        Y = {lab: w1.wide(train_long, lab, E, U) for lab in ("L_prob", "L_ge06")}
        Y["L_prob_noCE"] = Y["L_prob"].where(w1.ce_mask(E, U, line_ce))
        del train_long

        strata = {k: [l for l in E if G["lm"].stratum.get(l) in v] for k, v in ws.STRATA.items()}
        strata["all"] = list(E)
        A = strata["attribution"]

        # ---- per-line tables for every (parent, step, cell) on all 686 lines (tAP_line is per line)
        LT, gene_rows = {}, []
        for p in PARENTS:
            for s in steps:
                S = S_of(p, s)
                for cell in ("L_prob", "L_prob_noCE"):
                    LT[(p, s, cell)] = metrics.per_line(S, Y[cell])
                row, _, _ = w1.method_summary(S.loc[A], Y["L_prob"].loc[A], lineage, gene_unit=True)
                gene_rows.append({"parent": p, "step": s, **{k: row[k] for k in (
                    "tap_lin", "tap_lin_lo", "tap_lin_hi", "e_random_lin", "lift_lin", "m_delta_ap_gene",
                    "n_genes_estimable", "pooled_tap", "pooled_lift")}})
        print(f"per-line tables: {len(LT)}, {time.time() - t0:.0f} s", flush=True)
        gene_df = pd.DataFrame(gene_rows)

        def sub(lt, lines):
            return lt.loc[[l for l in lines if l in lt.index]]

        # ---- step table: Δ vs own M0 and vs H1 M0, per stratum, TG on every stratum (decision on attribution)
        rows = []
        for st_name, lines in strata.items():
            for p in PARENTS:
                for s in steps:
                    info = ws.parse_step(s)
                    lt = sub(LT[(p, s, "L_prob")], lines)
                    est = lt[lt.estimable]
                    r = {"stratum": st_name, "parent": p, "step": s, "family": info.family, "gated": info.gated,
                         "f": info.f, "eligible": info.eligible and p == "human1", "n_lines": len(lt),
                         "n_lineages": int(lineage.reindex(est.index).nunique()),
                         "tap_lin": metrics.lineage_equal_mean(est.tap, lineage)}
                    for ref_p, tag in ((p, "vs_M0"), ("human1", "vs_H1M0")):
                        if tag == "vs_H1M0" and p == "human1":
                            continue
                        ref = sub(LT[(ref_p, "M0", "L_prob")], lines)
                        if s == "M0" and ref_p == p:
                            continue
                        d = w1.delta_summary(lt, ref, lineage)
                        r.update({f"{k}_{tag}": v for k, v in d.items()})
                        if tag == "vs_M0":
                            dd = metrics.paired_line_delta(lt, ref)
                            r.update(ws.folds_positive(dd, lineage, fold_names))
                            r["tg_pass"] = ws.tg_pass(d["delta_lin"], d["delta_lo"], r["folds_positive"])
                            r["delta_noCE_vs_M0"] = w1.delta_summary(sub(LT[(p, s, "L_prob_noCE")], lines),
                                                                     sub(LT[(p, "M0", "L_prob_noCE")], lines),
                                                                     lineage)["delta_lin"]
                    if info.gated:  # gate effect: gated − gate-free twin (context lever)
                        tw = sub(LT[(p, info.gate_free_twin, "L_prob")], lines)
                        g = w1.delta_summary(lt, tw, lineage)
                        r.update({"gate_effect": g["delta_lin"], "gate_effect_lo": g["delta_lo"],
                                  "gate_effect_hi": g["delta_hi"]})
                    rows.append(r)
        step_df = pd.DataFrame(rows).merge(gene_df.assign(stratum="attribution"), on=["stratum", "parent", "step"],
                                           how="left", suffixes=("", "_attr"))
        TRDIR.mkdir(parents=True, exist_ok=True)
        print(f"step table: {len(step_df)} rows, {time.time() - t0:.0f} s", flush=True)

        # ---- shortlist (attribution stratum, H1 = G*, eligible gate-free documented steps)
        sa = step_df[(step_df.stratum == "attribution") & (step_df.parent == "human1") & step_df.eligible].copy()
        sa = sa.rename(columns={"delta_lin_vs_M0": "delta_lin"})
        best = sa.sort_values("delta_lin", ascending=False).iloc[0].step
        se_of = {}
        for s in sa.step:
            d = metrics.paired_line_delta(sub(LT[("human1", best, "L_prob")], A), sub(LT[("human1", s, "L_prob")], A))
            bt = metrics.bootstrap_mean(d, lineage.reindex(d.index), B=metrics.B_OTHER, seed=SEEDS[0])
            se_of[s] = (bt.estimate, float(np.std(bt.replicates, ddof=1)))
        decision = ws.pick_shortlist(sa[["step", "delta_lin", "tg_pass"]], se_of)
        h1_scores = step_df[(step_df.stratum == "attribution") & (step_df.parent == "human1") & (step_df.step != "M0")]
        sd = float(np.std(h1_scores.delta_lin_vs_M0.dropna(), ddof=1))
        m_screened = int(len(sa))
        decision.update({"eligible_steps": sa.step.tolist(), "n_candidates_screened_h1": int(len(h1_scores)),
                         "m_eligible": m_screened, "sigma_all_h1_steps": sd,
                         "winners_curse_sigma_sqrt_2lnm": sd * float(np.sqrt(2 * np.log(max(2, m_screened)))),
                         "tg_rule": "delta_lin >= 0.005, CI lower > 0, >= 9 of 12 training LOLO folds (fold < 5 lines = not positive)",
                         "eligibility": "documented gate-free M1, M2(f) on Human1 (G*); M3 excluded (A9i); gate = context lever"})
        cand = decision["shortlist"][0] if decision["shortlist"] else None

        # ---- attribution (attribution stratum; H1 primary, H2 secondary)
        attribution, flip_rows = {}, []
        for p in PARENTS:
            doc = step_df[(step_df.stratum == "attribution") & (step_df.parent == p) & step_df.step.isin(ws.ELIGIBLE)]
            best_doc = doc.sort_values("delta_lin_vs_M0", ascending=False).iloc[0]
            lt0 = sub(LT[(p, "M0", "L_prob")], A)
            lt_doc = sub(LT[(p, best_doc.step, "L_prob")], A)
            d_doc = metrics.paired_line_delta(lt_doc, lt0)
            taps = pd.concat({s: sub(LT[(p, s, "L_prob")], A).tap for s in ws.DOC_LADDER}, axis=1)
            est = lt0.estimable
            d_or = (taps.max(axis=1) - taps["M0"])[est]
            orc = ws.oracle_ratio(d_doc, d_or, lineage, B=metrics.B_OTHER, seed=SEEDS[0])
            orc_lin = metrics.bootstrap_mean(d_or, lineage.reindex(d_or.index), B=metrics.B_OTHER, seed=SEEDS[0])
            argbest = taps[est].idxmax(axis=1)
            # sensitivity oracle over every gate-free step incl. physiological, and over gated steps too
            taps_all = pd.concat({s: sub(LT[(p, s, "L_prob")], A).tap for s in steps if s != "M-open"}, axis=1)
            d_or_all = (taps_all.max(axis=1) - taps_all["M0"])[est]
            attribution[p] = {
                "medium_alone": {"best_gate_free_documented": best_doc.step, "delta_lin": best_doc.delta_lin_vs_M0,
                                 "ci": [best_doc.delta_lo_vs_M0, best_doc.delta_hi_vs_M0],
                                 "delta_line": best_doc.delta_line_vs_M0,
                                 "folds_positive": int(best_doc.folds_positive)},
                "oracle": {**orc, "oracle_minus_m0_ci": [orc_lin.ci_low, orc_lin.ci_high],
                           "oracle_choice_counts": argbest.value_counts().to_dict(),
                           "oracle_all_steps_minus_m0": metrics.lineage_equal_mean(d_or_all, lineage),
                           "ladder": list(ws.DOC_LADDER)},
            }
            P_genes = [g for g in U if G["gen"].loc[g, f"in_{SHORT[p]}"]]
            c0 = S_of(p, "M0", A)[P_genes] >= ws.CALL
            for s in [x for x in steps if x not in ("M0",)]:
                c1 = S_of(p, s, A)[P_genes] >= ws.CALL
                for lab in ("L_prob", "L_ge06"):
                    fc = ws.flip_counts(c0, c1, Y[lab].loc[A, P_genes])
                    lin_frac = (fc["per_line_fixed"] - fc["per_line_broken"]) / fc["per_line_errors"].replace(0, np.nan)
                    flip_rows.append({"parent": p, "step": s, "label": lab,
                                      **{k: v for k, v in fc.items() if not k.startswith("per_line")},
                                      "net_flip_fraction_lin": metrics.lineage_equal_mean(lin_frac, lineage)})
            fr = [r for r in flip_rows if r["parent"] == p and r["step"] == best_doc.step]
            attribution[p]["net_flip"] = {r["label"]: {k: r[k] for k in ("errors_ref", "fixed", "broken", "fixed_fp",
                                                                         "fixed_fn", "net_flip_fraction",
                                                                         "net_flip_fraction_lin")} for r in fr}
        flip_df = pd.DataFrame(flip_rows)
        # onboarded (sensitivity) and other stratum: medium alone for the same best step
        for st_name in ("onboarded", "other", "all"):
            for p in PARENTS:
                bs = attribution[p]["medium_alone"]["best_gate_free_documented"]
                r = step_df[(step_df.stratum == st_name) & (step_df.parent == p) & (step_df.step == bs)].iloc[0]
                attribution[p][f"medium_alone_{st_name}"] = {"step": bs, "delta_lin": r.delta_lin_vs_M0,
                                                             "ci": [r.delta_lo_vs_M0, r.delta_hi_vs_M0],
                                                             "n_lines": int(r.n_lines)}
        print(f"attribution done, {time.time() - t0:.0f} s", flush=True)

        # ---- negative control: shuffled medium within lineage x serum-free (SEEDS[0:5]), attribution lines
        sh = pd.read_csv(W2 / "media_library" / "shuffled_media.tsv", sep="\t").set_index("ModelID")
        A_sh = [l for l in A if l in sh.index]
        elig_lin = ws.shuffle_eligible_lineages(G["lm"].medium_key.reindex(A_sh), lineage)
        sh_rows = []
        for p in PARENTS:
            for s in [x for x in steps if x in ws.ELIGIBLE]:
                lt_doc = sub(LT[(p, s, "L_prob")], A_sh)
                lt0 = sub(LT[(p, "M0", "L_prob")], A_sh)
                seed_lts = []
                for seed in SEEDS[0:5]:
                    sw = ws.shuffled_sweeps(sh.loc[A_sh, f"shuffled_seed{seed}"], s, p, G["sl"])
                    Ssh = ws.line_matrix(sw, G["Esw"][p], A_sh)
                    seed_lts.append(metrics.per_line(Ssh, Y["L_prob"].loc[A_sh]))
                mean_lt = seed_lts[0].copy()
                mean_lt["tap"] = np.mean([x.tap.to_numpy() for x in seed_lts], axis=0)
                d = w1.delta_summary(lt_doc, mean_lt, lineage)
                per_seed = [w1.delta_summary(lt_doc, x, lineage)["delta_lin"] for x in seed_lts]
                sh_vs_m0 = [w1.delta_summary(x, lt0, lineage)["delta_lin"] for x in seed_lts]
                r = {"parent": p, "step": s, "scope": "pooled", "n_lines": len(lt_doc),
                     "doc_minus_shuffled": d["delta_lin"], "lo": d["delta_lo"], "hi": d["delta_hi"],
                     "per_seed_min": float(np.min(per_seed)), "per_seed_max": float(np.max(per_seed)),
                     "shuffled_minus_m0_mean": float(np.mean(sh_vs_m0)),
                     "doc_minus_m0": w1.delta_summary(lt_doc, lt0, lineage)["delta_lin"],
                     "beats_shuffled": bool(d["delta_lo"] > 0)}
                sh_rows.append(r)
                dd = metrics.paired_line_delta(lt_doc, mean_lt)
                for k in elig_lin:
                    v = dd[lineage.reindex(dd.index) == k]
                    sh_rows.append({"parent": p, "step": s, "scope": f"lineage:{k}", "n_lines": int(len(v)),
                                    "doc_minus_shuffled": float(v.mean()) if len(v) else np.nan})
        sh_df = pd.DataFrame(sh_rows)
        print(f"shuffled control done, {time.time() - t0:.0f} s", flush=True)

        # ---- positive control on training labels: ASNS by medium (Lagziel replication)
        asns_rows = []
        for p in PARENTS:
            g = G["sym"][p]
            gid = g.index[g == "ASNS"][0]
            ez = int(G["emap"][p][gid])
            for s in ("M0", "M1", "M2(f=0.1)", "M2(f=1)"):
                sw = sweep_of[(p, s)].reindex(A)
                r_asns = 1.0 - G["Esw"][p].loc[sw.to_numpy(), ez].to_numpy()
                grp = pd.Series(np.where(r_asns < 0.01, "asn_absent_lethal", "asn_present_viable"), index=A)
                y = Y["L_prob"].loc[A, ez]
                for k, idx in grp.groupby(grp).groups.items():
                    yy = y.loc[idx].dropna()
                    asns_rows.append({"parent": p, "step": s, "group": k, "n_lines": int(len(yy)),
                                      "n_lineages": int(lineage.reindex(yy.index).nunique()),
                                      "positive_rate": float(yy.mean()) if len(yy) else np.nan,
                                      "media": ";".join(sorted(set(G["lm"].medium_key.reindex(list(idx)))))[:300]})
                yv = y.dropna()
                if yv.nunique() == 2 and len(set(np.round(r_asns, 6))) > 1:
                    e = pd.Series(1.0 - r_asns, index=A).loc[yv.index]
                    asns_rows.append({"parent": p, "step": s, "group": "tAP_gene_over_lines",
                                      "n_lines": int(len(yv)), "positive_rate": float(yv.mean()),
                                      "tap_gene": metrics.tap(e.to_numpy(), yv.to_numpy()),
                                      "e_random": metrics.expected_ap_random(len(yv), int(yv.sum()))})
        asns_df = pd.DataFrame(asns_rows)

        # ---- label-free positive panel and mevalonate core at every step (all 835 attribution lines, all tiers)
        panel = pd.read_csv(W2 / "controls" / "panel_ratio_r.tsv", sep="\t")
        panel_genes = [c for c in panel.columns if c not in ("parent", "step", "medium_key", "role", "gate_open", "sweep_id")]
        A_all = G["lm"].index[G["lm"].stratum == "attribution"].tolist()
        pr, mev = [], []
        for p in PARENTS:
            ratio = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
            g = G["sym"][p]
            for s in steps:
                sw = sweep_of[(p, s)].reindex(A_all)
                usw = sw.unique()
                step_rows = []
                for symb in panel_genes:
                    gids = g.index[g == symb]
                    if len(gids) == 0 or gids[0] not in ratio.columns:
                        continue
                    rv = ratio.loc[usw, gids[0]]
                    per_line_r = rv.reindex(sw.to_numpy()).to_numpy()
                    step_rows.append({"parent": p, "step": s, "gene": symb, "n_media": len(usw),
                                      "n_media_lethal": int((rv < 0.01).sum()),
                                      "frac_lines_lethal": float((per_line_r < 0.01).mean()),
                                      "ratio_min": float(rv.min()), "ratio_max": float(rv.max())})
                pr += step_rows
                core = {r["gene"]: r for r in step_rows if r["gene"] in MEV_ALL}
                mev.append({"parent": p, "step": s, "n_media": len(usw),
                            "core_lethal_all_media": all(core[x]["n_media_lethal"] == core[x]["n_media"] for x in MEV_CORE if x in core),
                            **{f"{x}_max_ratio": core[x]["ratio_max"] for x in MEV_ALL if x in core}})
        panel_df, mev_df = pd.DataFrame(pr), pd.DataFrame(mev)

        # ---- binary endpoints (§6.1 fixed cutoff) vs L_ge06 (§6.8 label) and L_prob: M0, best doc, shortlist
        bin_rows = []
        bsteps = sorted({"M0", attribution["human1"]["medium_alone"]["best_gate_free_documented"]} | ({cand} if cand else set()))
        for p in PARENTS:
            P_genes = [g for g in U if G["gen"].loc[g, f"in_{SHORT[p]}"]]
            pls = {}
            for s in bsteps:
                C = S_of(p, s, A)[P_genes] >= ws.CALL
                for lab in ("L_ge06", "L_prob"):
                    pl = metrics.binary_per_line(C, Y[lab].loc[A, P_genes])
                    pls[(s, lab)] = pl
                    summ = metrics.summarise_binary(pl, lineage.reindex(pl.index))
                    summ = summ[summ.nan_policy == "exclude"]
                    for _, r in summ.iterrows():
                        v = pl[r.metric].dropna()
                        bt = metrics.bootstrap_mean(v, lineage.reindex(v.index), B=metrics.B_OTHER, seed=SEEDS[0])
                        row = {"parent": p, "step": s, "label": lab, "metric": r.metric,
                               "lineage_equal_mean": r.lineage_equal_mean, "lo": bt.ci_low, "hi": bt.ci_high,
                               "line_mean": r.line_mean, "median": r["median"], "n_nan": r.n_nan}
                        if s != "M0":
                            dv = (pl[r.metric] - pls[("M0", lab)][r.metric]).dropna()
                            bd = metrics.bootstrap_mean(dv, lineage.reindex(dv.index), B=metrics.B_OTHER, seed=SEEDS[0])
                            row.update({"delta_vs_M0": bd.estimate, "delta_lo": bd.ci_low, "delta_hi": bd.ci_high})
                        bin_rows.append(row)
        bin_df = pd.DataFrame(bin_rows)

        # ---- sensitivities for the shortlisted (or best documented) step: A9 SENS media keys
        sens_rows = []
        s_sens = cand or attribution["human1"]["medium_alone"]["best_gate_free_documented"]
        info = ws.parse_step(s_sens)
        for p in PARENTS:
            lab_tab = G["sl"][(G["sl"].parent == p) & (G["sl"].variant == ("M1" if s_sens == "M1" else "M2"))]
            if info.f is not None:
                lab_tab = lab_tab[np.isclose(lab_tab.f.astype(float), info.f)]
            lab_tab = lab_tab[lab_tab.gate_open.astype(str) == "True"].drop_duplicates("medium_key").set_index("medium_key").sweep_id
            base_sw = sweep_of[(p, s_sens)].reindex(A)
            lt_base = sub(LT[(p, s_sens, "L_prob")], A)
            for sk, bk in SENS_BASE.items():
                if sk not in lab_tab.index:
                    continue
                hit = [l for l in A if G["lm"].medium_key.get(l) == bk]
                sw = base_sw.copy()
                sw.loc[hit] = lab_tab[sk]
                lt_s = metrics.per_line(ws.line_matrix(sw, G["Esw"][p], A), Y["L_prob"].loc[A])
                d = w1.delta_summary(lt_s, lt_base, lineage)
                d0 = w1.delta_summary(lt_s, sub(LT[(p, "M0", "L_prob")], A), lineage)
                sens_rows.append({"parent": p, "step": s_sens, "sensitivity": sk, "n_lines_changed": len(hit),
                                  "n_sweeps_differ": int((sw != base_sw).sum()),
                                  "delta_vs_primary": d["delta_lin"], "lo": d["delta_lo"], "hi": d["delta_hi"],
                                  "delta_vs_M0": d0["delta_lin"], "delta_vs_M0_lo": d0["delta_lo"]})

            def add_alt(name, sw):
                lt_s = metrics.per_line(ws.line_matrix(sw, G["Esw"][p], A), Y["L_prob"].loc[A])
                d = w1.delta_summary(lt_base, lt_s, lineage)  # primary − alternative
                d0 = w1.delta_summary(lt_s, sub(LT[(p, "M0", "L_prob")], A), lineage)
                sens_rows.append({"parent": p, "step": s_sens, "sensitivity": name,
                                  "n_lines_changed": int((sw != base_sw).sum()), "n_sweeps_differ": int((sw != base_sw).sum()),
                                  "delta_vs_primary": -d["delta_lin"], "lo": -d["delta_hi"], "hi": -d["delta_lo"],
                                  "delta_vs_M0": d0["delta_lin"], "delta_vs_M0_lo": d0["delta_lo"]})

            # A9(f): serum components at f_eff = f x serum% / 10 (scripts/w2_serum_pct.py; training attribution lines)
            sp_dir = W2 / "sens_serum_pct"
            if (sp_dir / "line_sweeps.tsv").exists():
                spl = pd.read_csv(sp_dir / "line_sweeps.tsv", sep="\t")
                spl = spl[spl.parent == p].set_index("ModelID").sweep_id
                if set(A) <= set(spl.index):
                    spl = spl.reindex(A)
                    rs = pd.read_parquet(sp_dir / f"ratio_r_{p}.parquet")
                    Esp = pd.concat([G["Esw"][p], ws.sweep_e_matrix(rs, G["emap"][p], U)])
                    Esp = Esp[~Esp.index.duplicated()]
                    lt_s = metrics.per_line(ws.line_matrix(spl, Esp, A), Y["L_prob"].loc[A])
                    d = w1.delta_summary(lt_s, lt_base, lineage)
                    d0 = w1.delta_summary(lt_s, sub(LT[(p, "M0", "L_prob")], A), lineage)
                    sens_rows.append({"parent": p, "step": s_sens, "sensitivity": "A9f:serum_pct_scaling",
                                      "n_lines_changed": int((spl.reindex(A) != base_sw).sum()),
                                      "n_sweeps_differ": int((spl.reindex(A) != base_sw).sum()),
                                      "delta_vs_primary": d["delta_lin"], "lo": d["delta_lo"], "hi": d["delta_hi"],
                                      "delta_vs_M0": d0["delta_lin"], "delta_vs_M0_lo": d0["delta_lo"]})
            # POST HOC (added after the training aggregates were seen; descriptive only, enters no decision): the same
            # step with one line-constant documented medium for every line (classic RPMI + 10% FBS, the imputed medium)
            rp = lab_tab.get("RPMI")
            if rp is not None:
                add_alt("POSTHOC:constant_RPMI_all_lines", pd.Series(rp, index=A))
        sens_df = pd.DataFrame(sens_rows)

        # ---- write
        pl_out = pd.concat([lt.assign(parent=p, step=s, cell=c) for (p, s, c), lt in LT.items()])
        pl_out = pl_out.rename_axis("ModelID").reset_index()
        pl_out["lineage"] = pl_out.ModelID.map(lineage)
        pl_out["w2_stratum"] = pl_out.ModelID.map(G["lm"].stratum)
        pl_out = pl_out[["parent", "step", "cell", "ModelID", "lineage", "w2_stratum", "n", "n_pos", "estimable", "tap", "e_random"]]
        outs = [
            atomic_write_tsv(step_df, TRDIR / "step_scores.tsv", float_format="%.6g"),
            atomic_write_tsv(flip_df, TRDIR / "flip_counts.tsv", float_format="%.6g"),
            atomic_write_tsv(sh_df, TRDIR / "shuffled_control.tsv", float_format="%.6g"),
            atomic_write_tsv(panel_df, TRDIR / "panel_by_step.tsv", float_format="%.6g"),
            atomic_write_tsv(mev_df, TRDIR / "mevalonate_by_step.tsv", float_format="%.6g"),
            atomic_write_tsv(asns_df, TRDIR / "asns_by_medium.tsv", float_format="%.6g"),
            atomic_write_tsv(bin_df, TRDIR / "binary.tsv", float_format="%.6g"),
            atomic_write_tsv(sens_df, TRDIR / "sensitivity.tsv", float_format="%.6g"),
        ]
        pq = TRDIR / "per_line.parquet"
        pl_out.to_parquet(pq, index=False, compression="zstd")
        outs.append(pq)
        cand_flags = {}
        if cand:
            mrow = mev_df[(mev_df.parent == "human1") & (mev_df.step == cand)].iloc[0]
            cand_flags["mevalonate_core_lethal_all_media_generic_biomass"] = bool(mrow.core_lethal_all_media)
            shr = sh_df[(sh_df.parent == "human1") & (sh_df.step == cand) & (sh_df.scope == "pooled")].iloc[0]
            cand_flags["beats_shuffled_medium"] = bool(shr.beats_shuffled)
            mo = step_df[(step_df.stratum == "attribution") & (step_df.parent == "human1") & (step_df.step == "M-open")].iloc[0]
            cand_flags["m_open_delta_vs_M0"] = float(mo.delta_lin_vs_M0)
        sl_out = {**decision, "candidate_flags": cand_flags, "primary_parent": "human1 (G*)",
                  "candidate_row": (step_df[(step_df.stratum == "attribution") & (step_df.parent == "human1")
                                            & (step_df.step == cand)].iloc[0].to_dict() if cand else None)}
        (TRDIR / "shortlist.json").write_text(json.dumps(sl_out, indent=1, default=float))
        (TRDIR / "attribution.json").write_text(json.dumps(attribution, indent=1, default=float))
        summary = {"n_lines": len(E), "strata_n_lines": {k: len(v) for k, v in strata.items()}, "n_steps": len(steps),
                   "steps": steps, "checks": checks, "folds": fold_names, "shuffle_eligible_lineages": elig_lin,
                   "n_shuffle_eligible_lineages": len(elig_lin), "shortlist": decision["shortlist"],
                   "wall_s": round(time.time() - t0, 1)}
        (TRDIR / "train_summary.json").write_text(json.dumps(summary, indent=1, default=float))
        outs += [TRDIR / "shortlist.json", TRDIR / "attribution.json", TRDIR / "train_summary.json"]
        man.add(*outs)
        man.extra.update(n_lines=len(E), shortlist=decision["shortlist"])
        print(json.dumps({"checks": checks, "strata": summary["strata_n_lines"], "decision": decision["decision"],
                          "wall_s": summary["wall_s"]}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
