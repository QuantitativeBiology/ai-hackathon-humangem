"""W2 step 5: label-free controls and per-line sweep lookup from the cached KO sweeps (DESIGN W2 Controls, §12 W2).

Positive (label-free): ASNS DMEM vs RPMI; sterol branch serum-free vs serum (by f); GOT1/GOT2 and SLC1A3 with the
gate; SLC7A11 cysteine (M0) vs cystine (M1); mevalonate core essential at every step (diagnostic panel).
Negative: M-open collapse; shuffled medium (relabelling of cached sweeps within lineage x serum-free, SEEDS[0:5]):
only label-free magnitudes are computed here (how many lines change medium, mean per-line |delta ratio_r|).

Outputs (results/w2/controls/):
  panel_ratio_r.tsv      ladder label x panel gene ratio_r (mevalonate core, notes/apriori_lists_final.md §1.2)
  positive_controls.tsv  one row per control x parent with expectation, observed values and pass flag
  negative_controls.tsv  M-open collapse; shuffled-medium magnitudes
  line_sweeps.tsv        ModelID x parent x ladder step -> sweep_id (the per-line score lookup)
  ladder_summary.tsv     parent x variant x f x gate: lethal-count distribution over media

Usage: .venv/bin/python scripts/w2_controls.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from nextgem import media
from nextgem.config import SEEDS
from nextgem.w2_runtime import PARENTS, Manifest, atomic_write_tsv, w2_dir

MEVALONATE = ["HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "GGPS1"]
STEROL = ["FDFT1", "SQLE", "LSS", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "SC5D", "EBP", "DHCR7", "DHCR24"]
PANEL_12 = ["CAD", "DHODH", "UMPS", "TYMS", "TK1", "DHFR", "PPAT", "GART", "PFAS", "PAICS", "ATIC", "ASNS", "GLUL",
            "ALDH18A1", "PYCR1", "PYCR2", "PYCR3", "GOT1", "GOT2", "SLC1A3", "PC", "PCK2", *STEROL, "FASN", "ACACA",
            "SCD", "GCLC", "GCLM", "GSS", "SLC7A11", "NAMPT", "PNPO", "PDXK", "BBOX1", "TMLHE", "ALDH9A1", "SHMT1",
            "SHMT2", "MMAA", "MMAB", "MMUT", "PCCA", "PCCB", "SLC5A6", "LIAS", "LIPT1", "LIPT2", "ODC1", "SRM", "SMS",
            "AMD1", "ALAD", "HMBS", "UROS", "UROD", "CPOX", "ISYNA1", "SEPHS2", "GPX4", "PHGDH", "PSAT1", "PSPH",
            "ASS1", "ASL", "ADSL", "IMPDH2", "GMPS", "MDH2", "SDHA", "CS", "ATP5F1A"]
PURINE = ["PPAT", "GART", "PFAS", "PAICS", "ATIC"]
PYR_GENES = ["GOT1", "GOT2", "MDH2", "SDHA", "CS", "ATP5F1A"]
PYRIMIDINE = ["CAD", "DHODH", "UMPS"]


def step_name(r) -> str:
    v = str(r.variant)
    if v.startswith("M2"):
        return f"{v}(f={r.f:g})"
    return v


def main() -> None:
    out = w2_dir() / "controls"
    out.mkdir(parents=True, exist_ok=True)
    with Manifest("scripts/w2_controls.py", solver="none") as man:
        lab = pd.read_csv(w2_dir() / "ko" / "sweep_labels.tsv", sep="\t")
        lab["step"] = lab.apply(step_name, axis=1)
        sym = {p: pd.read_csv(w2_dir().parent / "baseline" / f"{p}_gene_mapping.tsv", sep="\t")
               .set_index("gene").symbol for p in PARENTS}
        mats = {p: pd.read_parquet(w2_dir() / "ko" / f"ratio_r_{p}.parquet") for p in PARENTS}
        for p in PARENTS:
            mats[p].columns = [sym[p].get(g, g) for g in mats[p].columns]
            mats[p] = mats[p].T.groupby(level=0).min().T  # symbol-level (duplicate symbols: most severe)

        # ---------------------------------------------------------------- panel table
        panel_genes = list(dict.fromkeys(MEVALONATE + PANEL_12))
        rows = []
        for r in lab[lab.precheck_pass].itertuples():
            m = mats[r.parent]
            if r.sweep_id not in m.index:
                continue
            vals = m.loc[r.sweep_id].reindex(panel_genes)
            rows.append(dict(parent=r.parent, step=r.step, medium_key=r.medium_key, role=r.role,
                             gate_open=r.gate_open, sweep_id=r.sweep_id, **vals.to_dict()))
        panel = pd.DataFrame(rows)

        def val(p, step, key, gene):
            s = panel[(panel.parent == p) & (panel.step == step) & (panel.medium_key == key)]
            return float(s[gene].iloc[0]) if len(s) and gene in s else np.nan

        # ---------------------------------------------------------------- positive controls
        pos = []
        for p in PARENTS:
            a_d, a_r = val(p, "M1", "DMEM", "ASNS"), val(p, "M1", "RPMI", "ASNS")
            pos.append(dict(control="ASNS DMEM vs RPMI (M1)", parent=p, expectation="DMEM < 0.01 (no Asn); RPMI = 1",
                            observed=f"DMEM {a_d}; RPMI {a_r}", passed=a_d < 0.01 and a_r == 1.0))
            for f in media.SERUM_F_LEVELS:
                sf = [val(p, "M1", "DMEM:F12|serum-free", g) for g in STEROL]
                se = [val(p, f"M2(f={f:g})", "DMEM:F12", g) for g in STEROL]
                pos.append(dict(control=f"sterol branch: serum-free vs serum M2(f={f:g})", parent=p,
                                expectation="serum-free: lethal; serum: ratio rises with f (~f when cholesterol "
                                            "uptake caps growth)",
                                observed=f"serum-free n_lethal {sum(x < 0.01 for x in sf)}/{len(STEROL)}; "
                                         f"serum median ratio {np.nanmedian(se):.4g}",
                                passed=sum(x < 0.01 for x in sf) >= 10 and np.nanmedian(se) >= min(f, 1) * 0.5))
            for gene in ("GOT1", "GOT2", "SLC1A3"):
                go, gc = val(p, "M2(f=0.01)", "RPMI", gene), val(p, "M2g(f=0.01)", "RPMI", gene)
                pos.append(dict(control=f"{gene} with the gate (RPMI M2 f=0.01)", parent=p,
                                expectation="ratio(gate closed) <= ratio(gate open)",
                                observed=f"open {go}; closed {gc}", passed=gc <= go))
            s0 = float(mats[p].loc[lab[(lab.parent == p) & (lab.variant == "M0")].sweep_id.iloc[0]].get("SLC7A11",
                                                                                                      np.nan))
            s1 = val(p, "M1", "RPMI", "SLC7A11")
            pos.append(dict(control="SLC7A11 cysteine (M0) vs cystine (M1 RPMI)", parent=p,
                            expectation="M1 (cystine only) ratio <= M0", observed=f"M0 {s0}; M1 {s1}", passed=s1 <= s0))
            for genes, name, key, ref, exp in (
                    (PURINE, "purines lose essentiality in HPLM (hypoxanthine)", "M-HPLM", "RPMI", "HPLM >= RPMI"),
                    (PYRIMIDINE, "pyrimidine core rescued by uridine in Plasmax", "M-Plasmax", "RPMI",
                     "Plasmax >= RPMI"),
                    (PURINE, "purines in Plasmax (hypoxanthine 5 uM)", "M-Plasmax", "RPMI", "Plasmax >= RPMI")):
                a = [val(p, "M2(f=0.01)", key, g) for g in genes]
                b = [val(p, "M2(f=0.01)", ref, g) for g in genes]
                pos.append(dict(control=f"physiological: {name} (M2 f=0.01)", parent=p, expectation=exp,
                                observed="; ".join(f"{g} {x:.3g}/{y:.3g}" for g, x, y in zip(genes, a, b)),
                                passed=all(x >= y for x, y in zip(a, b))))
            a = [val(p, "M1", "RPMI", g) for g in PYR_GENES]
            b = [val(p, "M1", "SENS:RPMI-no-pyruvate", g) for g in PYR_GENES]
            pos.append(dict(control="physiological: pyruvate relaxes GOT/MDH/OXPHOS (RPMI vs RPMI-no-pyruvate, M1)",
                            parent=p, expectation="with pyruvate >= without",
                            observed="; ".join(f"{g} {x:.3g}/{y:.3g}" for g, x, y in zip(PYR_GENES, a, b)),
                            passed=all(x >= y for x, y in zip(a, b))))
            mev = panel[(panel.parent == p) & panel.role.isin(["documented", "physiological"])]
            core = ["HMGCR", "MVK", "PMVK", "MVD"]
            n_loss = int((mev[core] >= 0.01).any(axis=1).sum())
            pos.append(dict(control="mevalonate core (HMGCR, MVK, PMVK, MVD) essential at every step", parent=p,
                            expectation="lethal in every documented/physiological label",
                            observed=f"{n_loss} of {len(mev)} labels lose essentiality of >= 1 core gene",
                            passed=n_loss == 0))
        pos = pd.DataFrame(pos)

        # ---------------------------------------------------------------- negative controls
        neg = []
        sw = pd.read_csv(w2_dir() / "ko" / "sweeps.tsv", sep="\t").set_index("sweep_id")
        for p in PARENTS:
            n0 = sw.loc[lab[(lab.parent == p) & (lab.variant == "M0")].sweep_id.iloc[0], "n_lethal"]
            no = sw.loc[lab[(lab.parent == p) & (lab.variant == "M-open")].sweep_id.iloc[0], "n_lethal"]
            neg.append(dict(control="M-open collapse", parent=p, metric="n_lethal M-open vs M0",
                            value=f"{no} vs {n0}", passed=no < n0))
        lines = pd.read_csv(w2_dir() / "media_library" / "line_media.tsv", sep="\t")
        sh = pd.read_csv(w2_dir() / "media_library" / "shuffled_media.tsv", sep="\t")
        for p in PARENTS:
            for step in ("M1", "M2(f=0.01)"):
                key2sweep = lab[(lab.parent == p) & (lab.step == step) & lab.gate_open].set_index("medium_key").sweep_id
                m = mats[p]
                doc = m.loc[sh.medium_key.map(key2sweep)].to_numpy()
                for s in SEEDS[0:5]:
                    col = f"shuffled_seed{s}"
                    shm = m.loc[sh[col].map(key2sweep)].to_numpy()
                    changed = float((sh[col] != sh.medium_key).mean())
                    d = np.abs(doc - shm).mean(axis=1)
                    neg.append(dict(control="shuffled medium (label-free magnitude)", parent=p, metric=step,
                                    seed=s, value=f"lines with changed medium {changed:.3f}; mean per-line "
                                                  f"|delta ratio_r| {d.mean():.2e}; max {d.max():.2e}", passed=None))
        neg = pd.DataFrame(neg)

        # ---------------------------------------------------------------- per-line sweep lookup
        lk = []
        base = lab[lab.role.isin(["documented", "baseline", "negative_control", "physiological"])]
        for p in PARENTS:
            lp = base[base.parent == p]
            glob = lp[lp.variant.isin(["M0", "M-open"])][["step", "sweep_id"]]
            doc = lp[lp.role == "documented"]
            phys = lp[lp.role == "physiological"]
            for r in lines.itertuples():
                for g in glob.itertuples():
                    lk.append((r.ModelID, p, g.step, g.sweep_id))
                for src, prefix in ((doc[doc.medium_key == r.medium_key], ""), (phys, "phys")):
                    for d in src[src.gate_open].itertuples():
                        name = d.step if not prefix else f"{d.medium_key}:{d.step}"
                        lk.append((r.ModelID, p, name, d.sweep_id))
                        # +g state: the gated sweep if the line's gate is closed and the medium has Asp/Glu,
                        # otherwise identical to the gate-free sweep
                        gated = src[~src.gate_open & (src.medium_key == d.medium_key) & (src.f.fillna(-1) ==
                                                                                          (-1 if pd.isna(d.f) else d.f))]
                        gname = name.replace("M1", "M1g").replace("M2(", "M2g(")
                        sid = gated.sweep_id.iloc[0] if (not r.gate_effective and len(gated)) else d.sweep_id
                        lk.append((r.ModelID, p, gname, sid))
        lk = pd.DataFrame(lk, columns=["ModelID", "parent", "ladder_step", "sweep_id"])
        lk = lk.merge(lines[["ModelID", "stratum", "medium_key", "gate_effective"]], on="ModelID")

        summ = (lab[lab.role == "documented"].groupby(["parent", "step"])
                .agg(n_media=("medium_key", "nunique"), lethal_min=("n_lethal", "min"),
                     lethal_median=("n_lethal", "median"), lethal_max=("n_lethal", "max"),
                     precheck_fail=("precheck_pass", lambda s: int((~s).sum()))).reset_index())
        for n, df in {"panel_ratio_r.tsv": panel, "positive_controls.tsv": pos, "negative_controls.tsv": neg,
                      "line_sweeps.tsv": lk, "ladder_summary.tsv": summ}.items():
            man.add(atomic_write_tsv(df, out / n))
        print(pos.to_string())
        print(neg.to_string())
        print(summ.to_string())


if __name__ == "__main__":
    main()
