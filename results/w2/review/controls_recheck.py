"""Re-check (A9) of results/w2/controls/positive_controls.tsv statuses, recomputed from the ratio_r matrices and
sweep_labels.tsv (label-free, no LP). Status rule (A9j): pass = expectation met with a strict change; non_informative =
both sides equal where a change is expected; model_structure_failure = expected rescue absent for a model reason.
Writes results/w2/review/controls_recheck.tsv.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
W2 = ROOT / "results" / "w2"
OUT = Path(__file__).resolve().parent / "controls_recheck.tsv"
STEROL = ["FDFT1", "SQLE", "LSS", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "SC5D", "EBP", "DHCR7", "DHCR24"]


def main() -> None:
    lab = pd.read_csv(W2 / "ko" / "sweep_labels.tsv", sep="\t")
    ctl = pd.read_csv(W2 / "controls" / "positive_controls.tsv", sep="\t")
    rows = []
    for p in ("human1", "human2"):
        R = pd.read_parquet(W2 / "ko" / f"ratio_r_{p}.parquet")
        sym = pd.read_csv(ROOT / "results" / "baseline" / f"{p}_gene_mapping.tsv", sep="\t").set_index("gene").symbol
        R.columns = [sym.get(g, g) for g in R.columns]
        R = R.T.groupby(level=0).min().T

        def v(key, variant, gene, f=None, gate=True):
            r = lab[(lab.parent == p) & (lab.medium_key == key) & (lab.variant == variant) & (lab.gate_open == gate)]
            r = r[r.f.isna()] if f is None else r[r.f == f]
            assert r.sweep_id.nunique() == 1, (key, variant, f, gate, r.sweep_id.unique())
            return float(R.loc[r.sweep_id.iloc[0], gene])

        def status(exp_change, a, b, strict_ok):
            if all(x == y for x, y in zip(a, b)):
                return "non_informative" if exp_change else "pass"
            return "pass" if strict_ok else "fail"

        mine = {}
        a, r_ = v("DMEM", "M1", "ASNS"), v("RPMI", "M1", "ASNS")
        mine["ASNS DMEM vs RPMI (M1)"] = ("pass" if a < 0.01 and r_ == 1 else "fail", f"DMEM {a} RPMI {r_}")
        for f in (0.001, 0.01, 0.1, 1.0):
            sf = [v("DMEM:F12|serum-free", "M1", g) for g in STEROL]
            se = [v("DMEM:F12", "M2", g, f) for g in STEROL]
            rescued = [g for g, x, y in zip(STEROL, sf, se) if x < 0.01 and y >= 0.01]
            lethal_both = [g for g, x, y in zip(STEROL, sf, se) if x < 0.01 and y < 0.01]
            mine[f"sterol branch: serum-free vs serum M2(f={f:g})"] = (
                "model_structure_failure (W3/W5)" if lethal_both else "pass",
                f"serum-free lethal {sum(x < 0.01 for x in sf)}; rescued by serum {rescued}; still lethal {lethal_both}")
        for g in ("GOT1", "GOT2", "SLC1A3"):
            o, c = v("RPMI", "M2", g, 0.01, True), v("RPMI", "M2g", g, 0.01, False)
            mine[f"{g} with the gate (RPMI M2 f=0.01)"] = (status(True, [c], [o], c < o), f"open {o} closed {c}")
        m0 = float(R.loc[lab[(lab.parent == p) & (lab.variant == "M0")].sweep_id.iloc[0], "SLC7A11"])
        m1 = v("RPMI", "M1", "SLC7A11")
        mine["SLC7A11 cysteine (M0) vs cystine (M1 RPMI)"] = (status(True, [m1], [m0], m1 < m0), f"M0 {m0} M1 {m1}")
        pur = ["PPAT", "GART", "PFAS", "PAICS", "ATIC"]
        for key, name in (("M-HPLM", "physiological: purines lose essentiality in HPLM (hypoxanthine) (M2 f=0.01)"),
                          ("M-Plasmax", "physiological: purines in Plasmax (hypoxanthine 5 uM) (M2 f=0.01)")):
            a = [v(key, "M2", g, 0.01) for g in pur]
            b = [v("RPMI", "M2", g, 0.01) for g in pur]
            mine[name] = (status(True, a, b, all(x > y for x, y in zip(a, b))), f"{a} vs RPMI {b}")
        a = [v("M-Plasmax", "M2", g, 0.01) for g in ("DHODH", "UMPS")]
        b = [v("RPMI", "M2", g, 0.01) for g in ("DHODH", "UMPS")]
        mine["physiological: DHODH/UMPS rescued by uridine, Plasmax only (M2 f=0.01)"] = (
            status(True, a, b, all(x > y for x, y in zip(a, b))), f"{a} vs {b}")
        a, b = v("M-Plasmax", "M2", "CAD", 0.01), v("RPMI", "M2", "CAD", 0.01)
        mine["physiological: CAD rescued by uridine, Plasmax only (M2 f=0.01)"] = (
            "model_structure_failure (W3/W5)" if a < 0.01 else "pass", f"Plasmax {a} RPMI {b}")
        gs = ["GOT1", "GOT2", "MDH2", "SDHA", "CS", "ATP5F1A"]
        a = [v("SENS:RPMI-ATCC", "M1", g) for g in gs]
        b = [v("RPMI", "M1", g) for g in gs]
        mine["physiological: pyruvate relaxes GOT/MDH/OXPHOS (RPMI-ATCC with pyruvate vs classic RPMI, M1)"] = (
            status(True, a, b, all(x >= y for x, y in zip(a, b)) and any(x > y for x, y in zip(a, b))), f"{a} vs {b}")
        docl = lab[(lab.parent == p) & lab.role.isin(["documented", "physiological"])]
        core = ["HMGCR", "MVK", "PMVK", "MVD"]
        n_loss = int((R.loc[docl.sweep_id.unique(), core] >= 0.01).any(axis=1).sum())
        mine["mevalonate core (HMGCR, MVK, PMVK, MVD) essential at every step"] = (
            "pass" if n_loss == 0 else "fail", f"{n_loss} of {docl.sweep_id.nunique()} sweeps lose a core gene")
        for ctrl, (st, obs) in mine.items():
            w2 = ctl[(ctl.parent == p) & (ctl.control == ctrl)]
            rows.append(dict(parent=p, control=ctrl, status_review=st,
                             status_w2=w2.status.iloc[0] if len(w2) else "MISSING", observed_review=obs,
                             agree=bool(len(w2)) and w2.status.iloc[0] == st))
        extra = set(ctl[ctl.parent == p].control) - set(mine)
        for c in extra:
            rows.append(dict(parent=p, control=c, status_review="not recomputed",
                             status_w2=ctl[(ctl.parent == p) & (ctl.control == c)].status.iloc[0], agree=None))
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print(df[["parent", "control", "status_review", "status_w2", "agree"]].to_string())
    print(df.groupby(["parent", "status_review"]).size())


if __name__ == "__main__":
    main()
