"""Re-check (A9): semantic consistency of results/w2/ko/label_changes_A9.tsv with DESIGN §17 A9 (label-free, no LP).

Expected under A9:
  (b) cysteine closed by default: DMEM:F12 defaults now point to the OLD 'SENS:DMEM:F12-cysteine-closed' sweep;
      the new 'SENS:DMEM:F12-cysteine-documented' points to the OLD DMEM:F12 default sweep; M-HPLM/M-Plasmax
      repointed (cysteine closed), their '|cys-documented' sensitivities added; RPMI+cys(2-ME) repointed (2-ME bug
      fix: cysteine now open);
  (d) classic RPMI: RPMI defaults now point to the OLD 'SENS:RPMI-no-pyruvate' sweep; the new 'SENS:RPMI-ATCC' points
      to the OLD RPMI default; RPMI texts that add pyruvate get their own keys (added);
  (h) the ATP key is removed and replaced by the adenosine key (added);
  unchanged: DMEM-based and IMDM-based keys and their sensitivities, M0, M-open.
Writes results/w2/review/label_check_A9.tsv (one row per check).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
L = ROOT / "results" / "w2" / "ko" / "label_changes_A9.tsv"
OUT = Path(__file__).resolve().parent / "label_check_A9.tsv"
K = ["parent", "variant", "f", "gate_open"]


def main() -> None:
    lc = pd.read_csv(L, sep="\t")
    lc["f"] = lc.f.fillna(-1.0)
    lc["medium_key"] = lc.medium_key.fillna("")
    rows = []

    def chk(name, ok, detail=""):
        rows.append(dict(check=name, ok=bool(ok), detail=detail))

    ct = lc.groupby("change").size().to_dict()
    chk("change_counts", ct == {"unchanged": 134, "repointed": 272, "label_added": 152, "label_removed": 60}, str(ct))
    chk("new_sweeps", lc[lc.new_sweep_existed_before_A9 == False].sweep_id_new.nunique() == 216,  # noqa: E712
        str(lc[lc.new_sweep_existed_before_A9 == False].sweep_id_new.nunique()))  # noqa: E712

    by = {c: set(g.medium_key) for c, g in lc.groupby("change")}
    exp_removed = {"SENS:RPMI-no-pyruvate", "SENS:DMEM:F12-cysteine-closed", "RPMI+ATP+adenine+hypoxanthine+thymidine"}
    chk("removed_keys", by.get("label_removed", set()) == exp_removed, str(sorted(by.get("label_removed", set()))))
    exp_added = {"SENS:RPMI-ATCC", "SENS:DMEM:F12-cysteine-documented", "SENS:DMEM:F12-no-standard-Gln",
                 "M-HPLM|cys-documented", "M-Plasmax|cys-documented", "RPMI+pyruvate", "RPMI+alanine+pyruvate",
                 "RPMI+adenine+adenosine+hypoxanthine+thymidine"}
    chk("added_keys", by.get("label_added", set()) == exp_added, str(sorted(by.get("label_added", set()))))
    unch = by.get("unchanged", set())
    chk("unchanged_keys_are_DMEM_IMDM_global",
        all(k == "" or k == "Hams_GR" or k.startswith(("DMEM+", "IMDM", "SENS:DMEM+", "SENS:DMEM-")) or k == "DMEM"
            for k in unch), str(sorted(unch)))
    rep = by.get("repointed", set())
    chk("repointed_keys_are_RPMI_or_DMEMF12_or_physiological",
        all(("RPMI" in k) or k.startswith("DMEM:F12") or k in ("M-HPLM", "M-Plasmax") for k in rep), str(sorted(rep)))

    # inheritance of old sweeps (bounds identical by construction)
    def mapping(key_new, key_old_src, col_new="sweep_id_new", col_old="sweep_id_old"):
        a = lc[lc.medium_key == key_new].set_index(K)[col_new]
        b = lc[lc.medium_key == key_old_src].set_index(K)[col_old]
        common = a.index.intersection(b.index)
        return len(common), int((a.loc[common] == b.loc[common]).sum())

    for new, old in (("RPMI", "SENS:RPMI-no-pyruvate"), ("SENS:RPMI-ATCC", "RPMI"),
                     ("DMEM:F12", "SENS:DMEM:F12-cysteine-closed"), ("SENS:DMEM:F12-cysteine-documented", "DMEM:F12")):
        n, eq = mapping(new, old)
        chk(f"inherits:{new}<-old {old}", n > 0 and n == eq, f"{eq}/{n} labels point to the old sweep")
    # physiological documented-cysteine sensitivities inherit the old physiological sweeps
    for ph in ("M-HPLM", "M-Plasmax"):
        n, eq = mapping(f"{ph}|cys-documented", ph)
        chk(f"inherits:{ph}|cys-documented<-old {ph}", n > 0 and n == eq, f"{eq}/{n}")
    # 2-ME: repointed to sweeps that did not exist before (cysteine now open)
    me = lc[lc.medium_key == "RPMI+cys(2-ME)"]
    chk("RPMI+cys(2-ME)_repointed_to_new_sweeps",
        (me.change == "repointed").all() and (me.new_sweep_existed_before_A9 == False).all(),  # noqa: E712
        f"{len(me)} labels; changes {me.change.value_counts().to_dict()}")
    # sweep ids hash the bounds, so DMEM:F12 old/new ids differ in both parents (cysteine exchange). The ratios should
    # reproduce the first review: Human1 identical, Human2 83-92 genes changed, 0 lethal flips.
    import numpy as np
    for parent in ("human1", "human2"):
        R = pd.read_parquet(ROOT / "results" / "w2" / "ko" / f"ratio_r_{parent}.parquet")
        d = lc[(lc.parent == parent) & (lc.medium_key == "DMEM:F12") & lc.gate_open & (lc.variant.isin(["M1", "M2"]))]
        nch, nflip = [], []
        for r in d.itertuples():
            a, b = R.loc[r.sweep_id_old], R.loc[r.sweep_id_new]
            nch.append(int((np.abs(a - b) > 1e-6).sum()))
            nflip.append(int(((a < 0.01) != (b < 0.01)).sum()))
        ok = (max(nch) == 0) if parent == "human1" else (max(nflip) == 0)
        chk(f"{parent}_DMEM:F12_old_vs_new_ratio_changes", ok, f"genes changed per label {nch}; lethal flips {nflip}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(df.to_string())


if __name__ == "__main__":
    main()
