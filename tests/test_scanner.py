"""§12 "Agent guards": the content scanner flags planted synthetic label copies (fingerprint and shape, in
many formats) and does not flag expression or CN matrices (toy and real); its output is coarse, so it
cannot be used as an oracle (review M6)."""

from __future__ import annotations

import os
import re

import numpy as np
import pandas as pd
import pytest

import toy


def _wide(df: pd.DataFrame, bare: bool = False) -> pd.DataFrame:
    out = df.copy()
    out.columns = [str(e) if bare else f"G{e - 1000} ({e})" for e in out.columns]
    out.index.name = "ModelID"
    return out


@pytest.fixture()
def planted(built_toy, tmp_path):
    t = built_toy
    d = tmp_path / "plant"
    d.mkdir()
    rng = np.random.default_rng(7)
    U = toy.U_ENTREZ
    ge, dp = t.ge[U], t.dp[U]
    held = t.ids_in("validation") + t.ids_in("test")
    # --- positives
    _wide(ge.iloc[:20]).to_csv(d / "copy_wide.csv")  # fingerprint + shape
    _wide(ge.iloc[:, :30]).T.reset_index().rename(columns={"index": "gene"}).to_csv(d / "copy_transposed.tsv", sep="\t", index=False)
    long = dp.stack(future_stack=True).rename("p").reset_index().rename(columns={"level_0": "ModelID", "level_1": "entrez"})
    long.to_parquet(d / "copy_long_prob.parquet")
    _wide(ge.loc[held]).to_pickle(d / "copy_heldout.pkl")
    _wide(ge.loc[held]).to_json(d / "copy_heldout.json")
    _wide(ge.loc[held]).to_csv(d / "copy_heldout.tsv.bz2", sep="\t")
    np.save(d / "copy_array.npy", ge.to_numpy())  # unkeyed array: shape rule
    bare = _wide(ge.loc[held], bare=True)
    bare.index = [int(m[4:]) for m in bare.index]  # ModelIDs without prefix and zero padding
    bare.index.name = "line"
    (-bare).to_csv(d / "copy_unprefixed_negated.csv")  # negated values, bare Entrez columns
    pad = pd.DataFrame(rng.normal(0, 0.3, (3000, 70)), index=[f"ACH-9{i:05d}" for i in range(3000)], columns=U)
    _wide(pd.concat([pad, ge.loc[held]])).to_csv(d / "copy_padded.csv")  # held-out rows after 3000 padding rows
    binm = (ge.loc[t.ids_in("test")] < -0.5).astype(int)
    _wide(binm).to_csv(d / "binary_test_keyed.csv")
    fake = pd.DataFrame(rng.normal(-0.4, 0.5, (12, 60)), index=t.model_ids[:12], columns=range(1000, 1060))
    _wide(fake).to_csv(d / "chronos_like_synthetic.csv")
    pd.DataFrame({"a": [1]}).to_csv(d / (toy.RAW_GE[:-4] + "_copy.csv"))  # protected (raw) name, wrong place
    # --- negatives
    expr = pd.DataFrame(rng.gamma(2, 1.5, (30, 70)), index=t.model_ids, columns=U)
    _wide(expr).to_csv(d / "expression_like.csv")
    cn = pd.DataFrame(rng.normal(1, 0.15, (30, 70)).clip(0), index=t.model_ids, columns=U)
    _wide(cn).to_csv(d / "cn_like.csv")
    _wide((ge.loc[t.ids_in("train")] < -0.5).astype(int)).to_csv(d / "binary_train_keyed.csv")
    _wide(ge + 0.37).to_csv(d / "shifted_chronos.csv")
    np.save(d / "positive_array.npy", rng.gamma(2, 1, (30, 70)))
    return d


def _scan(t, *args):
    r = t.run("scan_label_copies.py", *args)
    hits = {}
    for ln in r.stdout.splitlines():
        if ln.startswith("HIT\t"):
            _, p, reason = ln.split("\t", 2)
            hits[os.path.basename(p)] = reason
    return r, hits


def test_scanner_flags_planted_copies(built_toy, planted):
    r, hits = _scan(built_toy, "--roots", str(planted))
    assert r.returncode == 1
    fp = "fingerprint: matches label values"
    for f in ("copy_wide.csv", "copy_transposed.tsv", "copy_long_prob.parquet", "copy_heldout.pkl", "copy_heldout.json",
              "copy_heldout.tsv.bz2", "copy_unprefixed_negated.csv", "copy_padded.csv"):
        assert fp in hits.get(f, ""), (f, hits.get(f))
    assert "shape: Chronos-like" in hits["copy_wide.csv"]
    assert "shape: unkeyed Chronos-like array" in hits["copy_array.npy"]
    assert "shape: binary matrix keyed by held-out" in hits["binary_test_keyed.csv"]
    assert "shape: Chronos-like" in hits["chronos_like_synthetic.csv"] and "fingerprint" not in hits["chronos_like_synthetic.csv"]
    assert "name:" in hits[toy.RAW_GE[:-4] + "_copy.csv"]
    for f in ("expression_like.csv", "cn_like.csv", "binary_train_keyed.csv", "positive_array.npy"):
        assert f not in hits, (f, hits.get(f))
    assert "fingerprint" not in hits.get("shifted_chronos.csv", "")


def test_scanner_output_is_coarse_no_oracle(built_toy, tmp_path):
    """Reviewer's oracle: 30 known training cells + 1 held-out cell, with 60 candidate values in separate
    columns. The output must not depend on which column holds the true held-out value, and carries no counts."""
    t = built_toy
    tr = t.ids_in("train")
    held_mid, held_e = t.ids_in("test")[0], 1050
    true_v = round(float(t.ge.loc[held_mid, held_e]), 2)
    known = [(tr[i % len(tr)], 1040 + i, round(float(t.ge.loc[tr[i % len(tr)], 1040 + i]), 2)) for i in range(30)]
    outs = []
    for k_true in (3, 41):
        cand = np.round(np.linspace(-1.5, 0.5, 60), 2)
        cand[k_true] = true_v
        rows = [(m, e, *([v] * 60)) for m, e, v in known] + [(held_mid, held_e, *cand)]
        d = tmp_path / f"o{k_true}"
        d.mkdir()
        pd.DataFrame(rows, columns=["ModelID", "entrez", *[f"c{i}" for i in range(60)]]).to_csv(d / "probe.csv", index=False)
        r, hits = _scan(t, "--roots", str(d))
        outs.append(hits.get("probe.csv"))
        reasons = [ln.split("\t", 2)[2] for ln in r.stdout.splitlines() if ln.startswith("HIT\t")]
        assert reasons and not any(re.search(r"\d|cells|column", x) for x in reasons)  # no counts, no columns
    assert outs[0] == outs[1]


def test_scanner_delete_scratch_hits(built_toy, planted):
    r, hits = _scan(built_toy, "--roots", str(planted), "--delete-scratch-hits", "--scratch-glob", str(planted))
    assert not (planted / "copy_wide.csv").exists() and (planted / "expression_like.csv").exists()
    assert "DELETED" in r.stdout
    log = (built_toy.state / "logs" / "label_copy_deletions.tsv").read_text()
    assert "copy_wide.csv" in log


def test_scanner_negative_on_real_expression_and_cn(built_toy):
    """Real DepMap expression and CN matrices are not flagged (toy fingerprints loaded)."""
    dm = toy.PROJECT / "data" / "depmap"
    r, hits = _scan(built_toy, "--roots", str(dm / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
                    str(dm / "OmicsCNGeneWGS.csv"))
    assert r.returncode == 0 and hits == {}, hits


def test_fingerprints_unreadable_outside_scanner(built_toy, monkeypatch):
    from nextgem import labels

    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(built_toy.data))
    with pytest.raises(PermissionError):
        labels.load_fingerprints()


def test_canonical_split_exempt_but_copies_and_test_keyed_results_flagged(tmp_path_factory):
    """The canonical split (pinned by SPLIT_SHA256 and split_meta.json) is not a test-keyed hit; a copy
    elsewhere, a modified version in place, and a per-line results table keyed by test IDs still are."""
    import hashlib
    import json
    import shutil

    t = toy.make_toy(tmp_path_factory.mktemp("canon"))
    sp = t.state / "results" / "splits"
    split = pd.read_csv(sp / "lineage_split.tsv", sep="\t")
    split["n_screens"] = [1 + i % 3 for i in range(len(split))]  # integer covariate (parsed as a gene column)
    split["has_expression"] = True
    split.to_csv(sp / "lineage_split.tsv", sep="\t", index=False)
    sha = hashlib.sha256((sp / "lineage_split.tsv").read_bytes()).hexdigest()
    (sp / "SPLIT_SHA256").write_text(f"{sha}  lineage_split.tsv\n")
    (sp / "split_meta.json").write_text(json.dumps({"split_sha256": sha}))
    other = t.state / "results" / "other"
    other.mkdir()
    shutil.copy(sp / "lineage_split.tsv", other / "split_copy.tsv")
    test_ids = t.ids_in("test")
    pd.DataFrame({"ModelID": test_ids, "n_called": [3, 5, 2, 7, 4], "mcc": [0.1, 0.2, 0.3, 0.4, 0.5]}).to_csv(
        other / "per_line_results.tsv", sep="\t", index=False)

    def scan_hits():
        _, hits = _scan(t, "--test-keyed", "--roots", str(t.state / "results"))
        return hits

    hits = scan_hits()
    assert "lineage_split.tsv" not in hits, hits.get("lineage_split.tsv")
    assert "test-keyed" in hits.get("split_copy.tsv", "")
    assert "test-keyed" in hits.get("per_line_results.tsv", "")
    # modified in place -> no longer matches the pin -> flagged
    with open(sp / "lineage_split.tsv", "a") as f:
        f.write(f"{test_ids[0]}\ttest\tS3\tLinT\tAvana-containing\t9\tTrue\n")
    assert "test-keyed" in scan_hits().get("lineage_split.tsv", "")
    # pins that disagree -> no exemption
    shutil.copy(other / "split_copy.tsv", sp / "lineage_split.tsv")
    (sp / "split_meta.json").write_text(json.dumps({"split_sha256": "0" * 64}))
    assert "test-keyed" in scan_hits().get("lineage_split.tsv", "")


def test_small_integer_columns_are_not_bare_modelids():
    """Guard changelog 2026-10-27: a 2-row table whose integer 'group'/'rank' columns happen to equal known
    bare ModelIDs must not be read as ModelID-keyed; a long unique integer column still is."""
    import pandas as pd

    from nextgem import scan

    known = {f"ACH-{i:06d}" for i in range(1, 2000)}
    assert not scan._is_mid_col(pd.Series([76, 612]), known)
    assert scan._is_mid_col(pd.Series(list(range(100, 140))), known)
