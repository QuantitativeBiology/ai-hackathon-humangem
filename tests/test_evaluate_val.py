"""scripts/evaluate_val.py (W1 look) on synthetic toy data only. No real label file is touched: the validation
loader is a synthetic function and NEXTGEM_STATE_ROOT points to a toy git repository.

Covered: label-free preparation and the config JSON happen before the ledger ``started`` row; ``started`` is
committed before validation labels are read; scalar-only outputs and the hardened guard; a look id is never
repeated; the cap refuses before labels are read; a crashed look resumes with a ledger marker and is counted
once by the cap; the G-star sign rule (ties to Human2); both B6 branches; the pre-flight FROZEN/git check."""

from __future__ import annotations

import importlib.util
import json
import subprocess

import numpy as np
import pandas as pd
import pytest

import toy

from nextgem import paths

SPEC = importlib.util.spec_from_file_location("evaluate_val", paths.SCRIPTS_DIR / "evaluate_val.py")
ev = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ev)

GENES = list(range(1001, 1061))  # 60 toy genes
VAL_LINEAGES = {"VA": 6, "VB": 6, "VC": 5, "VD": 5}
TRAIN_LINEAGES = {"TA": 8, "TB": 8, "TC": 8}


def _ids(spec, start):
    out, lin = [], {}
    for name, n in spec.items():
        for _ in range(n):
            mid = f"ACH-{start + len(out):06d}"
            out.append(mid)
            lin[mid] = name
    return out, lin


def _long(lines, lineage, rng, essential):
    rows = []
    for l in lines:
        ge = rng.normal(0, 0.3, len(GENES)) - 1.2 * essential + rng.normal(0, 0.2, len(GENES))
        for g, v in zip(GENES, ge):
            dp = 1 / (1 + np.exp(4 * (v + 0.5)))
            rows.append((l, g, v, dp))
    df = pd.DataFrame(rows, columns=["ModelID", "entrez", "gene_effect", "dep_prob"])
    df["L_prob"] = (df.dep_prob > 0.5).astype(float)
    df["L_ge05"] = (df.gene_effect < -0.5).astype(float)
    df["L_ge06"] = (df.gene_effect < -0.6).astype(float)
    df["L_ge10"] = (df.gene_effect < -1.0).astype(float)
    df["lineage"] = df.ModelID.map(lineage)
    return df


@pytest.fixture()
def setup(tmp_path, monkeypatch):
    repo = tmp_path / "state"
    toy.init_repo(repo)
    (repo / "README").write_text("toy\n")
    toy.commit_all(repo, "init")
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(repo))
    rng = np.random.default_rng(5)
    essential = (rng.random(len(GENES)) < 0.25).astype(float)
    tr, lin_t = _ids(TRAIN_LINEAGES, 1)
    va, lin_v = _ids(VAL_LINEAGES, 500)
    lineage = pd.Series({**lin_t, **lin_v})
    train_long = _long(tr, lineage, rng, essential)
    val_long = _long(va, lineage, rng, essential).drop(columns="lineage")
    ratio = np.where((essential > 0) & (rng.random(len(GENES)) < 0.6), 0.0, 1.0)
    gen = pd.DataFrame(index=pd.Index(GENES, name="entrez"))
    for s in ("H1", "H2"):
        r = ratio if s == "H2" else np.where(rng.random(len(GENES)) < 0.1, 1 - ratio, ratio)
        gen[f"in_{s}"] = True
        gen[f"ratio_{s}"] = r
        gen[f"e_{s}"] = 1 - r
        gen[f"fba_call_{s}"] = r < 0.01
        gen[f"p57_{s}"] = r < 0.01
        gen[f"pgr_{s}"] = r < 0.01
        gen[f"et_{s}"] = 1 - r
    all_lines = tr + va
    expr = pd.DataFrame(rng.gamma(2, 1.5, (len(all_lines), len(GENES))), index=all_lines, columns=GENES).astype("float32")
    cn = pd.DataFrame(rng.normal(1, 0.2, (len(all_lines), len(GENES))), index=all_lines, columns=GENES).astype("float32")
    lf = pd.DataFrame({"base_medium_class": "sheet:RPMI", "library_class": "Avana-containing", "model_type": "Cell Line",
                       "has_expression_matrix": True, "serum": 1.0}, index=pd.Index(all_lines, name="ModelID"))
    ctx = {"genes": GENES, "ucap": GENES[:50], "gen": gen, "omics": {"expr": expr, "cn": cn}, "lf": lf, "lineage": lineage,
           "line_class": lf.base_medium_class, "val_lines": va, "val_expr_lines": va, "all_ids": all_lines,
           "depmap_ce": set(GENES[:3]), "strata": {"none": ("none", None), "expressed_pairs": ("pairs", "expressed")},
           "b6": None, "gene_symbols": {"GENEA", "GENEB"}}
    cfg = {"look_id": "TOY_W1", "workstream": "W1", "milestone": "M2", "labels": ["L_prob", "L_ge05", "L_ge10"],
           "binary_label": "L_ge06", "genesets": ["all", "noCE", "noDepMapCE"], "universes": ["U", "Ucap"],
           "strata": ["none", "expressed_pairs"], "B": 200, "seed": 11, "random_seeds": [11],
           "gstar_rule": "toy", "inputs": {}, "code": {}, "include_b6": False}
    calls = {"val": 0, "crash": 0}

    def load_val():
        # the started row must already be committed when validation labels are first read
        head = subprocess.run(["git", "-C", str(repo), "show", "HEAD:results/validation_ledger.tsv"],
                              capture_output=True, text=True).stdout
        rows = [ln.split("\t") for ln in head.splitlines()[1:]]
        assert any(r[3] == "started" and r[4] == cfg["look_id"] for r in rows), "labels read before 'started' was committed"
        calls["val"] += 1
        if calls["crash"]:
            calls["crash"] -= 1
            raise RuntimeError("simulated crash after the label read")
        return val_long

    return {"repo": repo, "cfg": cfg, "ctx": ctx, "train": train_long, "load_val": load_val, "calls": calls,
            "out": tmp_path / "out"}


def _run(s, include_b6=False, b6=None, write_cfg=True):
    cfg = dict(s["cfg"], include_b6=include_b6)
    ctx = dict(s["ctx"], b6=b6)
    prep = ev.prepare_look(cfg, ctx, s["train"])
    if write_cfg:
        ev.write_config(cfg, s["out"])
    return ev.run_look(cfg, prep=prep, load_val=s["load_val"], ctx=ctx, out_dir=s["out"], author="test")


def _b6_stub(calls):
    def b6(lines, fitted):
        calls["b6"] = calls.get("b6", 0) + 1
        rng = np.random.default_rng(3)
        return pd.DataFrame(np.round(rng.random((len(lines), len(GENES))), 8), index=lines, columns=GENES)
    return b6


def test_started_committed_before_labels_and_completed(setup):
    from nextgem import ledger

    res = _run(setup)
    assert setup["calls"]["val"] == 1
    rows = ledger.read("validation")
    assert [r["event"] for r in rows] == ["started", "output", "completed"]
    assert rows[1]["output_sha"] == res["output_sha"]
    assert rows[0]["config_sha"] == ev.config_sha(dict(setup["cfg"], include_b6=False))
    assert not ledger.verify("validation")


def test_config_json_required_before_started(setup):
    from nextgem import ledger

    with pytest.raises(RuntimeError, match="config JSON"):
        _run(setup, write_cfg=False)
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0
    cfg = dict(setup["cfg"], include_b6=False)
    doc = json.loads(ev.write_config(cfg, setup["out"]).read_text())
    assert doc["config_sha"] == ev.config_sha(cfg)


def test_label_free_failure_happens_before_started(setup):
    from nextgem import ledger

    def broken(lines, fitted):
        raise ValueError("label-free B6 failure")

    with pytest.raises(ValueError):
        _run(setup, include_b6=True, b6=broken)
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0


def test_outputs_are_scalar_only(setup):
    res = _run(setup)
    out = pd.read_csv(res["output"], sep="\t")
    ids = set(setup["ctx"]["all_ids"])
    for c in out.columns:
        assert not out[c].astype(str).isin(ids).any(), c
    assert not ({"ModelID", "entrez", "lineage", "gene"} & set(out.columns))
    assert set(out.section) == {"tAP_grid", "binary_6.8"}
    files = sorted(p.name for p in setup["out"].iterdir())
    assert files == ["val_look_TOY_W1.tsv", "val_look_TOY_W1_config.json", "val_look_TOY_W1_gstar.json"]


def test_b6_excluded_by_default(setup):
    calls = {}
    res = _run(setup, include_b6=False, b6=_b6_stub(calls))
    out = pd.read_csv(res["output"], sep="\t")
    assert not out.method.astype(str).str.startswith("B6").any()
    assert calls.get("b6", 0) == 0


def test_b6_branch_when_enabled(setup):
    calls = {}
    res = _run(setup, include_b6=True, b6=_b6_stub(calls))
    out = pd.read_csv(res["output"], sep="\t")
    b6 = out[out.method == "B6[L_prob]"]
    assert calls["b6"] == 1 and len(b6) == 3 * 2 * 2 and set(b6.label) == {"L_prob"}


def test_gstar_sign_rule():
    def grid(d):
        rows = []
        for m, dl in (("H1", np.nan), ("H2", d)):
            rows.append({"label": "L_prob", "geneset": "all", "universe": "U", "stratum": "none", "method": m,
                         "tap_lin": 0.3, "delta_lin_vs_H1": dl, "delta_lo_vs_H1": dl, "delta_hi_vs_H1": dl})
        return pd.DataFrame(rows)

    assert ev.decide_gstar(grid(-0.001))["G_star"] == "Human1"
    assert ev.decide_gstar(grid(0.0))["G_star"] == "Human2"  # tie -> Human2 (default parent)
    assert ev.decide_gstar(grid(0.002))["G_star"] == "Human2"


def test_gstar_matches_sign_in_look(setup):
    g = _run(setup)["gstar"]
    assert g["G_star"] == ("Human1" if g["delta_lin_H2_vs_H1"] < 0 else "Human2")


def test_repeat_look_refused_without_reading_labels(setup):
    from nextgem import ledger

    _run(setup)
    with pytest.raises(ledger.LedgerError, match="never repeated"):
        _run(setup)
    assert setup["calls"]["val"] == 1


def test_cap_refuses_before_labels(setup):
    from nextgem import ledger

    for i in range(40):
        ledger.append("validation", "started", look_id=f"v{i}", workstream="W2", config_sha="c", commit=False)
    toy.commit_all(setup["repo"], "40 looks")
    with pytest.raises(ledger.LedgerError, match="cap reached"):
        _run(setup)
    assert setup["calls"]["val"] == 0


def test_resumed_look_is_logged_and_counted_once(setup):
    from nextgem import ledger

    for i in range(39):
        ledger.append("validation", "started", look_id=f"v{i}", workstream="W2", config_sha="c", commit=False)
    toy.commit_all(setup["repo"], "39 looks")
    setup["calls"]["crash"] = 1
    with pytest.raises(RuntimeError, match="simulated crash"):
        _run(setup)
    rows = ledger.read("validation")
    assert [r["event"] for r in rows if r["look_id"] == "TOY_W1"] == ["started"]
    res = _run(setup)  # resume: same config, 40 looks in total -> allowed (the crashed look is already counted)
    assert res["resumed"]
    ev_rows = [r for r in ledger.read("validation") if r["look_id"] == "TOY_W1"]
    assert [r["event"] for r in ev_rows] == ["started", "output", "output", "completed"]
    assert json.loads(ev_rows[1]["detail"]) == {"event": "resumed", "attempt": 1}
    assert len(ledger.looks(ledger.read("validation"))) == 40
    assert not ledger.verify("validation")
    with pytest.raises(ledger.LedgerError, match="cap reached"):
        ledger.append("validation", "started", look_id="new", workstream="W2", config_sha="c", commit=False)


# ---- hardened scalar guard
def _valid_out(setup):
    res = _run(setup)
    return pd.read_csv(res["output"], sep="\t")


def _guard(setup, out):
    cfg = dict(setup["cfg"], include_b6=False)
    n = {lab: len([m for m in out[out.section == "tAP_grid"].query("label == @lab").method.unique()]) for lab in cfg["labels"]}
    ev.assert_scalar_only(out, cfg, setup["ctx"], n)


def test_guard_accepts_valid_output(setup):
    _guard(setup, _valid_out(setup))


def test_guard_rejects_duplicate_family_rows(setup):
    out = _valid_out(setup)
    with pytest.raises(RuntimeError, match="more than one row|rows >"):
        _guard(setup, pd.concat([out, out.iloc[[0]]], ignore_index=True))


def test_guard_rejects_too_many_rows(setup):
    out = _valid_out(setup)
    extra = out[out.section == "binary_6.8"].copy()
    extra["metric"] = extra["metric"] + "_x"
    with pytest.raises(RuntimeError, match="rows > expected"):
        _guard(setup, pd.concat([out, extra], ignore_index=True))


@pytest.mark.parametrize("value,why", [("ACH-000501", "line ID"), ("VA", "lineage name"), ("GENEA", "gene"), ("1005", "gene")])
def test_guard_rejects_unit_keys_in_cells(setup, value, why):
    out = _valid_out(setup)
    out.loc[out.index[0], "stratum"] = value
    with pytest.raises(RuntimeError, match=why):
        _guard(setup, out)


def test_guard_rejects_unit_key_columns(setup):
    out = _valid_out(setup)
    with pytest.raises(RuntimeError, match="line ID"):
        _guard(setup, out.assign(**{"ACH-000501": 0.1}))
    with pytest.raises(RuntimeError, match="gene"):
        _guard(setup, out.assign(**{"1001": 0.1}))


def test_guard_rejects_numeric_per_unit_columns(setup):
    out = _valid_out(setup)
    with pytest.raises(RuntimeError, match="not a known scalar"):
        _guard(setup, out.assign(tap_line_17=0.2))


# ---- pre-flight
def test_preflight_frozen_and_clean(tmp_path, monkeypatch):
    repo = tmp_path / "pf"
    toy.init_repo(repo)
    (repo / ".gitignore").write_text("ign/\n")
    (repo / "a.txt").write_text("a\n")
    (repo / "ign").mkdir()
    (repo / "ign" / "b.bin").write_text("b\n")
    (repo / "untracked.txt").write_text("u\n")
    sha = {f: ev._sha(repo / f) for f in ("a.txt", "ign/b.bin", "untracked.txt")}
    hdr = "item\tkind\tpath\tsha256_or_value\tfrozen_on\tamendment\tnote\n"
    (repo / "FROZEN.tsv").write_text(hdr + "".join(f"{f}\tfile\t{f}\t{sha[f]}\t2026-10-26\t\t\n" for f in sha))
    toy.git(repo, "add", "a.txt", ".gitignore", "FROZEN.tsv")
    toy.git(repo, "commit", "-q", "-m", "init")
    assert ev.preflight(["a.txt", "ign/b.bin"], root=repo) == []
    assert any("not tracked" in p for p in ev.preflight(["untracked.txt"], root=repo))
    (repo / "a.txt").write_text("changed\n")
    pr = ev.preflight(["a.txt"], root=repo)
    assert any("differs" in p for p in pr) and any("uncommitted" in p for p in pr)
    (repo / "ign" / "b.bin").write_text("changed\n")
    assert any("differs" in p for p in ev.preflight(["ign/b.bin"], root=repo))
    (repo / "c.txt").write_text("c\n")
    toy.git(repo, "add", "c.txt")
    toy.git(repo, "commit", "-q", "-m", "c")
    assert any("no file row" in p for p in ev.preflight(["c.txt"], root=repo))
