"""scripts/evaluate_val.py shortlist looks (a reviewed look spec, no harness change) on synthetic toy data only.
No real label file is touched: validation labels are a synthetic function, NEXTGEM_STATE_ROOT points to a toy git
repository, and the W1-config check reads only the committed W1 config record (no input is hashed).

Covered: the W1 look's config is unchanged apart from the harness's own SHA; look ids are registered from spec
files; a wrong scores SHA, a non-score or label-like scores column and bad spec fields are refused before any
label access; the reference must be the W1 G-star; generic H1/H2 in the scores file must equal the W1 look's;
a toy end-to-end look gives known-answer Δ values, a scalar-only output and the usual ledger rows; the §8 gate and
the §9 paired-difference one-SE rule and winner's-curse term; the dry run reads no validation label and writes no ledger row."""

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

LOOK = "TOY_W9_M4"
SPEC_REL = f"results/w9/val/look_spec_{LOOK}.json"
SCORES_REL = "results/w9/val/w9_val_scores.parquet"
GENES = list(range(1001, 1041))  # 40 toy genes
POS = GENES[:8]  # positive in every validation line
VAL_LINEAGES = {"VA": 5, "VB": 5, "VC": 5, "VD": 5}
TRAIN_LINEAGES = {"TA": 6, "TB": 6, "TC": 6}
PERFECT, SAME = "W9:perfect[H1]", "W9:same[H1]"
W1_RECORD = "results/w1/val_look_W1_M2_baseline_config.json"


def _ids(spec, start):
    out, lin = [], {}
    for name, n in spec.items():
        for _ in range(n):
            mid = f"ACH-{start + len(out):06d}"
            out.append(mid)
            lin[mid] = name
    return out, lin


def _ap_worst(n, p):
    """AP when all p positives are ranked below all n - p negatives (distinct scores)."""
    return sum(j / (n - p + j) for j in range(1, p + 1)) / p


def _ap_constant(n, p):
    """tAP of a constant scorer = E[AP_random] (DESIGN §6.1 single-block form)."""
    return sum((1 + (r - 1) * (p - 1) / (n - 1)) / r for r in range(1, n + 1)) / n


def _write_scores(root, df):
    p = root / SCORES_REL
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(p, index=False)
    return ev._sha(p)


def _write_spec(root, **over):
    spec = {"look_id": LOOK, "kind": "shortlist", "workstream": "W9", "milestone": "M4",
            "scores": {"path": SCORES_REL, "sha256": over.pop("sha256"), "layout": over.pop("layout", "gene")},
            "methods": ["H1", "H2", SAME, PERFECT], "candidates": [SAME, PERFECT], "reference": "H1",
            "labels": ["L_prob"], "genesets": ["all", "noCE"], "universes": ["U"], "strata": ["none", "expressed_pairs"],
            "B": 2000, "seed": 11, "decision_rule": "toy: §8 gate vs G*, then the §9 paired one-SE rule"}
    spec.update(over)
    p = root / SPEC_REL
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(spec, indent=1))
    return spec


@pytest.fixture()
def setup(tmp_path, monkeypatch):
    repo = tmp_path / "state"
    toy.init_repo(repo)
    for rel in ev.input_files().values():  # stand-ins for the W1 context inputs (hashed into the config only)
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(f"toy {rel}\n")
    (repo / ev.W1_GSTAR).write_text(json.dumps({"G_star": "Human1"}))
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(repo))
    rng = np.random.default_rng(7)
    tr, lin_t = _ids(TRAIN_LINEAGES, 1)
    va, lin_v = _ids(VAL_LINEAGES, 500)
    lineage = pd.Series({**lin_t, **lin_v})
    train = pd.DataFrame([(l, g) for l in tr for g in GENES], columns=["ModelID", "entrez"])
    train["gene_effect"] = rng.normal(0, 0.5, len(train))
    train["dep_prob"] = rng.random(len(train))
    for lab in ("L_prob", "L_ge05", "L_ge06", "L_ge10"):
        train[lab] = (rng.random(len(train)) < 0.2).astype(float)  # no gene reaches CE_fit (>= 90%)
    val_long = pd.DataFrame([(l, g, float(g in POS)) for l in va for g in GENES], columns=["ModelID", "entrez", "L_prob"])
    gen = pd.DataFrame(index=pd.Index(GENES, name="entrez"))
    # H1 ranks every positive below every negative (distinct scores); H2 is constant (one tied block)
    gen["e_H1"] = [0.01 * i if g in POS else 0.5 + 0.01 * i for i, g in enumerate(GENES)]
    gen["e_H2"] = 0.5
    all_lines = tr + va
    expr = pd.DataFrame(5.0, index=all_lines, columns=GENES).astype("float32")  # every pair expressed
    lf = pd.DataFrame({"base_medium_class": "sheet:RPMI", "has_expression_matrix": True},
                      index=pd.Index(all_lines, name="ModelID"))
    ctx = {"genes": GENES, "ucap": GENES, "gen": gen, "omics": {"expr": expr, "cn": expr}, "lf": lf,
           "lineage": lineage, "line_class": lf.base_medium_class, "val_lines": va, "val_expr_lines": va,
           "all_ids": all_lines, "depmap_ce": set(), "gene_symbols": {"GENEA"},
           "strata": {"none": ("none", None), "expressed_pairs": ("pairs", "expressed")}, "b6": None}
    scores = pd.DataFrame({"entrez": GENES, "H1": gen.e_H1.to_numpy(), "H2": gen.e_H2.to_numpy(), SAME: gen.e_H1.to_numpy(),
                           PERFECT: [0.9 if g in POS else 0.1 for g in GENES]})
    _write_spec(repo, sha256=_write_scores(repo, scores))
    toy.commit_all(repo, "init")
    calls = {"val": 0}

    def load_val():
        head = subprocess.run(["git", "-C", str(repo), "show", "HEAD:results/validation_ledger.tsv"],
                              capture_output=True, text=True).stdout
        rows = [ln.split("\t") for ln in head.splitlines()[1:]]
        assert any(r[3] == "started" and r[4] == LOOK for r in rows), "labels read before 'started' was committed"
        calls["val"] += 1
        return val_long

    return {"repo": repo, "ctx": ctx, "train": train, "scores": scores, "load_val": load_val, "calls": calls,
            "va": va}


def _run(s):
    cfg = ev.shortlist_config(LOOK)
    prep = ev.prepare_shortlist(cfg, s["ctx"], s["train"])
    out_dir = ev.look_out_dir(cfg)
    ev.write_config(cfg, out_dir)
    return ev.run_look(cfg, prep=prep, load_val=s["load_val"], ctx=s["ctx"], out_dir=out_dir, author="test")


def _forbid_labels(monkeypatch):
    from nextgem import labels

    hit = []

    def forbidden(name):
        def f(*a, **k):
            hit.append(name)
            raise AssertionError(f"labels.{name}() called")
        return f

    for name in ("train", "validation", "test", "load_split", "label_derived_list", "gene_universe"):
        monkeypatch.setattr(labels, name, forbidden(name))
    return hit


# ---- the W1 look is unchanged
def test_w1_config_unchanged(monkeypatch):
    """Every field of the W1 config except the harness's own SHA equals the committed record of the W1 look,
    and that record's SHA is reproduced by config_sha (so only code['evaluate_val.py'] moves)."""
    rec = json.loads((paths.PROJECT_ROOT / W1_RECORD).read_text())
    assert rec["config_sha"] == "331a2c8d1b25c509fb12ac6e1771fb792b9d48f13426f698ebd092b0801a60e5"
    monkeypatch.setattr(ev, "_sha", lambda p: "sha")  # no input is read
    cfg = ev.w1_config()
    keep = lambda d: {k: v for k, v in d.items() if k not in ("inputs", "code", "config_sha") and not k.startswith("_")}  # noqa: E731
    assert keep(cfg) == keep(rec)
    assert set(cfg["inputs"]) == set(rec["inputs"]) and set(cfg["code"]) == set(rec["code"])
    assert ev.config_sha({**cfg, "inputs": rec["inputs"], "code": rec["code"]}) == rec["config_sha"]
    assert ev.look_config(ev.W1_LOOK)["look_id"] == ev.W1_LOOK and "kind" not in cfg
    assert ev.look_out_dir(cfg) == paths.results_dir() / "w1"


def test_look_ids_registered_from_specs(setup, capsys):
    assert list(ev.look_specs()) == [LOOK]
    assert ev.main(["--look", LOOK, "--print-config"]) == 0
    doc = json.loads(capsys.readouterr().out)
    assert doc["config_sha"] == ev.config_sha(ev.shortlist_config(LOOK)) and doc["kind"] == "shortlist"
    assert doc["inputs"]["scores"] == ev._sha(setup["repo"] / SCORES_REL) and doc["G_star"] == "Human1"
    with pytest.raises(SystemExit):
        ev.main(["--look", "W9_unknown", "--print-config"])


# ---- refusals before any label access
def test_wrong_scores_sha_refused_before_label_access(setup, monkeypatch):
    _write_spec(setup["repo"], sha256="0" * 64)
    hit = _forbid_labels(monkeypatch)
    with pytest.raises(ev.PreflightError, match="SHA-256 differs"):
        ev.shortlist_config(LOOK)
    for mode in ("--print-config", "--preflight-only", "--dry-run"):
        assert ev.main(["--look", LOOK, mode]) == 2
    assert ev.main(["--look", LOOK]) == 2
    assert hit == [] and not (setup["repo"] / "results" / "validation_ledger.tsv").exists()


@pytest.mark.parametrize("change,match", [
    pytest.param(lambda d: d.assign(note="x"), "non-score column", id="string_note_column"),
    pytest.param(lambda d: d.assign(L_prob=0.0), "non-score column", id="label_column_L_prob"),
    pytest.param(lambda d: d.assign(dep_prob=0.3), "non-score column", id="label_column_dep_prob"),
    pytest.param(lambda d: d.assign(lineage="VA"), "non-score column", id="lineage_column"),
    pytest.param(lambda d: d.assign(ModelID="ACH-000500"), "non-score column", id="line_id_column"),
    pytest.param(lambda d: d.drop(columns=[PERFECT]), "no scores column", id="fewer_missing_candidate"),
    pytest.param(lambda d: d.drop(columns=["H1"]), "no scores column", id="fewer_missing_H1"),
    pytest.param(lambda d: d.drop(columns=["H1", "H2"]), "no scores column", id="fewer_missing_H1_H2"),
    pytest.param(lambda d: d[["entrez"]], "no scores column", id="fewer_keys_only"),
    pytest.param(lambda d: d.assign(**{PERFECT: d[PERFECT].astype(str)}), "finite floats", id="string_scores"),
    pytest.param(lambda d: d.assign(**{PERFECT: (d[PERFECT] > 0.5).astype(int)}), "finite floats", id="int_scores"),
    pytest.param(lambda d: d.assign(**{PERFECT: np.where(d.entrez == 1001, np.nan, d[PERFECT])}), "finite floats",
                 id="nan_score"),
    pytest.param(lambda d: pd.concat([d, d.iloc[[0]]]), "duplicate", id="duplicate_gene"),
])
def test_scores_non_score_or_label_like_columns_refused(setup, monkeypatch, change, match):
    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], change(setup["scores"])))  # SHA pin matches
    hit = _forbid_labels(monkeypatch)
    with pytest.raises(ev.PreflightError, match=match):
        ev.shortlist_config(LOOK)
    assert hit == []


def test_line_gene_scores_checked(setup):
    va = setup["va"]
    long = pd.DataFrame([(l, g) for l in va for g in GENES], columns=["ModelID", "entrez"])
    base = setup["scores"].set_index("entrez")
    for m in ("H1", "H2", SAME, PERFECT):
        long[m] = base[m].reindex(long.entrez).to_numpy()
    binary = long.assign(**{PERFECT: (long[PERFECT] > 0.5).astype(float)})
    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], binary), layout="line_gene")
    with pytest.raises(ev.PreflightError, match="label-shaped"):
        ev.shortlist_config(LOOK)
    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], long), layout="line_gene")
    cfg = ev.shortlist_config(LOOK)
    prep = ev.prepare_shortlist(cfg, setup["ctx"], setup["train"])
    assert np.array_equal(prep["S"][PERFECT].to_numpy(), np.tile(base[PERFECT].to_numpy(), (len(va), 1)))
    extra = pd.concat([long, long[long.ModelID == va[0]].assign(ModelID="ACH-000999")])  # a non-validation line
    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], extra), layout="line_gene")
    with pytest.raises(RuntimeError, match="exactly the validation expression lines"):
        ev.prepare_shortlist(ev.shortlist_config(LOOK), setup["ctx"], setup["train"])


@pytest.mark.parametrize("over,match", [
    pytest.param({"seed": 23}, "seed", id="seed_not_SEEDS0"),
    pytest.param({"B": 200}, "B must be", id="B_not_B_OTHER"),
    pytest.param({"candidates": ["H1"]}, "candidates", id="generic_candidate"),
    pytest.param({"candidates": [SAME, PERFECT, "H2", "x"]}, "candidates", id="four_candidates"),
    pytest.param({"methods": ["H1", "H2", SAME, PERFECT, "ACH-000001"]}, "refused", id="method_line_id"),
    pytest.param({"methods": ["H1", "H2", SAME, PERFECT, "1005"]}, "refused", id="method_entrez"),
    pytest.param({"methods": ["H1", "H2", SAME, PERFECT, "L_prob"]}, "refused", id="method_label_name"),
    pytest.param({"methods": [SAME, PERFECT]}, "H1 and H2", id="methods_without_generic"),
    pytest.param({"look_id": "OTHER"}, "file name", id="look_id_not_file_name"),
    pytest.param({"kind": "w1"}, "kind", id="kind_not_shortlist"),
    pytest.param({"strata": ["expressed_pairs"]}, "strata", id="strata_without_none"),
    pytest.param({"universes": ["U", "Uall"]}, "universes", id="unknown_universe"),
    pytest.param({"labels": ["L_ge06"]}, "labels", id="label_not_in_grid"),
    pytest.param({"extra_key": 1}, "keys", id="unknown_key"),
    pytest.param({"scores": {"path": "data/controls/x.parquet", "sha256": "0" * 64, "layout": "gene"}}, "scores.path",
                 id="scores_path_outside_results"),
    pytest.param({"scores": {"path": "results/w9/val/../../w1/x.parquet", "sha256": "0" * 64, "layout": "gene"}},
                 "scores.path", id="scores_path_dotdot"),
    pytest.param({"scores": {"path": "results/w3/val/x.parquet", "sha256": "0" * 64, "layout": "gene"}}, "scores.path",
                 id="scores_path_other_workstream"),
    pytest.param({"provenance": {"data/controls/x.tsv": "0" * 64}}, "provenance", id="provenance_outside_code_results"),
    pytest.param({"provenance": {"results/w9/x.json": "abc"}}, "provenance", id="provenance_bad_sha"),
])
def test_spec_fields_refused(setup, over, match):
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL), **over)
    with pytest.raises(ev.PreflightError, match=match):
        ev.shortlist_config(LOOK)


def test_reference_must_be_gstar(setup):
    sha = ev._sha(setup["repo"] / SCORES_REL)
    _write_spec(setup["repo"], sha256=sha, reference="H2")
    with pytest.raises(ev.PreflightError, match="not G-star"):
        ev.shortlist_config(LOOK)
    (setup["repo"] / ev.W1_GSTAR).write_text(json.dumps({"G_star": "Human2"}))
    assert ev.shortlist_config(LOOK)["reference"] == "H2"
    _write_spec(setup["repo"], sha256=sha, reference="H1")
    with pytest.raises(ev.PreflightError, match="not G-star"):
        ev.shortlist_config(LOOK)


def test_generic_scores_must_match_w1(setup):
    from nextgem import ledger

    bad = setup["scores"].assign(H1=setup["scores"].H1 + 0.001)
    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], bad))
    with pytest.raises(RuntimeError, match="differs from the generic H1"):
        ev.prepare_shortlist(ev.shortlist_config(LOOK), setup["ctx"], setup["train"])
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0


# ---- the look
def test_shortlist_known_answer_deltas(setup):
    from nextgem import ledger

    res = _run(setup)
    assert setup["calls"]["val"] == 1
    rows = ledger.read("validation")
    assert [r["event"] for r in rows] == ["started", "output", "completed"]
    assert rows[0]["config_sha"] == res["config_sha"] == ev.config_sha(ev.shortlist_config(LOOK))
    assert json.loads(rows[0]["detail"]) == {"kind": "shortlist", "milestone": "M4", "scalars_only": True}
    assert not ledger.verify("validation")
    out = pd.read_csv(res["output"], sep="\t").set_index(["geneset", "stratum", "method"]).sort_index()
    n, p = len(GENES), len(POS)
    worst, const = _ap_worst(n, p), _ap_constant(n, p)
    for gs in ("all", "noCE"):  # no CE_fit gene and every pair expressed: all four cells are identical
        for st in ("none", "expressed_pairs"):
            c = out.loc[(gs, st)]
            assert c.loc["H1", "tap_lin"] == pytest.approx(worst, abs=1e-6)
            assert c.loc["H2", "tap_lin"] == pytest.approx(const, abs=1e-6)
            assert c.loc[PERFECT, "tap_lin"] == pytest.approx(1.0, abs=1e-6)
            for k in ("delta_lin_vs_H1", "delta_lo_vs_H1", "delta_hi_vs_H1", "delta_line_vs_H1"):
                assert c.loc[PERFECT, k] == pytest.approx(1 - worst, abs=1e-6)
                assert c.loc[SAME, k] == pytest.approx(0.0, abs=1e-6)
            assert c.loc[PERFECT, "delta_lin_vs_H2"] == pytest.approx(1 - const, abs=1e-6)
            assert c.loc["H2", "delta_lin_vs_H1"] == pytest.approx(const - worst, abs=1e-6)
    dec = res["decision"]
    assert json.loads((setup["repo"] / "results/w9/val" / f"val_look_{LOOK}_decision.json").read_text()) == dec
    assert dec["gate_passed"] and dec["best_eligible"] == PERFECT and dec["selected"] == PERFECT
    assert [c["eligible"] for c in dec["candidates"]] == [False, True]
    assert dec["candidates"][1]["delta_lin"] == pytest.approx(1 - worst, abs=1e-12)
    wc = dec["winners_curse"]  # identical lines: σ = 0, so no optimism correction
    assert wc["method"] == PERFECT and wc["m"] == 2 and wc["term"] == pytest.approx(0.0, abs=1e-9)
    assert wc["best_delta_lin_minus_term"] == pytest.approx(1 - worst, abs=1e-9)


def test_shortlist_single_candidate(setup):
    """One gate candidate (m = 1): the other methods are descriptive only, the winner's-curse term is 0 and the
    gate/selection apply to the single candidate."""
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL), candidates=[PERFECT])
    res = _run(setup)
    dec, n, p = res["decision"], len(GENES), len(POS)
    assert [c["method"] for c in dec["candidates"]] == [PERFECT]
    assert dec["gate_passed"] and dec["best_eligible"] == dec["selected"] == PERFECT
    wc = dec["winners_curse"]
    assert wc["m"] == 1 and wc["term"] == 0.0 and wc["best_delta_lin_minus_term"] == wc["best_delta_lin"]
    assert wc["best_delta_lin"] == pytest.approx(1 - _ap_worst(n, p), abs=1e-9)
    out = pd.read_csv(res["output"], sep="\t")
    assert set(out.method) == {"H1", "H2", SAME, PERFECT}  # descriptive methods are still scored and reported


def test_shortlist_output_scalar_only(setup):
    res = _run(setup)
    out = pd.read_csv(res["output"], sep="\t")
    files = sorted(p.name for p in (setup["repo"] / "results/w9/val").iterdir() if p.name.startswith("val_look_"))
    assert files == [f"val_look_{LOOK}.tsv", f"val_look_{LOOK}_config.json", f"val_look_{LOOK}_decision.json"]
    assert set(out.section) == {"tAP_grid"} and len(out) == 4 * 1 * 2 * 1 * 2
    ids = set(setup["ctx"]["all_ids"]) | {str(g) for g in GENES} | {"VA", "VB", "VC", "VD"}
    for c in out.columns:
        assert str(c) not in ids and not out[c].astype(str).isin(ids).any(), c
    assert not ({"ModelID", "entrez", "lineage", "gene"} & set(out.columns))
    dec = json.dumps(res["decision"])
    assert not any(i in dec for i in ids if i.startswith("ACH-"))
    cfg, ctx = ev.shortlist_config(LOOK), setup["ctx"]
    n = {"L_prob": len(cfg["methods"])}
    ev.assert_scalar_only(out, cfg, ctx, n, n_bin=0)
    for bad, match in ((out.assign(tap_line_3=0.1), "not a known scalar"), (out.assign(**{"ACH-000500": 0.1}), "line ID"),
                       (pd.concat([out, out.iloc[[0]]], ignore_index=True), "rows > expected|more than one row"),
                       (pd.concat([out, out.iloc[[0]].assign(section="binary_6.8")], ignore_index=True), "rows > expected")):
        with pytest.raises(RuntimeError, match=match):
            ev.assert_scalar_only(bad, cfg, ctx, n, n_bin=0)


def _lt(taps, lines):
    return pd.DataFrame({"tap": taps, "n": 40, "n_pos": 8, "estimable": True}, index=lines)


def test_shortlist_one_se_rule(setup):
    """Δ̄(best − near) = 0.0025 with bootstrap SE ~ 0.0039 -> near (listed first = simpler) is selected; once near
    is uniformly 0.05 worse (SE 0) the best is selected; a candidate with CI lower bound <= 0 is not eligible."""
    va, lineage = setup["va"], setup["ctx"]["lineage"]
    cfg = {"look_id": LOOK, "labels": ["L_prob"], "reference": "H1", "G_star": "Human1", "B": 2000, "seed": 11,
           "candidates": ["near", "best", "weak"], "decision_rule": "toy"}
    grid = pd.DataFrame([{"label": "L_prob", "geneset": "all", "universe": "U", "stratum": "none", "method": m,
                          "delta_lin_vs_H1": d, "delta_lo_vs_H1": lo, "delta_hi_vs_H1": hi}
                         for m, d, lo, hi in (("near", 0.0975, 0.04, 0.15), ("best", 0.10, 0.05, 0.15),
                                              ("weak", 0.02, -0.01, 0.05))])
    shift = np.zeros(len(va))
    shift[[0, 1, 2]] = 0.05  # lineage VA: near worse on 3 lines
    shift[[5, 6]] = -0.05  # lineage VB: near better on 2 lines
    lts = {"H1": _lt(0.5, va), "best": _lt(0.6, va), "near": _lt(0.6 - shift, va), "weak": _lt(0.5, va)}
    dec = ev.decide_shortlist(grid, lts, lineage, cfg)
    c = {r["method"]: r for r in dec["candidates"]}
    assert [r["eligible"] for r in dec["candidates"]] == [True, True, False]
    assert c["near"]["diff_vs_best"] == pytest.approx(0.0025, abs=1e-12)
    assert 0.0025 < c["near"]["se_diff"] < 0.0055 and c["near"]["within_1se"]
    assert c["best"]["diff_vs_best"] == 0 and c["best"]["within_1se"] and "se_diff" not in c["weak"]
    assert dec["best_eligible"] == "best" and dec["selected"] == "near" and dec["gate_passed"]
    lts["near"] = _lt(0.55, va)
    dec = ev.decide_shortlist(grid, lts, lineage, cfg)
    assert dec["selected"] == "best" and not dec["candidates"][0]["within_1se"]


def test_shortlist_gate_none_eligible(setup):
    va = setup["va"]
    cfg = {"look_id": LOOK, "labels": ["L_prob"], "reference": "H1", "G_star": "Human1", "B": 2000, "seed": 11,
           "candidates": ["a", "b"], "decision_rule": "toy"}
    grid = pd.DataFrame([{"label": "L_prob", "geneset": "all", "universe": "U", "stratum": "none", "method": m,
                          "delta_lin_vs_H1": d, "delta_lo_vs_H1": lo, "delta_hi_vs_H1": 0.1}
                         for m, d, lo in (("a", 0.01, -0.02), ("b", 0.03, 0.0))])
    lts = {"H1": _lt(0.5, va), "a": _lt(0.5, va), "b": _lt(0.5, va)}
    dec = ev.decide_shortlist(grid, lts, setup["ctx"]["lineage"], cfg)
    assert not dec["gate_passed"] and dec["selected"] is None and dec["best_point_estimate"] == "b"
    assert dec["winners_curse"]["method"] == "b" and dec["winners_curse"]["m"] == 2


def test_shortlist_winners_curse_known_answer(setup):
    """§9 winner's curse (reporting only): best Δ = +0.1 on 3 of the 5 VA lines, 0 elsewhere -> Δ̄_lin = 0.015 and
    σ ≈ (0.1·√(0.6·0.4/5))/4 = 0.00548 (lineage-stratified bootstrap); term = σ·√(2 ln m), m = #candidates;
    m = 1 gives 0. The gate and the selection do not change."""
    from nextgem import metrics

    va, lineage = setup["va"], setup["ctx"]["lineage"]
    up = np.zeros(len(va))
    up[[0, 1, 2]] = 0.1
    lts = {"H1": _lt(0.5, va), "best": _lt(0.5 + up, va), "other": _lt(0.5, va), "third": _lt(0.5, va)}
    grid = pd.DataFrame([{"label": "L_prob", "geneset": "all", "universe": "U", "stratum": "none", "method": m,
                          "delta_lin_vs_H1": d, "delta_lo_vs_H1": lo, "delta_hi_vs_H1": 0.03}
                         for m, d, lo in (("other", 0.0, 0.0), ("best", 0.015, 0.004), ("third", 0.0, 0.0))])
    cfg = {"look_id": LOOK, "labels": ["L_prob"], "reference": "H1", "G_star": "Human1", "B": 2000, "seed": 11,
           "candidates": ["other", "best", "third"], "decision_rule": "toy"}
    wc = ev.decide_shortlist(grid, lts, lineage, cfg)["winners_curse"]
    rep = metrics.bootstrap_mean(pd.Series(up, index=va), lineage.reindex(va), B=2000, seed=11).replicates
    assert wc["method"] == "best" and wc["m"] == 3 and wc["best_delta_lin"] == 0.015
    assert wc["sigma"] == pytest.approx(float(np.std(rep, ddof=1)), rel=1e-12)
    assert wc["sigma"] == pytest.approx(0.1 * np.sqrt(0.6 * 0.4 / 5) / 4, rel=0.1)
    assert wc["term"] == pytest.approx(wc["sigma"] * np.sqrt(2 * np.log(3)), rel=1e-12)
    assert wc["best_delta_lin_minus_term"] == pytest.approx(0.015 - wc["term"], abs=1e-15)
    one = ev.decide_shortlist(grid, lts, lineage, dict(cfg, candidates=["best"]))
    assert one["winners_curse"]["term"] == 0 and one["winners_curse"]["best_delta_lin_minus_term"] == 0.015
    assert one["selected"] == "best" and one["gate_passed"]


def test_shortlist_dry_run_reads_no_validation_labels(setup, monkeypatch):
    from nextgem import labels, ledger

    def no_validation(*a, **k):
        raise AssertionError("validation labels read")

    monkeypatch.setattr(labels, "validation", no_validation)
    res = ev.dry_run(ev.shortlist_config(LOOK), setup["ctx"], setup["train"])
    assert res["scalar_guard"] == "pass" and res["n_rows"] == 16 and not res["validation_labels_read"]
    assert set(res["synthetic_decision"]) >= {"selected", "gate_passed", "candidates"}
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0
    assert not (setup["repo"] / "results/w9/val" / f"val_look_{LOOK}.tsv").exists()


def _patch_real(setup, monkeypatch):
    """main() on the toy: the real-data context and training labels replaced; validation labels forbidden."""
    from nextgem import labels

    def no_validation(*a, **k):
        raise AssertionError("validation labels read")

    monkeypatch.setattr(ev, "real_context", lambda cfg: setup["ctx"])
    monkeypatch.setattr(labels, "train", lambda *a, **k: setup["train"])
    monkeypatch.setattr(labels, "validation", no_validation)


def test_main_dry_run_prints_no_training_rates(setup, monkeypatch, capsys):
    from nextgem import ledger

    _patch_real(setup, monkeypatch)
    assert ev.main(["--look", LOOK, "--dry-run"]) == 0
    text = capsys.readouterr().out
    res = json.loads(text)
    assert res["scalar_guard"] == "pass" and res["ledger_rows_written"] == 0 and res["n_rows"] == 16
    rate = float(setup["train"].L_prob.mean())
    assert "bernoulli_rates" not in text and repr(rate) not in text and f"{rate:.6g}" not in text
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0


# ---- review fixes (notes/reviews/M4_harness_shortlist_review.md F2-F5)
def test_spec_workstream_must_match_directory(setup):
    """F2 (A11g): the W7 sub-cap counts looks by workstream, so a spec's workstream must be its directory's."""
    sha = ev._sha(setup["repo"] / SCORES_REL)
    for ws in ("W3", "w9", "W9a", "M9", "W09"):
        _write_spec(setup["repo"], sha256=sha, workstream=ws)
        with pytest.raises(ev.PreflightError, match="workstream"):
            ev.shortlist_config(LOOK)
    _write_spec(setup["repo"], sha256=sha, workstream="W9")
    assert ev.shortlist_config(LOOK)["workstream"] == "W9"


def test_spec_milestone_workstream_in_own_directory(setup):
    """A13: a milestone-level look (M<n>) lives in results/m<n>/val and is not a W7 look."""
    import shutil
    repo, look = setup["repo"], "TOY_M9_M6"
    scores = "results/m9/val/m9_val_scores.parquet"
    (repo / scores).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(repo / SCORES_REL, repo / scores)
    spec = _write_spec(repo, sha256=ev._sha(repo / scores), look_id=look, workstream="M9", milestone="M9")
    spec["scores"]["path"] = scores
    (repo / SPEC_REL).unlink()
    (repo / f"results/m9/val/look_spec_{look}.json").write_text(json.dumps(spec, indent=1))
    assert ev.shortlist_config(look)["workstream"] == "M9"
    for ws, ms in (("W9", "M9"), ("m9", "M9"), ("X9", "M9"), ("M9", "M6"), ("M09", "M09")):
        spec["workstream"], spec["milestone"] = ws, ms
        (repo / f"results/m9/val/look_spec_{look}.json").write_text(json.dumps(spec, indent=1))
        with pytest.raises(ev.PreflightError, match="workstream"):
            ev.shortlist_config(look)


def test_spec_primary_label_is_l_prob(setup):
    """F3: the gate cell is the first label, which must be L_prob (§8, A11g)."""
    sha = ev._sha(setup["repo"] / SCORES_REL)
    _write_spec(setup["repo"], sha256=sha, labels=["L_ge05", "L_prob"])
    with pytest.raises(ev.PreflightError, match="L_prob"):
        ev.shortlist_config(LOOK)
    _write_spec(setup["repo"], sha256=sha, labels=["L_prob", "L_ge05"])
    assert ev.shortlist_config(LOOK)["labels"][0] == "L_prob"


@pytest.mark.parametrize("over", [pytest.param({"B": 2000.0}, id="B_float"), pytest.param({"seed": 11.0}, id="seed_float"),
                                  pytest.param({"B": True}, id="B_bool")])
def test_non_integer_b_or_seed_refused(setup, over):
    """F4: 2000.0 == 2000 and 11.0 == 11, but they break the bootstrap after ``started``."""
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL), **over)
    with pytest.raises(ev.PreflightError, match="integer"):
        ev.shortlist_config(LOOK)
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL))
    cfg = dict(ev.shortlist_config(LOOK), **over)  # prepare_shortlist checks again (before started)
    with pytest.raises(ev.PreflightError, match="integers"):
        ev.prepare_shortlist(cfg, setup["ctx"], setup["train"])


@pytest.mark.parametrize("name,why", [pytest.param("GENEA", "gene", id="gene_symbol"),
                                      pytest.param("VA", "lineage name", id="lineage")])
def test_unit_key_method_name_refused_before_started(setup, monkeypatch, name, why):
    """F4: a method named like a gene symbol or lineage passes the static spec check but would fail the scalar
    guard after ``started``; prepare_shortlist (and so --dry-run) refuses it before any ledger row."""
    from nextgem import ledger

    _write_spec(setup["repo"], sha256=_write_scores(setup["repo"], setup["scores"].assign(**{name: 0.2})),
                methods=["H1", "H2", SAME, PERFECT, name])
    cfg = ev.shortlist_config(LOOK)
    with pytest.raises(ev.PreflightError, match=why):
        ev.prepare_shortlist(cfg, setup["ctx"], setup["train"])
    _patch_real(setup, monkeypatch)
    assert ev.main(["--look", LOOK, "--dry-run"]) == 2
    assert ledger.read("validation") == [] and setup["calls"]["val"] == 0


def test_missing_or_invalid_gstar_refused(setup, monkeypatch):
    """F5: no G-star file (or a bad one) is a pre-flight refusal, before any label access."""
    hit = _forbid_labels(monkeypatch)
    g = setup["repo"] / ev.W1_GSTAR
    for content in (None, "{not json", json.dumps({"other": 1}), json.dumps({"G_star": "Human3"})):
        if content is None:
            g.unlink()
        else:
            g.write_text(content)
        with pytest.raises(ev.PreflightError, match="G-star file"):
            ev.shortlist_config(LOOK)
        assert ev.main(["--look", LOOK, "--print-config"]) == 2
    assert hit == []


def test_duplicate_look_ids_refused(setup):
    """F5: the same look_id in two results/<ws>/val directories is refused, not shadowed."""
    other = setup["repo"] / f"results/w8/val/look_spec_{LOOK}.json"
    other.parent.mkdir(parents=True)
    other.write_text((setup["repo"] / SPEC_REL).read_text())
    with pytest.raises(ev.PreflightError, match="more than one look spec"):
        ev.look_specs()
    assert ev.main(["--look", LOOK, "--print-config"]) == 2


def test_pinned_problems_and_required_files(setup):
    """Provenance pins: a changed file is a pre-flight problem (a missing one is reported by preflight());
    required_files = code + W1 context inputs + spec + scores + G-star file + provenance."""
    a = setup["repo"] / "results/w9/train/a.json"
    a.parent.mkdir(parents=True)
    a.write_text("{}")
    prov = {"results/w9/train/a.json": ev._sha(a), "scripts/w9_missing.py": "0" * 64}
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL), provenance=prov)
    cfg = ev.shortlist_config(LOOK)
    assert cfg["provenance"] == prov and ev.pinned_problems(cfg) == []
    context = {v for k, v in ev.input_files(False).items() if k != "sesoi_c1b"}
    assert set(ev.required_files(cfg)) == set(ev.CODE_FILES) | context | {SPEC_REL, SCORES_REL, ev.W1_GSTAR} | set(prov)
    assert any("scripts/w9_missing.py: missing" in p for p in ev.preflight(ev.required_files(cfg), root=setup["repo"]))
    a.write_text('{"changed": 1}')
    assert ev.pinned_problems(cfg) == ["results/w9/train/a.json: SHA-256 differs from the look spec's provenance pin"]


def test_frozen_candidates_reports_missing_provenance(setup, capsys):
    """F5: --print-frozen-candidates lists a missing provenance file instead of crashing."""
    a = setup["repo"] / "results/w9/train/a.json"
    a.parent.mkdir(parents=True)
    a.write_text("{}")
    prov = {"results/w9/train/a.json": ev._sha(a), "scripts/w9_missing.py": "0" * 64}
    _write_spec(setup["repo"], sha256=ev._sha(setup["repo"] / SCORES_REL), provenance=prov)
    assert ev.main(["--look", LOOK, "--print-frozen-candidates"]) == 0
    rows = capsys.readouterr().out.splitlines()
    assert "# MISSING\tscripts/w9_missing.py\t(no row possible until the file exists)" in rows
    row = next(r.split("\t") for r in rows if r.startswith("results/w9/train/a.json\t"))
    assert row[1:4] == ["file", "results/w9/train/a.json", ev._sha(a)] and row[6] == "W9 validation look: input"
    assert not any(r.startswith("sesoi.C1b") for r in rows)
