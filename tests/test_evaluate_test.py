"""scripts/evaluate_test.py on synthetic toy data only. No real label file is touched: training, validation and test
labels are synthetic frames handed in by monkeypatching ``labels.train/validation/test`` (exactly as the
evaluate_val tests do), NEXTGEM_STATE_ROOT points to a toy git repository, NEXTGEM_DATA_ROOT to a toy directory,
and sealed writes go to a recorder that replaces ``labels._write_protected`` (nothing is written under data/).

Covered:
- known answers: the sign-flip p (all-zero, all-positive, exact enumeration and the replicated random stream), the
  lineage-stratified bootstrap CI, C3's permutation-null percentile, fixed-sequence gatekeeping that stops at the
  first failure (T1 and T2 states), claim-language tiers with each contrast's SESOI, the lineage t-interval with
  df = 6 and the exact sign flip over 2^7, k-of-7;
- spec checks before any label access; the look order (T1 needs the rehearsal; T2 needs the W7 stop row naming the
  P2 spec; unseal needs a completed T2);
- the rehearsal's validation stage is a validation look with its ledger rows; T1 writes sealed outputs only under
  data/locked/sealed/T1/ and only their SHA-256 into the ledger; T1 stdout is scalar-only;
- an end-to-end chain rehearsal → T1 → W7 stop → T2 → unseal with SHA-verified publication;
- no dry mode reads a validation or test label or writes a ledger row or a file.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import json
import math
import subprocess

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import toy

from nextgem import metrics, paths

SPEC = importlib.util.spec_from_file_location("evaluate_test", paths.SCRIPTS_DIR / "evaluate_test.py")
et = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(et)

GENES = list(range(1001, 1041))  # 40 toy genes
TRAIN = {"TA": 6, "TB": 6, "TC": 6}
VAL = {"VA": 4, "VB": 4, "VC": 4, "VD": 4}
TEST = {f"K{c}": 4 for c in "ABCDEFG"}  # 7 test lineages
METHODS = ["H1", "H2", "P1[H1]", "P1[H2]", "P1g[H1]"]
ARMS = {"P1[H1]": "H1", "P1[H2]": "H2"}
SPEC_ID = "TOYP1"
REH = et.REHEARSAL_PREFIX + SPEC_ID
FROZEN_HDR = "item\tkind\tpath\tsha256_or_value\tfrozen_on\tamendment\tnote\n"


# ----------------------------------------------------------------------------- toy world
def _ids(spec, start):
    out, lin = [], {}
    for name, n in spec.items():
        for _ in range(n):
            mid = f"ACH-{start + len(out):06d}"
            out.append(mid)
            lin[mid] = name
    return out, lin


def _long(lines, rng, essential, lineage=None):
    rows = []
    for l in lines:
        ge = rng.normal(0, 0.3, len(GENES)) - 1.2 * essential + rng.normal(0, 0.2, len(GENES))
        for g, v in zip(GENES, ge):
            rows.append((l, g, v, 1 / (1 + np.exp(4 * (v + 0.5)))))
    df = pd.DataFrame(rows, columns=["ModelID", "entrez", "gene_effect", "dep_prob"])
    df["L_prob"] = (df.dep_prob > 0.5).astype(float)
    df["L_ge05"] = (df.gene_effect < -0.5).astype(float)
    df["L_ge06"] = (df.gene_effect < -0.6).astype(float)
    df["L_ge10"] = (df.gene_effect < -1.0).astype(float)
    if lineage is not None:
        df["lineage"] = df.ModelID.map(lineage)
    return df


def _parquet(repo, rel, df):
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(p, index=False)
    return et._sha(p)


def _scores(lines, gen, essential, rng, layout):
    base = {"H1": gen.e_H1.to_numpy(), "H2": gen.e_H2.to_numpy()}
    p1 = np.round(gen.e_H1.to_numpy() + 0.3 * essential, 6)
    p2 = np.round(gen.e_H2.to_numpy() + 0.3 * essential, 6)
    if layout == "gene":
        return pd.DataFrame({"entrez": GENES, **base, "P1[H1]": p1, "P1[H2]": p2, "P1g[H1]": np.round(p1 * 0.9, 6)})
    rows = []
    for l in lines:
        jit = np.round(rng.random(len(GENES)) * 0.01, 6)
        rows.append(pd.DataFrame({"ModelID": l, "entrez": GENES, **base, "P1[H1]": p1 + jit, "P1[H2]": p2 + jit,
                                  "P1g[H1]": np.round(p1 * 0.9, 6) + jit}))
    return pd.concat(rows, ignore_index=True)


def _calls(lines, gen, essential, layout):
    parts = []
    for arm, s in ARMS.items():
        call = (gen[f"e_{s}"].to_numpy() + 0.3 * essential) > 0.6
        d = pd.DataFrame({"entrez": GENES, "method": arm, "p57": call, "pgr": call, "fba": call})
        parts.append(d if layout == "gene" else pd.concat([d.assign(ModelID=l) for l in lines], ignore_index=True))
    out = pd.concat(parts, ignore_index=True)
    return out[(["ModelID"] if layout == "line_gene" else []) + ["entrez", "method", "p57", "pgr", "fba"]]


def p1_spec(sha, **over):
    spec = {"spec_id": SPEC_ID, "kind": "p1", "milestone": "M6", "primary": "P1[H1]", "methods": list(METHODS),
            "arms": dict(ARMS),
            "scores": {"validation": {"path": "results/m6/val_scores.parquet", "sha256": sha["sv"], "layout": "gene"},
                       "test": {"path": "results/m6/tmp/test_scores.parquet", "sha256": sha["st"], "layout": "line_gene"}},
            "calls": {"validation": {"path": "results/m6/val_calls.parquet", "sha256": sha["cv"], "layout": "gene"},
                      "test": {"path": "results/m6/tmp/test_calls.parquet", "sha256": sha["ct"], "layout": "line_gene"}},
            "failed_builds": {"validation": None, "test": {"path": "results/m6/tmp/test_failed.parquet", "sha256": sha["ft"]}},
            "attribution": {"i": ["P1[H1]", "G*"], "ii": et.LEVER_NOT_SELECTED, "iii": et.LEVER_NOT_SELECTED,
                            "iv": ["P1g[H1]", "G*"], "v": ["P1[H1]", "P1g[H1]"]},
            "task_calls": {"path": "results/m6/task_calls.parquet", "sha256": sha["tc"]},
            "gene_annotations": {"path": "results/m6/gene_annotations.parquet", "sha256": sha["ga"]},
            "bootstrap_method": "percentile", "fitted_on_labels": False, "notes": "toy"}
    spec.update(over)
    return spec


def p2_spec(sha, **over):
    spec = {"spec_id": "TOYP2", "kind": "p2", "milestone": "M7", "primary": "P2", "methods": ["H1", "H2", "P2"],
            "scores": {"test": {"path": "results/m6/p2_test_scores.parquet", "sha256": sha, "layout": "gene"}},
            "failed_builds": {"test": None}, "labels": ["L_prob"], "genesets": ["all"], "universes": ["U"],
            "strata": ["none"], "b6": {"features": "w1", "hgb": {"max_iter": 20, "learning_rate": 0.1,
                                                                 "max_leaf_nodes": 7, "early_stopping": False}},
            "bootstrap_method": "percentile", "fitted_on_labels": False}
    spec.update(over)
    return spec


def _b6_stub(fit_long, lines, genes, b6, seed):
    """B6 stand-in. Like w1_b6.stacked_priors it needs a lineage on every fitting row (training ∪ validation)."""
    assert fit_long["lineage"].notna().all(), "fitting rows without lineage (validation frames carry none)"
    return pd.DataFrame(np.round(np.random.default_rng(seed).random((len(lines), len(genes))), 8),
                        index=list(lines), columns=list(genes))


def write_json(repo, rel, doc):
    p = repo / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(doc, indent=1))
    return et._sha(p)


@pytest.fixture()
def world(tmp_path, monkeypatch):
    repo = tmp_path / "state"
    toy.init_repo(repo)
    data = tmp_path / "data"
    data.mkdir()
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(repo))
    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(data))
    # small resampling counts on toy data (the real constants are checked by check_constants in main)
    for k, v in (("B_CONFIRMATORY", 300), ("N_SIGN_FLIPS", 400), ("N_C3_PERMUTATIONS", 50), ("B_OTHER", 200)):
        monkeypatch.setattr(metrics, k, v)
    rng = np.random.default_rng(17)
    essential = (rng.random(len(GENES)) < 0.3).astype(float)
    tr, lin_t = _ids(TRAIN, 1)
    va, lin_v = _ids(VAL, 300)
    te, lin_e = _ids(TEST, 600)
    noexpr, lin_n = _ids({"KA": 2}, 900)
    lineage = pd.Series({**lin_t, **lin_v, **lin_e, **lin_n})
    train = _long(tr, rng, essential, lineage)
    val = _long(va, rng, essential)
    test = _long(te + noexpr, rng, essential)
    ratio = np.where((essential > 0) & (rng.random(len(GENES)) < 0.6), 0.0, rng.choice([0.4, 0.8, 1.0], len(GENES)))
    gen = pd.DataFrame(index=pd.Index(GENES, name="entrez"))
    for s in ("H1", "H2"):
        r = ratio if s == "H1" else np.where(rng.random(len(GENES)) < 0.15, 1 - ratio, ratio)
        gen[f"in_{s}"] = [g != (1040 if s == "H1" else 1039) for g in GENES]
        gen[f"ratio_{s}"] = r
        gen[f"e_{s}"] = 1 - r
        gen[f"fba_call_{s}"] = r < 0.01
        gen[f"p57_{s}"] = r < 0.01
        gen[f"pgr_{s}"] = r < 0.01
        gen[f"et_{s}"] = 1 - r
    expr_lines = tr + va + te
    expr = pd.DataFrame(rng.gamma(2, 1.5, (len(expr_lines), len(GENES))), index=expr_lines, columns=GENES).astype("float32")
    cn = pd.DataFrame(rng.normal(1, 0.8, (len(expr_lines), len(GENES))), index=expr_lines, columns=GENES).astype("float32")
    all_lines = expr_lines + noexpr
    lf = pd.DataFrame({"base_medium_class": "sheet:RPMI", "library_class": "Avana-containing", "model_type": "Cell Line",
                       "has_expression_matrix": [l in set(expr_lines) for l in all_lines], "serum": 1.0},
                      index=pd.Index(all_lines, name="ModelID"))
    lf.loc[te[:3], "model_type"] = "Organoid"
    gf = pd.DataFrame({"addressable_human1": [g % 2 == 0 for g in GENES], "addressable_human2": [g % 3 == 0 for g in GENES]},
                      index=pd.Index(GENES, name="entrez"))
    base_ctx = {"genes": GENES, "ucap": GENES[:35], "gen": gen, "omics": {"expr": expr, "cn": cn}, "lf": lf,
                "lineage": lineage, "line_class": lf.base_medium_class, "all_ids": all_lines, "depmap_ce": set(GENES[:2]),
                "gene_symbols": {"GENEA", "GENEB"}, "gene_features": gf,
                "strata": {"none": ("none", None), "expressed_pairs": ("pairs", "expressed"),
                           "organoid": ("rows", ("model_type", "Organoid")), "lib_KY": ("rows", ("library_class", "KY-only")),
                           "multicopy_in": ("genes_in", set(GENES[:5])), "cn_amp4_in": ("pairs", "cn_amp4_in")},
                "b6_fit": _b6_stub}

    def ctx_for(tier):
        tl = te + noexpr if tier == "test" else va
        return {**base_ctx, "tier": tier, "tier_lines": sorted(tl), "eval_lines": sorted(l for l in tl if l in set(expr_lines))}

    # committed state: FROZEN values, the W1 G-star file, stand-ins for the W1 context inputs, spec and its files
    (repo / "FROZEN.tsv").write_text(FROZEN_HDR + "label.primary\tvalue\t\tL_prob\t2026-10-26\t\t\n"
                                     "sesoi.C1b\tvalue\t\t0.0057228089\t2026-10-27\t\t\n")
    for rel in [*et.ev.input_files(False).values(), "results/w1/tmp/mut_dam.parquet", "results/w1/tmp/mut_hot.parquet",
                "results/w1/tmp/gene_neighbours.json"]:
        (repo / rel).parent.mkdir(parents=True, exist_ok=True)
        (repo / rel).write_text(f"toy {rel}\n")
    write_json(repo, et.ev.W1_GSTAR, {"G_star": "Human1"})
    te_sorted = sorted(te)
    sha = {"sv": _parquet(repo, "results/m6/val_scores.parquet", _scores(va, gen, essential, rng, "gene")),
           "st": _parquet(repo, "results/m6/tmp/test_scores.parquet", _scores(te_sorted, gen, essential, rng, "line_gene")),
           "cv": _parquet(repo, "results/m6/val_calls.parquet", _calls(va, gen, essential, "gene")),
           "ct": _parquet(repo, "results/m6/tmp/test_calls.parquet", _calls(te_sorted, gen, essential, "line_gene")),
           "ft": _parquet(repo, "results/m6/tmp/test_failed.parquet", pd.DataFrame({"ModelID": [te_sorted[5]], "method": ["P1[H2]"]})),
           "tc": _parquet(repo, "results/m6/task_calls.parquet", pd.DataFrame(
               {"entrez": GENES * 2, "task": ["ESS001"] * len(GENES) + ["ESS002"] * len(GENES),
                "fails": [bool(e) for e in essential] + [g % 5 == 0 for g in GENES]})),
           "ga": _parquet(repo, "results/m6/gene_annotations.parquet", pd.DataFrame(
               {"entrez": GENES, "paralog_buffered": [g % 4 == 0 for g in GENES]}))}
    write_json(repo, et.P1_SPEC, p1_spec(sha))
    p2_sha = _parquet(repo, "results/m6/p2_test_scores.parquet",
                      pd.DataFrame({"entrez": GENES, "H1": gen.e_H1.to_numpy(), "H2": gen.e_H2.to_numpy(),
                                    "P2": np.round(gen.e_H1.to_numpy() + 0.4 * essential, 6)}))
    write_json(repo, et.P2_SPEC, p2_spec(p2_sha))
    (repo / ".gitignore").write_text("results/**/tmp/\n")
    toy.commit_all(repo, "toy world")
    calls = {"train": 0, "validation": 0, "test": 0}
    sealed = []

    def recorder(path, data, category):
        sealed.append((str(path), category, hashlib.sha256(data).hexdigest()))

    from nextgem import labels

    monkeypatch.setattr(labels, "_write_protected", recorder)
    return {"repo": repo, "data": data, "ctx_for": ctx_for, "train": train, "val": val, "test": test, "calls": calls,
            "sealed": sealed, "lineage": lineage, "te": te_sorted, "p2_sha": p2_sha, "sha": sha}


def _committed_rows(repo, name):
    head = subprocess.run(["git", "-C", str(repo), "show", f"HEAD:results/{name}_ledger.tsv"], capture_output=True, text=True).stdout
    return [ln.split("\t") for ln in head.splitlines()[1:]]


def _patch_labels(w, monkeypatch, *, allow=("train", "validation", "test"), check_started=None):
    """labels.* return the toy frames; a tier not in ``allow`` raises. ``check_started`` = (ledger, look_id) that
    must already be committed as 'started' when that ledger's evaluation tier is first read."""
    from nextgem import labels

    def make(tier, frame):
        def f(*a, **k):
            if tier not in allow:
                raise AssertionError(f"labels.{tier}() called in a dry mode")
            if check_started and tier == {"validation": "validation", "test": "test"}.get(check_started[0]):
                rows = _committed_rows(w["repo"], check_started[0])
                assert any(r[3] == "started" and r[4] == check_started[1] for r in rows), "labels read before 'started'"
            w["calls"][tier] += 1
            if tier == "train" and k.get("columns"):
                return frame[["ModelID", "entrez", "lineage", *k["columns"]]]
            return frame.copy()
        return f

    monkeypatch.setattr(labels, "train", make("train", w["train"]))
    monkeypatch.setattr(labels, "validation", make("validation", w["val"]))
    monkeypatch.setattr(labels, "test", make("test", w["test"]))


def _patch_main(w, monkeypatch):
    from nextgem import run

    monkeypatch.setattr(et, "real_context", lambda cfg, tier: w["ctx_for"](tier))
    monkeypatch.setattr(et, "preflight_problems", lambda cfg: [])
    monkeypatch.setattr(run, "git_dirty", lambda *a, **k: False)


def _rehearse(w, monkeypatch):
    """The rehearsal's validation stage through run_p1_look (validation labels read after 'started')."""
    from nextgem import labels

    _patch_labels(w, monkeypatch, check_started=("validation", REH))
    cfg = et.look_config("REHEARSAL")
    ctx = w["ctx_for"]("validation")
    prep = et.prepare(cfg, ctx, labels.train(columns=et.ev.LABEL_COLS), "validation")
    out_dir = w["repo"] / et.REHEARSAL_DIR
    et.write_config(cfg, out_dir, "val_look")
    res = et.run_p1_look(cfg, ctx=ctx, prep=prep, load_labels=labels.validation, out_dir=out_dir,
                         ledger_name="validation", author="test", workstream=et.REHEARSAL_WS)
    toy.commit_all(w["repo"], "rehearsal outputs")
    return cfg, res


# ============================================================================= known answers: statistics
def _cfg(**over):
    from nextgem.config import SEEDS

    c = {"seed": SEEDS[0], "B_conf": 10_000, "n_flip": 10_000, "n_perm": 2_000, "bootstrap_method": "percentile",
         "alpha": 0.025}
    c.update(over)
    return c


def _lineage_series(spec):
    ids, lin = _ids(spec, 1)
    return ids, pd.Series(lin)


def test_sign_flip_p_known_answers():
    """p = (1 + #{T* >= T_obs}) / (1 + 10,000): all-zero Δ -> every flip ties -> p = 1; 24 equal positive Δ ->
    no flip reaches T_obs except all '+' (prob 2^-24) -> p = 1/10,001."""
    ids, lin = _lineage_series({"A": 8, "B": 8, "C": 8})
    r0 = et.confirmatory(pd.Series(0.0, index=ids), lin, "C1", 0.01, _cfg())
    assert r0["p_sign_flip"] == 1.0 and not r0["passed"]
    r1 = et.confirmatory(pd.Series(0.02, index=ids), lin, "C1", 0.01, _cfg())
    assert r1["p_sign_flip"] == pytest.approx(1 / 10_001, abs=1e-15) and r1["passed"]
    assert r1["estimate"] == pytest.approx(0.02) and r1["ci_low"] == pytest.approx(0.02) and r1["ci_high"] == pytest.approx(0.02)


def test_sign_flip_p_matches_exact_enumeration_and_stream():
    """4 lines in 2 lineages: the exact one-sided p over the 16 sign patterns of the lineage-equal mean; the
    10,000-flip Monte-Carlo p is within 4 SE of it and equals the count replicated from the spawned SEEDS[0] stream."""
    ids, lin = _lineage_series({"A": 2, "B": 2})
    d = pd.Series([0.05, 0.01, 0.03, -0.02], index=ids)
    w = metrics.lineage_weights(lin.reindex(ids)).to_numpy()
    obs = float(w @ d.to_numpy())
    exact = np.mean([float(w @ (np.array(s) * d.to_numpy())) >= obs - 1e-12 for s in itertools.product((1, -1), repeat=4)])
    assert exact == pytest.approx(3 / 16)
    r = et.confirmatory(d, lin, "C1", 0.01, _cfg())
    assert abs(r["p_sign_flip"] - exact) < 4 * math.sqrt(exact * (1 - exact) / 10_000)
    _, ss_flip = np.random.SeedSequence(11).spawn(2)
    rng = np.random.default_rng(ss_flip)
    count, done = 0, 0
    while done < 10_000:
        b = min(2_000, 10_000 - done)
        signs = rng.integers(0, 2, size=(b, 4), dtype=np.int8) * 2 - 1
        count += int(np.sum(signs @ (w * d.to_numpy()) >= obs - 1e-12 * max(1.0, abs(obs))))
        done += b
    assert r["p_sign_flip"] == (1 + count) / 10_001


def test_bootstrap_ci_known_answer():
    """Lineage A = {0, 1}, B = {0, 0}: Δ̄_lin = 0.25; a replicate is mean(A*)/2 ∈ {0, 0.25, 0.5} with
    probabilities 1/4, 1/2, 1/4, so the 95% percentile CI over B = 10,000 is [0, 0.5]."""
    ids, lin = _lineage_series({"A": 2, "B": 2})
    a = pd.DataFrame({"tap": [1.0, 0.0, 0.0, 0.0], "n": 10, "n_pos": 5, "estimable": True}, index=ids)
    b = pd.DataFrame({"tap": 0.0, "n": 10, "n_pos": 5, "estimable": True}, index=ids)
    r = et.contrast_est(a, b, lin, B=10_000, seed=11)
    assert r["estimate"] == pytest.approx(0.25) and r["ci_lo"] == 0.0 and r["ci_hi"] == 0.5
    assert r["line_mean"] == pytest.approx(0.25) and r["n_lines"] == 4 and r["n_lineages"] == 2
    c = et.confirmatory(pd.Series([1.0, 0.0, 0.0, 0.0], index=ids), lin, "C1", 0.01, _cfg())
    assert (c["estimate"], c["ci_low"], c["ci_high"]) == (pytest.approx(0.25), 0.0, 0.5)


def _c3_toy(perfect: bool):
    ids, lin = _lineage_series({"A": 12, "B": 12})
    rng = np.random.default_rng(3)
    genes = list(range(1, 31))
    Y = pd.DataFrame((rng.random((len(ids), len(genes))) < 0.4).astype(float), index=ids, columns=genes)
    Y.iloc[:6, :] = 1.0  # >= 5 positive and >= 5 negative lines per gene
    Y.iloc[6:12, :] = 0.0
    S = Y + 0.0 if perfect else pd.DataFrame(np.tile(rng.random(len(genes)), (len(ids), 1)), index=ids, columns=genes)
    B3 = pd.DataFrame(rng.random(Y.shape), index=ids, columns=genes)
    return S, B3, Y, lin


def test_c3_permutation_null_percentile():
    """C3's null threshold is the 97.5th percentile of the within-lineage permutation null; a line-constant P1
    has mΔAP_gene = 0 = every null draw (fails); P1 = the labels has tAP_gene = 1 per gene (passes)."""
    S, B3, Y, lin = _c3_toy(perfect=False)
    r = et.c3(S, B3, Y, lin, _cfg(B_conf=500, n_perm=200))
    assert r["m_delta_ap_gene_p1"] == pytest.approx(0.0, abs=1e-12) and r["null_threshold"] == pytest.approx(0.0, abs=1e-12)
    assert not r["permutation_passed"] and not r["passed"] and r["permutation_p"] == 1.0
    S, B3, Y, lin = _c3_toy(perfect=True)
    r = et.c3(S, B3, Y, lin, _cfg(B_conf=500, n_perm=200))
    _, ss_perm = np.random.SeedSequence(11).spawn(2)
    null = metrics.c3_permutation_null(S, Y, lin, n_perm=200, seed=ss_perm)
    assert r["null_threshold"] == pytest.approx(float(np.percentile(null["null"], 97.5)), abs=1e-15)
    g = metrics.per_gene(S, Y)
    assert r["m_delta_ap_gene_p1"] == pytest.approx(float((1.0 - g.e_random).mean()), abs=1e-12)
    assert r["permutation_passed"] and r["m_delta_ap_gene_p1"] > r["null_threshold"]
    assert r["permutation_p"] == pytest.approx(1 / 201) and r["bootstrap_passed"] and r["passed"]


def _res(passed, est=0.02, lo=0.012, sesoi=0.01, name="C1"):
    return {"contrast": name, "passed": passed, "estimate": est, "ci_low": lo, "sesoi": sesoi}


@pytest.mark.parametrize("pattern,statuses", [
    pytest.param((True, False, True), ["confirmed", "not confirmed", "estimate only", "estimate only", "estimate only"],
                 id="C1b_fails_stops_chain"),
    pytest.param((False, True, True), ["not confirmed", "estimate only", "estimate only", "estimate only", "estimate only"],
                 id="C1_fails_stops_chain"),
    pytest.param((True, True, True), ["confirmed", "confirmed", "confirmed", "pending", "pending"], id="T1_all_pass_C4_C5_pending"),
    pytest.param((True, True, False), ["confirmed", "confirmed", "not confirmed", "estimate only", "estimate only"],
                 id="C3_fails_C4_C5_estimates"),
])
def test_gatekeeping_stops_at_first_failure_t1(pattern, statuses):
    res = {"C1": _res(pattern[0]), "C1b": _res(pattern[1], sesoi=0.0057, name="C1b"),
           "C3": {"contrast": "C3", "passed": pattern[2], "estimate": 0.01, "ci_low": 0.002, "sesoi": None}}
    rows = et.gate_table(res)
    assert [r["contrast"] for r in rows] == ["C1", "C1b", "C3", "C4", "C5"]
    assert [r["status"] for r in rows] == statuses
    for r in rows:
        if r["status"] == "estimate only":
            assert r["claim"] == "estimate only (earlier gate failed)"
    if statuses[3] == "pending":
        assert rows[3]["claim"] == rows[4]["claim"] == "pending (T2)"


def test_gatekeeping_t2_chain_after_c3_failure():
    """T2 adds C4 and C5; a C3 failure at T1 makes both estimates even if their own tests pass."""
    res = {"C1": _res(True), "C1b": _res(True, sesoi=0.0057, name="C1b"),
           "C3": {"contrast": "C3", "passed": False, "estimate": -0.01, "ci_low": -0.02, "sesoi": None},
           "C4": _res(True, name="C4"), "C5": _res(True, name="C5")}
    rows = et.gate_table(res)
    assert [r["status"] for r in rows] == ["confirmed", "confirmed", "not confirmed", "estimate only", "estimate only"]
    assert rows[3]["claim"] == rows[4]["claim"] == "estimate only (earlier gate failed)"
    with pytest.raises(ValueError):
        et.gate_table({"C1": _res(True), "C3": res["C3"]})  # not a prefix of the frozen order


@pytest.mark.parametrize("name,sesoi,est,lo,passed,claim", [
    pytest.param("C1", 0.01, 0.015, 0.012, True, "meaningfully improves", id="C1_meaningful"),
    pytest.param("C1", 0.01, 0.015, 0.004, True, "improves", id="C1_improves"),
    pytest.param("C1", 0.01, 0.008, 0.002, True, "detectable improvement", id="C1_detectable"),
    pytest.param("C1", 0.01, 0.015, 0.012, False, "no detectable improvement", id="C1_not_passed"),
    pytest.param("C1b", 0.0057228089, 0.006, 0.003, True, "improves", id="C1b_lift_scaled_improves"),
    pytest.param("C1b", 0.0057228089, 0.0058, 0.0058, True, "meaningfully improves", id="C1b_lift_scaled_meaningful"),
    pytest.param("C3", None, 0.5, 0.4, True, "detectable improvement", id="C3_detectable_only"),
])
def test_claim_language_tiers(name, sesoi, est, lo, passed, claim):
    order = ["C1", "C1b", "C3"]
    res = {}
    for c in order[:order.index(name) + 1]:
        res[c] = _res(True, name=c, sesoi=0.0057228089 if c == "C1b" else (None if c == "C3" else 0.01))
    res[name] = {"contrast": name, "passed": passed, "estimate": est, "ci_low": lo, "sesoi": sesoi}
    rows = {r["contrast"]: r for r in et.gate_table(res)}
    assert rows[name]["claim"] == claim


def test_transfer_t_interval_df6_and_exact_sign_flip():
    """7 lineage means 0.01..0.07: mean 0.04, t-interval with t(0.975, 6) = 2.446912, exact one-sided sign-flip p
    over 2^7 = 1/128 (only all '+' reaches the observed mean); 7 of 7 positive >= ceil(0.7·7) = 5 -> generalises.
    Lines within a lineage differ, so the line deltas average to the lineage means."""
    lin_names = list("ABCDEFG")
    means = np.arange(1, 8) / 100
    ids, lin, d = [], {}, []
    for k, m in zip(lin_names, means):
        for off in (-0.004, 0.004):
            mid = f"ACH-{len(ids) + 1:06d}"
            ids.append(mid)
            lin[mid] = k
            d.append(m + off)
    t = et.transfer(pd.Series(d, index=ids), pd.Series(lin))
    half = stats.t.ppf(0.975, 6) * np.std(means, ddof=1) / math.sqrt(7)
    assert stats.t.ppf(0.975, 6) == pytest.approx(2.446912, abs=1e-6)
    assert t["t_mean"] == pytest.approx(0.04) and t["t_low"] == pytest.approx(0.04 - half) and t["t_high"] == pytest.approx(0.04 + half)
    assert t["exact_sign_flip_p"] == pytest.approx(1 / 128) and t["n_lineages"] == 7
    assert t["k_positive"] == 7 and t["k_required"] == 5 and t["generalises"]
    mixed = np.array([0.05, 0.04, 0.03, 0.02, -0.01, -0.02, -0.03])
    t2 = et.transfer(pd.Series(np.repeat(mixed, 2), index=ids), pd.Series(lin))
    signs = np.array(list(itertools.product((1.0, -1.0), repeat=7)))
    assert t2["exact_sign_flip_p"] == pytest.approx(float(np.mean(signs @ mixed / 7 >= mixed.mean() - 1e-12)))
    assert t2["k_positive"] == 4 and not t2["generalises"]


# ============================================================================= specs and look order
@pytest.mark.parametrize("over,match", [
    pytest.param({"attribution": {"i": ["P1[H2]", "G*"], "ii": "lever not selected", "iii": "lever not selected",
                                  "iv": "lever not selected", "v": "lever not selected"}}, "attribution.i", id="attribution_i_not_primary"),
    pytest.param({"attribution": {"i": ["P1[H1]", "G*"], "ii": "not chosen", "iii": "lever not selected",
                                  "iv": "lever not selected", "v": "lever not selected"}}, "attribution.ii", id="attribution_bad_token"),
    pytest.param({"fitted_on_labels": True}, "fitted_on_labels", id="label_fitted_p1_refused"),
    pytest.param({"bootstrap_method": "basic"}, "bootstrap_method", id="bootstrap_method"),
    pytest.param({"primary": "H1"}, "primary", id="primary_generic"),
    pytest.param({"arms": {"P1[H1]": "H1", "P1g[H1]": "H1"}}, "arms", id="two_arms_one_parent"),
    pytest.param({"methods": METHODS + ["L_prob"]}, "refused", id="label_like_method"),
    pytest.param({"methods": METHODS + ["ACH-000001"]}, "refused", id="line_id_method"),
    pytest.param({"spec_id": "bad id"}, "spec_id", id="spec_id_chars"),
    pytest.param({"extra": 1}, "keys", id="unknown_key"),
])
def test_p1_spec_fields_refused(world, monkeypatch, over, match):
    _patch_labels(world, monkeypatch, allow=())
    write_json(world["repo"], et.P1_SPEC, p1_spec(world["sha"], **over))
    with pytest.raises(et.PreflightError, match=match):
        et.look_config("T1")
    assert world["calls"] == {"train": 0, "validation": 0, "test": 0}


@pytest.mark.parametrize("tier,key,entry,match", [
    pytest.param("test", "scores", {"path": "results/m6/tmp/test_scores.parquet", "sha256": "0" * 64, "layout": "line_gene"},
                 "SHA-256 differs", id="test_scores_wrong_sha"),
    pytest.param("validation", "scores", {"path": "results/w3/val/x.parquet", "sha256": "0" * 64, "layout": "gene"},
                 "below results/m6", id="scores_outside_m6"),
    pytest.param("test", "calls", {"path": "results/m6/tmp/test_calls.parquet", "sha256": "0" * 64, "layout": "matrix"},
                 "layout", id="calls_bad_layout"),
])
def test_p1_spec_files_refused_before_label_access(world, monkeypatch, tier, key, entry, match):
    _patch_labels(world, monkeypatch, allow=())
    spec = p1_spec(world["sha"])
    spec[key][tier] = entry
    write_json(world["repo"], et.P1_SPEC, spec)
    with pytest.raises(et.PreflightError, match=match):
        et.look_config("T1")
    assert world["calls"] == {"train": 0, "validation": 0, "test": 0}


def test_calls_file_must_be_boolean_and_cover_arms(world, monkeypatch):
    rel = "results/m6/tmp/test_calls.parquet"
    df = pd.read_parquet(world["repo"] / rel)
    spec = p1_spec(world["sha"])
    for bad, match in ((df.assign(p57=df.p57.astype(int)), "boolean"), (df[df.method == "P1[H1]"], "exactly the arms"),
                       (df.assign(score=0.1), "columns")):
        spec["calls"]["test"]["sha256"] = _parquet(world["repo"], rel, bad)
        write_json(world["repo"], et.P1_SPEC, spec)
        with pytest.raises(et.PreflightError, match=match):
            et.look_config("T1")


def test_t1_refused_before_rehearsal(world, monkeypatch, capsys):
    _patch_labels(world, monkeypatch, allow=())
    _patch_main(world, monkeypatch)
    assert et.main(["--look", "T1"]) == 2
    assert "rehearsal" in capsys.readouterr().err
    assert not (world["repo"] / "results" / "test_ledger.tsv").exists()


def _complete_t1_rows():
    from nextgem import ledger

    ledger.append("test", "started", look_id="T1", workstream="M6", config_sha="c", commit=False)
    ledger.append("test", "completed", look_id="T1", commit=False)


def test_t2_refused_before_w7_stop_row(world, monkeypatch, capsys):
    from nextgem import ledger

    _patch_labels(world, monkeypatch, allow=())
    _patch_main(world, monkeypatch)
    _complete_t1_rows()
    toy.commit_all(world["repo"], "T1 rows")
    assert et.main(["--look", "T2"]) == 2
    err = capsys.readouterr().err
    assert "W7 stop" in err and "T2 refused" in err
    assert any("W7 stop" in p for p in et.state_problems("T2"))
    ledger.append("validation", "non_look", look_id=et.W7_STOP, workstream="W7",
                  detail=json.dumps({"p2_spec_sha256": "0" * 64}), commit=False)
    toy.commit_all(world["repo"], "W7 stop naming another P2")
    assert any("does not name the current P2 spec" in p for p in et.state_problems("T2"))
    ledger.append("validation", "non_look", look_id=et.W7_STOP, workstream="W7",
                  detail=json.dumps({"p2_spec_sha256": et._sha(world["repo"] / et.P2_SPEC)}), commit=False)
    toy.commit_all(world["repo"], "W7 stop")
    assert et.state_problems("T2") == []
    assert [r[3] for r in _committed_rows(world["repo"], "test")] == ["started", "completed"]
    assert world["calls"] == {"train": 0, "validation": 0, "test": 0}


def test_unseal_refused_before_t2(world, monkeypatch, capsys):
    from nextgem import ledger

    _patch_labels(world, monkeypatch, allow=())
    _patch_main(world, monkeypatch)
    assert et.main(["--unseal"]) == 2
    assert "T2 is not completed" in capsys.readouterr().err
    _complete_t1_rows()
    toy.commit_all(world["repo"], "T1 rows")
    assert et.main(["--unseal"]) == 2
    with pytest.raises(ledger.LedgerError, match="completed T2"):
        ledger.append("test", "unseal", look_id="unseal", commit=False)
    assert not any(r["event"] == "unseal" for r in ledger.read("test"))
    assert world["calls"] == {"train": 0, "validation": 0, "test": 0}


def test_third_look_needs_authorisation(world, monkeypatch):
    _complete_t1_rows()
    toy.commit_all(world["repo"], "T1 rows")
    probs = et.state_problems("T3")
    assert any("test_third_look.md" in p for p in probs)


# ============================================================================= rehearsal, T1, sealing, stdout
def test_rehearsal_validation_stage_writes_validation_ledger_rows(world, monkeypatch):
    from nextgem import ledger

    cfg, res = _rehearse(world, monkeypatch)
    rows = [r for r in ledger.read("validation") if r["look_id"] == REH]
    assert [r["event"] for r in rows] == ["started", "output", "output", "completed"]
    assert rows[0]["workstream"] == "M6" and rows[0]["config_sha"] == et.config_sha(cfg)
    d = json.loads(rows[0]["detail"])
    assert d["tier"] == "validation" and d["fit_tiers"] == ["train"] and d["spec_sha"] == et._sha(world["repo"] / et.P1_SPEC)
    assert [json.loads(r["detail"])["kind"] for r in rows[1:3]] == ["revealed", "sealed"]
    assert not ledger.verify("validation") and not ledger.read("test")
    assert world["calls"]["validation"] == 1 and world["calls"]["test"] == 0
    assert all(p.startswith(str(paths.sealed_dir(REH))) and c == "sealed" for p, c, _ in world["sealed"])
    out = pd.read_csv(res["revealed"]["files"][0], sep="\t")
    assert {"confirmatory", "transfer", "estimates", "binary_6.8", "failed_builds", "counts"} == set(out.section)
    tr = out[out.section == "transfer"].iloc[0]
    assert tr.n_lineages == 4 and tr.k_required == 3  # k-of-K on the 4 validation lineages
    assert et.state_problems("T1") == []


def test_t1_seals_only_under_sealed_dir_and_ledger_has_hashes_only(world, monkeypatch):
    from nextgem import labels, ledger

    _rehearse(world, monkeypatch)
    world["sealed"].clear()
    _patch_labels(world, monkeypatch, check_started=("test", "T1"))
    cfg = et.look_config("T1")
    ctx = world["ctx_for"]("test")
    fit = pd.concat([labels.train(columns=et.ev.LABEL_COLS), labels.validation()], ignore_index=True)
    prep = et.prepare(cfg, ctx, fit, "test")
    out_dir = world["repo"] / et.REVEAL_DIR
    et.write_config(cfg, out_dir, "test_look")
    before = {p for p in world["repo"].rglob("*") if p.is_file() and ".git" not in p.parts}
    res = et.run_p1_look(cfg, ctx=ctx, prep=prep, load_labels=labels.test, out_dir=out_dir, ledger_name="test", author="test")
    assert world["calls"]["test"] == 1
    sealed_root = str(paths.sealed_dir("T1"))
    assert sealed_root.startswith(str(world["data"] / "locked" / "sealed"))
    assert world["sealed"] and all(p.startswith(sealed_root + "/") and c == "sealed" for p, c, _ in world["sealed"])
    names = {p.rsplit("/", 1)[1] for p, _, _ in world["sealed"]}
    assert names == {f"{k}.tsv" for k in ("per_lineage", "per_gene", "w1_benchmark", "generic_all_lines", "strata", "attribution",
                                           "binary_threshold_grid", "per_task", "failure_account", "worst_cells", "sensitivities")}
    new = {p for p in world["repo"].rglob("*") if p.is_file() and ".git" not in p.parts} - before
    assert {p.name for p in new} == {"test_look_T1.tsv", "test_look_T1_gates.json", "test_ledger.tsv"}
    assert not any(n in {p.name for p in new} for n in names)
    rows = ledger.read("test")
    assert [r["event"] for r in rows] == ["started", "output", "output", "completed"]
    sealed_row = json.loads(rows[2]["detail"])
    assert sealed_row["kind"] == "sealed" and sealed_row["files"] == {p.rsplit("/", 1)[1]: s for p, _, s in world["sealed"]}
    assert rows[2]["output_sha"] == et.manifest_sha(sealed_row["files"]) == res["sealed"]["manifest_sha256"]
    assert all(et.SHA_RE.match(v) for v in sealed_row["files"].values())
    ledger_text = (world["repo"] / "results/test_ledger.tsv").read_text()
    ids = set(world["lineage"].index)
    assert not any(i in ledger_text for i in ids) and not any(f"\t{g}\t" in ledger_text for g in GENES)
    assert not any(k in ledger_text for k in TEST)  # no test lineage name in the ledger
    # revealed outputs: scalar-only; lineage names only in failed_builds
    out = pd.read_csv(out_dir / "test_look_T1.tsv", sep="\t")
    et.assert_revealed_scalar_only(out.assign(**{c: out[c].astype(object) for c in out.columns if c in et.STR_COLS}), ctx)
    assert set(out.loc[out.lineage.notna(), "section"]) == {"failed_builds"}
    fb = out[(out.section == "failed_builds") & (out.method == "P1[H2]") & (out.lineage == "all")].iloc[0]
    assert fb.n_failed == 1 and fb.n_lines == len(ctx["tier_lines"])
    gates = json.loads((out_dir / "test_look_T1_gates.json").read_text())
    assert [g["contrast"] for g in gates["gates"]] == ["C1", "C1b", "C3", "C4", "C5"]
    assert gates["headline"]["contrast"] == "C1" and "7 test lineages" in gates["headline"]["conditional_on"]
    assert json.loads(rows[1]["detail"])["headline"] == gates["headline"]


@pytest.mark.parametrize("bad,match", [
    pytest.param({"name": "ACH-000600"}, "line ID", id="line_id_cell"),
    pytest.param({"method": "1005"}, "gene", id="entrez_cell"),
    pytest.param({"lineage": "KA"}, "only in the failed_builds", id="lineage_outside_failed_builds"),
    pytest.param({"tap_line_3": 0.1}, "not a known scalar", id="per_unit_column"),
])
def test_revealed_guard_rejects_unit_keys(world, bad, match):
    ctx = world["ctx_for"]("test")
    out = pd.DataFrame([{"section": "estimates", "name": "C2:B2", "method": "P1[H1]", "estimate": 0.1}])
    out = out.astype({"section": object, "name": object, "method": object})
    et.assert_revealed_scalar_only(out, ctx)
    with pytest.raises(RuntimeError, match=match):
        et.assert_revealed_scalar_only(out.assign(**bad), ctx)


@pytest.mark.parametrize("name,kind", [pytest.param("p1_spec.template.json", "p1", id="p1_template"),
                                       pytest.param("p2_spec.template.json", "p2", id="p2_template")])
def test_spec_templates_valid_after_removing_doc(name, kind):
    """The committed templates document every field in ``_doc``; the checker refuses ``_doc`` and accepts the rest."""
    check = et.check_p1_spec if kind == "p1" else et.check_p2_spec
    doc = json.loads((paths.PROJECT_ROOT / "results" / "m6" / name).read_text())
    with pytest.raises(et.PreflightError, match="keys"):
        check(doc)
    check({k: v for k, v in doc.items() if k != "_doc"})
    keys = et.P1_KEYS | et.P1_OPTIONAL if kind == "p1" else et.P2_KEYS | et.P2_OPTIONAL
    assert set(doc) - {"_doc"} <= keys
    if kind == "p1":
        assert keys <= set(doc["_doc"])


def test_t1_stdout_is_scalar_only(world, monkeypatch, capsys):
    _rehearse(world, monkeypatch)
    _patch_labels(world, monkeypatch)
    _patch_main(world, monkeypatch)
    world["sealed"].clear()
    capsys.readouterr()
    assert et.main(["--look", "T1"]) == 0
    text = capsys.readouterr().out
    doc = json.loads(text)
    et.assert_stdout_scalar_only(text, world["ctx_for"]("test"))
    assert "ACH-" not in text and not any(f'"{k}"' in text for k in TEST)
    assert not any(f'"{g}"' in text or f" {g}," in text for g in GENES)
    assert doc["look"] == "T1" and set(doc["sealed"]["files"]) == {p.rsplit("/", 1)[1] for p, _, _ in world["sealed"]}
    assert all(et.SHA_RE.match(v) for v in doc["sealed"]["files"].values())
    assert {g["contrast"] for g in doc["gates"]} == {"C1", "C1b", "C3", "C4", "C5"}
    with pytest.raises(RuntimeError):
        et.assert_stdout_scalar_only(json.dumps({"x": "ACH-000600"}), world["ctx_for"]("test"))
    with pytest.raises(RuntimeError):
        et.assert_stdout_scalar_only(json.dumps({"x": "KA"}), world["ctx_for"]("test"))


def test_chain_t1_w7stop_t2_unseal(world, monkeypatch, capsys):
    """Rehearsal → T1 → W7 stop → T2 (C4/C5 with the T1 chain) → unseal: every sealed file re-derived with the same
    harness and published only if its SHA-256 equals the ledger's; the ledger has exactly one unseal row."""
    from nextgem import ledger

    _rehearse(world, monkeypatch)
    _patch_labels(world, monkeypatch)
    _patch_main(world, monkeypatch)
    assert et.main(["--look", "T1"]) == 0
    toy.commit_all(world["repo"], "T1 outputs")
    assert et.main(["--look", "T2"]) == 2  # no W7 stop row yet
    ledger.append("validation", "non_look", look_id=et.W7_STOP, workstream="W7",
                  detail=json.dumps({"p2_spec_sha256": et._sha(world["repo"] / et.P2_SPEC)}))
    capsys.readouterr()
    assert et.main(["--look", "T2"]) == 0
    t2 = json.loads(capsys.readouterr().out)
    t1g = json.loads((world["repo"] / et.REVEAL_DIR / "test_look_T1_gates.json").read_text())
    assert t2["headline"] == t1g["headline"]
    assert [g["contrast"] for g in t2["gates"]] == ["C1", "C1b", "C3", "C4", "C5"]
    assert [g["status"] for g in t2["gates"]][:3] == [g["status"] for g in t1g["gates"]][:3]
    assert "pending" not in {g["status"] for g in t2["gates"]}
    grid = pd.read_csv(world["repo"] / et.REVEAL_DIR / "test_look_T2_grid.tsv", sep="\t")
    assert set(grid.method) == {"H1", "H2", "P2", "B6[L_prob]"}
    toy.commit_all(world["repo"], "T2 outputs")
    assert et.main(["--unseal"]) == 0
    rep = json.loads(capsys.readouterr().out)["unseal"]
    assert set(rep) == {"T1"} and not rep["T1"]["sha_mismatch_not_published"] and len(rep["T1"]["published"]) == 11
    pub = world["repo"] / et.UNSEALED_DIR / "T1"
    want = json.loads([r for r in ledger.read("test") if r["look_id"] == "T1" and r["event"] == "output"][1]["detail"])["files"]
    assert {p.name: et._sha(p) for p in pub.iterdir()} == want
    assert et.main(["--unseal"]) == 0  # idempotent: the row is not repeated
    assert sum(r["event"] == "unseal" for r in ledger.read("test")) == 1
    assert not ledger.verify("test")


# ============================================================================= dry modes read no evaluation label
@pytest.mark.parametrize("argv", [
    pytest.param(["--print-config"], id="print_config_T1"),
    pytest.param(["--look", "T2", "--print-config"], id="print_config_T2"),
    pytest.param(["--rehearsal", "--print-config"], id="print_config_rehearsal"),
    pytest.param(["--unseal", "--print-config"], id="print_config_unseal"),
    pytest.param(["--print-frozen-candidates"], id="print_frozen_candidates"),
    pytest.param(["--look", "T1", "--preflight-only"], id="preflight_only_T1"),
    pytest.param(["--rehearsal", "--preflight-only"], id="preflight_only_rehearsal"),
    pytest.param(["--dry-run"], id="dry_run_T1_test_tier"),
    pytest.param(["--dry-run", "--tier", "validation"], id="dry_run_T1_validation_tier"),
    pytest.param(["--look", "T2", "--dry-run"], id="dry_run_T2"),
])
def test_dry_modes_read_no_validation_or_test_labels(world, monkeypatch, capsys, argv):
    from nextgem import ledger

    _patch_labels(world, monkeypatch, allow=("train",))
    monkeypatch.setattr(et, "real_context", lambda cfg, tier: world["ctx_for"](tier))
    before = {p for p in world["repo"].rglob("*") if p.is_file() and ".git" not in p.parts}
    rc = et.main(argv)
    out = capsys.readouterr().out
    assert rc in (0, 2)
    if "--dry-run" in argv:
        assert rc == 0
        doc = json.loads(out)
        assert doc["scalar_guard"] == "pass" and not doc["test_labels_read"] and not doc["validation_labels_read"]
        assert doc["ledger_rows_written"] == 0 and doc["files_written"] == 0
        et.assert_stdout_scalar_only(out, world["ctx_for"]("test"))
    assert world["calls"]["validation"] == 0 and world["calls"]["test"] == 0
    assert world["sealed"] == [] and ledger.read("test") == [] and ledger.read("validation") == []
    after = {p for p in world["repo"].rglob("*") if p.is_file() and ".git" not in p.parts}
    assert {p.relative_to(world["repo"]).parts[0] for p in after - before} <= {"logs"}


def test_spec_override_refused_in_real_modes(world, monkeypatch, capsys):
    _patch_labels(world, monkeypatch, allow=())
    assert et.main(["--look", "T1", "--spec", et.P1_SPEC]) == 2
    assert "--spec" in capsys.readouterr().err


def test_real_constants_are_enforced(world):
    cfg = et.look_config("T1")
    assert et.check_constants(cfg) == []  # the patched constants are self-consistent
    probs = et.check_constants(dict(cfg, B_conf=999, seed=12))
    assert any(p.startswith("B_conf") for p in probs) and any(p.startswith("seed") for p in probs)
