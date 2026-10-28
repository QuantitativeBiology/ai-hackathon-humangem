"""Toy / known-answer tests for the W3 library (src/nextgem/w3.py). Label-free, GLPK, fast."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

cobra = pytest.importorskip("cobra")


def _toy():
    """A_e -> A_c (no GPR) ; A_c -> B_c by g1 ; A_c -> B_c by (g2 or g3) bounded 2 ; B_c -> biomass ; P_c from B_c by g4."""
    cobra.Configuration().solver = "glpk"
    m = cobra.Model("toy")
    A_e, A, B, P, X = (cobra.Metabolite(i, name=n, compartment=c) for i, n, c in
                       (("A_e", "A", "e"), ("A_c", "A", "c"), ("B_c", "B", "c"), ("P_c", "P", "c"), ("X_c", "X", "c")))
    ex = cobra.Reaction("EX_A", lower_bound=-10, upper_bound=1000)
    ex.add_metabolites({A_e: -1})
    tr = cobra.Reaction("TR_A", lower_bound=0, upper_bound=1000)
    tr.add_metabolites({A_e: -1, A: 1})
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=6)
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "g1"
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=2)
    r2.add_metabolites({A: -1, B: 1})
    r2.gene_reaction_rule = "g2 or g3"
    r3 = cobra.Reaction("R3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({B: -1, P: 1})
    r3.gene_reaction_rule = "g4"
    rx = cobra.Reaction("RX", lower_bound=0, upper_bound=1000)
    rx.add_metabolites({A: -1, X: 1})
    rx.gene_reaction_rule = "g5"
    bio = cobra.Reaction("MAR13082", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({B: -1, X: -0.1})
    m.add_reactions([ex, tr, r1, r2, r3, rx, bio])
    m.objective = "MAR13082"
    return m


def _task(outputs, inputs=(("A[e]", 0.0, 1000.0),), equations=()):
    from raven_toolbox.tasks.tasklist import Task

    return Task(id="BS", description="toy", should_fail=False, print_fluxes=False, comments="", inputs=list(inputs),
                outputs=list(outputs), equations=list(equations), changed=[])


def _setup(m):
    from raven_toolbox.tasks.check import task_name_maps

    from nextgem import w3

    base = w3.task_base(m)
    maps = task_name_maps(base)
    saved = {x.id: (base.constraints[x.id].lb, base.constraints[x.id].ub) for x in base.metabolites}
    gd = w3.gene_disabled_reactions(base)
    return base, maps, saved, gd, {r.id for r in base.reactions}


def test_w3_task_feasibility_and_capacity_toy():
    from nextgem import w3

    base, maps, saved, gd, orig = _setup(_toy())
    kos = list(w3.ko_groups(gd))
    r = w3.task_capacity(base, _task([("B[c]", 1.0, 1000.0)]), maps, saved, kos, orig)
    assert r["status"] == "ok"
    lam = r["lam"]
    assert r["lam_wt"] == pytest.approx(8.0)  # R1 (6) + R2 (2), capped by uptake 10
    assert r["feasible"][frozenset({"R1"})] is True  # R2 still makes B
    assert w3.capacity_ratio(lam[frozenset({"R1"})], r["lam_wt"]) == pytest.approx(0.25)
    assert all(frozenset({"R2"}) != k for k in kos)  # g2 / g3 are isozymes: no single-gene KO closes R2
    # P needs R3 (g4): KO infeasible, capacity 0
    r = w3.task_capacity(base, _task([("P[c]", 1.0, 1000.0)]), maps, saved, kos, orig)
    assert r["feasible"][frozenset({"R3"})] is False
    assert w3.capacity_ratio(r["lam"][frozenset({"R3"})], r["lam_wt"]) == 0.0
    # constraint bounds restored after every task
    assert all((base.constraints[k].lb, base.constraints[k].ub) == v for k, v in saved.items())


def test_w3_task_essential_reactions_toy():
    from nextgem import w3

    base, maps, saved, gd, orig = _setup(_toy())
    r = w3.task_capacity(base, _task([("P[c]", 1.0, 1000.0)]), maps, saved, [], orig, essential_reactions=True)
    assert "R3" in r["ess_rxns"] and "TR_A" in r["ess_rxns"]
    assert "R1" not in r["ess_rxns"]  # R2 can replace R1 at lb 1


def test_w3_variant_apply_is_reverted():
    from nextgem import w3

    m = _toy()
    before = (m.reactions.MAR13082.reaction, len(m.reactions), m.reactions.R3.bounds)
    vs = [w3.Variant("rm", "remove", {"mets": ["X_c"]}), w3.Variant("add", "add", {"stoich": {"P_c": -1e-4}}),
          w3.Variant("fix", "fix", {"fluxes": {"R3": 0.5}}),
          w3.Variant("dm", "demands", {"demands": [("DM_P", {"P_c": -1.0}, 0.1)]})]
    for v in vs + [w3.Variant("combo", "combo", {"variants": vs[:2]})]:
        with m:
            v.apply(m)
        assert (m.reactions.MAR13082.reaction, len(m.reactions), m.reactions.R3.bounds) == before
        assert w3.Variant.from_json(v.to_json()).to_json() == v.to_json()


def test_w3_ko_sweep_shortcut_equals_brute_force():
    from nextgem import w3

    m = _toy()
    genes = sorted(g.id for g in m.genes)
    clos = w3.gene_closures(m)
    bounds = {"EX_A": 10.0}
    for v in (w3.Variant("T0"), w3.Variant("rm", "remove", {"mets": ["X_c"]}),
              w3.Variant("dm", "demands", {"demands": [("DM_P", {"P_c": -1.0}, 0.5)]})):
        r = w3.ko_sweep(m, genes, clos, bounds, v, verify_n=5)
        with m:
            m.medium = bounds
            v.apply(m)
            wt = m.slim_optimize()
            brute = []
            for g in genes:
                with m:
                    m.genes.get_by_id(g).knock_out()
                    x = m.slim_optimize(error_value=0.0)
                brute.append(np.clip(np.round(x / wt, 6), 0, 1))
        assert np.allclose(r["ratio_r"], brute), (v.name, r["ratio_r"], brute)
    # X is needed for biomass (g5 lethal) unless X is left out of the biomass
    r = w3.ko_sweep(m, genes, clos, bounds, w3.Variant("T0"))
    assert r["ratio_r"][genes.index("g5")] == 0.0
    r = w3.ko_sweep(m, genes, clos, bounds, w3.Variant("rm", "remove", {"mets": ["X_c"]}))
    assert r["ratio_r"][genes.index("g5")] == 1.0


def test_w3_pseudo_task_is_deterministic_and_producible():
    from nextgem import w3

    base, maps, saved, gd, orig = _setup(_toy())
    comp = w3.compartment_tokens(base)
    t = _task([("B[c]", 1.0, 1000.0)])
    a = w3.make_pseudo_task(t, base, maps, saved, np.random.default_rng([11, 0]), comp)
    b = w3.make_pseudo_task(t, base, maps, saved, np.random.default_rng([11, 0]), comp)
    assert a is not None and a[1] == b[1]
    new, swaps, _ = a
    assert swaps[0][0] == "B[c]" and swaps[0][1] != "B[c]" and swaps[0][1].endswith("[c]")
    r = w3.task_capacity(base, new, maps, saved, [], orig)
    assert r["status"] == "ok"


def test_w3_cellfie_rule_and_gpr_score():
    from nextgem import w3

    expr = pd.DataFrame(np.log2(np.array([[1.0, 10.0, 100.0], [3.0, 30.0, 300.0]]) + 1), index=["l1", "l2"],
                        columns=[1, 2, 3])
    thr = w3.cellfie_thresholds(expr, ["l1", "l2"])
    lin = np.array([[1.0, 10.0, 100.0], [3.0, 30.0, 300.0]]).ravel()
    lo, hi = np.quantile(lin, 0.25), np.quantile(lin, 0.75)
    assert thr.min() >= lo - 1e-9 and thr.max() <= hi + 1e-9
    gs = w3.cellfie_gene_scores(expr, ["l1"], thr)
    assert gs.loc["l1", 2] == pytest.approx(5 * np.log1p(10.0 / thr[2]))
    assert w3.eval_gpr_score(("and", ("a", ("or", ("b", "c")))), {"a": 2.0, "b": 1.0, "c": 5.0}) == 2.0
    assert w3.eval_gpr_score(("or", ("a", "b")), {"a": 2.0}) == 2.0
    assert np.isnan(w3.eval_gpr_score("z", {}))


def test_w3_mode_e_rules():
    from nextgem import w3

    r = pd.Series({"a": 1.0, "b": 0.5, "c": 0.0, "d": np.nan})
    e = w3.mode_e(r, pd.Series({"a": True}))
    assert e.to_dict() == {"a": 1.0, "b": 0.5, "c": 1.0, "d": 0.0}


def test_w3_essential_task_groups():
    from nextgem import w3

    for p in ("human1", "human2"):
        t = w3.load_tasks(p, "essential")
        assert len(t) == 57
        assert t.group.nunique() == 11  # DESIGN: 11 functional groups, one of which is "Growth" (GR only)
        assert t.group[t.category != "GR"].nunique() == 10
        assert t.category.value_counts().to_dict() == {"BS": 18, "SU": 17, "IC": 12, "ER": 9, "GR": 1}


def test_w3_variant_pool_remove():
    """pool_remove drops species from the named pool reaction (default MAR10065), is reverted, and changes the KO."""
    from nextgem import w3

    m = _toy()
    Y = cobra.Metabolite("Y_c", name="Y", compartment="c")
    pool_y = cobra.Metabolite("pool_c", name="pool", compartment="c")
    ry = cobra.Reaction("RY", lower_bound=0, upper_bound=1000)
    ry.add_metabolites({m.metabolites.A_c: -1, Y: 1})
    ry.gene_reaction_rule = "g6"
    pool = cobra.Reaction("MAR10065", lower_bound=0, upper_bound=1000)
    pool.add_metabolites({Y: -1, m.metabolites.X_c: -1, pool_y: 1})
    m.add_reactions([ry, pool])
    m.reactions.MAR13082.add_metabolites({pool_y: -0.1})
    before = m.reactions.MAR10065.reaction
    v = w3.Variant("B-loo:Y", "pool_remove", {"mets": ["Y_c"]})
    with m:
        v.apply(m)
        assert "Y_c" not in {x.id for x in m.reactions.MAR10065.metabolites}
        assert "X_c" in {x.id for x in m.reactions.MAR10065.metabolites}
    assert m.reactions.MAR10065.reaction == before
    genes = sorted(g.id for g in m.genes)
    clos = w3.gene_closures(m)
    r0 = w3.ko_sweep(m, genes, clos, {"EX_A": 10.0}, w3.Variant("T0"))
    r1 = w3.ko_sweep(m, genes, clos, {"EX_A": 10.0}, v)
    assert r0["ratio_r"][genes.index("g6")] == 0.0 and r1["ratio_r"][genes.index("g6")] == 1.0
    # explicit pool id (as the MAR10064 dihydroorotate component) and JSON round trip
    v2 = w3.Variant("B-loo:Y2", "pool_remove", {"mets": ["Y_c"], "pool": "MAR10065"})
    assert w3.Variant.from_json(v2.to_json()).params["pool"] == "MAR10065"
    with m:
        v2.apply(m)
        assert "Y_c" not in {x.id for x in m.reactions.MAR10065.metabolites}


def _score_module():
    import importlib
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    return importlib.import_module("w3_score_train")


def _toy_G(st):
    """Two folds (lineages k1, k2), 2 lines each, 6 genes; T0 calls genes 1 and 2; mode adds 2 calls and removes 1."""
    U = [1, 2, 3, 4, 5, 6]
    lines = ["ACH-000001", "ACH-000002", "ACH-000003", "ACH-000004"]
    t0 = pd.Series([1.0, 1.0, 0.5, 0.0, 0.0, 0.2], index=U)
    mode = pd.Series([1.0, 0.0, 1.0, 1.0, 0.0, 0.2], index=U)  # adds 3, 4; removes 2
    st.G.clear()
    st.G.update({"U": U, "E": lines, "modes": {("human1", "T0"): t0, ("human1", "M"): mode},
                 "extra": {"human1": {"in_p": pd.Series([True, True, True, True, True, False], index=U)}},
                 "fold_of": {"ACH-000001": "k1", "ACH-000002": "k1", "ACH-000003": "k2", "ACH-000004": "k2"},
                 "prior": {"k1": pd.Series({1: 0.9, 2: 0.1, 3: 0.2, 4: 0.3, 5: 0.8, 6: 0.99}),
                           "k2": pd.Series({1: 0.05, 2: 0.9, 3: 0.7, 4: 0.1, 5: 0.6, 6: 0.99})},
                 "lineage": pd.Series({"ACH-000001": "k1", "ACH-000002": "k1", "ACH-000003": "k2", "ACH-000004": "k2"}),
                 "folds": ["k1", "k2"]})
    return U, lines


def test_w3_mcr_matrix_literal_adds_top_b2_only():
    st = _score_module()
    U, lines = _toy_G(st)
    M = st.mcr_matrix("human1", "M", symmetric=False)
    # 2 added calls go to the uncalled in-parent genes with the highest fold-scoped B2; gene 6 (not in parent) never;
    # removed calls are NOT mirrored (literal DESIGN control)
    assert M.loc["ACH-000001"].tolist() == [1.0, 1.0, 0.5, 1.0, 1.0, 0.2]  # k1: B2 5 (0.8), 4 (0.3)
    assert M.loc["ACH-000003"].tolist() == [1.0, 1.0, 1.0, 0.0, 1.0, 0.2]  # k2: B2 3 (0.7), 5 (0.6)
    # a mode that only removes calls gives the identity (= T0): the literal control is vacuous for removals (A11a)
    st.G["modes"][("human1", "R")] = pd.Series([1.0, 0.0, 0.5, 0.0, 0.0, 0.2], index=U)
    R = st.mcr_matrix("human1", "R", symmetric=False)
    assert (R.to_numpy() == np.tile(st.G["modes"][("human1", "T0")].to_numpy(), (4, 1))).all()


def test_w3_mcr_matrix_symmetric_removes_lowest_b2():
    st = _score_module()
    U, lines = _toy_G(st)
    M = st.mcr_matrix("human1", "M", symmetric=True)
    # k1: one removal from the called genes {1, 2} with the lowest B2 -> gene 2 (0.1); k2 -> gene 1 (0.05)
    assert M.loc["ACH-000002"].tolist() == [1.0, 0.0, 0.5, 1.0, 1.0, 0.2]
    assert M.loc["ACH-000004"].tolist() == [0.0, 1.0, 1.0, 0.0, 1.0, 0.2]


def test_w3_fold_stats_min_lines_rule():
    """A fold counts as positive only with >= 5 estimable lines and mean delta > 0 (DESIGN §9 TG)."""
    st = _score_module()
    lines = [f"ACH-{i:06d}" for i in range(1, 12)]
    lin = pd.Series(["A"] * 6 + ["B"] * 4 + ["C"], index=lines)
    st.G.clear()
    st.G.update({"lineage": lin, "folds": ["A", "B", "C", "D"]})

    def table(taps):
        return pd.DataFrame({"n": 10, "n_pos": 5, "estimable": True, "tap": taps}, index=lines)

    a = table([0.6] * 6 + [0.9] * 4 + [0.9])
    b = table([0.5] * 6 + [0.1] * 4 + [0.1])
    r = st.fold_stats(a, b)
    assert r["folds_positive"] == 1 and r["n_folds"] == 4  # B (4 lines) and C (1) are "not positive"; D is absent
    a2 = table([0.4] * 6 + [0.9] * 5)
    a2.loc[lines[-1], "tap"] = 0.9
    assert st.fold_stats(a2, b)["folds_positive"] == 0  # A negative; B, C too small
