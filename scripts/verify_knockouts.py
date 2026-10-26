"""Verify the single-gene knockout pipeline (nextgem.knockouts) before trusting genome-scale results.

Checks (each prints PASS/FAIL; exit code 1 if any fails; table -> results/baseline/verification.tsv):
  toy_*      toy network with hand-derived answers: isozymes (OR), complex (AND), nested rules,
             redundant routes (exact ratio 1/3), dead-end and blocked-branch genes, a KO that makes
             the LP infeasible, swapped medium and objective, blocked objective, serial == parallel,
             status policy (non-optimal solves raise; NaN never silently non-essential), ratio_r
  gpr_*      on Human1 and Human2, reactions closed by each KO == reactions whose GPR is False under
             an independent evaluator (Python ast), for all genes
  human*_*   positive controls (WT growth, CAD/UMPS lethal on Ham's), FVA consistency for genes
             whose KO closes exactly one reaction (ratio == 1 iff that reaction's FVA range at
             optimal growth contains 0), a named non-lethal gene, negative controls (unknown IDs
             rejected), serial == parallel (bitwise on ratio_r), task file derived from the model's
             own release, genome-scale infeasible KO (biomass forced >= 1, lethal KO -> status
             infeasible, ratio 0)
  solver_*   Gurobi vs GLPK on 200 random Human2 KOs (seed 0)

Usage: .venv/bin/python scripts/verify_knockouts.py [--skip-glpk] [--processes N]
"""

import argparse
import ast
import math
import random
import sys
import time
import warnings

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

import nextgem.knockouts as ko_module
from nextgem import GEMS, NonOptimalError, essential_genes, gem_of, load_model, single_gene_knockouts
from nextgem.genes import read_genes_tsv
from nextgem.knockouts import check_genes, configure, rounded_ratio
from nextgem.model import PROJECT_ROOT, tasks_file_for

OUT = PROJECT_ROOT / "results" / "baseline" / "verification.tsv"
SEED = 0
BIOMASS = "MAR13082"
RESULTS: list[tuple[str, bool, str, str]] = []

# Control failures that are properties of a model, not of the pipeline. The check still runs and is
# reported as failed ("XFAIL" if it fails, "XPASS" if it unexpectedly passes); it does not fail the script.
KNOWN_MODEL_ISSUES = {
    "human1_lethal_UMPS": (
        "Human1 v1.19.0 has thymidylate synthase MAR04644 reversible, so thymidine (in Ham's) -> dTMP -> dUMP"
        " -> deoxyuridine -> uracil -> UMP (UPRT) bypasses UMPS; Human2 has MAR04644 irreversible."
        " See check human1_UMPS_lethal_if_TYMS_irreversible."
    ),
}


def check(name: str, passed: bool, detail: str = "") -> None:
    issue = KNOWN_MODEL_ISSUES.get(name, "")
    label = ("XPASS" if passed else "XFAIL") if issue else ("PASS" if passed else "FAIL")
    RESULTS.append((name, bool(passed), detail, issue))
    print(f"{label}  {name}  {detail}" + (f"  [known model issue: {issue}]" if issue else ""), flush=True)


# ---------------------------------------------------------------- toy network
# a_e -T_a-> a_c -R1-> b_c -R2-> c_c -R3-> d_c -R4-> e_c ; BIOMASS: e_c + f_c ->
# f_c from uptake (T_f) or synthesis 2 a_c -> f_c (R_fsyn); forced maintenance a_c -> m_c (lb 1)
# whose product must leave via R_mx; R_side: a_c -> w_c (secreted); dead end R_dead: c_c -> z_c;
# blocked branch R_blk: q_c -> e_c (q_c never produced). Uptake: a <= 10, f <= 10.
# WT growth: e needs 1 a/unit, maintenance 1 a -> growth 9. Without f uptake: x + 2x + 1 = 10 -> 3.
TOY_REACTIONS = {
    # id: (stoichiometry, lb, ub, GPR)
    "EX_a": ({"a_e": -1}, -10, 1000, ""),
    "EX_f": ({"f_e": -1}, -10, 1000, ""),
    "EX_w": ({"w_e": -1}, 0, 1000, ""),
    "EX_m": ({"m_e": -1}, 0, 1000, ""),
    "T_a": ({"a_e": -1, "a_c": 1}, 0, 1000, "g_ta"),
    "T_f": ({"f_e": -1, "f_c": 1}, 0, 1000, "g_fup"),
    "R1": ({"a_c": -1, "b_c": 1}, 0, 1000, "g1a or g1b"),
    "R2": ({"b_c": -1, "c_c": 1}, 0, 1000, "g2a and g2b"),
    "R3": ({"c_c": -1, "d_c": 1}, 0, 1000, "(g3a and g3b) or g3c"),
    "R4": ({"d_c": -1, "e_c": 1}, 0, 1000, "(g4a or g4b) and g4c"),
    "R_fsyn": ({"a_c": -2, "f_c": 1}, 0, 1000, "g_fsyn"),
    "R_maint": ({"a_c": -1, "m_c": 1}, 1, 1000, ""),
    "R_mx": ({"m_c": -1, "m_e": 1}, 0, 1000, "g_mx"),
    "R_side": ({"a_c": -1, "w_c": 1}, 0, 1000, "g_side"),
    "R_wx": ({"w_c": -1, "w_e": 1}, 0, 1000, ""),
    "R_dead": ({"c_c": -1, "z_c": 1}, 0, 1000, "g_dead"),
    "R_blk": ({"q_c": -1, "e_c": 1}, 0, 1000, "g_blk"),
    "BIOMASS": ({"e_c": -1, "f_c": -1}, 0, 1000, ""),
}
TOY_MEDIUM = {"EX_a": 10.0, "EX_f": 10.0}
TOY_WT = 9.0
TOY_ESSENTIAL = {"g_ta", "g2a", "g2b", "g4c", "g_mx"}
TOY_RATIOS = {"g_fup": 3.0 / 9.0, "g_fsyn": 1.0, "g1a": 1.0, "g3c": 1.0, "g_dead": 1.0, "g_blk": 1.0}
TOY_STATUS = {"g_mx": "infeasible"}  # KO closes the only outlet of a forced flux
TOY_CLOSED = {  # reactions closed by each single-gene KO
    "g_ta": {"T_a"}, "g_fup": {"T_f"}, "g1a": set(), "g1b": set(), "g2a": {"R2"}, "g2b": {"R2"},
    "g3a": set(), "g3b": set(), "g3c": set(), "g4a": set(), "g4b": set(), "g4c": {"R4"},
    "g_fsyn": {"R_fsyn"}, "g_mx": {"R_mx"}, "g_side": {"R_side"}, "g_dead": {"R_dead"}, "g_blk": {"R_blk"},
}


def toy_model(solver: str = "glpk") -> cobra.Model:
    model = cobra.Model("toy")
    model.solver = solver
    mets = {m: cobra.Metabolite(m, compartment=m[-1]) for s, *_ in TOY_REACTIONS.values() for m in s}
    for rid, (stoich, lb, ub, gpr) in TOY_REACTIONS.items():
        rxn = cobra.Reaction(rid, lower_bound=lb, upper_bound=ub)
        model.add_reactions([rxn])
        rxn.add_metabolites({mets[m]: c for m, c in stoich.items()})
        rxn.gene_reaction_rule = gpr
    model.objective = "BIOMASS"
    return model


def closed_reactions(model: cobra.Model, gene_id: str) -> set[str]:
    """Reactions whose bounds a KO of ``gene_id`` changes to (0, 0) (already-closed ones excluded)."""
    gene = model.genes.get_by_id(gene_id)
    before = {r.id: r.bounds for r in gene.reactions}
    with model:
        gene.knock_out()
        return {r.id for r in gene.reactions if r.bounds == (0, 0) != before[r.id]}


def verify_toy(processes: int) -> None:
    for solver in ["glpk", "gurobi"]:
        model = toy_model(solver)
        res = single_gene_knockouts(model, medium=TOY_MEDIUM, objective="BIOMASS").set_index("gene")
        got = essential_genes(res.reset_index())
        check(f"toy_wt_growth[{solver}]", math.isclose(res.wt_growth.iloc[0], TOY_WT, abs_tol=1e-6),
              f"wt={res.wt_growth.iloc[0]:.6g} expected {TOY_WT}")
        check(f"toy_essential_set[{solver}]", got == TOY_ESSENTIAL,
              f"got {sorted(got)}; missing {sorted(TOY_ESSENTIAL - got)}; extra {sorted(got - TOY_ESSENTIAL)}")
        bad = {g: res.ratio[g] for g, r in TOY_RATIOS.items() if not math.isclose(res.ratio[g], r, abs_tol=1e-6)}
        check(f"toy_exact_ratios[{solver}]", not bad, f"mismatches {bad}" if bad else f"{len(TOY_RATIOS)} ratios exact")
        bad = {g: res.status[g] for g, s in TOY_STATUS.items() if res.status[g] != s}
        check(f"toy_infeasible_ko[{solver}]", not bad and res.ko_growth["g_mx"] == 0.0,
              f"g_mx status={res.status['g_mx']} growth={res.ko_growth['g_mx']}")
        check(f"toy_model_unchanged[{solver}]",
              model.medium == {"EX_a": 10, "EX_f": 10} and all(not r.bounds == (0, 0) for r in model.reactions)
              and all(g.functional for g in model.genes), "bounds/medium/gene states restored after run")

    model = toy_model()
    bad = {g: (closed_reactions(model, g), exp) for g, exp in TOY_CLOSED.items() if closed_reactions(model, g) != exp}
    check("toy_gpr_closure", not bad and set(TOY_CLOSED) == {g.id for g in model.genes},
          f"mismatches {bad}" if bad else f"{len(TOY_CLOSED)} genes: closed reactions as expected")

    # Medium swap: without f uptake, f must be synthesised -> g_fsyn essential, WT 3, g_fup ratio 1.
    res = single_gene_knockouts(model, medium={"EX_a": 10.0}, objective="BIOMASS")
    exp = TOY_ESSENTIAL | {"g_fsyn"}
    check("toy_medium_swap", essential_genes(res) == exp and math.isclose(res.wt_growth[0], 3.0, abs_tol=1e-6),
          f"wt={res.wt_growth[0]:.6g}; essential {sorted(essential_genes(res))}")
    # Objective swap: maximise w secretion -> needs T_a, R_side, and the maintenance outlet.
    res = single_gene_knockouts(model, medium=TOY_MEDIUM, objective="EX_w")
    exp = {"g_ta", "g_side", "g_mx"}
    check("toy_objective_swap", essential_genes(res) == exp and math.isclose(res.wt_growth[0], 9.0, abs_tol=1e-6),
          f"wt={res.wt_growth[0]:.6g}; essential {sorted(essential_genes(res))}")
    # Blocked objective: no carbon uptake -> no WT growth; ratios undefined (NaN), nothing called essential.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        res = single_gene_knockouts(model, medium={"EX_f": 10.0}, objective="BIOMASS")
    try:
        essential_genes(res)
        nan_raises = False
    except ValueError:
        nan_raises = True
    check("toy_blocked_objective",
          res.ratio.isna().all() and res.ratio_r.isna().all() and (res.status == "wt_no_growth").all()
          and nan_raises and not essential_genes(res, nan="exclude") and caught,
          f"wt={res.wt_growth[0]:.3g}; all ratios NaN, status wt_no_growth, warning raised, "
          "essential_genes raises on NaN (nan='exclude' -> empty)")
    verify_toy_status_policy(model)
    # Serial == parallel (spawn workers, small chunks).
    serial = single_gene_knockouts(model, medium=TOY_MEDIUM, objective="BIOMASS")
    par = single_gene_knockouts(model, medium=TOY_MEDIUM, objective="BIOMASS", processes=min(processes, 3), chunksize=4)
    diff = (serial.ratio - par.ratio).abs().max()
    check("toy_serial_eq_parallel", serial.gene.equals(par.gene) and serial.status.equals(par.status) and diff < 1e-9
          and serial.ratio_r.equals(par.ratio_r),
          f"{len(par)} genes, {min(processes, 3)} processes, max |ratio diff|={diff:.2g}, ratio_r identical")
    # Negative control: unknown gene ids are rejected.
    try:
        single_gene_knockouts(model, genes=["g_ta", "not_a_gene"], medium=TOY_MEDIUM, objective="BIOMASS")
        check("toy_unknown_gene_rejected", False, "no error raised")
    except KeyError as err:
        check("toy_unknown_gene_rejected", True, str(err))


def verify_toy_status_policy(model: cobra.Model) -> None:
    # ratio_r: rounded to 6 decimals, clipped to [0, 1], no negative zero.
    raw = pd.Series([1 + 2e-13, 1 - 1e-13, -1e-13, 0.3333333333, 1.5, math.nan])
    r = rounded_ratio(raw)
    check("toy_ratio_r_rounding",
          r[:5].tolist() == [1.0, 1.0, 0.0, 0.333333, 1.0] and math.isnan(r[5]) and math.copysign(1, r[2]) == 1,
          f"{raw.tolist()} -> {r.tolist()}")
    # Unbounded wild type (unbounded source and drain of a_c, source maximised) must raise, not
    # return NaN ratios.
    with model:
        src = cobra.Reaction("SRC", lower_bound=0, upper_bound=math.inf)
        snk = cobra.Reaction("SNK", lower_bound=0, upper_bound=math.inf)
        model.add_reactions([src, snk])
        src.add_metabolites({model.metabolites.a_c: 1})
        snk.add_metabolites({model.metabolites.a_c: -1})
        try:
            single_gene_knockouts(model, medium=TOY_MEDIUM, objective="SRC")
            check("toy_wt_nonoptimal_raises", False, "no error raised")
        except NonOptimalError as err:
            check("toy_wt_nonoptimal_raises", True, str(err))
    # A KO ending in a non-optimal status other than infeasible (simulated time limit) must raise;
    # with strict=False it is kept as NaN with its status, and essential_genes refuses NaN.
    original = ko_module._solve
    ko_module._solve = lambda m: (math.nan, "time_limit") if not m.genes.g2a.functional else original(m)
    try:
        try:
            single_gene_knockouts(model, medium=TOY_MEDIUM, objective="BIOMASS")
            raised = False
        except NonOptimalError:
            raised = True
        res = single_gene_knockouts(model, medium=TOY_MEDIUM, objective="BIOMASS", strict=False).set_index("gene")
    finally:
        ko_module._solve = original
    try:
        essential_genes(res.reset_index())
        refused = False
    except ValueError:
        refused = True
    check("toy_ko_nonoptimal_raises",
          raised and res.status["g2a"] == "time_limit" and math.isnan(res.ratio_r["g2a"]) and refused
          and (res.status.drop("g2a") != "time_limit").all(),
          "strict run raised; strict=False kept g2a as time_limit/NaN; essential_genes refused NaN")
    # Tasks file must come from the model's own release; untagged models without a medium are refused.
    try:
        configure(model)
        check("toy_untagged_tasks_file_refused", False, "no error raised")
    except ValueError as err:
        check("toy_untagged_tasks_file_refused", True, str(err)[:90])


# ---------------------------------------------------------------- GPR logic
def eval_rule(rule: str, knocked: str) -> bool:
    """Independent GPR evaluator: True if the reaction stays active after knocking out one gene."""
    def ev(node) -> bool:
        if isinstance(node, ast.BoolOp):
            vals = [ev(v) for v in node.values]
            return any(vals) if isinstance(node.op, ast.Or) else all(vals)
        if isinstance(node, ast.Name):
            return node.id != knocked
        raise TypeError(type(node))
    return True if not rule.strip() else ev(ast.parse(rule, mode="eval").body)


def verify_gpr_closure(name: str, model: cobra.Model) -> None:
    mismatches = []
    for gene in model.genes:
        expected = {r.id for r in gene.reactions
                    if not eval_rule(r.gene_reaction_rule, gene.id) and r.bounds != (0, 0)}
        if closed_reactions(model, gene.id) != expected:
            mismatches.append(gene.id)
    n_buffered = sum(1 for g in model.genes if not closed_reactions(model, g.id))
    check(f"gpr_closure_{name}", not mismatches,
          f"{len(model.genes)} genes; mismatches {mismatches[:5]}; {n_buffered} genes close no reaction (fully buffered)")


# ---------------------------------------------------------------- genome-scale controls
def symbol_to_gene(gem) -> dict[str, str]:
    table = read_genes_tsv(gem.genes_tsv)
    return dict(zip(table.geneSymbols, table.index))


def verify_gem(name: str, model: cobra.Model, processes: int, expected_wt: float | None) -> None:
    gem = GEMS[name]
    sym = symbol_to_gene(gem)
    derived = tasks_file_for(model)
    check(f"{name}_tasks_file_derived", gem_of(model) is gem and derived == gem.essential_tasks,
          f"{derived.relative_to(PROJECT_ROOT)}")
    controls = ["CAD", "UMPS"]
    res = single_gene_knockouts(model, genes=[sym[s] for s in controls], tasks_file=gem.essential_tasks)
    wt = res.wt_growth[0]
    if expected_wt is not None:
        check(f"{name}_wt_growth", math.isclose(wt, expected_wt, rel_tol=1e-4), f"wt={wt:.4f} expected {expected_wt}")
    else:
        check(f"{name}_wt_growth", wt > 1, f"wt={wt:.4f} (> 1 required)")
    for s, r in zip(controls, res.ratio):
        check(f"{name}_lethal_{s}", r < 0.01, f"ratio={r:.3g}")
    if name == "human1":  # diagnostic for the known UMPS bypass: make TYMS irreversible as in Human2
        with model:
            model.reactions.MAR04644.lower_bound = 0
            r = single_gene_knockouts(model, genes=[sym["UMPS"]], tasks_file=gem.essential_tasks).ratio[0]
        check("human1_UMPS_lethal_if_TYMS_irreversible", r < 0.01, f"MAR04644 lb=0: UMPS ratio={r:.3g}")

    # FVA consistency: genes whose KO closes exactly one reaction (and no other) -> ratio < 1 iff the
    # reaction must carry flux at optimal growth. Sample 60 such genes (seed 0).
    with model:
        configure(model, tasks_file=gem.essential_tasks)
        single = {g.id: next(iter(closed)) for g in model.genes if len(closed := closed_reactions(model, g.id)) == 1}
        sample = sorted(random.Random(SEED).sample(sorted(single), 60))
        fva = flux_variability_analysis(model, sorted({single[g] for g in sample}), fraction_of_optimum=1.0,
                                        processes=1)
    ko = single_gene_knockouts(model, genes=sample, tasks_file=gem.essential_tasks).set_index("gene")
    tol = 1e-6
    required = {g: not (fva.minimum[single[g]] <= tol and fva.maximum[single[g]] >= -tol) for g in sample}
    bad = [g for g in sample if required[g] != (ko.ratio[g] < 1 - 1e-6)]
    check(f"{name}_fva_consistency", not bad,
          f"{len(sample)} single-reaction KOs; {sum(required.values())} required; mismatches {bad[:5]}")

    # Named non-lethal control: a gene whose only reaction is a blocked/unrelated one, picked from the
    # FVA sample as the first gene whose reaction has an FVA range of exactly [0, 0] at optimal growth.
    zero = [g for g in sample if abs(fva.minimum[single[g]]) < tol and abs(fva.maximum[single[g]]) < tol]
    if zero:
        g = zero[0]
        symbol = read_genes_tsv(gem.genes_tsv).geneSymbols.get(g, "?")
        check(f"{name}_nonlethal_control", math.isclose(ko.ratio[g], 1.0, abs_tol=1e-6),
              f"{symbol} ({g}, only closes {single[g]}, flux fixed at 0 at optimum): ratio={ko.ratio[g]:.6g}")
    else:
        check(f"{name}_nonlethal_control", False, "no zero-flux single-reaction gene in sample")

    # Negative controls: non-model IDs are rejected (TP53 Ensembl, a random non-model DepMap Entrez id,
    # and a metabolic gene symbol, which is not an Ensembl id).
    model_entrez = set(read_genes_tsv(gem.genes_tsv).geneEntrezID.dropna().str.split(";").explode())
    depmap = pd.read_csv(PROJECT_ROOT / "data" / "depmap" / "CRISPRGeneEffect.csv", nrows=0).columns[1:]
    non_model = sorted(c for c in depmap if c.rsplit("(", 1)[-1].rstrip(")") not in model_entrez)
    random_label = random.Random(SEED).choice(non_model)
    random_entrez = random_label.rsplit("(", 1)[-1].rstrip(")")
    for label, bad_id in [("TP53_ensembl", "ENSG00000141510"), (f"random_entrez_{random_label}", random_entrez),
                          ("symbol_CAD", "CAD")]:
        try:
            check_genes(model, [bad_id])
            check(f"{name}_reject_{label}", False, f"{bad_id} accepted")
        except KeyError:
            check(f"{name}_reject_{label}", True, f"{bad_id} rejected")

    # Serial == parallel on a random subset of 48 genes.
    sub = sorted(random.Random(SEED + 1).sample(sorted(g.id for g in model.genes), 48))
    serial = single_gene_knockouts(model, genes=sub, tasks_file=gem.essential_tasks)
    par = single_gene_knockouts(model, genes=sub, tasks_file=gem.essential_tasks, processes=min(processes, 4), chunksize=12)
    diff = (serial.ratio - par.ratio).abs().max()
    check(f"{name}_serial_eq_parallel",
          diff < 1e-6 and (serial.status == par.status).all() and serial.ratio_r.equals(par.ratio_r),
          f"48 genes; max |ratio diff|={diff:.2g}; ratio_r bitwise identical")

    # Genome-scale infeasible KO: force biomass >= 1; a lethal KO (CAD) then makes the LP infeasible,
    # which must be recorded as status infeasible, growth 0, ratio 0 (essential); a neutral KO stays optimal.
    neutral = zero[0] if zero else sample[0]
    with model:
        model.reactions.get_by_id(BIOMASS).lower_bound = 1.0
        res = single_gene_knockouts(model, genes=[sym["CAD"], neutral]).set_index("gene")
    cad = res.loc[sym["CAD"]]
    check(f"{name}_infeasible_ko",
          cad.status == "infeasible" and cad.ko_growth == 0.0 and cad.ratio_r == 0.0
          and essential_genes(res.reset_index()) == {sym["CAD"]} and res.status[neutral] == "optimal"
          and res.ratio_r[neutral] == 1.0,
          f"biomass lb=1: CAD status={cad.status} ratio_r={cad.ratio_r}; {neutral} status={res.status[neutral]} "
          f"ratio_r={res.ratio_r[neutral]}")


def verify_solvers(model: cobra.Model) -> None:
    gem = GEMS["human2"]
    sub = sorted(random.Random(SEED).sample(sorted(g.id for g in model.genes), 200))
    t = time.perf_counter()
    grb = single_gene_knockouts(model, genes=sub, tasks_file=gem.essential_tasks, solver="gurobi")
    t_grb = time.perf_counter() - t
    t = time.perf_counter()
    glpk = single_gene_knockouts(model, genes=sub, tasks_file=gem.essential_tasks, solver="glpk")
    t_glpk = time.perf_counter() - t
    diff = (grb.ratio - glpk.ratio).abs()
    agree = ((grb.ratio < 0.01) == (glpk.ratio < 0.01)) & ((grb.ratio < 0.5) == (glpk.ratio < 0.5))
    check("solver_gurobi_vs_glpk_200", diff.max() < 1e-4 and agree.all(),
          f"WT {grb.wt_growth[0]:.6g} vs {glpk.wt_growth[0]:.6g}; max |ratio diff|={diff.max():.2g}; "
          f"lethal calls agree {agree.sum()}/200; {(grb.ratio < 0.01).sum()} lethal; "
          f"time gurobi {t_grb:.0f}s, glpk {t_glpk:.0f}s (incl. solver switch)")
    if not agree.all():
        print(pd.concat([grb[~agree], glpk[~agree]]).to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--skip-glpk", action="store_true")
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()

    t0 = time.perf_counter()
    verify_toy(args.processes)
    human2 = load_model(GEMS["human2"].sbml)
    verify_gpr_closure("human2", human2)
    verify_gem("human2", human2, args.processes, expected_wt=74.7568)
    if not args.skip_glpk:
        verify_solvers(human2)
    del human2
    human1 = load_model(GEMS["human1"].sbml)
    verify_gpr_closure("human1", human1)
    verify_gem("human1", human1, args.processes, expected_wt=None)

    df = pd.DataFrame(RESULTS, columns=["check", "passed", "detail", "known_model_issue"])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, sep="\t", index=False)
    unexpected = df[~df.passed & (df.known_model_issue == "")]
    print(f"\n{df.passed.sum()}/{len(df)} checks passed ({(~df.passed & (df.known_model_issue != '')).sum()} known"
          f" model issues, {len(unexpected)} unexpected failures) in {time.perf_counter() - t0:.0f}s -> {OUT}")
    sys.exit(0 if unexpected.empty else 1)


if __name__ == "__main__":
    main()
