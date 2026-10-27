"""W4 review, part 2: GPR evaluation and the A4 (a) KO gene-state convention on hand-made toy GPRs (GLPK only).

1. Random nested AND/OR GPRs with on/off/unknown genes: author's presence_score vs the reviewer's own
   evaluator vs cobra GPR.eval; checks s_r < 0 <=> GPR False with the off genes removed.
2. Hand-made toys through the author's closure_map / ko_sweep: an off isozyme does not buffer (presence
   families); nested AND/OR; complexes; a reaction kept despite a False pruned GPR keeps its full GPR;
   E-Flux KO zeroes the gene's capacity (and closes a reaction whose GPR becomes False).
"""

import ast
import json
import math
import sys
from pathlib import Path

import cobra

cobra.Configuration().solver = "glpk"
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import own  # noqa: E402

from nextgem import context as cx  # noqa: E402

OUT = Path(__file__).parent
rng = np.random.default_rng(11)
res = {}


def rand_gpr(depth, genes):
    if depth == 0 or rng.random() < 0.3:
        return str(rng.choice(genes))
    op = " and " if rng.random() < 0.5 else " or "
    k = int(rng.integers(2, 4))
    return "(" + op.join(rand_gpr(depth - 1, genes) for _ in range(k)) + ")"


genes = [f"g{i}" for i in range(8)]
n, bad_auth_vs_own, bad_sign = 0, 0, 0
for _ in range(5000):
    rule = rand_gpr(3, genes)
    gpr = cobra.core.GPR.from_string(rule)
    s = {}
    for g in genes:
        u = rng.random()
        s[g] = math.nan if u < 0.25 else float(rng.normal())
    sd = {g: v for g, v in s.items() if not math.isnan(v)}
    node = cx._compile(gpr.body)
    a_val = cx.presence_score(node, sd)
    o_val = own.presence(gpr.body, s)[0]
    if not ((math.isnan(a_val) and math.isnan(o_val)) or a_val == o_val):
        bad_auth_vs_own += 1
    off = {g for g, v in sd.items() if v < 0}
    gpr_false = not gpr.eval(knockouts=off)
    if (not math.isnan(a_val) and a_val < 0) != gpr_false:
        bad_sign += 1
    n += 1
res["random_gprs"] = {"n": n, "author_vs_own_mismatch": bad_auth_vs_own, "sign_rule_violations": bad_sign}


# ------------------------------------------------------------------ toy network for closure_map / ko_sweep
def toy(spec):
    m = cobra.Model("toy")
    for rid, eq, gpr, lb, ub in spec:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.build_reaction_from_string(eq, verbose=False)
        r.gene_reaction_rule = gpr
        r.bounds = (lb, ub)
    m.objective = cx.BIOMASS_REACTION
    return cx.parent_from_model(m)


BIO = cx.BIOMASS_REACTION
SPEC = [
    ("EX_a", "--> a", "", 0, 10),
    ("Riso", "a --> b", "gI1 or gI2", 0, 1000),                 # isozymes; gI2 off
    ("Rnest", "b --> c", "(gA and gB) or (gC and gD)", 0, 1000),  # nested; gC off
    ("Rcplx", "c --> d", "gX and gY", 0, 1000),                 # complex; gY off -> F0 closes, GIMME/repair may keep
    ("Rbyp", "c --> d", "gZ", 0, 1000),                         # bypass of the complex, gZ on
    ("Runk", "d --> e", "gU and (gV or gW)", -1000, 1000),      # gV off, gW unknown
    (BIO, "e --> ", "", 0, 1000),
]
p = toy(SPEC)
s = pd.Series({"gI1": 2.0, "gI2": -1.0, "gA": 1.0, "gB": 1.5, "gC": -0.5, "gD": 3.0, "gX": 1.0, "gY": -2.0,
               "gZ": 1.0, "gU": 1.0, "gV": -1.0}).reindex(p.genes)  # gW: no evidence (NaN)
off = set(s.index[(s < 0).to_numpy()])

# presence family with Rcplx kept open despite a False pruned GPR (as GIMME / repair would)
lm = cx.LineModel(line="toy", family="F2", closed=set(), off=off)
close, rebound = cx.closure_map(p, lm)
close = {g: sorted(v) for g, v in close.items()}
exp = {"gI1": ["Riso"], "gA": ["Rnest"], "gB": ["Rnest"], "gX": ["Rcplx"], "gY": ["Rcplx"], "gZ": ["Rbyp"],
       "gU": ["Runk"], "gW": ["Runk"]}
res["closure_presence"] = {"got": close, "expected": exp, "pass": close == exp}
ko = cx.ko_sweep(p, lm).set_index("gene").ratio_r.to_dict()
wt, own_ko = own.ko_presence(p.model, set(), off)
res["ko_presence"] = {"author": ko, "own_bruteforce": {g: round(v, 6) for g, v in own_ko.items()},
                      "equal": all(abs(ko[g] - own_ko[g]) < 1e-9 for g in ko),
                      "off_isozyme_does_not_buffer": ko["gI1"] == 0.0,
                      "generic_parent_buffers": None}
# generic parent (no off genes): gI1 KO is buffered by gI2
lm0 = cx.LineModel(line="toy", family="F0", closed=set(), off=set())
res["ko_presence"]["generic_parent_buffers"] = cx.ko_sweep(p, lm0).set_index("gene").ratio_r["gI1"] == 1.0
# F0 on the same scores: Rcplx closed (GPR False under off); bypass carries flux; gZ KO lethal, gX KO no effect
lmF0 = cx.extract(p, "F0", s)
res["f0"] = {"closed": sorted(lmF0.closed), "expected_closed": ["Rcplx"]}

# E-Flux: lin = 2**s (on) / 0 (off); KO zeroes the gene's linear value; GPR-False reactions are closed
lmE = cx.extract(p, "F4", s)
lin = lmE.lin
bnd = {k: tuple(round(x, 6) for x in v) for k, v in lmE.eflux.items()}
close_e, rebound_e = cx.closure_map(p, lmE)
res["eflux"] = {"lin": lin, "cap_max": lmE.cap_max, "bounds": bnd,
                "closures": {g: sorted(v) for g, v in close_e.items()},
                "rebound": {g: {r: tuple(round(x, 6) for x in b) for r, b in v.items()} for g, v in rebound_e.items()}}
# expected capacities: Riso = 4 + 0 = 4; Rnest = min(2, 2.83) + min(0, 8) = 2; Rcplx = min(2, 0) = 0; Rbyp = 2;
# Runk: OR(gV=0, gW unknown) -> unknown -> AND ignores -> min(2) = 2.  cap_max = 4.
exp_b = {"Riso": (0, 1000.0), "Rnest": (0, 500.0), "Rcplx": (0, 10.0), "Rbyp": (0, 500.0), "Runk": (-500.0, 500.0)}
res["eflux"]["expected_bounds"] = exp_b
res["eflux"]["bounds_pass"] = all(abs(bnd[k][0] - v[0]) < 1e-6 and abs(bnd[k][1] - v[1]) < 1e-6 for k, v in exp_b.items())
# KO gI1: Riso cap 0 but GPR still True (gI2 present though off) -> floor 10, not closed
res["eflux"]["ko_gI1_rebound_floor"] = rebound_e.get("gI1", {}).get("Riso") == (0, 10.0)
# KO gZ: Rbyp GPR False -> closed
res["eflux"]["ko_gZ_closes"] = "Rbyp" in close_e.get("gZ", set())
koE = cx.ko_sweep(p, lmE).set_index("gene").ratio_r.to_dict()
b_own, cm = own.eflux_bounds(p.model, own.eflux_lin(s.to_dict()))
_, own_koE = own.ko_eflux(p.model, b_own, own.eflux_lin(s.to_dict()), cm)
res["eflux"]["ko_author"] = koE
res["eflux"]["ko_equal_own"] = all(abs(koE[g] - own_koE[g]) < 1e-6 for g in koE)

# genomic loss in a presence family: lost gene non-functional everywhere; reopened reaction keeps full GPR
lmL = cx.LineModel(line="toy", family="F0", closed={"Rcplx"}, off=off, lost={"gI2"})
lmL.lost_closed = cx.loss_closure(p, lmL.lost)
lmL.reopened = {"Rcplx"}
close_l, _ = cx.closure_map(p, lmL)
res["genomic"] = {"closures": {g: sorted(v) for g, v in close_l.items()},
                  "gY_ko_closes_reopened_complex": "Rcplx" in close_l.get("gY", set())}
(OUT / "r2_gpr_toys.json").write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
