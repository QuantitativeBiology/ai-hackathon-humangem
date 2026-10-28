"""Assemble final/index.html (standalone, with doctype) and final/index.artifact.html (no html/head/body wrapper)
from scripts/template.html, the figure SVGs (md/figures) and the worked-example JSON (data/worked_example.json).

Usage: python build_html.py   (paths are relative to this file's parent folder, final/)
"""
import html
import json
import re
from pathlib import Path

FINAL = Path(__file__).resolve().parents[1]
FIG = FINAL / "md" / "figures"
DATA = FINAL / "data"


# ------------------------------------------------------------------ figures
def inline_svg(name: str) -> str:
    s = (FIG / f"{name}.svg").read_text()
    s = re.sub(r"<!--FX-STYLE-->.*?<!--/FX-STYLE-->", "", s, flags=re.S)
    s = re.sub(r'(<svg[^>]*?) width="[^"]+" height="[^"]+"', r"\1", s, count=1)
    s = re.sub(r"var\(--(fx-[a-z0-9]+), #[0-9a-f]{6}\)", r"var(--\1)", s)
    return s


# ------------------------------------------------------------------ pipeline diagram
def stage(x, y, w, h, target, num, lines, sub, mono, kind="", label=""):
    cls = "pd-stage" + (f" pd-{kind}" if kind else "")
    out = [f'<g class="{cls}" data-target="{target}" role="link" tabindex="0" aria-label="{html.escape(label)}">',
           f'<rect class="pd-box" x="{x}" y="{y}" width="{w}" height="{h}" rx="5"/>']
    ty = y + 24
    if num:
        out.append(f'<circle cx="{x + 16}" cy="{y + 18}" r="10" style="fill:none;stroke:var(--accent);stroke-width:1.2"/>'
                   f'<text class="pd-n" x="{x + 16}" y="{y + 22}" text-anchor="middle">{num}</text>')
        tx = x + 32
    else:
        tx = x + 12
    for i, ln in enumerate(lines):
        out.append(f'<text class="pd-t" x="{tx}" y="{ty + i * 17}">{html.escape(ln)}</text>')
    yy = ty + len(lines) * 17 + 4
    for i, ln in enumerate(sub):
        out.append(f'<text class="pd-s" x="{x + 12}" y="{yy + i * 15}">{html.escape(ln)}</text>')
    if mono:
        out.append(f'<text class="pd-m" x="{x + 12}" y="{y + h - 10}">{html.escape(mono)}</text>')
    out.append("</g>")
    return "".join(out)


def inp(x, y, w, h, title, sub, mono):
    return (f'<g class="pd-input"><rect class="pd-box" x="{x}" y="{y}" width="{w}" height="{h}" rx="5"/>'
            f'<text class="pd-s" x="{x + 10}" y="{y + 17}" style="text-transform:uppercase;letter-spacing:.06em;font-size:10px">input</text>'
            f'<text class="pd-t" x="{x + 10}" y="{y + 36}">{html.escape(title)}</text>'
            f'<text class="pd-s" x="{x + 10}" y="{y + 52}">{html.escape(sub)}</text>'
            f'<text class="pd-m" x="{x + 10}" y="{y + h - 9}">{html.escape(mono)}</text></g>')


def arrow(d, m="pd-ah"):
    return f'<path class="pd-arrow" d="{d}" marker-end="url(#{m})"/>'


DEFS = ('<defs><marker id="{mid}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
        '<path class="pd-head" d="M0,0 L10,5 L0,10 z"/></marker></defs>')


def defs(mid):
    return DEFS.replace("{mid}", mid)

LABELS = {
    "step-1": "Step 1: generic Human1 network",
    "step-2": "Step 2: remove the cofactor pool from biomass",
    "step-3": "Step 3: medium and serum",
    "step-4": "Step 4: knock out each gene",
    "step-5": "Step 5: read as a tied ranking against baselines",
    "worked-example": "Output: worked example",
}


def pipeline_svg() -> str:
    W = ['<svg class="pd-wide" viewBox="0 0 1000 262" role="img" aria-labelledby="pd-t1"><title id="pd-t1">P1 pipeline: '
         'inputs, four stages, output and how it is read</title>', defs("pd-ah")]
    W.append(inp(8, 26, 150, 82, "Human1 v1.19.0", "SBML, 2,887 genes", "Human-GEM"))
    W.append(inp(8, 132, 150, 82, "Line's medium", "DepMap annotation", '"RPMI + 10% FBS"'))
    sw, sy, sh = 150, 40, 140
    xs = [196, 362, 528, 694]
    W.append(stage(xs[0], sy, sw, sh, "step-1", "1", ["Generic", "network"], ["gene rules as", "released; no omics"], "Human1 as is", label=LABELS["step-1"]))
    W.append(stage(xs[1], sy, sw, sh, "step-2", "2", ["Remove", "cofactor pool"], ["from the biomass", "reaction (all lines)"], "MAR13082 − MAM10012c", "fix", LABELS["step-2"]))
    W.append(stage(xs[2], sy, sw, sh, "step-3", "3", ["Medium", "+ serum"], ["uptake limits per", "line; RPMI if unknown"], "not line-specific", "fix", LABELS["step-3"]))
    W.append(stage(xs[3], sy, sw, sh, "step-4", "4", ["Knock out", "each gene"], ["FBA: e = 1 − KO /", "normal growth"], "2,866 genes", label=LABELS["step-4"]))
    W.append(stage(862, sy, 130, sh, "worked-example", "", ["Score e", "per gene"], ["a ranking, mostly", "tied at 0"], "NCI-H460: 106 lethal", "out", LABELS["worked-example"]))
    for a, b in ((xs[0] + sw, xs[1]), (xs[1] + sw, xs[2]), (xs[2] + sw, xs[3]), (xs[3] + sw, 862)):
        W.append(arrow(f"M{a + 2},{sy + sh / 2} L{b - 3},{sy + sh / 2}"))
    W.append(arrow(f"M160,67 C178,67 176,{sy + 40} {xs[0] - 3},{sy + 40}"))
    W.append(arrow(f"M160,173 C300,244 {xs[2] + 30},244 {xs[2] + 60},{sy + sh + 3}"))
    W.append(stage(694, 202, 298, 52, "step-5", "5", ["Read as a tied ranking (tAP)"], [], "vs gene prior and expression baselines", label=LABELS["step-5"]))
    W.append(arrow(f"M927,{sy + sh + 2} L927,199"))
    W.append("</svg>")

    T = ['<svg class="pd-tall" viewBox="0 0 360 760" role="img" aria-labelledby="pd-t2"><title id="pd-t2">P1 pipeline, '
         'vertical layout</title>', defs("pd-ah2")]
    T.append(inp(8, 8, 168, 82, "Human1 v1.19.0", "SBML, 2,887 genes", "Human-GEM"))
    T.append(inp(184, 8, 168, 82, "Line's medium", "DepMap annotation", '"RPMI + 10% FBS"'))
    x, w, h = 30, 300, 92
    ys = [118, 232, 346, 460]
    T.append(stage(x, ys[0], w, h, "step-1", "1", ["Generic network"], ["gene rules as released; no omics"], "Human1 as is", label=LABELS["step-1"]))
    T.append(stage(x, ys[1], w, h, "step-2", "2", ["Remove cofactor pool"], ["from the biomass reaction (all lines)"], "MAR13082 − MAM10012c", "fix", LABELS["step-2"]))
    T.append(stage(x, ys[2], w, h, "step-3", "3", ["Medium + serum"], ["uptake limits per line; RPMI if unknown"], "not line-specific", "fix", LABELS["step-3"]))
    T.append(stage(x, ys[3], w, h, "step-4", "4", ["Knock out each gene"], ["FBA: e = 1 − KO / normal growth"], "2,866 genes", label=LABELS["step-4"]))
    T.append(stage(x, 574, w, 82, "worked-example", "", ["Score e per gene"], ["a ranking, mostly tied at 0"], "106 lethal · 64 partial · 2,696 tied", "out", LABELS["worked-example"]))
    T.append(stage(x, 678, w, 74, "step-5", "5", ["Read as a tied ranking"], [], "tAP vs prior and expression", label=LABELS["step-5"]))
    T.append(arrow("M92,92 L92,115", "pd-ah2"))
    for a in ys:
        T.append(arrow(f"M180,{a + h + 2} L180,{a + h + 19}", "pd-ah2"))
    T.append(arrow("M180,658 L180,675", "pd-ah2"))
    T.append(arrow("M352,60 C358,250 352,380 333,390", "pd-ah2"))
    T.append("</svg>")
    return "".join(W) + "".join(T)


# ------------------------------------------------------------------ turning points
TPS = [
    dict(id="tp-1", n="1", when="Day 1 · 14:46–17:05 · commits b105f88 → 61b0beb", cls="",
         title="Measure honestly before measuring anything",
         gist="Ties, common essentials and missing fair baselines would have made any 'better model' claim unreadable; the design fixed the metric, the split and the baselines before any scoring.",
         rows=[("Situation", "The PI drafted a design (DESIGN.md v0) and sent it to seven reviewer agents. What metric and split make a 'better metabolic model' claim believable?"),
               ("Evidence", "Generic knockout scores are nearly binary (about 2,650 of 2,870 genes at ratio 1); on data of that shape, random tie order alone moved standard AP by 0.025–0.053 per line. Common essentials are 86.5% of positives. Published GEM-vs-CRISPR comparisons report best-of-settings MCC without tissue hold-out. <span class=\"src\">notes/design_review_r1/statistics.md; LITERATURE.md §0</span>"),
               ("Options", "Plain AUPRC or MCC with random line splits; or tie-aware AP, tie-break baselines and whole-lineage hold-out."),
               ("Outcome", "Tie-aware AP as primary metric; expression tie-break baselines; a second endpoint without common essentials; a fixed order of confirmatory tests; ≥ 7 test lineages. Two review rounds produced 51 decisions; an unexposed verifier added 17 fixes before the freeze (61b0beb). A hook blocks raw label reads (user: \"yes, add the hook\"); a retired seed was replaced by the user's seed 38, revealed after the freeze."),
               ("Human input", "\"5–6 days compute is too much\"; ftINIT excluded (\"poor performance from our own analysis\"); compare only with the generic model, add MCC/F1. Budget ≈ 2,950 → ≈ 1,020 CPU-h. <span class=\"src\">RESEARCH_TREE.md D2; a8ba7b8</span>"),
               ("Produced", '<a href="#step-5">Step 5</a> (tied rankings vs baselines); the rounding in <a href="#step-4">step 4</a> came from the first baseline review (float noise moved AP by 0.018).')]),
    dict(id="tp-2", n="2", when="Day 2 · 00:50 · commits 67f5f4f, 09b69a2", cls="",
         title="The baseline: Human1 beats Human2, and a gene prior is three times the model",
         gist="The newer Human2 lost to Human1 everywhere, and a label-only gene-frequency prior scored 0.90 against the model's 0.30.",
         rows=[("Situation", "First contact with labels (training only). Which generic model is the reference, and how much room is there?"),
               ("Options", "Human2 by default (its paper reports improvement), or Human1."),
               ("Evidence", "Human1 0.300 vs Human2 0.292 on training, ahead in 12/12 lineages; validation Human2 − Human1 = −0.007. Gene prior 0.900; expression alone 0.326; reference with an expression tie-break 0.424; ~45% of essential pairs invisible to any knockout. <span class=\"src\">results/w1/README.md; train_grid.tsv; val_look_W1_M2_baseline_gstar.json; train_addressable.tsv</span>"),
               ("Outcome", "Human1 is P1's base; every result is reported against the gene prior and the expression tie-break; the ML search was scoped to gene-prior hybrids."),
               ("Produced", '<a href="#step-1">Step 1</a> (Human1) and the baselines of <a href="#step-5">step 5</a>.')]),
    dict(id="tp-3", n="3", when="Day 2 · 19:46 · commit 75c0805", cls="neg",
         title="Every expression-tailored model loses",
         gist="Five context-specific families and CORDA, 1,380 builds and 108.7 CPU-h: all lowered accuracy, and less pruning was always better.",
         rows=[("Situation", "The brief centres on per-line models from omics. Stage A tested five families on 30 training lines × 3 thresholds × RNA/protein × genomic constraint, plus CORDA on 10 lines (kept at the user's request)."),
               ("Evidence", "Change vs reference: thresholding −0.008, GIMME −0.007, FASTCORE −0.010, E-Flux −0.003, label-calibrated −0.011, CORDA −0.013. Method explained more knockout variance than tissue (R² 0.199 vs 0.085). <span class=\"src\">results/w4/stage_a_analysis/train_selection.json, h4_tests.json; independent recompute</span>"),
               ("Outcome", "On top of the biomass fix, GIMME and E-Flux added +0.003 and +0.004, below the +0.005 bar; group models were never run. Post hoc, E-Flux on full P1 gave −0.0007. <span class=\"src\">results/m5/train/train_decision.json; results/f4/README.md</span>"),
               ("Produced", '<a href="#step-1">Step 1</a>: no extraction.')]),
    dict(id="tp-4", n="4", when="Day 2 · 20:42 · commit b280526", cls="neg",
         title="GPR curation backfires",
         gist="Evidence-based gene-rule edits lost 0.029: a paralog rule merged RNA polymerase subunits and made 15 essential genes 'viable'. A side finding pointed at the biomass.",
         rows=[("Situation", "Wrong isozyme rules are a known error source (the Human2 paper credits its gains to fixing them). Do pre-registered paralog and complex edits help?"),
               ("Evidence", "−0.029 on both models, 0/12 folds, worse than random edits. Paralogs of 21–25% identity were merged inside the lumped RNA-synthesis reaction. A ≥ 30% identity rule still lost 0.002. <span class=\"src\">results/w5/README.md; train_summary.json</span>"),
               ("Outcome", "GPRs kept as released. A one-CPU-hour scan for lumped reactions would have predicted the failure. The reviewer also found that CAD stays lethal under uridine only because dihydroorotate sits in a biomass pool: a biomass clue."),
               ("Produced", '<a href="#step-1">Step 1</a> (rules as released) and a clue behind <a href="#tp-5">TP5</a>.')]),
    dict(id="tp-5", n="5", when="Day 2 19:06 – Day 3 00:47 · commits efc4138, 5537c83", cls="fixes",
         title="The cofactor pool: the one structural fix that works",
         gist="Of 134 ways to change the model's objective, only removing the biomass cofactor pool passed; it was the first fix to beat the reference on unseen tissues (+0.032).",
         rows=[("Situation", "134 pre-registered modes per model: metabolic tasks, alternative objectives (MOMA, ATP maintenance), demands, biomass additions and leave-outs. Are false positives caused by the objective?"),
               ("Where the idea came from", "Leave-one-out of biomass reactants was a <em>diagnostic</em> in the first design. It became a candidate <em>fix</em> through label-free clues: a design reviewer saw expression-pruned models die of 'biomass-artefact lethals: vitamin/xenobiotic CYPs and carnitine synthesis'; the medium pilot found serum cholesterol could not rescue sterol genes because of a vitamin-D catabolite demand; the GPR reviewer found the dihydroorotate entry."),
               ("Evidence", "Pool removal +0.034, 12/12 folds (Human2 +0.037); task modes −0.015; MOMA −0.002. <span class=\"src\">results/w3/README.md</span>"),
               ("The control problem", "The agent's own control (drop the calls with the lowest gene prior) disqualified everything; after seeing that, it switched to the literal pre-registered control, which did not apply to leave-outs. The reviewer ruled the literal reading correct, showed the agent's control was near-oracle, and added a label-free control: random removal loses 0.028, so the targeted removal wins by +0.062. Every P1 claim carries this flag. <span class=\"src\">notes/reviews/M4_W3_review.md; amendment A11</span>"),
               ("Outcome", "Over 134 modes, best-of-many inflation (0.049) exceeded every training effect, so only held-out data could decide: validation +0.032 [0.031, 0.033]. A weekly usage limit had delayed this result by ~6.6 h. <span class=\"src\">results/w3/val/val_look_W3_M4_shortlist_decision.json</span>"),
               ("Produced", '<a href="#step-2">Step 2</a>, and growth-only calls in <a href="#step-4">step 4</a>.')]),
    dict(id="tp-6", n="6", when="Day 3 · 01:46 · commits 80c992d, 1d1c51c", cls="fixes",
         title="The medium helps, but any realistic medium does",
         gist="Documented medium plus serum gained +0.016, but a shuffled medium kept 97% of that and a constant RPMI did slightly better.",
         rows=[("Situation", "Medium is a known confounder of metabolic dependency, and DepMap annotates each screen's medium. Is the gain line-specific?"),
               ("Options", "Ham's; documented base medium; plus serum at several doses; calibrated rates (infeasible); HPLM, Plasmax. Controls: all-open and shuffled media."),
               ("Evidence", "+0.016 with serum, +0.009 without; HPLM +0.012, Plasmax +0.010. Documented − shuffled = +0.0005; constant RPMI + FBS +0.0015 above documented (post hoc); all-open +0.020 overall but −0.017 without common essentials. <span class=\"src\">results/w2/README.md; train/shuffled_control.tsv; step_scores.tsv</span>"),
               ("Outcome", "Kept as 'a realistic medium instead of Ham's', flagged not line-specific; its validation came from the factorial. W2's scoring code was committed only after its first label read."),
               ("Produced", '<a href="#step-3">Step 3</a>.')]),
    dict(id="tp-7", n="7", when="Day 3 · 01:47–07:02 · commits 36b0529, 93d6544", cls="fixes",
         title="The combination wins validation and holds on test",
         gist="Biomass + medium was the only combination within one SE of the best on validation (+0.052), and gave +0.050 on the 7 sealed test tissues, 7/7 positive.",
         rows=[("Situation", "Two fixes passed training, one passed validation. How to combine, and does it hold on test?"),
               ("Options", "Three pre-declared cells: medium only, biomass only, both. A reviewer caught a candidate order that broke the 'simplest first' rule; it was corrected before scoring."),
               ("Evidence", "Training +0.016 / +0.034 / +0.053 (slightly super-additive). Validation +0.017 / +0.032 / +0.052. Test +0.0497 [0.0482, 0.0513], 7/7 lineages; +0.037 without common essentials; line specificity failed (−0.011); P1 loses to the expression tie-break (−0.071). <span class=\"src\">results/m6/train/summary.json; results/m6/val/val_look_M6_factorial_decision.json; results/m6/test/test_look_T1_gates.json</span>"),
               ("Outcome", "P1 = Human1 + pool removed + documented medium with serum; headline confirmed."),
               ("Produced", 'The combination of <a href="#step-2">steps 2</a> and <a href="#step-3">3</a>; the <a href="#results">results</a>.')]),
    dict(id="tp-8", n="8", when="Day 3 · 03:13–05:31 · commits 3153ef0, ea8823f", cls="neg",
         title="Negatives that bound what the model adds",
         gist="Simulated fluxes added nothing beyond omics, and model outputs added nothing to gene-prior machine learning.",
         rows=[("Situation", "Do flux states carry information beyond omics? Do model outputs help ML that already has gene priors?"),
               ("Evidence", "Flux added nothing for subtypes, drug response or metabolomics (all Holm p ≥ 0.38), because generic-model flux depends only on the medium (77/106 drug units had a change of exactly 0). Model features in ML −0.003, stacking +0.000; the one promoted ML idea has no model output, and its lead shrank from +0.014 to +0.0007 on validation. <span class=\"src\">results/w6/; autoresearch/; notes/reviews/M7_*</span>"),
               ("Outcome", "The best ML model uses no GEM. P1 ships as 48 models, one per medium, not one per line."),
               ("Produced", 'The <a href="#limitations">limitations and when-not-to-use</a> list.')]),
    dict(id="epilogue", n="E", when="Day 3 · 09:31–11:36 · post hoc · commits 0da7565, b62bd80, 101aed1", cls="",
         title="Follow-ups F1–F4 leave P1 unchanged",
         gist="A stricter GPR rule failed, a safer biomass variant costs 0.021, the worst errors lie outside the model, and line-specific constraints add nothing.",
         rows=[("Situation", "After the report, the user asked to \"revisit the GPR curation with a stricter paralog rule\", to try an artefact-only biomass variant, to trace the worst errors and to add line-specific constraints. Rules fixed before scoring; training and validation only."),
               ("Evidence", "F1: −0.0006, 0/12 (RRM2 merged with RRM2B). F2: +0.031 vs reference on validation vs +0.052 for P1 (cost ≈ 0.021); decomposition: pool genes 62% of the gain, mevalonate artefact −0.0014. F3: the 20 worst errors are housekeeping genes outside Human1's scope. F4: genomic loss −0.0001; E-Flux −0.0007. <span class=\"src\">results/f1–f4/README.md; results/f2/decomposition.tsv; results/m10/val/</span>"),
               ("Outcome", "P1 stands. The follow-ups explain each step: rules as released, whole pool removed with a small known cost, no line-specific constraints, a recall ceiling set by Human1's scope."),
               ("Produced", '<a href="#follow-ups">Follow-ups</a>, <a href="#step-1">step 1</a>, <a href="#step-2">step 2</a>, <a href="#limitations">limitations</a>.')]),
]


def tp_html() -> str:
    out = []
    for t in TPS:
        rows = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in t["rows"])
        out.append(f'<li class="tp {t["cls"]}" id="{t["id"]}"><span class="tp-mark" aria-hidden="true">{t["n"]}</span>'
                   f'<details><summary><span class="when">{t["when"]}</span><span class="ttl">{html.escape(t["title"])}</span>'
                   f'<span class="gist">{html.escape(t["gist"])}</span><span class="more">Open turning point</span></summary>'
                   f'<div class="tp-body"><dl class="wgw">{rows}</dl></div></details></li>')
    return "\n".join(out)


# ------------------------------------------------------------------ worked example data
def we_data() -> str:
    d = json.loads((DATA / "worked_example.json").read_text())
    st, tap, tr = d["steps"], d["tAP"], d["transitions"]
    panel = [dict(symbol=r["symbol"], pathway=r["pathway"], e0=r["e_step0"], e1=r["e_step1"], em=r["e_medium_only"],
                  e2=r["e_step2"], frac=r.get("test_frac_lines_essential"), rank=r["P1_rank_block"]) for r in d["panel"]]
    cols = [dict(key="e0", label="Step 0", prev=None), dict(key="e1", label="Step 1", prev=0),
            dict(key="em", label="Medium only", prev=0), dict(key="e2", label="Step 2 = P1", prev=1)]

    def mets(k, ref=None):
        s = st[k]
        m = [dict(k="normal growth (a.u.)", v=f"{s['wt_growth_au']:.3f}"),
             dict(k="genes called lethal", v=str(s["n_lethal"]), d=(f"{s['n_lethal'] - st[ref]['n_lethal']:+d} vs previous" if ref else "")),
             dict(k="partial (0 < e < 0.99)", v=str(s["n_partial"])),
             dict(k="distinct scores", v=str(s["n_distinct_scores"])),
             dict(k="NCI-H460 tAP", v=f"{tap[k]['tAP_all_genes']:.3f}",
                  d=(f"{tap[k]['tAP_all_genes'] - tap[ref]['tAP_all_genes']:+.3f} vs previous" if ref else f"random ranking {tap['random_ranking_expectation']:.3f}"))]
        return m

    def moved(key):
        return dict(down=tr[key]["uncalled (lethal -> not)"], up=tr[key]["newly called (not -> lethal)"])

    p1 = {r["symbol"]: r for r in panel}
    top = [s for s in p1 if p1[s]["e2"] >= 0.99]
    mid = [s for s in p1 if 0 < p1[s]["e2"] < 0.99]
    low = [s for s in p1 if p1[s]["e2"] == 0]
    stepper = [
        dict(short="0 · Reference", col=0, title="Step 0. The reference: generic Human1 on Ham's medium",
             text="Human1 as released, with the 44 inputs of Ham's medium, which include hypoxanthine and thymidine. The biomass reaction demands the cofactor pool, so every gene that makes one of its 38 species is lethal. Re-run: 197 linear programs, 42 CPU-s; scores equal the committed table exactly.",
             metrics=mets("step0")),
        dict(short="1 · Remove pool", col=1, title="Step 1. Remove the cofactor pool from the biomass reaction",
             text="One reactant (<code>MAM10012c</code>, coefficient 0.0012) is deleted from <code>MAR13082</code>. Normal growth is unchanged. 35 genes leave the lethal set and none join: retinoid and vitamin-D catabolism (CYP26A1, CYP24A1), carnitine synthesis (BBOX1, TMLHE), and also real cofactors (heme: HMBS, UROD; CoQ: COQ2; B6: PDXK). Across 380 test lines these 35 genes are essential in 14% of lines on average.",
             metrics=mets("step1", "step0"), moved=moved("step0->step1 (pool removed)")),
        dict(short="2 · Medium + serum", col=3, title="Step 2. Add the line's medium and serum: the P1 model",
             text="Uptake switches from Ham's to RPMI-1640 + 10% FBS: hypoxanthine, thymidine, alanine and cysteine disappear; ions, glutathione, cystine, carnitine, serum lipids, cholesterol and ethanolamine appear (58 uptakes, 12 capped at the model's requirement). De novo nucleotide genes become lethal; serum cholesterol rescues sterol genes, and with the pool gone it also rescues the mevalonate core (HMGCR), which is the known cost. Side cell, medium without the biomass fix: 142 lethal, tAP 0.291.",
             metrics=mets("step2", "step1"), moved=moved("step1->step2 (medium + serum added)")),
        dict(short="3 · Rank and score", col=3, nochg=True, title="Step 3. Read the scores as a ranking with ties",
             text="P1's scores for NCI-H460 form three tiers: 106 genes tied at the top, 64 partial genes in 15 levels, and 2,696 genes tied at zero. Tie-aware AP averages precision over random orders within each tied block. The line's tAP comes from the committed per-line table; across all 686 training lines the median gain over the reference is +0.053 (IQR 0.044–0.063), positive in all 686.",
             tapbars=[dict(k="random ranking", v=tap["random_ranking_expectation"], cls="rnd"),
                      dict(k="step 0 reference", v=tap["step0"]["tAP_all_genes"]),
                      dict(k="medium only", v=tap["medium_only"]["tAP_all_genes"]),
                      dict(k="step 1 − pool", v=tap["step1"]["tAP_all_genes"]),
                      dict(k="step 2 = P1", v=tap["step2"]["tAP_all_genes"], cls="p1")],
             blocks=[dict(r="1–106", e="1 (lethal)", n=106, g=", ".join(sorted(top))),
                     dict(r="107–170", e="0 < e < 1 (15 levels)", n=64, g=", ".join(sorted(mid)) + " (ranks 150–160)"),
                     dict(r="171–2,866", e="0 (no effect)", n=2696, g=", ".join(sorted(low)))]),
    ]
    order = ["step0", "step1", "medium_only", "step2"]
    summary = [dict(k="Uptake", v=["Ham's (44)", "Ham's (44)", "RPMI + serum (58)", "RPMI + serum (58)"]),
               dict(k="Biomass", v=["with pool", "pool removed", "with pool", "pool removed"]),
               dict(k="Normal growth (a.u.)", v=[f"{st[k]['wt_growth_au']:.3f}" for k in order]),
               dict(k="Genes called lethal", v=[str(st[k]["n_lethal"]) for k in order]),
               dict(k="Partial", v=[str(st[k]["n_partial"]) for k in order]),
               dict(k="Tied at e = 0", v=[f"{st[k]['n_no_effect']:,}" for k in order]),
               dict(k="NCI-H460 tAP", v=[f"{tap[k]['tAP_all_genes']:.3f}" for k in order]),
               dict(k="Max difference vs committed scores", v=[f"{st[k]['max_abs_dev_vs_committed_sweep_e']:.0f}" for k in order])]
    return json.dumps(dict(panel=panel, cols=cols, stepper=stepper, summary=summary)).replace("</", "<\\/")


def main():
    t = (FINAL / "scripts" / "template.html").read_text()
    t = t.replace("{{PIPELINE_SVG}}", pipeline_svg())
    for i, n in enumerate(["fig1_tap_tiers", "fig2_lineage_gain", "fig3_factorial", "fig4_decomposition", "fig5_timeline"], 1):
        t = t.replace("{{FIG%d}}" % i, inline_svg(n))
    t = t.replace("{{TURNING_POINTS}}", tp_html())
    t = t.replace("{{WE_DATA}}", we_data())
    assert "{{" not in t, re.findall(r"\{\{[A-Z0-9_]+\}\}", t)
    (FINAL / "index.artifact.html").write_text(t)
    head_end = t.index("</style>") + len("</style>")
    full = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            + t[:head_end] + "\n</head>\n<body>\n" + t[head_end:] + "\n</body>\n</html>\n")
    (FINAL / "index.html").write_text(full)
    print("index.html", len(full), "bytes; index.artifact.html", len(t), "bytes")


if __name__ == "__main__":
    main()
