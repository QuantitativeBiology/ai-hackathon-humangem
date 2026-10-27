"""W2 media library tests (DESIGN §12 row W2): toy nutrient removal, every formulation string parsed, sheet facts,
line-assignment counts equal to the frozen design numbers, bound scheme, gate, shuffle, req_c != 0. Label-free.

Run: NEXTGEM_NO_TEST_RECORD=1 .venv/bin/python -m pytest tests/test_w2_media.py -q
"""

from __future__ import annotations

import json
import re

import cobra
import pandas as pd
import pytest

from nextgem import media, paths
from nextgem.knockouts import single_gene_knockouts

ON_DISK = set(media.SHEET_FILES)


# ------------------------------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def lines():
    return media.assign_line_media(ON_DISK)


@pytest.fixture(scope="module")
def sheets():
    return media.load_sheets()


def _all_texts() -> list[str]:
    dm = paths.depmap_dir()
    mc = pd.read_csv(dm / "ModelCondition.csv", usecols=["FormulationID"]).FormulationID
    mo = pd.read_csv(dm / "Model.csv", usecols=["FormulationID"]).FormulationID
    return sorted(set(mc.dropna()) | set(mo.dropna()))


ALL_TEXTS = _all_texts()


# ------------------------------------------------------------------------------------------------- strings
@pytest.mark.parametrize("text", ALL_TEXTS)
def test_every_formulation_string_parses(text):
    """Every distinct DepMap medium string (ModelCondition and Model.csv) parses with no unclassified token."""
    p = media.parse_formulation(text)
    assert p.base_class == media.base_class(text)
    if re.search(r"FBS|FCS|serum", text, re.I):
        assert p.serum_pct > 0
    else:
        assert p.serum_free
    for s in p.supplements:
        assert s.kind in {"nutrient", "non_nutrient", "cysteine_rule", "unrepresentable"}
        for c in s.components:
            assert c in media.EXCHANGES or c in media.NO_EXCHANGE, (text, c)
    if p.base_class in ("mixture_of_sheets",) or p.base_class.startswith("sheet:"):
        assert abs(sum(f for _, f in p.parts) - 1) < 1e-12


@pytest.mark.parametrize("text,expect", [
    ("RPMI + 10% FBS + 2mM Glutamine + 25mM HEPES + 25mM Sodium bicarbonate",
     dict(serum=10, comps={"glutamine": 2000.0, "HCO3-": 25000.0})),
    ("DMEM + 10% FBS+ 1% Glutamax", dict(serum=10, comps={"glutamine": None})),
    ("IMDM:RPMI (4:1) + 20% FBS + 2mM Glutamine", dict(serum=20, parts=[("IMDM", 0.8), ("RPMI", 0.2)])),
    ("MEM + 10% FBS + 4mM Glutamine + 1mM Sodium pyruvate + 2x MEM NEAA + EAA",
     dict(serum=10, comps={"glutamine": 4000.0, "pyruvate": 1000.0, "alanine": 200.0, "valine": None})),
    ("RPMI + 10% FBS + 50uM 2-Mercaptoethanol", dict(serum=10, cysrule=True)),
    ("DMEM:F12", dict(serum=0)),
    ("DMEM:F12 + 20% horse serum", dict(serum=20, flag="non_FBS_serum_treated_as_serum")),
    ("RPMI + 10% FBS + 30 mg/L Proline + 35 mg/L Cystine + 3.57 g/L HEPES +15 mg/L Hypoxanthine +1 mg/L ATP "
     "+10 mg/L Adenine + 1 mg/L Thymidine", dict(serum=10, comps={"hypoxanthine": pytest.approx(110.2, rel=1e-3),
                                                                "thymidine": pytest.approx(4.129, rel=1e-3)})),
])
def test_known_parses(text, expect):
    p = media.parse_formulation(text)
    assert p.serum_pct == expect["serum"]
    comps = {c: v for s in p.supplements for c, v in s.components.items()}
    for c, v in expect.get("comps", {}).items():
        assert comps[c] == v, (c, comps.get(c))
    if "parts" in expect:
        assert [(b, pytest.approx(f)) for b, f in expect["parts"]] == p.parts
    if expect.get("cysrule"):
        assert any(s.kind == "cysteine_rule" for s in p.supplements)
    if "flag" in expect:
        assert expect["flag"] in p.flags


def test_unknown_token_raises():
    with pytest.raises(ValueError):
        media.parse_formulation("RPMI + 10% FBS + 5mM Unobtainium")


# ------------------------------------------------------------------------------------------------- sheets
def test_sheet_facts(sheets):
    r, d, df, i = sheets["RPMI"], sheets["DMEM"], sheets["DMEM:F12"], sheets["IMDM"]
    assert r["asparagine"] == pytest.approx(378.788, rel=1e-4) and "asparagine" not in d
    assert "proline" not in d and "biotin" not in d and "cobalamin" not in d
    assert "pyridoxal" in i and "pyridoxine" not in i and "pyridoxine" in r
    for c in ("thymidine", "hypoxanthine", "lipoate", "linoleate", "putrescine", "cysteine"):
        assert c in df and c not in r
    assert r["GSH"] > 0 and "GSH" not in d
    assert all("cystine" in s for s in (r, d, df, i))
    assert "selenite" not in media.EXCHANGES and i.get("selenite", 0) > 0


def test_physiological_media(sheets):
    h, p = sheets["HPLM"], sheets["Plasmax"]
    assert "uridine" not in h  # Cantor 2017 Table S1 HPLM has uric acid but no uridine
    assert h["urate"] == 350 and p["uridine"] == 3 and p["urate"] == 270
    assert h["hypoxanthine"] == 10 and p["hypoxanthine"] == 5
    assert len(h) >= 60 and len(p) >= 70


def test_every_raw_component_classified():
    rows = media.read_compositions()
    comp = media.canonical_composition(rows)
    st = comp.canonical.map(media.mapping_status)
    assert set(st) <= {"mapped", "no_exchange", "non_nutrient"}


# ------------------------------------------------------------------------------------------------- assignment
def test_assignment_equals_design_numbers(lines):
    dn = json.loads((paths.results_dir() / "design" / "design_numbers.json").read_text())["medium"]
    assert lines.tier.value_counts().to_dict() == dn["tier_counts"]
    assert lines.base_class.value_counts().to_dict() == dn["base_class_counts"]
    att = lines[lines.attribution_set]
    assert len(att) == dn["attribution_set_screen_tier_on_disk"]
    assert att.base_class.value_counts().to_dict() == dn["attribution_set_by_base"]
    sheet_missing = lines.base_class.str.startswith("sheet:") & ~lines.base_class.isin({"sheet:" + b for b in ON_DISK})
    assert int(sheet_missing.sum()) == dn["fallback_lines_sheet_missing_any_tier"]


def test_strata_cover_all_lines(lines):
    assert lines.ModelID.is_unique and lines.stratum.notna().all()
    assert (lines[lines.stratum.isin(["imputed", "fallback_sheet_missing", "unmappable"])].medium_text_used
            == media.IMPUTED_TEXT).all()
    assert (lines.attribution_set == (lines.stratum == "attribution")).all()


def test_shuffle_preserves_strata(lines):
    lines = lines.copy()
    lines["medium_key"] = lines.medium_text_used
    sh = media.shuffled_assignment(lines, [11, 23])
    for col in ("shuffled_seed11", "shuffled_seed23"):
        for _, g in sh.groupby(["lineage", "serum_free"]):
            assert sorted(g[col]) == sorted(g.medium_key)
    again = media.shuffled_assignment(lines, [11])
    assert (again.shuffled_seed11.values == sh.shuffled_seed11.values).all()


# ------------------------------------------------------------------------------------------------- bounds (toy)
def _toy():
    """Toy: glutamine (G) or glucose-derived synthesis (gene SYN) -> P; floor lipoate L -> biomass; Asp import."""
    m = cobra.Model("toy")
    mets = {k: cobra.Metabolite(k, compartment="c") for k in ("glc", "gln", "P", "lip", "asp", "o2", "h2o", "bio")}
    ex = {}
    for canon, met in (("glucose", "glc"), ("glutamine", "gln"), ("lipoate", "lip"), ("aspartate", "asp"),
                       ("O2", "o2"), ("H2O", "h2o")):
        r = cobra.Reaction(media.EXCHANGES[canon][0], lower_bound=-1000, upper_bound=1000)
        r.add_metabolites({mets[met]: -1})
        ex[canon] = r
    syn = cobra.Reaction("SYN", upper_bound=1000)
    syn.add_metabolites({mets["glc"]: -1, mets["o2"]: -1, mets["P"]: 1})
    syn.gene_reaction_rule = "SYNG"
    use = cobra.Reaction("USE", upper_bound=1000)
    use.add_metabolites({mets["gln"]: -1, mets["P"]: 1})
    use.gene_reaction_rule = "TRANS"
    bio = cobra.Reaction("BIO", upper_bound=1000)
    bio.add_metabolites({mets["P"]: -1, mets["lip"]: -0.01, mets["bio"]: 1})
    sink = cobra.Reaction("EX_bio", lower_bound=0, upper_bound=1000)
    sink.add_metabolites({mets["bio"]: -1})
    m.add_reactions(list(ex.values()) + [syn, use, bio, sink])
    m.objective = "BIO"
    return m


def test_toy_nutrient_removal_flips_essentiality():
    m = _toy()
    pc = media.ParentConstraints(floor=("lipoate",), req={"lipoate": 0.1})
    with_gln = media.Medium("A", {"glucose": 1.0, "glutamine": 1.0}, serum=False)
    no_gln = media.Medium("B", {"glucose": 1.0}, serum=False)
    b1 = media.uptake_bounds(m, "M1", with_gln, pc)
    b0 = media.uptake_bounds(m, "M1", no_gln, pc)
    assert b1[media.EXCHANGES["lipoate"][0]] == pytest.approx(media.FLOOR_MULT * 0.1)
    assert media.EXCHANGES["glutamine"][0] not in b0
    r1 = single_gene_knockouts(m, medium=b1, objective="BIO").set_index("gene").ratio_r
    r0 = single_gene_knockouts(m, medium=b0, objective="BIO").set_index("gene").ratio_r
    assert r1["SYNG"] > 0.01 and r0["SYNG"] == 0.0  # removing glutamine makes the synthesis gene essential
    assert r0["TRANS"] == 1.0


def test_floor_binding_detected_on_toy():
    m = _toy()
    med = media.Medium("A", {"glucose": 1000.0, "glutamine": 1000.0}, serum=False)
    req = 0.1
    pc = media.ParentConstraints(floor=("lipoate",), req={"lipoate": req})
    with m:
        m.medium = media.uptake_bounds(m, "M1", med, pc)
        g10 = m.slim_optimize()
        m.medium = media.uptake_bounds(m, "M1", med, pc, floor_override=1000)
        gop = m.slim_optimize()
    assert g10 == pytest.approx(100.0) and gop > 2 * g10  # 10 x req binds here: the check must see it
    assert abs(g10 - gop) / gop > media.FLOOR_NONBINDING_RTOL


def test_bound_scheme_variants():
    m = _toy()
    pc = media.ParentConstraints(floor=("lipoate",), req={"lipoate": 0.1, **{c: 2.0 for c in media.SERUM_F_SCALED}})
    med = media.Medium("A", {"glucose": 1.0, "aspartate": 1.0, "carnitine": 5.0}, serum=True)
    sf = media.Medium("B", {"glucose": 1.0}, serum=False)
    b2 = media.uptake_bounds(m, "M2", med, pc, f=0.01)
    assert b2[media.EXCHANGES["cholesterol"][0]] == pytest.approx(0.02)
    assert b2[media.EXCHANGES["carnitine"][0]] == media.UPTAKE_OPEN  # supplied by the medium itself
    assert b2[media.EXCHANGES["Fe3+"][0]] == media.UPTAKE_OPEN
    assert media.uptake_bounds(m, "M2", sf, pc, f=1.0) == media.uptake_bounds(m, "M1", sf, pc)  # serum-free
    g = media.uptake_bounds(m, "M2", med, pc, f=0.01, gate_open=False)
    assert media.EXCHANGES["aspartate"][0] not in g and media.EXCHANGES["aspartate"][0] in b2
    b3 = media.uptake_bounds(m, "M3", med, pc, f=0.01, calibration={"glucose": 0.3, "glutamine": 0.05, "O2": 1.0})
    assert b3[media.EXCHANGES["glucose"][0]] == 0.3
    assert len(media.uptake_bounds(m, "M-open")) == len(m.exchanges)
    cys = media.Medium("C", {"cysteine": 1.0, "cystine": 1.0}, serum=False)
    assert media.EXCHANGES["cysteine"][0] not in media.uptake_bounds(m, "M1", cys, pc)
    cys.allow_cysteine = True
    assert media.EXCHANGES["cysteine"][0] in media.uptake_bounds(m, "M1", cys, pc)


def test_dmem_defaults_and_mixture(sheets):
    p = media.parse_formulation("DMEM + 10% FBS")
    med = media.build_medium(p, sheets)
    assert "glutamine" in med.present() and "pyruvate" not in med.present()
    assert "dmem_pyruvate_unstated_default_none" in med.flags
    mix = media.build_medium(media.parse_formulation("DMEM:RPMI (1:1) + 10% FBS"), sheets)
    assert mix.components["asparagine"] == pytest.approx(sheets["RPMI"]["asparagine"] / 2)
    df12 = media.build_medium(media.parse_formulation("DMEM:F12 + 10% FBS"), sheets)
    assert df12.allow_cysteine  # documented free cysteine (CYSTEINE_RULE = 'documented')


# ------------------------------------------------------------------------------------------------- outputs
def test_req_nonzero_unless_flagged():
    """DESIGN W2: a unit test fails if any serum component has req_c = 0 unless it is flagged."""
    f = paths.results_dir() / "w2" / "floor" / "requirements.tsv"
    if not f.exists():
        pytest.skip("run scripts/w2_medium_floor.py first")
    rq = pd.read_csv(f, sep="\t")
    flagged = rq.get("flagged", pd.Series(False, index=rq.index)).astype("boolean").fillna(False)
    assert ((rq.req_c > 0) | flagged).all(), rq[~((rq.req_c > 0) | flagged)]
    assert set(media.SERUM_F_SCALED) <= set(rq.component)
    assert "cholesterol" not in set(rq[rq.kind == "floor"].component)
