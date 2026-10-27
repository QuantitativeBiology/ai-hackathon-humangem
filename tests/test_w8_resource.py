"""W8 tests (DESIGN.md section 12 row W8): toy SBML round-trip, GR passes, corrupted SBML rejected; plus
report_build wording rules. GLPK only; label-free."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import cobra
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from nextgem import resource as rs  # noqa: E402

import report_build as rb  # noqa: E402


@pytest.fixture(autouse=True)
def _glpk():
    old = cobra.Configuration().solver
    cobra.Configuration().solver = "glpk"
    yield
    cobra.Configuration().solver = old


def toy_model() -> cobra.Model:
    m = cobra.Model("toy")
    mets = {n: cobra.Metabolite(n, compartment="c") for n in "ABC"}
    def rxn(i, stoich, lb, ub, gpr=""):
        r = cobra.Reaction(i, lower_bound=lb, upper_bound=ub)
        m.add_reactions([r])
        r.add_metabolites({mets[k]: v for k, v in stoich.items()})
        if gpr:
            r.gene_reaction_rule = gpr
        return r
    rxn("EX_A", {"A": -1}, -10, 0)
    rxn("R1", {"A": -1, "B": 1}, 0, 1000, "g1")
    rxn("R2", {"A": -1, "B": 1}, 0, 1000, "g2")
    rxn("R3", {"B": -1, "C": 1}, 0, 1000, "g3 and g4")
    rxn("BIO", {"C": -1}, 0, 1000)
    m.objective = "BIO"
    return m


META = rs.Meta(model_id="toy", parent="toyparent", parent_version="0.1", pipeline_id="none",
               git_commit="deadbeef", input_shas={"x": "00"}, repaired_reactions=["R2"],
               medium_source_tier="M0", release_status="demo, not a release")


def test_w8_toy_sbml_roundtrip(tmp_path):
    m = toy_model()
    p = rs.export_model(m, META, tmp_path)
    assert p.name == "toy_toyparent-0.1_none.xml.gz"
    assert gzip.open(p).read(5) == b"<?xml"
    b = rs.load_checked(p)
    assert rs.roundtrip_diff(m, b) == {}
    notes = rs.read_notes(b)
    assert notes["git_commit"] == "deadbeef" and notes["repaired_reactions"] == "R2"
    assert notes["medium_source_tier"] == "M0" and notes["pipeline_id"] == "none"
    assert '"x": "00"' in notes["input_shas"] and notes["release_status"] == "demo, not a release"
    v = rs.validate_model(m, p, None, None)
    assert {k: g[0] for k, g in v.gates.items()} == {
        "libsbml_no_new_errors": True, "roundtrip_identical": True, "wt_growth_match": True,
        "random_ko_match": True, "viability_tasks": "skipped"}
    # export leaves the source model untouched
    assert "release_status" not in m.notes and m.id == "toy"


def test_w8_roundtrip_detects_bound_change(tmp_path):
    m = toy_model()
    p = rs.export_model(m, META, tmp_path)
    b = rs.load_checked(p)
    b.reactions.R1.upper_bound = 5
    assert rs.roundtrip_diff(m, b)["bound_mismatch"] == 1
    v = rs.validate_model(b, p, None, None)  # "source" no longer equals the file
    assert v.status == "fail" and "roundtrip_identical" in v.detail


def test_w8_corrupted_sbml_rejected(tmp_path):
    m = toy_model()
    p = rs.export_model(m, META, tmp_path)
    raw = gzip.decompress(p.read_bytes())
    trunc = tmp_path / "trunc.xml.gz"
    trunc.write_bytes(gzip.compress(raw[: len(raw) // 2]))
    garbage = tmp_path / "garbage.xml.gz"
    garbage.write_bytes(b"not gzip at all")
    badgz = tmp_path / "badgz.xml.gz"
    badgz.write_bytes(p.read_bytes()[:-40])
    for q in (trunc, garbage, badgz):
        with pytest.raises(rs.ResourceError):
            rs.load_checked(q)
        v = rs.validate_model(m, q, None, None)
        assert v.status == "fail" and "load" in v.gates


def test_w8_libsbml_new_errors_vs_parent():
    from collections import Counter
    assert rs.new_errors(Counter({1: 2}), Counter({1: 2})) == {}
    assert rs.new_errors(Counter({1: 3, 5: 1}), Counter({1: 2})) == {1: 1, 5: 1}


def test_w8_random_ko_detects_difference():
    a, b = toy_model(), toy_model()
    b.reactions.R3.bounds = (0, 0)  # g3/g4 KO matters in a; b already lethal -> WT differs
    out = rs.random_ko_compare(a, b, n=4)
    assert out["mismatches"]


def test_w8_line_model_from_mask():
    import numpy as np
    m = toy_model()
    bits = np.array([1, 1, 0, 1, 1], dtype=bool)
    lm, n = rs.line_model_from_mask(m, np.packbits(bits).tobytes().hex())
    assert n == 1 and lm.reactions.R2.bounds == (0.0, 0.0) and m.reactions.R2.upper_bound == 1000


@pytest.fixture(scope="module")
def human2_ham():
    old = cobra.Configuration().solver
    cobra.Configuration().solver = "glpk"
    from nextgem.knockouts import configure
    from nextgem.model import load_model
    m = load_model("human2")
    configure(m)
    yield m
    cobra.Configuration().solver = old


def test_w8_gr_passes(human2_ham):
    from nextgem.model import GEMS
    res = rs.viability_tasks(human2_ham, GEMS["human2"].essential_tasks, types=("GR",))
    assert res == {"GR": True}


# ---- report_build wording rules --------------------------------------------------------------------------
def test_w8_claim_tiers():
    f = rb.claim_tier
    assert f("C1", passes_test=False, point=0.05, ci_lower=0.02) == "no detectable improvement"
    assert f("C1", passes_test=True, point=0.005, ci_lower=0.001) == "detectable improvement"
    assert f("C1", passes_test=True, point=0.02, ci_lower=0.005) == "improves"
    assert f("C1", passes_test=True, point=0.02, ci_lower=0.012) == "meaningfully improves"
    assert f("C3", passes_test=True, point=0.5, ci_lower=0.4) == "detectable improvement"  # no SESOI
    assert f("C1b", passes_test=True, point=0.004, ci_lower=0.0031, sesoi=0.003) == "meaningfully improves"


def test_w8_sesoi_c1b():
    assert rb.sesoi_c1b(0.1, 0.05) == pytest.approx(0.005)
    with pytest.raises(ValueError):
        rb.sesoi_c1b(0.0, 0.0)


def test_w8_transfer_wording():
    assert rb.transfer_wording(0.01, 5).startswith("generalises")
    assert rb.transfer_wording(0.01, 4).startswith("no claim")
    assert rb.transfer_wording(-0.001, 7).startswith("no claim")


def test_w8_headroom_rule():
    assert rb.headroom_closed(0.01, 0.1, (0.05, 0.15)) == pytest.approx(0.1)
    assert rb.headroom_closed(0.01, 0.01, (-0.01, 0.03)) is None


def test_w8_flags():
    assert rb.derive_flags(gate_passed=False, exploratory=True) == ["gate not passed", "exploratory"]
    assert rb.apply_flags("x.", ["exploratory", "gate not passed"]) == "x. [gate not passed; exploratory]"
    assert rb.apply_flags("x.", []) == "x."
    with pytest.raises(ValueError):
        rb.apply_flags("x", ["bogus"])
    s = rb.render_claim("C1", "L-prob", 0.02, (0.012, 0.03), passes_test=True, flags=["label-calibrated"])
    assert "meaningfully improves" in s and s.endswith("[label-calibrated]")


def test_w8_collect_numbers_stub(tmp_path):
    out = rb.collect_numbers(tmp_path)
    assert set(out) == set(rb.NUMBER_SOURCES) and not any(out.values())
    with pytest.raises(NotImplementedError):
        rb.render("", {})


def test_w8_rel_equal_noise_floor():
    assert rs._rel_equal(1.7e-12, 0.0) and rs._rel_equal(46.69420500301919, 46.69420500301934)
    assert not rs._rel_equal(1.0, 1.0 + 1e-5) and not rs._rel_equal(1e-6, 0.0)
