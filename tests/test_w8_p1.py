"""W8 P1 resource tests (DESIGN §9 W8; RESEARCH_TREE D11): toy medium + biomass-variant build, notes, reload gate,
line manifest, per-line materialisation, task keys; plus label-free checks on the real line assignment and the
Human2 viability tasks. GLPK only; no label is read."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path
from types import SimpleNamespace

import cobra
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from nextgem import resource as rs  # noqa: E402
from nextgem.run import sha256_file  # noqa: E402
from nextgem.w3 import Variant  # noqa: E402


@pytest.fixture(autouse=True)
def _glpk():
    old = cobra.Configuration().solver
    cobra.Configuration().solver = "glpk"
    yield
    cobra.Configuration().solver = old


def toy_parent() -> cobra.Model:
    """A_e, D_e taken up; biomass MAR13082 = B + P; B from A (isozymes g1/g2), P from D (g5)."""
    m = cobra.Model("toyparent")
    mets = {n: cobra.Metabolite(n, compartment=n.split("_")[1]) for n in ("A_e", "D_e", "A_c", "D_c", "B_c", "P_c")}

    def rxn(i, stoich, lb, ub, gpr=""):
        r = cobra.Reaction(i, lower_bound=lb, upper_bound=ub)
        m.add_reactions([r])
        r.add_metabolites({mets[k]: v for k, v in stoich.items()})
        if gpr:
            r.gene_reaction_rule = gpr
        return r

    rxn("EX_A", {"A_e": -1}, -1000, 1000)
    rxn("EX_D", {"D_e": -1}, -1000, 1000)
    rxn("TA", {"A_e": -1, "A_c": 1}, 0, 1000, "t1")
    rxn("TD", {"D_e": -1, "D_c": 1}, 0, 1000, "t2")
    rxn("R1", {"A_c": -1, "B_c": 1}, 0, 1000, "g1")
    rxn("R2", {"A_c": -1, "B_c": 1}, 0, 1000, "g2")
    rxn("R5", {"D_c": -1, "P_c": 1}, 0, 1000, "g5")
    rxn("CLOSED", {"A_c": -1, "D_c": 1}, 0, 0, "g6")
    rxn("MAR13082", {"B_c": -1, "P_c": -1}, 0, 1000)
    m.objective = "MAR13082"
    return m


BLOO = Variant("B-loo:toy_pool", "remove", {"mets": ["P_c"]})
MEDIUM = {"EX_A": 10.0}  # D not in the medium


def p1_meta(**kw) -> rs.Meta:
    base = dict(model_id="medium-abc123", parent="toyparent", parent_version="0.1", pipeline_id=rs.PIPELINE_P1,
                git_commit="deadbeef", input_shas={"x": "00"}, n_repaired_reactions=0,
                repaired_note="none (generic parent)", flags=rs.p1_flags(False), release_status="pilot, not a release",
                extra={"pipeline": rs.P1_DESCRIPTION, "medium_signature": "abc123", "w3_variant": BLOO.to_json(),
                       "human_gem_release": "toy 0.1"})
    base.update(kw)
    return rs.Meta(**base)


def test_w8_p1_flags():
    assert rs.p1_flags(True) == "W3-D1;W2-SHUF;W2-MOPEN"
    assert rs.p1_flags(False) == "W3-D1;W3-MEV;W2-SHUF;W2-MOPEN"


def test_w8_p1_build_medium_variant():
    parent = toy_parent()
    with parent:  # the plain parent on the medium cannot make P: no growth
        parent.medium = MEDIUM
        assert parent.slim_optimize() == pytest.approx(0.0, abs=1e-9)
    m = rs.build_medium_variant_model(parent, MEDIUM, BLOO)
    assert m.slim_optimize() == pytest.approx(10.0)
    assert m.reactions.EX_D.bounds == (0.0, 1000.0) and m.reactions.EX_A.lower_bound == -10.0
    assert "P_c" not in {x.id for x in m.reactions.MAR13082.metabolites}
    # parent untouched; the closed reaction is kept with bounds (0, 0)
    assert parent.reactions.EX_D.lower_bound == -1000 and "P_c" in {x.id for x in parent.reactions.MAR13082.metabolites}
    assert rs.n_closed(m) == 1 and m.reactions.CLOSED.bounds == (0.0, 0.0)


def test_w8_p1_notes_roundtrip(tmp_path):
    m = rs.build_medium_variant_model(toy_parent(), MEDIUM, BLOO)
    p = rs.export_model(m, p1_meta(), tmp_path)
    assert p.name == "medium-abc123_toyparent-0.1_P1.xml.gz"
    notes = rs.read_notes(rs.load_checked(p))
    assert notes["pipeline_id"] == "P1" and notes["pipeline"].startswith("P1 (RESEARCH_TREE D11)")
    assert notes["flags"] == "W3-D1;W3-MEV;W2-SHUF;W2-MOPEN" and notes["medium_signature"] == "abc123"
    assert notes["repaired_reactions"] == "none (generic parent)" and notes["n_repaired_reactions"] == "0"
    assert '"remove"' in notes["w3_variant"] and notes["git_commit"] == "deadbeef"
    with pytest.raises(ValueError):
        rs.set_notes(m, p1_meta(extra={"flags": "x"}))


def test_w8_reload_sample_deterministic():
    ref = {f"g{i:02d}": (0.5 if i < 10 else 1.0) for i in range(40)}
    s1, s2 = rs.reload_sample(ref, 4, 3), rs.reload_sample(ref, 4, 3)
    assert s1 == s2 and len(s1) == 7
    assert sum(ref[g] < 1 for g in s1) == 4 and sum(ref[g] == 1 for g in s1) == 3
    assert len(rs.reload_sample({"a": 0.0}, 25, 15)) == 1


def test_w8_reload_check(tmp_path):
    m = rs.build_medium_variant_model(toy_parent(), MEDIUM, BLOO)
    genes = sorted(g.id for g in m.genes)
    wt, r = rs.ko_ratio_vector(m, genes)
    assert wt == pytest.approx(10.0)
    ref = dict(zip(genes, r))
    assert ref["t1"] == 0.0 and ref["g1"] == 1.0 and ref["g5"] == 1.0  # isozyme g1/g2; P no longer required
    p = rs.export_model(m, p1_meta(), tmp_path)
    b = rs.load_checked(p)
    ok, d = rs.reload_check(m, b, ref, round(wt, 4), genes)
    assert ok and d["max_abs_diff_reloaded_m6"] == 0.0
    bad = dict(ref, g1=0.5)
    ok, d = rs.reload_check(m, b, bad, wt, genes)
    assert not ok and d["mismatches"][0][0] == "g1"
    ok, _ = rs.reload_check(m, b, ref, wt * (1 + 1e-3), genes)  # WT off the stored precision
    assert not ok


def test_w8_validate_extra_gate_and_cells(tmp_path):
    m = rs.build_medium_variant_model(toy_parent(), MEDIUM, BLOO)
    p = rs.export_model(m, p1_meta(), tmp_path)
    v = rs.validate_model(m, p, None, None, extra_gates={"reload_vs_m6": lambda a, b: (False, {"why": "toy"})})
    assert v.status == "fail" and "reload_vs_m6" in v.detail and "viability_tasks" in v.detail
    cells = rs.gate_cells(v)
    assert cells["gate_reload_m6"] == "FAIL" and cells["gate_viability"] == "skipped"
    assert cells["gate_libsbml"] == cells["gate_roundtrip"] == cells["gate_wt_growth"] == cells["gate_random_ko"] == "pass"


def test_w8_task_keys():
    ts = [SimpleNamespace(id=i) for i in ("ER", "ER", "GR", "ER")]
    assert rs.task_keys(ts) == ["ER1", "ER2", "GR", "ER3"]


def _lines():
    return pd.DataFrame({
        "ModelID": ["ACH-1", "ACH-2", "ACH-3", "ACH-1"], "parent": ["human1", "human1", "human1", "human2"],
        "w2_sweep_id": ["human1:s1", "human1:s1", "human1:s2", "human2:s1"],
        "medium_key": ["RPMI", "RPMI", "DMEM", "RPMI"],
        "medium_stratum": ["attribution", "imputed", "attribution", "attribution"]})


def _models():
    return pd.DataFrame({
        "parent": ["human1", "human1", "human2"], "w2_sweep_id": ["human1:s1", "human1:s2", "human2:s1"],
        "file": ["v1/a.xml.gz", "v1/b.xml.gz", "v1/c.xml.gz"], "sha256": ["sa", "sb", "sc"],
        "version": ["1.19.0", "1.19.0", "2.1.1"], "pipeline_id": ["P1"] * 3, "medium_signature": ["s1", "s2", "s1"],
        "flags": ["f"] * 3, "validation_status": ["pass"] * 3, "release_status": ["r"] * 3})


def test_w8_line_manifest():
    man = rs.line_manifest(_lines(), _models())
    assert list(man.columns) == rs.LINE_COLUMNS and len(man) == 4
    assert man.parent.tolist() == ["human2", "human1", "human1", "human1"]  # Human2 first
    h1 = man[man.parent == "human1"].set_index("ModelID")
    assert h1.loc["ACH-2", "model_file"] == "v1/a.xml.gz" and h1.loc["ACH-3", "model_sha256"] == "sb"
    assert h1.loc["ACH-2", "medium_fallback"].startswith("classic RPMI") and h1.loc["ACH-1", "medium_fallback"] == ""
    with pytest.raises(rs.ResourceError):  # a line without a model file
        rs.line_manifest(_lines(), _models().iloc[:2])
    with pytest.raises(rs.ResourceError):  # two files for one signature
        rs.line_manifest(_lines(), pd.concat([_models(), _models().iloc[:1]]))
    with pytest.raises(rs.ResourceError):  # duplicated (line, parent)
        rs.line_manifest(pd.concat([_lines(), _lines().iloc[:1]]), _models())


def test_w8_load_line_and_materialise(tmp_path):
    res = tmp_path / "resource"
    m = rs.build_medium_variant_model(toy_parent(), MEDIUM, BLOO)
    p = rs.export_model(m, p1_meta(), res / "v9.9.9")
    models = pd.DataFrame({"parent": ["toyparent"], "w2_sweep_id": ["toyparent:abc123"],
                           "file": [str(p.relative_to(res))], "sha256": [sha256_file(p)], "version": ["0.1"],
                           "pipeline_id": ["P1"], "medium_signature": ["abc123"], "flags": ["f"],
                           "validation_status": ["pass"], "release_status": ["pilot"]})
    lines = pd.DataFrame({"ModelID": ["ACH-000009"], "parent": ["toyparent"], "w2_sweep_id": ["toyparent:abc123"],
                          "medium_key": ["RPMI"], "medium_stratum": ["unmappable"], "lineage": ["Lung"]})
    rs.line_manifest(lines, models).to_csv(res / "manifest.tsv", sep="\t", index=False)
    lm = rs.load_line(res, "ACH-000009", "toyparent")
    notes = rs.read_notes(lm)
    assert notes["model_id"] == "ACH-000009" and notes["shared_model_id"] == "medium-abc123"
    assert notes["line_medium_stratum"] == "unmappable" and notes["line_medium_fallback"].startswith("classic RPMI")
    assert lm.id == "ACH_000009_toyparent_P1" and rs.roundtrip_diff(m, lm) == {}
    out = rs.materialise_line(res, "ACH-000009", "toyparent", tmp_path / "lines")
    assert out.name == "ACH-000009_toyparent-0.1_P1.xml.gz" and gzip.open(out).read(5) == b"<?xml"
    b = rs.load_checked(out)
    assert rs.roundtrip_diff(m, b) == {} and rs.read_notes(b)["line_lineage"] == "Lung"
    assert b.slim_optimize() == pytest.approx(m.slim_optimize())
    p.write_bytes(p.read_bytes()[:-10] + b"0" * 10)  # corrupt the shared file: sha256 check rejects it
    with pytest.raises(rs.ResourceError):
        rs.load_line(res, "ACH-000009", "toyparent")


# ---- label-free checks on the real inputs ---------------------------------------------------------------
def test_w8_line_table_covers_all_lines():
    import w8_build_p1 as wb

    t, _ = wb.line_table()
    assert t.groupby("parent").ModelID.nunique().to_dict() == {"human1": 1361, "human2": 1361}
    assert not t.duplicated(["ModelID", "parent"]).any()
    assert t.groupby("parent").w2_sweep_id.nunique().to_dict() == {"human1": 20, "human2": 20}
    assert t.medium_key.notna().all() and t.medium_stratum.notna().all()
    assert set(t.cancer) <= {"true", "false", "unknown"}
    assert "tier" not in t.columns  # no split tier in the resource tables


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


def test_w8_viability_all_er_tasks(human2_ham):
    """All 9 ER tasks and GR are checked and kept (raven ids are categories; the stub kept only the last ER)."""
    from nextgem.model import GEMS
    res = rs.viability_tasks(human2_ham, GEMS["human2"].essential_tasks)
    assert sorted(res) == sorted(["GR"] + [f"ER{i}" for i in range(1, 10)])
    assert all(res.values())


def test_w8_canonical_group_members(tmp_path):
    """Group members are written in a stable order (cobra keeps them in a set): identical models give identical
    bytes whatever the member insertion order; membership survives the round trip."""
    from cobra.core import Group

    files = []
    for order in (["R1", "R2", "R5", "TA"], ["TA", "R5", "R2", "R1"]):
        m = rs.build_medium_variant_model(toy_parent(), MEDIUM, BLOO)
        m.add_groups([Group("grp1", name="toy subsystem", members=[m.reactions.get_by_id(r) for r in order])])
        files.append(rs.export_model(m, p1_meta(), tmp_path / "".join(order)))
    assert files[0].read_bytes() == files[1].read_bytes()
    raw = gzip.decompress(files[0].read_bytes())
    i = [raw.index(f'idRef="R_{r}"'.encode()) for r in ("R1", "R2", "R5", "TA")]
    assert i == sorted(i)
    b = rs.load_checked(files[0])
    assert {r.id for r in b.groups.get_by_id("grp1").members} == {"R1", "R2", "R5", "TA"}
    assert rs.canonical_group_members(b"<x/>") == b"<x/>"


def test_w8_mevalonate_rule_matches_m6():
    """The W3-MEV rule used for new signatures reproduces results/m6/scores/w3_mev.tsv on every M6 P1 sweep."""
    import m6_common as mc
    import m6_sweeps
    import w8_extend_unscreened as we

    mev = pd.read_csv(mc.SCDIR / "w3_mev.tsv", sep="\t").set_index("key")
    cat = pd.read_csv(mc.SWDIR / "catalog.tsv", sep="\t")
    n = 0
    for p in ("human1", "human2"):
        ids = cat[(cat.parent == p) & (cat.variant == mc.BLOO) & (cat.step == mc.M2)].w2_sweep_id
        new = pd.read_parquet(mc.SWDIR / f"ratio_r_{p}.parquet").set_index("sweep")
        reused = m6_sweeps.reused_rows()[p]
        for sid in ids:
            k = mc.sweep_key(mc.BLOO, sid)
            row = new.loc[k] if k in new.index else reused.loc[k]
            assert we.mevalonate_core_lethal(p, row.to_dict()) == bool(mev.loc[k, "core_lethal"]), k
            n += 1
    assert n == 40


def test_w8_unscreened_assignment_tables():
    """The committed W8 extension tables: 753 in-scope unscreened lines, W2 tiers only (onboarded/imputed), one
    signature per (parent, medium key), no overlap with the screened cohort."""
    out = ROOT / "results" / "w8" / "unscreened"
    if not (out / "line_media_unscreened.tsv").exists():
        pytest.skip("extension not built")
    lines = pd.read_csv(out / "line_media_unscreened.tsv", sep="\t")
    sigs = pd.read_csv(out / "signatures.tsv", sep="\t")
    la = pd.read_csv(ROOT / "results" / "m6" / "scores" / "line_assignment.tsv", sep="\t")
    assert lines.ModelID.is_unique and len(lines) == 753
    assert not set(lines.ModelID) & set(la.ModelID)
    assert set(lines.tier) <= {"onboarded", "imputed"} and (lines.screened.astype(str).str.lower() == "false").all()
    assert set(lines[lines.tier == "imputed"].medium_text_used) == {"RPMI + 10% FBS"}
    assert not sigs.duplicated(["parent", "medium_key"]).any()
    assert set(sigs.medium_key) == set(lines.medium_key)


def test_w8_gate_exit_code(capsys):
    ok = pd.DataFrame({"parent": ["human2"] * 2, "w2_sweep_id": ["a", "b"], "job_status": ["ok", "ok"],
                       "validation_status": ["pass", "pass"], "validation_detail": ["", ""]})
    assert rs.gate_exit_code(ok) == 0
    one = ok.assign(validation_status=["pass", "fail"], validation_detail=["", "random_ko_match"])
    assert rs.gate_exit_code(one) == 3 and "random_ko_match" in capsys.readouterr().err
    assert rs.gate_exit_code(ok.assign(job_status=["ok", "failed-timeout"])) == 3
    assert rs.gate_exit_code(ok.iloc[:0]) == 3


def test_w8_build_main_propagates_exit_code(monkeypatch):
    """main() returns run()'s code through the provenance block, and the script exits with it."""
    import contextlib

    import w8_build_p1 as wb

    @contextlib.contextmanager
    def fake_prov(*a, **k):
        yield None

    monkeypatch.setattr("nextgem.run.provenance", fake_prov)
    monkeypatch.setattr(wb, "run", lambda *a, **k: 3)
    assert wb.main(["--pilot", "1", "--allow-dirty", "--parents", "human2"]) == 3
    monkeypatch.setattr(wb, "run", lambda *a, **k: 0)
    assert wb.main(["--pilot", "1", "--allow-dirty", "--parents", "human2"]) == 0


def test_w8_extension_build_idempotent():
    """A re-run of the extension build drops its earlier rows first (models of the rebuilt signatures, unscreened
    manifest rows) and keeps the screened cohort untouched."""
    import w8_extend_unscreened as we

    models = pd.DataFrame({"parent": ["human2", "human2", "human1"], "w2_sweep_id": ["human2:a", "human2:new", "human1:b"]})
    man = pd.DataFrame({"ModelID": ["ACH-1", "ACH-2", "ACH-2"], "parent": ["human2", "human2", "human1"],
                        "screened": ["true", "false", "false"]})
    new = pd.DataFrame({"parent": ["human2"], "w2_sweep_id": ["human2:new"]})
    m2, man2 = we.drop_previous_extension(models, man, new)
    assert m2.w2_sweep_id.tolist() == ["human2:a", "human1:b"] and man2.ModelID.tolist() == ["ACH-1"]
