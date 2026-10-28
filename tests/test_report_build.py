"""Known-answer tests for the REPORT.md number checker (scripts/report_build.py).

All tests use toy files under tmp_path; none reads a results table, a label file or the real REPORT.md, so their
outcome cannot change when results are added. Names are stable (results/test_runs.tsv tracks them).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import report_build as rb  # noqa: E402

M = "−"


def _toy(tmp_path: Path) -> rb.Source:
    (tmp_path / "res").mkdir()
    (tmp_path / "res" / "a.json").write_text('{"delta": 0.05174, "ci": [0.04988, 0.05376], "n": 1361}')
    (tmp_path / "res" / "t.tsv").write_text("method\ttap\nH1\t0.301676\nB2\t0.891696\n")
    return rb.Source(tmp_path, ref=None)


def _numbers():
    return {
        "a.delta": rb.Number(("res/a.json",), lambda s: rb.sci(rb.J(s, "res/a.json", "delta"),
                                                               *rb.J(s, "res/a.json", "ci"))),
        "a.n": rb.Number(("res/a.json",), lambda s: rb.ints(rb.J(s, "res/a.json", "n"))),
        "t.h1": rb.Number(("res/t.tsv",), lambda s: rb.fx(rb.T(s, "res/t.tsv", "tap", method="H1"))),
        "missing": rb.Number(("res/none.json",), lambda s: rb.fx(rb.J(s, "res/none.json", "x"))),
    }


def test_report_build_formatters():
    assert rb.fx(-0.0071, 4) == f"{M}0.0071"
    assert rb.sx(0.05174) == "+0.0517" and rb.sx(-0.00708) == f"{M}0.0071"
    assert rb.ci(-0.0079, -0.0063) == f"[{M}0.0079, {M}0.0063]"
    assert rb.sci(0.0322, 0.0308, 0.0334) == "+0.0322 [0.0308, 0.0334]"
    assert rb.ints(1361) == "1,361" and rb.pct(0.4232) == "42%" and rb.pct(0.02724, 1) == "2.7%"


def test_report_build_shape_regex_matches_changed_numbers():
    rx = rb.shape_regex("+0.0517 [0.0499, 0.0538]")
    assert rx.match(f"{M}0.0100 [{M}0.0200, 0.0010] rest").group(0) == f"{M}0.0100 [{M}0.0200, 0.0010]"
    assert rx.match("+0.05 [0.04, 0.06]") is not None
    assert rb.shape_regex("Human1").match("Human2") is None  # words are literal


def test_report_build_check_passes_on_match(tmp_path):
    src = _toy(tmp_path)
    text = ("Lines: <!--n:a.n-->1,361 and P1 <!--n:a.delta-->+0.0517 [0.0499, 0.0538] (res/a.json). "
            "G* <!--n:t.h1-->0.302 (res/t.tsv).\n")
    issues, n = rb.check_text(text, src, _numbers())
    assert n == 3 and issues == []


def test_report_build_check_flags_mismatch(tmp_path):
    src = _toy(tmp_path)
    issues, _ = rb.check_text("P1 <!--n:a.delta-->+0.0520 [0.0499, 0.0538] (res/a.json).\n", src, _numbers())
    assert [i.kind for i in issues] == ["mismatch"]
    assert "+0.0520 [0.0499, 0.0538]" in issues[0].detail and "+0.0517" in issues[0].detail


def test_report_build_check_rejects_longer_number(tmp_path):
    """0.302 must not pass when the text says 0.3021."""
    src = _toy(tmp_path)
    issues, _ = rb.check_text("G* <!--n:t.h1-->0.3021 (res/t.tsv)\n", src, _numbers())
    assert [i.kind for i in issues] == ["mismatch"]


def test_report_build_sentence_period_is_not_part_of_number(tmp_path):
    src = _toy(tmp_path)
    issues, _ = rb.check_text("G* is <!--n:t.h1-->0.302. Next sentence (res/t.tsv).\n", src, _numbers())
    assert issues == []


def test_report_build_check_flags_unknown_key(tmp_path):
    issues, n = rb.check_text("x <!--n:nope-->1\n", _toy(tmp_path), _numbers())
    assert n == 1 and [i.kind for i in issues] == ["unknown-key"]


def test_report_build_check_flags_uncited_source(tmp_path):
    src = _toy(tmp_path)
    text = "G* <!--n:t.h1-->0.302.\n\nSource: res/t.tsv\n"  # cited in another paragraph only
    issues, _ = rb.check_text(text, src, _numbers())
    assert [i.kind for i in issues] == ["uncited"]


def test_report_build_table_caption_counts_as_citation(tmp_path):
    src = _toy(tmp_path)
    text = ("Values from res/t.tsv.\n\n| m | tAP |\n|---|---|\n| H1 | <!--n:t.h1-->0.302 |\n\n"
            "Note: none.\n")
    issues, _ = rb.check_text(text, src, _numbers())
    assert issues == []


def test_report_build_indented_table_in_list_uses_caption(tmp_path):
    src = _toy(tmp_path)
    text = ("- Estimates (res/t.tsv).\n\n  | m | tAP |\n  |---|---|\n  | H1 | <!--n:t.h1-->0.302 |\n\n- next item\n")
    issues, _ = rb.check_text(text, src, _numbers())
    assert issues == []


def test_report_build_source_error_is_an_issue(tmp_path):
    issues, _ = rb.check_text("x <!--n:missing-->? (res/none.json)\n", _toy(tmp_path), _numbers())
    assert [i.kind for i in issues] == ["source-error"]


def test_report_build_update_fills_slots_and_rewrites_values(tmp_path):
    src = _toy(tmp_path)
    text = "a <!--n:a.delta-->? and <!--n:t.h1-->0.999. (res/a.json, res/t.tsv)\n"
    new, issues = rb.update_text(text, src, _numbers())
    assert issues == []
    assert new == "a <!--n:a.delta-->+0.0517 [0.0499, 0.0538] and <!--n:t.h1-->0.302. (res/a.json, res/t.tsv)\n"
    assert rb.check_text(new, src, _numbers())[0] == []


def test_report_build_update_keeps_comma_after_integer_value(tmp_path):
    """Regression: a value like 19/20 followed by ', ' must not swallow the comma (commas are thousands only)."""
    (tmp_path / "f.json").write_text('{"p": 19, "n": 20, "big": 2722}')
    src = rb.Source(tmp_path, ref=None)
    nums = {"r": rb.Number(("f.json",), lambda s: f"{rb.J(s, 'f.json', 'p')}/{rb.J(s, 'f.json', 'n')}"),
            "b": rb.Number(("f.json",), lambda s: rb.ints(rb.J(s, "f.json", "big")))}
    text = "x <!--n:r-->19/20, y <!--n:b-->2,722, z (f.json)\n"
    new, issues = rb.update_text(text, src, nums)
    assert issues == [] and new == text
    assert rb.check_text(new, src, nums)[0] == []


def test_report_build_update_refuses_reshaped_prose(tmp_path):
    """A value whose shape no longer matches is reported, never overwritten (prose stays untouched)."""
    src = _toy(tmp_path)
    text = "a <!--n:a.delta-->about five percent (res/a.json)\n"
    new, issues = rb.update_text(text, src, _numbers())
    assert new == text and [i.kind for i in issues] == ["unmatched"]


def test_report_build_placeholders_counted():
    text = "x [T1: C1 estimate] y [T2: C4] z [T2 unseal: per-lineage] w [not a placeholder]"
    assert len(rb.placeholders(text)) == 3


def test_report_build_runtime_table_sums_manifests(tmp_path):
    base = tmp_path / "results"
    (base / "w9" / "sub").mkdir(parents=True)
    (base / "w9" / "MANIFEST.tsv").write_text("script\twall_s\tcpu_s\textra\na\t10\t3600\t{}\nb\t20\t1800\t{}\n")
    # an early-format row with extra trailing fields must still parse
    (base / "w9" / "sub" / "PROVENANCE.tsv").write_text("script\twall_s\tcpu_s\nc\t5\t1800\textra\tfields\n")
    (base / "w9" / "notes.tsv").write_text("script\twall_s\tcpu_s\nz\t1\t99999\n")  # not a manifest
    t = rb.runtime_table(rb.Source(tmp_path, ref=None))
    assert t["w9"]["runs"] == 3
    assert t["w9"]["cpu_s"] == pytest.approx(7200) and t["w9"]["wall_s"] == pytest.approx(35)


def test_report_build_source_reads_committed_ref(tmp_path):
    """The default source is the committed tree: an uncommitted edit must not reach the report."""
    def git(*a):
        subprocess.run(["git", *a], cwd=tmp_path, check=True, capture_output=True)
    git("init", "-q")
    git("config", "user.email", "toy@example.org")
    git("config", "user.name", "toy")
    (tmp_path / "v.json").write_text('{"x": 1}')
    git("add", "v.json")
    git("commit", "-q", "-m", "toy")
    (tmp_path / "v.json").write_text('{"x": 2}')
    assert rb.Source(tmp_path, ref="HEAD").json("v.json")["x"] == 1
    assert rb.Source(tmp_path, ref=None).json("v.json")["x"] == 2
    with pytest.raises(FileNotFoundError):
        rb.Source(tmp_path, ref="HEAD").text("absent.json")


def test_report_build_registry_entries_are_well_formed():
    assert len(rb.NUMBERS) > 100
    for key, num in rb.NUMBERS.items():
        assert rb.MARKER.fullmatch(f"<!--n:{key}-->"), key
        assert num.sources and all(isinstance(p, str) and p for p in num.sources), key
        assert callable(num.fn), key


def test_report_build_registry_never_names_label_files():
    import re
    protected = re.compile(r"CRISPRGene|ScreenGeneEffect|CommonEssential|Nonessential|labels_test|labels_val|"
                           r"data/locked|data/depmap")
    for key, num in rb.NUMBERS.items():
        assert not any(protected.search(p) for p in num.sources), key


def test_report_build_cli_exit_codes(tmp_path):
    ok = tmp_path / "ok.md"
    ok.write_text("No markers here.\n")
    assert rb.main(["--report", str(ok), "--check", "--worktree"]) == 0
    bad = tmp_path / "bad.md"
    bad.write_text("x <!--n:not.a.key-->1\n")
    assert rb.main(["--report", str(bad), "--check", "--worktree"]) == 1
    ph = tmp_path / "ph.md"
    ph.write_text("C1: [T1: estimate]\n")
    assert rb.main(["--report", str(ph), "--check", "--worktree"]) == 0
    assert rb.main(["--report", str(ph), "--check", "--worktree", "--require-filled"]) == 1
