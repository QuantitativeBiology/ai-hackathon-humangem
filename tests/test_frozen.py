"""§11.3 FROZEN.tsv and test-run checks, in a toy git repository."""

from __future__ import annotations

import re

import pytest

import toy


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    r = tmp_path / "state"
    toy.init_repo(r)
    (r / "src").mkdir()
    (r / "src" / "a.py").write_text("a = 1\n")
    (r / "DESIGN.md").write_text(toy.DESIGN_TOY)
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(r))
    from nextgem import frozen

    frozen.upsert_rows([frozen.file_row("src/a.py"), {"item": "label.primary", "kind": "value", "sha256_or_value": "L_prob"}])
    toy.commit_all(r, "freeze")
    return r


def test_clean(repo):
    from nextgem import frozen

    assert frozen.check() == []
    assert frozen.check(staged=True) == []


def test_modified_frozen_file(repo):
    from nextgem import frozen

    (repo / "src" / "a.py").write_text("a = 2\n")
    assert any("differs from its frozen" in p for p in frozen.check())
    # staged mode looks at the index, not the working tree
    assert frozen.check(staged=True) == []
    toy.git(repo, "add", "src/a.py")
    assert any("differs from its frozen" in p for p in frozen.check(staged=True))


def test_git_ignored_frozen_file_checked_in_working_tree(repo):
    """A5: an ignored input is never in the index, so staged mode checks it in the working tree."""
    from nextgem import frozen

    (repo / ".gitignore").write_text("tmp/\n")
    (repo / "tmp").mkdir()
    (repo / "tmp" / "x.parquet").write_bytes(b"x1")
    frozen.upsert_rows([frozen.file_row("tmp/x.parquet")])
    toy.commit_all(repo, "freeze ignored input")
    assert frozen.check(staged=True) == []
    (repo / "tmp" / "x.parquet").write_bytes(b"x2")
    assert any("differs from its frozen" in p for p in frozen.check(staged=True))
    (repo / "tmp" / "x.parquet").unlink()
    assert any("is missing" in p for p in frozen.check(staged=True))


def test_changed_row_needs_logged_amendment(repo):
    from nextgem import frozen

    (repo / "src" / "a.py").write_text("a = 2\n")
    frozen.upsert_rows([frozen.file_row("src/a.py")])
    assert any("without an amendment" in p for p in frozen.check())
    frozen.upsert_rows([{**frozen.file_row("src/a.py"), "amendment": "A9"}])
    assert any("not in the DESIGN.md" in p for p in frozen.check())
    frozen.upsert_rows([{**frozen.file_row("src/a.py"), "amendment": "A1"}])
    assert frozen.check() == []


def test_removed_row_and_changed_value(repo):
    from nextgem import frozen

    rows = [r for r in frozen.read_rows() if r["item"] != "label.primary"]
    frozen.write_rows(rows)
    assert any("removed" in p for p in frozen.check())
    frozen.upsert_rows([{"item": "label.primary", "kind": "value", "sha256_or_value": "L_ge05"}])
    assert any("without an amendment" in p for p in frozen.check())


def test_new_rows_are_free(repo):
    from nextgem import frozen

    (repo / "src" / "b.py").write_text("b = 1\n")
    frozen.upsert_rows([frozen.file_row("src/b.py")])
    assert frozen.check() == []


def test_amendment_ids_only_from_table_data_rows(repo):
    """M5: header words, the separator and the placeholder row are not amendment ids."""
    from nextgem import frozen

    ids = frozen.amendment_ids(toy.DESIGN_TOY)
    assert ids == {"A1"}
    assert not ({"Change", "Date", "Section", "Reason", "analysis", "the"} & ids)
    (repo / "src" / "a.py").write_text("a = 9\n")
    frozen.upsert_rows([{**frozen.file_row("src/a.py"), "amendment": "Change"}])
    assert any("not in the DESIGN.md" in p for p in frozen.check())
    real = (toy.PROJECT / "DESIGN.md").read_text()
    real_ids = frozen.amendment_ids(real)  # only explicit A<n> ids from data rows of the real §17 table
    assert all(re.fullmatch(r"A\d+", i) for i in real_ids) and not ({"Change", "Date", "Section"} & real_ids)


def test_test_runs_union_baseline_and_retirement(repo):
    from nextgem import frozen

    frozen.record_test_run({"t::a": "passed", "t::b": "passed"}, True, "c1")
    frozen.record_test_run({"t::a": "passed", "t::b": "passed"}, True, "c2")
    assert frozen.check_test_runs() == []
    frozen.record_test_run({"t::a": "passed"}, False, "c3")  # partial runs do not count
    assert frozen.check_test_runs() == []
    frozen.record_test_run({"t::a": "xfailed", "t::c": "passed"}, True, "c4")
    probs = frozen.check_test_runs()
    assert any("t::a" in p and "xfailed" in p for p in probs) and any("t::b" in p and "missing" in p for p in probs)
    # not forgotten after one more full run (review S1)
    frozen.record_test_run({"t::a": "passed", "t::c": "passed"}, True, "c5")
    assert any("t::b" in p for p in frozen.check_test_runs())
    # ... unless retired by a logged amendment
    (repo / "results" / "test_retirements.tsv").write_text("nodeid\tamendment\nt::b\tA1\n")
    assert frozen.check_test_runs() == []
    (repo / "results" / "test_retirements.tsv").write_text("nodeid\tamendment\nt::b\tA7\n")
    assert frozen.check_test_runs()
    assert frozen.check(test_runs=False) == []  # record-only commits skip the test-run check
