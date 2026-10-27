"""nextgem.run.provenance (DESIGN §11.4): one manifest row per successful run with full SHA-256 of inputs and
outputs, commit, dirty flag, versions, seeds, solver, wall/CPU time; no row if the body fails; a dirty tree
is refused unless explicitly allowed. Toy files only."""

from __future__ import annotations

import json
import re

import pytest

from nextgem import run


def _rows(path):
    lines = path.read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in lines[1:]]


def test_provenance_row_full_shas_and_versions(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "git_dirty", lambda *a, **k: False)
    inp = tmp_path / "in.txt"
    inp.write_text("input\n")
    out = tmp_path / "out.txt"
    man = tmp_path / "MANIFEST.tsv"
    with run.provenance(man, inputs=[inp], seeds=[11, 23], solver="gurobi", workers=3, author="test",
                        argv=["toy.py", "--x"]) as prov:
        out.write_text("output\n")
        prov.output(out)
        prov.note(look_id="L1")
    (row,) = _rows(man)
    assert list(row) == run.MANIFEST_COLUMNS
    ins, outs = json.loads(row["inputs"]), json.loads(row["outputs"])
    assert list(ins.values()) == [run.sha256_file(inp)] and list(outs.values()) == [run.sha256_file(out)]
    assert all(re.fullmatch(r"[0-9a-f]{64}", v) for v in [*ins.values(), *outs.values()])
    assert row["script"] == "toy.py" and row["argv"] == "--x" and row["dirty"] == "False"
    assert re.fullmatch(r"[0-9a-f]{40}", row["commit"])
    assert json.loads(row["seeds"]) == [11, 23] and row["workers"] == "3" and row["solver"] == "gurobi"
    v = json.loads(row["versions"])
    assert {"python", "cobra", "optlang", "gurobipy"} <= set(v)
    assert row["solver_version"] == (v["gurobipy"] or "")
    assert float(row["wall_s"]) >= 0 and float(row["cpu_s"]) >= 0
    assert json.loads(row["extra"]) == {"look_id": "L1"}


def test_no_row_when_body_fails(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "git_dirty", lambda *a, **k: False)
    man = tmp_path / "MANIFEST.tsv"
    with pytest.raises(ValueError):
        with run.provenance(man):
            raise ValueError("boom")
    assert not man.exists()


def test_dirty_tree_refused_unless_allowed(tmp_path, monkeypatch):
    monkeypatch.setattr(run, "git_dirty", lambda *a, **k: True)
    man = tmp_path / "MANIFEST.tsv"
    with pytest.raises(run.DirtyTreeError):
        with run.provenance(man):
            pass
    assert not man.exists()
    with run.provenance(man, allow_dirty=True):
        pass
    (row,) = _rows(man)
    assert row["dirty"] == "True"
