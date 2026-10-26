"""§11.3 ledgers: hash chain, commit on append (rollback if the commit fails), additions only in git
history, look caps in code (40 validation, <= 20 W7, raised only by a committed authorisation file),
idempotent test looks, split rows never twice, unseal only after T2. All in toy git repositories."""

from __future__ import annotations

import pytest

import toy


@pytest.fixture()
def repo(tmp_path, monkeypatch):
    r = tmp_path / "state"
    toy.init_repo(r)
    (r / "README").write_text("toy\n")
    toy.commit_all(r, "init")
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(r))
    return r


SHA = "a" * 64


def test_append_chain_and_commit(repo):
    from nextgem import ledger

    n0 = int(toy.git(repo, "rev-list", "--count", "HEAD").stdout)
    ledger.append("validation", "non_look", look_id="kappa_rule", detail='{"x": 1}')
    ledger.append("validation", "started", look_id="v1", workstream="W2", config_sha="c1")
    ledger.append("validation", "output", look_id="v1", output_sha=SHA)
    ledger.append("validation", "completed", look_id="v1")
    rows = ledger.read("validation")
    assert [r["seq"] for r in rows] == ["1", "2", "3", "4"]
    assert rows[0]["prev_sha"] == ledger.genesis("validation") and rows[1]["prev_sha"] == rows[0]["row_sha"]
    assert ledger.verify("validation") == []
    assert int(toy.git(repo, "rev-list", "--count", "HEAD").stdout) == n0 + 4  # one commit per append
    assert toy.git(repo, "status", "--porcelain").stdout.strip() == ""


def test_commit_only_leaves_other_staged_changes(repo):
    from nextgem import ledger

    (repo / "pi_work.txt").write_text("staged by the PI\n")
    toy.git(repo, "add", "pi_work.txt")
    ledger.append("test", "non_look", look_id="note")
    assert toy.git(repo, "show", "--name-only", "--format=", "HEAD").stdout.split() == ["results/test_ledger.tsv"]
    assert "A  pi_work.txt" in toy.git(repo, "status", "--porcelain").stdout


def test_tampering_detected_and_append_refused(repo):
    from nextgem import ledger

    for i in range(3):
        ledger.append("validation", "non_look", look_id=f"n{i}")
    p = repo / "results" / "validation_ledger.tsv"
    p.write_text(p.read_text().replace("\tn1\t", "\tn9\t"))
    assert any("row_sha" in x for x in ledger.verify("validation"))
    with pytest.raises(ledger.LedgerError):
        ledger.append("validation", "non_look", look_id="n3")


def test_history_rewrite_detected(repo):
    from nextgem import ledger

    for i in range(3):
        ledger.append("validation", "non_look", look_id=f"n{i}")
    rows = ledger.read("validation")
    # rewrite history with a *consistent* chain: drop row 2, re-chain row 3, commit
    keep = [rows[0], dict(rows[2])]
    keep[1]["seq"], keep[1]["prev_sha"] = "2", keep[0]["row_sha"]
    keep[1]["row_sha"] = ledger.row_hash(keep[1])
    text = "\t".join(ledger.COLUMNS) + "\n" + "".join("\t".join(r[c] for c in ledger.COLUMNS) + "\n" for r in keep)
    (repo / "results" / "validation_ledger.tsv").write_text(text)
    toy.commit_all(repo, "rewrite")
    assert ledger.verify_chain("validation", ledger.read("validation")) == []  # chain alone looks fine
    assert any("additions only" in x for x in ledger.verify_history("validation"))


def test_validation_cap_and_authorisation(repo):
    from nextgem import ledger

    for i in range(40):
        ledger.append("validation", "started", look_id=f"v{i}", workstream="W2", config_sha="c", commit=False)
    with pytest.raises(ledger.LedgerError, match="cap reached"):
        ledger.append("validation", "started", look_id="v40", workstream="W2", config_sha="c", commit=False)
    auth = repo / "results" / "authorisations"
    auth.mkdir(parents=True)
    (auth / "val_cap_raise_42.md").write_text("PI authorisation\n")
    with pytest.raises(ledger.LedgerError):  # not committed yet
        ledger.append("validation", "started", look_id="v40", workstream="W2", config_sha="c", commit=False)
    toy.git(repo, "add", "results/authorisations/val_cap_raise_42.md")
    toy.git(repo, "commit", "-q", "-m", "auth")
    ledger.append("validation", "started", look_id="v40", workstream="W2", config_sha="c", commit=False)
    ledger.append("validation", "started", look_id="v41", workstream="W2", config_sha="c", commit=False)
    with pytest.raises(ledger.LedgerError, match="cap reached"):
        ledger.append("validation", "started", look_id="v42", workstream="W2", config_sha="c", commit=False)


def test_w7_subcap(repo):
    from nextgem import ledger

    for i in range(20):
        ledger.append("validation", "started", look_id=f"w{i}", workstream="W7", config_sha="c", commit=False)
    with pytest.raises(ledger.LedgerError, match="W7"):
        ledger.append("validation", "started", look_id="w20", workstream="W7", config_sha="c", commit=False)
    ledger.append("validation", "started", look_id="x0", workstream="W3", config_sha="c", commit=False)


def test_test_looks_idempotent(repo):
    from nextgem import ledger

    with ledger.look("test", "T1", workstream="M6", config_sha="cfg1") as lk:
        lk.output(SHA)
    with pytest.raises(ledger.LedgerError, match="already exists"):
        ledger.append("test", "started", look_id="T1", config_sha="cfg1")
    # crashed look: started without completed
    with pytest.raises(RuntimeError):
        with ledger.look("test", "T2", workstream="M7", config_sha="cfg2"):
            raise RuntimeError("crash")
    n = len(ledger.read("test"))
    resumed = ledger.append("test", "started", look_id="T2", config_sha="cfg2")
    assert len(ledger.read("test")) == n and resumed["look_id"] == "T2"  # still one look
    with pytest.raises(ledger.LedgerError, match="same config"):
        ledger.append("test", "started", look_id="T2", config_sha="other")
    with pytest.raises(ledger.LedgerError):
        ledger.append("test", "output", look_id="T3", output_sha=SHA)  # never started
    with pytest.raises(ledger.LedgerError):
        ledger.append("test", "output", look_id="T2", output_sha="not-a-sha")


def test_split_rows_once_and_unseal_after_t2(repo):
    from nextgem import ledger

    ledger.append("test", "split_started", look_id="split")
    ledger.append("test", "split_completed", look_id="split", output_sha=SHA)
    with pytest.raises(ledger.LedgerError, match="never redrawn"):
        ledger.append("test", "split_started", look_id="split")
    with pytest.raises(ledger.LedgerError):
        ledger.append("validation", "split_started", look_id="split")
    with pytest.raises(ledger.LedgerError, match="T2"):
        ledger.append("test", "unseal")
    with ledger.look("test", "T2", workstream="M7", config_sha="c"):
        pass
    ledger.append("test", "unseal")
    assert ledger.verify("test") == []


def test_commit_failure_rolls_back(repo):
    from nextgem import ledger

    ledger.append("validation", "non_look", look_id="ok")
    before = (repo / "results" / "validation_ledger.tsv").read_text()
    hook = repo / ".git" / "hooks" / "pre-commit"
    hook.write_text("#!/bin/sh\nexit 1\n")
    hook.chmod(0o755)
    with pytest.raises(ledger.LedgerError, match="rolled back"):
        ledger.append("validation", "started", look_id="v1", workstream="W2", config_sha="c")
    assert (repo / "results" / "validation_ledger.tsv").read_text() == before
    # a look whose started row cannot be committed never runs its body
    ran = []
    with pytest.raises(ledger.LedgerError):
        with ledger.look("validation", "v2", workstream="W2", config_sha="c"):
            ran.append(1)
    assert ran == []
