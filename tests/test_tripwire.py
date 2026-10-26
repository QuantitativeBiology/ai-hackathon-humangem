"""§12 "Agent guards" and "Labels" negative controls: the tripwire (DESIGN.md §11.1)."""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import toy

RAW = toy.RAW_GE


def tw():
    return sys.modules["nextgem_tripwire"]


@pytest.fixture()
def raw_file(tmp_path):
    """A toy raw-named label file outside the data roots (written by the test: allowed)."""
    d = tmp_path / "depmap"
    d.mkdir()
    f = d / RAW
    pd.DataFrame({"ModelID": ["ACH-000001"], "G1 (1)": [-0.1]}).to_csv(f, index=False)
    return f


def test_sitecustomize_installed_and_current():
    r = subprocess.run([toy.PY, str(toy.SCRIPTS / "install_tripwire.py"), "--check"], capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert tw().installed()


def test_tripwire_active_without_importing_nextgem():
    code = "import sys; print('nextgem' in sys.modules, sys.modules['nextgem_tripwire'].installed())"
    r = subprocess.run([toy.PY, "-c", code], capture_output=True, text=True)
    assert r.stdout.split() == ["False", "True"]


def test_open_and_pandas_read_of_raw_file_denied(raw_file):
    with pytest.raises(PermissionError):
        open(raw_file)
    with pytest.raises(PermissionError):
        pd.read_csv(raw_file)
    with pytest.raises(PermissionError):
        os.open(raw_file, os.O_RDONLY)


def test_pyarrow_readers_denied(tmp_path):
    import pyarrow as pa
    import pyarrow.csv as pacsv
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq

    f = tmp_path / "labels_train.parquet"
    # writing a tier-named file from a non-entry point is itself denied
    with pytest.raises(PermissionError):
        pq.write_table(pa.table({"a": [1]}), f)
    with pytest.raises(PermissionError):
        pd.DataFrame({"a": [1]}).to_parquet(f)
    # create it under a neutral name, then rename -> the rename is denied too
    g = tmp_path / "neutral.parquet"
    pq.write_table(pa.table({"a": [1]}), g)
    with pytest.raises(PermissionError):
        os.rename(g, f)
    # a raw-named csv (toy, outside the data roots) can be written; every reader is denied
    raw = tmp_path / RAW
    raw.write_text("ModelID,G1 (1)\nACH-000001,-0.2\n")
    for fn in (lambda: pq.read_table(raw), lambda: pq.ParquetFile(raw), lambda: ds.dataset(raw),
               lambda: pacsv.read_csv(raw), lambda: pa.memory_map(str(raw)), lambda: pd.read_parquet(raw)):
        with pytest.raises(PermissionError):
            fn()


def test_symlink_to_label_file_denied(raw_file, tmp_path):
    link = tmp_path / "innocuous.csv"
    os.symlink(raw_file, link)  # toy target outside data roots: link creation allowed
    with pytest.raises(PermissionError):
        open(link)


def test_links_to_real_raw_files_denied(tmp_path):
    from nextgem import paths

    real = paths.PROJECT_ROOT / "data" / "depmap" / RAW
    with pytest.raises(PermissionError):
        os.link(real, tmp_path / "x.csv")
    with pytest.raises(PermissionError):
        os.symlink(real, tmp_path / "y.csv")


def test_locked_dir_denied_by_path():
    from nextgem import paths

    with pytest.raises(PermissionError):
        open(paths.PROJECT_ROOT / "data" / "locked" / "HOOK_PROBE")
    assert tw().classify("/anything/data/locked/x.bin") == "locked"
    assert tw().classify("/anything/data/locked/sealed/T1/a.tsv") == "sealed"
    assert tw().classify("/anything/data/locked/labels_test.parquet") == "test"
    assert tw().classify("/x/Model.csv") is None


def test_process_arguments_naming_label_files_denied():
    with pytest.raises(PermissionError):
        subprocess.run(["echo", "data/depmap/" + RAW])
    with pytest.raises(PermissionError):
        os.system("echo labels_test")
    assert subprocess.run(["echo", "fine"], capture_output=True, text=True).stdout.strip() == "fine"


def test_private_labels_reader_from_non_entry_script_denied(built_toy, tmp_path):
    """§12: tripwire raises on a private labels.py reader called from a non-entry-point script."""
    script = tmp_path / "sneaky.py"
    script.write_text(textwrap.dedent("""
        from nextgem import labels
        try:
            labels._read_raw("gene_effect", ["G1 (1001)"])
            print("READ")
        except PermissionError:
            print("DENIED")
    """))
    r = subprocess.run([toy.PY, str(script)], capture_output=True, text=True, env={**os.environ, **built_toy.env})
    assert r.stdout.strip() == "DENIED", r.stderr[-1000:]


def test_direct_pyarrow_read_of_training_labels_denied_in_subprocess(built_toy):
    """§12: tripwire raises on pyarrow.parquet.read_table of a label file (not through labels.py)."""
    # (the path goes through the environment: naming it in argv would itself be denied by the Popen check)
    code = ("import pyarrow.parquet as pq, os\n"
            "try:\n    pq.read_table(os.environ['TARGET']); print('READ')\n"
            "except PermissionError:\n    print('DENIED')\n")
    target = str(built_toy.data / "labels" / "labels_train.parquet")
    r = subprocess.run([toy.PY, "-c", code], capture_output=True, text=True,
                       env={**os.environ, **built_toy.env, "TARGET": target})
    assert r.stdout.strip() == "DENIED"


def test_fold_api_allowed_from_any_entry_point(built_toy, monkeypatch):
    from nextgem import labels

    monkeypatch.setenv("NEXTGEM_DATA_ROOT", str(built_toy.data))
    monkeypatch.setenv("NEXTGEM_STATE_ROOT", str(built_toy.state))
    df = labels.train()
    assert set(df.ModelID) == set(built_toy.ids_in("train"))


def test_entry_point_policy_table():
    """The policy table alone: evaluate_val.py reads validation only; evaluate_test.py validation and test;
    check_frozen.py/verify_ledgers.py hash only; build_labels.py everything; other scripts nothing."""
    pol = tw().policy
    assert pol("val", False, "evaluate_val.py", {"val"})[0] and not pol("test", False, "evaluate_val.py", {"test"})[0]
    assert pol("test", False, "evaluate_test.py", {"test"})[0] and pol("val", False, "evaluate_test.py", {"val"})[0]
    assert pol("sealed", True, "evaluate_test.py", {"sealed"})[0] and not pol("test", True, "evaluate_test.py", {"test"})[0]
    assert pol("test", False, "check_frozen.py", {"hash"})[0] and not pol("test", False, "check_frozen.py", {"test"})[0]
    assert pol("test", False, "check_frozen.py", {"hash"})[2] is True  # hashing needs a frozen entry point
    assert pol("raw", False, "build_labels.py", {"raw"})[0] and not pol("raw", False, "sneaky.py", {"raw"})[0]
    assert pol("train", False, None, {"train"}) == (True, "fold-scoped training API", False)


_ENTRY_CODE = textwrap.dedent("""
    import sys
    sys.modules['__main__'].__file__ = sys.argv[1]
    from nextgem import labels
    out = []
    for name, fn in [('val', labels.validation), ('test', labels.test)]:
        try:
            fn(); out.append(name + ':ok')
        except PermissionError:
            out.append(name + ':denied')
    try:
        labels.sha256_protected(labels.tier_file('test')); out.append('hash:ok')
    except PermissionError:
        out.append('hash:denied')
    print(' '.join(out))
""")


def test_entry_points_must_be_frozen(built_toy):
    """M4: an entry point is trusted only if its SHA-256 is in FROZEN.tsv. A spoofed name of a script that
    does not exist (not frozen) gets nothing; the frozen check_frozen.py may hash."""
    res = {}
    for entry in ("evaluate_val.py", "evaluate_test.py", "check_frozen.py", "not_an_entry.py"):
        r = subprocess.run([toy.PY, "-c", _ENTRY_CODE, str(toy.SCRIPTS / entry)], capture_output=True, text=True,
                           env={**os.environ, **built_toy.env})
        res[entry] = r.stdout.strip()
    assert res["evaluate_val.py"] == "val:denied test:denied hash:denied"
    assert res["evaluate_test.py"] == "val:denied test:denied hash:denied"
    assert res["check_frozen.py"] == "val:denied test:denied hash:ok"
    assert res["not_an_entry.py"] == "val:denied test:denied hash:denied"


def test_copied_checkout_with_fake_entry_point_denied(built_toy, tmp_path):
    """M4 (reviewer's bypass): copy src/ elsewhere and add scripts/evaluate_test.py next to it."""
    import shutil

    co = tmp_path / "co"
    shutil.copytree(toy.PROJECT / "src", co / "src", ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
    (co / "scripts").mkdir()
    fake = co / "scripts" / "evaluate_test.py"
    fake.write_text("from nextgem import labels\ntry:\n    labels.test(); print('READ')\n"
                    "except PermissionError:\n    print('DENIED')\n")
    r = subprocess.run([toy.PY, str(fake)], capture_output=True, text=True,
                       env={**os.environ, **built_toy.env, "PYTHONPATH": str(co / "src")})
    assert r.stdout.strip() == "DENIED", r.stderr[-800:]
    # and the copied labels.py really was the one on the stack
    r = subprocess.run([toy.PY, "-c", "import nextgem.labels as l; print(l.__file__)"], capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": str(co / "src")})
    assert r.stdout.strip().startswith(str(co))


def test_glob_arguments_denied():
    with pytest.raises(PermissionError):
        subprocess.run(["ls", "/nonexistent/lock*/x"])
    with pytest.raises(PermissionError):
        subprocess.run(["ls", "/nonexistent/dep?ap/*.csv"])


def test_access_log_records_denials(raw_file):
    from nextgem import paths

    with pytest.raises(PermissionError):
        open(raw_file)
    log = (paths.PROJECT_ROOT / "logs" / "label_access.log").read_text().splitlines()
    assert any("DENY" in ln and str(raw_file) in ln for ln in log[-50:])


def test_fail_closed_fallback(tmp_path, raw_file):
    """If tripwire.py cannot be loaded, sitecustomize installs a hook that denies every label open."""
    sys.path.insert(0, str(toy.SCRIPTS))
    try:
        import install_tripwire
    finally:
        sys.path.pop(0)
    site = tmp_path / "site"
    site.mkdir()
    (site / "sitecustomize.py").write_text(install_tripwire.TEMPLATE.format(mark=install_tripwire.MARK, path="/nonexistent/tripwire.py"))
    code = f"import sys, os; sys.path.insert(0, {str(site)!r}); import sitecustomize\n" \
           "try:\n    open(os.environ['TARGET']); print('READ')\nexcept PermissionError:\n    print('DENIED')\n"
    r = subprocess.run([toy.PY, "-S", "-c", code], capture_output=True, text=True, env={**os.environ, "TARGET": str(raw_file)})
    assert r.stdout.strip() == "DENIED"
    assert "fail-closed" in r.stderr


def test_unprotected_io_unaffected(tmp_path):
    f = tmp_path / "expression.parquet"
    pd.DataFrame({"x": np.arange(3)}).to_parquet(f)
    assert pd.read_parquet(f).x.sum() == 3
