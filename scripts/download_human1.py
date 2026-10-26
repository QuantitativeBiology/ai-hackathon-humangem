"""Download Human1 = the last Human-GEM v1.x release (v1.19.0) into data/Human1/.

Why v1.19.0: the GitHub releases/tags of SysBioChalmers/Human-GEM (checked 2026-10-26) go
... v1.18.0 (2023-12-08), v1.19.0 (2024-06-24), then v2.0.0 (2026-03-30), v2.0.1, v2.1.0, v2.1.1.
There is no v1.19.x patch release, so v1.19.0 is the latest 1.x. The script checks this again
against the GitHub API (warns if a newer 1.x tag appears) and pins the files to the tag's commit.

Every file is checked against the pinned upstream SHA-256 in ``FILES`` (downloads before they are
written, existing files before they are kept); any mismatch raises and nothing unverified is written.

Layout mirrors the Human-GEM repo so the same relative paths work for both versions:
  data/Human1/model/{Human-GEM.xml, genes.tsv, reactions.tsv, metabolites.tsv, README.md}
  data/Human1/data/metabolicTasks/metabolicTasks_{Essential,Full}.txt
  data/Human1/LICENSE.md, data/Human1/SOURCE.txt (tag, commit, URLs, date)

Usage: .venv/bin/python scripts/download_human1.py [--force]
"""

import argparse
import hashlib
import json
import re
import urllib.request
from datetime import date
from pathlib import Path

REPO = "SysBioChalmers/Human-GEM"
TAG = "v1.19.0"
COMMIT = "022ed5c2c971e27f739e149c4bbe1328de38337f"  # refs/tags/v1.19.0
# path in repo -> SHA-256 of the file at COMMIT (raw.githubusercontent.com, checked 2026-10-26;
# the SBML hash was independently re-derived by the baseline reviewer from both tag and commit URLs)
FILES = {
    "model/Human-GEM.xml": "5aa161fbd52e79051dc445d87965db0716219d95ef1a6f0a8c5bd51fd60e075c",
    "model/genes.tsv": "690b580ba322b423ef24b6f69e9211ee1c215e4defc48179b67378d7ab20646a",
    "model/reactions.tsv": "b45c542eaeb0ed52d31eb05858c005f1fa5bf7a90807d8811497a34fddd836d0",
    "model/metabolites.tsv": "7b971bdbc4ccc71a34e500f26ddf06f6be1e21c30f3d50dbe512b0389ca52ed1",
    "model/README.md": "8bfa0bfd39ce7783df8b7191184695a8c6d38a928c71fb6181070d22007ff9f8",
    "data/metabolicTasks/metabolicTasks_Essential.txt": "d0c58503fc44d1eb49ea1b89f4cc90930a06bcd2c3972083c0986cfb126a0f18",
    "data/metabolicTasks/metabolicTasks_Full.txt": "4896f6505d284725a30aadeea6548f3a6653f32a9c387d219860033a735ecc42",
    "LICENSE.md": "bb5e0179a1e9bdb55634a4303784cf73c297a1dcee27b92cd546bf398075531c",
}
OUT = Path(__file__).resolve().parents[1] / "data" / "Human1"


def github_json(url: str):
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def check_tag() -> None:
    """Confirm TAG points at COMMIT and that it is still the highest v1.x tag."""
    ref = github_json(f"https://api.github.com/repos/{REPO}/git/refs/tags/{TAG}")
    if ref["object"]["sha"] != COMMIT:
        raise RuntimeError(f"{TAG} now points to {ref['object']['sha']}, expected {COMMIT}")
    tags = [t["name"] for t in github_json(f"https://api.github.com/repos/{REPO}/tags?per_page=100")]
    v1 = sorted(
        (tuple(map(int, m.groups())), t) for t in tags if (m := re.fullmatch(r"v1\.(\d+)\.(\d+)", t))
    )
    if v1[-1][1] != TAG:
        print(f"WARNING: newer 1.x tag exists: {v1[-1][1]} (pinned {TAG})")
    print(f"{TAG} -> {COMMIT}; highest v1.x tag: {v1[-1][1]}")


def verify(rel: str, content: bytes) -> None:
    digest = hashlib.sha256(content).hexdigest()
    if digest != FILES[rel]:
        raise RuntimeError(f"SHA-256 mismatch for {rel}: got {digest}, expected {FILES[rel]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    args = parser.parse_args()

    check_tag()
    for rel in FILES:
        dest = OUT / rel
        if dest.exists() and not args.force:
            verify(rel, dest.read_bytes())
            print(f"exists, SHA-256 ok: {dest.relative_to(OUT.parent)}")
            continue
        url = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/{rel}"
        with urllib.request.urlopen(url, timeout=300) as r:
            content = r.read()
        verify(rel, content)
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(content)
        print(f"{dest.relative_to(OUT.parent)}: {len(content):,} bytes, SHA-256 ok")
    (OUT / "SOURCE.txt").write_text(
        f"repo\thttps://github.com/{REPO}\ntag\t{TAG}\ncommit\t{COMMIT}\n"
        f"url_template\thttps://raw.githubusercontent.com/{REPO}/{COMMIT}/<path>\n"
        f"downloaded\t{date.today().isoformat()}\n"
    )


if __name__ == "__main__":
    main()
