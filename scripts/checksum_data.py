"""Write SHA-256 checksums of every file under data/ (except data/README.md) to data/checksums.tsv.

Usage: .venv/bin/python scripts/checksum_data.py [subdir ...]
"""

import hashlib
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
OUT = DATA / "checksums.tsv"


def sha256(path: Path) -> tuple[str, str]:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return str(path.relative_to(DATA)), h.hexdigest()


def main():
    roots = [DATA / s for s in sys.argv[1:]] or [DATA]
    # data/locked and data/labels hold the built label files: their hashes live in FROZEN.tsv, and the
    # tripwire (DESIGN.md §11.1) denies this script from opening them.
    skip = (DATA / "locked", DATA / "labels")
    files = sorted(p for r in roots for p in r.rglob("*") if p.is_file() and p not in (OUT, DATA / "README.md")
                   and not any(s in p.parents for s in skip))
    # Keep entries for files outside the requested subdirs; path -> digest
    existing = {}
    if OUT.exists():
        rows = (line.split("\t") for line in OUT.read_text().splitlines()[1:])
        existing = {path: digest for digest, path in rows}
    with ProcessPoolExecutor(8) as ex:
        existing.update({rel: digest for rel, digest in ex.map(sha256, files)})
    OUT.write_text("sha256\tpath\n" + "".join(f"{d}\t{p}\n" for p, d in sorted(existing.items())))
    print(f"{len(files)} files hashed -> {OUT}")


if __name__ == "__main__":
    main()
