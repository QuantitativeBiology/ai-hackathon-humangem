"""Shared helpers for download_*.py scripts (streamed, resumable-by-skip, size-checked downloads)."""
import sys
import time
from pathlib import Path

import requests

DATA = Path(__file__).resolve().parents[1] / "data"
UA = {"User-Agent": "nextgem-data-fetch/1.0 (academic research; contact: project PI)"}
MAX_BYTES = 15 * 1024**3


def fetch(url: str, dest: Path, overwrite: bool = False, retries: int = 3, **kw) -> Path:
    """Stream `url` to `dest`. Skips if present with nonzero size. Refuses files > 15 GB."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        print(f"[skip] {dest.name} exists ({dest.stat().st_size:,} B)")
        return dest
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=(30, 300), headers=UA, allow_redirects=True, **kw) as r:
                r.raise_for_status()
                n = int(r.headers.get("content-length", 0))
                if n > MAX_BYTES:
                    sys.exit(f"{url} is {n/1e9:.1f} GB > 15 GB limit; ask the PI")
                ctype = r.headers.get("content-type", "")
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                if "text/html" in ctype and dest.suffix not in (".html", ".htm"):
                    tmp.unlink()
                    raise RuntimeError(f"{url} returned HTML (likely a login/bot-challenge page), not data")
                tmp.rename(dest)
                print(f"[ok] {dest.name} {dest.stat().st_size:,} B  <- {url}")
                return dest
        except requests.RequestException as e:
            print(f"[retry {attempt+1}] {url}: {e}")
            time.sleep(5 * (attempt + 1))
    sys.exit(f"failed: {url}")
