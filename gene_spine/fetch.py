"""Download every source with a url or api entry and record it in the manifest.

    python -m gene_spine.fetch            # fetch all
    python -m gene_spine.fetch clingen    # one source

Manual sources are listed at the end with where to get them. Nothing is
overwritten without being re-hashed; the manifest is the provenance record
for a release.
"""
from __future__ import annotations

import hashlib
import json
import sys
from datetime import date
from pathlib import Path

import requests
import yaml

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "sources"
MANIFEST = SRC / "manifest.json"
HEADERS = {"User-Agent": "gene-spine/0.1 (research; contact via repo)"}


def load_config() -> dict:
    return yaml.safe_load((ROOT / "config" / "sources.yaml").read_text())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def snapshot_hash(path: Path, cfg: dict) -> tuple[str, str, int]:
    """Hash of the payload that matters, the scope it covers, and its size.

    Some services wrap the data in an envelope that changes between identical responses.
    SysNDD reports its own `meta.executionTime`, so two fetches of the same 3,271 rows differ
    in bytes and release.json stops being comparable between runs. Where `hash_canonical` names
    a path, the hash covers that node serialised with sorted keys.

    The size travels with the hash for the same reason. Recording the file size next to a
    canonical hash would leave one volatile field behind: "0.63 secs" and "0.6 secs" are a byte
    apart, which is exactly what made release.json churn after the hash was already stable.
    """
    scope = cfg.get("hash_canonical")
    if not scope:
        return sha256(path), "file bytes", path.stat().st_size
    payload = json.loads(path.read_text(encoding="utf-8"))
    for key in str(scope).split("."):
        payload = payload[key]
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return (hashlib.sha256(blob).hexdigest(),
            f"canonicalised JSON at '{scope}', keys sorted", len(blob))


def record(manifest: dict, sid: str, path: Path, cfg: dict, how: str):
    digest, scope, size = snapshot_hash(path, cfg)
    entry = {
        # as_posix, not str: a manifest written on Windows would otherwise carry backslashes
        # and flip separators on every alternate local and CI build.
        "file": path.relative_to(ROOT).as_posix(),
        "fetched": date.today().isoformat(),
        "how": how,
        "url": cfg.get("url") or cfg.get("api") or cfg.get("landing"),
        "version": cfg.get("version"),
        "bytes": size,
        "sha256": digest,
        "sha256_scope": scope,
    }
    # Provenance a hand export carries that the file name cannot: what the site called it, when
    # the data was released, when it was pulled. Declared in config, so it travels into
    # release.json with the hash and a reader can tell which release a figure came from.
    for key in ("original_filename", "release_date", "exported", "unscored_note"):
        if cfg.get(key):
            entry[key] = cfg[key]
    manifest[sid] = entry


def fetch_quarterly(sid: str, cfg: dict, manifest: dict):
    """Try the newest quarterly snapshot first (1 Jan/Apr/Jul/Oct), walking back up to two years."""
    out = SRC / sid / cfg["file"]
    out.parent.mkdir(parents=True, exist_ok=True)
    today = date.today()
    quarters = []
    y, m = today.year, ((today.month - 1) // 3) * 3 + 1
    for _ in range(8):
        quarters.append(f"{y}-{m:02d}-01")
        m -= 3
        if m < 1:
            m += 12
            y -= 1
    last_err = None
    for q in quarters:
        url = cfg["quarterly_url"].format(date=q)
        try:
            r = requests.get(url, headers=HEADERS, timeout=180)
            if r.status_code == 200 and r.content[:7] == b"hgnc_id":
                out.write_bytes(r.content)
                record(manifest, sid, out, dict(cfg, url=url, version=q), "quarterly")
                print(f"  {sid}: snapshot {q}, {out.stat().st_size:,} bytes")
                return
            last_err = f"HTTP {r.status_code}, starts {r.content[:20]!r}"
        except Exception as e:
            last_err = str(e)
    raise RuntimeError(f"{sid}: no quarterly snapshot found; last error {last_err}")


def fetch_url(sid: str, cfg: dict, manifest: dict):
    out = SRC / sid / cfg["file"]
    out.parent.mkdir(parents=True, exist_ok=True)
    # Streamed, with a per-source timeout: mondo.obo is 53 MB and does not finish inside the
    # 120 seconds that suits the small tables.
    r = requests.get(cfg["url"], headers=HEADERS, timeout=cfg.get("timeout", 120), stream=True)
    r.raise_for_status()
    with out.open("wb") as fh:
        for chunk in r.iter_content(1 << 20):
            fh.write(chunk)
    record(manifest, sid, out, cfg, "url")
    print(f"  {sid}: {out.stat().st_size:,} bytes")


def fetch_panelapp(sid: str, cfg: dict, manifest: dict):
    """Use panel_id when given, otherwise find the panel by name; then pull all genes."""
    base = cfg["api"].rstrip("/") + "/"
    want = cfg.get("panel_name_match", "").lower()
    panel = None
    if cfg.get("panel_id"):
        r = requests.get(f"{base}{cfg['panel_id']}/", headers=HEADERS, timeout=120)
        r.raise_for_status()
        panel = r.json()
    url = base
    while url and panel is None:
        r = requests.get(url, headers=HEADERS, timeout=120)
        r.raise_for_status()
        data = r.json()
        for p in data.get("results", []):
            if want in p.get("name", "").lower():
                panel = p
                break
        url = data.get("next")
    if panel is None:
        raise RuntimeError(f"{sid}: no panel matching '{cfg['panel_name_match']}' at {base}")
    genes, url = [], f"{base}{panel['id']}/genes/?page_size=500"
    while url:
        r = requests.get(url, headers=HEADERS, timeout=120)
        r.raise_for_status()
        data = r.json()
        genes.extend(data.get("results", []))
        url = data.get("next")
    out = SRC / sid / cfg["file"]
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"id": panel["id"], "name": panel["name"],
                               "version": panel.get("version"), "genes": genes}, indent=1),
                   encoding="utf-8", newline="\n")
    cfg = dict(cfg, version=str(panel.get("version")))
    record(manifest, sid, out, cfg, f"panelapp:{panel['id']}")
    print(f"  {sid}: panel {panel['id']} '{panel['name']}' v{panel.get('version')}, {len(genes)} genes")


def main(argv: list[str]):
    cfg_all = load_config()
    SRC.mkdir(exist_ok=True)
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    wanted = set(argv) or set(cfg_all)
    manual, failed = [], []
    for sid, cfg in cfg_all.items():
        if sid not in wanted:
            continue
        if cfg.get("manual"):
            existing = next((c for c in (SRC / sid / cfg["file"], SRC / cfg["file"], ROOT / cfg["file"]) if c.exists()), SRC / sid / cfg["file"])
            if existing.exists():
                record(manifest, sid, existing, cfg, "manual")
            else:
                manual.append((sid, cfg))
            continue
        try:
            if "url" in cfg and "quarterly_url" in cfg:
                try:
                    fetch_url(sid, cfg, manifest)
                    ok = (SRC / sid / cfg["file"]).read_bytes()[:7] == b"hgnc_id"
                except Exception as e:
                    print(f"  {sid}: current-release URL failed ({e}); trying quarterly archive")
                    ok = False
                if not ok:
                    fetch_quarterly(sid, cfg, manifest)
            elif "quarterly_url" in cfg:
                fetch_quarterly(sid, cfg, manifest)
            elif "api" in cfg:
                fetch_panelapp(sid, cfg, manifest)
            elif "url" in cfg:
                fetch_url(sid, cfg, manifest)
        except Exception as e:  # keep going; a release can still be built from what fetched
            print(f"  {sid}: FAILED {e}")
            failed.append(sid)
    MANIFEST.write_text(json.dumps(manifest, indent=1), encoding="utf-8", newline="\n")
    if failed:
        print("\nFailed to fetch (build will skip these):", ", ".join(failed))
    if manual:
        print("\nManual sources still missing. Place each file at sources/<id>/<file>:")
        for sid, cfg in manual:
            print(f"  {sid:28s} {cfg['file']:34s} {cfg.get('landing','')}")
            if cfg.get("note"):
                print(f"  {'':28s} {cfg['note']}")


if __name__ == "__main__":
    main(sys.argv[1:])
