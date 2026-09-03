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


def record(manifest: dict, sid: str, path: Path, cfg: dict, how: str):
    manifest[sid] = {
        "file": str(path.relative_to(ROOT)),
        "fetched": date.today().isoformat(),
        "how": how,
        "url": cfg.get("url") or cfg.get("api") or cfg.get("landing"),
        "version": cfg.get("version"),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def fetch_url(sid: str, cfg: dict, manifest: dict):
    out = SRC / sid / cfg["file"]
    out.parent.mkdir(parents=True, exist_ok=True)
    r = requests.get(cfg["url"], headers=HEADERS, timeout=120)
    r.raise_for_status()
    out.write_bytes(r.content)
    record(manifest, sid, out, cfg, "url")
    print(f"  {sid}: {out.stat().st_size:,} bytes")


def fetch_panelapp(sid: str, cfg: dict, manifest: dict):
    """Find the panel by name (ids drift between deployments), then pull all genes."""
    base = cfg["api"].rstrip("/") + "/"
    want = cfg["panel_name_match"].lower()
    url = base
    panel = None
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
        raise SystemExit(f"{sid}: no panel matching '{cfg['panel_name_match']}' at {base}")
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
                               "version": panel.get("version"), "genes": genes}, indent=1))
    cfg = dict(cfg, version=str(panel.get("version")))
    record(manifest, sid, out, cfg, f"panelapp:{panel['id']}")
    print(f"  {sid}: panel {panel['id']} '{panel['name']}' v{panel.get('version')}, {len(genes)} genes")


def main(argv: list[str]):
    cfg_all = load_config()
    SRC.mkdir(exist_ok=True)
    manifest = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}
    wanted = set(argv) or set(cfg_all)
    manual = []
    for sid, cfg in cfg_all.items():
        if sid not in wanted:
            continue
        if cfg.get("manual"):
            existing = SRC / sid / cfg["file"]
            if existing.exists():
                record(manifest, sid, existing, cfg, "manual")
            else:
                manual.append((sid, cfg))
            continue
        try:
            if "api" in cfg:
                fetch_panelapp(sid, cfg, manifest)
            elif "url" in cfg:
                fetch_url(sid, cfg, manifest)
        except Exception as e:  # keep going; a release can still be built from what fetched
            print(f"  {sid}: FAILED {e}")
    MANIFEST.write_text(json.dumps(manifest, indent=1))
    if manual:
        print("\nManual sources still missing. Place each file at sources/<id>/<file>:")
        for sid, cfg in manual:
            print(f"  {sid:28s} {cfg['file']:34s} {cfg.get('landing','')}")
            if cfg.get("note"):
                print(f"  {'':28s} {cfg['note']}")


if __name__ == "__main__":
    main(sys.argv[1:])
