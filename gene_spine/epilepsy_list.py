"""Derived lists: the epilepsy and NDD gene list, and the organisation list.

    python -m gene_spine.epilepsy_list          # after gene_spine.build

Runs on the spine the build wrote, so it has to run after it. `make lists` and the workflow both
call it in that order. Keeping it in the package rather than in a scratch script is the point: a
published CSV that no code in the repository can reproduce would go stale against the spine at
the next scheduled build, and its hash in the README would quietly stop matching.

outputs/epilepsy_ndd_gene_list.csv is the union of three branches:

  (a) epilepsy        on Genes4Epilepsy
  (b) ndd_only        on an NDD or ID source, not in (a)
  (c) epilepsy_other  on an epilepsy panel, or on ClinGen or Orphadata naming epilepsy,
                      seizures or encephalopathy, and in neither branch above

`both` marks a gene on Genes4Epilepsy and on an NDD source. `tier` and `domain` are the spine's
own values, carried as description: a gene is on this list because a source names it, not
because it reached a tier.

outputs/epilepsy_orgs.csv is one row per organisation and gene, deduplicated on normalised
website host, every row status=candidate. curated/organizations.csv is never written here.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import urllib.parse
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
SRC = ROOT / "sources"

NDD_FLAGS = ["on_panelapp_au_id", "on_panelapp_ge_id", "on_sysndd", "on_g2p_dd", "on_sfari"]
EPI_PANEL_FLAGS = ["on_panelapp_au_epilepsy", "on_panelapp_ge_epilepsy"]
# \bDEE\b rather than a bare substring: "dee" turns up inside ordinary words.
EPI_TEXT = re.compile(r"epilep|seizure|encephalopath|\bDEE\b", re.I)

# Organisation sources come from the config, every entry whose role is org_candidates, so this
# list cannot drift from the declarations the way a hardcoded one did: the first version pointed
# at epilepsylive_orgs.csv in the repository root, which is a deleted duplicate, and silently
# lost the largest source.
ORG_ROLE = "org_candidates"
CURATED_ORGS = "curated/organizations.csv"


def _truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def _flagset(spine: pd.DataFrame, cols: list[str]) -> pd.Series:
    out = pd.Series(False, index=spine.index)
    for c in cols:
        if c in spine.columns:
            out = out | _truthy(spine[c])
    return out


def _codes(cell: object) -> list[str]:
    if pd.isna(cell):
        return []
    return [c.strip() for c in str(cell).replace("|", ";").split(";") if c.strip()]


def mondo_xrefs(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]], str]:
    """ORPHA and OMIM to MONDO, from the ontology's xref lines."""
    raw = path.read_bytes()
    orpha: dict[str, set[str]] = {}
    omim: dict[str, set[str]] = {}
    term = None
    for line in raw.decode("utf-8", errors="replace").splitlines():
        if line == "[Term]":
            term = None
        elif line.startswith("id: MONDO:"):
            term = line[4:].strip()
        elif line.startswith("xref: ") and term:
            x = line[6:].split("{")[0].split("!")[0].strip()
            if x.startswith("Orphanet:"):
                orpha.setdefault(x.split(":", 1)[1], set()).add(term)
            elif x.startswith("OMIM:"):
                omim.setdefault(x.split(":", 1)[1], set()).add(term)
    return orpha, omim, hashlib.sha256(raw).hexdigest()


def _host(url: object) -> str:
    """Normalised website host, the identity of an organisation for deduplication."""
    if pd.isna(url) or not str(url).strip():
        return ""
    u = str(url).strip()
    if "//" not in u:
        u = "https://" + u
    host = urllib.parse.urlparse(u).netloc.lower()
    return host[4:] if host.startswith("www.") else host


def _source_path(sid: str, cfg: dict) -> Path | None:
    """Same candidate order the spine build uses: sources/<id>/, sources/, then the root."""
    f = cfg.get("file")
    if not f:
        return None
    for cand in (SRC / sid / f, SRC / f, ROOT / f):
        if cand.exists():
            return cand
    return None


def _pick(df: pd.DataFrame, names: list[str] | str | None) -> str | None:
    for n in ([names] if isinstance(names, str) else (names or [])):
        if n in df.columns:
            return n
    return None


def build_org_list(spine: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    frames = []
    sources = [(sid, c) for sid, c in cfg.items()
               if isinstance(c, dict) and c.get("role") == ORG_ROLE]
    sources.append(("organizations_curated", {"file": CURATED_ORGS,
                                              "columns": {"gene": ["symbol"],
                                                          "org_name": ["org_name"],
                                                          "url": ["url"]}}))
    for sid, scfg in sources:
        p = (ROOT / CURATED_ORGS) if sid == "organizations_curated" else _source_path(sid, scfg)
        if p is None or not p.exists():
            print(f"  absent : {sid}")
            continue
        cols = scfg.get("columns") or {}
        df = pd.read_csv(p, low_memory=False)
        gene_col = _pick(df, cols.get("gene")) or _pick(df, ["gene", "symbol"])
        name_col = _pick(df, cols.get("org_name")) or _pick(df, ["org_name"])
        if df.empty or not gene_col or not name_col:
            print(f"  skipped: {sid} ({len(df)} rows, gene={gene_col}, org_name={name_col})")
            continue
        df = df.rename(columns={name_col: "org_name"})
        # fda_pfdd rows without a host describe an FDA-led session; the config names the column
        # to fall back to so those rows keep a usable name instead of an empty one.
        fb_col = _pick(df, cols.get("fallback_name"))
        if fb_col:
            blank = df["org_name"].isna() | (df["org_name"].astype(str).str.strip() == "")
            df.loc[blank, "org_name"] = df.loc[blank, fb_col]
        df = df[df["org_name"].notna() & (df["org_name"].astype(str).str.strip() != "")]
        url_col = _pick(df, cols.get("url"))
        df["url"] = df[url_col] if url_col else ""
        df = df.dropna(subset=[gene_col]).copy()
        df["symbol"] = df[gene_col].astype(str).str.split(";")
        df = df.explode("symbol")
        df["symbol"] = df["symbol"].str.strip().str.upper()
        df = df[df["symbol"] != ""]
        df["source_files"] = p.relative_to(ROOT).as_posix()
        for col in ("org_name", "url", "country", "org_type", "registry_platform",
                    "coalitions"):
            if col not in df.columns:
                df[col] = ""
        frames.append(df[["org_name", "symbol", "url", "country", "org_type",
                          "registry_platform", "coalitions", "source_files"]])
        print(f"  {sid:28} {len(df):>5} org-gene rows  "
              f"({p.relative_to(ROOT).as_posix()})")
    if not frames:
        return pd.DataFrame(columns=["org_name", "symbol", "hgnc_id", "url", "country",
                                     "org_type", "registry_platform", "coalitions",
                                     "source_files", "status", "org_key"])
    orgs = pd.concat(frames, ignore_index=True)
    orgs["host"] = orgs["url"].map(_host)
    # The host identifies an organisation where there is one, the name otherwise, so two
    # spellings of one charity collapse into a single row per gene.
    orgs["org_key"] = orgs["host"].where(
        orgs["host"] != "", orgs["org_name"].astype(str).str.strip().str.lower())
    before = len(orgs)
    orgs = (orgs.sort_values(["org_key", "symbol", "source_files"])
                .groupby(["org_key", "symbol"], as_index=False)
                .agg({"org_name": "first", "url": "first", "country": "first",
                      "org_type": "first", "registry_platform": "first",
                      "coalitions": "first",
                      "source_files": lambda s: ";".join(sorted(set(s)))}))
    print(f"  deduplicated on normalised host: {before:,} -> {len(orgs):,} (org, gene) rows")
    sym_to_hgnc = dict(zip(spine["symbol"].astype(str).str.upper(), spine["hgnc_id"]))
    orgs["hgnc_id"] = orgs["symbol"].map(sym_to_hgnc).fillna("")
    unresolved = orgs.loc[orgs["hgnc_id"] == "", "symbol"].nunique()
    if unresolved:
        print(f"  {int((orgs['hgnc_id'] == '').sum())} rows name something that is not a spine "
              f"symbol ({unresolved} distinct labels), kept with a blank hgnc_id")
    orgs["status"] = "candidate"
    return orgs


def main() -> int:
    spine = pd.read_csv(OUT / "gene_spine.csv", low_memory=False)
    long = pd.read_csv(OUT / "gene_source_long.csv", low_memory=False).fillna("")
    cfg = yaml.safe_load((ROOT / "config" / "sources.yaml").read_text(encoding="utf-8"))
    print(f"spine: {len(spine):,} genes")

    print("\norganisation sources:")
    orgs = build_org_list(spine, cfg)
    orgs_out = orgs[["org_name", "symbol", "hgnc_id", "url", "country", "org_type",
                     "registry_platform", "coalitions", "source_files", "status"]]
    orgs_out = orgs_out.sort_values(["org_name", "symbol"])
    orgs_path = OUT / "epilepsy_orgs.csv"
    with orgs_path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("# One row per (organisation, gene) from the organisation sources, "
                 "deduplicated on normalised website host.\n")
        fh.write("# status=candidate throughout: nothing here has been reviewed, and "
                 "curated/organizations.csv is untouched.\n")
        orgs_out.to_csv(fh, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"wrote {orgs_path.relative_to(ROOT).as_posix()}: {len(orgs_out):,} rows, "
          f"{orgs_out['org_name'].nunique():,} organisations")

    in_a = _truthy(spine["on_genes4epilepsy"])
    ndd_present = [c for c in NDD_FLAGS if c in spine.columns]
    ndd_absent = [c for c in NDD_FLAGS if c not in spine.columns]
    ndd = _flagset(spine, NDD_FLAGS)
    in_b = ndd & ~in_a

    epi_panel = _flagset(spine, EPI_PANEL_FLAGS)
    clingen_epi = pd.Series(False, index=spine.index)
    if "pheno_clingen" in spine.columns:
        clingen_epi = spine["pheno_clingen"].astype(str).str.contains(EPI_TEXT, na=False)
    orpha_hits = long[(long["source_id"] == "orphadata_genes")
                      & long["phenotype"].astype(str).str.contains(EPI_TEXT, na=False)]
    orpha_epi = spine["hgnc_id"].isin(set(orpha_hits["hgnc_id"]))
    in_c = (epi_panel | clingen_epi | orpha_epi) & ~in_a & ~ndd
    union = in_a | ndd | in_c

    print(f"\nNDD/ID sources used   : {', '.join(c[3:] for c in ndd_present)}")
    print(f"NDD/ID sources absent : {', '.join(c[3:] for c in ndd_absent) or 'none'}")
    print(f"(a) epilepsy            : {int(in_a.sum()):>6,}")
    print(f"(b) ndd_only            : {int(in_b.sum()):>6,}")
    print(f"(c) epilepsy_other      : {int(in_c.sum()):>6,}"
          f"   (panel {int((epi_panel & ~in_a & ~ndd).sum())},"
          f" ClinGen {int((clingen_epi & ~in_a & ~ndd).sum())},"
          f" Orphadata {int((orpha_epi & ~in_a & ~ndd).sum())})")
    print(f"union                   : {int(union.sum()):>6,}")

    sel = spine[union].copy()
    sel["list_membership"] = [
        "both" if a and n else ("epilepsy" if a else ("ndd_only" if n else "epilepsy_other"))
        for a, n in zip(in_a[union], ndd[union])
    ]

    # The Genes4Epilepsy phenotype, from the snapshot the build already read.
    g4e_cfg = cfg.get("genes4epilepsy", {})
    g4e_path = SRC / "genes4epilepsy" / str(g4e_cfg.get("file", ""))
    pheno, pheno_by_symbol = {}, {}
    if g4e_path.exists():
        g4e = pd.read_csv(g4e_path, sep="\t", dtype=str).fillna("")
        col = next((c for c in g4e.columns if c.lower().startswith("phenotype")), None)
        if col:
            pheno = dict(zip(g4e["HGNC_ID"].str.strip(), g4e[col].str.strip()))
            pheno_by_symbol = dict(zip(g4e["Gene"].str.strip().str.upper(),
                                       g4e[col].str.strip()))
    sel["genes4epilepsy_phenotype"] = [
        pheno.get(str(h).strip()) or pheno_by_symbol.get(str(s).strip().upper(), "")
        for h, s in zip(sel["hgnc_id"], sel["symbol"])
    ]
    sel.loc[~in_a[union].values, "genes4epilepsy_phenotype"] = ""

    mondo_note = "not resolved: no mondo.obo snapshot"
    sel["mondo_ids"] = ""
    mondo_path = SRC / "mondo" / str((cfg.get("mondo") or {}).get("file", "mondo.obo"))
    if mondo_path.exists():
        orpha_x, omim_x, sha = mondo_xrefs(mondo_path)
        mondo_note = (f"{mondo_path.relative_to(ROOT).as_posix()} sha256 {sha[:16]}, "
                      f"{len(orpha_x):,} ORPHA and {len(omim_x):,} OMIM xrefs")

        def resolve(row) -> str:
            ids: set[str] = set()
            for c in _codes(row["orphacodes"]):
                ids |= orpha_x.get(c.replace("ORPHA:", "").strip(), set())
            for c in _codes(row["omim_id"]):
                ids |= omim_x.get(c.replace("OMIM:", "").strip(), set())
            return ";".join(sorted(ids))

        sel["mondo_ids"] = sel.apply(resolve, axis=1)
    print(f"mondo_ids: {int((sel['mondo_ids'] != '').sum()):,} of {len(sel):,} resolved "
          f"[{mondo_note}]")

    counts = orgs.groupby("symbol")["org_key"].nunique() if len(orgs) else pd.Series(dtype=int)
    names = (orgs.groupby("symbol")["org_name"]
             .apply(lambda s: "; ".join(sorted({str(v) for v in s.dropna()})))
             if len(orgs) else pd.Series(dtype=str))
    upper = sel["symbol"].astype(str).str.upper()
    sel["org_count"] = upper.map(counts).fillna(0).astype(int)
    sel["org_names"] = upper.map(names).fillna("")

    cols = ["hgnc_id", "symbol", "gene_name", "list_membership", "genes4epilepsy_phenotype",
            "tier", "domain", "n_sources", "sources", "orphacodes", "omim_id", "mondo_ids",
            "org_count", "org_names"]
    out = sel[cols].sort_values(["list_membership", "symbol"])

    sfari = cfg.get("sfari") or {}
    header = [
        "# Epilepsy and NDD gene list, written by gene_spine.epilepsy_list after the spine",
        "# build. Union of three branches:",
        "#   epilepsy / both : on Genes4Epilepsy ('both' also on an NDD source)",
        "#   ndd_only        : on an NDD or ID source, not on Genes4Epilepsy",
        "#   epilepsy_other  : on an epilepsy panel, or ClinGen or Orphadata naming epilepsy,",
        "#                     seizures or encephalopathy, and in neither branch above",
        "# tier and domain are the spine's own values, carried as description, not selection.",
        f"# NDD/ID sources: {', '.join(c[3:] for c in ndd_present)}.",
        f"# Declared but absent: {', '.join(c[3:] for c in ndd_absent) or 'none'}.",
        f"# SFARI: {sfari.get('version', 'n/a')}, exported {sfari.get('exported', 'n/a')}."
        " Its unscored rows count as presence and lift no tier.",
        "# genes4epilepsy_phenotype is blank for genes that are not on Genes4Epilepsy.",
        f"# mondo_ids from {mondo_note}.",
        "# org_count and org_names come from epilepsy_orgs.csv, blank where there is no match.",
    ]
    path = OUT / "epilepsy_ndd_gene_list.csv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        for line in header:
            fh.write(line + "\n")
        out.to_csv(fh, index=False, quoting=csv.QUOTE_MINIMAL)

    print(f"wrote {path.relative_to(ROOT).as_posix()}: {len(out):,} rows")
    print("  membership: " + ", ".join(
        f"{k} {v:,}" for k, v in out["list_membership"].value_counts().sort_index().items()))
    print(f"  with at least one organisation: {int((out['org_count'] > 0).sum()):,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
