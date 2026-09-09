"""Build the gene spine and the organization layer from whatever is in sources/.

    python -m gene_spine.build

Outputs (all in outputs/):
  gene_spine.csv / .json     one row per HGNC gene, wide presence + evidence columns
  gene_source_long.csv       one row per (gene, source) with the raw evidence label
  gene_orpha.csv             gene -> ORPHAcode bridge for the organization join
  cnv_regions.csv            CNV entries that are not genes (Simons Searchlight)
  unresolved_symbols.csv     symbols that did not map to HGNC, by source
  org_layer.csv              curated organizations + candidates, keyed to HGNC / ORPHAcode
  org_gaps.csv               genes with no organization row at all, by tier
  release.json               dates, versions, hashes, tier rule text, row counts
"""
from __future__ import annotations

import json
import re
from datetime import date
from pathlib import Path

import pandas as pd
import yaml

from . import parsers, tiers
from .hgnc import HgncResolver

ROOT = Path(__file__).resolve().parents[1]
SRC, OUT, CUR = ROOT / "sources", ROOT / "outputs", ROOT / "curated"


def source_path(sid: str, cfg: dict) -> Path:
    """Manual files may be uploaded to sources/<id>/, to sources/, or to the repo root
    (GitHub's uploader flattens folders). The first match wins."""
    for cand in (SRC / sid / cfg["file"], SRC / cfg["file"], ROOT / cfg["file"]):
        if cand.exists():
            return cand
    return SRC / sid / cfg["file"]


def load_sources(cfg_all: dict, resolver: HgncResolver):
    long_frames, cnv_frames, missing = [], [], []
    known = set(resolver.by_symbol)
    for sid, cfg in cfg_all.items():
        if sid == "hgnc" or cfg.get("role") in ("org_candidates", "gene_annotation", "gene_registries"):
            continue
        path = source_path(sid, cfg)
        if not path.exists():
            missing.append(sid)
            continue
        fmt = cfg.get("format")
        if "api" in cfg:
            df = parsers.parse_panelapp(sid, cfg, path)
        elif fmt == "pdf":
            df, cnv = parsers.parse_simons_pdf(sid, cfg, path, known)
            cnv_frames.append(cnv)
        elif fmt == "orphadata_xml":
            df = parsers.parse_orphadata_genes(sid, cfg, path)
        else:
            df = parsers.parse_table(sid, cfg, path)
            if sid == "wang_omim_2023":
                df = _apply_wang_map(df, cfg.get("category_map", {}))
        long_frames.append(df)
        print(f"  {sid}: {len(df)} rows")
    if missing:
        print("  skipped (no file):", ", ".join(missing))
    long = pd.concat(long_frames, ignore_index=True) if long_frames else pd.DataFrame(columns=parsers.COLS)
    cnv = pd.concat(cnv_frames, ignore_index=True) if cnv_frames else pd.DataFrame(columns=["source_id", "region_text"])
    return long, cnv


def _apply_wang_map(df: pd.DataFrame, cmap: dict) -> pd.DataFrame:
    def m(ev):
        low = (ev or "").lower()
        for frag, dom in cmap.items():
            if frag in low:
                return dom
        return ""
    df = df.copy()
    df["extra"] = [json.dumps({**(json.loads(e) if e else {}), "domain": m(ev)}) for e, ev in zip(df["extra"], df["evidence"])]
    return df


def resolve_all(long: pd.DataFrame, resolver: HgncResolver):
    ids, hows = [], []
    for sym, hid in zip(long["symbol_in_source"], long["hgnc_id_in_source"]):
        h, how = resolver.resolve(sym, hid)
        ids.append(h or "")
        hows.append(how)
    long = long.assign(hgnc_id=ids, resolved_by=hows)
    unresolved = long[long["hgnc_id"] == ""][["source_id", "symbol_in_source", "hgnc_id_in_source", "resolved_by"]]
    stats = long.assign(ok=long["hgnc_id"] != "").groupby("source_id")["ok"].agg(["sum", "count"])
    for sid, r in stats.iterrows():
        print(f"  resolved {sid}: {int(r['sum'])}/{int(r['count'])}")
    resolved = long[long["hgnc_id"] != ""].copy()
    if resolved.empty:
        sample = long.head(5)[["source_id", "symbol_in_source", "hgnc_id_in_source", "resolved_by"]].to_string()
        raise SystemExit("No symbol from any source resolved to an HGNC ID.\n"
                         f"HGNC table: {len(resolver.table)} rows, first symbols {list(resolver.table['symbol'].head(5))}\n"
                         f"Sample source rows:\n{sample}")
    return resolved, unresolved.drop_duplicates()


ENRICH_ONLY = {"orphadata_genes"}   # identifier and disorder bridges, never a reason to include a gene


def in_scope(long: pd.DataFrame, cfg_all: dict) -> pd.DataFrame:
    """Keep genes with at least one in-scope assertion.

    A gene enters the spine because an epilepsy or NDD source lists it, or because
    a broad source (ClinGen) records a relationship whose disease text matches the
    epilepsy or NDD terms in the config. Orphadata and non-matching ClinGen rows
    still enrich a gene that is already in scope; they never pull one in. Without
    this, a gene like CFTR arrives through the disorder bridge and is displayed in
    an epilepsy directory at full strength.
    """
    epi = cfg_all.get("clingen", {}).get("epilepsy_terms", [])
    ndd = cfg_all.get("clingen", {}).get("ndd_terms", [])
    terms = [t.lower() for t in epi + ndd]

    def qualifies(row) -> bool:
        sid = row.source_id
        if sid in ENRICH_ONLY:
            return False
        if sid == "clingen":
            text = f"{row.phenotype}".lower()
            return any(t in text for t in terms)
        return True

    keep = {h for h, g in long.groupby("hgnc_id") if any(qualifies(r) for r in g.itertuples())}
    dropped = long[~long["hgnc_id"].isin(keep)]
    if len(dropped):
        print(f"  out of scope: {dropped['hgnc_id'].nunique()} genes with no epilepsy or NDD assertion")
    return long[long["hgnc_id"].isin(keep)].copy()


def widen(long: pd.DataFrame, resolver: HgncResolver, cfg_all: dict) -> pd.DataFrame:
    epi_terms = cfg_all.get("clingen", {}).get("epilepsy_terms", [])
    ndd_terms = cfg_all.get("clingen", {}).get("ndd_terms", [])
    rows = []
    for hid, g in long.groupby("hgnc_id", sort=True):
        base = resolver.info(hid)
        rec = dict(base)
        rec["sources"] = ";".join(sorted(g["source_id"].unique()))
        rec["n_sources"] = g["source_id"].nunique()
        for sid, gg in g.groupby("source_id"):
            rec[f"on_{sid}"] = True
            ev = ";".join(sorted({e for e in gg["evidence"] if e}))
            if ev:
                rec[f"ev_{sid}"] = ev
            ph = "; ".join(sorted({p for p in gg["phenotype"] if p}))
            if ph and sid != "orphadata_genes":
                rec[f"pheno_{sid}"] = ph[:500]
        tier, disputed, disputed_note = tiers.tier_for(zip(g["source_id"], g["evidence"]))
        extras = [json.loads(e) if e else {} for e in g["extra"]]
        rec["tier"] = tier
        rec["disputed"] = disputed
        rec["disputed_note"] = disputed_note
        rec["domain"] = tiers.domain_for(zip(g["source_id"], g["evidence"], g["phenotype"], extras), epi_terms, ndd_terms)
        orpha = sorted({x.get("orphacode") for x in extras if x.get("orphacode")})
        rec["orphacodes"] = ";".join(orpha)
        rows.append(rec)
    wide = pd.DataFrame(rows)
    on_cols = [c for c in wide.columns if c.startswith("on_")]
    wide[on_cols] = wide[on_cols].fillna(False)
    lead = ["hgnc_id", "symbol", "gene_name", "tier", "domain", "disputed", "disputed_note", "n_sources", "sources", "orphacodes", "synonyms", "omim_id", "ensembl_gene_id", "locus_group"]
    rest = sorted(c for c in wide.columns if c not in lead)
    return wide[lead + rest].sort_values(["tier", "symbol"]).reset_index(drop=True)


# ---------------------------------------------------------------- organization layer

ORG_COLS = ["org_id", "org_name", "url", "hgnc_id", "symbol", "orphacode", "country", "org_type",
            "registry_platform", "registry_url", "natural_history_study", "data_elements_documented",
            "irb_status", "maturity_stage", "coalitions", "contact", "status", "source", "reviewed_by",
            "review_date", "notes"]


def build_org_layer(wide: pd.DataFrame, resolver: HgncResolver, cfg_all: dict, gene_orpha: pd.DataFrame):
    cur_path = CUR / "organizations.csv"
    curated = pd.read_csv(cur_path, dtype=str).fillna("") if cur_path.exists() else pd.DataFrame(columns=ORG_COLS)
    curated["status"] = curated["status"].replace("", "curated")
    # A row may only read "verified" when someone signed it: reviewer, date, and a website
    # to check against. Anything else is a candidate however it was labelled.
    unsigned = (curated["status"].str.lower() == "curated") & (
        (curated["reviewed_by"].str.strip() == "") | (curated["review_date"].str.strip() == "") | (curated["url"].str.strip() == ""))
    if unsigned.any():
        print(f"  organizations: {int(unsigned.sum())} curated rows lack reviewer, date, or website; shown as awaiting review")
        curated.loc[unsigned, "status"] = "candidate"
    n_rejected = int((curated["status"].str.lower() == "rejected").sum())
    curated = curated[curated["status"].str.lower() != "rejected"]
    if n_rejected:
        print(f"  organizations: {n_rejected} rejected rows withheld from the public output")
    cands = []
    for sid, cfg in cfg_all.items():
        if cfg.get("role") != "org_candidates":
            continue
        path = source_path(sid, cfg)
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str).fillna("")
        cols = cfg["columns"]
        name_c = parsers._pick(df, cols["org_name"])
        fb_c = parsers._pick(df, cols.get("fallback_name", []), required=False)
        gene_c = parsers._pick(df, cols.get("gene", []), required=False)
        url_c = parsers._pick(df, cols.get("url", []), required=False)
        note_cs = [c for c in (parsers._pick(df, [n], required=False) for n in cols.get("note_fields", [])) if c]
        passthrough = {k: parsers._pick(df, [k] + cols.get(k, []), required=False)
                       for k in ORG_COLS if k not in ("org_name", "url", "hgnc_id", "symbol", "status", "source")}
        defaults = cfg.get("defaults", {})
        for _, r in df.iterrows():
            name = r[name_c].strip() or (r[fb_c].strip() if fb_c else "")
            if not name:
                continue
            genes = [g.strip() for g in str(r[gene_c]).replace(",", ";").split(";") if g.strip()] if gene_c else []
            base = {k: (r[c].strip() if c else "") for k, c in passthrough.items()}
            base.update({k: v for k, v in defaults.items() if not base.get(k)})
            if note_cs:
                extra = "; ".join(f"{c}: {r[c]}" for c in note_cs if r[c])
                base["notes"] = (base.get("notes", "") + ("; " if base.get("notes") else "") + extra).strip()
            base.update({"org_name": name, "url": (r[url_c].strip() if url_c else ""), "status": "candidate", "source": sid})
            for gsym in (genes or [""]):
                hid, _ = resolver.resolve(gsym) if gsym else (None, "")
                cands.append({**base, "hgnc_id": hid or "", "symbol": resolver.current_symbol(hid) if hid else gsym})
    cand = pd.DataFrame(cands, columns=ORG_COLS) if cands else pd.DataFrame(columns=ORG_COLS, dtype=str)
    # A candidate is dropped when a curated row already has the same HGNC ID and same normalized URL host.
    # Hosts that many unrelated groups share. Two Facebook groups are two
    # organizations, so identity on these needs the path, not just the domain.
    SHARED_HOSTS = ("facebook.com", "groups.io", "sites.google.com", "wixsite.com",
                    "wordpress.com", "blogspot.com", "linktr.ee", "instagram.com",
                    "x.com", "twitter.com", "iamrare.org", "citizen.health", "redcap.link")

    def host(u):
        u = (u or "").lower().strip()
        u = u.replace("https://", "").replace("http://", "").replace("www.", "")
        if not u:
            return ""
        domain = u.split("/")[0]
        if any(domain == h or domain.endswith("." + h) for h in SHARED_HOSTS):
            return u.rstrip("/")          # keep the path: the group is the path
        return domain
    def norm(n):
        return re.sub(r"[^a-z0-9]", "", (n or "").lower())
    cur_keys = {(h, host(u)) for h, u in zip(curated["hgnc_id"], curated["url"]) if u}
    cur_names = {(h, norm(n)) for h, n in zip(curated["hgnc_id"], curated["org_name"])}
    cand = cand[[(h, host(u)) not in cur_keys and (h, norm(n)) not in cur_names
                 for h, u, n in zip(cand["hgnc_id"], cand["url"], cand["org_name"])]]
    # Merge candidates that describe the same organization for the same gene: same website host,
    # or same normalized name. Fields fill from the first non-empty value; the origin ids are kept.
    cand = cand.copy().reset_index(drop=True)
    if cand.empty:
        org = curated.copy().fillna("")
        org = org.merge(gene_orpha.rename(columns={"orphacodes": "_orpha"}), on="hgnc_id", how="left")
        org["orphacode"] = org["orphacode"].where(org["orphacode"] != "", org["_orpha"].fillna(""))
        org = org.drop(columns="_orpha")
        covered = set(org["hgnc_id"])
        gaps = wide.loc[~wide["hgnc_id"].isin(covered), ["hgnc_id", "symbol", "tier", "domain", "sources"]]
        return org[ORG_COLS], gaps
    # union-find over (gene, name) and (gene, host) so a row with a URL and a row without one still merge
    parent = list(range(len(cand)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    seen = {}
    for i, (h, u, n) in enumerate(zip(cand["hgnc_id"], cand["url"], cand["org_name"])):
        for key in ((h, "n:" + norm(n)), (h, "h:" + host(u)) if host(u) else None):
            if key is None:
                continue
            if key in seen:
                parent[find(i)] = find(seen[key])
            else:
                seen[key] = i
    cand["_key"] = [find(i) for i in range(len(cand))]
    merged = []
    for key, grp in cand.groupby("_key", sort=False):
        rec = {}
        for col in ORG_COLS:
            vals = [v for v in grp[col] if v]
            rec[col] = vals[0] if vals else ""
        rec["org_name"] = max(grp["org_name"], key=len)  # fullest name variant
        rec["source"] = ";".join(sorted({s for s in grp["source"] if s}))
        rec["coalitions"] = ";".join(sorted({c for cs in grp["coalitions"] for c in cs.split(";") if c}))
        rec["status"] = "candidate"
        merged.append(rec)
    cand = pd.DataFrame(merged, columns=ORG_COLS) if merged else pd.DataFrame(columns=ORG_COLS)
    private = {sid for sid, c in cfg_all.items() if c.get("private")}
    def public_source(v):
        parts = [p for p in str(v).split(";") if p]
        out = sorted({("editor curation" if p in private else p) for p in parts})
        return ";".join(out)
    cand["source"] = cand["source"].map(public_source)
    org = pd.concat([curated, cand], ignore_index=True).fillna("")
    org = org.merge(gene_orpha.rename(columns={"orphacodes": "_orpha"}), on="hgnc_id", how="left")
    org["orphacode"] = org["orphacode"].where(org["orphacode"] != "", org["_orpha"].fillna(""))
    org = org.drop(columns="_orpha")
    covered = set(org.loc[org["status"] == "curated", "hgnc_id"]) | set(org.loc[org["status"] == "candidate", "hgnc_id"])
    gaps = wide.loc[~wide["hgnc_id"].isin(covered), ["hgnc_id", "symbol", "tier", "domain", "sources"]]
    return org[ORG_COLS], gaps


# ---------------------------------------------------------------- gene annotations

def apply_annotations(wide: pd.DataFrame, cnv: pd.DataFrame, resolver: HgncResolver, cfg_all: dict):
    """Join platform pages (e.g. Simons Searchlight) onto gene rows by symbol; CNV rows go to cnv_regions."""
    unresolved_ann = []
    for sid, cfg in cfg_all.items():
        if cfg.get("role") != "gene_annotation":
            continue
        path = source_path(sid, cfg)
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str).fillna("")
        cols = cfg["columns"]
        gene_c = parsers._pick(df, cols["gene"])
        hid_c = parsers._pick(df, cols.get("hgnc_id", []), required=False)
        group_c = parsers._pick(df, cols.get("group", []), required=False)
        ann = cfg.get("annotate", {})
        for v in ann.values():
            parsers._pick(df, [v])  # fail early with the real header if a column is missing
        gene_rows, cnv_rows = {}, []
        for _, r in df.iterrows():
            is_cnv = group_c and r[group_c].strip().lower() == "cnv"
            hid, _ = (None, "") if is_cnv else resolver.resolve(r[gene_c], r[hid_c] if hid_c else None)
            if hid:
                gene_rows[hid] = {k: r[parsers._pick(df, [v])] for k, v in ann.items()}
            elif is_cnv:
                cnv_rows.append({"source_id": sid, "region_text": r[gene_c], **{k: r[parsers._pick(df, [v])] for k, v in ann.items()}})
            else:
                unresolved_ann.append({"source_id": sid, "symbol_in_source": r[gene_c], "hgnc_id_in_source": r[hid_c] if hid_c else "", "resolved_by": "unresolved"})
        missing = [h for h in gene_rows if h not in set(wide["hgnc_id"])]
        if missing and cfg.get("adds_genes"):
            extra = pd.DataFrame([{**resolver.info(h), "tier": "T3", "domain": "unclassified", "disputed": False,
                                   "n_sources": 0, "sources": ""} for h in missing])
            wide = pd.concat([wide, extra], ignore_index=True)
        for k in ann:
            wide[k] = wide["hgnc_id"].map(lambda h: gene_rows.get(h, {}).get(k, ""))
        wide["on_" + sid] = wide["hgnc_id"].isin(gene_rows)
        if cnv_rows:
            cnv = pd.concat([cnv, pd.DataFrame(cnv_rows)], ignore_index=True)
        print(f"  {sid}: {len(gene_rows)} genes annotated, {len(cnv_rows)} CNV rows, {sum(1 for u in unresolved_ann if u['source_id']==sid)} unresolved")
    on_cols = [c for c in wide.columns if c.startswith("on_")]
    wide[on_cols] = wide[on_cols].fillna(False).astype(bool)
    return wide.fillna(""), cnv, pd.DataFrame(unresolved_ann)


# ---------------------------------------------------------------- registries

def apply_registries(wide: pd.DataFrame, resolver: HgncResolver, cfg_all: dict) -> pd.DataFrame:
    """Attach the registry rows for each gene as a JSON list, so a card can render them
    without re-deriving anything. Genes with no registry get an empty list."""
    by_gene: dict[str, list] = {}
    for sid, cfg in cfg_all.items():
        if cfg.get("role") != "gene_registries":
            continue
        path = source_path(sid, cfg)
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str).fillna("")
        for _, r in df.iterrows():
            hid, _ = resolver.resolve(r.get("gene", ""))
            if not hid:
                continue
            by_gene.setdefault(hid, []).append({k: r.get(k, "") for k in
                ("platform", "registry_name", "registry_url", "platform_url", "platform_basis", "notes")})
        print(f"  {sid}: {len(df)} registry rows over {len(by_gene)} genes")
    wide["registries"] = wide["hgnc_id"].map(lambda h: json.dumps(by_gene.get(h, [])))
    return wide


# ---------------------------------------------------------------- shards and search index

SHARDS = 64


def shard_of(hgnc_id: str) -> int:
    return int(hgnc_id.split(":")[1]) % SHARDS


def write_shards(wide: pd.DataFrame, org: pd.DataFrame, out: Path):
    """The card page needs one gene, not the whole spine. Genes are split into
    SHARDS files by HGNC number, and each gene record carries its organization rows,
    so a card is two small fetches: the index, then its shard."""
    gdir = out / "genes"
    gdir.mkdir(exist_ok=True)
    org_by = {h: g.to_dict(orient="records") for h, g in org.groupby("hgnc_id") if h}
    buckets: dict[int, list] = {i: [] for i in range(SHARDS)}
    for rec in wide.to_dict(orient="records"):
        rec = {k: v for k, v in rec.items() if v not in ("", None, False)}
        rec["organizations"] = org_by.get(rec["hgnc_id"], [])
        buckets[shard_of(rec["hgnc_id"])].append(rec)
    for i, recs in buckets.items():
        (gdir / f"{i:02d}.json").write_text(json.dumps(recs, separators=(",", ":")))
    print(f"  wrote {SHARDS} gene shards")


def write_search_index(wide: pd.DataFrame, org: pd.DataFrame, cnv: pd.DataFrame, out: Path):
    """A small index for search: genes with synonyms and disorder text, organizations,
    conditions (from disorder and phenotype text), and CNV regions."""
    genes = []
    for r in wide.itertuples(index=False):
        pheno = " | ".join(str(getattr(r, c, "") or "") for c in wide.columns if c.startswith("pheno_"))
        genes.append({"s": r.symbol, "h": r.hgnc_id, "n": r.gene_name, "t": r.tier, "d": r.domain,
                      "y": getattr(r, "synonyms", "") or "", "c": getattr(r, "disorder", "") or "",
                      "p": pheno[:300], "k": shard_of(r.hgnc_id)})
    orgs = [{"o": r.org_name, "u": r.url, "s": r.symbol, "h": r.hgnc_id, "st": r.status,
             "ty": r.org_type, "co": r.country}
            for r in org.itertuples(index=False)]
    cnvs = [{"r": r.region_text, "src": r.source_id,
             "u": getattr(r, "simons_page_url", "") or ""} for r in cnv.itertuples(index=False)] if len(cnv) else []
    (out / "search_index.json").write_text(json.dumps({"genes": genes, "orgs": orgs, "cnvs": cnvs}, separators=(",", ":")))
    print(f"  search index: {len(genes)} genes, {len(orgs)} organization rows, {len(cnvs)} CNV regions")


# ---------------------------------------------------------------- main

def main():
    cfg_all = yaml.safe_load((ROOT / "config" / "sources.yaml").read_text())
    OUT.mkdir(exist_ok=True)
    hgnc_path = source_path("hgnc", cfg_all["hgnc"])
    if not hgnc_path.exists():
        raise SystemExit("HGNC complete set missing: python -m gene_spine.fetch hgnc")
    print("Loading HGNC")
    resolver = HgncResolver.from_file(hgnc_path)
    print(f"  HGNC: {len(resolver.table)} genes, e.g. {list(resolver.table['symbol'].head(3))}")
    print("Parsing sources")
    long, cnv = load_sources(cfg_all, resolver)
    long, unresolved = resolve_all(long, resolver)
    long = in_scope(long, cfg_all)
    wide = widen(long, resolver, cfg_all)
    wide, cnv, unresolved_ann = apply_annotations(wide, cnv, resolver, cfg_all)
    if len(unresolved_ann):
        unresolved = pd.concat([unresolved, unresolved_ann], ignore_index=True)
    wide = apply_registries(wide, resolver, cfg_all)
    gene_orpha = wide[["hgnc_id", "orphacodes"]]
    org, gaps = build_org_layer(wide, resolver, cfg_all, gene_orpha)

    write_shards(wide, org, OUT)
    write_search_index(wide, org, cnv, OUT)
    wide.to_csv(OUT / "gene_spine.csv", index=False)
    wide.to_json(OUT / "gene_spine.json", orient="records", indent=1)
    long.to_csv(OUT / "gene_source_long.csv", index=False)
    gene_orpha[gene_orpha["orphacodes"] != ""].to_csv(OUT / "gene_orpha.csv", index=False)
    cnv.to_csv(OUT / "cnv_regions.csv", index=False)
    unresolved.to_csv(OUT / "unresolved_symbols.csv", index=False)
    org.to_csv(OUT / "org_layer.csv", index=False)
    gaps.to_csv(OUT / "org_gaps.csv", index=False)

    manifest_path = SRC / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    # The published source list covers only sources in the current config that are not marked
    # private. Private entries are the editor's own working files: they are used in the build and
    # recorded in the manifest, but they are not public documents, so listing them would point
    # readers at something they cannot open. Sources dropped from the config are excluded too,
    # so a stale manifest entry does not linger on the page.
    public_sources = {sid: rec for sid, rec in manifest.items()
                      if sid in cfg_all and not cfg_all[sid].get("private")}
    release = {
        "built": date.today().isoformat(),
        "sources": public_sources,
        "tier_rule": tiers.__doc__,
        "genes_total": int(len(wide)),
        "by_tier": wide["tier"].value_counts().to_dict(),
        "by_domain": wide["domain"].value_counts().to_dict(),
        "unresolved_symbols": int(len(unresolved)),
        "org_rows": {"curated": int((org["status"] == "curated").sum()), "candidate": int((org["status"] == "candidate").sum())},
        "genes_without_org": int(len(gaps)),
    }
    (OUT / "release.json").write_text(json.dumps(release, indent=1))
    print("\nGenes:", release["genes_total"], "| by tier:", release["by_tier"])
    print("By domain:", release["by_domain"])
    print("Unresolved symbols:", release["unresolved_symbols"], "| genes without an org row:", release["genes_without_org"])


if __name__ == "__main__":
    main()
