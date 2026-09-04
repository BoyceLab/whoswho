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


def load_sources(cfg_all: dict, resolver: HgncResolver):
    long_frames, cnv_frames, missing = [], [], []
    known = set(resolver.by_symbol)
    for sid, cfg in cfg_all.items():
        if sid == "hgnc" or cfg.get("role") in ("org_candidates", "gene_annotation"):
            continue
        path = SRC / sid / cfg["file"]
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
        tier, disputed = tiers.tier_for(zip(g["source_id"], g["evidence"]))
        extras = [json.loads(e) if e else {} for e in g["extra"]]
        rec["tier"] = tier
        rec["disputed"] = disputed
        rec["domain"] = tiers.domain_for(zip(g["source_id"], g["evidence"], g["phenotype"], extras), epi_terms, ndd_terms)
        orpha = sorted({x.get("orphacode") for x in extras if x.get("orphacode")})
        rec["orphacodes"] = ";".join(orpha)
        rows.append(rec)
    wide = pd.DataFrame(rows)
    on_cols = [c for c in wide.columns if c.startswith("on_")]
    wide[on_cols] = wide[on_cols].fillna(False)
    lead = ["hgnc_id", "symbol", "gene_name", "tier", "domain", "disputed", "n_sources", "sources", "orphacodes", "omim_id", "ensembl_gene_id", "locus_group"]
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
    cands = []
    for sid, cfg in cfg_all.items():
        if cfg.get("role") != "org_candidates":
            continue
        path = SRC / sid / cfg["file"]
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
    cand = pd.DataFrame(cands, columns=ORG_COLS) if cands else pd.DataFrame(columns=ORG_COLS)
    # A candidate is dropped when a curated row already has the same HGNC ID and same normalized URL host.
    def host(u):
        u = (u or "").lower().replace("https://", "").replace("http://", "").replace("www.", "")
        return u.split("/")[0]
    def norm(n):
        return re.sub(r"[^a-z0-9]", "", (n or "").lower())
    cur_keys = {(h, host(u)) for h, u in zip(curated["hgnc_id"], curated["url"]) if u}
    cur_names = {(h, norm(n)) for h, n in zip(curated["hgnc_id"], curated["org_name"])}
    cand = cand[[(h, host(u)) not in cur_keys and (h, norm(n)) not in cur_names
                 for h, u, n in zip(cand["hgnc_id"], cand["url"], cand["org_name"])]]
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
    for sid, cfg in cfg_all.items():
        if cfg.get("role") != "gene_annotation":
            continue
        path = SRC / sid / cfg["file"]
        if not path.exists():
            continue
        df = pd.read_csv(path, dtype=str).fillna("")
        cols = cfg["columns"]
        gene_c = parsers._pick(df, cols["gene"])
        group_c = parsers._pick(df, cols.get("group", []), required=False)
        ann = cfg.get("annotate", {})
        gene_rows, cnv_rows = {}, []
        for _, r in df.iterrows():
            is_cnv = group_c and r[group_c].strip().lower() == "cnv"
            hid, _ = (None, "") if is_cnv else resolver.resolve(r[gene_c])
            if hid:
                gene_rows[hid] = {k: r[parsers._pick(df, [v])] for k, v in ann.items()}
            else:
                cnv_rows.append({"source_id": sid, "region_text": r[gene_c], **{k: r[parsers._pick(df, [v])] for k, v in ann.items()}})
        for k in ann:
            wide[k] = wide["hgnc_id"].map(lambda h: gene_rows.get(h, {}).get(k, ""))
        wide["on_" + sid] = wide["hgnc_id"].isin(gene_rows)
        if cnv_rows:
            cnv = pd.concat([cnv, pd.DataFrame(cnv_rows)], ignore_index=True)
        print(f"  {sid}: {len(gene_rows)} genes annotated, {len(cnv_rows)} CNV rows")
    return wide, cnv


# ---------------------------------------------------------------- main

def main():
    cfg_all = yaml.safe_load((ROOT / "config" / "sources.yaml").read_text())
    OUT.mkdir(exist_ok=True)
    hgnc_path = SRC / "hgnc" / cfg_all["hgnc"]["file"]
    if not hgnc_path.exists():
        raise SystemExit("HGNC complete set missing: python -m gene_spine.fetch hgnc")
    print("Loading HGNC")
    resolver = HgncResolver.from_file(hgnc_path)
    print(f"  HGNC: {len(resolver.table)} genes, e.g. {list(resolver.table['symbol'].head(3))}")
    print("Parsing sources")
    long, cnv = load_sources(cfg_all, resolver)
    long, unresolved = resolve_all(long, resolver)
    wide = widen(long, resolver, cfg_all)
    wide, cnv = apply_annotations(wide, cnv, resolver, cfg_all)
    gene_orpha = wide[["hgnc_id", "orphacodes"]]
    org, gaps = build_org_layer(wide, resolver, cfg_all, gene_orpha)

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
    release = {
        "built": date.today().isoformat(),
        "sources": manifest,
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
