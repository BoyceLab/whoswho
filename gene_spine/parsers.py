"""Source parsers.

Every parser returns a long-format DataFrame with these columns:

    source_id, symbol_in_source, hgnc_id_in_source, evidence, phenotype, extra

evidence is the source's own label, untouched (ClinGen "Definitive", PanelApp
"Green", SFARI "1", SysNDD "Definitive", Wang category text, SAGAS score).
Tiering happens later, in tiers.py, so the raw label stays auditable.

A parser never resolves HGNC IDs itself; build.py does that once, uniformly.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

COLS = ["source_id", "symbol_in_source", "hgnc_id_in_source", "evidence", "phenotype", "extra"]


def _pick(df: pd.DataFrame, candidates: list[str], required: bool = True) -> str | None:
    lower = {c.lower().strip(): c for c in df.columns}
    for cand in candidates or []:
        if cand.lower().strip() in lower:
            return lower[cand.lower().strip()]
    if required:
        raise KeyError(f"None of {candidates} found. Actual header: {list(df.columns)}")
    return None


def _frame(source_id, symbols, hgnc_ids=None, evidence=None, phenotype=None, extra=None):
    n = len(symbols)
    return pd.DataFrame({
        "source_id": [source_id] * n,
        "symbol_in_source": list(symbols),
        "hgnc_id_in_source": list(hgnc_ids) if hgnc_ids is not None else [""] * n,
        "evidence": list(evidence) if evidence is not None else [""] * n,
        "phenotype": list(phenotype) if phenotype is not None else [""] * n,
        "extra": list(extra) if extra is not None else [""] * n,
    })


# ---------------------------------------------------------------- tabular sources

def parse_table(source_id: str, cfg: dict, path: Path) -> pd.DataFrame:
    fmt = cfg.get("format", "csv")
    if fmt == "tsv":
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    elif fmt == "xlsx":
        df = pd.read_excel(path, sheet_name=cfg.get("sheet") or 0, dtype=str).fillna("")
    elif fmt == "json":
        # An API that answers with records rather than a file. `records_path` names the key
        # holding the list, so the envelope a service wraps its rows in stays in config.
        payload = json.loads(path.read_text(encoding="utf-8"))
        for key in str(cfg.get("records_path") or "data").split("."):
            if key and isinstance(payload, dict):
                payload = payload[key]
        if not isinstance(payload, list):
            raise ValueError(
                f"{source_id}: records_path {cfg.get('records_path')!r} did not reach a list")
        df = pd.DataFrame(payload).fillna("").astype(str)
    else:
        df = _read_csv_with_preamble(path, cfg.get("skiprows_until"))
    cols = cfg.get("columns", {})
    sym = _pick(df, cols.get("symbol", ["symbol"]))
    hid = _pick(df, cols.get("hgnc_id", []), required=False)
    ev_col = None
    for key in ("classification", "category", "confidence", "score", "evidence"):
        if key in cols:
            ev_col = _pick(df, cols[key], required=False)
            if ev_col:
                break
    ph_col = None
    for key in ("phenotype", "disease"):
        if key in cols:
            ph_col = _pick(df, cols[key], required=False)
            if ph_col:
                break
    extras = {}
    for key in ("inheritance", "moi", "omim", "mondo", "gcep", "syndromic", "monogenic"):
        if key in cols:
            c = _pick(df, cols[key], required=False)
            if c:
                extras[key] = df[c]
    extra = [json.dumps({k: str(v.iloc[i]) for k, v in extras.items()}) for i in range(len(df))] if extras else None

    if not ph_col and cols.get("phenotype"):
        print(f"  note: {source_id} has no phenotype column; header is {list(df.columns)}")
    out = _frame(
        source_id,
        df[sym].str.strip(),
        df[hid].str.strip() if hid else None,
        df[ev_col].str.strip() if ev_col else None,
        df[ph_col].str.strip() if ph_col else None,
        extra,
    )
    return out[out["symbol_in_source"] != ""]


def _read_csv_with_preamble(path: Path, header_token: str | None) -> pd.DataFrame:
    if not header_token:
        return pd.read_csv(path, dtype=str).fillna("")
    text = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    tok = header_token.lower()
    start = next((i for i, line in enumerate(text)
                  if line.lstrip().lstrip('"').lower().startswith(tok)), None)
    if start is None:
        raise KeyError(f"No header line starting with '{header_token}'. First lines: {text[:6]}")
    from io import StringIO
    body = "\n".join(line for line in text[start:] if not line.lstrip().lstrip('"').startswith("+"))
    return pd.read_csv(StringIO(body), dtype=str).fillna("")


# ---------------------------------------------------------------- PanelApp

def parse_panelapp(source_id: str, cfg: dict, path: Path) -> pd.DataFrame:
    data = json.loads(Path(path).read_text())
    genes = data.get("genes", data)
    rows = []
    for g in genes:
        gd = g.get("gene_data", {})
        rows.append({
            "symbol": gd.get("gene_symbol") or g.get("gene_symbol", ""),
            "hgnc": gd.get("hgnc_id", ""),
            "conf": g.get("confidence_level", ""),
            "pheno": "; ".join(g.get("phenotypes", []) or []),
            "extra": json.dumps({"moi": g.get("mode_of_inheritance", ""),
                                 "panel": data.get("name", ""),
                                 "panel_version": data.get("version", "")}),
        })
    df = pd.DataFrame(rows)
    conf_map = {"3": "Green", "2": "Amber", "1": "Red", "0": "None"}
    df["conf"] = df["conf"].astype(str).map(lambda x: conf_map.get(x, x))
    return _frame(source_id, df["symbol"], df["hgnc"], df["conf"], df["pheno"], df["extra"])


# ---------------------------------------------------------------- Simons Searchlight PDF

CNV_RE = re.compile(r"^\d{1,2}[pq]\d")


def parse_simons_pdf(source_id: str, cfg: dict, path: Path, known_symbols: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (gene rows, cnv rows). Tokens are checked against HGNC symbols so
    heading words in the PDF are not mistaken for genes."""
    try:
        import pdfplumber
    except ImportError as e:
        raise SystemExit("pip install pdfplumber") from e
    tokens = []
    with pdfplumber.open(path) as pdf:
        for page in pdf.pages:
            txt = page.extract_text() or ""
            tokens.extend(re.split(r"[\s,]+", txt))
    genes, cnvs = [], []
    for t in tokens:
        t = t.strip().strip("·•")
        if not t:
            continue
        if CNV_RE.match(t) or "deletion" in t.lower() or "duplication" in t.lower():
            cnvs.append(t)
        elif t.upper() in known_symbols:
            genes.append(t.upper())
    genes = sorted(set(genes))
    cnv_df = pd.DataFrame({"source_id": source_id, "region_text": sorted(set(cnvs))})
    return _frame(source_id, genes, evidence=["studied"] * len(genes)), cnv_df


# ---------------------------------------------------------------- Orphadata product 6

def parse_orphadata_genes(source_id: str, cfg: dict, path: Path) -> pd.DataFrame:
    """Gene -> ORPHAcode bridge. One row per (gene, disorder) association."""
    import xml.etree.ElementTree as ET
    root = ET.parse(path).getroot()
    rows = []
    for dis in root.iter("Disorder"):
        orpha = (dis.findtext("OrphaCode") or "").strip()
        dname = (dis.findtext("Name") or "").strip()
        for assoc in dis.iter("DisorderGeneAssociation"):
            gene = assoc.find("Gene")
            if gene is None:
                continue
            sym = (gene.findtext("Symbol") or "").strip()
            hgnc = ""
            for ref in gene.iter("ExternalReference"):
                if (ref.findtext("Source") or "") == "HGNC":
                    hgnc = "HGNC:" + (ref.findtext("Reference") or "").strip()
            atype = assoc.findtext("DisorderGeneAssociationType/Name") or ""
            status = assoc.findtext("DisorderGeneAssociationStatus/Name") or ""
            rows.append({"sym": sym, "hgnc": hgnc, "ev": status.strip(), "ph": dname,
                         "extra": json.dumps({"orphacode": orpha, "association_type": atype.strip()})})
    df = pd.DataFrame(rows)
    return _frame(source_id, df["sym"], df["hgnc"], df["ev"], df["ph"], df["extra"])
