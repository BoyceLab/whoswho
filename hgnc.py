"""HGNC identifier backbone.

Every list is keyed on HGNC ID before anything is joined. Symbols are only
used to get to an HGNC ID, and a symbol that does not resolve is written to
outputs/unresolved_symbols.csv rather than dropped silently.

Resolution order for a bare symbol:
  1. exact current symbol
  2. previous symbol (renamed genes)
  3. alias symbol (only if the alias maps to exactly one gene)
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import pandas as pd

HGNC_ID_RE = re.compile(r"^(?:HGNC:)?(\d+)$", re.I)


@dataclass
class HgncResolver:
    table: pd.DataFrame
    by_symbol: dict = field(default_factory=dict)
    by_prev: dict = field(default_factory=dict)
    by_alias: dict = field(default_factory=dict)
    ambiguous_alias: set = field(default_factory=set)

    @classmethod
    def from_file(cls, path) -> "HgncResolver":
        t = _read_hgnc(path)
        keep = ["hgnc_id", "symbol", "name", "locus_group", "status",
                "prev_symbol", "alias_symbol", "omim_id", "ensembl_gene_id"]
        missing = [c for c in ("hgnc_id", "symbol") if c not in t.columns]
        if missing:
            raise SystemExit(f"HGNC file lacks {missing}. Header found: {list(t.columns)[:20]}")
        keep = [c for c in keep if c in t.columns]
        t = t[keep].fillna("")
        t["hgnc_id"] = t["hgnc_id"].str.strip().str.strip('"')
        t["symbol"] = t["symbol"].str.strip().str.strip('"')
        t = t[t["hgnc_id"].str.match(r"^HGNC:\d+$", na=False)]
        if t.empty:
            raise SystemExit("HGNC table has no rows with an HGNC:n id; the download is not the complete set.")
        if len(t) < 10000:
            print(f"  warning: HGNC table has only {len(t)} rows")
        r = cls(table=t)
        for row in t.itertuples(index=False):
            r.by_symbol[row.symbol.upper()] = row.hgnc_id
            for p in _split(row.prev_symbol):
                r.by_prev.setdefault(p.upper(), row.hgnc_id)
            for a in _split(row.alias_symbol):
                a = a.upper()
                if a in r.by_alias and r.by_alias[a] != row.hgnc_id:
                    r.ambiguous_alias.add(a)
                r.by_alias.setdefault(a, row.hgnc_id)
        return r

    def resolve(self, symbol: str | None, hgnc_id: str | None = None) -> tuple[str | None, str]:
        """Return (HGNC ID or None, how it was resolved)."""
        if hgnc_id:
            m = HGNC_ID_RE.match(str(hgnc_id).strip())
            if m:
                hid = f"HGNC:{m.group(1)}"
                if (self.table["hgnc_id"] == hid).any():
                    return hid, "id"
        if not symbol or not isinstance(symbol, str):
            return None, "no_symbol"
        s = symbol.strip().upper()
        if s in self.by_symbol:
            return self.by_symbol[s], "symbol"
        if s in self.by_prev:
            return self.by_prev[s], "prev_symbol"
        if s in self.by_alias and s not in self.ambiguous_alias:
            return self.by_alias[s], "alias"
        if s in self.ambiguous_alias:
            return None, "ambiguous_alias"
        return None, "unresolved"

    def current_symbol(self, hgnc_id: str) -> str:
        hit = self.table.loc[self.table["hgnc_id"] == hgnc_id, "symbol"]
        return hit.iloc[0] if len(hit) else ""

    def info(self, hgnc_id: str) -> dict:
        hit = self.table.loc[self.table["hgnc_id"] == hgnc_id]
        if not len(hit):
            return {}
        row = hit.iloc[0]
        return {
            "hgnc_id": hgnc_id,
            "symbol": row.get("symbol", ""),
            "gene_name": row.get("name", ""),
            "locus_group": row.get("locus_group", ""),
            "omim_id": row.get("omim_id", ""),
            "ensembl_gene_id": row.get("ensembl_gene_id", ""),
        }


def _read_hgnc(path) -> pd.DataFrame:
    """Read the HGNC complete set, tolerating a BOM, header casing, and either quoting style."""
    import csv
    best = None
    for quoting in (csv.QUOTE_NONE, csv.QUOTE_MINIMAL):
        try:
            t = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, quoting=quoting,
                            encoding="utf-8-sig", on_bad_lines="skip")
        except Exception as e:  # try the other quoting style
            print(f"  HGNC read with quoting={quoting} failed: {e}")
            continue
        t.columns = [str(c).strip().strip('"').lower() for c in t.columns]
        if "hgnc_id" in t.columns and (best is None or len(t) > len(best)):
            best = t
    if best is None:
        raise SystemExit(f"Could not read HGNC file at {path}")
    return best


def _split(cell: str) -> list[str]:
    if not cell:
        return []
    return [x.strip() for x in str(cell).split("|") if x.strip()]
