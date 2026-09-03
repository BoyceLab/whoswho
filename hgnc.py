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
        t = pd.read_csv(path, sep="\t", dtype=str, low_memory=False, quoting=3)
        keep = ["hgnc_id", "symbol", "name", "locus_group", "status",
                "prev_symbol", "alias_symbol", "omim_id", "ensembl_gene_id"]
        keep = [c for c in keep if c in t.columns]
        t = t[keep].fillna("")
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


def _split(cell: str) -> list[str]:
    if not cell:
        return []
    return [x.strip() for x in str(cell).split("|") if x.strip()]
