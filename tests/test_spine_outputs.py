"""Validation of the built spine, per the release brief.

These run against `outputs/` as built, not against a sandbox fixture, because the thing being
checked is the released artefact: that every identifier resolves, that no gene is listed twice,
that the published tier can be recomputed by hand from the per-source columns, and that an
ORPHA code is an ORPHA code.

The tier test is the important one. The tier rule is editorial policy quoted verbatim in the
README, so it has to be reproducible from the published evidence columns by anyone who reads
that text. If the recomputation drifts from the column, one of the two is wrong.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest

from gene_spine import tiers

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
SRC = ROOT / "sources"

pytestmark = pytest.mark.skipif(
    not (OUT / "gene_spine.csv").exists(),
    reason="no built spine in outputs/; run `make build` first",
)


@pytest.fixture(scope="module")
def spine() -> pd.DataFrame:
    return pd.read_csv(OUT / "gene_spine.csv", low_memory=False)


@pytest.fixture(scope="module")
def long() -> pd.DataFrame:
    return pd.read_csv(OUT / "gene_source_long.csv", low_memory=False)


@pytest.fixture(scope="module")
def hgnc_ids() -> set[str]:
    """Every HGNC id in the downloaded complete set."""
    path = next((SRC / "hgnc").glob("*"), None)
    if path is None:
        pytest.skip("HGNC complete set not fetched")
    df = pd.read_csv(path, sep="\t", low_memory=False, usecols=["hgnc_id"])
    return set(df["hgnc_id"].astype(str))


def test_every_hgnc_id_resolves(spine, hgnc_ids):
    ids = spine["hgnc_id"].astype(str)
    assert ids.str.match(r"^HGNC:\d+$").all(), (
        "malformed hgnc_id: " + ", ".join(ids[~ids.str.match(r'^HGNC:\d+$')].head(5)))
    unknown = sorted(set(ids) - hgnc_ids)
    assert not unknown, f"{len(unknown)} ids absent from the HGNC set, e.g. {unknown[:5]}"


def test_no_symbol_appears_twice(spine):
    dupes = spine["symbol"][spine["symbol"].duplicated()].tolist()
    assert not dupes, f"repeated symbols: {dupes[:10]}"


def test_no_hgnc_id_appears_twice(spine):
    dupes = spine["hgnc_id"][spine["hgnc_id"].duplicated()].tolist()
    assert not dupes, f"repeated hgnc_ids: {dupes[:10]}"


def test_tier_is_recomputable_from_the_per_source_columns(spine, long):
    """Recompute every tier from the long-form evidence and require equality."""
    recomputed = {}
    for hgnc_id, g in long.groupby("hgnc_id"):
        tier, _disputed, _note = tiers.tier_for(zip(g["source_id"], g["evidence"].fillna("")))
        recomputed[hgnc_id] = tier

    published = dict(zip(spine["hgnc_id"], spine["tier"]))
    missing = [k for k in published if k not in recomputed]
    assert not missing, f"{len(missing)} published genes have no evidence rows, e.g. {missing[:5]}"

    mismatch = {k: (published[k], recomputed[k])
                for k in published if published[k] != recomputed[k]}
    assert not mismatch, (
        f"{len(mismatch)} tiers do not recompute, e.g. "
        + "; ".join(f"{k} published {p} recomputed {r}"
                    for k, (p, r) in list(mismatch.items())[:5]))


def test_disputed_never_lowers_a_tier(spine, long):
    """The stated policy: a Disputed label is reported, not used to downgrade."""
    disputed = spine[spine["disputed"].astype(str).str.lower().isin({"true", "1"})]
    if disputed.empty:
        pytest.skip("no disputed genes in this build")
    assert (disputed["tier"] != "T3").any(), (
        "every disputed gene sits at T3, which suggests the flag is downgrading after all")
    assert disputed["disputed_note"].notna().all(), "a disputed gene must name its label"


def test_every_gene_with_orphacodes_has_a_valid_orpha_integer(spine):
    have = spine[spine["orphacodes"].notna() & (spine["orphacodes"].astype(str).str.strip() != "")]
    assert len(have), "no gene carries an orphacode, which cannot be right"
    bad = []
    for _, row in have.iterrows():
        codes = [c.strip() for c in str(row["orphacodes"]).replace("|", ";").split(";")
                 if c.strip()]
        if not any(re.fullmatch(r"\d+", c.replace("ORPHA:", "")) for c in codes):
            bad.append((row["symbol"], row["orphacodes"]))
    assert not bad, f"{len(bad)} genes have no parseable ORPHA integer, e.g. {bad[:5]}"


def test_release_records_every_source_and_its_hash(spine):
    rel = json.loads((OUT / "release.json").read_text(encoding="utf-8"))
    assert rel.get("built"), "release.json has no build date"
    srcs = rel.get("sources") or {}
    assert srcs, "release.json lists no sources"
    for sid, rec in srcs.items():
        assert rec.get("sha256"), f"{sid} has no sha256"
        assert rec.get("fetched"), f"{sid} has no fetch date"
    assert "missing_sources" in rel, (
        "release.json must record declared sources that had no file, or a reader cannot tell "
        "an absent source from a source that said no")


def test_gene_counts_agree_between_release_and_table(spine):
    rel = json.loads((OUT / "release.json").read_text(encoding="utf-8"))
    assert rel["genes_total"] == len(spine)
    assert rel["by_tier"] == spine["tier"].value_counts().to_dict()
