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


def test_every_source_states_what_its_hash_covers():
    """A hash over a canonicalised node is not a hash of the file, so it has to say so.

    SysNDD's response carries its own query timing, which would make the recorded hash change
    between two fetches of identical data.
    """
    rel = json.loads((OUT / "release.json").read_text(encoding="utf-8"))
    for sid, rec in rel["sources"].items():
        assert rec.get("sha256_scope"), f"{sid} does not say what its hash covers"
    sysndd = rel["sources"].get("sysndd")
    if sysndd:
        assert "canonicalised" in sysndd["sha256_scope"]


# --------------------------------------------------------- the derived lists

@pytest.fixture(scope="module")
def gene_list() -> pd.DataFrame:
    p = OUT / "epilepsy_ndd_gene_list.csv"
    if not p.exists():
        pytest.skip("no gene list; run `make lists` after `make build`")
    return pd.read_csv(p, comment="#", low_memory=False)


def test_gene_list_membership_values_are_the_four_branches(gene_list):
    assert set(gene_list["list_membership"]) <= {"epilepsy", "both", "ndd_only",
                                                 "epilepsy_other"}
    assert set(gene_list["list_membership"]) >= {"both", "ndd_only"}


def test_gene_list_is_a_subset_of_the_spine(gene_list, spine):
    extra = set(gene_list["hgnc_id"]) - set(spine["hgnc_id"])
    assert not extra, f"gene list holds ids absent from the spine: {sorted(extra)[:5]}"
    assert len(gene_list) == gene_list["hgnc_id"].nunique(), "a gene is listed twice"


def test_only_genes4epilepsy_members_carry_its_phenotype(gene_list):
    """The column is that source's own text, so it is blank for a gene the source does not list."""
    off_list = gene_list[~gene_list["list_membership"].isin(["epilepsy", "both"])]
    assert off_list["genes4epilepsy_phenotype"].fillna("").eq("").all()
    on_list = gene_list[gene_list["list_membership"].isin(["epilepsy", "both"])]
    assert on_list["genes4epilepsy_phenotype"].fillna("").ne("").all()


def test_org_list_is_candidate_throughout():
    p = OUT / "epilepsy_orgs.csv"
    if not p.exists():
        pytest.skip("no org list; run `make lists`")
    orgs = pd.read_csv(p, comment="#", low_memory=False)
    assert (orgs["status"] == "candidate").all(), (
        "every row is a candidate until a curator reviews it")
    assert orgs["org_name"].notna().all(), "an organisation row with no name is not usable"
    # One row per organisation and gene.
    assert not orgs.duplicated(subset=["org_name", "symbol"]).any()


def test_org_counts_in_the_gene_list_match_the_org_list(gene_list):
    """Every row, not a sample: the two files are written from one frame and must agree."""
    p = OUT / "epilepsy_orgs.csv"
    if not p.exists():
        pytest.skip("no org list; run `make lists`")
    orgs = pd.read_csv(p, comment="#", low_memory=False)
    counts = orgs.groupby(orgs["symbol"].str.upper())["org_name"].nunique()
    wrong = []
    for _, row in gene_list.iterrows():
        expected = int(counts.get(str(row["symbol"]).upper(), 0))
        if int(row["org_count"]) != expected:
            wrong.append(f"{row['symbol']}: list {row['org_count']}, org file {expected}")
    assert not wrong, f"{len(wrong)} disagreements, e.g. {wrong[:5]}"


def test_org_names_match_the_org_list(gene_list):
    p = OUT / "epilepsy_orgs.csv"
    if not p.exists():
        pytest.skip("no org list; run `make lists`")
    orgs = pd.read_csv(p, comment="#", low_memory=False)
    joined = (orgs.groupby(orgs["symbol"].str.upper())["org_name"]
              .apply(lambda s: "; ".join(sorted({str(v) for v in s.dropna()}))))
    listed = gene_list[gene_list["org_count"] > 0]
    for _, row in listed.iterrows():
        assert row["org_names"] == joined.get(str(row["symbol"]).upper(), ""), (
            f"{row['symbol']}: names differ between the two files")


def test_the_emitter_is_deterministic(tmp_path):
    """Run it twice and require identical bytes.

    Two runs that differ mean something in the pipeline depends on dictionary or filesystem
    order, which would show up as a spurious diff on every scheduled build and would make the
    hashes published in the README wrong at random.
    """
    from gene_spine import epilepsy_list

    if not (OUT / "gene_spine.csv").exists():
        pytest.skip("no built spine")
    targets = [OUT / "epilepsy_ndd_gene_list.csv", OUT / "epilepsy_orgs.csv"]

    epilepsy_list.main()
    first = {p.name: p.read_bytes() for p in targets}
    epilepsy_list.main()
    second = {p.name: p.read_bytes() for p in targets}

    for name in first:
        assert first[name] == second[name], f"{name} differs between two runs of the emitter"
