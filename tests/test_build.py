"""Runs the full build against tiny fixtures so the join, tiering, and org
layer can be checked without network access. Fixtures mimic each source's
real layout (ClinGen preamble, PanelApp JSON, Orphadata XML, HGNC pipes)."""
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest
import yaml

import gene_spine.build as build

ROOT = Path(__file__).resolve().parents[1]

HGNC = """hgnc_id\tsymbol\tname\tlocus_group\tstatus\tprev_symbol\talias_symbol\tomim_id\tensembl_gene_id
HGNC:10588\tSCN2A\tsodium voltage-gated channel alpha subunit 2\tprotein-coding gene\tApproved\tSCN2A1|SCN2A2\tNAC2|BFIS3\t182390\tENSG00000136531
HGNC:11444\tSTXBP1\tsyntaxin binding protein 1\tprotein-coding gene\tApproved\t\tMUNC18-1|UNC18\t602926\tENSG00000136854
HGNC:11497\tSYNGAP1\tsynaptic Ras GTPase activating protein 1\tprotein-coding gene\tApproved\t\tMRD5|RASA1\t603384\tENSG00000197283
HGNC:6294\tKCNT1\tpotassium sodium-activated channel subfamily T member 1\tprotein-coding gene\tApproved\t\tSLACK|KCa4.1\t608167\tENSG00000107147
HGNC:9948\tRAI1\tretinoic acid induced 1\tprotein-coding gene\tApproved\t\tSMCR\t607642\tENSG00000108557
HGNC:2231\tCOL4A1\tcollagen type IV alpha 1 chain\tprotein-coding gene\tApproved\t\t\t120130\tENSG00000187498
"""

G4E = """Gene\tHGNC_ID\tEpilepsy_phenotype\tInheritance
SCN2A\tHGNC:10588\tDEE\tAD
STXBP1\tHGNC:11444\tDEE\tAD
KCNT1\tHGNC:6294\tDEE; Focal\tAD
COL4A1\tHGNC:2231\tMalformation\tAD
SCN2A1\t\tDEE\tAD
"""

CLINGEN = """CLINGEN GENE VALIDITY CURATIONS
FILE CREATED: 2026-09-01
WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity
+++++++++++++
GENE SYMBOL,GENE ID (HGNC),DISEASE LABEL,DISEASE ID (MONDO),MOI,SOP,CLASSIFICATION,ONLINE REPORT,CLASSIFICATION DATE,GCEP
+++++++++++++
SCN2A,HGNC:10588,developmental and epileptic encephalopathy 11,MONDO:0013388,AD,SOP7,Definitive,url,2020-01-01,Epilepsy
SYNGAP1,HGNC:11497,intellectual disability autosomal dominant 5,MONDO:0013250,AD,SOP7,Definitive,url,2020-01-01,ID/Autism
RAI1,HGNC:9948,Smith-Magenis syndrome,MONDO:0008525,AD,SOP7,Moderate,url,2020-01-01,ID/Autism
KCNT1,HGNC:6294,some phenotype,MONDO:0000001,AD,SOP7,Disputed,url,2020-01-01,Epilepsy
"""

PANELAPP = {
    "id": 202, "name": "Genetic Epilepsy", "version": "1.20",
    "genes": [
        {"gene_data": {"gene_symbol": "STXBP1", "hgnc_id": "HGNC:11444"}, "confidence_level": "3",
         "phenotypes": ["Developmental and epileptic encephalopathy 4"], "mode_of_inheritance": "MONOALLELIC"},
        {"gene_data": {"gene_symbol": "COL4A1", "hgnc_id": "HGNC:2231"}, "confidence_level": "2",
         "phenotypes": ["Porencephaly"], "mode_of_inheritance": "MONOALLELIC"},
    ],
}

SFARI = """status,gene-symbol,gene-name,gene-score,syndromic
9,SYNGAP1,synaptic Ras GTPase activating protein 1,1,1
9,SCN2A,sodium channel 2,1,1
9,RAI1,retinoic acid induced 1,3,1
"""

ORPHA = """<?xml version="1.0"?>
<JDBOR><DisorderList>
<Disorder><OrphaCode>1934</OrphaCode><Name lang="en">Early infantile epileptic encephalopathy</Name>
<DisorderGeneAssociationList>
<DisorderGeneAssociation><Gene><Symbol>STXBP1</Symbol><ExternalReferenceList>
<ExternalReference><Source>HGNC</Source><Reference>11444</Reference></ExternalReference></ExternalReferenceList></Gene>
<DisorderGeneAssociationType><Name lang="en">Disease-causing germline mutation(s) in</Name></DisorderGeneAssociationType>
<DisorderGeneAssociationStatus><Name lang="en">Assessed</Name></DisorderGeneAssociationStatus>
</DisorderGeneAssociation></DisorderGeneAssociationList></Disorder>
</DisorderList></JDBOR>
"""

REN = """Organization,Gene,Website
STXBP1 Foundation,STXBP1,https://www.stxbp1disorders.org
SynGAP Research Fund,SYNGAP1,https://curesyngap1.org
FamilieSCN2A Foundation,SCN2A,https://www.scn2a.org
"""

CURATED = """org_id,org_name,url,hgnc_id,symbol,orphacode,country,org_type,registry_platform,registry_url,natural_history_study,data_elements_documented,irb_status,maturity_stage,coalitions,contact,status,source,reviewed_by,review_date,notes
stxbp1-foundation,STXBP1 Foundation,https://stxbp1disorders.org,HGNC:11444,STXBP1,,US,foundation,Ciitizen,,yes,yes,approved,4,REN;COMBINEDBrain,,curated,hand,DB,2026-08-01,
"""


@pytest.fixture
def sandbox(tmp_path, monkeypatch):
    src, out, cur = tmp_path / "sources", tmp_path / "outputs", tmp_path / "curated"
    for d in ("hgnc", "genes4epilepsy", "clingen", "panelapp_au_epilepsy", "sfari", "orphadata_genes", "ren_members"):
        (src / d).mkdir(parents=True)
    cur.mkdir()
    (src / "hgnc" / "hgnc_complete_set.txt").write_text(HGNC)
    (src / "genes4epilepsy" / "EpilepsyGenes_v2026-03.tsv").write_text(G4E)
    (src / "clingen" / "clingen_gene_validity.csv").write_text(CLINGEN)
    (src / "panelapp_au_epilepsy" / "panelapp_au_epilepsy.json").write_text(json.dumps(PANELAPP))
    (src / "sfari" / "SFARI-Gene_genes.csv").write_text(SFARI)
    (src / "orphadata_genes" / "en_product6.xml").write_text(ORPHA)
    (src / "ren_members" / "ren_members.csv").write_text(REN)
    (cur / "organizations.csv").write_text(CURATED)
    shutil.copytree(ROOT / "config", tmp_path / "config")
    monkeypatch.setattr(build, "ROOT", tmp_path)
    monkeypatch.setattr(build, "SRC", src)
    monkeypatch.setattr(build, "OUT", out)
    monkeypatch.setattr(build, "CUR", cur)
    return tmp_path


def test_end_to_end(sandbox):
    build.main()
    out = sandbox / "outputs"
    spine = pd.read_csv(out / "gene_spine.csv").set_index("symbol")

    # STXBP1: ClinGen absent, one PanelApp Green, G4E, Orphadata -> T2, epilepsy with DEE text
    # -> ndd_with_epilepsy. T2 and not T1 because the tier rule was revised: a single Green is
    # T2, and T1 needs an expert-panel definitive or strong call, or two or more independent
    # panels Green. This fixture supplies only panelapp_au_epilepsy, so one Green is all there
    # is. The assertion said T1 until the rule changed under it.
    assert spine.loc["STXBP1", "tier"] == "T2"
    assert spine.loc["STXBP1", "domain"] == "ndd_with_epilepsy"
    assert str(spine.loc["STXBP1", "orphacodes"]).startswith("1934")

    # SCN2A: on G4E + ClinGen Definitive + SFARI 1 -> T1, both domains -> ndd_with_epilepsy
    assert spine.loc["SCN2A", "tier"] == "T1"
    assert spine.loc["SCN2A", "domain"] == "ndd_with_epilepsy"

    # SYNGAP1: ClinGen Definitive + SFARI 1, no epilepsy source -> T1, ndd_no_documented_epilepsy
    assert spine.loc["SYNGAP1", "tier"] == "T1"
    assert spine.loc["SYNGAP1", "domain"] == "ndd_no_documented_epilepsy"

    # RAI1: ClinGen Moderate + SFARI 3 -> T2, NDD only
    assert spine.loc["RAI1", "tier"] == "T2"

    # KCNT1: ClinGen Disputed plus G4E -> T2 with the disputed flag set. The rule was revised:
    # a Disputed label is reported and never downgrades, because the classification applies to
    # one gene-disease relationship and a gene can be definitive for one condition and disputed
    # for another. This assertion read T3 under the earlier rule.
    assert spine.loc["KCNT1", "tier"] == "T2"
    assert bool(spine.loc["KCNT1", "disputed"]) is True
    assert str(spine.loc["KCNT1", "disputed_note"]).strip() != "", (
        "a disputed gene has to name the label that made it disputed")

    # COL4A1: PanelApp Amber + G4E malformation -> T2, systemic_with_seizures
    assert spine.loc["COL4A1", "tier"] == "T2"
    assert spine.loc["COL4A1", "domain"] == "systemic_with_seizures"

    # SCN2A1 (previous symbol) resolves to SCN2A rather than creating a new row
    assert "SCN2A1" not in spine.index
    long = pd.read_csv(out / "gene_source_long.csv")
    assert (long["resolved_by"] == "prev_symbol").sum() == 1
    assert pd.read_csv(out / "unresolved_symbols.csv").empty

    org = pd.read_csv(out / "org_layer.csv", dtype=str).fillna("")
    # curated STXBP1 row kept, REN STXBP1 candidate dropped (same host), others become candidates
    stx = org[org["symbol"] == "STXBP1"]
    assert len(stx) == 1 and stx.iloc[0]["status"] == "curated"
    assert set(org.loc[org["status"] == "candidate", "symbol"]) == {"SYNGAP1", "SCN2A"}
    # orphacode filled from bridge on the curated row
    assert str(stx.iloc[0]["orphacode"]) == "1934"

    gaps = pd.read_csv(out / "org_gaps.csv")
    assert set(gaps["symbol"]) == {"KCNT1", "RAI1", "COL4A1"}

    rel = json.loads((out / "release.json").read_text())
    assert rel["genes_total"] == 6
