"""Derived columns: evidence tier and phenotype domain.

These rules are the editorial policy of the spine. They are deliberately
simple and stated in the README verbatim, so a tier can be recomputed by
hand from the per-source columns in the output.

Tier (strength of gene-disease evidence, any domain):
  T1  ClinGen Definitive or Strong
      OR any PanelApp Green
      OR SysNDD Definitive
      OR G2P definitive / strong
  T2  Genes4Epilepsy (its inclusion is itself a manual curation)
      OR SFARI score 1 or 2
      OR SysNDD Moderate  OR G2P moderate  OR PanelApp Amber
      OR ClinGen Moderate
      OR Simons Searchlight studied gene
  T3  everything else that appears on at least one list
      (Wang OMIM-derived, SAGAS, PanelApp Red, ClinGen Limited, SFARI 3, SysNDD Limited)
  Refuted/Disputed labels from ClinGen or SysNDD are recorded in a flag column
  and the gene is held at T3 with disputed=True regardless of other lists.

Domain (what phenotype the evidence speaks to):
  epilepsy_primary          epilepsy source says the epilepsy is the core phenotype
  ndd_with_epilepsy         gene appears on both an epilepsy source and an NDD source,
                            or Wang labels it NDD-with-epilepsy
  ndd_no_documented_epilepsy gene appears on NDD sources only
  systemic_with_seizures    Wang systemic category, or Genes4Epilepsy phenotype in its
                            malformation/metabolic groupings and no NDD source
  epilepsy_only_unclassified epilepsy sources only, no category text to go on
"""
from __future__ import annotations

import re

STRONG = {
    "clingen": re.compile(r"^(definitive|strong)$", re.I),
    "panelapp": re.compile(r"^green$", re.I),
    "sysndd": re.compile(r"^definitive$", re.I),
    "g2p_dd": re.compile(r"(definitive|strong)", re.I),
}
MODERATE = {
    "clingen": re.compile(r"^moderate$", re.I),
    "panelapp": re.compile(r"^amber$", re.I),
    "sysndd": re.compile(r"^moderate$", re.I),
    "g2p_dd": re.compile(r"moderate", re.I),
    "sfari": re.compile(r"^(1|2)(\.0)?$"),
}
DISPUTED = re.compile(r"(refuted|disputed)", re.I)

EPILEPSY_SOURCES = {"genes4epilepsy", "sagas", "panelapp_au_epilepsy", "panelapp_ge_epilepsy"}
NDD_SOURCES = {"sysndd", "sfari", "g2p_dd", "panelapp_au_id", "panelapp_ge_id", "simons_searchlight"}

G4E_SYSTEMIC = re.compile(r"(malformation|metabolic|mitochondrial|storage|syndrom)", re.I)
G4E_NDD = re.compile(r"(DEE|encephalopath|developmental)", re.I)


def _family(source_id: str) -> str:
    if source_id.startswith("panelapp"):
        return "panelapp"
    return source_id


def tier_for(rows) -> tuple[str, bool]:
    """rows: iterable of (source_id, evidence) for one gene."""
    disputed = False
    best = 3
    on_any = False
    for sid, ev in rows:
        on_any = True
        ev = (ev or "").strip()
        if DISPUTED.search(ev):
            disputed = True
            continue
        fam = _family(sid)
        if fam in STRONG and STRONG[fam].match(ev):
            best = min(best, 1)
        elif fam in MODERATE and MODERATE[fam].match(ev):
            best = min(best, 2)
        elif sid in ("genes4epilepsy", "simons_searchlight"):
            best = min(best, 2)
        elif sid == "clingen" and ev:
            best = min(best, 3)
    if not on_any:
        return "", disputed
    if disputed:
        return "T3", True
    return f"T{best}", False


def domain_for(rows, clingen_epi_terms, clingen_ndd_terms) -> str:
    """rows: iterable of (source_id, evidence, phenotype, extra_dict)."""
    epi = ndd = False
    wang = None
    g4e_pheno = ""
    for sid, ev, ph, extra in rows:
        if sid in EPILEPSY_SOURCES:
            epi = True
        if sid in NDD_SOURCES:
            ndd = True
        if sid == "genes4epilepsy":
            g4e_pheno = ph or ""
        if sid == "wang_omim_2023":
            wang = (extra or {}).get("domain") or ev
        if sid == "clingen":
            text = (ph or "").lower()
            if any(t in text for t in clingen_epi_terms):
                epi = True
            if any(t in text for t in clingen_ndd_terms):
                ndd = True
    if wang in ("epilepsy_primary", "ndd_with_epilepsy", "systemic_with_seizures"):
        if wang == "epilepsy_primary" and ndd:
            return "ndd_with_epilepsy"
        return wang
    if epi and ndd:
        return "ndd_with_epilepsy"
    if ndd and not epi:
        return "ndd_no_documented_epilepsy"
    if epi and G4E_NDD.search(g4e_pheno):
        return "ndd_with_epilepsy"
    if epi and G4E_SYSTEMIC.search(g4e_pheno):
        return "systemic_with_seizures"
    if epi and g4e_pheno:
        return "epilepsy_primary"
    if epi:
        return "epilepsy_only_unclassified"
    return "unclassified"
