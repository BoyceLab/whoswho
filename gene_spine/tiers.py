"""Derived columns: evidence tier and phenotype domain.

These rules are the editorial policy of the spine. They are deliberately
simple and stated in the README verbatim, so a tier can be recomputed by
hand from the per-source columns in the output.

Tier (strength of the gene-disease evidence recorded by the sources):
  T1  An expert panel curation says definitive or strong: ClinGen Definitive or
      Strong, SysNDD Definitive, or Gene2Phenotype definitive or strong.
      Also T1 when two or more independent clinical panels rate the gene Green,
      since agreement across panels is the point, not any single listing.
  T2  A single PanelApp Green, ClinGen Moderate, SysNDD Moderate, G2P moderate,
      SFARI 1 or 2, PanelApp Amber, Genes4Epilepsy inclusion, or a Simons
      Searchlight studied gene.
  T3  Present on at least one list and none of the above.

  A Refuted or Disputed label is recorded in the `disputed` flag and named in
  `disputed_note`. It does NOT lower the tier, because a classification applies
  to one gene-disease relationship: a gene can be definitive for one condition
  and disputed for another, and collapsing that into one score misstates both.

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
NDD_SOURCES = {"sysndd", "sfari", "g2p_dd", "panelapp_au_id", "panelapp_ge_id", "simons_searchlight", "simons_gene_list"}

# MCD is how Genes4Epilepsy writes malformation of cortical development. The full words below
# never match its abbreviated vocabulary, so without MCD the systemic branch could not fire at
# all. G4E_NDD is tested first, so "DEE, MCD" still reads as ndd_with_epilepsy; only an
# MCD-without-DEE phenotype reaches here. PME stays epilepsy by decision.
G4E_SYSTEMIC = re.compile(r"(MCD|malformation|metabolic|mitochondrial|storage|syndrom)", re.I)
G4E_NDD = re.compile(r"(DEE|encephalopath|developmental)", re.I)


def _family(source_id: str) -> str:
    if source_id.startswith("panelapp"):
        return "panelapp"
    return source_id


def tier_for(rows) -> tuple[str, bool, str]:
    """rows: iterable of (source_id, evidence) for one gene.

    Returns (tier, disputed, disputed_note). Disputed is reported, never used to
    downgrade, since the label belongs to one gene-disease relationship.
    """
    disputed_labels = []
    best = 3
    green_panels = set()
    on_any = False
    for sid, ev in rows:
        on_any = True
        ev = (ev or "").strip()
        fam = _family(sid)
        for label in re.split(r"[;,]", ev):
            label = label.strip()
            if DISPUTED.search(label):
                disputed_labels.append(f"{sid}: {label}")
                continue
            if fam in STRONG and STRONG[fam].match(label):
                if fam == "panelapp":
                    green_panels.add(sid)      # one green is not agreement
                else:
                    best = min(best, 1)        # an expert panel curation is
            elif fam in MODERATE and MODERATE[fam].match(label):
                best = min(best, 2)
            elif sid in ("genes4epilepsy", "simons_searchlight"):
                best = min(best, 2)
    if len(green_panels) >= 2:
        best = min(best, 1)
    elif green_panels:
        best = min(best, 2)
    if not on_any:
        return "", False, ""
    return f"T{best}", bool(disputed_labels), "; ".join(sorted(set(disputed_labels)))
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
