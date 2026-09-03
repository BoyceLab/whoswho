# Gene Spine


A harmonized, HGNC-keyed list of genes associated with epilepsy and
neurodevelopmental disorders, with a curated organization layer joined to it.
Standalone; no dependency on any institutional platform.

The spine has two halves that are built and versioned differently:

1. The gene half is mechanical. Each source list is parsed, every symbol is
   resolved to an HGNC ID, and the lists are joined into one row per gene with
   the source's own evidence label kept verbatim. A tier and a domain are
   derived from those labels by a rule stated in `gene_spine/tiers.py`.
2. The organization half is curated. `curated/organizations.csv` is the
   authoritative table. Membership rosters and scraped pages come in only as
   candidate rows, and a candidate is dropped the moment a curated row exists
   for the same gene and website.

## Sources

Declared in `config/sources.yaml`. Each entry records the URL or landing page,
the release version where the source has one, the license, and the column
names the parser should look for. The parser prints the real header on a
mismatch so the fix is a config edit.

| id | list | scope | fetch |
|---|---|---|---|
| hgnc | HGNC complete set | identifier backbone | automatic |
| genes4epilepsy | Genes4Epilepsy v2026-03 | epilepsy | automatic |
| clingen | ClinGen gene-disease validity | both | automatic |
| panelapp_au_epilepsy, panelapp_au_id | PanelApp Australia | epilepsy, ID | automatic (panel found by name) |
| panelapp_ge_epilepsy, panelapp_ge_id | Genomics England PanelApp | epilepsy, ID | automatic (panel found by name) |
| simons_searchlight | Simons Searchlight gene list PDF | NDD | automatic |
| orphadata_genes | Orphadata product 6 | gene to ORPHAcode bridge | automatic |
| wang_omim_2023 | Wang et al. 2023 OMIM/HGMD supplement | both (three categories) | manual |
| sagas | SAGAS | epilepsy, cross-species | manual |
| sysndd | SysNDD | NDD | manual export |
| sfari | SFARI Gene human genes | autism/NDD | manual export |
| g2p_dd | Gene2Phenotype DD panel | NDD | manual download |
| ren_members, combinedbrain_members, epilepsy_foundation_genetic, simons_support_links | organization rosters | org candidates | manual CSV |

Orphanet does not publish patient organizations as a bulk file, so the
ORPHAcode is carried on every gene and organization row as the join key, and
the organizations themselves enter through the curated table.

## Running it

```
pip install -r requirements.txt
python -m gene_spine.fetch        # downloads automatic sources, lists what is still manual
python -m gene_spine.build        # writes outputs/
python -m pytest tests            # offline fixture test of the whole pipeline
```

## Outputs

- `gene_spine.csv`, `gene_spine.json`: one row per HGNC gene. Lead columns are
  `hgnc_id, symbol, gene_name, tier, domain, disputed, n_sources, sources,
  orphacodes, omim_id, ensembl_gene_id`, followed by `on_<source>`,
  `ev_<source>` (the source's own label), and `pheno_<source>` columns.
- `gene_source_long.csv`: one row per gene per source, with how the symbol was
  resolved (id, symbol, prev_symbol, alias).
- `unresolved_symbols.csv`: anything that did not map to HGNC. Empty is the goal.
- `cnv_regions.csv`: CNV entries that are not genes.
- `gene_orpha.csv`: gene to ORPHAcode bridge.
- `org_layer.csv`: curated rows plus surviving candidates, keyed to HGNC and ORPHAcode.
- `org_gaps.csv`: genes with no organization row at all, with tier and domain,
  in the order you would check them.
- `release.json`: build date, the manifest of every source file with its fetch
  date and sha256, the tier rule text, and row counts by tier and domain.

## Tier rule

T1: ClinGen Definitive or Strong, or any PanelApp Green, or SysNDD Definitive,
or G2P definitive/strong.
T2: Genes4Epilepsy inclusion, SFARI 1 or 2, SysNDD Moderate, G2P moderate,
PanelApp Amber, ClinGen Moderate, Simons Searchlight studied gene.
T3: on at least one list and none of the above.
A Refuted or Disputed label from ClinGen or SysNDD pins the gene at T3 with
`disputed = true` regardless of other lists.

## Domain rule

`epilepsy_primary`, `ndd_with_epilepsy`, `ndd_no_documented_epilepsy`,
`systemic_with_seizures`, `epilepsy_only_unclassified`. Wang's three
categories take precedence where present; otherwise the domain follows which
scope of sources the gene appears on, with Genes4Epilepsy phenotype text
deciding between DEE and malformation/metabolic groupings.

## Release discipline

A release is the `outputs/` directory plus `sources/manifest.json`. Every
number in a release is reproducible from the manifest: the same source
snapshots and the same tier rule give the same spine. Bump the Genes4Epilepsy
version in the config when its half-yearly update lands, re-fetch, rebuild,
and diff `gene_spine.csv` against the previous release.

## Layout

```
config/sources.yaml     source registry
gene_spine/hgnc.py      HGNC resolver
gene_spine/parsers.py   one parser per source layout
gene_spine/tiers.py     tier and domain rules
gene_spine/fetch.py     downloader and manifest
gene_spine/build.py     join and outputs
curated/                organizations.csv (authoritative) and its column guide
schema/                 JSON schema for a spine row
tests/                  offline fixture test
```
