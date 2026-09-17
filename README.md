# whoswho: a gene spine for rare epilepsy and neurodevelopmental disorders

One table of genes, assembled from the public curation sources, with an evidence tier, a
phenotype domain, and the organisations that serve each gene. The build is deterministic: given
the same source snapshots it produces the same outputs, and every snapshot's URL, fetch date and
SHA-256 is recorded in `outputs/release.json`.

This file documents the project. The curated organisation table has its own notes in
[curated/README.md](curated/README.md).

## Published files

Served from GitHub Pages at the repository root, so a consumer can fetch them by URL without
cloning. Replace `<pages-base>` with the Pages address once Pages is enabled; it is left blank
here rather than guessed.

| path | bytes | sha256 |
|---|---|---|
| `<pages-base>/outputs/epilepsy_ndd_gene_list.csv` | 874,057 | `e8f8c8148d435d64e0bf1ef73d2158ed123658cd139cd65cebc9e2def669e456` |
| `<pages-base>/outputs/epilepsy_orgs.csv` | 55,294 | `4349b56c524c67962ad9b94e6df569a450980d26068146f5dba8c9ba94e14edf` |
| `<pages-base>/outputs/gene_spine.json` | 9,550,535 | `99af53934d055d63ef978f4b05e6c69f72f798d9d519429404c71cdf7f3fb9f6` |
| `<pages-base>/outputs/gene_spine.csv` | 2,799,406 | `416ac7a1e76b4a945bc943401c0b79dc81d46fe783822a445a97ea8a5073cb75` |
| `<pages-base>/outputs/release.json` | 8,744 | see the file itself |

Those hashes are the 2026-09-16 build. **The build rewrites them, so treat the table as a
snapshot and read `outputs/release.json` for the current values**; it carries a hash per output
and per source. The paths are stable, the hashes are not.

## The two gene lists

`epilepsy_ndd_gene_list.csv` is the union of three branches, with `list_membership` naming which
one a gene came from:

| membership | rule | count |
|---|---|---:|
| `both` | on Genes4Epilepsy and on an NDD or ID source | 1,059 |
| `epilepsy` | on Genes4Epilepsy only | 19 |
| `ndd_only` | on an NDD or ID source, not on Genes4Epilepsy | 3,732 |
| `epilepsy_other` | on an epilepsy panel, or ClinGen or Orphadata naming epilepsy, seizures or encephalopathy, and in neither branch above | 31 |
| | **union** | **4,841** |

`tier` and `domain` travel with each row as description. They are not selection criteria for this
file: a gene is on the list because a source names it, not because it reached a tier.

`epilepsy_orgs.csv` is one row per organisation and gene, deduplicated on normalised website
host, every row `status = candidate` until reviewed. 527 rows over 420 organisations. The curated
layer in `curated/organizations.csv` is never written by the build.

## Tier and domain rules

Editorial policy, quoted verbatim from `gene_spine/tiers.py`, which is the only place they are
implemented.

Tier, the strength of the gene-disease evidence recorded by the sources:

- **T1** An expert panel curation says definitive or strong: ClinGen Definitive or Strong, SysNDD
  Definitive, or Gene2Phenotype definitive or strong. Also T1 when two or more independent
  clinical panels rate the gene Green, since agreement across panels is the point, not any single
  listing.
- **T2** A single PanelApp Green, ClinGen Moderate, SysNDD Moderate, G2P moderate, SFARI 1 or 2,
  PanelApp Amber, Genes4Epilepsy inclusion, or a Simons Searchlight studied gene.
- **T3** Present on at least one list and none of the above.

A Refuted or Disputed label is recorded in the `disputed` flag and named in `disputed_note`. It
does **not** lower the tier, because a classification applies to one gene-disease relationship: a
gene can be definitive for one condition and disputed for another, and collapsing that into one
score misstates both.

Domain, what phenotype the evidence speaks to:

- `epilepsy_primary` an epilepsy source says the epilepsy is the core phenotype
- `ndd_with_epilepsy` on both an epilepsy source and an NDD source
- `ndd_no_documented_epilepsy` on NDD sources only
- `systemic_with_seizures` a systemic category, or a Genes4Epilepsy malformation or metabolic
  phenotype with no NDD source
- `epilepsy_only_unclassified` epilepsy sources only, with no category text to go on

Current distribution over 4,842 genes: `ndd_no_documented_epilepsy` 3,497,
`ndd_with_epilepsy` 1,302, `epilepsy_only_unclassified` 30, `epilepsy_primary` 9,
`systemic_with_seizures` 4.

## Running it

```
pip install -r requirements.txt
make fetch     # python -m gene_spine.fetch, downloads the automatic sources
make build     # python -m gene_spine.build, writes outputs/
make test      # python -m pytest -q tests
```

`make fetch` prints any declared source whose file is absent, and the build continues without it.
An absent source's evidence columns are **missing from the outputs rather than False**, so a
reader cannot mistake "we did not look" for "the source said no"; `outputs/release.json` records
the absences under `missing_sources`.

## Sources

Automatic, downloaded by `make fetch`: HGNC, Genes4Epilepsy, ClinGen, PanelApp Australia
(Genetic Epilepsy, Intellectual disability), PanelApp Genomics England (Early onset or syndromic
epilepsy, Intellectual disability), SysNDD, Gene2Phenotype DD, Orphadata.

By hand, placed under `sources/<id>/`: SFARI Gene (the site requires its download button), the
Wang 2023 supplement, SAGAS, and the organisation and registry files. Each hand export records
its original filename, release date and export date in the manifest, because the filename the
build reads is not the name the site produced.

Two remain absent and are reported on every build: `wang_omim_2023` and `sagas`.

## Layout

```
gene_spine/     the build: fetch, parsers, hgnc resolver, tiers, build
config/         sources.yaml, the only source declaration
curated/        organizations.csv, the authoritative organisation table
sources/        source snapshots, plus manifest.json
outputs/        everything the build writes and Pages serves
tests/          pytest, run before any release
schema/         gene_spine.schema.json
*.html          the site, served from the repository root
```
