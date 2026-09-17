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
| `<pages-base>/outputs/epilepsy_ndd_gene_list.csv` | 877,594 | `bf7807cdfc48edd6892211f557f4cfab469db40e082dd96d3188d01227ef71b7` |
| `<pages-base>/outputs/epilepsy_orgs.csv` | 95,920 | `5711354b82623619f7b073d00c68c098f04c497fef57227106c40d9edca61823` |
| `<pages-base>/outputs/gene_spine.json` | 9,550,535 | `99af53934d055d63ef978f4b05e6c69f72f798d9d519429404c71cdf7f3fb9f6` |
| `<pages-base>/outputs/gene_spine.csv` | 2,794,563 | `c3b00968293567f2a851c14264fec492429acf5b55d54744bcedd793c81b69f6` |
| `<pages-base>/outputs/release.json` | 9,447 | `22e0f5824c556d4c842779ed586134b5fb9a893adf19279664df3816faaaa8d1` |

Those hashes are the 2026-09-16 build. **A build against fresher sources rewrites them, so read
`outputs/release.json` for the current values**; it carries a hash per source. The paths are
stable.

The hashes do not depend on where the build ran. Every output is written with LF line endings
whatever the platform, `.gitattributes` keeps the checkout consistent, and a source whose
response carries its own timing is hashed over its canonicalised records rather than its bytes,
with `sha256_scope` saying so. Two consecutive fetch-and-build cycles reproduce `release.json`
and `sources/manifest.json` byte for byte, and `tests/test_spine_outputs.py` runs the list
emitter twice and requires identical output.

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

`epilepsy_orgs.csv` is one row per organisation and gene, every row `status = candidate` until
reviewed. 770 rows over 586 organisations, drawn from every source the config marks
`role: org_candidates` that has a file, which is seven of ten. The curated layer in
`curated/organizations.csv` is never written by the build.

An organisation is identified by its normalised website host where it has one, and by its
normalised name otherwise, lowercased with punctuation stripped. Rows are then collapsed on the
name as well, because one organisation can arrive with two URLs, so `(org_name, symbol)` is
unique and the gene list's `org_count` cannot disagree with this file. Both counts are taken
from the same frame.

Both files are written by `gene_spine/epilepsy_list.py`, which runs after the spine build, so
they cannot drift away from it.

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
make build     # python -m gene_spine.build, writes the spine into outputs/
make lists     # python -m gene_spine.epilepsy_list, the two derived CSVs; run after build
make test      # python -m pytest -q tests
make release   # fetch, build, lists
```

`make lists` reads the spine the build wrote, so the order matters. The workflow runs them in
that order too.

`make fetch` prints any declared source whose file is absent, and the build continues without it.
An absent source's evidence columns are **missing from the outputs rather than False**, so a
reader cannot mistake "we did not look" for "the source said no"; `outputs/release.json` records
the absences under `missing_sources`.

## Sources

Automatic, downloaded by `make fetch`: HGNC, Genes4Epilepsy, ClinGen, PanelApp Australia
(Genetic Epilepsy, Intellectual disability), PanelApp Genomics England (Early onset or syndromic
epilepsy, Intellectual disability), SysNDD, Gene2Phenotype DD, Orphadata, and MONDO.

MONDO is a crosswalk rather than evidence: the list builder reads its `xref` lines to turn ORPHA
codes and OMIM ids into MONDO ids, and the spine build skips it. Its own API is no use here,
because OLS4's search index does not expose cross-references, so one hashed snapshot of the
ontology is both cheaper and reproducible.

Two endpoints needed care. SysNDD's browse call applies a hidden `max_category='Definitive'`
filter unless all four categories are named, which silently halves the table. Gene2Phenotype's
documented `DDG2P.csv.gz` path now answers 200 with the single-page app's HTML, so a status-code
check would take a web page for a dataset; the API download path serves the CSV.

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
