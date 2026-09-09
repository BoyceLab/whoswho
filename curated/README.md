# Curated layer

`organizations.csv` is the authoritative organization table. Rows here are never
overwritten by the build; candidate rows from scraped or membership sources are
appended with `status = candidate` and dropped automatically once a curated row
exists for the same gene and website host.

Column meanings:

- org_id: stable slug you assign (e.g. `stxbp1-foundation`)
- hgnc_id / symbol: the gene the organization serves; one row per gene for multi-gene groups
- orphacode: filled from the Orphadata bridge if blank
- org_type: foundation, family group, coalition, registry-only, professional society
- registry_platform: REDCap, Simons Searchlight, RARE-X, Ciitizen, COMBINEDBrain, Across Healthcare, other
- natural_history_study: yes / planned / no
- data_elements_documented: yes / partial / no
- irb_status: approved / pending / none / unknown
- maturity_stage: your five-stage self-assessment value
- coalitions: REN; COMBINEDBrain; AGENDA; ELC (semicolon-separated)
- status: curated | candidate | rejected. A row reads "verified" on the site only when status is
  curated AND reviewed_by, review_date, and url are all filled; otherwise it is shown as awaiting
  review. Rejected rows never appear in the public outputs.
- source: where the row came from
- reviewed_by / review_date: who confirmed it and when
