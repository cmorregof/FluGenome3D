# FluGenome3D Fase 0: data inventory report

Audit date: 2026-05-09
Scope: local read-only audit of the parent repository data used by the AntigenLM replication. No raw sequences, accessions, isolate names, or restricted metadata are reproduced here.

## Data roots inspected

The FluGenome3D child repository currently has empty, gitignored data directories:

- `data/raw/.gitkeep`
- `data/interim/.gitkeep`
- `data/processed/.gitkeep`

For this Fase 0, the parent repository is the local data source:

- `../data/gisaid/`: GISAID FASTA files.
- `../data/gisaid_metadata_private/`: GISAID metadata, accession/EPI_SET files, and dedup caches.
- `../data/processed_gisaid/`: AntigenLM replication JSON datasets.

No BLAST index files (`*.nhr`, `*.nin`, `*.nsq`, `*.phr`, `*.pin`, `*.psq`, etc.) were found under the current workspace during this audit. No Parquet files were found.

## Relevant files found

### FASTA

| File | Size | Records | Subtype | HA records | NA records |
|---|---:|---:|---|---:|---:|
| `../data/gisaid/GISAID_H1N1_2000_2015.fasta` | 66,036,501 B | 39,523 | H1N1 | 19,758 | 19,765 |
| `../data/gisaid/GISAID_H1N1_2016_2019.fasta` | 65,889,960 B | 39,215 | H1N1 | 19,608 | 19,607 |
| `../data/gisaid/GISAID_H1N1_2019_2022.fasta` | 57,781,722 B | 34,334 | H1N1 | 17,166 | 17,168 |
| `../data/gisaid/GISAID_H3N2_2000_2014.fasta` | 60,435,974 B | 35,976 | H3N2 | 19,287 | 16,689 |
| `../data/gisaid/GISAID_H3N2_2015_2017.fasta` | 64,563,005 B | 38,523 | H3N2 | 19,263 | 19,260 |
| `../data/gisaid/GISAID_H3N2_2017_2019.fasta` | 64,770,376 B | 38,495 | H3N2 | 19,249 | 19,246 |
| `../data/gisaid/GISAID_H3N2_2019_2021.fasta` | 66,842,346 B | 39,715 | H3N2 | 19,855 | 19,860 |
| `../data/gisaid/GISAID_H3N2_2021_2022.fasta` | 15,329,640 B | 9,090 | H3N2 | 4,542 | 4,548 |

FASTA header structure observed from local parser conventions:

`EPI_ISL|isolate_name|subtype|collection_date|segment`

Segment appears as text (`HA`, `NA`) or segment numbers (`4`, `6`). The audit normalized `4` to HA and `6` to NA.

### Metadata and accession tables

| File | Rows | Columns | Notes |
|---|---:|---:|---|
| `../data/gisaid_metadata_private/gisaid_epiflu_isolates_2000_2022_epi_set_260506bu_combined.csv` | 81,943 | 55 | Rich isolate metadata with host, location, date, subtype, clade, segment IDs. |
| `../data/gisaid_metadata_private/gisaid_epiflu_isolates_2000_2022_epi_set_260506bu_joined_dedup_cache.csv` | 82,306 | 12 | Dedup-level join cache with subtype, year/month, clade, lineage, genotype fields. |
| `../data/gisaid_metadata_private/gisaid_episet_accessions_full_all_111756.csv` | 111,756 | 1 | EPI_ISL only; restricted accession-level table. |
| `../data/gisaid_metadata_private/gisaid_episet_accessions_dedup_ha_na_all.csv` | 82,306 | 1 | EPI_ISL only; restricted accession-level table. |
| `../data/gisaid_metadata_private/gisaid_epiflu_isolates*.xls` | not parsed | not parsed | Present locally, but `xlrd` is not installed; the combined CSV appears to be the usable export. |

### Processed JSON

| File | Paired strains | Top-level keys | Record keys |
|---|---:|---|---|
| `../data/processed_gisaid/dataset_H1N1.json` | 46,125 | `subtype`, `paired_strains`, `monthly_groups`, `windows`, `stats` | `epi_isl`, `strain_name`, `subtype`, `subtype_token`, `year`, `month`, `day`, `ha_sequence`, `na_sequence` |
| `../data/processed_gisaid/dataset_H3N2.json` | 65,631 | `subtype`, `paired_strains`, `monthly_groups`, `windows`, `stats` | `epi_isl`, `strain_name`, `subtype`, `subtype_token`, `year`, `month`, `day`, `ha_sequence`, `na_sequence` |

These JSON files contain raw nucleotide sequences and must remain local/gitignored.

## Metadata columns detected

Rich combined metadata (`gisaid_epiflu_isolates_2000_2022_epi_set_260506bu_combined.csv`) includes:

- subtype: `Subtype`
- segment/protein: `HA Segment_Id`, `NA Segment_Id`, plus other segment ID columns
- host: `Host`
- country/region: `Location`
- collection date/year: `Collection_Date`
- accession/isolate: `Isolate_Id`, `Isolate_Name`
- clade: `Clade`
- additional fields: `Genotype`, `Lineage`, `Pathogenicity`, passage, submitter, publication, resistance fields, demographics, source fields, source file/sheet.

Dedup cache (`gisaid_epiflu_isolates_2000_2022_epi_set_260506bu_joined_dedup_cache.csv`) includes:

- accession/isolate: `epi_isl`
- subtype: `subtype`
- collection date/year: `year`, `month`
- clade: `clade_raw`, `clade`, `major_clade`
- lineage/genotype: `lineage_raw`, `lineage`, `genotype_raw`, `genotype`
- match flag: `matched`

## Sequence counts

### Raw FASTA records

| Category | Count |
|---|---:|
| Total FASTA records | 274,871 |
| HA records | 138,728 |
| NA records | 136,143 |
| H1N1 records | 113,072 |
| H3N2 records | 161,799 |
| H1N1 HA | 56,532 |
| H1N1 NA | 56,540 |
| H3N2 HA | 82,196 |
| H3N2 NA | 79,603 |

All parsed FASTA records had valid 5-part headers and parseable collection dates.

### Valid paired processed records

The AntigenLM replication processed data contains paired HA/NA strains:

| Subtype | Valid paired strains | HA sequences | NA sequences | Total sequence records |
|---|---:|---:|---:|---:|
| H1N1 | 46,125 | 46,125 | 46,125 | 92,250 |
| H3N2 | 65,631 | 65,631 | 65,631 | 131,262 |
| Total | 111,756 | 111,756 | 111,756 | 223,512 |

The processed `stats.txt` reports:

- H1N1: 46,125 valid strains, 198 months, 195 temporal windows.
- H3N2: 65,631 valid strains, 242 months, 239 temporal windows.
- Total: 111,756 valid paired strains.

### Host metadata

Host is not encoded in the FASTA headers or processed JSON records. The rich metadata CSV joins to 81,943 processed EPI_ISL records. Within that joined rich metadata subset:

| Host | Isolate rows | Implied paired HA/NA sequence records |
|---|---:|---:|
| Human | 81,943 | 163,886 |

Interpretation: all rich metadata rows with host information are `Human`, but the full 111,756 paired-strain panel has host metadata coverage for 81,943 records in the current rich metadata export. The remaining full-panel records should not be described as independently host-verified from the files audited here.

### Complete or near-complete sequence records

Near-complete was audited using the current thesis filter windows:

- HA: length 1650-1800 nt, ambiguous fraction <= 1%, alphabet limited to A/C/G/T/N.
- NA: length 1200-1600 nt, ambiguous fraction <= 1%, alphabet limited to A/C/G/T/N.

| Segment/subtype | Raw records | Near-complete records |
|---|---:|---:|
| H1N1 HA | 56,532 | 48,631 |
| H1N1 NA | 56,540 | 52,203 |
| H3N2 HA | 82,196 | 74,671 |
| H3N2 NA | 79,603 | 70,447 |
| HA total | 138,728 | 123,302 |
| NA total | 136,143 | 122,650 |
| Total | 274,871 | 245,952 |

The processed paired data is stricter: each retained strain has both HA and NA passing the local quality filters.

## Quality audit

### Length distributions

| Group | n | Min | p05 | Median | Mean | p95 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|
| H1N1 HA raw | 56,532 | 284 | 1701 | 1726 | 1727.85 | 1777 | 1900 |
| H1N1 NA raw | 56,540 | 839 | 1410 | 1420 | 1423.45 | 1458 | 2278 |
| H3N2 HA raw | 82,196 | 586 | 1701 | 1731 | 1722.45 | 1762 | 1901 |
| H3N2 NA raw | 79,603 | 246 | 1410 | 1432 | 1428.55 | 1466 | 1666 |
| H1N1 HA processed | 46,125 | 1695 | 1701 | 1734 | 1729.01 | 1777 | 1799 |
| H1N1 NA processed | 46,125 | 1350 | 1410 | 1420 | 1423.94 | 1458 | 1493 |
| H3N2 HA processed | 65,631 | 1672 | 1701 | 1735 | 1723.94 | 1762 | 1800 |
| H3N2 NA processed | 65,631 | 1251 | 1410 | 1436 | 1429.42 | 1466 | 1590 |

Raw FASTA contains short/long outliers. The processed JSONs remove those outliers for both segments.

### Ambiguous or non-ACGT characters

| Group | Records with any non-ACGT character | Records with >1% non-ACGT | Records with non-ACGTN character |
|---|---:|---:|---:|
| H1N1 HA raw | 8,205 | 238 | 7,609 |
| H1N1 NA raw | 4,725 | 210 | 4,133 |
| H3N2 HA raw | 8,412 | 459 | 6,997 |
| H3N2 NA raw | 10,808 | 479 | 8,656 |
| HA total raw | 16,617 | 697 | 14,606 |
| NA total raw | 15,533 | 689 | 12,789 |

Processed paired JSON still contains some `N` ambiguity but no non-ACGTN characters:

| Subtype | HA ambiguous records | HA non-ACGTN | NA ambiguous records | NA non-ACGTN |
|---|---:|---:|---:|---:|
| H1N1 | 363 | 0 | 362 | 0 |
| H3N2 | 940 | 0 | 1,664 | 0 |

### Exact duplicates

Raw FASTA exact duplicate counts were computed from sequence hashes only; no sequence strings are exported.

| Group | Records | Unique exact sequences | Duplicate records |
|---|---:|---:|---:|
| H1N1 HA raw | 56,532 | 39,479 | 17,053 |
| H1N1 NA raw | 56,540 | 33,453 | 23,087 |
| H3N2 HA raw | 82,196 | 48,159 | 34,037 |
| H3N2 NA raw | 79,603 | 42,994 | 36,609 |
| HA total raw | 138,728 | 87,636 | 51,092 |
| NA total raw | 136,143 | 76,445 | 59,698 |

Processed paired exact duplicates:

| Dataset | Records | Unique exact units | Duplicate records |
|---|---:|---:|---:|
| H1N1 HA processed | 46,125 | 30,626 | 15,499 |
| H1N1 NA processed | 46,125 | 26,643 | 19,482 |
| H3N2 HA processed | 65,631 | 36,287 | 29,344 |
| H3N2 NA processed | 65,631 | 33,176 | 32,455 |
| HA+NA exact paired processed | 111,756 | 82,306 | 29,450 |

The deduplicated HA+NA pair universe of 82,306 records is therefore the strongest candidate for representation and distance-style descriptive analyses.

### Missing metadata

Rich metadata CSV:

| Field | Missing rows |
|---|---:|
| `Isolate_Id` | 0 |
| `Isolate_Name` | 0 |
| `Subtype` | 0 |
| `HA Segment_Id` | 0 |
| `NA Segment_Id` | 0 |
| `Host` | 0 |
| `Location` | 0 |
| `Collection_Date` | 0 |
| `Clade` | 0 |
| `Lineage` | 50,638 |
| `Genotype` | 0 |

Dedup cache:

| Field | Missing rows |
|---|---:|
| `epi_isl` | 0 |
| `subtype` | 0 |
| `year` | 0 |
| `month` | 0 |
| `clade_raw` | 363 |
| `clade` | 6,068 |
| `major_clade` | 6,068 |
| `lineage_raw` | 51,001 |
| `lineage` | 51,001 |
| `genotype_raw` | 363 |
| `genotype` | 82,306 |

Rich metadata geography coverage is broad but uneven. Top countries in the rich metadata are United States (24,585), Australia (5,590), China (3,985), United Kingdom (3,397), Russian Federation (2,720), New Zealand (2,348), Brazil (2,074), Singapore (2,070), France (1,902), and Netherlands (1,618).

## Audit conclusions

1. The local dataset is real and large enough for the proposed FluGenome3D MVP.
2. The core usable unit is a paired HA/NA nucleotide record, not a protein-only record.
3. The strongest full local panel is the 111,756 valid paired-strain JSON set.
4. The strongest duplicate-controlled panel is the 82,306 exact HA+NA deduplicated set.
5. Host/country/clade metadata should be treated as rich but not full-panel universal: host/country coverage is confirmed for 81,943 rich metadata rows; major clade coverage is available for 76,238 of 82,306 dedup rows.
6. Raw FASTA and processed JSON files contain restricted sequence data and must not be copied into FluGenome3D or committed.
