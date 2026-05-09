# FluGenome3D Phase 1 dataset report

Generated from local restricted parent-repository data. This report contains only aggregate counts and no raw nucleotide sequences.

## Detected datasets

- FASTA files detected: 8
- Paired HA/NA JSON datasets detected: 2
- Rich metadata CSV detected: True
- Dedup metadata CSV detected: True

## Columns used

- Paired JSON: `epi_isl`, `subtype`, `year`, `month`, `day`, `ha_sequence`, `na_sequence`.
- Rich metadata: `Isolate_Id`, `Subtype`, `Host`, `Location`, `Collection_Date`, `Clade`.
- Dedup metadata: `epi_isl`, `subtype`, `year`, `month`, `matched`, `clade`, `major_clade`.

## Filtering and labels

- Subtypes retained: H1N1 and H3N2.
- Segment labels: paired HA and paired NA.
- Date availability: year, month, and day required for `paired_pass_qc`.
- HA expected length: 1650-1800 nt.
- NA expected length: 1200-1600 nt.
- Maximum ambiguous fraction: 0.01.
- Alphabet filter: A/C/G/T/N only.
- Duplicate definition: exact SHA256 identity of the HA+NA sequence-hash pair.

## Pair counts by stage

| stage | n_pairs | n_ha_sequences | n_na_sequences | n_sequence_records | n_unique_pair_hashes | n_duplicate_pair_records | n_h1n1 | n_h3n2 | n_rich_metadata | n_human_host | year_min | year_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_valid_paired | 111756 | 111756 | 111756 | 223512 | 82306 | 29450 | 46125 | 65631 | 81943 | 81943 | 2000 | 2022 |
| exact_pair_deduplicated | 82306 | 82306 | 82306 | 164612 | 82306 | 0 | 36753 | 45553 | 81943 | 81943 | 2000 | 2022 |
| metadata_rich_subset | 81943 | 81943 | 81943 | 163886 | 81943 | 0 | 36723 | 45220 | 81943 | 81943 | 2000 | 2022 |
| human_metadata_subset | 81943 | 81943 | 81943 | 163886 | 81943 | 0 | 36723 | 45220 | 81943 | 81943 | 2000 | 2022 |
| smoke_panel | 400 | 400 | 400 | 800 | 400 | 0 | 200 | 200 | 400 | 400 | 2000 | 2022 |
| mvp_panel | 10000 | 10000 | 10000 | 20000 | 10000 | 0 | 5000 | 5000 | 10000 | 10000 | 2000 | 2022 |
| full_panel | 82306 | 82306 | 82306 | 164612 | 82306 | 0 | 36753 | 45553 | 81943 | 81943 | 2000 | 2022 |

## Panel definitions

- `smoke_panel`: small deduplicated local panel for fast tests.
- `mvp_panel`: 10,000-pair deduplicated panel, balanced by subtype and spread across years where possible, preferring rich human metadata.
- `full_panel`: all 82,306 exact HA+NA deduplicated pairs passing Phase 1 QC.

| panel | subtype | n_pairs | n_sequence_records | n_unique_pair_hashes | n_duplicate_pair_records | n_rich_metadata | n_human_host | year_min | year_max |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| smoke_panel | H1N1 | 200 | 400 | 200 | 0 | 200 | 200 | 2000 | 2022 |
| smoke_panel | H3N2 | 200 | 400 | 200 | 0 | 200 | 200 | 2000 | 2022 |
| mvp_panel | H1N1 | 5000 | 10000 | 5000 | 0 | 5000 | 5000 | 2000 | 2022 |
| mvp_panel | H3N2 | 5000 | 10000 | 5000 | 0 | 5000 | 5000 | 2000 | 2022 |
| full_panel | H1N1 | 36753 | 73506 | 36753 | 0 | 36723 | 36723 | 2000 | 2022 |
| full_panel | H3N2 | 45553 | 91106 | 45553 | 0 | 45220 | 45220 | 2000 | 2022 |

## Duplicate accounting

| panel | n_pairs | n_unique_pair_hashes | n_duplicate_pair_records |
| --- | --- | --- | --- |
| all_valid | 111756 | 82306 | 29450 |
| all_valid:H1N1 | 46125 | 36753 | 9372 |
| all_valid:H3N2 | 65631 | 45553 | 20078 |
| smoke_panel | 400 | 400 | 0 |
| smoke_panel:H1N1 | 200 | 200 | 0 |
| smoke_panel:H3N2 | 200 | 200 | 0 |
| mvp_panel | 10000 | 10000 | 0 |
| mvp_panel:H1N1 | 5000 | 5000 | 0 |
| mvp_panel:H3N2 | 5000 | 5000 | 0 |
| full_panel | 82306 | 82306 | 0 |
| full_panel:H1N1 | 36753 | 36753 | 0 |
| full_panel:H3N2 | 45553 | 45553 | 0 |

## Metadata missingness

| field | n_rows | missing | missing_fraction |
| --- | --- | --- | --- |
| subtype | 111756 | 0 | 0 |
| year | 111756 | 0 | 0 |
| month | 111756 | 0 | 0 |
| day | 111756 | 0 | 0 |
| host | 111756 | 29813 | 0.266769 |
| region | 111756 | 29813 | 0.266769 |
| country | 111756 | 29828 | 0.266903 |
| metadata_collection_date | 111756 | 29813 | 0.266769 |
| major_clade | 111756 | 29813 | 0.266769 |
| dedup_matched | 111756 | 29450 | 0.263521 |

## Year distribution

The full panel spans 2000-2022. Counts by year/subtype are stored in `results/tables/phase1_year_subtype_counts.csv`.

## Length and ambiguity QC

| panel | subtype | segment | n | min | p05 | median | mean | p95 | max | n_any_ambiguity | mean_ambiguous_fraction | n_non_acgtn | n_pass_qc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_valid | H1N1 | HA | 46125 | 1695 | 1701 | 1734 | 1729.01 | 1777 | 1799 | 363 | 1.40458e-05 | 0 | 46125 |
| all_valid | H1N1 | NA | 46125 | 1350 | 1410 | 1420 | 1423.94 | 1458 | 1493 | 362 | 1.27732e-05 | 0 | 46125 |
| all_valid | H3N2 | HA | 65631 | 1672 | 1701 | 1735 | 1723.94 | 1762 | 1800 | 940 | 2.33979e-05 | 0 | 65631 |
| all_valid | H3N2 | NA | 65631 | 1251 | 1410 | 1436 | 1429.42 | 1466 | 1590 | 1664 | 4.28717e-05 | 0 | 65631 |
| full_panel | H1N1 | HA | 36753 | 1695 | 1701 | 1734 | 1729.18 | 1777 | 1799 | 350 | 1.69746e-05 | 0 | 36753 |
| full_panel | H1N1 | NA | 36753 | 1350 | 1410 | 1420 | 1424.08 | 1458 | 1493 | 333 | 1.51374e-05 | 0 | 36753 |
| full_panel | H3N2 | HA | 45553 | 1672 | 1701 | 1732 | 1723.55 | 1762 | 1800 | 897 | 3.27243e-05 | 0 | 45553 |
| full_panel | H3N2 | NA | 45553 | 1251 | 1410 | 1432 | 1428.76 | 1466 | 1590 | 1255 | 5.26804e-05 | 0 | 45553 |
| mvp_panel | H1N1 | HA | 5000 | 1695 | 1699.9 | 1724 | 1726.44 | 1777 | 1797 | 39 | 9.80244e-06 | 0 | 5000 |
| mvp_panel | H1N1 | NA | 5000 | 1350 | 1410 | 1418 | 1423.27 | 1458 | 1478 | 35 | 9.23204e-06 | 0 | 5000 |
| mvp_panel | H3N2 | HA | 5000 | 1701 | 1701 | 1721 | 1721.38 | 1762 | 1790 | 45 | 1.76337e-05 | 0 | 5000 |
| mvp_panel | H3N2 | NA | 5000 | 1398 | 1410 | 1431 | 1428.86 | 1466 | 1510 | 79 | 3.23446e-05 | 0 | 5000 |
| smoke_panel | H1N1 | HA | 200 | 1698 | 1698 | 1732.5 | 1728.61 | 1777 | 1778 | 1 | 8.66051e-06 | 0 | 200 |
| smoke_panel | H1N1 | NA | 200 | 1410 | 1410 | 1424 | 1424.72 | 1458 | 1463 | 4 | 3.14737e-05 | 0 | 200 |
| smoke_panel | H3N2 | HA | 200 | 1701 | 1701 | 1718 | 1718.86 | 1751 | 1774 | 5 | 5.74496e-05 | 0 | 200 |
| smoke_panel | H3N2 | NA | 200 | 1407 | 1410 | 1431 | 1428.38 | 1466 | 1467 | 7 | 6.91425e-05 | 0 | 200 |

## Local-only outputs

- `data/processed/panels/smoke_panel.parquet`
- `data/processed/panels/mvp_panel.parquet`
- `data/processed/panels/full_panel.parquet`

These Parquet files may contain raw sequences and accession-level fields. They must remain gitignored.

## GitHub-safe aggregate outputs

- `results/tables/phase1_dataset_summary.csv`
- `results/tables/phase1_panel_summary.csv`
- `results/tables/phase1_missingness_summary.csv`
- `results/tables/phase1_year_subtype_counts.csv`
- `results/tables/phase1_length_qc_summary.csv`
- `results/tables/phase1_duplicate_summary.csv`
- `results/figures/fig1_dataset_overview.png`
- `results/figures/fig2_year_subtype_distribution.png`
- `results/figures/fig3_length_qc_distribution.png`

## Limitations

- This phase constructs descriptive local panels only.
- Host is confirmed from the rich metadata subset; not every valid paired record has rich metadata coverage.
- Exact deduplication removes identical HA+NA pairs but does not collapse near-duplicates.
- No GROVER tokenizer, 3D structure mapping, predictive modeling, antigenicity, vaccine, escape, fitness, or optimization analysis is implemented in Phase 1.
