# FluGenome3D Fase 0: dataset schema

This schema separates local restricted data from GitHub-safe derived outputs. Fields marked `local only` must not be committed when derived from GISAID data.

## Source file inventory

Suggested table: `results/tables/data_inventory_aggregate.csv`

GitHub status: allowed if aggregate only.

| Field | Type | Status | Description |
|---|---|---|---|
| `source_root` | string | aggregate | Local root label, e.g. `parent_gisaid`. |
| `relative_path` | string | aggregate | Relative path without copying file contents. |
| `file_type` | enum | aggregate | `fasta`, `csv`, `xls`, `json`, `txt`, `other`. |
| `size_bytes` | integer | aggregate | File size. |
| `n_records` | integer | aggregate | Record count if parseable. |
| `n_columns` | integer | aggregate | Column count for tabular files. |
| `contains_sequences` | boolean | aggregate | True for FASTA/processed JSON with sequence strings. |
| `contains_accessions` | boolean | aggregate | True if row-level accessions/EPI_ISL are present. |
| `redistribution_class` | enum | aggregate | `restricted_local`, `public`, `unknown`. |

## Raw FASTA record schema

Suggested local table: `data/interim/raw_fasta_records.parquet`

GitHub status: local only.

| Field | Type | Status | Description |
|---|---|---|---|
| `record_id_local` | string | local only | Internal non-public ID. |
| `epi_isl` | string | local only | GISAID isolate identifier. |
| `isolate_name` | string | local only | Isolate name from FASTA header. |
| `subtype` | enum | local only | Normalized `H1N1` or `H3N2`. |
| `segment` | enum | local only | Normalized `HA` or `NA`. |
| `collection_date` | date/string | local only | Raw parsed date. |
| `year` | integer | aggregate allowed if binned | Parsed collection year. |
| `month` | integer | aggregate allowed if binned | Parsed collection month. |
| `sequence` | string | local only | Nucleotide sequence. |
| `sequence_sha256` | string | local only | Hash for duplicate detection. |
| `source_file` | string | local only | FASTA source file. |

## Rich metadata isolate schema

Suggested local table: `data/interim/gisaid_metadata_rich.parquet`

GitHub status: local only at row level; aggregate summaries only.

| Field | Type | Source column | Status |
|---|---|---|---|
| `epi_isl` | string | `Isolate_Id` | local only |
| `isolate_name` | string | `Isolate_Name` | local only |
| `subtype` | enum | `Subtype` | aggregate allowed |
| `ha_segment_id` | string | `HA Segment_Id` | local only |
| `na_segment_id` | string | `NA Segment_Id` | local only |
| `host` | string | `Host` | aggregate allowed |
| `location_raw` | string | `Location` | local only |
| `region` | string | parsed from `Location` | aggregate allowed |
| `country` | string | parsed from `Location` | aggregate allowed if coarse enough |
| `collection_date` | date/string | `Collection_Date` | local only |
| `year` | integer | parsed from `Collection_Date` | aggregate allowed |
| `clade_raw` | string | `Clade` | aggregate allowed |
| `major_clade` | string | normalized from `Clade` | aggregate allowed |
| `lineage_raw` | string | `Lineage` | aggregate allowed |
| `genotype_raw` | string | `Genotype` | aggregate allowed |

## Dedup metadata schema

Suggested local table: `data/interim/gisaid_dedup_metadata.parquet`

GitHub status: local only at row level; aggregate summaries only.

| Field | Type | Source column | Status |
|---|---|---|---|
| `epi_isl` | string | `epi_isl` | local only |
| `subtype` | enum | `subtype` | aggregate allowed |
| `year` | integer | `year` | aggregate allowed |
| `month` | integer | `month` | aggregate allowed if binned |
| `matched` | boolean | `matched` | aggregate allowed |
| `clade_raw` | string | `clade_raw` | aggregate allowed |
| `clade` | string | `clade` | aggregate allowed |
| `major_clade` | string | `major_clade` | aggregate allowed |
| `lineage` | string | `lineage` | aggregate allowed |
| `genotype` | string | `genotype` | aggregate allowed |

## Paired HA/NA record schema

Suggested local table: `data/processed/paired_ha_na_records.parquet`

GitHub status: local only.

| Field | Type | Status | Description |
|---|---|---|---|
| `pair_id_local` | string | local only | Stable local ID. |
| `epi_isl` | string | local only | GISAID isolate identifier. |
| `subtype` | enum | aggregate allowed | `H1N1` or `H3N2`. |
| `year` | integer | aggregate allowed | Collection year. |
| `month` | integer | aggregate allowed if binned | Collection month. |
| `day` | integer | local only | Collection day. |
| `ha_sequence` | string | local only | HA nucleotide sequence. |
| `na_sequence` | string | local only | NA nucleotide sequence. |
| `ha_sha256` | string | local only | HA sequence hash. |
| `na_sha256` | string | local only | NA sequence hash. |
| `pair_sha256` | string | local only | HA+NA exact-pair hash. |
| `is_exact_pair_duplicate` | boolean | local only | Duplicate flag. |
| `metadata_join_status` | enum | aggregate allowed | `rich_metadata`, `dedup_only`, `missing`. |

## Sequence QC metrics schema

Suggested local row-level table: `data/processed/sequence_qc_metrics.parquet`
Suggested GitHub-safe aggregate table: `results/tables/sequence_qc_summary.csv`

| Field | Type | Row-level status | Aggregate status |
|---|---|---|---|
| `seq_id_local` | string | local only | omit |
| `subtype` | enum | local only | allowed |
| `segment` | enum | local only | allowed |
| `year_bin` | string | local only | allowed |
| `length_nt` | integer | local only | summarize only |
| `gc_fraction` | float | local only | summarize only |
| `n_fraction` | float | local only | summarize only |
| `non_acgtn_count` | integer | local only | summarize only |
| `passes_length_filter` | boolean | local only | count only |
| `passes_ambiguity_filter` | boolean | local only | count only |
| `passes_qc` | boolean | local only | count only |

## Codon usage schema

Suggested local row-level table: `data/processed/codon_usage_by_sequence.parquet`
Suggested aggregate table: `results/tables/codon_usage_summary_by_group.csv`

| Field | Type | Row-level status | Aggregate status |
|---|---|---|---|
| `seq_id_local` | string | local only | omit |
| `subtype` | enum | local only | allowed |
| `segment` | enum | local only | allowed |
| `year_bin` | string | local only | allowed |
| `codon` | string | local only | allowed |
| `amino_acid` | string | local only | allowed |
| `codon_count` | integer | local only | summarize only |
| `codon_frequency` | float | local only | summarize only |
| `rscu` | float | local only | summarize only |
| `codon_position_frame_status` | enum | local only | count only |

## Dinucleotide metrics schema

Suggested local row-level table: `data/processed/dinucleotide_metrics_by_sequence.parquet`
Suggested aggregate table: `results/tables/dinucleotide_bias_summary_by_group.csv`

| Field | Type | Row-level status | Aggregate status |
|---|---|---|---|
| `seq_id_local` | string | local only | omit |
| `subtype` | enum | local only | allowed |
| `segment` | enum | local only | allowed |
| `year_bin` | string | local only | allowed |
| `dinucleotide` | enum | local only | allowed |
| `observed_count` | integer | local only | summarize only |
| `expected_count` | float | local only | summarize only |
| `odds_ratio` | float | local only | summarize only |
| `log_odds_ratio` | float | local only | summarize only |
| `cpg_odds_ratio` | float | local only | summarize only |
| `upa_odds_ratio` | float | local only | summarize only |

## Tokenization metrics schema

Suggested local row-level table: `data/processed/tokenization_metrics_by_sequence.parquet`
Suggested aggregate table: `results/tables/tokenization_metrics_summary.csv`

| Field | Type | Row-level status | Aggregate status |
|---|---|---|---|
| `seq_id_local` | string | local only | omit |
| `subtype` | enum | local only | allowed |
| `segment` | enum | local only | allowed |
| `tokenizer_name` | string | local only | allowed |
| `tokenizer_version` | string | local only | allowed |
| `sequence_length` | integer | local only | summarize only |
| `n_tokens` | integer | local only | summarize only |
| `tokens_per_kb` | float | local only | summarize only |
| `mean_token_length` | float | local only | summarize only |
| `vocab_size_observed` | integer | local only | summarize only |
| `vocab_entropy` | float | local only | summarize only |
| `cross_codon_boundary_fraction` | float | local only | summarize only |
| `cpg_token_fraction` | float | local only | summarize only |
| `upa_token_fraction` | float | local only | summarize only |

## Structure mapping schema

Suggested local table: `data/processed/structure_mapping_metrics.parquet`
Suggested GitHub-safe output: aggregate static figure/HTML without accession-level data.

| Field | Type | Status | Description |
|---|---|---|---|
| `pdb_id` | string | public | PDB identifier. |
| `protein` | enum | public | `HA` or `NA`. |
| `subtype` | enum | public | `H1N1`/`H3N2` or `N1`/`N2` context. |
| `chain_id` | string | public | PDB chain used. |
| `pdb_residue_number` | string/integer | public | Residue identifier. |
| `alignment_status` | enum | aggregate allowed | `mapped`, `gap`, `ambiguous`, `unmapped`. |
| `codon_index` | integer | local only | Consensus/reference codon coordinate. |
| `metric_name` | string | aggregate allowed | e.g. `windowed_cpg_log_odds`. |
| `metric_value_aggregate` | float | aggregate allowed | Aggregated sequence-context signal. |
| `n_sequences_contributing` | integer | aggregate allowed | Group count only. |

## Required privacy rules

1. Never commit `sequence`, `ha_sequence`, `na_sequence`, isolate names, or accession/EPI_ISL row tables.
2. Never commit row-level GISAID-derived metadata, even without sequences.
3. Aggregate by subtype/segment/year bin/clade only when the table cannot reconstruct restricted records.
4. Keep all local GISAID-derived Parquet/CSV/JSON files under gitignored `data/interim` or `data/processed`.
5. Public demo datasets must be independently redistributable and documented separately from local GISAID mode.
