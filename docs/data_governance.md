# Data governance

## Purpose

FluGenome3D is designed for a GISAID-sensitive workflow: local thesis data can be analyzed privately, while the deployed app exposes only safe derived artifacts.

## Local full mode

Raw GISAID FASTA, metadata, accessions, isolate names, BLAST databases, sequence-containing tables and restricted per-record Parquet files remain local.

Local-only paths include:

- `data/raw/`
- `data/interim/`
- `data/processed/`
- `data/processed/tokenization/`
- `data/processed/tokenization_stability/`
- `config/local_paths.yml`
- `app/data-local/`
- `app/**/*.local.json`
- `app/**/*.restricted.json`

These paths must not be committed or deployed.

## Public safe mode

The deployed app uses only safe derived artifacts in:

```text
app/data/*.safe.json
```

The deployed app may include:

- aggregate summaries;
- binned temporal/geographic summaries;
- reduced PCA/t-SNE coordinates;
- short token summaries;
- hash-based visual identifiers;
- public PDB identifiers;
- alignment-QC summaries;
- captions, guide chunks and limitations.

The deployed app excludes:

- raw HA/NA sequences;
- FASTA files;
- accessions;
- source record identifiers;
- isolate names;
- source sequence hashes;
- unrestricted sample-level metadata;
- long tokens or sequence-containing fields;
- any table that could reasonably reconstruct restricted records.

## Safe artifact rules

Safe artifacts should be derived, aggregate, binned, reduced or hash-identified. When in doubt, exclude a column from the public export.

Coordinate views are allowed only as reduced PCA/t-SNE artifacts with hash-based IDs and minimal metadata. They should not contain source identifiers, sequence hashes, accessions or exact sample-level locations.

## Pre-publication checklist

Run:

```bash
grep -R -E "[ACGTN]{80,}" app/ reports/*.md results/tables/*.csv
find app -type f -size +5M
find app -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.ffn" \)
```

Expected result:

- no long nucleotide-like strings in public app/docs/tables;
- no FASTA-like files in `app/`;
- no large files in `app/` without review.
