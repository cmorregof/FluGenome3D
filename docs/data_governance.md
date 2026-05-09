# Data governance

## Local full mode

Local full mode may use non-redistributable GISAID-backed data already available on the developer's machine. Raw FASTA, metadata, BLAST databases, and accession-level restricted tables must not be committed.

Recommended policy:

- Commit scripts, configs, documentation and tests.
- Commit only aggregate derived tables if allowed.
- Do not commit raw sequences or raw metadata from restricted sources.
- Store dataset DOI/acknowledgement information in `data_manifest/dataset_doi.txt`.
- Keep `data/raw`, `data/interim`, and `data/processed` gitignored by default.

## Public demo mode

Public demo mode should use redistributable data or scripts that instruct users how to fetch public data themselves, for example through NCBI Datasets/NCBI Virus.
