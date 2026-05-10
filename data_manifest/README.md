# Data manifest

This repository is designed to support two modes:

1. **Local full mode**: uses locally available non-redistributable Influenza A HA/NA data, for example GISAID-backed thesis datasets. Raw FASTA and metadata must not be committed.
2. **Cryptographic derived-data layer**: uses only real derived exports under `app/data/*.safe.json`, with aggregate views and hash-based IDs.

No simulated dataset is used as the primary app data source. Only scripts, derived non-sensitive tables, aggregate figures, safe JSON exports, documentation, and reproducibility metadata should be committed.
