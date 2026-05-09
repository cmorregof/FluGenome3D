# Data manifest

This repository is designed to support two modes:

1. **Local full mode**: uses locally available non-redistributable Influenza A HA/NA data, for example GISAID-backed thesis datasets. Raw FASTA and metadata must not be committed.
2. **Vercel safe mode**: uses only real derived safe exports under `app/data/*.safe.json`.

No simulated dataset is used as the primary app data source. Only scripts, derived non-sensitive tables, aggregate figures, safe JSON exports, documentation, and reproducibility metadata should be committed.
