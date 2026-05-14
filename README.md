# FluGenome3D

FluGenome3D is a geometric audit and research-visualization lab for thesis-derived Influenza A HA/NA artifacts.

It converts local analyses from a master's thesis workflow on Influenza A, AntigenLM, latent representation geometry, deterministic tokenization, sequence-context metrics, and structure-aware QC into a deployed visual interface built from safe derived artifacts.

## Live app

Live app: https://flugenome3d.vercel.app/

Recommended first path:

1. Dataset Atlas
2. Sequence/Token Inspector
3. AntigenLM Latent Atlas
4. Representation Projector
5. 3D Molecular Viewer
6. Bridge View

For Vercel deployment, set the project root to `app/`.

## Scientific framing

FluGenome3D asks a conservative visualization question:

Can thesis-derived HA/NA artifacts be audited visually across sequence context, deterministic tokenization, learned representation geometry, and public structure references while preserving strict data-governance boundaries?

The project should be read as a geometric audit, latent representation analysis, reproducibility companion, and research visualization lab. It is not a forecasting system and does not make biological-causality claims.

## What it is

- A reproducible audit of Influenza A HA/NA sequence-context metrics.
- Deterministic tokenization summaries for codons, k-mers, overlapping/non-overlapping variants, and frame-aware variants.
- AntigenLM-derived latent-coordinate views from the parent thesis repository.
- PCA/t-SNE visualizations for descriptive inspection.
- Aggregate dataset, QC, temporal, geographic, subtype and segment/protein summaries.
- Public RCSB HA/NA structure viewers with conservative alignment-QC status.
- A deployed Next.js interface built from derived safe JSON files.
- A local-vs-public data boundary designed for GISAID-sensitive workflows.
- A grounded in-app guide for explaining formulas, views, governance and limitations from safe artifacts.

## What it is not

- Not antigenic drift prediction.
- Not immune escape detection.
- Not vaccine recommendation.
- Not viral fitness, pathogenicity or transmissibility prediction.
- Not viral sequence design or optimization.
- Not redistribution of raw HA/NA sequences.
- Not exposure of accessions, source record identifiers, isolate names, or restricted per-record metadata.
- Not biological validation of AntigenLM, GROVER, BPE or any tokenizer.

## Data governance

FluGenome3D separates private local analysis from the deployed safe derived-data layer.

### Local full mode

Local full mode may include restricted and gitignored artifacts such as:

- `data/raw/`
- `data/interim/`
- `data/processed/`
- `data/processed/tokenization/`
- `data/processed/tokenization_stability/`
- `app/data-local/`

These may contain sequence-bearing or sensitive derived data and must not be pushed to GitHub or deployed to Vercel.

### Public safe mode

The deployed app uses:

```text
app/data/*.safe.json
```

Safe files may include:

- aggregate summaries;
- reduced PCA/t-SNE coordinates;
- binned temporal/geographic summaries;
- short token summaries;
- hash-based point identifiers;
- public PDB identifiers;
- alignment-QC summaries;
- explanatory guide chunks.

Safe files exclude:

- raw HA/NA sequences;
- FASTA files;
- accessions, source record identifiers or accession tables;
- isolate names;
- source sequence hashes;
- unrestricted sample-level metadata;
- long tokens;
- any table that could reasonably reconstruct restricted records.

## Repository structure

```text
FluGenome3D/
  app/                  Next.js visual lab
  app/data/             Safe derived JSON exports for deployment
  config/               Example configs and local path templates
  data/                 Gitignored local data roots
  data_export/          Safe-bundle export scripts
  data_manifest/        Dataset provenance and governance notes
  docs/                 Governance, claims, reproducibility and methods docs
  reports/              Phase reports generated from local analysis
  results/              Aggregate figures and derived outputs
  scripts/              Phase-based analysis pipeline
  src/flugenome3d/      Python package
  tests/                Unit and regression tests
```

## Pipeline overview

The Makefile exposes the current analysis and app targets:

```bash
make setup
make phase1
make phase2-smoke
make phase2-mvp
make phase3
make phase4
make phase5
make phase6
make phase7-9
make app-export
make app-install
make app-dev
make app-build
make test
```

`phase7-9` runs the AntigenLM bridge, latent atlas, structure mapping, and safe export targets. The earlier phases require local restricted paths and are not runnable from the public repository alone.

## Exporting the safe data layer

The safe export script is:

```bash
python3 data_export/export_vercel_safe_bundle.py
```

It writes:

```text
app/data/dataset_overview.safe.json
app/data/representation_maps.safe.json
app/data/metric_summaries.safe.json
app/data/tokenization_summaries.safe.json
app/data/stability_summaries.safe.json
app/data/antigenlm_latent_atlas.safe.json
app/data/structure_catalog.safe.json
app/data/structure_mapping.safe.json
app/data/lab_guide.safe.json
app/data/claims_and_limits.safe.json
app/data/data_governance.safe.json
```

These files are derived from real Phase 0-9 outputs and project documentation. They are not simulated substitutes.

## Running the app locally

```bash
cd app
npm install
npm run dev
```

Build and typecheck:

```bash
cd app
npm run typecheck
npm run build
```

## Python setup

```bash
pip install -e .
pip install -e ".[dev]"
```

Run tests:

```bash
make test
```

## Safety checks before publishing

```bash
grep -R -E "[ACGTN]{80,}" app/ reports/*.md results/tables/*.csv
find app -type f -size +5M
find app -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.ffn" \)
```

The sequence and FASTA searches should return no matches. Any file over 5 MB in `app/` should be reviewed before deployment.

## Claims and limitations

Allowed claims:

- FluGenome3D provides descriptive exploration of derived HA/NA artifacts.
- It audits sequence-context summaries, deterministic tokenization behavior, and representation geometry.
- PCA and t-SNE maps are safe reduced-coordinate views for inspection.
- Structure views use public PDB entries and conservative alignment-QC status.
- CDS-dependent summaries are limited to records passing documented CDS/QC filters.
- The deployed app is designed to avoid raw sequence, accession, or isolate-name redistribution.

Disallowed claims:

- FluGenome3D predicts antigenic drift.
- FluGenome3D identifies immune escape mutations.
- FluGenome3D recommends vaccine candidates.
- FluGenome3D predicts fitness, pathogenicity, transmissibility, or selection.
- FluGenome3D validates AntigenLM, GROVER, BPE, or any tokenizer as biologically causal.
- FluGenome3D provides viral design or sequence optimization guidance.

## Suggested citation / portfolio description

FluGenome3D is a deployed research-visualization lab for thesis-derived Influenza A HA/NA artifacts. It audits sequence-context summaries, deterministic tokenization behavior, AntigenLM-derived latent geometry, and public structure-alignment QC through a GISAID-safe derived-data layer.

## Suggested figures and screenshots

1. Landing page with governance badges.
2. Dataset Atlas geographic/subtype coverage view.
3. Sequence/Token Inspector composition and tokenization diagnostics.
4. AntigenLM Latent Atlas PCA/t-SNE map.
5. Public structure viewer with alignment-QC panel.
6. Bridge View connecting sequence context, latent geometry, and structure context.
7. Data-governance boundary diagram.
8. Pipeline/reproducibility diagram.

## License

No license file is currently included in this repository. Add one before reuse or redistribution outside the intended portfolio/research context.
