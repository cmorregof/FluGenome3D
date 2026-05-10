# FluGenome3D

**FluGenome3D is a deployable visual lab for real derived Influenza A HA/NA research artifacts.**

It grew out of my AntigenLM/AntigenSDE thesis work as a satellite project: the Python pipeline audits sequence context, codon/CDS behavior, simple representations and deterministic tokenization; the web app turns derived artifacts into a shareable interface with aggregate views and hash-based identifiers.

## What It Is

- A reproducible audit of Influenza A HA/NA sequence-context metrics.
- A deterministic tokenization and representation explorer.
- A Vercel-ready visual interface built from real derived outputs.
- A governance-first app that separates restricted local analysis from a cryptographic derived-data layer.

## What It Is Not

- Not a viral design tool.
- Not a vaccine, antigenicity, escape, fitness, pathogenicity or evolution predictor.
- Not a raw sequence redistribution site.
- Not a replacement for AntigenLM, AntigenSDE, GROVER or biological validation.

## Data Modes

### Local Full Mode

Local full mode can use restricted, gitignored artifacts under `data/processed/`, `data/processed/tokenization/`, `data/processed/tokenization_stability/` and optional `app/data-local/`.

This mode is for private analysis only. It must not be pushed to GitHub or deployed to Vercel.

```bash
cp config/local_paths.example.yml config/local_paths.yml
make phase1
make phase2-mvp
make phase3
make phase4
make phase5
make phase6
```

To let the app prefer local-only JSON files during development:

```bash
cd app
FLUGENOME3D_DATA_MODE=local npm run dev
```

### Cryptographic Data Layer

The deployable app uses only `app/data/*.safe.json`, generated from real derived artifacts. These exports use aggregate summaries, reduced coordinates, short tokens and hash-based IDs rather than raw records.

It excludes:

- raw HA/NA sequences;
- FASTA;
- restricted Parquet panels;
- accessions and isolate names;
- source sequence hashes;
- tokens longer than 6 nt;
- per-sample sensitive tables.

Generate the derived data layer:

```bash
python3 data_export/export_vercel_safe_bundle.py
```

Run the app locally:

```bash
cd app
npm install
npm run dev
```

Build for Vercel:

```bash
cd app
npm run typecheck
npm run build
```

## Web App Views

The app has seven views:

- **Home / Overview**: cinematic entry point for the visual lab.
- **Project Guide**: plain-language overview of project intent, formulas, data layers and current models.
- **Dataset Atlas**: panel counts, subtype/protein balance, temporal distribution, CDS reliability and deduplication summaries.
- **Representation Projector**: TensorFlow Projector-style scatter maps from real reduced-coordinate artifacts with safe hashed IDs.
- **Sequence/Token Inspector**: aggregate GC/CpG/UpA metrics, token entropy, effective vocabulary and stability summaries.
- **3D Molecular Viewer**: public RCSB structures `3LZG`, `3VUN`, `3NSS`, `6BR6` using 3Dmol.js. Metric-to-structure mapping is marked pending until validated.
- **Bridge View**: integrated group-level view connecting sequence context, representation maps and structure catalog entries.

The visual language is a dark minimal research lab style inspired by The Velveteen Project and TensorFlow Projector, without copying assets or layouts.

## Safe Export Files

`data_export/export_vercel_safe_bundle.py` writes:

```text
app/data/dataset_overview.safe.json
app/data/representation_maps.safe.json
app/data/metric_summaries.safe.json
app/data/tokenization_summaries.safe.json
app/data/stability_summaries.safe.json
app/data/structure_catalog.safe.json
app/data/claims_and_limits.safe.json
app/data/data_governance.safe.json
```

These files are derived from real Phase 0-6 outputs. They are not simulated substitutes.

## Data Governance

Raw GISAID-derived data, private metadata, FASTA files, restricted Parquet panels and local thesis datasets must not be committed.

Gitignored local paths include:

```text
data/raw/
data/interim/
data/processed/
config/local_paths.yml
app/data-local/
app/**/*.local.json
app/**/*.restricted.json
```

Before publishing, run:

```bash
grep -E "[ACGTN]{80,}" app/ reports/*.md results/tables/*.csv
find app -type f -size +5M
find app -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.ffn" \)
```

The first and third commands should return no matches. Any file over 5 MB in `app/` must be reviewed before deployment.

## Claims

Allowed:

- This app provides descriptive exploration of real derived FluGenome3D artifacts.
- Deterministic tokenization metrics are compared under bootstrap and temporal summaries.
- CDS-dependent views are restricted to refined CDS subsets.
- Public PDB structures are shown for structure-aware context.

Prohibited:

- Predicts antigenic drift.
- Identifies escape mutations.
- Predicts vaccine candidates.
- Explains fitness, pathogenicity, transmissibility or selection.
- Validates GROVER, BPE or any learned tokenizer.

## Developer Commands

```bash
make app-export
make app-install
make app-build
make test
```

For Vercel, set the project root to `app/`. No restricted local data is required for build.
