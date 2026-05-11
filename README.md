# FluGenome3D

**FluGenome3D is a deployable visual lab for real derived Influenza A HA/NA research artifacts.**

It grew out of my AntigenLM/AntigenSDE thesis work as a satellite project: the Python pipeline audits sequence context, codon/CDS behavior, deterministic tokenization, learned AntigenLM geometry and public structure alignment QC; the web app turns derived artifacts into a shareable interface with aggregate views and hash-based identifiers.

## What It Is

- A reproducible audit of Influenza A HA/NA sequence-context metrics.
- A deterministic tokenization and representation explorer.
- A learned AntigenLM latent-geometry atlas connected to the parent thesis repository.
- A structure-aware alignment-QC layer for public HA/NA PDB entries.
- A grounded in-app guide that explains formulas, views and limits from safe project artifacts.
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
make phase7-9
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

The app has nine views:

- **Home / Overview**: cinematic entry point for the visual lab.
- **Project Guide**: plain-language overview of project intent, formulas, data layers and current models.
- **Ask FluGenome3D**: grounded explanatory guide using safe docs, reports and manifests; no raw sequence access and no external model call.
- **Dataset Atlas**: panel counts, subtype/protein balance, temporal distribution, CDS reliability and deduplication summaries.
- **AntigenLM Latent Atlas**: learned HA+NA embedding geometry from the parent thesis repo, shown as hash-based PCA coordinates with aggregate geometry diagnostics.
- **Representation Projector**: TensorFlow Projector-style scatter maps from real reduced-coordinate artifacts with safe hashed IDs.
- **Sequence/Token Inspector**: aggregate GC/CpG/UpA metrics, token entropy, effective vocabulary and stability summaries.
- **3D Molecular Viewer**: public RCSB structures `3LZG`, `3VUN`, `3NSS`, `6BR6` using 3Dmol.js, with alignment QC and aggregate residue-signal summaries.
- **Bridge View**: integrated group-level view connecting sequence context, representation maps and structure catalog entries.

The visual language is a dark minimal research lab style inspired by The Velveteen Project and TensorFlow Projector, without copying assets or layouts.

## Ask FluGenome3D

The `Ask FluGenome3D` view is a local, grounded guide for visitors. It retrieves short explanation chunks from `app/data/lab_guide.safe.json`, cites the source docs/reports, and keeps answers inside the descriptive scope of the project.

It does not call an external LLM, does not access raw sequence files, and cannot answer from restricted local panels.

## Safe Export Files

`data_export/export_vercel_safe_bundle.py` writes:

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

These files are derived from real Phase 0-9 outputs and public project documentation. They are not simulated substitutes.

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
- AntigenLM latent coordinates are shown as a descriptive learned-representation layer.
- CDS-dependent views are restricted to refined CDS subsets.
- Public PDB structures are shown with alignment QC and aggregate residue-signal context.

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
