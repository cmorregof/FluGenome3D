# Reproducibility

FluGenome3D is reproducible in two layers: a local full analysis layer that depends on restricted local data, and a deployable safe app layer built from derived JSON artifacts.

## Python setup

```bash
pip install -e .
pip install -e ".[dev]"
```

The package is defined in `pyproject.toml` and requires Python 3.10 or newer.

## Makefile targets

Available targets include:

```bash
make setup
make inventory
make build
make phase1
make phase2-smoke
make phase2-mvp
make phase3
make phase4
make phase5
make phase6
make phase7
make phase8
make phase9
make phase7-9
make app-export
make app-install
make app-dev
make app-build
make test
```

The targets from `phase1` through `phase9` assume local data and prior phase outputs. They should be run only in a local full environment with `config/local_paths.yml` configured.

## Safe export

The deployable derived-data layer is generated with:

```bash
python3 data_export/export_vercel_safe_bundle.py
```

or:

```bash
make app-export
```

This writes `app/data/*.safe.json` files for the Next.js app.

## Running the app

```bash
cd app
npm install
npm run dev
```

Build checks:

```bash
cd app
npm run typecheck
npm run build
```

## Tests

```bash
make test
```

## Local-data caveat

The public repository cannot reproduce restricted GISAID-derived raw inputs. It can reproduce code, safe exports, docs, app build behavior and tests. Full scientific regeneration requires the private local thesis data environment.

## Publication safety

Before pushing or deploying, run:

```bash
grep -R -E "[ACGTN]{80,}" app/ reports/*.md results/tables/*.csv
find app -type f -size +5M
find app -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" -o -name "*.ffn" \)
```
