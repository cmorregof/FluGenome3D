# Reproducibility

## Environment

Use:

```bash
conda env create -f environment.yml
conda activate flugenome3d
pip install -e .
```

## Minimal pipeline

```bash
make inventory CONFIG=config/local_paths.yml
make build CONFIG=config/local_paths.yml
make metrics
make tokenization
make structure
```

## Reproducibility philosophy

The repository should reproduce code, filters, figures and aggregate descriptive outputs. It should not redistribute raw restricted biological datasets.
