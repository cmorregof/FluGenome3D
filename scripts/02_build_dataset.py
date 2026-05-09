#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from flugenome3d.filters import filter_records, summarize_records
from flugenome3d.io import write_fasta
from flugenome3d.utils import load_yaml


def find_fasta_files(roots: list[str]) -> list[Path]:
    exts = {".fa", ".fasta", ".fna", ".ffn"}
    files = []
    for root in roots:
        p = Path(root).expanduser()
        if p.is_file() and p.suffix.lower() in exts:
            files.append(p)
        elif p.exists():
            files.extend([x for x in p.rglob("*") if x.is_file() and x.suffix.lower() in exts])
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description="Build local non-redistributable HA/NA working dataset.")
    parser.add_argument("--config", default="config/local_paths.yml")
    parser.add_argument("--filters", default="config/filters.yml")
    parser.add_argument("--out-fasta", default="data/processed/ha_na_sequences.fasta")
    parser.add_argument("--out-meta", default="data/processed/ha_na_metadata.parquet")
    parser.add_argument("--summary", default="results/tables/dataset_summary.csv")
    args = parser.parse_args()

    cfg = load_yaml(args.config)
    filters = load_yaml(args.filters)
    roots = cfg.get("local_data_roots", [])
    fasta_files = find_fasta_files(roots)
    if not fasta_files:
        raise SystemExit("No FASTA files discovered in local_data_roots. Run scripts/01_inventory_data.py first.")

    records = []
    for fasta in fasta_files:
        records.extend(list(SeqIO.parse(str(fasta), "fasta")))
    kept = filter_records(records, filters)

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)
    write_fasta(kept, args.out_fasta)
    meta = summarize_records(kept)
    Path(args.out_meta).parent.mkdir(parents=True, exist_ok=True)
    meta.to_parquet(args.out_meta, index=False)

    summary = meta.groupby(["segment_inferred", "subtype_inferred"], dropna=False).size().reset_index(name="n_sequences")
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.summary, index=False)

    print(f"Input FASTA files: {len(fasta_files)}")
    print(f"Input records: {len(records)}")
    print(f"Kept records: {len(kept)}")
    print(f"Wrote {args.out_fasta}, {args.out_meta}, {args.summary}")


if __name__ == "__main__":
    main()
