#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO

from flugenome3d.plotting import save_tokenization_summary
from flugenome3d.tokenization import audit_records, train_basic_bpe_from_records


def main() -> None:
    parser = argparse.ArgumentParser(description="Run lightweight tokenization audit inspired by DNA language models.")
    parser.add_argument("--fasta", default="data/processed/ha_na_sequences.fasta")
    parser.add_argument("--out", default="results/tables/tokenization_metrics.csv")
    parser.add_argument("--summary", default="results/tables/tokenization_summary.csv")
    parser.add_argument("--train-bpe", action="store_true")
    parser.add_argument("--bpe-out", default="results/tables/bpe_tokenizer.json")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    metrics = audit_records(records)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.out, index=False)
    save_tokenization_summary(metrics, args.summary)
    if args.train_bpe:
        train_basic_bpe_from_records(records, out_path=args.bpe_out)
        print(f"Wrote exploratory BPE tokenizer: {args.bpe_out}")
    print(f"Wrote tokenization metrics for {len(records)} records x methods.")


if __name__ == "__main__":
    main()
