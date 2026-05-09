#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from flugenome3d.codon_usage import codon_usage_table, rscu_table
from flugenome3d.plotting import save_metric_histograms
from flugenome3d.sequence_metrics import sequence_record_metrics, translate_qc


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute descriptive sequence-context metrics.")
    parser.add_argument("--fasta", default="data/processed/ha_na_sequences.fasta")
    parser.add_argument("--metadata", default="data/processed/ha_na_metadata.parquet")
    parser.add_argument("--out", default="results/tables/sequence_metrics.csv")
    parser.add_argument("--codon-out", default="results/tables/codon_usage.csv")
    parser.add_argument("--rscu-out", default="results/tables/rscu.csv")
    parser.add_argument("--figures", default="results/figures")
    args = parser.parse_args()

    records = list(SeqIO.parse(args.fasta, "fasta"))
    metrics = sequence_record_metrics(records)
    qc = pd.DataFrame([{"seq_id": r.id, **translate_qc(str(r.seq))} for r in records])
    metrics = metrics.merge(qc, on="seq_id", how="left")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(args.out, index=False)
    codon_usage_table(records).to_csv(args.codon_out, index=False)
    rscu_table(records).to_csv(args.rscu_out, index=False)
    save_metric_histograms(metrics, args.figures)
    print(f"Wrote sequence metrics for {len(records)} records.")


if __name__ == "__main__":
    main()
