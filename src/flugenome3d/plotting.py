from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def save_metric_histograms(metrics: pd.DataFrame, outdir: str | Path) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    for col in ["length", "gc_fraction", "cpg_oe", "upa_oe", "ambiguous_fraction"]:
        if col not in metrics.columns:
            continue
        plt.figure(figsize=(7, 4))
        metrics[col].dropna().hist(bins=30)
        plt.title(col)
        plt.xlabel(col)
        plt.ylabel("count")
        plt.tight_layout()
        plt.savefig(out / f"{col}_histogram.png", dpi=200)
        plt.close()


def save_tokenization_summary(metrics: pd.DataFrame, outpath: str | Path) -> None:
    summary = metrics.groupby("method", as_index=False).agg(
        mean_tokens_per_kb=("tokens_per_kb", "mean"),
        mean_token_length=("mean_token_length", "mean"),
        mean_boundary_crossing=("cross_codon_boundary_fraction", "mean"),
        mean_cpg_fraction=("cpg_token_fraction", "mean"),
        mean_upa_fraction=("upa_token_fraction", "mean"),
    )
    p = Path(outpath)
    p.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(p, index=False)
