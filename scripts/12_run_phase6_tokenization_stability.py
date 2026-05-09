#!/usr/bin/env python
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import jensenshannon

from flugenome3d.token_stability import (
    bootstrap_js_distances,
    bootstrap_token_metrics,
    bootstrap_top_tokens,
    summarize_bootstrap_ci,
    temporal_bins,
    token_metrics_by_time_window,
    top_token_jaccard_stability,
)
from flugenome3d.token_stats import token_enrichment_by_group
from flugenome3d.tokenization import (
    codon_tokenize_with_positions,
    fixed_kmer_tokenize,
    frame_aware_kmer_tokenize_with_positions,
    non_overlapping_kmer_tokenize_with_positions,
)


RAW_PANEL = Path("data/processed/panels/mvp_panel.parquet")
CDS_PANEL = Path("data/processed/panels/mvp_cds_refined_panel.parquet")
TOKEN_DIR = Path("data/processed/tokenization")
STABILITY_DIR = Path("data/processed/tokenization_stability")
TABLE_DIR = Path("results/tables")
FIGURE_DIR = Path("results/figures")
REPORT_PATH = Path("reports/phase6_tokenization_stability_report.md")
DEFAULT_CONFIG = Path("config/phase6.yml")
GROUP_ORDER = ["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"]
TARGET_COMPARISONS = {
    "HA-H1N1 vs NA-H1N1",
    "HA-H3N2 vs NA-H3N2",
    "HA-H1N1 vs HA-H3N2",
    "NA-H1N1 vs NA-H3N2",
}


def load_config(path: Path) -> dict[str, Any]:
    if path.exists():
        with path.open("r", encoding="utf-8") as handle:
            return yaml.safe_load(handle) or {}
    return {}


def priority_tokenizers(config: dict[str, Any]) -> list[str]:
    tokenizers = config.get("tokenizers", {})
    raw = tokenizers.get("raw", ["raw_overlap_k3", "raw_overlap_k6", "raw_nonoverlap_k6"])
    cds = tokenizers.get(
        "cds",
        [
            "cds_codon",
            "cds_frame_k6",
            "cds_nonoverlap_k3_offset0",
            "cds_nonoverlap_k3_offset1",
            "cds_nonoverlap_k3_offset2",
        ],
    )
    return list(raw) + list(cds)


def mvp_sequence_rows(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in panel.itertuples(index=False):
        for protein, seq_col, hash_col in [("HA", "ha_sequence", "ha_sha256"), ("NA", "na_sequence", "na_sha256")]:
            rows.append(
                {
                    "internal_sequence_id": f"{row.internal_strain_id}_{protein}",
                    "sequence_sha256": getattr(row, hash_col),
                    "subtype": row.subtype,
                    "protein": protein,
                    "protein_subtype": f"{protein}-{row.subtype}",
                    "year": row.year,
                    "sequence": getattr(row, seq_col),
                }
            )
    return pd.DataFrame(rows)


def cds_sequence_rows(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in panel.itertuples(index=False):
        rows.append(
            {
                "internal_sequence_id": row.internal_sequence_id,
                "sequence_sha256": row.sequence_sha256,
                "subtype": row.subtype,
                "protein": row.protein,
                "protein_subtype": f"{row.protein}-{row.subtype}",
                "year": row.year,
                "sequence": row.refined_sequence,
            }
        )
    return pd.DataFrame(rows)


def tokenizer_spec(tokenizer: str) -> dict[str, Any]:
    if tokenizer == "raw_overlap_k3":
        return {"tokenizer": tokenizer, "dataset": "raw", "family": "overlapping_kmer", "k": 3, "offset": 0}
    if tokenizer == "raw_overlap_k6":
        return {"tokenizer": tokenizer, "dataset": "raw", "family": "overlapping_kmer", "k": 6, "offset": 0}
    if tokenizer == "raw_nonoverlap_k6":
        return {"tokenizer": tokenizer, "dataset": "raw", "family": "non_overlapping_kmer", "k": 6, "offset": 0}
    if tokenizer == "cds_codon":
        return {"tokenizer": tokenizer, "dataset": "cds", "family": "codon", "k": 3, "offset": 0}
    if tokenizer == "cds_frame_k6":
        return {"tokenizer": tokenizer, "dataset": "cds", "family": "frame_aware_kmer", "k": 6, "offset": 0}
    if tokenizer.startswith("cds_nonoverlap_k3_offset"):
        offset = int(tokenizer.rsplit("offset", 1)[1])
        return {"tokenizer": tokenizer, "dataset": "cds", "family": "non_overlapping_kmer", "k": 3, "offset": offset}
    raise ValueError(f"Unsupported Phase 6 tokenizer: {tokenizer}")


def tokens_for_spec(seq: str, spec: dict[str, Any]) -> list[str]:
    family = spec["family"]
    k = int(spec["k"])
    offset = int(spec["offset"])
    if family == "codon":
        return [token for token, _, _ in codon_tokenize_with_positions(seq)]
    if family == "overlapping_kmer":
        return [token for token, _, _ in fixed_kmer_tokenize(seq, k=k, step=1)]
    if family == "non_overlapping_kmer":
        return [token for token, _, _ in non_overlapping_kmer_tokenize_with_positions(seq, k=k, offset=offset)]
    if family == "frame_aware_kmer":
        return [token for token, _, _ in frame_aware_kmer_tokenize_with_positions(seq, k=k, frame=offset)]
    raise ValueError(f"Unknown tokenizer family: {family}")


def summarize_bootstrap_metrics(bootstrap_metrics: pd.DataFrame) -> pd.DataFrame:
    metric_map = {
        "mean_token_entropy_bits": "token_entropy_bits",
        "mean_effective_vocab_size": "effective_vocab_size",
        "mean_cpg_token_fraction": "cpg_token_fraction",
        "mean_upa_token_fraction": "upa_token_fraction",
    }
    rows = []
    for key, group in bootstrap_metrics.groupby(["tokenizer", "protein_subtype"], dropna=False):
        tokenizer, protein_subtype = key
        for col, metric in metric_map.items():
            summary = summarize_bootstrap_ci(group[col])
            rows.append(
                {
                    "tokenizer": tokenizer,
                    "protein_subtype": protein_subtype,
                    "metric": metric,
                    **summary,
                }
            )
    return pd.DataFrame(rows)


def summarize_js(js_boot: pd.DataFrame) -> pd.DataFrame:
    non_diag = js_boot[js_boot["group_a"] != js_boot["group_b"]].copy()
    non_diag["comparison_key"] = non_diag.apply(lambda r: " vs ".join(sorted([r["group_a"], r["group_b"]])), axis=1)
    non_diag = non_diag.drop_duplicates(["bootstrap_id", "tokenizer", "comparison_key"])
    rows = []
    target = non_diag[non_diag["comparison"].isin(TARGET_COMPARISONS)].copy()
    for key, group in target.groupby(["tokenizer", "comparison"], dropna=False):
        tokenizer, comparison = key
        rows.append({"tokenizer": tokenizer, "comparison": comparison, **summarize_bootstrap_ci(group["js_distance"])})
    for key, group in non_diag.groupby(["bootstrap_id", "tokenizer"], dropna=False):
        bootstrap_id, tokenizer = key
        rows.append(
            {
                "tokenizer": tokenizer,
                "comparison": "mean_pairwise_nonidentical",
                "bootstrap_id": bootstrap_id,
                "value": float(group["js_distance"].mean()),
            }
        )
    overall = pd.DataFrame([row for row in rows if "value" in row])
    summary_rows = [row for row in rows if "value" not in row]
    if not overall.empty:
        for tokenizer, group in overall.groupby("tokenizer", dropna=False):
            summary_rows.append({"tokenizer": tokenizer, "comparison": "mean_pairwise_nonidentical", **summarize_bootstrap_ci(group["value"])})
    return pd.DataFrame(summary_rows)


def summarize_jaccard(jaccard: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, group in jaccard.groupby(["tokenizer", "protein_subtype"], dropna=False):
        tokenizer, protein_subtype = key
        rows.append({"tokenizer": tokenizer, "protein_subtype": protein_subtype, **summarize_bootstrap_ci(group["jaccard_vs_global_top"])})
    return pd.DataFrame(rows)


def temporal_window_counts(stats: pd.DataFrame, strategy: str) -> tuple[pd.DataFrame, pd.Series]:
    bins = temporal_bins(stats, "year", strategy=strategy)
    temp = stats.copy()
    temp["time_window"] = bins
    counts = (
        temp.groupby(["tokenizer", "time_window", "protein_subtype"], dropna=False)["internal_sequence_id"]
        .nunique()
        .reset_index(name="n_sequences")
    )
    return counts, bins


def usable_temporal_windows(counts: pd.DataFrame, min_per_group: int) -> set[tuple[str, str]]:
    usable = set()
    for key, group in counts.groupby(["tokenizer", "time_window"], dropna=False):
        if group["protein_subtype"].nunique() == 4 and group["n_sequences"].min() >= min_per_group:
            usable.add(key)
    return usable


def build_temporal_token_distributions(
    raw_rows: pd.DataFrame,
    cds_rows: pd.DataFrame,
    tokenizers: list[str],
    time_bins_by_id: pd.DataFrame,
    usable: set[tuple[str, str]],
) -> pd.DataFrame:
    rows = []
    time_lookup = dict(zip(time_bins_by_id["internal_sequence_id"], time_bins_by_id["time_window"], strict=False))
    row_sets = {"raw": raw_rows, "cds": cds_rows}
    for tokenizer in tokenizers:
        spec = tokenizer_spec(tokenizer)
        seq_rows = row_sets[spec["dataset"]]
        counters: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
        print(f"Building temporal token distributions for {tokenizer}")
        for row in seq_rows.itertuples(index=False):
            time_window = time_lookup.get(row.internal_sequence_id, "unknown")
            if (tokenizer, time_window) not in usable:
                continue
            tokens = tokens_for_spec(row.sequence, spec)
            counters[(tokenizer, time_window, row.protein, row.subtype)].update(tokens)
        for (tokenizer_name, time_window, protein, subtype), counter in counters.items():
            total = sum(counter.values())
            for token, count in counter.items():
                rows.append(
                    {
                        "tokenizer": tokenizer_name,
                        "time_window": time_window,
                        "protein": protein,
                        "subtype": subtype,
                        "protein_subtype": f"{protein}-{subtype}",
                        "token": token,
                        "count": int(count),
                        "frequency": count / total if total else 0.0,
                    }
                )
    return pd.DataFrame(rows)


def temporal_js_summary(temporal_dist: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if temporal_dist.empty:
        return pd.DataFrame(rows)
    for key, group in temporal_dist.groupby(["tokenizer", "time_window"], dropna=False):
        tokenizer, time_window = key
        pivot = group.pivot_table(index="protein_subtype", columns="token", values="frequency", aggfunc="sum", fill_value=0.0)
        pivot = pivot.reindex([g for g in GROUP_ORDER if g in pivot.index])
        groups = list(pivot.index)
        values = pivot.to_numpy(dtype=float)
        distances = []
        for i, group_a in enumerate(groups):
            for j, group_b in enumerate(groups):
                if i >= j:
                    continue
                comparison = f"{group_a} vs {group_b}"
                distance = float(jensenshannon(values[i], values[j], base=2.0))
                distances.append(distance)
                rows.append(
                    {
                        "row_type": "js_distance",
                        "tokenizer": tokenizer,
                        "time_window": time_window,
                        "protein_subtype": "",
                        "comparison": comparison,
                        "n_sequences": np.nan,
                        "mean_token_entropy_bits": np.nan,
                        "mean_effective_vocab_size": np.nan,
                        "mean_cpg_token_fraction": np.nan,
                        "mean_upa_token_fraction": np.nan,
                        "js_distance": distance,
                    }
                )
        if distances:
            rows.append(
                {
                    "row_type": "js_mean_pairwise",
                    "tokenizer": tokenizer,
                    "time_window": time_window,
                    "protein_subtype": "",
                    "comparison": "mean_pairwise_nonidentical",
                    "n_sequences": np.nan,
                    "mean_token_entropy_bits": np.nan,
                    "mean_effective_vocab_size": np.nan,
                    "mean_cpg_token_fraction": np.nan,
                    "mean_upa_token_fraction": np.nan,
                    "js_distance": float(np.mean(distances)),
                }
            )
    return pd.DataFrame(rows)


def build_temporal_summary(
    stats: pd.DataFrame,
    strategy: str,
    min_per_group: int,
    raw_rows: pd.DataFrame,
    cds_rows: pd.DataFrame,
    tokenizers: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, set[tuple[str, str]]]:
    counts, bins = temporal_window_counts(stats, strategy=strategy)
    usable = usable_temporal_windows(counts, min_per_group=min_per_group)
    metrics = token_metrics_by_time_window(stats, None, bins)
    metrics = metrics[metrics.apply(lambda r: (r["tokenizer"], r["time_window"]) in usable, axis=1)].copy()
    metrics["row_type"] = "group_metric"
    metrics["comparison"] = ""
    metrics["js_distance"] = np.nan

    time_bins_by_id = stats[["internal_sequence_id"]].copy()
    time_bins_by_id["time_window"] = bins.to_numpy()
    time_bins_by_id = time_bins_by_id.drop_duplicates("internal_sequence_id")
    temporal_dist = build_temporal_token_distributions(raw_rows, cds_rows, tokenizers, time_bins_by_id, usable)
    temporal_js = temporal_js_summary(temporal_dist)

    cols = [
        "row_type",
        "tokenizer",
        "time_window",
        "protein_subtype",
        "comparison",
        "n_sequences",
        "mean_token_entropy_bits",
        "mean_effective_vocab_size",
        "mean_cpg_token_fraction",
        "mean_upa_token_fraction",
        "js_distance",
    ]
    temporal_summary = pd.concat([metrics[cols], temporal_js[cols]], ignore_index=True) if not temporal_js.empty else metrics[cols]
    return temporal_summary, temporal_dist, usable


def minmax_scale(series: pd.Series, higher_is_better: bool = True) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce").astype(float)
    lo = values.min()
    hi = values.max()
    if not np.isfinite(lo) or not np.isfinite(hi) or hi == lo:
        scaled = pd.Series(1.0, index=series.index)
    else:
        scaled = (values - lo) / (hi - lo)
    return scaled if higher_is_better else 1.0 - scaled


def build_robustness_ranking(
    js_summary: pd.DataFrame,
    jaccard_summary: pd.DataFrame,
    metric_summary: pd.DataFrame,
    tokenizer_summary: pd.DataFrame,
) -> pd.DataFrame:
    js = js_summary[js_summary["comparison"] == "mean_pairwise_nonidentical"].copy()
    js = js.rename(columns={"mean": "mean_js_distance", "ci_width": "js_ci_width", "std": "js_std"})
    jac = jaccard_summary.groupby("tokenizer", dropna=False).agg(mean_top_token_jaccard=("mean", "mean"), top_token_jaccard_ci_width=("ci_width", "mean")).reset_index()
    entropy = metric_summary[metric_summary["metric"] == "token_entropy_bits"].groupby("tokenizer", dropna=False).agg(entropy_ci_width=("ci_width", "mean")).reset_index()
    coverage = tokenizer_summary[["tokenizer", "n_sequences"]].drop_duplicates()
    ranking = js[["tokenizer", "mean_js_distance", "js_ci_width", "js_std"]].merge(jac, on="tokenizer", how="left").merge(entropy, on="tokenizer", how="left").merge(coverage, on="tokenizer", how="left")
    ranking["score_js_distance"] = minmax_scale(ranking["mean_js_distance"], higher_is_better=True)
    ranking["score_js_stability"] = minmax_scale(ranking["js_ci_width"], higher_is_better=False)
    ranking["score_top_token_stability"] = minmax_scale(ranking["mean_top_token_jaccard"], higher_is_better=True)
    ranking["score_entropy_stability"] = minmax_scale(ranking["entropy_ci_width"], higher_is_better=False)
    ranking["score_coverage"] = ranking["n_sequences"] / ranking["n_sequences"].max()
    ranking["robustness_score"] = ranking[
        [
            "score_js_distance",
            "score_js_stability",
            "score_top_token_stability",
            "score_entropy_stability",
            "score_coverage",
        ]
    ].mean(axis=1)
    ranking = ranking.sort_values("robustness_score", ascending=False).reset_index(drop=True)
    ranking.insert(0, "rank", np.arange(1, len(ranking) + 1))
    return ranking


def write_tables(
    metric_summary: pd.DataFrame,
    js_summary: pd.DataFrame,
    jaccard_summary: pd.DataFrame,
    temporal_summary: pd.DataFrame,
    ranking: pd.DataFrame,
) -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    metric_summary.to_csv(TABLE_DIR / "phase6_bootstrap_metric_summary.csv", index=False)
    js_summary.to_csv(TABLE_DIR / "phase6_js_distance_stability.csv", index=False)
    jaccard_summary.to_csv(TABLE_DIR / "phase6_top_token_jaccard_stability.csv", index=False)
    temporal_summary.to_csv(TABLE_DIR / "phase6_temporal_token_summary.csv", index=False)
    ranking.to_csv(TABLE_DIR / "phase6_tokenizer_robustness_ranking.csv", index=False)


def plot_js_ci(js_summary: pd.DataFrame) -> None:
    df = js_summary[js_summary["comparison"] == "mean_pairwise_nonidentical"].sort_values("mean", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5.3))
    err = [df["mean"] - df["ci_lower"], df["ci_upper"] - df["mean"]]
    ax.errorbar(df["mean"], np.arange(len(df)), xerr=err, fmt="o", color="#3b6ea8", ecolor="#9ebad8", capsize=3)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["tokenizer"])
    ax.set_xlabel("mean pairwise JS distance with 95% CI")
    ax.set_title("Figure 24. Bootstrap JS distance stability")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fig24_bootstrap_js_distance_ci.png", dpi=220)
    plt.close(fig)


def plot_entropy_stability(metric_summary: pd.DataFrame) -> None:
    df = metric_summary[metric_summary["metric"] == "token_entropy_bits"].copy()
    tokenizers = df["tokenizer"].drop_duplicates().tolist()
    groups = [g for g in GROUP_ORDER if g in set(df["protein_subtype"])]
    fig, ax = plt.subplots(figsize=(12.5, 5.5))
    x = np.arange(len(tokenizers))
    width = min(0.18, 0.8 / len(groups))
    cmap = plt.get_cmap("tab10")
    for idx, group in enumerate(groups):
        part = df[df["protein_subtype"] == group].set_index("tokenizer")
        means = [part.loc[tokenizer, "mean"] if tokenizer in part.index else np.nan for tokenizer in tokenizers]
        lows = [part.loc[tokenizer, "ci_lower"] if tokenizer in part.index else np.nan for tokenizer in tokenizers]
        highs = [part.loc[tokenizer, "ci_upper"] if tokenizer in part.index else np.nan for tokenizer in tokenizers]
        yerr = [np.asarray(means) - np.asarray(lows), np.asarray(highs) - np.asarray(means)]
        ax.bar(x + (idx - (len(groups) - 1) / 2) * width, means, width=width, yerr=yerr, capsize=2, label=group, color=cmap(idx))
    ax.set_xticks(x)
    ax.set_xticklabels(tokenizers, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel("bootstrap mean entropy (bits)")
    ax.set_title("Figure 25. Token entropy stability by group")
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fig25_token_entropy_stability.png", dpi=220)
    plt.close(fig)


def plot_jaccard(jaccard_summary: pd.DataFrame) -> None:
    df = jaccard_summary.groupby("tokenizer", dropna=False).agg(mean=("mean", "mean"), ci_lower=("ci_lower", "mean"), ci_upper=("ci_upper", "mean")).reset_index().sort_values("mean", ascending=False)
    fig, ax = plt.subplots(figsize=(10, 5.3))
    err = [df["mean"] - df["ci_lower"], df["ci_upper"] - df["mean"]]
    ax.errorbar(df["mean"], np.arange(len(df)), xerr=err, fmt="o", color="#4c8b57", ecolor="#a9c9ad", capsize=3)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["tokenizer"])
    ax.set_xlabel("mean Jaccard vs global top 20")
    ax.set_title("Figure 26. Top-token Jaccard stability")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fig26_top_token_jaccard_stability.png", dpi=220)
    plt.close(fig)


def plot_temporal_entropy(temporal_summary: pd.DataFrame) -> None:
    df = temporal_summary[temporal_summary["row_type"] == "group_metric"].copy()
    if df.empty:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No temporal windows passed the minimum group threshold.", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(FIGURE_DIR / "fig27_temporal_token_entropy.png", dpi=220)
        plt.close(fig)
        return
    order = ["pre-2009", "2009-2014", "2015-2019", "2020+"]
    df["time_window"] = pd.Categorical(df["time_window"], categories=order, ordered=True)
    plot_df = df.groupby(["tokenizer", "time_window"], dropna=False)["mean_token_entropy_bits"].mean().reset_index()
    fig, ax = plt.subplots(figsize=(11, 5.5))
    for tokenizer, group in plot_df.groupby("tokenizer", dropna=False):
        group = group.sort_values("time_window")
        ax.plot(group["time_window"].astype(str), group["mean_token_entropy_bits"], marker="o", linewidth=1.4, label=tokenizer)
    ax.set_ylabel("mean token entropy (bits)")
    ax.set_xlabel("time window")
    ax.set_title("Figure 27. Temporal token entropy")
    ax.legend(fontsize=7, ncols=2)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fig27_temporal_token_entropy.png", dpi=220)
    plt.close(fig)


def plot_ranking(ranking: pd.DataFrame) -> None:
    df = ranking.sort_values("robustness_score", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 5.3))
    ax.barh(df["tokenizer"], df["robustness_score"], color="#6f5aa7")
    ax.set_xlabel("composite robustness score")
    ax.set_title("Figure 28. Tokenizer robustness ranking")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "fig28_tokenizer_robustness_ranking.png", dpi=220)
    plt.close(fig)


def write_figures(
    js_summary: pd.DataFrame,
    metric_summary: pd.DataFrame,
    jaccard_summary: pd.DataFrame,
    temporal_summary: pd.DataFrame,
    ranking: pd.DataFrame,
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plot_js_ci(js_summary)
    plot_entropy_stability(metric_summary)
    plot_jaccard(jaccard_summary)
    plot_temporal_entropy(temporal_summary)
    plot_ranking(ranking)


def markdown_table(df: pd.DataFrame, max_rows: int = 12) -> str:
    shown = df.head(max_rows)
    if shown.empty:
        return "_No rows._"
    cols = list(shown.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in shown.iterrows():
        vals = []
        for col in cols:
            val = row[col]
            if pd.isna(val):
                vals.append("")
            elif isinstance(val, float):
                vals.append(f"{val:.6g}")
            else:
                vals.append(str(val))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Table truncated to {max_rows} of {len(df)} rows._")
    return "\n".join(lines)


def write_report(
    tokenizers: list[str],
    n_bootstraps: int,
    max_per_group: int,
    usable_windows: set[tuple[str, str]],
    js_summary: pd.DataFrame,
    jaccard_summary: pd.DataFrame,
    temporal_summary: pd.DataFrame,
    ranking: pd.DataFrame,
) -> None:
    js_overall = js_summary[js_summary["comparison"] == "mean_pairwise_nonidentical"].sort_values("mean", ascending=False)
    jac_overall = jaccard_summary.groupby("tokenizer", dropna=False)["mean"].mean().reset_index(name="mean_top_token_jaccard").sort_values("mean_top_token_jaccard", ascending=False)
    temporal_js = temporal_summary[temporal_summary["row_type"] == "js_mean_pairwise"][
        ["tokenizer", "time_window", "comparison", "js_distance"]
    ].sort_values(["tokenizer", "time_window"])
    temporal_windows = sorted({window for _, window in usable_windows})
    if temporal_windows:
        temporal_note = ", ".join(temporal_windows)
    else:
        temporal_note = "No temporal window passed the minimum per-group threshold."

    report = f"""# FluGenome3D Phase 6 tokenization stability report

This phase evaluates the stability of deterministic tokenization metrics under stratified bootstrap and temporal grouping. It is a descriptive robustness audit before learned BPE or GROVER tokenizers.

## Why robustness before learned tokenizers

Phase 5 compared deterministic tokenization choices. Phase 6 checks whether those descriptive differences are stable under resampling, whether group distances are consistent through time windows, and whether top-token enrichment is robust to sampling variation. This avoids treating a single full-panel estimate as if it were automatically stable.

## Tokenizers evaluated

{", ".join(tokenizers)}

## Bootstrap design

- Bootstrap replicates: {n_bootstraps}
- Sequence metric bootstrap: stratified by HA-H1N1, NA-H1N1, HA-H3N2 and NA-H3N2.
- Maximum sampled rows per group per bootstrap: {max_per_group}
- JS and top-token bootstrap: aggregate group token distributions were resampled within each tokenizer/group because Phase 5 does not store ordered token lists per sequence.

## JS distance stability

Mean pairwise Jensen-Shannon distance is summarized with bootstrap confidence intervals. Higher values indicate stronger descriptive separation of aggregate token distributions, not predictive performance.

{markdown_table(js_overall[["tokenizer", "mean", "std", "ci_lower", "ci_upper", "ci_width"]])}

## Top-token stability

Top-token stability is measured as Jaccard overlap between bootstrap top-20 enriched tokens and the full-panel top-20 reference for the same tokenizer/group.

{markdown_table(jac_overall)}

## Temporal stability

Temporal bins used: {temporal_note}

Temporal metrics are reported only for tokenizer/window combinations where all four protein-subtype groups pass the minimum group count threshold.

{markdown_table(temporal_summary.head(12))}

Temporal mean pairwise JS distances:

{markdown_table(temporal_js.head(16), max_rows=16)}

## Robustness ranking

The ranking is a composite descriptive score combining higher JS distance, narrower JS confidence intervals, higher top-token Jaccard stability, narrower entropy intervals and coverage.

{markdown_table(ranking[["rank", "tokenizer", "robustness_score", "mean_js_distance", "js_ci_width", "mean_top_token_jaccard", "n_sequences"]])}

## Figures

- `results/figures/fig24_bootstrap_js_distance_ci.png`
- `results/figures/fig25_token_entropy_stability.png`
- `results/figures/fig26_top_token_jaccard_stability.png`
- `results/figures/fig27_temporal_token_entropy.png`
- `results/figures/fig28_tokenizer_robustness_ranking.png`

## What this phase can support

- We evaluate the stability of deterministic tokenization metrics under stratified bootstrap.
- k=6 tokenizers can be described as more or less stable based on bootstrap confidence intervals and top-token overlap.
- Top-token enrichment is assessed descriptively through Jaccard stability.
- Temporal stability is evaluated only where metadata supports it.
- This phase identifies robust deterministic baselines before learned tokenizers.

## What this phase does not support

- Stable tokens are not antigenic markers.
- Token stability does not predict evolution, escape, vaccine relevance or fitness.
- Temporal token changes do not imply selection.
- This phase does not validate GROVER or any learned tokenizer.

## Recommendation for Phase 7

Proceed to a learned-tokenizer audit only if the goal is methodological comparison. A conservative Phase 7 would introduce local BPE as an optional learned baseline and compare it against the robust deterministic baselines identified here, still without GROVER, prediction or biological efficacy claims.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Phase 6 tokenization stability audit.")
    parser.add_argument("--config", default=str(DEFAULT_CONFIG))
    args = parser.parse_args()

    config = load_config(Path(args.config))
    tokenizers = priority_tokenizers(config)
    n_bootstraps = int(config.get("n_bootstraps", 100))
    seed = int(config.get("random_seed", 42))
    max_per_group = int(config.get("max_per_group_per_bootstrap", 1000))
    max_tokens = int(config.get("aggregate_bootstrap_max_tokens_per_group", 200_000))
    top_n = int(config.get("top_n_tokens", 20))
    min_top_count = int(config.get("min_top_token_count", 100))
    min_top_freq = float(config.get("min_top_token_frequency", 1e-5))
    temporal_config = config.get("temporal", {})
    temporal_strategy = str(temporal_config.get("strategy", "fixed"))
    min_per_window = int(temporal_config.get("min_per_group_per_window", 50))

    STABILITY_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    raw_stats = pd.read_parquet(TOKEN_DIR / "mvp_token_stats.parquet")
    cds_stats = pd.read_parquet(TOKEN_DIR / "mvp_cds_token_stats.parquet")
    stats = pd.concat([raw_stats, cds_stats], ignore_index=True)
    stats = stats[stats["tokenizer"].isin(tokenizers)].reset_index(drop=True)
    distributions = pd.read_parquet(TOKEN_DIR / "group_token_distributions.parquet")
    distributions = distributions[distributions["tokenizer"].isin(tokenizers)].copy()
    distributions = distributions[distributions["token"].astype(str).str.len() <= 6].copy()

    bootstrap_frames = []
    for offset, tokenizer in enumerate(tokenizers):
        print(f"Bootstrapping sequence metrics for {tokenizer}")
        bootstrap_frames.append(
            bootstrap_token_metrics(
                stats,
                None,
                tokenizer,
                n_bootstraps=n_bootstraps,
                group_col="protein_subtype",
                n_per_group=max_per_group,
                random_state=seed + offset,
            )
        )
    bootstrap_metrics = pd.concat(bootstrap_frames, ignore_index=True)
    bootstrap_metrics.to_parquet(STABILITY_DIR / "bootstrap_metrics.parquet", index=False)
    metric_summary = summarize_bootstrap_metrics(bootstrap_metrics)

    print("Bootstrapping JS distances")
    js_boot = bootstrap_js_distances(distributions, n_bootstraps=n_bootstraps, random_state=seed, max_tokens_per_group=max_tokens)
    js_boot.to_parquet(STABILITY_DIR / "bootstrap_js_distances.parquet", index=False)
    js_summary = summarize_js(js_boot)

    print("Bootstrapping top tokens")
    global_top = token_enrichment_by_group(distributions, top_n=top_n, min_count=min_top_count, min_group_frequency=min_top_freq)
    global_top = global_top[global_top["token"].astype(str).str.len() <= 6].copy()
    top_boot = bootstrap_top_tokens(
        distributions,
        n_bootstraps=n_bootstraps,
        random_state=seed + 99,
        max_tokens_per_group=max_tokens,
        top_n=top_n,
        min_count=min_top_count,
        min_group_frequency=min_top_freq,
    )
    top_boot = top_boot[top_boot["token"].astype(str).str.len() <= 6].copy()
    top_boot.to_parquet(STABILITY_DIR / "bootstrap_top_tokens.parquet", index=False)
    global_top.to_parquet(STABILITY_DIR / "global_top_tokens_reference.parquet", index=False)
    jaccard = top_token_jaccard_stability(top_boot, global_top)
    jaccard.to_parquet(STABILITY_DIR / "top_token_jaccard_by_bootstrap.parquet", index=False)
    jaccard_summary = summarize_jaccard(jaccard)

    print("Computing temporal stability")
    raw_rows = mvp_sequence_rows(pd.read_parquet(RAW_PANEL))
    cds_rows = cds_sequence_rows(pd.read_parquet(CDS_PANEL))
    temporal_summary, temporal_dist, usable_windows = build_temporal_summary(
        stats,
        strategy=temporal_strategy,
        min_per_group=min_per_window,
        raw_rows=raw_rows,
        cds_rows=cds_rows,
        tokenizers=tokenizers,
    )
    temporal_summary.to_parquet(STABILITY_DIR / "temporal_metrics.parquet", index=False)
    temporal_dist.to_parquet(STABILITY_DIR / "temporal_group_token_distributions.parquet", index=False)

    tokenizer_summary = pd.read_csv(TABLE_DIR / "phase5_tokenizer_summary.csv")
    tokenizer_summary = tokenizer_summary[tokenizer_summary["tokenizer"].isin(tokenizers)].copy()
    ranking = build_robustness_ranking(js_summary, jaccard_summary, metric_summary, tokenizer_summary)

    write_tables(metric_summary, js_summary, jaccard_summary, temporal_summary, ranking)
    write_figures(js_summary, metric_summary, jaccard_summary, temporal_summary, ranking)
    write_report(tokenizers, n_bootstraps, max_per_group, usable_windows, js_summary, jaccard_summary, temporal_summary, ranking)

    print(f"Wrote local Phase 6 outputs to {STABILITY_DIR}")
    print(f"Wrote Phase 6 tables to {TABLE_DIR}")
    print(f"Wrote Phase 6 figures to {FIGURE_DIR}")
    print(f"Wrote report to {REPORT_PATH}")


if __name__ == "__main__":
    main()
