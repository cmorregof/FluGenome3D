#!/usr/bin/env python
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from flugenome3d.token_stats import (
    group_token_distribution,
    jensen_shannon_distance_between_groups,
    token_enrichment_by_group,
)
from flugenome3d.tokenization import (
    codon_tokenize_with_positions,
    fixed_kmer_tokenize,
    frame_aware_kmer_tokenize_with_positions,
    non_overlapping_kmer_tokenize_with_positions,
    token_contains_cpg,
    token_contains_upa,
    token_crosses_codon_boundary,
)


MVP_PANEL = Path("data/processed/panels/mvp_panel.parquet")
CDS_REFINED_PANEL = Path("data/processed/panels/mvp_cds_refined_panel.parquet")
TOKEN_DIR = Path("data/processed/tokenization")
TABLE_DIR = Path("results/tables")
FIGURE_DIR = Path("results/figures")
REPORT_PATH = Path("reports/phase5_tokenization_audit_report.md")
GROUP_ORDER = ["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"]


def mvp_sequence_rows(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in panel.itertuples(index=False):
        for protein, seq_col, hash_col in [("HA", "ha_sequence", "ha_sha256"), ("NA", "na_sequence", "na_sha256")]:
            sequence = getattr(row, seq_col)
            rows.append(
                {
                    "internal_strain_id": row.internal_strain_id,
                    "internal_sequence_id": f"{row.internal_strain_id}_{protein}",
                    "sequence_sha256": getattr(row, hash_col),
                    "pair_sha256": row.pair_sha256,
                    "subtype": row.subtype,
                    "protein": protein,
                    "protein_subtype": f"{protein}-{row.subtype}",
                    "year": getattr(row, "year", np.nan),
                    "sequence": sequence,
                }
            )
    return pd.DataFrame(rows)


def cds_sequence_rows(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in panel.itertuples(index=False):
        sequence = getattr(row, "refined_sequence")
        rows.append(
            {
                "internal_sequence_id": row.internal_sequence_id,
                "sequence_sha256": getattr(row, "sequence_sha256", ""),
                "pair_sha256": getattr(row, "pair_sha256", ""),
                "subtype": row.subtype,
                "protein": row.protein,
                "protein_subtype": f"{row.protein}-{row.subtype}",
                "year": getattr(row, "year", np.nan),
                "cds_status": getattr(row, "status", getattr(row, "cds_status", "")),
                "rescue_method": getattr(row, "rescue_method", ""),
                "sequence": sequence,
            }
        )
    return pd.DataFrame(rows)


def tokenizer_specs() -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for k in [3, 4, 5, 6]:
        specs.append(
            {
                "tokenizer": f"raw_overlap_k{k}",
                "dataset": "mvp_panel",
                "source": "raw_nucleotide",
                "family": "overlapping_kmer",
                "mode": "overlapping",
                "k": k,
                "offset": 0,
                "frame": 0,
                "boundary_frame": 0,
                "boundary_context": "position0_proxy_not_cds",
            }
        )
    for k in [3, 4, 5, 6]:
        specs.append(
            {
                "tokenizer": f"raw_nonoverlap_k{k}",
                "dataset": "mvp_panel",
                "source": "raw_nucleotide",
                "family": "non_overlapping_kmer",
                "mode": "non_overlapping",
                "k": k,
                "offset": 0,
                "frame": 0,
                "boundary_frame": 0,
                "boundary_context": "position0_proxy_not_cds",
            }
        )
    specs.append(
        {
            "tokenizer": "cds_codon",
            "dataset": "mvp_cds_refined_panel",
            "source": "cds_refined",
            "family": "codon",
            "mode": "codon",
            "k": 3,
            "offset": 0,
            "frame": 0,
            "boundary_frame": 0,
            "boundary_context": "cds_refined_codon_aligned",
        }
    )
    for k in [3, 6]:
        specs.append(
            {
                "tokenizer": f"cds_frame_k{k}",
                "dataset": "mvp_cds_refined_panel",
                "source": "cds_refined",
                "family": "frame_aware_kmer",
                "mode": "frame_aware",
                "k": k,
                "offset": 0,
                "frame": 0,
                "boundary_frame": 0,
                "boundary_context": "cds_refined_frame0",
            }
        )
    for offset in [0, 1, 2]:
        specs.append(
            {
                "tokenizer": f"cds_nonoverlap_k3_offset{offset}",
                "dataset": "mvp_cds_refined_panel",
                "source": "cds_refined",
                "family": "non_overlapping_kmer",
                "mode": "non_overlapping",
                "k": 3,
                "offset": offset,
                "frame": offset,
                "boundary_frame": 0,
                "boundary_context": "cds_refined_offset_sensitivity",
            }
        )
    return specs


def tokens_with_positions(seq: str, spec: dict[str, Any]) -> list[tuple[str, int, int]]:
    family = spec["family"]
    k = int(spec["k"])
    offset = int(spec.get("offset", 0))
    if family == "codon":
        return codon_tokenize_with_positions(seq)
    if family == "overlapping_kmer":
        return fixed_kmer_tokenize(seq, k=k, step=1)
    if family == "non_overlapping_kmer":
        return non_overlapping_kmer_tokenize_with_positions(seq, k=k, offset=offset)
    if family == "frame_aware_kmer":
        return frame_aware_kmer_tokenize_with_positions(seq, k=k, frame=int(spec.get("frame", offset)))
    raise ValueError(f"Unknown tokenizer family: {family}")


def entropy_from_counter(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total <= 0:
        return 0.0
    return float(-sum((count / total) * np.log2(count / total) for count in counts.values() if count > 0))


def summarize_sequence_tokens(row: Any, spec: dict[str, Any], toks_pos: list[tuple[str, int, int]]) -> dict[str, Any]:
    tokens = [token for token, _, _ in toks_pos]
    counts = Counter(tokens)
    n_tokens = len(tokens)
    seq_len = len(row.sequence) if isinstance(row.sequence, str) else 0
    boundary_frame = int(spec.get("boundary_frame", 0))
    n_cpg = sum(token_contains_cpg(token) for token in tokens)
    n_upa = sum(token_contains_upa(token) for token in tokens)
    n_crossing = sum(token_crosses_codon_boundary(start, end, frame=boundary_frame) for _, start, end in toks_pos)
    entropy_bits = entropy_from_counter(counts)
    return {
        "internal_sequence_id": row.internal_sequence_id,
        "sequence_sha256": getattr(row, "sequence_sha256", ""),
        "tokenizer": spec["tokenizer"],
        "dataset": spec["dataset"],
        "source": spec["source"],
        "family": spec["family"],
        "mode": spec["mode"],
        "k": spec["k"],
        "offset": spec["offset"],
        "frame": spec["frame"],
        "boundary_frame": spec["boundary_frame"],
        "boundary_context": spec["boundary_context"],
        "protein": row.protein,
        "subtype": row.subtype,
        "protein_subtype": row.protein_subtype,
        "year": getattr(row, "year", np.nan),
        "sequence_length": seq_len,
        "n_tokens": n_tokens,
        "unique_tokens": len(counts),
        "token_entropy_bits": entropy_bits,
        "effective_vocab_size": float(2**entropy_bits),
        "tokens_per_kb": n_tokens / (seq_len / 1000) if seq_len else np.nan,
        "cpg_token_fraction": n_cpg / n_tokens if n_tokens else np.nan,
        "upa_token_fraction": n_upa / n_tokens if n_tokens else np.nan,
        "codon_boundary_crossing_fraction": n_crossing / n_tokens if n_tokens else np.nan,
    }


def run_tokenizers(seq_rows: pd.DataFrame, specs: list[dict[str, Any]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    stats_rows = []
    group_counters: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for spec in specs:
        print(f"Auditing {spec['tokenizer']} on {len(seq_rows)} sequences")
        for row in seq_rows.itertuples(index=False):
            toks_pos = tokens_with_positions(row.sequence, spec)
            stats_rows.append(summarize_sequence_tokens(row, spec, toks_pos))
            tokens = [token for token, _, _ in toks_pos]
            group_counters[(spec["tokenizer"], row.protein, row.subtype)].update(tokens)

    token_rows = []
    for (tokenizer, protein, subtype), counter in group_counters.items():
        for token, count in counter.items():
            token_rows.append({"tokenizer": tokenizer, "protein": protein, "subtype": subtype, "token": token, "count": count})
    distributions = group_token_distribution(pd.DataFrame(token_rows))
    return pd.DataFrame(stats_rows), distributions


def write_public_tables(stats: pd.DataFrame, distributions: pd.DataFrame) -> dict[str, pd.DataFrame]:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    observed_vocab = distributions.groupby("tokenizer")["token"].nunique().rename("observed_vocab_size")
    summary = (
        stats.groupby(["tokenizer", "dataset", "source", "family", "mode", "k", "offset", "frame", "boundary_context"], dropna=False)
        .agg(
            n_sequences=("internal_sequence_id", "nunique"),
            total_tokens=("n_tokens", "sum"),
            mean_sequence_length=("sequence_length", "mean"),
            mean_tokens_per_sequence=("n_tokens", "mean"),
            mean_tokens_per_kb=("tokens_per_kb", "mean"),
            mean_unique_tokens=("unique_tokens", "mean"),
            mean_entropy_bits=("token_entropy_bits", "mean"),
            mean_effective_vocab_size=("effective_vocab_size", "mean"),
            mean_cpg_token_fraction=("cpg_token_fraction", "mean"),
            mean_upa_token_fraction=("upa_token_fraction", "mean"),
            mean_boundary_crossing_fraction=("codon_boundary_crossing_fraction", "mean"),
        )
        .reset_index()
    )
    summary = summary.merge(observed_vocab.reset_index(), on="tokenizer", how="left")
    summary.to_csv(TABLE_DIR / "phase5_tokenizer_summary.csv", index=False)

    entropy_by_group = (
        stats.groupby(["tokenizer", "protein", "subtype", "protein_subtype"], dropna=False)
        .agg(
            n_sequences=("internal_sequence_id", "nunique"),
            mean_entropy_bits=("token_entropy_bits", "mean"),
            median_entropy_bits=("token_entropy_bits", "median"),
            sd_entropy_bits=("token_entropy_bits", "std"),
            mean_effective_vocab_size=("effective_vocab_size", "mean"),
        )
        .reset_index()
    )
    entropy_by_group.to_csv(TABLE_DIR / "phase5_token_entropy_by_group.csv", index=False)

    vocab_rows = []
    for key, group in distributions.groupby(["tokenizer", "protein", "subtype", "protein_subtype"], dropna=False):
        tokenizer, protein, subtype, protein_subtype = key
        counts = Counter(dict(zip(group["token"], group["count"], strict=False)))
        entropy_bits = entropy_from_counter(counts)
        vocab_rows.append(
            {
                "tokenizer": tokenizer,
                "protein": protein,
                "subtype": subtype,
                "protein_subtype": protein_subtype,
                "total_tokens": int(group["count"].sum()),
                "observed_vocab_size": int(group["token"].nunique()),
                "group_token_entropy_bits": entropy_bits,
                "group_effective_vocab_size": float(2**entropy_bits),
            }
        )
    effective_vocab_by_group = pd.DataFrame(vocab_rows)
    effective_vocab_by_group.to_csv(TABLE_DIR / "phase5_effective_vocab_by_group.csv", index=False)

    cpg_upa = (
        stats.groupby(["tokenizer", "protein", "subtype", "protein_subtype"], dropna=False)
        .agg(
            n_sequences=("internal_sequence_id", "nunique"),
            mean_cpg_token_fraction=("cpg_token_fraction", "mean"),
            mean_upa_token_fraction=("upa_token_fraction", "mean"),
            median_cpg_token_fraction=("cpg_token_fraction", "median"),
            median_upa_token_fraction=("upa_token_fraction", "median"),
        )
        .reset_index()
    )
    cpg_upa.to_csv(TABLE_DIR / "phase5_cpg_upa_token_summary.csv", index=False)

    boundary = (
        stats.groupby(["tokenizer", "source", "family", "boundary_context", "protein", "subtype", "protein_subtype"], dropna=False)
        .agg(
            n_sequences=("internal_sequence_id", "nunique"),
            mean_codon_boundary_crossing_fraction=("codon_boundary_crossing_fraction", "mean"),
            median_codon_boundary_crossing_fraction=("codon_boundary_crossing_fraction", "median"),
        )
        .reset_index()
    )
    boundary.to_csv(TABLE_DIR / "phase5_codon_boundary_crossing_summary.csv", index=False)

    top_tokens = token_enrichment_by_group(distributions, top_n=10, min_count=100, min_group_frequency=1e-5)
    top_tokens["token_length"] = top_tokens["token"].astype(str).str.len()
    top_tokens["contains_cpg"] = top_tokens["token"].map(token_contains_cpg)
    top_tokens["contains_upa"] = top_tokens["token"].map(token_contains_upa)
    top_tokens = top_tokens[top_tokens["token_length"] <= 6].copy()
    top_tokens.to_csv(TABLE_DIR / "phase5_top_tokens_by_group.csv", index=False)

    js_distances = jensen_shannon_distance_between_groups(distributions)
    js_distances.to_csv(TABLE_DIR / "phase5_group_js_distances.csv", index=False)

    return {
        "summary": summary,
        "entropy": entropy_by_group,
        "effective_vocab": effective_vocab_by_group,
        "cpg_upa": cpg_upa,
        "boundary": boundary,
        "top_tokens": top_tokens,
        "js": js_distances,
    }


def grouped_bars(
    df: pd.DataFrame,
    value_col: str,
    title: str,
    ylabel: str,
    outpath: Path,
    tokenizers: list[str] | None = None,
) -> None:
    plot_df = df.copy()
    if tokenizers is not None:
        plot_df = plot_df[plot_df["tokenizer"].isin(tokenizers)].copy()
    tokenizers = list(dict.fromkeys(plot_df["tokenizer"].tolist())) if tokenizers is None else tokenizers
    groups = [group for group in GROUP_ORDER if group in set(plot_df["protein_subtype"])]
    fig, ax = plt.subplots(figsize=(max(11, 0.62 * len(tokenizers) + 4), 5.4))
    x = np.arange(len(tokenizers))
    width = min(0.18, 0.8 / max(1, len(groups)))
    cmap = plt.get_cmap("tab10")
    for idx, group in enumerate(groups):
        vals = [
            plot_df[(plot_df["tokenizer"] == tokenizer) & (plot_df["protein_subtype"] == group)][value_col].mean()
            for tokenizer in tokenizers
        ]
        ax.bar(x + (idx - (len(groups) - 1) / 2) * width, vals, width=width, label=group, color=cmap(idx))
    ax.set_xticks(x)
    ax.set_xticklabels(tokenizers, rotation=35, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_cpg_upa(cpg_upa: pd.DataFrame, outpath: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    for ax, value_col, label in [
        (axes[0], "mean_cpg_token_fraction", "CpG-containing token fraction"),
        (axes[1], "mean_upa_token_fraction", "UpA/TA-containing token fraction"),
    ]:
        tokenizers = list(dict.fromkeys(cpg_upa["tokenizer"].tolist()))
        groups = [group for group in GROUP_ORDER if group in set(cpg_upa["protein_subtype"])]
        x = np.arange(len(tokenizers))
        width = min(0.18, 0.8 / max(1, len(groups)))
        cmap = plt.get_cmap("tab10")
        for idx, group in enumerate(groups):
            vals = [
                cpg_upa[(cpg_upa["tokenizer"] == tokenizer) & (cpg_upa["protein_subtype"] == group)][value_col].mean()
                for tokenizer in tokenizers
            ]
            ax.bar(x + (idx - (len(groups) - 1) / 2) * width, vals, width=width, label=group, color=cmap(idx))
        ax.set_ylabel(label)
        ax.set_title(label)
    axes[1].set_xticks(np.arange(len(tokenizers)))
    axes[1].set_xticklabels(tokenizers, rotation=35, ha="right", fontsize=8)
    axes[0].legend(fontsize=8, ncols=2)
    fig.suptitle("Figure 20. CpG and UpA/TA token fractions")
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_js_heatmaps(js: pd.DataFrame, outpath: Path) -> None:
    tokenizers = list(dict.fromkeys(js["tokenizer"].tolist()))
    ncols = 4
    nrows = int(np.ceil(len(tokenizers) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, max(7, 3.4 * nrows)))
    axes_arr = np.atleast_1d(axes).ravel()
    vmax = max(0.001, float(js["js_distance"].max()))
    for ax, tokenizer in zip(axes_arr, tokenizers, strict=False):
        part = js[js["tokenizer"] == tokenizer]
        matrix = part.pivot_table(index="group_a", columns="group_b", values="js_distance", aggfunc="mean").reindex(
            index=GROUP_ORDER, columns=GROUP_ORDER
        )
        image = ax.imshow(matrix.values, cmap="magma", vmin=0, vmax=vmax)
        ax.set_title(tokenizer, fontsize=9)
        ax.set_xticks(range(len(GROUP_ORDER)))
        ax.set_xticklabels(GROUP_ORDER, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(GROUP_ORDER)))
        ax.set_yticklabels(GROUP_ORDER, fontsize=7)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    for ax in axes_arr[len(tokenizers) :]:
        ax.axis("off")
    fig.suptitle("Figure 22. Jensen-Shannon distances between groups")
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_top_tokens(top_tokens: pd.DataFrame, outpath: Path) -> None:
    selected = ["raw_overlap_k3", "raw_overlap_k6", "cds_codon", "cds_frame_k6"]
    plot_df = top_tokens[(top_tokens["tokenizer"].isin(selected)) & (top_tokens["rank"] <= 3)].copy()
    plot_df["label"] = plot_df["protein_subtype"] + " " + plot_df["token"]
    fig, axes = plt.subplots(2, 2, figsize=(14, 8))
    axes_arr = axes.ravel()
    for ax, tokenizer in zip(axes_arr, selected, strict=False):
        part = plot_df[plot_df["tokenizer"] == tokenizer].sort_values("enrichment_ratio", ascending=True)
        ax.barh(part["label"], part["enrichment_ratio"], color="#4c78a8")
        ax.set_title(tokenizer)
        ax.set_xlabel("enrichment ratio")
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle("Figure 23. Representative top enriched tokens by group")
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_figures(tables: dict[str, pd.DataFrame]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    grouped_bars(
        tables["entropy"],
        "mean_entropy_bits",
        "Figure 18. Token entropy by tokenizer and group",
        "mean entropy (bits)",
        FIGURE_DIR / "fig18_token_entropy_by_tokenizer.png",
    )
    grouped_bars(
        tables["effective_vocab"],
        "group_effective_vocab_size",
        "Figure 19. Effective vocabulary size by group",
        "effective vocabulary size",
        FIGURE_DIR / "fig19_effective_vocab_by_group.png",
    )
    plot_cpg_upa(tables["cpg_upa"], FIGURE_DIR / "fig20_cpg_upa_token_fraction.png")
    grouped_bars(
        tables["boundary"],
        "mean_codon_boundary_crossing_fraction",
        "Figure 21. Codon-boundary crossing fraction",
        "mean crossing fraction",
        FIGURE_DIR / "fig21_codon_boundary_crossing.png",
    )
    plot_js_heatmaps(tables["js"], FIGURE_DIR / "fig22_token_js_distance_heatmap.png")
    plot_top_tokens(tables["top_tokens"], FIGURE_DIR / "fig23_top_tokens_by_group.png")


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
            if isinstance(val, float):
                vals.append(f"{val:.6g}")
            else:
                vals.append(str(val))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Table truncated to {max_rows} of {len(df)} rows._")
    return "\n".join(lines)


def write_report(tables: dict[str, pd.DataFrame], raw_n: int, cds_n: int) -> None:
    summary = tables["summary"]
    entropy = tables["entropy"]
    effective_vocab = tables["effective_vocab"]
    cpg_upa = tables["cpg_upa"]
    boundary = tables["boundary"]
    js = tables["js"]
    top_tokens = tables["top_tokens"]

    non_diag_js = js[js["group_a"] != js["group_b"]].copy()
    best_js = non_diag_js.groupby("tokenizer", dropna=False)["js_distance"].mean().reset_index().sort_values("js_distance", ascending=False)
    best_entropy = entropy.groupby("tokenizer", dropna=False)["mean_entropy_bits"].mean().reset_index().sort_values("mean_entropy_bits", ascending=False)
    boundary_display = (
        boundary.groupby(["tokenizer", "boundary_context"], dropna=False)["mean_codon_boundary_crossing_fraction"]
        .mean()
        .reset_index()
        .sort_values(["tokenizer", "boundary_context"])
    )

    report = f"""# FluGenome3D Phase 5 tokenization audit report

This phase audits how deterministic tokenization choices segment Influenza A HA/NA sequences. It is a descriptive tokenization audit before learned BPE or GROVER tokenizers.

## Representation vs. tokenization

Phase 4 treated sequences as vector representations. Phase 5 treats tokenization itself as the object of study: vocabulary size, token entropy, CpG/UpA-containing token fractions, codon-boundary crossing and Jensen-Shannon distances between HA/NA and subtype groups.

## Datasets used

- Raw nucleotide tokenizers use `mvp_panel`: {raw_n:,} HA/NA sequences.
- CDS-aware tokenizers use `mvp_cds_refined_panel`: {cds_n:,} refined CDS sequences.
- Codon and frame-aware claims are restricted to the CDS-refined panel.

## Tokenizers audited

{markdown_table(summary[["tokenizer", "dataset", "source", "family", "mode", "k", "offset", "n_sequences", "observed_vocab_size", "total_tokens"]])}

## Entropy and effective vocabulary

Token entropy and effective vocabulary vary by tokenizer family, k, protein and subtype. Larger k generally increases the observed vocabulary available to the audit, while non-overlapping tokenizers produce fewer tokens per sequence.

{markdown_table(best_entropy.rename(columns={"mean_entropy_bits": "mean_group_entropy_bits"}))}

{markdown_table(effective_vocab[["tokenizer", "protein_subtype", "observed_vocab_size", "group_effective_vocab_size"]].head(16), max_rows=16)}

## CpG and UpA/TA token analysis

CpG-containing and UpA/TA-containing token fractions are summarized descriptively by tokenizer and group. Because records are stored as DNA alphabet, TA is used as the DNA proxy for UpA.

{markdown_table(cpg_upa[["tokenizer", "protein_subtype", "mean_cpg_token_fraction", "mean_upa_token_fraction"]].head(16), max_rows=16)}

## Codon-boundary crossing

Codon-boundary crossing is interpreted as CDS-aware only for tokenizers built on `mvp_cds_refined_panel`. For raw nucleotide tokenizers, boundary fractions are reported as position-0 segmentation diagnostics and should not be treated as coding-frame evidence.

{markdown_table(boundary_display, max_rows=18)}

## Jensen-Shannon distances

Jensen-Shannon distances compare group-level token distributions across HA-H1N1, NA-H1N1, HA-H3N2 and NA-H3N2. Higher values indicate stronger descriptive separation of aggregate token distributions, not predictive performance.

{markdown_table(best_js.rename(columns={"js_distance": "mean_pairwise_js_distance"}))}

## Top enriched tokens

Top token tables contain only short public-safe tokens: k-mers of length <= 6 or codons. These tokens are descriptive enrichment summaries, not antigenic, escape, fitness or vaccine markers.

{markdown_table(top_tokens[["tokenizer", "protein_subtype", "rank", "token", "group_frequency", "background_frequency", "enrichment_ratio", "contains_cpg", "contains_upa"]].head(16), max_rows=16)}

## Figures

- `results/figures/fig18_token_entropy_by_tokenizer.png`
- `results/figures/fig19_effective_vocab_by_group.png`
- `results/figures/fig20_cpg_upa_token_fraction.png`
- `results/figures/fig21_codon_boundary_crossing.png`
- `results/figures/fig22_token_js_distance_heatmap.png`
- `results/figures/fig23_top_tokens_by_group.png`

## What this phase can support

- We audit how deterministic tokenization choices segment HA/NA sequences.
- Token entropy and effective vocabulary vary by protein/subtype group.
- CpG/UpA-containing token fractions are summarized descriptively.
- Codon-boundary crossing is evaluated as CDS-aware only in CDS-refined sequences.
- This phase establishes transparent tokenization baselines before learned BPE/GROVER tokenizers.

## What this phase does not support

- It does not identify antigenic sites.
- It does not predict escape, viral evolution, vaccine relevance or fitness.
- It does not replace biological language models.
- CpG/UpA token fractions should not be interpreted as pathogenicity or antigenicity evidence.

## Limitations

- Raw nucleotide tokenizers do not require a validated coding frame and therefore cannot support codon-frame claims.
- CDS-aware tokenizers inherit the Phase 3 rescue criteria and should be interpreted together with CDS refinement QC.
- Token enrichment is sensitive to k, overlap mode and group composition.
- Figures summarize aggregate token distributions; they do not inspect individual biological variants.

## Recommendation for Phase 6

Use this tokenization baseline to decide whether Phase 6 should add a learned tokenizer audit. If approved, BPE should be introduced as an explicitly optional learned baseline, still separate from GROVER and still without prediction, antigenicity or vaccine claims.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    TOKEN_DIR.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    raw_panel = pd.read_parquet(MVP_PANEL)
    cds_panel = pd.read_parquet(CDS_REFINED_PANEL)
    raw_rows = mvp_sequence_rows(raw_panel)
    cds_rows = cds_sequence_rows(cds_panel)

    specs = tokenizer_specs()
    raw_specs = [spec for spec in specs if spec["dataset"] == "mvp_panel"]
    cds_specs = [spec for spec in specs if spec["dataset"] == "mvp_cds_refined_panel"]

    raw_stats, raw_dist = run_tokenizers(raw_rows, raw_specs)
    cds_stats, cds_dist = run_tokenizers(cds_rows, cds_specs)
    stats = pd.concat([raw_stats, cds_stats], ignore_index=True)
    distributions = pd.concat([raw_dist, cds_dist], ignore_index=True)

    raw_stats.to_parquet(TOKEN_DIR / "mvp_token_stats.parquet", index=False)
    cds_stats.to_parquet(TOKEN_DIR / "mvp_cds_token_stats.parquet", index=False)
    distributions.to_parquet(TOKEN_DIR / "group_token_distributions.parquet", index=False)

    tables = write_public_tables(stats, distributions)
    write_figures(tables)
    write_report(tables, raw_n=len(raw_rows), cds_n=len(cds_rows))

    print(f"Wrote local token stats to {TOKEN_DIR}")
    print(f"Wrote public Phase 5 tables to {TABLE_DIR}")
    print(f"Wrote Phase 5 figures to {FIGURE_DIR}")
    print(f"Wrote report to {REPORT_PATH}")


if __name__ == "__main__":
    main()
