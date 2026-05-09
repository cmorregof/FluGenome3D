#!/usr/bin/env python
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from flugenome3d.codon_usage import (
    ALL_CODONS,
    CODON_TO_AA,
    codon_counts,
    codon_frequencies_from_counts,
    rscu_from_counts,
    translation_qc,
)
from flugenome3d.kmer_profiles import kmer_counts, sequence_kmer_metrics
from flugenome3d.sequence_metrics import (
    DINUCLEOTIDES,
    ambiguous_fraction,
    dinucleotide_odds_ratios,
    gc_content,
    sequence_length,
)


PANEL_PATHS = {
    "smoke": Path("data/processed/panels/smoke_panel.parquet"),
    "mvp": Path("data/processed/panels/mvp_panel.parquet"),
}
METRICS_DIR = Path("data/processed/metrics")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")
REPORT_PATH = Path("reports/phase2_sequence_context_report.md")
GROUP_COLS = ["panel", "subtype", "protein"]


def panel_to_sequence_rows(panel: str, df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        for protein, seq_col, hash_col in [("HA", "ha_sequence", "ha_sha256"), ("NA", "na_sequence", "na_sha256")]:
            rows.append(
                {
                    "panel": panel,
                    "internal_strain_id": row["internal_strain_id"],
                    "internal_sequence_id": f"{row['internal_strain_id']}_{protein}",
                    "pair_sha256": row["pair_sha256"],
                    "sequence_sha256": row[hash_col],
                    "subtype": row["subtype"],
                    "protein": protein,
                    "year": row["year"],
                    "sequence": row[seq_col],
                }
            )
    return pd.DataFrame(rows)


def compute_sequence_metrics(seq_rows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    base_cols = ["panel", "internal_sequence_id", "pair_sha256", "sequence_sha256", "subtype", "protein", "year", "sequence"]
    for row in seq_rows[base_cols].itertuples(index=False):
        seq = row.sequence
        out = {
            "panel": row.panel,
            "internal_sequence_id": row.internal_sequence_id,
            "pair_sha256": row.pair_sha256,
            "sequence_sha256": row.sequence_sha256,
            "subtype": row.subtype,
            "protein": row.protein,
            "year": row.year,
        }
        ratios = dinucleotide_odds_ratios(seq)
        out.update(
            {
                "sequence_length": sequence_length(seq),
                "ambiguous_fraction": ambiguous_fraction(seq),
                "gc_content": gc_content(seq),
                "cpg_oe": ratios["CG"],
                "upa_oe": ratios["TA"],
            }
        )
        for dinuc in DINUCLEOTIDES:
            out[f"dinuc_oe_{dinuc}"] = ratios[dinuc]
        out.update(sequence_kmer_metrics(seq, (3, 4, 5)))
        rows.append(out)
    return pd.DataFrame(rows)


def compute_codon_metrics(seq_rows: pd.DataFrame, max_ambiguous_fraction: float = 0.01) -> pd.DataFrame:
    rows = []
    base_cols = ["panel", "internal_sequence_id", "pair_sha256", "sequence_sha256", "subtype", "protein", "year", "sequence"]
    for row in seq_rows[base_cols].itertuples(index=False):
        seq = row.sequence
        qc = translation_qc(seq, max_ambiguous_fraction=max_ambiguous_fraction)
        out = {
            "panel": row.panel,
            "internal_sequence_id": row.internal_sequence_id,
            "pair_sha256": row.pair_sha256,
            "sequence_sha256": row.sequence_sha256,
            "subtype": row.subtype,
            "protein": row.protein,
            "year": row.year,
        }
        out.update(qc)
        counts = codon_counts(seq) if qc["codon_reliable"] else Counter()
        freqs = codon_frequencies_from_counts(counts) if qc["codon_reliable"] else {codon: np.nan for codon in ALL_CODONS}
        rscu_values = rscu_from_counts(counts) if qc["codon_reliable"] else {codon: np.nan for codon in sorted(CODON_TO_AA) if CODON_TO_AA[codon] != "*"}
        for codon in ALL_CODONS:
            out[f"codon_count_{codon}"] = counts[codon]
            out[f"codon_freq_{codon}"] = freqs[codon]
        for codon, value in rscu_values.items():
            out[f"rscu_{codon}"] = value
        rows.append(out)
    return pd.DataFrame(rows)


def write_local_metrics(panel: str, sequence_metrics: pd.DataFrame, codon_metrics: pd.DataFrame) -> dict[str, Path]:
    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    seq_path = METRICS_DIR / f"{panel}_sequence_metrics.parquet"
    codon_path = METRICS_DIR / f"{panel}_codon_metrics.parquet"
    sequence_metrics.to_parquet(seq_path, index=False)
    codon_metrics.to_parquet(codon_path, index=False)
    return {"sequence": seq_path, "codon": codon_path}


def load_available_metrics() -> tuple[pd.DataFrame, pd.DataFrame]:
    seq_frames = [pd.read_parquet(path) for path in sorted(METRICS_DIR.glob("*_sequence_metrics.parquet"))]
    codon_frames = [pd.read_parquet(path) for path in sorted(METRICS_DIR.glob("*_codon_metrics.parquet"))]
    if not seq_frames or not codon_frames:
        raise FileNotFoundError("No Phase 2 local metrics were found.")
    return pd.concat(seq_frames, ignore_index=True), pd.concat(codon_frames, ignore_index=True)


def _summary_stats(df: pd.DataFrame, value: str) -> pd.DataFrame:
    return df.groupby(GROUP_COLS, dropna=False)[value].agg(
        n="count",
        mean="mean",
        median="median",
        std="std",
        q05=lambda x: x.quantile(0.05),
        q95=lambda x: x.quantile(0.95),
    ).reset_index()


def build_gc_cpg_upa_summary(seq_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for metric in ["sequence_length", "gc_content", "cpg_oe", "upa_oe", "ambiguous_fraction"]:
        part = _summary_stats(seq_metrics, metric)
        part.insert(3, "metric", metric)
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def build_dinucleotide_summary(seq_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dinuc in DINUCLEOTIDES:
        col = f"dinuc_oe_{dinuc}"
        part = _summary_stats(seq_metrics, col)
        part.insert(3, "dinucleotide", dinuc)
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def build_translation_qc_summary(codon_metrics: pd.DataFrame) -> pd.DataFrame:
    grouped = codon_metrics.groupby(GROUP_COLS, dropna=False)
    rows = []
    for key, group in grouped:
        panel, subtype, protein = key
        n = len(group)
        reliable = int(group["codon_reliable"].sum())
        rows.append(
            {
                "panel": panel,
                "subtype": subtype,
                "protein": protein,
                "n_sequences": n,
                "n_frame_fail": int(group["frame_fail"].sum()),
                "n_ambiguous_fail": int(group["ambiguous_fail"].sum()),
                "n_internal_stop_fail": int(group["internal_stop_fail"].sum()),
                "n_translation_fail": int(group["translation_fail"].sum()),
                "n_codon_reliable": reliable,
                "reliable_fraction": reliable / n if n else np.nan,
                "mean_internal_stop_count": float(group["internal_stop_count"].mean()) if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def build_codon_usage_summary(codon_metrics: pd.DataFrame) -> pd.DataFrame:
    reliable = codon_metrics[codon_metrics["codon_reliable"]].copy()
    rows = []
    for key, group in reliable.groupby(GROUP_COLS, dropna=False):
        panel, subtype, protein = key
        for codon in ALL_CODONS:
            rows.append(
                {
                    "panel": panel,
                    "subtype": subtype,
                    "protein": protein,
                    "codon": codon,
                    "amino_acid": CODON_TO_AA[codon],
                    "n_reliable_sequences": len(group),
                    "total_count": int(group[f"codon_count_{codon}"].sum()),
                    "mean_count": float(group[f"codon_count_{codon}"].mean()),
                    "mean_frequency": float(group[f"codon_freq_{codon}"].mean()),
                    "median_frequency": float(group[f"codon_freq_{codon}"].median()),
                }
            )
    return pd.DataFrame(rows)


def build_rscu_summary(codon_metrics: pd.DataFrame) -> pd.DataFrame:
    reliable = codon_metrics[codon_metrics["codon_reliable"]].copy()
    sense_codons = [codon for codon in ALL_CODONS if CODON_TO_AA[codon] != "*"]
    rows = []
    for key, group in reliable.groupby(GROUP_COLS, dropna=False):
        panel, subtype, protein = key
        for codon in sense_codons:
            values = group[f"rscu_{codon}"].dropna()
            rows.append(
                {
                    "panel": panel,
                    "subtype": subtype,
                    "protein": protein,
                    "codon": codon,
                    "amino_acid": CODON_TO_AA[codon],
                    "n_reliable_sequences": len(values),
                    "mean_rscu": float(values.mean()) if not values.empty else np.nan,
                    "median_rscu": float(values.median()) if not values.empty else np.nan,
                }
            )
    return pd.DataFrame(rows)


def build_kmer_entropy_summary(seq_rows_with_sequence: pd.DataFrame, seq_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, metrics_group in seq_metrics.groupby(GROUP_COLS, dropna=False):
        panel, subtype, protein = key
        sequence_group = seq_rows_with_sequence[
            (seq_rows_with_sequence["panel"] == panel)
            & (seq_rows_with_sequence["subtype"] == subtype)
            & (seq_rows_with_sequence["protein"] == protein)
        ]
        for k in [3, 4, 5]:
            counts: Counter[str] = Counter()
            for seq in sequence_group["sequence"]:
                counts.update(kmer_counts(seq, k))
            total = sum(counts.values())
            top = counts.most_common(3)
            row = {
                "panel": panel,
                "subtype": subtype,
                "protein": protein,
                "k": k,
                "n_sequences": len(metrics_group),
                "mean_entropy": float(metrics_group[f"k{k}_entropy"].mean()),
                "median_entropy": float(metrics_group[f"k{k}_entropy"].median()),
                "mean_valid_kmers": float(metrics_group[f"k{k}_valid_count"].mean()),
            }
            for rank in range(1, 4):
                kmer, count = top[rank - 1] if len(top) >= rank else ("", 0)
                row[f"top{rank}_kmer"] = kmer
                row[f"top{rank}_frequency"] = count / total if total else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def write_aggregate_tables(seq_rows_with_sequence: pd.DataFrame, seq_metrics: pd.DataFrame, codon_metrics: pd.DataFrame) -> dict[str, Path]:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    tables = {
        "gc_cpg_upa": (build_gc_cpg_upa_summary(seq_metrics), TABLES_DIR / "phase2_gc_cpg_upa_summary.csv"),
        "dinucleotide": (build_dinucleotide_summary(seq_metrics), TABLES_DIR / "phase2_dinucleotide_odds_summary.csv"),
        "codon_usage": (build_codon_usage_summary(codon_metrics), TABLES_DIR / "phase2_codon_usage_summary.csv"),
        "rscu": (build_rscu_summary(codon_metrics), TABLES_DIR / "phase2_rscu_summary.csv"),
        "kmer_entropy": (build_kmer_entropy_summary(seq_rows_with_sequence, seq_metrics), TABLES_DIR / "phase2_kmer_entropy_summary.csv"),
        "translation_qc": (build_translation_qc_summary(codon_metrics), TABLES_DIR / "phase2_translation_qc_summary.csv"),
    }
    paths = {}
    for name, (df, path) in tables.items():
        df.to_csv(path, index=False)
        paths[name] = path
    return paths


def _group_label(df: pd.DataFrame) -> pd.Series:
    return df["protein"] + "-" + df["subtype"]


def plot_gc_cpg_upa(seq_metrics: pd.DataFrame, outpath: Path) -> None:
    mvp = seq_metrics[seq_metrics["panel"] == "mvp"].copy()
    if mvp.empty:
        mvp = seq_metrics.copy()
    mvp["group"] = _group_label(mvp)
    groups = sorted(mvp["group"].unique())
    metrics = [("gc_content", "GC fraction"), ("cpg_oe", "CpG O/E"), ("upa_oe", "UpA O/E")]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    for ax, (metric, label) in zip(axes, metrics, strict=False):
        data = [mvp.loc[mvp["group"] == group, metric].dropna() for group in groups]
        ax.boxplot(data, tick_labels=groups, showfliers=False)
        ax.set_title(label)
        ax.tick_params(axis="x", rotation=25)
        ax.set_ylabel(label)
    fig.suptitle("Figure 4. Descriptive composition metrics by subtype and protein", y=1.02)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_dinucleotide_heatmap(dinuc_summary: pd.DataFrame, outpath: Path) -> None:
    mvp = dinuc_summary[dinuc_summary["panel"] == "mvp"].copy()
    if mvp.empty:
        mvp = dinuc_summary.copy()
    mvp["group"] = mvp["protein"] + "-" + mvp["subtype"]
    matrix = mvp.pivot_table(index="group", columns="dinucleotide", values="mean", aggfunc="mean").reindex(columns=list(DINUCLEOTIDES))
    fig, ax = plt.subplots(figsize=(12, 4.8))
    image = ax.imshow(matrix.values, aspect="auto", cmap="viridis")
    ax.set_title("Figure 5. Mean dinucleotide odds ratios")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    fig.colorbar(image, ax=ax, label="odds ratio")
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_rscu_or_qc(rscu_summary: pd.DataFrame, qc_summary: pd.DataFrame, outpath: Path) -> None:
    mvp_qc = qc_summary[qc_summary["panel"] == "mvp"].copy()
    acceptable = not mvp_qc.empty and bool((mvp_qc["reliable_fraction"] >= 0.5).all()) and bool((mvp_qc["n_codon_reliable"] >= 100).all())
    fig, ax = plt.subplots(figsize=(14, 5.5))
    if acceptable and not rscu_summary.empty:
        mvp = rscu_summary[rscu_summary["panel"] == "mvp"].copy()
        mvp["group"] = mvp["protein"] + "-" + mvp["subtype"]
        matrix = mvp.pivot_table(index="group", columns="codon", values="mean_rscu", aggfunc="mean")
        image = ax.imshow(matrix.values, aspect="auto", cmap="magma")
        ax.set_title("Figure 6. Mean RSCU for translation-QC reliable sequences")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index)
        fig.colorbar(image, ax=ax, label="RSCU")
    else:
        mvp_qc["group"] = mvp_qc["protein"] + "-" + mvp_qc["subtype"]
        ax.bar(mvp_qc["group"], mvp_qc["reliable_fraction"], color="#4C78A8")
        ax.set_ylim(0, 1)
        ax.set_ylabel("reliable fraction")
        ax.set_title("Figure 6. Codon/translation QC summary")
        ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_kmer_entropy(kmer_summary: pd.DataFrame, outpath: Path) -> None:
    mvp = kmer_summary[kmer_summary["panel"] == "mvp"].copy()
    if mvp.empty:
        mvp = kmer_summary.copy()
    mvp["group"] = mvp["protein"] + "-" + mvp["subtype"]
    pivot = mvp.pivot_table(index="group", columns="k", values="mean_entropy", aggfunc="mean")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    pivot.plot(kind="bar", ax=ax)
    ax.set_title("Figure 7. Mean k-mer entropy by group")
    ax.set_ylabel("Shannon entropy")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_figures(seq_metrics: pd.DataFrame, table_paths: dict[str, Path]) -> dict[str, Path]:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    dinuc = pd.read_csv(table_paths["dinucleotide"], keep_default_na=False)
    rscu_table = pd.read_csv(table_paths["rscu"], keep_default_na=False)
    qc = pd.read_csv(table_paths["translation_qc"], keep_default_na=False)
    kmer = pd.read_csv(table_paths["kmer_entropy"], keep_default_na=False)
    paths = {
        "gc_cpg_upa": FIGURES_DIR / "fig4_gc_cpg_upa_by_group.png",
        "dinucleotide": FIGURES_DIR / "fig5_dinucleotide_odds_heatmap.png",
        "rscu": FIGURES_DIR / "fig6_codon_usage_rscu_heatmap.png",
        "kmer": FIGURES_DIR / "fig7_kmer_entropy_by_group.png",
    }
    plot_gc_cpg_upa(seq_metrics, paths["gc_cpg_upa"])
    plot_dinucleotide_heatmap(dinuc, paths["dinucleotide"])
    plot_rscu_or_qc(rscu_table, qc, paths["rscu"])
    plot_kmer_entropy(kmer, paths["kmer"])
    return paths


def _markdown_table(df: pd.DataFrame, max_rows: int = 16) -> str:
    shown = df.head(max_rows).copy()
    if shown.empty:
        return "_No rows._"
    columns = list(shown.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for _, row in shown.iterrows():
        values = []
        for col in columns:
            value = row[col]
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Table truncated to {max_rows} of {len(df)} rows._")
    return "\n".join(lines)


def write_report(table_paths: dict[str, Path], figure_paths: dict[str, Path]) -> None:
    gc = pd.read_csv(table_paths["gc_cpg_upa"], keep_default_na=False)
    qc = pd.read_csv(table_paths["translation_qc"], keep_default_na=False)
    dinuc = pd.read_csv(table_paths["dinucleotide"], keep_default_na=False)
    codon = pd.read_csv(table_paths["codon_usage"], keep_default_na=False)
    rscu_table = pd.read_csv(table_paths["rscu"], keep_default_na=False)
    kmer = pd.read_csv(table_paths["kmer_entropy"], keep_default_na=False)

    mvp_gc = gc[(gc["panel"] == "mvp") & (gc["metric"].isin(["sequence_length", "gc_content", "cpg_oe", "upa_oe"]))]
    mvp_qc = qc[qc["panel"] == "mvp"]
    mvp_dinuc = dinuc[dinuc["panel"] == "mvp"].sort_values("mean", ascending=False)
    mvp_kmer = kmer[kmer["panel"] == "mvp"]
    reliable = int(mvp_qc["n_codon_reliable"].sum()) if not mvp_qc.empty else 0
    total = int(mvp_qc["n_sequences"].sum()) if not mvp_qc.empty else 0
    figure6_is_rscu = not mvp_qc.empty and bool((mvp_qc["reliable_fraction"] >= 0.5).all()) and bool((mvp_qc["n_codon_reliable"] >= 100).all())
    figure6_mode = (
        "RSCU heatmap, because each MVP group passed the codon-QC threshold."
        if figure6_is_rscu
        else "codon/translation QC summary, because at least one MVP group was below the codon-QC threshold."
    )
    codon_preview = codon[(codon["panel"] == "mvp") & (codon["amino_acid"] != "*")].sort_values("mean_frequency", ascending=False)
    rscu_preview = rscu_table[rscu_table["panel"] == "mvp"].sort_values("mean_rscu", ascending=False)

    report = f"""# FluGenome3D Phase 2 sequence-context report

This phase provides a descriptive sequence-context audit over the local Phase 1 panels. It compares compositional and codon-level patterns across HA/NA and H1N1/H3N2. It does not implement GROVER, BPE tokenization, 3D structure mapping, prediction, antigenicity, vaccine, escape, fitness, or sequence optimization.

## Panel coverage

- Smoke panel local metrics: `data/processed/metrics/smoke_sequence_metrics.parquet` and `data/processed/metrics/smoke_codon_metrics.parquet`.
- MVP panel local metrics: `data/processed/metrics/mvp_sequence_metrics.parquet` and `data/processed/metrics/mvp_codon_metrics.parquet`.
- MVP sequences analyzed: {total}.
- MVP codon/translation-QC reliable sequences: {reliable}.

## Composition summary

GC fraction, CpG observed/expected and UpA observed/expected are summarized as descriptive features. UpA is measured as TA in DNA alphabet.

{_markdown_table(mvp_gc)}

## Dinucleotide odds ratios

The table below shows the highest mean dinucleotide odds ratios in the MVP panel by group. Full 16-dinucleotide summaries are in `results/tables/phase2_dinucleotide_odds_summary.csv`.

{_markdown_table(mvp_dinuc[["subtype", "protein", "dinucleotide", "n", "mean", "median", "q05", "q95"]])}

## Translation and codon QC

Codon usage and RSCU are reported only for sequences that pass the naive frame check, ambiguity check, translation check and internal-stop check. These checks use the available nucleotide strings as provided; they do not infer CDS boundaries.

{_markdown_table(mvp_qc)}

## Codon usage and RSCU

Codon usage summaries use translation-QC reliable sequences only. Stop codons are included in codon usage counts but excluded from RSCU.

Top codons by mean frequency:

{_markdown_table(codon_preview[["subtype", "protein", "codon", "amino_acid", "n_reliable_sequences", "mean_frequency", "median_frequency"]])}

Top codons by mean RSCU:

{_markdown_table(rscu_preview[["subtype", "protein", "codon", "amino_acid", "n_reliable_sequences", "mean_rscu", "median_rscu"]])}

## K-mer entropy

K-mer profiles were summarized for k=3,4,5. Versionable outputs keep entropy and top-kmer aggregate summaries only, not row-level sequence strings.

{_markdown_table(mvp_kmer)}

## Figures

- `{figure_paths["gc_cpg_upa"]}`: GC, CpG O/E and UpA O/E by subtype/protein group.
- `{figure_paths["dinucleotide"]}`: 16-dinucleotide mean odds-ratio heatmap.
- `{figure_paths["rscu"]}`: {figure6_mode}
- `{figure_paths["kmer"]}`: k-mer entropy by group for k=3,4,5.

## What this phase can say

- This phase provides a descriptive sequence-context audit.
- We compare compositional and codon-level patterns across HA/NA and H1N1/H3N2.
- CpG/UpA and dinucleotide odds ratios are summarized as descriptive features.
- Codon usage is reported only after translation/frame QC.

## What this phase cannot say

- It does not predict antigenic drift.
- It does not identify escape mutations.
- It does not explain pathogenicity.
- It does not predict vaccine candidates.
- It does not claim codon usage explains fitness.
- It does not claim CpG/UpA determines antigenicity.
- It does not make causal claims.

## Limitations

- The translation QC is a naive check on the available nucleotide strings and does not infer curated CDS start/end coordinates.
- Ambiguous bases reduce valid k-mer and codon counts.
- The full panel remains a later extension; Phase 2 figures are anchored on the MVP panel.
- Aggregated summaries are public-safe, while row-level metrics remain local/gitignored.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")


def run_panel(panel: str) -> dict[str, Path]:
    panel_path = PANEL_PATHS[panel]
    if not panel_path.exists():
        raise FileNotFoundError(f"Missing panel: {panel_path}")
    panel_df = pd.read_parquet(panel_path)
    seq_rows = panel_to_sequence_rows(panel, panel_df)
    seq_metrics = compute_sequence_metrics(seq_rows)
    codon_metrics = compute_codon_metrics(seq_rows)
    return write_local_metrics(panel, seq_metrics, codon_metrics)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute FluGenome3D Phase 2 sequence-context metrics.")
    parser.add_argument("--panel", choices=sorted(PANEL_PATHS), required=True)
    args = parser.parse_args()

    paths = run_panel(args.panel)
    seq_metrics, codon_metrics = load_available_metrics()
    # Recreate sequence rows only for panels available locally so top-kmer summaries can remain aggregate-only.
    seq_rows = []
    for panel in sorted(seq_metrics["panel"].unique()):
        panel_df = pd.read_parquet(PANEL_PATHS[panel])
        seq_rows.append(panel_to_sequence_rows(panel, panel_df))
    seq_rows_df = pd.concat(seq_rows, ignore_index=True)
    table_paths = write_aggregate_tables(seq_rows_df, seq_metrics, codon_metrics)
    figure_paths = write_figures(seq_metrics, table_paths)
    write_report(table_paths, figure_paths)

    print(f"Wrote local metrics for {args.panel}: {paths['sequence']}, {paths['codon']}")
    print(f"Wrote aggregate Phase 2 tables: {', '.join(str(p) for p in table_paths.values())}")
    print(f"Wrote Phase 2 report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
