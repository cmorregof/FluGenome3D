#!/usr/bin/env python
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from flugenome3d.cds_refinement import rescue_cds
from flugenome3d.codon_usage import ALL_CODONS, CODON_TO_AA, codon_counts, codon_frequencies_from_counts, rscu_from_counts, translation_qc
from flugenome3d.utils import load_yaml


PANEL_PATH = Path("data/processed/panels/mvp_panel.parquet")
PHASE2_CODON_PATH = Path("data/processed/metrics/mvp_codon_metrics.parquet")
PANEL_OUTDIR = Path("data/processed/panels")
METRICS_OUTDIR = Path("data/processed/metrics")
TABLES_DIR = Path("results/tables")
FIGURES_DIR = Path("results/figures")
REPORT_PATH = Path("reports/phase3_cds_refinement_report.md")
GROUP_COLS = ["subtype", "protein"]


def panel_to_sequence_rows(panel_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in panel_df.iterrows():
        for protein, seq_col, hash_col in [("HA", "ha_sequence", "ha_sha256"), ("NA", "na_sequence", "na_sha256")]:
            rows.append(
                {
                    "internal_strain_id": row["internal_strain_id"],
                    "internal_sequence_id": f"{row['internal_strain_id']}_{protein}",
                    "pair_sha256": row["pair_sha256"],
                    "sequence_sha256": row[hash_col],
                    "subtype": row["subtype"],
                    "protein": protein,
                    "year": row["year"],
                    "raw_sequence": row[seq_col],
                }
            )
    return pd.DataFrame(rows)


def load_phase3_inputs() -> pd.DataFrame:
    panel = pd.read_parquet(PANEL_PATH)
    seq_rows = panel_to_sequence_rows(panel)
    codon = pd.read_parquet(PHASE2_CODON_PATH)
    qc_cols = [
        "internal_sequence_id",
        "frame_fail",
        "ambiguous_fail",
        "internal_stop_fail",
        "translation_fail",
        "internal_stop_count",
        "aa_length",
        "codon_total",
        "codon_reliable",
    ]
    return seq_rows.merge(codon[qc_cols], on="internal_sequence_id", how="left", validate="one_to_one")


def failure_breakdown(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for key, group in df.groupby(GROUP_COLS, dropna=False):
        subtype, protein = key
        n = len(group)
        rows.append(
            {
                "subtype": subtype,
                "protein": protein,
                "n_sequences": n,
                "n_naive_qc_pass": int(group["codon_reliable"].sum()),
                "n_any_qc_fail": int((~group["codon_reliable"]).sum()),
                "n_frame_fail": int(group["frame_fail"].sum()),
                "n_internal_stop_fail": int(group["internal_stop_fail"].sum()),
                "n_ambiguous_fail": int(group["ambiguous_fail"].sum()),
                "n_translation_fail": int(group["translation_fail"].sum()),
                "mean_length_nt": float(group["raw_sequence"].str.len().mean()),
                "median_length_nt": float(group["raw_sequence"].str.len().median()),
                "mean_internal_stop_count": float(group["internal_stop_count"].mean()),
            }
        )
    return pd.DataFrame(rows)


def length_mod3_by_group(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work["length_nt"] = work["raw_sequence"].str.len()
    work["length_mod3"] = work["length_nt"] % 3
    return work.groupby([*GROUP_COLS, "length_mod3"], dropna=False).size().reset_index(name="n_sequences")


def internal_stop_distribution(df: pd.DataFrame) -> pd.DataFrame:
    work = df.copy()
    work["internal_stop_count"] = work["internal_stop_count"].fillna(-1).astype(int)
    return work.groupby([*GROUP_COLS, "internal_stop_count"], dropna=False).size().reset_index(name="n_sequences")


def refine_sequences(df: pd.DataFrame, config: dict) -> pd.DataFrame:
    rows = []
    for row in df.itertuples(index=False):
        result = rescue_cds(row.raw_sequence, protein=row.protein, subtype=row.subtype, config=config)
        out = {
            "internal_strain_id": row.internal_strain_id,
            "internal_sequence_id": row.internal_sequence_id,
            "pair_sha256": row.pair_sha256,
            "sequence_sha256": row.sequence_sha256,
            "subtype": row.subtype,
            "protein": row.protein,
            "year": row.year,
            "raw_sequence": row.raw_sequence,
            "refined_sequence": result.refined_sequence,
            "naive_codon_reliable": bool(row.codon_reliable),
            "naive_frame_fail": bool(row.frame_fail),
            "naive_internal_stop_fail": bool(row.internal_stop_fail),
            "naive_ambiguous_fail": bool(row.ambiguous_fail),
            "naive_translation_fail": bool(row.translation_fail),
            "naive_internal_stop_count": row.internal_stop_count,
        }
        out.update(result.public_dict())
        rows.append(out)
    return pd.DataFrame(rows)


def write_cds_panels(refined: pd.DataFrame) -> dict[str, Path]:
    PANEL_OUTDIR.mkdir(parents=True, exist_ok=True)
    strict = refined[refined["status"] == "strict_pass"].copy()
    strict_plus_rescued = refined[refined["status"].isin(["strict_pass", "rescued"])].copy()
    final = strict_plus_rescued.copy()
    paths = {
        "strict": PANEL_OUTDIR / "mvp_cds_strict_panel.parquet",
        "rescued": PANEL_OUTDIR / "mvp_cds_rescued_panel.parquet",
        "refined": PANEL_OUTDIR / "mvp_cds_refined_panel.parquet",
    }
    strict.to_parquet(paths["strict"], index=False)
    strict_plus_rescued.to_parquet(paths["rescued"], index=False)
    final.to_parquet(paths["refined"], index=False)
    return paths


def compute_refined_codon_metrics(refined_panel: pd.DataFrame, cds_config: dict) -> pd.DataFrame:
    rows = []
    max_ambig = float(cds_config.get("max_ambiguous_fraction", 0.01))
    accepted = refined_panel[refined_panel["status"].isin(["strict_pass", "rescued"])].copy()
    for row in accepted.itertuples(index=False):
        seq = row.refined_sequence
        qc = translation_qc(seq, max_ambiguous_fraction=max_ambig)
        counts = codon_counts(seq) if qc["codon_reliable"] else Counter()
        freqs = codon_frequencies_from_counts(counts) if qc["codon_reliable"] else {codon: np.nan for codon in ALL_CODONS}
        rscu_values = rscu_from_counts(counts) if qc["codon_reliable"] else {codon: np.nan for codon in ALL_CODONS if CODON_TO_AA[codon] != "*"}
        out = {
            "internal_sequence_id": row.internal_sequence_id,
            "pair_sha256": row.pair_sha256,
            "sequence_sha256": row.sequence_sha256,
            "subtype": row.subtype,
            "protein": row.protein,
            "year": row.year,
            "cds_status": row.status,
            "rescue_method": row.rescue_method,
            "chosen_frame": row.chosen_frame,
            "trim_left": row.trim_left,
            "trim_right": row.trim_right,
        }
        out.update(qc)
        for codon in ALL_CODONS:
            out[f"codon_count_{codon}"] = counts[codon]
            out[f"codon_freq_{codon}"] = freqs[codon]
        for codon, value in rscu_values.items():
            out[f"rscu_{codon}"] = value
        rows.append(out)
    metrics = pd.DataFrame(rows)
    METRICS_OUTDIR.mkdir(parents=True, exist_ok=True)
    metrics.to_parquet(METRICS_OUTDIR / "mvp_cds_refined_codon_metrics.parquet", index=False)
    return metrics


def refined_translation_qc_summary(metrics: pd.DataFrame, refined: pd.DataFrame) -> pd.DataFrame:
    status_counts = refined.groupby([*GROUP_COLS, "status"], dropna=False).size().reset_index(name="n_sequences")
    rows = []
    for key, group in metrics.groupby(GROUP_COLS, dropna=False):
        subtype, protein = key
        n = len(group)
        rows.append(
            {
                "subtype": subtype,
                "protein": protein,
                "n_refined_sequences": n,
                "n_codon_reliable_after_refinement": int(group["codon_reliable"].sum()),
                "reliable_fraction_after_refinement": float(group["codon_reliable"].mean()) if n else np.nan,
                "n_frame_fail_after_refinement": int(group["frame_fail"].sum()),
                "n_internal_stop_fail_after_refinement": int(group["internal_stop_fail"].sum()),
                "n_translation_fail_after_refinement": int(group["translation_fail"].sum()),
                "n_strict_pass": int(status_counts[(status_counts["subtype"] == subtype) & (status_counts["protein"] == protein) & (status_counts["status"] == "strict_pass")]["n_sequences"].sum()),
                "n_rescued": int(status_counts[(status_counts["subtype"] == subtype) & (status_counts["protein"] == protein) & (status_counts["status"] == "rescued")]["n_sequences"].sum()),
                "n_unrescued": int(status_counts[(status_counts["subtype"] == subtype) & (status_counts["protein"] == protein) & (status_counts["status"] == "unrescued")]["n_sequences"].sum()),
            }
        )
    return pd.DataFrame(rows)


def refined_codon_usage_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    reliable = metrics[metrics["codon_reliable"]].copy()
    rows = []
    for key, group in reliable.groupby(GROUP_COLS, dropna=False):
        subtype, protein = key
        for codon in ALL_CODONS:
            rows.append(
                {
                    "subtype": subtype,
                    "protein": protein,
                    "codon": codon,
                    "amino_acid": CODON_TO_AA[codon],
                    "n_reliable_sequences": len(group),
                    "total_count": int(group[f"codon_count_{codon}"].sum()),
                    "mean_frequency": float(group[f"codon_freq_{codon}"].mean()),
                    "median_frequency": float(group[f"codon_freq_{codon}"].median()),
                }
            )
    return pd.DataFrame(rows)


def refined_rscu_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    reliable = metrics[metrics["codon_reliable"]].copy()
    sense_codons = [codon for codon in ALL_CODONS if CODON_TO_AA[codon] != "*"]
    rows = []
    for key, group in reliable.groupby(GROUP_COLS, dropna=False):
        subtype, protein = key
        for codon in sense_codons:
            values = group[f"rscu_{codon}"].dropna()
            rows.append(
                {
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


def rescue_status_summary(refined: pd.DataFrame) -> pd.DataFrame:
    return refined.groupby([*GROUP_COLS, "status"], dropna=False).size().reset_index(name="n_sequences")


def write_tables(diagnosis: pd.DataFrame, refined: pd.DataFrame, metrics: pd.DataFrame) -> dict[str, Path]:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    tables = {
        "failure_breakdown": (failure_breakdown(diagnosis), TABLES_DIR / "phase3_cds_qc_failure_breakdown.csv"),
        "length_mod3": (length_mod3_by_group(diagnosis), TABLES_DIR / "phase3_length_mod3_by_group.csv"),
        "internal_stop_distribution": (internal_stop_distribution(diagnosis), TABLES_DIR / "phase3_internal_stop_distribution.csv"),
        "rescue_status": (rescue_status_summary(refined), TABLES_DIR / "phase3_rescue_status_by_group.csv"),
        "refined_codon_usage": (refined_codon_usage_summary(metrics), TABLES_DIR / "phase3_refined_codon_usage_summary.csv"),
        "refined_rscu": (refined_rscu_summary(metrics), TABLES_DIR / "phase3_refined_rscu_summary.csv"),
        "refined_translation_qc": (refined_translation_qc_summary(metrics, refined), TABLES_DIR / "phase3_refined_translation_qc_summary.csv"),
    }
    paths = {}
    for name, (df, path) in tables.items():
        df.to_csv(path, index=False)
        paths[name] = path
    return paths


def plot_failure_breakdown(path: Path, outpath: Path) -> None:
    df = pd.read_csv(path, keep_default_na=False)
    labels = df["protein"] + "-" + df["subtype"]
    fig, ax = plt.subplots(figsize=(10, 5))
    bottom = np.zeros(len(df))
    for col, color in [
        ("n_frame_fail", "#4C78A8"),
        ("n_internal_stop_fail", "#F58518"),
        ("n_ambiguous_fail", "#B279A2"),
        ("n_translation_fail", "#E45756"),
    ]:
        ax.bar(labels, df[col], bottom=bottom, label=col.replace("n_", ""), color=color)
        bottom += df[col].to_numpy()
    ax.set_title("Figure 8. CDS QC failure flags by group")
    ax.set_ylabel("flag count; flags can overlap")
    ax.tick_params(axis="x", rotation=25)
    ax.legend()
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_length_mod3(path: Path, outpath: Path) -> None:
    df = pd.read_csv(path, keep_default_na=False)
    df["group"] = df["protein"] + "-" + df["subtype"]
    pivot = df.pivot_table(index="group", columns="length_mod3", values="n_sequences", aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 5))
    pivot.plot(kind="bar", stacked=True, ax=ax)
    ax.set_title("Figure 9. Nucleotide length modulo 3 by group")
    ax.set_ylabel("sequences")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_refined_rscu_or_qc(rscu_path: Path, qc_path: Path, outpath: Path) -> None:
    rscu_df = pd.read_csv(rscu_path, keep_default_na=False)
    qc = pd.read_csv(qc_path, keep_default_na=False)
    enough = not qc.empty and bool((qc["n_codon_reliable_after_refinement"] >= 100).all())
    fig, ax = plt.subplots(figsize=(14, 5.5))
    if enough and not rscu_df.empty:
        rscu_df["group"] = rscu_df["protein"] + "-" + rscu_df["subtype"]
        matrix = rscu_df.pivot_table(index="group", columns="codon", values="mean_rscu", aggfunc="mean")
        image = ax.imshow(matrix.values, aspect="auto", cmap="magma")
        ax.set_title("Figure 10. Refined CDS mean RSCU")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index)
        fig.colorbar(image, ax=ax, label="RSCU")
    else:
        qc["group"] = qc["protein"] + "-" + qc["subtype"]
        ax.bar(qc["group"], qc["n_codon_reliable_after_refinement"], color="#4C78A8")
        ax.set_title("Figure 10. Refined CDS QC count")
        ax.set_ylabel("reliable sequences")
        ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_rescue_status(path: Path, outpath: Path) -> None:
    df = pd.read_csv(path, keep_default_na=False)
    df["group"] = df["protein"] + "-" + df["subtype"]
    pivot = df.pivot_table(index="group", columns="status", values="n_sequences", aggfunc="sum", fill_value=0)
    fig, ax = plt.subplots(figsize=(9, 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, color=["#54A24B", "#F58518", "#E45756"])
    ax.set_title("Figure 11. CDS rescue status by group")
    ax.set_ylabel("sequences")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_figures(table_paths: dict[str, Path]) -> dict[str, Path]:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    paths = {
        "failure_breakdown": FIGURES_DIR / "fig8_cds_qc_failure_breakdown.png",
        "length_mod3": FIGURES_DIR / "fig9_length_mod3_by_group.png",
        "refined_rscu": FIGURES_DIR / "fig10_refined_rscu_heatmap.png",
        "rescue_status": FIGURES_DIR / "fig11_rescue_status_by_group.png",
    }
    plot_failure_breakdown(table_paths["failure_breakdown"], paths["failure_breakdown"])
    plot_length_mod3(table_paths["length_mod3"], paths["length_mod3"])
    plot_refined_rscu_or_qc(table_paths["refined_rscu"], table_paths["refined_translation_qc"], paths["refined_rscu"])
    plot_rescue_status(table_paths["rescue_status"], paths["rescue_status"])
    return paths


def _markdown_table(df: pd.DataFrame, max_rows: int = 16) -> str:
    shown = df.head(max_rows).copy()
    if shown.empty:
        return "_No rows._"
    cols = list(shown.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in shown.iterrows():
        values = []
        for col in cols:
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
    failure = pd.read_csv(table_paths["failure_breakdown"], keep_default_na=False)
    status = pd.read_csv(table_paths["rescue_status"], keep_default_na=False)
    refined_qc = pd.read_csv(table_paths["refined_translation_qc"], keep_default_na=False)
    length_mod = pd.read_csv(table_paths["length_mod3"], keep_default_na=False)

    strict_total = int(status.loc[status["status"] == "strict_pass", "n_sequences"].sum())
    rescued_total = int(status.loc[status["status"] == "rescued", "n_sequences"].sum())
    unrescued_total = int(status.loc[status["status"] == "unrescued", "n_sequences"].sum())
    refined_total = strict_total + rescued_total
    codon_reliable_after = int(refined_qc["n_codon_reliable_after_refinement"].sum()) if not refined_qc.empty else 0

    report = f"""# FluGenome3D Phase 3 CDS refinement report

This phase was necessary because Phase 2 found that naive codon-frame assumptions were not reliable for all HA/NA nucleotide records. Non-coding sequence-context metrics remain valid on raw nucleotide context, while codon usage and RSCU require a separate CDS-aware subset.

## Representations

- `raw_nucleotide_context`: original normalized nucleotide strings for GC, CpG/UpA, dinucleotide odds and k-mer entropy. No coding frame is required.
- `cds_strict`: sequences that pass conservative frame, ambiguity, translation, internal-stop and protein-length checks without trimming.
- `cds_rescued`: sequences that fail naive QC but pass transparent 0-2 nt trim/frame rescue rules.

## Why naive QC failed

Likely contributors include UTR or non-coding flanks, partial records, gaps or ambiguous characters, nonzero CDS frame offsets, aligned-vs-CDS nucleotide representations, non-multiple-of-3 lengths, and records not directly translatable from position 0. This phase audits those possibilities without assuming a single cause.

{_markdown_table(failure)}

## Length modulo 3

{_markdown_table(length_mod)}

## Rescue rules

- Normalize to uppercase DNA alphabet and convert U to T.
- Remove gaps only when configured, while recording `gaps_removed` and `gap_count_removed`.
- Try deterministic trims of 0-2 nt at the left and 0-2 nt at the right.
- Accept only candidates with no internal stops, ambiguity under threshold, no ambiguous amino acids, and expected HA/NA protein-length range.
- Leave all other sequences as `unrescued`.
- No reverse-complement rescue, alignment rescue, reference-guided repair, prediction, optimization or biological validation is performed.

## Rescue outcome

- Strict pass sequences: {strict_total}
- Rescued sequences: {rescued_total}
- Unrescued sequences: {unrescued_total}
- Final refined CDS panel size: {refined_total}
- Codon-reliable sequences after refinement: {codon_reliable_after}

{_markdown_table(status)}

## Refined codon/RSCU QC

Codon usage and RSCU are now reported on the strict/rescued CDS subset with explicit rescue status and QC flags. These are more reliable than Phase 2 naive codon summaries, but they are still not biological validation of CDS boundaries.

{_markdown_table(refined_qc)}

## Outputs

Local gitignored panels:

- `data/processed/panels/mvp_cds_strict_panel.parquet`
- `data/processed/panels/mvp_cds_rescued_panel.parquet`
- `data/processed/panels/mvp_cds_refined_panel.parquet`
- `data/processed/metrics/mvp_cds_refined_codon_metrics.parquet`

Public aggregate tables:

- `results/tables/phase3_cds_qc_failure_breakdown.csv`
- `results/tables/phase3_length_mod3_by_group.csv`
- `results/tables/phase3_internal_stop_distribution.csv`
- `results/tables/phase3_rescue_status_by_group.csv`
- `results/tables/phase3_refined_codon_usage_summary.csv`
- `results/tables/phase3_refined_rscu_summary.csv`
- `results/tables/phase3_refined_translation_qc_summary.csv`

Figures:

- `{figure_paths["failure_breakdown"]}`
- `{figure_paths["length_mod3"]}`
- `{figure_paths["refined_rscu"]}`
- `{figure_paths["rescue_status"]}`

## Permitted claims

- We identified limitations of naive codon-frame assumptions in HA/NA sequence records.
- We separated non-coding sequence-context metrics from CDS-dependent codon analyses.
- Codon usage and RSCU are reported only on strict/rescued CDS subsets with explicit QC.
- This phase improves reliability of downstream tokenization and codon-level analysis.

## Prohibited claims

- Rescued CDS are biologically validated.
- Codon usage explains antigenic drift.
- RSCU predicts fitness.
- Frame rescue identifies functional variants.
- Translation QC implies antigenic relevance.

## Recommendation for Phase 4

Use `raw_nucleotide_context` for tokenizer-independent composition and k-mer analyses, and use `mvp_cds_refined_panel` only for CDS-dependent codon/RSCU summaries. The next phase should focus on tokenizer-ready representations over the MVP/refined panels without GROVER/BPE unless explicitly approved.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run FluGenome3D Phase 3 CDS-aware refinement on the MVP panel.")
    parser.add_argument("--filters", default="config/filters.yml")
    args = parser.parse_args()

    filters = load_yaml(args.filters)
    cds_config = filters.get("cds_refinement", {})
    inputs = load_phase3_inputs()
    refined = refine_sequences(inputs, cds_config)
    panel_paths = write_cds_panels(refined)
    refined_panel = pd.read_parquet(panel_paths["refined"])
    metrics = compute_refined_codon_metrics(refined_panel, cds_config)
    table_paths = write_tables(inputs, refined, metrics)
    figure_paths = write_figures(table_paths)
    write_report(table_paths, figure_paths)

    status = pd.read_csv(table_paths["rescue_status"], keep_default_na=False)
    print(status.to_string(index=False))
    print(f"Wrote refined CDS panel: {panel_paths['refined']} ({len(refined_panel)} sequence rows)")
    print(f"Wrote refined codon metrics: {METRICS_OUTDIR / 'mvp_cds_refined_codon_metrics.parquet'}")
    print(f"Wrote report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
