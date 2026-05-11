#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

from flugenome3d.antigenlm_bridge import (
    build_latent_pca_points,
    build_latent_tsne_points,
    load_json,
    safe_point_records,
    write_public_tables,
)


PROJECT = Path(__file__).resolve().parents[1]


def markdown_table(df: pd.DataFrame, max_rows: int | None = None) -> str:
    if df.empty:
        return ""
    use = df.head(max_rows).copy() if max_rows else df.copy()
    columns = list(use.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in use.to_dict(orient="records"):
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def read_config(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def write_report(summary: dict, tables: dict[str, pd.DataFrame]) -> None:
    out = PROJECT / "reports" / "phase7_antigenlm_bridge_report.md"
    cache = summary["cache"]
    spearman = tables["phase7_antigenlm_spearman_summary.csv"]
    pca = tables["phase7_antigenlm_pca_summary.csv"]
    clade = tables["phase7_antigenlm_clade_enrichment_summary.csv"]
    lines = [
        "# Phase 7 AntigenLM bridge report",
        "",
        "Phase 7 imports the parent-repository AntigenLM geometry as a derived representation layer. It does not export raw sequences, accessions, isolate names, sequence hashes, FASTA, checkpoints, or restricted Parquet records.",
        "",
        "## Cache summary",
        "",
        f"- Records represented: {cache.get('n_records'):,}",
        f"- Embedding dimension: {cache.get('embedding_dim')}",
        f"- Year range: {cache.get('year_min')}-{cache.get('year_max')}",
        f"- Exported latent atlas points: {summary.get('n_exported_points'):,}",
        f"- Projection: {summary.get('projection')}",
        "",
        "## Why this matters",
        "",
        "The child repo already compares nucleotide, codon and token-level baselines. Phase 7 adds the learned AntigenLM layer from the thesis repo so the app can show how simple biological features relate to a learned influenza representation.",
        "",
        "## Molecular geometry summary",
        "",
        markdown_table(spearman) if not spearman.empty else "No Spearman table available.",
        "",
        "## PCA effective dimension",
        "",
        markdown_table(pca) if not pca.empty else "No PCA table available.",
        "",
        "## Clade enrichment summary",
        "",
        markdown_table(clade, max_rows=12) if not clade.empty else "No clade enrichment table available.",
        "",
        "## App boundary",
        "",
        "- The app may show reduced coordinates with hash-based IDs and minimal metadata.",
        "- The app may show aggregate Spearman, PCA, TwoNN, temporal-locality and clade-enrichment summaries.",
        "- The app must not show raw sequences, source accessions, isolate names, sequence hashes or checkpoints.",
        "",
        "## Claims allowed",
        "",
        "- AntigenLM embeddings are summarized as a learned representation layer.",
        "- Latent distances can be compared descriptively with molecular and temporal proxies.",
        "- The learned layer is compared against deterministic k-mer/codon/RSCU baselines.",
        "",
        "## Claims not allowed",
        "",
        "- This phase does not predict antigenicity, vaccine relevance, escape, pathogenicity, fitness or evolution.",
        "- This phase does not validate sequence generation or optimize sequences.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_figures(points: pd.DataFrame, tables: dict[str, pd.DataFrame], summary: dict) -> None:
    fig_dir = PROJECT / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    colors = {"H1N1": "#79d99c", "H3N2": "#5cdce2"}
    plt.figure(figsize=(8, 6), facecolor="#080807")
    ax = plt.gca()
    ax.set_facecolor("#01050a")
    for subtype, group in points.groupby("subtype"):
        ax.scatter(group["axis1"], group["axis2"], s=3, alpha=0.45, c=colors.get(subtype, "#edf7f4"), label=subtype)
    ax.set_title("AntigenLM PCA atlas sample", color="#edf7f4")
    ax.set_xlabel("PC1", color="#edf7f4")
    ax.set_ylabel("PC2", color="#edf7f4")
    ax.tick_params(colors="#9fb3ae")
    ax.legend(facecolor="#071017", edgecolor="#1b5561", labelcolor="#edf7f4")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig29_antigenlm_pca_by_subtype.png", dpi=220)
    plt.close()

    spearman = tables["phase7_antigenlm_spearman_summary.csv"]
    if not spearman.empty:
        pivot = spearman.pivot(index="metric", columns="subtype", values="rho_mean")
        ax = pivot.plot(kind="bar", figsize=(8, 5), color=["#79d99c", "#5cdce2"])
        ax.figure.set_facecolor("#080807")
        ax.set_facecolor("#01050a")
        ax.set_title("AntigenLM latent distance correlations", color="#edf7f4")
        ax.set_ylabel("Mean Spearman rho", color="#edf7f4")
        ax.tick_params(colors="#9fb3ae")
        ax.legend(facecolor="#071017", edgecolor="#1b5561", labelcolor="#edf7f4")
        plt.xticks(rotation=20, ha="right")
        plt.tight_layout()
        plt.savefig(fig_dir / "fig30_antigenlm_geometry_summary.png", dpi=220)
        plt.close()

    explained = summary.get("pca_explained_variance", [])
    if explained:
        plt.figure(figsize=(6, 4), facecolor="#080807")
        ax = plt.gca()
        ax.set_facecolor("#01050a")
        ax.bar([f"PC{i+1}" for i in range(len(explained))], explained, color="#5cdce2")
        ax.set_title("Export PCA variance", color="#edf7f4")
        ax.set_ylabel("Explained variance ratio", color="#edf7f4")
        ax.tick_params(colors="#9fb3ae")
        plt.tight_layout()
        plt.savefig(fig_dir / "fig31_antigenlm_export_pca_variance.png", dpi=220)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Phase 7 AntigenLM bridge artifacts.")
    parser.add_argument("--config", default=str(PROJECT / "config" / "phase7_9.yml"))
    args = parser.parse_args()
    cfg = read_config(args.config).get("phase7_antigenlm_bridge", {})

    cache_path = PROJECT / cfg.get("cache_path", "../results/embeddings_cache_full_all_available.pkl")
    max_points = int(cfg.get("max_export_points", 30000))
    tsne_max_points = int(cfg.get("tsne_max_export_points", 10000))
    random_seed = int(cfg.get("random_seed", 42))

    points, summary = build_latent_pca_points(cache_path, max_points=max_points, random_state=random_seed)
    tsne_2d_points, tsne_2d_summary = build_latent_tsne_points(
        cache_path,
        max_points=tsne_max_points,
        random_state=random_seed,
        n_components=2,
    )
    tsne_3d_points, tsne_3d_summary = build_latent_tsne_points(
        cache_path,
        max_points=tsne_max_points,
        random_state=random_seed,
        n_components=3,
    )
    local_dir = PROJECT / "data" / "processed" / "antigenlm"
    local_dir.mkdir(parents=True, exist_ok=True)
    points.to_parquet(local_dir / "antigenlm_full_pca_points.parquet", index=False)
    tsne_2d_points.to_parquet(local_dir / "antigenlm_tsne_2d_points.parquet", index=False)
    tsne_3d_points.to_parquet(local_dir / "antigenlm_tsne_3d_points.parquet", index=False)

    results_dir = PROJECT / "results" / "tables"
    tables = write_public_tables(
        PROJECT / cfg.get("latent_metrics_path", "../results/latent_geometry_full_metrics.json"),
        PROJECT / cfg.get("clade_enrichment_path", "../results/gisaid_clade_enrichment_results.json"),
        PROJECT / cfg.get("random_baseline_path", "../paper_revision_outputs/random_embedding_baseline_results.json"),
        results_dir,
    )
    pd.DataFrame([summary["cache"]]).to_csv(results_dir / "phase7_antigenlm_cache_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "projection": summary["projection"],
                "n_source_points": summary["n_source_points"],
                "n_exported_points": summary["n_exported_points"],
                "pca_explained_variance": ";".join(str(x) for x in summary["pca_explained_variance"]),
            }
        ]
    ).to_csv(results_dir / "phase7_antigenlm_projection_summary.csv", index=False)

    # Keep a compact local JSON manifest to make safe export deterministic without reloading the large pickle.
    manifest = {
        **summary,
        "point_schema": ["id", "x", "y", "z", "subtype", "year_bin", "representation", "source"],
        "points": safe_point_records(points, representation_label="AntigenLM PCA"),
        "additional_projections": [
            {
                "id": "antigenlm_tsne_2d",
                "label": "AntigenLM t-SNE 2D",
                "description": "Nonlinear t-SNE map of sampled AntigenLM HA+NA embeddings. Coordinates are derived artifacts with hash-based IDs.",
                "projection": "tsne_2d",
                "axis_labels": ["t-SNE 1", "t-SNE 2"],
                "point_schema": ["id", "x", "y", "z", "subtype", "year_bin", "representation", "source"],
                "n_source_points": tsne_2d_summary["n_source_points"],
                "n_exported_points": tsne_2d_summary["n_exported_points"],
                "sampling": tsne_2d_summary["sampling"],
                "tsne_parameters": tsne_2d_summary["tsne_parameters"],
                "privacy": "hash-based point IDs and coarse metadata only; no sequences, accessions, isolate names, sequence hashes, or checkpoint weights",
                "points": safe_point_records(tsne_2d_points, representation_label="AntigenLM t-SNE 2D"),
            },
            {
                "id": "antigenlm_tsne_3d",
                "label": "AntigenLM t-SNE 3D",
                "description": "Three-dimensional t-SNE map of sampled AntigenLM HA+NA embeddings. Coordinates are derived artifacts with hash-based IDs.",
                "projection": "tsne_3d",
                "axis_labels": ["t-SNE 1", "t-SNE 2", "t-SNE 3"],
                "point_schema": ["id", "x", "y", "z", "subtype", "year_bin", "representation", "source"],
                "n_source_points": tsne_3d_summary["n_source_points"],
                "n_exported_points": tsne_3d_summary["n_exported_points"],
                "sampling": tsne_3d_summary["sampling"],
                "tsne_parameters": tsne_3d_summary["tsne_parameters"],
                "privacy": "hash-based point IDs and coarse metadata only; no sequences, accessions, isolate names, sequence hashes, or checkpoint weights",
                "points": safe_point_records(tsne_3d_points, representation_label="AntigenLM t-SNE 3D"),
            },
        ],
    }
    (local_dir / "antigenlm_latent_atlas.local.json").write_text(
        __import__("json").dumps(manifest, separators=(",", ":"), ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    write_report(summary, tables)
    make_figures(points, tables, summary)
    print(f"Phase 7 complete: {len(points):,} safe latent points prepared locally.")


if __name__ == "__main__":
    main()
