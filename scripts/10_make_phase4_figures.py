#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from flugenome3d.feature_analysis import centroid_distances_for_representation, representation_summary_row, top_features_for_representation
from flugenome3d.reduction import compute_silhouette_safe, run_pca, run_umap
from flugenome3d.representations import load_feature_names


REP_DIR = Path("data/processed/representations")
TABLE_DIR = Path("results/tables")
FIGURE_DIR = Path("results/figures")
REPORT_PATH = Path("reports/phase4_representation_audit_report.md")
RAW_META = REP_DIR / "mvp_raw_representation_metadata.parquet"
CDS_META = REP_DIR / "mvp_cds_representation_metadata.parquet"

REPRESENTATIONS = {
    "kmer3_freq": {"path": REP_DIR / "mvp_kmer3_freq.npz", "features": REP_DIR / "feature_names_mvp_kmer3_freq.txt", "metadata": RAW_META, "source": "raw_nucleotide"},
    "kmer4_freq": {"path": REP_DIR / "mvp_kmer4_freq.npz", "features": REP_DIR / "feature_names_mvp_kmer4_freq.txt", "metadata": RAW_META, "source": "raw_nucleotide"},
    "kmer5_freq": {"path": REP_DIR / "mvp_kmer5_freq.npz", "features": REP_DIR / "feature_names_mvp_kmer5_freq.txt", "metadata": RAW_META, "source": "raw_nucleotide"},
    "kmer3_tfidf": {"path": REP_DIR / "mvp_kmer3_tfidf.npz", "features": REP_DIR / "feature_names_mvp_kmer3_tfidf.txt", "metadata": RAW_META, "source": "raw_nucleotide"},
    "kmer4_tfidf": {"path": REP_DIR / "mvp_kmer4_tfidf.npz", "features": REP_DIR / "feature_names_mvp_kmer4_tfidf.txt", "metadata": RAW_META, "source": "raw_nucleotide"},
    "codon_freq": {"path": REP_DIR / "mvp_codon_freq.npz", "features": REP_DIR / "feature_names_mvp_codon_freq.txt", "metadata": CDS_META, "source": "cds_refined"},
    "rscu": {"path": REP_DIR / "mvp_rscu.npz", "features": REP_DIR / "feature_names_mvp_rscu.txt", "metadata": CDS_META, "source": "cds_refined"},
}


def load_representation(name: str):
    spec = REPRESENTATIONS[name]
    matrix = sparse.load_npz(spec["path"])
    metadata = pd.read_parquet(spec["metadata"])
    feature_names = load_feature_names(spec["features"])
    return matrix, metadata, feature_names, spec


def stratified_sample(metadata: pd.DataFrame, max_total: int = 4000, random_state: int = 42) -> np.ndarray:
    if len(metadata) <= max_total:
        return np.arange(len(metadata))
    rng = np.random.default_rng(random_state)
    selected = []
    groups = sorted(metadata["protein_subtype"].unique())
    per_group = max(1, max_total // len(groups))
    for group in groups:
        idx = np.flatnonzero(metadata["protein_subtype"].to_numpy() == group)
        take = min(per_group, idx.size)
        selected.extend(rng.choice(idx, size=take, replace=False).tolist())
    return np.asarray(sorted(selected), dtype=int)


def write_embedding_table(name: str, embedding: np.ndarray, metadata: pd.DataFrame, kind: str) -> None:
    out = metadata[["internal_sequence_id", "subtype", "protein", "protein_subtype", "year"]].copy()
    out["representation"] = name
    out["axis1"] = embedding[:, 0]
    out["axis2"] = embedding[:, 1]
    out.to_parquet(REP_DIR / f"{name}_{kind}_embedding.parquet", index=False)


def compute_tables() -> tuple[dict[str, np.ndarray], dict[str, pd.DataFrame], dict[str, dict[str, object]]]:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    embeddings: dict[str, np.ndarray] = {}
    metadata_by_rep: dict[str, pd.DataFrame] = {}
    pca_info_by_rep: dict[str, dict[str, object]] = {}
    summary_rows = []
    silhouette_rows = []
    distance_frames = []
    top_kmer_frames = []
    top_codon_frames = []
    top_rscu_frames = []
    umap_info = {}

    for name in REPRESENTATIONS:
        matrix, metadata, feature_names, spec = load_representation(name)
        embedding, pca_info = run_pca(matrix, n_components=2)
        embeddings[name] = embedding
        metadata_by_rep[name] = metadata
        pca_info_by_rep[name] = pca_info
        write_embedding_table(name, embedding, metadata, "pca")

        if name == "kmer4_tfidf":
            sample_idx = stratified_sample(metadata)
            umap_embedding, info = run_umap(matrix[sample_idx], n_neighbors=30, min_dist=0.1)
            umap_meta = metadata.iloc[sample_idx].reset_index(drop=True)
            write_embedding_table(name, umap_embedding, umap_meta, "umap")
            umap_info = {**info, "sample_n": len(sample_idx)}

        summary_rows.append(representation_summary_row(name, spec["source"], matrix, pca_info, umap_info if name == "kmer4_tfidf" else None))

        for label_type in ["protein", "subtype", "protein_subtype"]:
            silhouette_rows.append(
                {
                    "representation": name,
                    "label_type": label_type,
                    "silhouette": compute_silhouette_safe(embedding, metadata[label_type].to_numpy()),
                    "space": "pca_2d",
                }
            )

        distance_frames.append(centroid_distances_for_representation(name, matrix, metadata, ["protein_subtype"]))

        if name.startswith("kmer") and name.endswith("freq"):
            top = top_features_for_representation(name, matrix, metadata, feature_names, ["protein", "subtype"], top_n=10)
            top["k"] = top["feature"].str.len()
            top_kmer_frames.append(top)
        elif name == "codon_freq":
            top_codon_frames.append(top_features_for_representation(name, matrix, metadata, feature_names, ["protein", "subtype"], top_n=10))
        elif name == "rscu":
            top_rscu_frames.append(top_features_for_representation(name, matrix, metadata, feature_names, ["protein", "subtype"], top_n=10))

    pd.DataFrame(summary_rows).to_csv(TABLE_DIR / "phase4_representation_summary.csv", index=False)
    pd.DataFrame(silhouette_rows).to_csv(TABLE_DIR / "phase4_silhouette_scores.csv", index=False)
    pd.concat(distance_frames, ignore_index=True).to_csv(TABLE_DIR / "phase4_group_centroid_distances.csv", index=False)
    pd.concat(top_kmer_frames, ignore_index=True).to_csv(TABLE_DIR / "phase4_top_kmers_by_group.csv", index=False)
    pd.concat(top_codon_frames, ignore_index=True).to_csv(TABLE_DIR / "phase4_top_codons_by_group.csv", index=False)
    pd.concat(top_rscu_frames, ignore_index=True).to_csv(TABLE_DIR / "phase4_top_rscu_by_group.csv", index=False)
    return embeddings, metadata_by_rep, pca_info_by_rep


def _scatter_embedding(embedding: np.ndarray, metadata: pd.DataFrame, title: str, outpath: Path) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    groups = sorted(metadata["protein_subtype"].unique())
    cmap = plt.get_cmap("tab10")
    for idx, group in enumerate(groups):
        mask = metadata["protein_subtype"].to_numpy() == group
        ax.scatter(embedding[mask, 0], embedding[mask, 1], s=6, alpha=0.55, label=group, color=cmap(idx))
    ax.set_title(title)
    ax.set_xlabel("axis 1")
    ax.set_ylabel("axis 2")
    ax.legend(markerscale=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_distance_heatmaps(outpath: Path) -> None:
    distances = pd.read_csv(TABLE_DIR / "phase4_group_centroid_distances.csv", keep_default_na=False)
    reps = ["kmer4_tfidf", "codon_freq", "rscu"]
    fig, axes = plt.subplots(1, len(reps), figsize=(14, 4.2))
    for ax, rep in zip(axes, reps, strict=False):
        part = distances[distances["representation"] == rep]
        matrix = part.pivot_table(index="group_a", columns="group_b", values="distance", aggfunc="mean")
        image = ax.imshow(matrix.values, cmap="viridis")
        ax.set_title(rep)
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels(matrix.index, fontsize=8)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Figure 16. Group centroid distances")
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def plot_silhouette(outpath: Path) -> None:
    scores = pd.read_csv(TABLE_DIR / "phase4_silhouette_scores.csv")
    keep = ["kmer3_freq", "kmer4_freq", "kmer5_freq", "codon_freq", "rscu"]
    scores = scores[scores["representation"].isin(keep)].copy()
    fig, ax = plt.subplots(figsize=(10, 5))
    labels = sorted(scores["label_type"].unique())
    x = np.arange(len(keep))
    width = 0.25
    for idx, label in enumerate(labels):
        vals = [scores[(scores["representation"] == rep) & (scores["label_type"] == label)]["silhouette"].mean() for rep in keep]
        ax.bar(x + (idx - 1) * width, vals, width=width, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(keep, rotation=25, ha="right")
    ax.set_ylabel("silhouette on PCA 2D")
    ax.set_title("Figure 17. Representation silhouette comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_figures(embeddings: dict[str, np.ndarray], metadata_by_rep: dict[str, pd.DataFrame]) -> None:
    _scatter_embedding(embeddings["kmer4_tfidf"], metadata_by_rep["kmer4_tfidf"], "Figure 12. k-mer TF-IDF PCA by group", FIGURE_DIR / "fig12_kmer_pca_by_group.png")

    umap_path = REP_DIR / "kmer4_tfidf_umap_embedding.parquet"
    if umap_path.exists():
        umap_df = pd.read_parquet(umap_path)
        embedding = umap_df[["axis1", "axis2"]].to_numpy()
        _scatter_embedding(embedding, umap_df, "Figure 13. k-mer TF-IDF UMAP/fallback by group", FIGURE_DIR / "fig13_kmer_umap_by_group.png")

    _scatter_embedding(embeddings["codon_freq"], metadata_by_rep["codon_freq"], "Figure 14. Codon-frequency PCA by group", FIGURE_DIR / "fig14_codon_pca_by_group.png")
    _scatter_embedding(embeddings["rscu"], metadata_by_rep["rscu"], "Figure 15. RSCU PCA by group", FIGURE_DIR / "fig15_rscu_pca_by_group.png")
    plot_distance_heatmaps(FIGURE_DIR / "fig16_group_centroid_distance_heatmap.png")
    plot_silhouette(FIGURE_DIR / "fig17_representation_silhouette_comparison.png")


def _markdown_table(df: pd.DataFrame, max_rows: int = 14) -> str:
    shown = df.head(max_rows)
    if shown.empty:
        return "_No rows._"
    cols = list(shown.columns)
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in shown.iterrows():
        vals = []
        for col in cols:
            val = row[col]
            vals.append(f"{val:.6g}" if isinstance(val, float) else str(val))
        lines.append("| " + " | ".join(vals) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Table truncated to {max_rows} of {len(df)} rows._")
    return "\n".join(lines)


def write_report() -> None:
    summary = pd.read_csv(TABLE_DIR / "phase4_representation_summary.csv", keep_default_na=False)
    silhouettes = pd.read_csv(TABLE_DIR / "phase4_silhouette_scores.csv", keep_default_na=False)
    distances = pd.read_csv(TABLE_DIR / "phase4_group_centroid_distances.csv", keep_default_na=False)
    top_kmers = pd.read_csv(TABLE_DIR / "phase4_top_kmers_by_group.csv", keep_default_na=False)
    best_group = silhouettes[silhouettes["label_type"] == "protein_subtype"].sort_values("silhouette", ascending=False).head(5)

    report = f"""# FluGenome3D Phase 4 representation audit report

This phase compares simple nucleotide and codon-level representations of Influenza A HA/NA before any GROVER or BPE tokenization. The analysis is descriptive and representation-focused.

## Dataset separation

- Raw nucleotide representations use `mvp_panel`, because k-mer and compositional features do not require a validated CDS frame.
- CDS/codon representations use `mvp_cds_refined_panel` and `mvp_cds_refined_codon_metrics`, because codon-frequency and RSCU features require explicit CDS QC.

## Representations built

{_markdown_table(summary)}

## PCA and UMAP summary

PCA was run for all representations. UMAP was attempted for `kmer4_tfidf` on a stratified sample; if UMAP is unavailable, the script records and plots a deterministic PCA fallback under the UMAP figure filename.

## Silhouette scores

Silhouette scores are computed on 2D PCA embeddings for protein, subtype and protein-subtype labels. These scores are descriptive separation diagnostics, not predictive performance metrics.

{_markdown_table(best_group)}

## Centroid distances

Pairwise centroid distances are stored in `results/tables/phase4_group_centroid_distances.csv`. Figure 16 visualizes group distances for k-mer TF-IDF, codon frequency and RSCU representations.

{_markdown_table(distances.head(12))}

## Top features

Top k-mer, codon-frequency and RSCU features are aggregated by group. These are high-weight descriptive features, not mutation, fitness, antigenicity or escape markers.

{_markdown_table(top_kmers.head(12))}

## Figures

- `results/figures/fig12_kmer_pca_by_group.png`
- `results/figures/fig13_kmer_umap_by_group.png`
- `results/figures/fig14_codon_pca_by_group.png`
- `results/figures/fig15_rscu_pca_by_group.png`
- `results/figures/fig16_group_centroid_distance_heatmap.png`
- `results/figures/fig17_representation_silhouette_comparison.png`

## Limitations

- UMAP is optional and may fall back to PCA if the dependency is not installed.
- CDS-dependent representations cover only the refined CDS subset.
- Feature matrices are local/gitignored and should not be treated as public sequence redistribution artifacts.
- This phase does not make claims about antigenic drift, vaccine relevance, escape, fitness, causality or prediction.

## Recommendation for Phase 5

Use these simple representations as baselines for tokenizer audits. The next phase can compare codon/k-mer tokenization behavior against these matrices before considering GROVER or BPE.
"""
    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    embeddings, metadata_by_rep, _ = compute_tables()
    write_figures(embeddings, metadata_by_rep)
    write_report()
    print(f"Wrote Phase 4 tables to {TABLE_DIR}")
    print(f"Wrote Phase 4 figures to {FIGURE_DIR}")
    print(f"Wrote Phase 4 report: {REPORT_PATH}")


if __name__ == "__main__":
    main()
