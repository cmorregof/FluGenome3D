from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from .reduction import compute_pairwise_group_distances
from .representations import aggregate_group_centroids, compute_feature_entropy, top_features_by_group


def matrix_sparsity(feature_matrix: sparse.spmatrix) -> float:
    total = feature_matrix.shape[0] * feature_matrix.shape[1]
    if total == 0:
        return float("nan")
    nnz = feature_matrix.nnz if sparse.issparse(feature_matrix) else int(np.count_nonzero(feature_matrix))
    return 1.0 - (nnz / total)


def representation_summary_row(name: str, source: str, feature_matrix, pca_info: dict[str, object], umap_info: dict[str, object] | None = None) -> dict[str, object]:
    explained = pca_info.get("explained_variance_ratio", [])
    ev1 = explained[0] if len(explained) > 0 else np.nan
    ev2 = explained[1] if len(explained) > 1 else np.nan
    nnz = feature_matrix.nnz if sparse.issparse(feature_matrix) else int(np.count_nonzero(feature_matrix))
    return {
        "representation": name,
        "source": source,
        "n_sequences": int(feature_matrix.shape[0]),
        "n_features": int(feature_matrix.shape[1]),
        "nnz": int(nnz),
        "sparsity": matrix_sparsity(feature_matrix),
        "pca_explained_variance_pc1": ev1,
        "pca_explained_variance_pc2": ev2,
        "pca_explained_variance_total_2pc": float(np.nansum(explained[:2])),
        "umap_method": (umap_info or {}).get("method", ""),
    }


def centroid_distances_for_representation(name: str, feature_matrix, metadata: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    centroids = aggregate_group_centroids(feature_matrix, metadata, group_cols)
    distances = compute_pairwise_group_distances(centroids)
    distances.insert(0, "representation", name)
    return distances


def top_features_for_representation(name: str, feature_matrix, metadata: pd.DataFrame, feature_names: list[str], group_cols: list[str], top_n: int = 10) -> pd.DataFrame:
    top = top_features_by_group(feature_matrix, metadata, feature_names, group_cols, top_n=top_n)
    top.insert(0, "representation", name)
    return top


def entropy_summary_by_group(name: str, feature_matrix, metadata: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    entropy = compute_feature_entropy(feature_matrix)
    work = metadata[group_cols].copy()
    work["entropy"] = entropy
    summary = work.groupby(group_cols, dropna=False)["entropy"].agg(n="count", mean="mean", median="median").reset_index()
    summary.insert(0, "representation", name)
    return summary
