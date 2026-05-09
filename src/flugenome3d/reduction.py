from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import svds


def _to_dense_float(feature_matrix) -> np.ndarray:
    if sparse.issparse(feature_matrix):
        return feature_matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(feature_matrix, dtype=np.float32)


def run_pca(feature_matrix, n_components: int = 2) -> tuple[np.ndarray, dict[str, object]]:
    x = _to_dense_float(feature_matrix)
    if x.ndim != 2:
        raise ValueError("feature_matrix must be 2-dimensional")
    n_samples, n_features = x.shape
    if n_samples == 0:
        return np.empty((0, n_components), dtype=np.float32), {"explained_variance_ratio": [np.nan] * n_components}
    x = x - x.mean(axis=0, keepdims=True)
    max_components = max(1, min(n_components, n_samples - 1, n_features - 1 if n_features > 1 else 1))
    if max_components == 1:
        u, s, _ = np.linalg.svd(x, full_matrices=False)
        s = s[:1]
        embedding = u[:, :1] * s
    else:
        u, s, _ = svds(x, k=max_components)
        order = np.argsort(s)[::-1]
        s = s[order]
        u = u[:, order]
        embedding = u * s
    total_var = float(np.sum(x * x) / max(n_samples - 1, 1))
    explained = (s**2) / max(n_samples - 1, 1)
    ratio = (explained / total_var).tolist() if total_var else [np.nan] * len(s)
    if embedding.shape[1] < n_components:
        pad = np.zeros((n_samples, n_components - embedding.shape[1]), dtype=np.float32)
        embedding = np.hstack([embedding, pad])
        ratio.extend([np.nan] * (n_components - len(ratio)))
    return embedding[:, :n_components].astype(np.float32), {"explained_variance_ratio": ratio[:n_components], "method": "pca"}


def run_umap(feature_matrix, n_neighbors: int = 30, min_dist: float = 0.1) -> tuple[np.ndarray, dict[str, object]]:
    try:
        import umap  # type: ignore

        reducer = umap.UMAP(n_components=2, n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
        embedding = reducer.fit_transform(feature_matrix)
        return np.asarray(embedding, dtype=np.float32), {"method": "umap", "n_neighbors": n_neighbors, "min_dist": min_dist}
    except Exception as exc:
        embedding, info = run_pca(feature_matrix, n_components=2)
        info = {**info, "method": "pca_fallback_for_umap", "fallback_reason": type(exc).__name__}
        return embedding, info


def _pairwise_distances_dense(x: np.ndarray) -> np.ndarray:
    sq = np.sum(x * x, axis=1, keepdims=True)
    d2 = np.maximum(sq + sq.T - 2 * (x @ x.T), 0)
    return np.sqrt(d2).astype(np.float32)


def compute_silhouette_safe(embedding_or_features, labels, max_samples: int = 2000, random_state: int = 42) -> float:
    labels = np.asarray(labels)
    valid = labels.astype(str)
    unique, counts = np.unique(valid, return_counts=True)
    if unique.size < 2 or np.any(counts < 2):
        return float("nan")
    x = _to_dense_float(embedding_or_features)
    if x.shape[0] != labels.shape[0]:
        raise ValueError("labels length must match matrix rows")
    if x.shape[0] > max_samples:
        rng = np.random.default_rng(random_state)
        selected: list[int] = []
        per_group = max(2, max_samples // unique.size)
        for group in unique:
            idx = np.flatnonzero(valid == group)
            take = min(per_group, idx.size)
            selected.extend(rng.choice(idx, size=take, replace=False).tolist())
        if len(selected) > max_samples:
            selected = rng.choice(np.asarray(selected), size=max_samples, replace=False).tolist()
        selected = sorted(selected)
        x = x[selected]
        valid = valid[selected]
        unique, counts = np.unique(valid, return_counts=True)
        if unique.size < 2 or np.any(counts < 2):
            return float("nan")
    distances = _pairwise_distances_dense(x)
    scores = []
    for i, label in enumerate(valid):
        same = valid == label
        same[i] = False
        a = float(distances[i, same].mean()) if np.any(same) else 0.0
        b_values = [float(distances[i, valid == other].mean()) for other in unique if other != label and np.any(valid == other)]
        if not b_values:
            continue
        b = min(b_values)
        denom = max(a, b)
        scores.append((b - a) / denom if denom else 0.0)
    return float(np.mean(scores)) if scores else float("nan")


def compute_pairwise_group_distances(group_centroids: pd.DataFrame) -> pd.DataFrame:
    if group_centroids.empty:
        return pd.DataFrame(columns=["group_a", "group_b", "distance"])
    labels = group_centroids["group"].astype(str).tolist()
    matrix = np.vstack(group_centroids["centroid"].to_list()).astype(np.float32)
    dist = squareform(pdist(matrix, metric="euclidean"))
    rows = []
    for i, group_a in enumerate(labels):
        for j, group_b in enumerate(labels):
            rows.append({"group_a": group_a, "group_b": group_b, "distance": float(dist[i, j])})
    return pd.DataFrame(rows)
