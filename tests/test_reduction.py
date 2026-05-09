from __future__ import annotations

import numpy as np
import pandas as pd

from flugenome3d.reduction import compute_pairwise_group_distances, compute_silhouette_safe, run_pca


def test_run_pca_produces_requested_dimensions() -> None:
    x = np.array([[0, 1, 0], [0, 2, 0], [3, 0, 1], [4, 0, 1]], dtype=float)
    embedding, info = run_pca(x, n_components=2)
    assert embedding.shape == (4, 2)
    assert len(info["explained_variance_ratio"]) == 2


def test_silhouette_safe_handles_degenerate_labels() -> None:
    x = np.array([[0, 0], [1, 1], [2, 2]], dtype=float)
    assert np.isnan(compute_silhouette_safe(x, ["A", "A", "A"]))


def test_pairwise_group_distances_runs() -> None:
    centroids = pd.DataFrame({"group": ["A", "B"], "centroid": [np.array([0.0, 0.0]), np.array([3.0, 4.0])]})
    distances = compute_pairwise_group_distances(centroids)
    assert set(distances.columns) == {"group_a", "group_b", "distance"}
    assert distances[(distances["group_a"] == "A") & (distances["group_b"] == "B")]["distance"].iloc[0] == 5.0
