from __future__ import annotations

from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

from .codon_usage import ALL_CODONS, CODON_TO_AA
from .sequence_metrics import clean_sequence


def all_kmers(k: int) -> list[str]:
    return ["".join(chars) for chars in product("ACGT", repeat=k)]


def build_kmer_count_matrix(sequences: Iterable[object], k: int, mode: str = "overlapping") -> tuple[sparse.csr_matrix, list[str]]:
    if mode not in {"overlapping", "nonoverlapping"}:
        raise ValueError("mode must be 'overlapping' or 'nonoverlapping'")
    feature_names = all_kmers(k)
    index = {kmer: idx for idx, kmer in enumerate(feature_names)}
    step = 1 if mode == "overlapping" else k
    data: list[int] = []
    rows: list[int] = []
    cols: list[int] = []
    n_rows = 0
    for row_idx, seq in enumerate(sequences):
        n_rows += 1
        s = clean_sequence(seq)
        counts: dict[int, int] = {}
        for start in range(0, len(s) - k + 1, step):
            kmer = s[start : start + k]
            col = index.get(kmer)
            if col is not None:
                counts[col] = counts.get(col, 0) + 1
        for col, value in counts.items():
            rows.append(row_idx)
            cols.append(col)
            data.append(value)
    matrix = sparse.csr_matrix((data, (rows, cols)), shape=(n_rows, len(feature_names)), dtype=np.float32)
    return matrix, feature_names


def row_normalize(matrix: sparse.spmatrix) -> sparse.csr_matrix:
    csr = matrix.tocsr(copy=True).astype(np.float32)
    row_sums = np.asarray(csr.sum(axis=1)).ravel()
    nonzero = row_sums > 0
    inv = np.zeros_like(row_sums, dtype=np.float32)
    inv[nonzero] = 1.0 / row_sums[nonzero]
    return sparse.diags(inv).dot(csr).tocsr()


def build_kmer_frequency_matrix(sequences: Iterable[object], k: int, mode: str = "overlapping") -> tuple[sparse.csr_matrix, list[str]]:
    counts, names = build_kmer_count_matrix(sequences, k, mode=mode)
    return row_normalize(counts), names


def build_kmer_tfidf_matrix(sequences: Iterable[object], k: int, mode: str = "overlapping") -> tuple[sparse.csr_matrix, list[str]]:
    counts, names = build_kmer_count_matrix(sequences, k, mode=mode)
    tf = row_normalize(counts)
    n_docs = counts.shape[0]
    df = np.asarray((counts > 0).sum(axis=0)).ravel()
    idf = np.log((1 + n_docs) / (1 + df)) + 1.0
    tfidf = tf.dot(sparse.diags(idf.astype(np.float32))).tocsr()
    norms = np.sqrt(np.asarray(tfidf.multiply(tfidf).sum(axis=1)).ravel())
    nonzero = norms > 0
    inv = np.zeros_like(norms, dtype=np.float32)
    inv[nonzero] = 1.0 / norms[nonzero]
    return sparse.diags(inv).dot(tfidf).tocsr(), names


def build_codon_frequency_matrix(codon_metrics: pd.DataFrame) -> tuple[sparse.csr_matrix, list[str]]:
    feature_names = list(ALL_CODONS)
    cols = [f"codon_freq_{codon}" for codon in feature_names]
    missing = [col for col in cols if col not in codon_metrics.columns]
    if missing:
        raise ValueError(f"Missing codon frequency columns: {', '.join(missing[:5])}")
    values = codon_metrics[cols].fillna(0.0).to_numpy(dtype=np.float32)
    return sparse.csr_matrix(values), feature_names


def build_rscu_matrix(codon_metrics: pd.DataFrame) -> tuple[sparse.csr_matrix, list[str]]:
    feature_names = [codon for codon in ALL_CODONS if CODON_TO_AA[codon] != "*"]
    cols = [f"rscu_{codon}" for codon in feature_names]
    missing = [col for col in cols if col not in codon_metrics.columns]
    if missing:
        raise ValueError(f"Missing RSCU columns: {', '.join(missing[:5])}")
    values = codon_metrics[cols].fillna(0.0).to_numpy(dtype=np.float32)
    return sparse.csr_matrix(values), feature_names


def aggregate_group_centroids(feature_matrix: sparse.spmatrix, metadata: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    matrix = feature_matrix.tocsr()
    groups = metadata[group_cols].astype(str).agg("|".join, axis=1)
    rows = []
    for group in sorted(groups.unique()):
        idx = np.flatnonzero(groups.to_numpy() == group)
        if idx.size == 0:
            continue
        centroid = np.asarray(matrix[idx].mean(axis=0)).ravel()
        row = {col: value for col, value in zip(group_cols, group.split("|"), strict=False)}
        row["group"] = group
        row["n_sequences"] = int(idx.size)
        row["centroid"] = centroid
        rows.append(row)
    return pd.DataFrame(rows)


def compute_feature_entropy(feature_matrix: sparse.spmatrix) -> np.ndarray:
    matrix = row_normalize(feature_matrix).tocsr()
    entropies = np.zeros(matrix.shape[0], dtype=np.float32)
    for i in range(matrix.shape[0]):
        row = matrix.getrow(i)
        vals = row.data[row.data > 0]
        entropies[i] = float(-np.sum(vals * np.log2(vals))) if vals.size else np.nan
    return entropies


def top_features_by_group(
    feature_matrix: sparse.spmatrix,
    metadata: pd.DataFrame,
    feature_names: list[str],
    group_col: str | list[str],
    top_n: int = 10,
) -> pd.DataFrame:
    group_cols = [group_col] if isinstance(group_col, str) else list(group_col)
    matrix = feature_matrix.tocsr()
    groups = metadata[group_cols].astype(str).agg("|".join, axis=1)
    rows = []
    for group in sorted(groups.unique()):
        idx = np.flatnonzero(groups.to_numpy() == group)
        if idx.size == 0:
            continue
        means = np.asarray(matrix[idx].mean(axis=0)).ravel()
        top_idx = np.argsort(means)[::-1][:top_n]
        group_parts = group.split("|")
        base = {col: value for col, value in zip(group_cols, group_parts, strict=False)}
        for rank, feat_idx in enumerate(top_idx, start=1):
            rows.append({**base, "group": group, "rank": rank, "feature": feature_names[feat_idx], "mean_value": float(means[feat_idx])})
    return pd.DataFrame(rows)


def save_feature_names(feature_names: list[str], path: str | Path) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("\n".join(feature_names) + "\n", encoding="utf-8")


def load_feature_names(path: str | Path) -> list[str]:
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
