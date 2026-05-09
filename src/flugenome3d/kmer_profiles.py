from __future__ import annotations

from collections import Counter
from math import log2
from typing import Iterable

import numpy as np
import pandas as pd

from .sequence_metrics import clean_sequence


def kmer_counts(seq: object, k: int) -> Counter[str]:
    s = clean_sequence(seq)
    if k <= 0 or len(s) < k:
        return Counter()
    return Counter(s[i : i + k] for i in range(len(s) - k + 1) if set(s[i : i + k]) <= set("ACGT"))


def kmer_frequencies(seq: object, k: int) -> dict[str, float]:
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    if total == 0:
        return {}
    return {kmer: count / total for kmer, count in counts.items()}


def kmer_entropy_from_counts(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total == 0:
        return np.nan
    return -sum((count / total) * log2(count / total) for count in counts.values() if count > 0)


def kmer_entropy(seq: object, k: int) -> float:
    return kmer_entropy_from_counts(kmer_counts(seq, k))


def sequence_kmer_metrics(seq: object, k_values: Iterable[int] = (3, 4, 5)) -> dict[str, object]:
    row: dict[str, object] = {}
    for k in k_values:
        counts = kmer_counts(seq, k)
        total = sum(counts.values())
        row[f"k{k}_valid_count"] = total
        row[f"k{k}_unique"] = len(counts)
        row[f"k{k}_entropy"] = kmer_entropy_from_counts(counts)
    return row


def top_kmers_by_group(df: pd.DataFrame, group_cols: list[str], k: int, top_n: int = 5) -> pd.DataFrame:
    rows = []
    for group_key, group in df.groupby(group_cols, dropna=False):
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        counts: Counter[str] = Counter()
        for seq in group["sequence"]:
            counts.update(kmer_counts(seq, k))
        total = sum(counts.values())
        for rank, (kmer, count) in enumerate(counts.most_common(top_n), start=1):
            row = dict(zip(group_cols, group_key, strict=False))
            row.update({"k": k, "rank": rank, "kmer": kmer, "count": count, "frequency": count / total if total else np.nan})
            rows.append(row)
    return pd.DataFrame(rows)
