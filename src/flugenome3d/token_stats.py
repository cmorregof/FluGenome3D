from __future__ import annotations

from collections import Counter
from math import log2

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


def token_counts(tokens: list[str]) -> Counter[str]:
    return Counter(tokens)


def token_frequencies(tokens: list[str]) -> dict[str, float]:
    counts = token_counts(tokens)
    total = sum(counts.values())
    if not total:
        return {}
    return {token: count / total for token, count in counts.items()}


def token_entropy(tokens: list[str]) -> float:
    counts = token_counts(tokens)
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum((count / total) * log2(count / total) for count in counts.values() if count > 0)


def effective_vocab_size(tokens: list[str]) -> float:
    return float(2 ** token_entropy(tokens))


def tokens_per_kb(tokens: list[str], seq_length: int) -> float:
    return len(tokens) / (seq_length / 1000) if seq_length else np.nan


def group_token_distribution(token_table: pd.DataFrame, metadata: pd.DataFrame | None = None) -> pd.DataFrame:
    del metadata
    required = {"tokenizer", "protein", "subtype", "token", "count"}
    missing = required - set(token_table.columns)
    if missing:
        raise ValueError(f"Missing token_table columns: {', '.join(sorted(missing))}")
    grouped = token_table.groupby(["tokenizer", "protein", "subtype", "token"], dropna=False)["count"].sum().reset_index()
    totals = grouped.groupby(["tokenizer", "protein", "subtype"], dropna=False)["count"].transform("sum")
    grouped["frequency"] = grouped["count"] / totals
    grouped["protein_subtype"] = grouped["protein"] + "-" + grouped["subtype"]
    return grouped


def top_tokens_by_group(distribution: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    rows = []
    for key, group in distribution.groupby(["tokenizer", "protein", "subtype"], dropna=False):
        tokenizer, protein, subtype = key
        top = group.sort_values("frequency", ascending=False).head(top_n)
        for rank, row in enumerate(top.itertuples(index=False), start=1):
            rows.append(
                {
                    "tokenizer": tokenizer,
                    "protein": protein,
                    "subtype": subtype,
                    "protein_subtype": f"{protein}-{subtype}",
                    "rank": rank,
                    "token": row.token,
                    "frequency": float(row.frequency),
                    "count": int(row.count),
                }
            )
    return pd.DataFrame(rows)


def token_enrichment_by_group(
    distribution: pd.DataFrame,
    top_n: int = 10,
    pseudocount: float = 1e-9,
    min_count: int = 1,
    min_group_frequency: float = 0.0,
) -> pd.DataFrame:
    rows = []
    for tokenizer, tok_df in distribution.groupby("tokenizer", dropna=False):
        total_by_token = tok_df.groupby("token")["count"].sum()
        total = float(total_by_token.sum())
        for key, group in tok_df.groupby(["protein", "subtype"], dropna=False):
            protein, subtype = key
            group_total = float(group["count"].sum())
            merged = group.set_index("token")
            values = []
            for token, group_count in merged["count"].items():
                background_count = total_by_token[token] - group_count
                background_total = total - group_total
                group_freq = group_count / group_total if group_total else 0.0
                bg_freq = background_count / background_total if background_total else 0.0
                if group_count < min_count or group_freq < min_group_frequency:
                    continue
                values.append((token, group_freq, bg_freq, (group_freq + pseudocount) / (bg_freq + pseudocount), int(group_count)))
            for rank, (token, group_freq, bg_freq, enrichment, count) in enumerate(sorted(values, key=lambda x: x[3], reverse=True)[:top_n], start=1):
                rows.append(
                    {
                        "tokenizer": tokenizer,
                        "protein": protein,
                        "subtype": subtype,
                        "protein_subtype": f"{protein}-{subtype}",
                        "rank": rank,
                        "token": token,
                        "group_frequency": group_freq,
                        "background_frequency": bg_freq,
                        "enrichment_ratio": enrichment,
                        "count": count,
                    }
                )
    return pd.DataFrame(rows)


def jensen_shannon_distance_between_groups(distribution: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tokenizer, tok_df in distribution.groupby("tokenizer", dropna=False):
        pivot = tok_df.pivot_table(index="protein_subtype", columns="token", values="frequency", aggfunc="sum", fill_value=0.0)
        groups = list(pivot.index)
        values = pivot.to_numpy(dtype=float)
        for i, group_a in enumerate(groups):
            for j, group_b in enumerate(groups):
                distance = float(jensenshannon(values[i], values[j], base=2.0))
                rows.append({"tokenizer": tokenizer, "group_a": group_a, "group_b": group_b, "js_distance": distance})
    return pd.DataFrame(rows)
