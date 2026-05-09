from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy.spatial.distance import jensenshannon


GROUP_ORDER = ["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"]


def stratified_bootstrap_indices(
    metadata: pd.DataFrame,
    group_col: str,
    n_per_group: int | None = None,
    random_state: int | np.random.Generator | None = None,
) -> np.ndarray:
    rng = random_state if isinstance(random_state, np.random.Generator) else np.random.default_rng(random_state)
    indices: list[np.ndarray] = []
    for _, group in metadata.groupby(group_col, dropna=False):
        values = group.index.to_numpy()
        if values.size == 0:
            continue
        n = values.size if n_per_group is None else min(int(n_per_group), values.size)
        indices.append(rng.choice(values, size=n, replace=True))
    if not indices:
        return np.asarray([], dtype=int)
    return np.concatenate(indices).astype(int)


def summarize_bootstrap_ci(values: Iterable[float], ci: float = 0.95) -> dict[str, float]:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {"n": 0, "mean": 0.0, "std": 0.0, "ci_lower": 0.0, "ci_upper": 0.0, "ci_width": 0.0}
    alpha = (1.0 - ci) / 2.0
    lower = float(np.quantile(arr, alpha))
    upper = float(np.quantile(arr, 1.0 - alpha))
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "std": float(arr.std(ddof=1)) if arr.size > 1 else 0.0,
        "ci_lower": lower,
        "ci_upper": upper,
        "ci_width": upper - lower,
    }


def bootstrap_token_metrics(
    token_stats: pd.DataFrame,
    metadata: pd.DataFrame | None,
    tokenizer: str,
    n_bootstraps: int = 100,
    group_col: str = "protein_subtype",
    n_per_group: int | None = None,
    random_state: int | None = None,
) -> pd.DataFrame:
    del metadata
    stats = token_stats[token_stats["tokenizer"] == tokenizer].reset_index(drop=True)
    rng = np.random.default_rng(random_state)
    rows = []
    metric_cols = [
        "token_entropy_bits",
        "effective_vocab_size",
        "cpg_token_fraction",
        "upa_token_fraction",
    ]
    for bootstrap_id in range(n_bootstraps):
        sample_idx = stratified_bootstrap_indices(stats, group_col, n_per_group=n_per_group, random_state=rng)
        sampled = stats.loc[sample_idx]
        for group_value, group in sampled.groupby(group_col, dropna=False):
            row = {
                "bootstrap_id": bootstrap_id,
                "tokenizer": tokenizer,
                group_col: group_value,
                "n_sequences": int(group["internal_sequence_id"].nunique()) if "internal_sequence_id" in group else int(len(group)),
                "sample_rows": int(len(group)),
            }
            for col in metric_cols:
                row[f"mean_{col}"] = float(group[col].mean())
            rows.append(row)
    return pd.DataFrame(rows)


def _bootstrap_group_counts(
    tok_df: pd.DataFrame,
    rng: np.random.Generator,
    max_tokens_per_group: int | None = 200_000,
) -> tuple[dict[str, np.ndarray], list[str]]:
    pivot = tok_df.pivot_table(index="protein_subtype", columns="token", values="count", aggfunc="sum", fill_value=0)
    pivot = pivot.reindex([g for g in GROUP_ORDER if g in pivot.index])
    tokens = list(pivot.columns)
    sampled: dict[str, np.ndarray] = {}
    for group_name, row in pivot.iterrows():
        counts = row.to_numpy(dtype=float)
        total = counts.sum()
        if total <= 0:
            sampled[group_name] = np.zeros_like(counts, dtype=float)
            continue
        probs = counts / total
        n = int(total if max_tokens_per_group is None else min(total, max_tokens_per_group))
        sampled_counts = rng.multinomial(n, probs)
        sampled[group_name] = sampled_counts.astype(float)
    return sampled, tokens


def _counts_to_frequencies(sampled_counts: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    freqs = {}
    for group, counts in sampled_counts.items():
        denom = counts.sum()
        freqs[group] = counts / denom if denom else counts.astype(float)
    return freqs


def bootstrap_js_distances(
    group_token_distributions: pd.DataFrame,
    n_bootstraps: int = 100,
    random_state: int | None = None,
    max_tokens_per_group: int | None = 200_000,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_state)
    rows = []
    for tokenizer, tok_df in group_token_distributions.groupby("tokenizer", dropna=False):
        for bootstrap_id in range(n_bootstraps):
            sampled_counts, _ = _bootstrap_group_counts(tok_df, rng, max_tokens_per_group=max_tokens_per_group)
            sampled = _counts_to_frequencies(sampled_counts)
            groups = list(sampled)
            for i, group_a in enumerate(groups):
                for j, group_b in enumerate(groups):
                    distance = float(jensenshannon(sampled[group_a], sampled[group_b], base=2.0))
                    rows.append(
                        {
                            "bootstrap_id": bootstrap_id,
                            "tokenizer": tokenizer,
                            "group_a": group_a,
                            "group_b": group_b,
                            "comparison": f"{group_a} vs {group_b}",
                            "js_distance": distance,
                        }
                    )
    return pd.DataFrame(rows)


def _top_enriched_from_sampled_counts(
    sampled_counts: dict[str, np.ndarray],
    tokens: list[str],
    tokenizer: str,
    bootstrap_id: int,
    top_n: int = 20,
    min_count: int = 100,
    min_group_frequency: float = 1e-5,
    pseudocount: float = 1e-9,
) -> list[dict[str, object]]:
    count_matrix = np.vstack([sampled_counts[group] for group in sampled_counts])
    group_names = list(sampled_counts)
    total_by_token = count_matrix.sum(axis=0)
    total = float(total_by_token.sum())
    rows = []
    for group_idx, group_name in enumerate(group_names):
        group_counts = count_matrix[group_idx]
        group_total = float(group_counts.sum())
        background_counts = total_by_token - group_counts
        background_total = total - group_total
        values = []
        for token_idx, token in enumerate(tokens):
            group_count = float(group_counts[token_idx])
            group_freq = group_count / group_total if group_total else 0.0
            if group_count < min_count or group_freq < min_group_frequency:
                continue
            bg_freq = float(background_counts[token_idx] / background_total) if background_total else 0.0
            enrichment = (group_freq + pseudocount) / (bg_freq + pseudocount)
            values.append((token, int(group_count), group_freq, bg_freq, enrichment))
        protein, subtype = group_name.split("-", 1)
        for rank, (token, count, group_freq, bg_freq, enrichment) in enumerate(sorted(values, key=lambda x: x[4], reverse=True)[:top_n], start=1):
            rows.append(
                {
                    "bootstrap_id": bootstrap_id,
                    "tokenizer": tokenizer,
                    "protein": protein,
                    "subtype": subtype,
                    "protein_subtype": group_name,
                    "rank": rank,
                    "token": token,
                    "count": count,
                    "group_frequency": group_freq,
                    "background_frequency": bg_freq,
                    "enrichment_ratio": enrichment,
                }
            )
    return rows


def bootstrap_top_tokens(
    group_token_distributions: pd.DataFrame,
    n_bootstraps: int = 100,
    random_state: int | None = None,
    max_tokens_per_group: int | None = 200_000,
    top_n: int = 20,
    min_count: int = 100,
    min_group_frequency: float = 1e-5,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_state)
    rows = []
    for tokenizer, tok_df in group_token_distributions.groupby("tokenizer", dropna=False):
        for bootstrap_id in range(n_bootstraps):
            sampled_counts, tokens = _bootstrap_group_counts(tok_df, rng, max_tokens_per_group=max_tokens_per_group)
            rows.extend(
                _top_enriched_from_sampled_counts(
                    sampled_counts,
                    tokens,
                    tokenizer=tokenizer,
                    bootstrap_id=bootstrap_id,
                    top_n=top_n,
                    min_count=min_count,
                    min_group_frequency=min_group_frequency,
                )
            )
    return pd.DataFrame(rows)


def jaccard_index(left: Iterable[str], right: Iterable[str]) -> float:
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    if not union:
        return 1.0
    return len(left_set & right_set) / len(union)


def top_token_jaccard_stability(top_tokens_by_bootstrap: pd.DataFrame, reference_top_tokens: pd.DataFrame | None = None) -> pd.DataFrame:
    required = {"bootstrap_id", "tokenizer", "protein_subtype", "token"}
    missing = required - set(top_tokens_by_bootstrap.columns)
    if missing:
        raise ValueError(f"Missing top token columns: {', '.join(sorted(missing))}")
    if reference_top_tokens is None:
        reference_top_tokens = top_tokens_by_bootstrap.groupby(["tokenizer", "protein_subtype", "token"], dropna=False).size().reset_index(name="n")
    rows = []
    ref_sets = {
        key: set(group["token"].astype(str))
        for key, group in reference_top_tokens.groupby(["tokenizer", "protein_subtype"], dropna=False)
    }
    for key, group in top_tokens_by_bootstrap.groupby(["bootstrap_id", "tokenizer", "protein_subtype"], dropna=False):
        bootstrap_id, tokenizer, protein_subtype = key
        reference = ref_sets.get((tokenizer, protein_subtype), set())
        rows.append(
            {
                "bootstrap_id": bootstrap_id,
                "tokenizer": tokenizer,
                "protein_subtype": protein_subtype,
                "jaccard_vs_global_top": jaccard_index(group["token"].astype(str), reference),
                "n_bootstrap_tokens": int(group["token"].nunique()),
                "n_reference_tokens": int(len(reference)),
            }
        )
    return pd.DataFrame(rows)


def temporal_bins(
    metadata: pd.DataFrame,
    year_col: str,
    strategy: str = "fixed",
    n_bins: int = 4,
) -> pd.Series:
    years = pd.to_numeric(metadata[year_col], errors="coerce")
    labels = pd.Series("unknown", index=metadata.index, dtype="object")
    valid = years.notna()
    if not valid.any():
        return labels
    if strategy == "fixed":
        bins = [-np.inf, 2008, 2014, 2019, np.inf]
        names = ["pre-2009", "2009-2014", "2015-2019", "2020+"]
        labels.loc[valid] = pd.cut(years.loc[valid], bins=bins, labels=names, right=True).astype(str)
        return labels
    if strategy == "quantile":
        ranks = years.loc[valid].rank(method="first")
        q = min(n_bins, int(valid.sum()))
        labels.loc[valid] = pd.qcut(ranks, q=q, labels=[f"q{i + 1}" for i in range(q)]).astype(str)
        return labels
    raise ValueError(f"Unknown temporal bin strategy: {strategy}")


def token_metrics_by_time_window(token_stats: pd.DataFrame, metadata: pd.DataFrame | None, time_bins: pd.Series) -> pd.DataFrame:
    del metadata
    stats = token_stats.copy()
    stats["time_window"] = time_bins.reindex(stats.index).fillna("unknown").to_numpy()
    rows = []
    for key, group in stats.groupby(["tokenizer", "time_window", "protein_subtype"], dropna=False):
        tokenizer, time_window, protein_subtype = key
        rows.append(
            {
                "tokenizer": tokenizer,
                "time_window": time_window,
                "protein_subtype": protein_subtype,
                "n_sequences": int(group["internal_sequence_id"].nunique()) if "internal_sequence_id" in group else int(len(group)),
                "mean_token_entropy_bits": float(group["token_entropy_bits"].mean()),
                "mean_effective_vocab_size": float(group["effective_vocab_size"].mean()),
                "mean_cpg_token_fraction": float(group["cpg_token_fraction"].mean()),
                "mean_upa_token_fraction": float(group["upa_token_fraction"].mean()),
            }
        )
    return pd.DataFrame(rows)
