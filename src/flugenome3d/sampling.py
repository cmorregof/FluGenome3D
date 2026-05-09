from __future__ import annotations

from math import ceil

import pandas as pd


def _sample_without_replacement(df: pd.DataFrame, n: int, random_state: int) -> pd.DataFrame:
    if n <= 0 or df.empty:
        return df.iloc[0:0].copy()
    if len(df) <= n:
        return df.copy()
    return df.sample(n=n, random_state=random_state)


def sample_by_subtype_year(
    df: pd.DataFrame,
    target_total: int,
    random_seed: int = 42,
    subtype_col: str = "subtype",
    year_col: str = "year",
) -> pd.DataFrame:
    if target_total <= 0:
        return df.iloc[0:0].copy()
    if len(df) <= target_total:
        return df.copy().reset_index(drop=True)

    subtypes = sorted(df[subtype_col].dropna().unique())
    if not subtypes:
        return _sample_without_replacement(df, target_total, random_seed).reset_index(drop=True)

    base_quota = target_total // len(subtypes)
    remainder = target_total % len(subtypes)
    selected_parts: list[pd.DataFrame] = []

    for idx, subtype in enumerate(subtypes):
        subtype_target = base_quota + (1 if idx < remainder else 0)
        group = df[df[subtype_col] == subtype]
        if group.empty:
            continue
        years = sorted(group[year_col].dropna().unique())
        if not years:
            selected_parts.append(_sample_without_replacement(group, subtype_target, random_seed + idx))
            continue

        cap_per_year = max(1, ceil(subtype_target / len(years)))
        subtype_parts: list[pd.DataFrame] = []
        selected_index: set[int] = set()
        for y_idx, year in enumerate(years):
            year_group = group[group[year_col] == year]
            take = min(cap_per_year, len(year_group))
            sampled = _sample_without_replacement(year_group, take, random_seed + idx * 1000 + y_idx)
            subtype_parts.append(sampled)
            selected_index.update(sampled.index.tolist())

        subtype_selected = pd.concat(subtype_parts, ignore_index=False) if subtype_parts else group.iloc[0:0]
        if len(subtype_selected) > subtype_target:
            subtype_selected = subtype_selected.sample(n=subtype_target, random_state=random_seed + idx)
        elif len(subtype_selected) < subtype_target:
            remaining = group.drop(index=list(selected_index), errors="ignore")
            fill = _sample_without_replacement(remaining, subtype_target - len(subtype_selected), random_seed + idx + 5000)
            subtype_selected = pd.concat([subtype_selected, fill], ignore_index=False)
        selected_parts.append(subtype_selected)

    selected = pd.concat(selected_parts, ignore_index=False) if selected_parts else df.iloc[0:0]
    if len(selected) > target_total:
        selected = selected.sample(n=target_total, random_state=random_seed)
    elif len(selected) < target_total:
        remaining = df.drop(index=selected.index, errors="ignore")
        fill = _sample_without_replacement(remaining, target_total - len(selected), random_seed + 9000)
        selected = pd.concat([selected, fill], ignore_index=False)
    return selected.sort_values([subtype_col, year_col, "internal_strain_id"]).reset_index(drop=True)
