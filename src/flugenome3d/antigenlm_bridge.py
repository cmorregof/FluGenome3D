from __future__ import annotations

import hashlib
import json
import pickle
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


LATENT_SAFE_SALT = "flugenome3d-antigenlm-latent-v1"


def safe_latent_id(value: str, salt: str = LATENT_SAFE_SALT) -> str:
    digest = hashlib.sha256(f"{salt}::{value}".encode("utf-8")).hexdigest()
    return f"lm_{digest[:18]}"


def year_bin(year: Any) -> str:
    try:
        value = int(year)
    except (TypeError, ValueError):
        return "unknown"
    if value <= 2008:
        return "pre-2009"
    if value <= 2014:
        return "2009-2014"
    if value <= 2019:
        return "2015-2019"
    return "2020+"


def finite_number(value: Any, digits: int = 6) -> float | int | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    if number.is_integer():
        return int(number)
    return round(number, digits)


def load_json(path: str | Path) -> dict[str, Any]:
    full = Path(path)
    if not full.exists():
        return {}
    return json.loads(full.read_text(encoding="utf-8"))


def load_cache(path: str | Path) -> dict[str, Any]:
    with Path(path).open("rb") as handle:
        payload = pickle.load(handle)
    required = {"embeddings", "years", "months", "types", "records"}
    missing = required - set(payload)
    if missing:
        raise ValueError(f"Invalid AntigenLM cache, missing keys: {sorted(missing)}")
    n = len(payload["embeddings"])
    for key in ("years", "months", "types", "records"):
        if len(payload[key]) != n:
            raise ValueError(f"AntigenLM cache is not aligned: {key} has {len(payload[key])}, embeddings has {n}")
    return payload


def cache_summary(payload: dict[str, Any]) -> dict[str, Any]:
    embeddings = np.asarray(payload["embeddings"])
    types = np.asarray(payload["types"]).astype(str)
    years = np.asarray(payload["years"])
    metadata = payload.get("metadata", {})
    checkpoint = metadata.get("checkpoint", {}) if isinstance(metadata, dict) else {}
    subtype_counts = {subtype: int(np.sum(types == subtype)) for subtype in sorted(set(types))}
    return {
        "n_records": int(embeddings.shape[0]),
        "embedding_dim": int(embeddings.shape[1]) if embeddings.ndim == 2 else None,
        "subtype_counts": subtype_counts,
        "year_min": int(np.nanmin(years)) if len(years) else None,
        "year_max": int(np.nanmax(years)) if len(years) else None,
        "checkpoint_sha256": checkpoint.get("sha256", "not_recorded"),
        "max_seq_length": metadata.get("max_seq_length") if isinstance(metadata, dict) else None,
        "embedding_batch_size": metadata.get("embedding_batch_size") if isinstance(metadata, dict) else None,
        "source": "parent AntigenLM embedding cache; sequences are not exported",
    }


def stratified_sample_indices(metadata: pd.DataFrame, max_points: int, random_state: int = 42) -> np.ndarray:
    if len(metadata) <= max_points:
        return metadata.index.to_numpy()
    rng = np.random.default_rng(random_state)
    selected: list[int] = []
    strata = metadata.assign(stratum=metadata["subtype"].astype(str) + "|" + metadata["year_bin"].astype(str))
    groups = sorted(strata["stratum"].dropna().unique())
    per_group = max(1, max_points // max(1, len(groups)))
    for group in groups:
        idx = strata.index[strata["stratum"] == group].to_numpy()
        take = min(per_group, len(idx))
        selected.extend(rng.choice(idx, size=take, replace=False).tolist())
    if len(selected) < max_points:
        remaining = np.setdiff1d(metadata.index.to_numpy(), np.asarray(selected), assume_unique=False)
        take = min(max_points - len(selected), len(remaining))
        if take:
            selected.extend(rng.choice(remaining, size=take, replace=False).tolist())
    return np.asarray(sorted(selected), dtype=int)


def build_latent_pca_points(
    cache_path: str | Path,
    max_points: int = 30000,
    random_state: int = 42,
    n_components: int = 3,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    payload = load_cache(cache_path)
    embeddings = np.asarray(payload["embeddings"], dtype=np.float32)
    types = np.asarray(payload["types"]).astype(str)
    years = np.asarray(payload["years"])
    months = np.asarray(payload["months"])
    records = payload["records"]

    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=random_state)
    coords = pca.fit_transform(embeddings)

    meta = pd.DataFrame(
        {
            "cache_index": np.arange(len(embeddings), dtype=int),
            "subtype": types,
            "year": years,
            "month": months,
        }
    )
    meta["year_bin"] = meta["year"].map(year_bin)
    meta["source"] = "AntigenLM HA+NA embedding"
    meta["representation"] = "antigenlm_full_pca"
    meta["safe_id"] = [
        safe_latent_id(f"{idx}|{record.get('epi_isl', '')}|{record.get('subtype', subtype)}")
        for idx, record, subtype in zip(meta["cache_index"], records, types, strict=True)
    ]
    for axis in range(n_components):
        meta[f"axis{axis + 1}"] = coords[:, axis]

    sample_idx = stratified_sample_indices(meta, max_points=max_points, random_state=random_state)
    points = meta.loc[sample_idx].copy()
    summary = {
        "cache": cache_summary(payload),
        "projection": "pca_3d",
        "pca_explained_variance": [finite_number(value, 6) for value in pca.explained_variance_ratio_.tolist()],
        "n_source_points": int(len(meta)),
        "n_exported_points": int(len(points)),
        "sampling": {
            "method": "subtype_year_bin_stratified",
            "max_points": int(max_points),
            "random_seed": int(random_state),
        },
    }
    return points, summary


def summarize_spearman(metrics: dict[str, Any]) -> pd.DataFrame:
    rows = pd.DataFrame(metrics.get("spearman", []))
    if rows.empty:
        return pd.DataFrame(columns=["metric", "subtype", "rho_mean", "rho_sd", "valid_pairs_mean", "n_runs"])
    return (
        rows.groupby(["metric", "subtype"], dropna=False)
        .agg(
            rho_mean=("rho", "mean"),
            rho_sd=("rho", "std"),
            valid_pairs_mean=("valid_pairs", "mean"),
            omitted_pairs_mean=("omitted_pairs", "mean"),
            n_runs=("rho", "size"),
        )
        .reset_index()
    )


def summarize_pca(metrics: dict[str, Any]) -> pd.DataFrame:
    pca = metrics.get("pca", {})
    rows = []
    for group, payload in pca.items():
        rows.append(
            {
                "group": group,
                "n": payload.get("n"),
                "n80": payload.get("n80"),
                "n90": payload.get("n90"),
                "n95": payload.get("n95"),
                "n99": payload.get("n99"),
                "participation_ratio": payload.get("participation_ratio"),
                "top1_evr": (payload.get("explained_variance_ratio") or [None])[0],
                "top2_evr": (payload.get("explained_variance_ratio") or [None, None])[1],
            }
        )
    return pd.DataFrame(rows)


def summarize_twonn(metrics: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(metrics.get("twonn", {}).get("rows", []))


def summarize_temporal_locality(metrics: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(metrics.get("temporal_locality", {}).get("rows", []))


def summarize_clade_enrichment(clade: dict[str, Any]) -> pd.DataFrame:
    rows = pd.DataFrame(clade.get("clade_enrichment", []))
    if rows.empty:
        return pd.DataFrame()
    keep = [
        "subtype",
        "label",
        "k",
        "classes",
        "n_labeled",
        "mean_precision",
        "random_baseline",
        "enrichment_vs_random",
        "permutation_p05",
        "permutation_p95",
    ]
    return rows[[column for column in keep if column in rows.columns]].copy()


def summarize_random_baseline(random_baseline: dict[str, Any]) -> pd.DataFrame:
    rows = pd.DataFrame(random_baseline.get("aggregate", []))
    if rows.empty:
        return pd.DataFrame()
    return rows.copy()


def write_public_tables(
    metrics_path: str | Path,
    clade_path: str | Path,
    random_path: str | Path,
    out_dir: str | Path,
) -> dict[str, pd.DataFrame]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    metrics = load_json(metrics_path)
    clade = load_json(clade_path)
    random_baseline = load_json(random_path)
    tables = {
        "phase7_antigenlm_spearman_summary.csv": summarize_spearman(metrics),
        "phase7_antigenlm_pca_summary.csv": summarize_pca(metrics),
        "phase7_antigenlm_twonn_summary.csv": summarize_twonn(metrics),
        "phase7_antigenlm_temporal_locality_summary.csv": summarize_temporal_locality(metrics),
        "phase7_antigenlm_clade_enrichment_summary.csv": summarize_clade_enrichment(clade),
        "phase7_random_embedding_baseline_summary.csv": summarize_random_baseline(random_baseline),
    }
    for filename, table in tables.items():
        table.to_csv(out / filename, index=False)
    return tables


def safe_point_records(points: pd.DataFrame) -> list[list[Any]]:
    rows = []
    for row in points.itertuples(index=False):
        rows.append(
            [
                row.safe_id,
                finite_number(row.axis1, 7),
                finite_number(row.axis2, 7),
                finite_number(row.axis3, 7),
                row.subtype,
                year_bin(row.year),
                "AntigenLM HA+NA",
                "learned_latent",
            ]
        )
    return rows
