from __future__ import annotations

import numpy as np
import pandas as pd

from flugenome3d.token_stability import (
    bootstrap_token_metrics,
    jaccard_index,
    stratified_bootstrap_indices,
    summarize_bootstrap_ci,
    temporal_bins,
    top_token_jaccard_stability,
)


def test_stratified_bootstrap_indices_preserves_groups() -> None:
    metadata = pd.DataFrame({"group": ["A", "A", "B", "B", "C", "C"]})
    idx = stratified_bootstrap_indices(metadata, "group", n_per_group=2, random_state=1)
    sampled = metadata.loc[idx]
    assert set(sampled["group"]) == {"A", "B", "C"}
    assert sampled.groupby("group").size().to_dict() == {"A": 2, "B": 2, "C": 2}


def test_summarize_bootstrap_ci_is_finite() -> None:
    summary = summarize_bootstrap_ci([1.0, 2.0, 3.0])
    assert summary["n"] == 3
    assert np.isfinite(summary["mean"])
    assert np.isfinite(summary["ci_lower"])
    assert np.isfinite(summary["ci_upper"])


def test_jaccard_index_cases() -> None:
    assert jaccard_index({"AAA", "CCC"}, {"CCC", "AAA"}) == 1.0
    assert jaccard_index({"AAA"}, {"CCC"}) == 0.0


def test_top_token_jaccard_stability_against_reference() -> None:
    boot = pd.DataFrame(
        [
            {"bootstrap_id": 0, "tokenizer": "tok", "protein_subtype": "HA-H1N1", "token": "AAA"},
            {"bootstrap_id": 0, "tokenizer": "tok", "protein_subtype": "HA-H1N1", "token": "CCC"},
        ]
    )
    ref = pd.DataFrame(
        [
            {"tokenizer": "tok", "protein_subtype": "HA-H1N1", "token": "AAA"},
            {"tokenizer": "tok", "protein_subtype": "HA-H1N1", "token": "CCC"},
        ]
    )
    out = top_token_jaccard_stability(boot, ref)
    assert out["jaccard_vs_global_top"].iloc[0] == 1.0


def test_temporal_bins_handles_missing_years() -> None:
    metadata = pd.DataFrame({"year": [2007, 2010, 2018, 2021, None]})
    bins = temporal_bins(metadata, "year", strategy="fixed")
    assert bins.tolist() == ["pre-2009", "2009-2014", "2015-2019", "2020+", "unknown"]


def test_bootstrap_token_metrics_runs_on_small_table() -> None:
    token_stats = pd.DataFrame(
        {
            "internal_sequence_id": ["a1", "a2", "b1", "b2"],
            "tokenizer": ["tok", "tok", "tok", "tok"],
            "protein_subtype": ["A", "A", "B", "B"],
            "token_entropy_bits": [1.0, 1.2, 2.0, 2.2],
            "effective_vocab_size": [2.0, 2.3, 4.0, 4.5],
            "cpg_token_fraction": [0.1, 0.2, 0.3, 0.4],
            "upa_token_fraction": [0.2, 0.2, 0.1, 0.1],
        }
    )
    out = bootstrap_token_metrics(token_stats, None, "tok", n_bootstraps=3, n_per_group=2, random_state=2)
    assert set(out["protein_subtype"]) == {"A", "B"}
    assert out["bootstrap_id"].nunique() == 3
