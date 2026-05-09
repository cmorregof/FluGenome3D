from __future__ import annotations

import numpy as np
import pandas as pd

from flugenome3d.token_stats import (
    effective_vocab_size,
    group_token_distribution,
    jensen_shannon_distance_between_groups,
    token_counts,
    token_entropy,
    token_frequencies,
    tokens_per_kb,
)


def test_token_counts_and_frequencies() -> None:
    counts = token_counts(["AAA", "AAA", "CCC"])
    freqs = token_frequencies(["AAA", "AAA", "CCC"])
    assert counts["AAA"] == 2
    assert freqs["AAA"] == 2 / 3


def test_entropy_and_effective_vocab_are_finite() -> None:
    assert token_entropy([]) == 0.0
    entropy = token_entropy(["AAA", "CCC", "CCC", "GGG"])
    assert np.isfinite(entropy)
    assert effective_vocab_size(["AAA", "CCC", "CCC", "GGG"]) > 0


def test_tokens_per_kb() -> None:
    assert tokens_per_kb(["AAA", "CCC"], 1000) == 2.0


def test_group_distribution_and_js_symmetry() -> None:
    token_table = pd.DataFrame(
        [
            {"tokenizer": "tok", "protein": "HA", "subtype": "H1N1", "token": "AAA", "count": 8},
            {"tokenizer": "tok", "protein": "HA", "subtype": "H1N1", "token": "CCC", "count": 2},
            {"tokenizer": "tok", "protein": "NA", "subtype": "H1N1", "token": "AAA", "count": 1},
            {"tokenizer": "tok", "protein": "NA", "subtype": "H1N1", "token": "CCC", "count": 9},
        ]
    )
    dist = group_token_distribution(token_table)
    js = jensen_shannon_distance_between_groups(dist)
    ab = js[(js["group_a"] == "HA-H1N1") & (js["group_b"] == "NA-H1N1")]["js_distance"].iloc[0]
    ba = js[(js["group_a"] == "NA-H1N1") & (js["group_b"] == "HA-H1N1")]["js_distance"].iloc[0]
    aa = js[(js["group_a"] == "HA-H1N1") & (js["group_b"] == "HA-H1N1")]["js_distance"].iloc[0]
    assert ab == ba
    assert aa == 0.0
