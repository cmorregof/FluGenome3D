from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
STABILITY_DIR = PROJECT / "data" / "processed" / "tokenization_stability"
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"


def test_phase6_local_outputs_exist() -> None:
    expected = [
        "bootstrap_metrics.parquet",
        "bootstrap_js_distances.parquet",
        "bootstrap_top_tokens.parquet",
        "temporal_metrics.parquet",
    ]
    assert [name for name in expected if not (STABILITY_DIR / name).exists()] == []
    boot = pd.read_parquet(STABILITY_DIR / "bootstrap_metrics.parquet")
    assert boot["bootstrap_id"].nunique() == 100
    assert boot["tokenizer"].nunique() == 8


def test_phase6_public_outputs_exist() -> None:
    tables = [
        "phase6_bootstrap_metric_summary.csv",
        "phase6_js_distance_stability.csv",
        "phase6_top_token_jaccard_stability.csv",
        "phase6_temporal_token_summary.csv",
        "phase6_tokenizer_robustness_ranking.csv",
    ]
    figures = [
        "fig24_bootstrap_js_distance_ci.png",
        "fig25_token_entropy_stability.png",
        "fig26_top_token_jaccard_stability.png",
        "fig27_temporal_token_entropy.png",
        "fig28_tokenizer_robustness_ranking.png",
    ]
    assert [name for name in tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in figures if not (FIGURE_DIR / name).exists()] == []
    assert (REPORT_DIR / "phase6_tokenization_stability_report.md").exists()


def test_phase6_ranking_exists_and_is_nonempty() -> None:
    ranking = pd.read_csv(TABLE_DIR / "phase6_tokenizer_robustness_ranking.csv")
    assert not ranking.empty
    assert {"rank", "tokenizer", "robustness_score"}.issubset(ranking.columns)
    assert ranking["robustness_score"].between(0, 1).all()


def test_phase6_public_tables_do_not_include_sequence_columns() -> None:
    forbidden = {"sequence", "raw_sequence", "refined_sequence", "ha_sequence", "na_sequence"}
    for path in TABLE_DIR.glob("phase6_*.csv"):
        df = pd.read_csv(path, nrows=1, keep_default_na=False)
        assert forbidden.isdisjoint(set(df.columns))


def test_phase6_public_tokens_are_short_if_present() -> None:
    for path in TABLE_DIR.glob("phase6_*.csv"):
        df = pd.read_csv(path, keep_default_na=False)
        if "token" in df.columns:
            assert df["token"].astype(str).str.len().max() <= 6


def test_phase6_reports_and_tables_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = [REPORT_DIR / "phase6_tokenization_stability_report.md"] + list(TABLE_DIR.glob("phase6_*.csv"))
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
