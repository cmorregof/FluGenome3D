from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TOKEN_DIR = PROJECT / "data" / "processed" / "tokenization"
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"


def test_phase5_local_outputs_exist() -> None:
    expected = [
        "mvp_token_stats.parquet",
        "mvp_cds_token_stats.parquet",
        "group_token_distributions.parquet",
    ]
    assert [name for name in expected if not (TOKEN_DIR / name).exists()] == []
    raw = pd.read_parquet(TOKEN_DIR / "mvp_token_stats.parquet")
    cds = pd.read_parquet(TOKEN_DIR / "mvp_cds_token_stats.parquet")
    assert raw["internal_sequence_id"].nunique() > 0
    assert cds["internal_sequence_id"].nunique() > 0
    assert raw["tokenizer"].nunique() == 8
    assert cds["tokenizer"].nunique() == 6


def test_phase5_public_outputs_exist() -> None:
    tables = [
        "phase5_tokenizer_summary.csv",
        "phase5_token_entropy_by_group.csv",
        "phase5_effective_vocab_by_group.csv",
        "phase5_cpg_upa_token_summary.csv",
        "phase5_codon_boundary_crossing_summary.csv",
        "phase5_top_tokens_by_group.csv",
        "phase5_group_js_distances.csv",
    ]
    figures = [
        "fig18_token_entropy_by_tokenizer.png",
        "fig19_effective_vocab_by_group.png",
        "fig20_cpg_upa_token_fraction.png",
        "fig21_codon_boundary_crossing.png",
        "fig22_token_js_distance_heatmap.png",
        "fig23_top_tokens_by_group.png",
    ]
    assert [name for name in tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in figures if not (FIGURE_DIR / name).exists()] == []
    assert (REPORT_DIR / "phase5_tokenization_audit_report.md").exists()


def test_phase5_public_tables_do_not_include_sequence_columns() -> None:
    forbidden = {"sequence", "raw_sequence", "refined_sequence", "ha_sequence", "na_sequence"}
    for path in TABLE_DIR.glob("phase5_*.csv"):
        df = pd.read_csv(path, nrows=1, keep_default_na=False)
        assert forbidden.isdisjoint(set(df.columns))


def test_phase5_top_tokens_are_short() -> None:
    top = pd.read_csv(TABLE_DIR / "phase5_top_tokens_by_group.csv", keep_default_na=False)
    assert top["token"].astype(str).str.len().max() <= 6


def test_phase5_reports_and_tables_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = [REPORT_DIR / "phase5_tokenization_audit_report.md"] + list(TABLE_DIR.glob("phase5_*.csv"))
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
