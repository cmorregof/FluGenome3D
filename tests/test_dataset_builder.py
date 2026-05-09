from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest


PROJECT = Path(__file__).resolve().parents[1]
PANEL_DIR = PROJECT / "data" / "processed" / "panels"
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"


def _read_panel(name: str) -> pd.DataFrame:
    path = PANEL_DIR / f"{name}_panel.parquet"
    if not path.exists():
        pytest.skip(f"{path} has not been generated")
    return pd.read_parquet(path)


def test_phase1_panels_are_non_empty_and_have_required_columns() -> None:
    required = {
        "internal_strain_id",
        "subtype",
        "year",
        "ha_sequence",
        "na_sequence",
        "ha_sha256",
        "na_sha256",
        "pair_sha256",
        "paired_pass_qc",
        "rich_metadata_available",
        "host_is_human",
    }
    for name in ["smoke", "mvp", "full"]:
        df = _read_panel(name)
        assert len(df) > 0
        assert required <= set(df.columns)


def test_deduplicated_panels_have_no_exact_ha_na_pair_duplicates() -> None:
    for name in ["smoke", "mvp", "full"]:
        df = _read_panel(name)
        assert df["pair_sha256"].is_unique


def test_phase1_aggregate_reports_and_figures_exist() -> None:
    expected_tables = [
        "phase1_dataset_summary.csv",
        "phase1_panel_summary.csv",
        "phase1_missingness_summary.csv",
        "phase1_year_subtype_counts.csv",
        "phase1_length_qc_summary.csv",
        "phase1_duplicate_summary.csv",
    ]
    expected_figures = [
        "fig1_dataset_overview.png",
        "fig2_year_subtype_distribution.png",
        "fig3_length_qc_distribution.png",
    ]
    missing_tables = [name for name in expected_tables if not (TABLE_DIR / name).exists()]
    missing_figures = [name for name in expected_figures if not (FIGURE_DIR / name).exists()]
    assert missing_tables == []
    assert missing_figures == []
    assert (REPORT_DIR / "phase1_dataset_report.md").exists()
