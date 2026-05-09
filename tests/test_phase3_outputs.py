from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
PANEL_DIR = PROJECT / "data" / "processed" / "panels"
METRICS_DIR = PROJECT / "data" / "processed" / "metrics"
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"


def test_phase3_outputs_exist() -> None:
    expected_panels = [
        "mvp_cds_strict_panel.parquet",
        "mvp_cds_rescued_panel.parquet",
        "mvp_cds_refined_panel.parquet",
    ]
    expected_tables = [
        "phase3_cds_qc_failure_breakdown.csv",
        "phase3_length_mod3_by_group.csv",
        "phase3_internal_stop_distribution.csv",
        "phase3_refined_codon_usage_summary.csv",
        "phase3_refined_rscu_summary.csv",
        "phase3_refined_translation_qc_summary.csv",
    ]
    expected_figures = [
        "fig8_cds_qc_failure_breakdown.png",
        "fig9_length_mod3_by_group.png",
        "fig10_refined_rscu_heatmap.png",
    ]
    assert [name for name in expected_panels if not (PANEL_DIR / name).exists()] == []
    assert [name for name in expected_tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in expected_figures if not (FIGURE_DIR / name).exists()] == []
    assert (METRICS_DIR / "mvp_cds_refined_codon_metrics.parquet").exists()
    assert (REPORT_DIR / "phase3_cds_refinement_report.md").exists()


def test_phase3_refined_panel_non_empty_and_local_metrics_do_not_store_sequences() -> None:
    refined = pd.read_parquet(PANEL_DIR / "mvp_cds_refined_panel.parquet")
    metrics = pd.read_parquet(METRICS_DIR / "mvp_cds_refined_codon_metrics.parquet")
    assert len(refined) > 0
    assert len(metrics) > 0
    assert "refined_sequence" in refined.columns
    assert "raw_sequence" in refined.columns
    assert "refined_sequence" not in metrics.columns
    assert "raw_sequence" not in metrics.columns


def test_phase3_public_tables_do_not_include_sequence_columns() -> None:
    forbidden = {"sequence", "raw_sequence", "refined_sequence", "ha_sequence", "na_sequence"}
    for path in TABLE_DIR.glob("phase3_*.csv"):
        df = pd.read_csv(path, nrows=1, keep_default_na=False)
        assert forbidden.isdisjoint(set(df.columns))


def test_phase3_public_text_outputs_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = list(REPORT_DIR.glob("*.md")) + list(TABLE_DIR.glob("*.csv"))
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
