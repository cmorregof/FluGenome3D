from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"
METRICS_DIR = PROJECT / "data" / "processed" / "metrics"


def test_phase2_aggregate_outputs_exist() -> None:
    expected_tables = [
        "phase2_gc_cpg_upa_summary.csv",
        "phase2_dinucleotide_odds_summary.csv",
        "phase2_codon_usage_summary.csv",
        "phase2_rscu_summary.csv",
        "phase2_kmer_entropy_summary.csv",
        "phase2_translation_qc_summary.csv",
    ]
    expected_figures = [
        "fig4_gc_cpg_upa_by_group.png",
        "fig5_dinucleotide_odds_heatmap.png",
        "fig6_codon_usage_rscu_heatmap.png",
        "fig7_kmer_entropy_by_group.png",
    ]
    assert [name for name in expected_tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in expected_figures if not (FIGURE_DIR / name).exists()] == []
    assert (REPORT_DIR / "phase2_sequence_context_report.md").exists()


def test_phase2_local_metric_outputs_exist_and_do_not_store_sequences() -> None:
    for panel in ["smoke", "mvp"]:
        seq_path = METRICS_DIR / f"{panel}_sequence_metrics.parquet"
        codon_path = METRICS_DIR / f"{panel}_codon_metrics.parquet"
        assert seq_path.exists()
        assert codon_path.exists()
        seq_df = pd.read_parquet(seq_path)
        codon_df = pd.read_parquet(codon_path)
        assert "sequence" not in seq_df.columns
        assert "sequence" not in codon_df.columns
        assert len(seq_df) > 0
        assert len(codon_df) > 0


def test_phase2_reports_and_tables_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = list(REPORT_DIR.glob("*.md")) + list(TABLE_DIR.glob("*.csv"))
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
