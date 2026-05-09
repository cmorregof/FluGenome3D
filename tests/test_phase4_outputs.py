from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from scipy import sparse


PROJECT = Path(__file__).resolve().parents[1]
REP_DIR = PROJECT / "data" / "processed" / "representations"
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"


def test_phase4_local_matrices_exist() -> None:
    expected = [
        "mvp_kmer3_freq.npz",
        "mvp_kmer4_freq.npz",
        "mvp_kmer5_freq.npz",
        "mvp_kmer3_tfidf.npz",
        "mvp_kmer4_tfidf.npz",
        "mvp_codon_freq.npz",
        "mvp_rscu.npz",
    ]
    missing = [name for name in expected if not (REP_DIR / name).exists()]
    assert missing == []
    matrix = sparse.load_npz(REP_DIR / "mvp_kmer3_freq.npz")
    assert matrix.shape[0] > 0
    assert matrix.shape[1] == 64


def test_phase4_public_outputs_exist() -> None:
    tables = [
        "phase4_representation_summary.csv",
        "phase4_group_centroid_distances.csv",
        "phase4_silhouette_scores.csv",
        "phase4_top_kmers_by_group.csv",
        "phase4_top_codons_by_group.csv",
        "phase4_top_rscu_by_group.csv",
    ]
    figures = [
        "fig12_kmer_pca_by_group.png",
        "fig13_kmer_umap_by_group.png",
        "fig14_codon_pca_by_group.png",
        "fig15_rscu_pca_by_group.png",
        "fig16_group_centroid_distance_heatmap.png",
        "fig17_representation_silhouette_comparison.png",
    ]
    assert [name for name in tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in figures if not (FIGURE_DIR / name).exists()] == []
    assert (REPORT_DIR / "phase4_representation_audit_report.md").exists()


def test_phase4_public_tables_do_not_include_sequence_columns() -> None:
    forbidden = {"sequence", "raw_sequence", "refined_sequence", "ha_sequence", "na_sequence"}
    for path in TABLE_DIR.glob("phase4_*.csv"):
        df = pd.read_csv(path, nrows=1, keep_default_na=False)
        assert forbidden.isdisjoint(set(df.columns))


def test_phase4_reports_and_tables_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = list(REPORT_DIR.glob("*.md")) + list(TABLE_DIR.glob("*.csv"))
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
