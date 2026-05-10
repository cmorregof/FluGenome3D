from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIGURE_DIR = PROJECT / "results" / "figures"
REPORT_DIR = PROJECT / "reports"
APP_DATA = PROJECT / "app" / "data"


def test_phase7_9_public_outputs_exist() -> None:
    tables = [
        "phase7_antigenlm_spearman_summary.csv",
        "phase7_antigenlm_pca_summary.csv",
        "phase7_antigenlm_temporal_locality_summary.csv",
        "phase7_antigenlm_clade_enrichment_summary.csv",
        "phase8_representation_family_comparison.csv",
        "phase9_structure_mapping_qc.csv",
        "phase9_structure_signal_catalog.csv",
    ]
    figures = [
        "fig29_antigenlm_pca_by_subtype.png",
        "fig30_antigenlm_geometry_summary.png",
        "fig31_antigenlm_export_pca_variance.png",
        "fig32_representation_atlas_coverage.png",
        "fig33_structure_mapping_qc.png",
        "fig34_structure_signal_catalog.png",
    ]
    reports = [
        "phase7_antigenlm_bridge_report.md",
        "phase8_latent_atlas_report.md",
        "phase9_structure_mapping_report.md",
    ]
    assert [name for name in tables if not (TABLE_DIR / name).exists()] == []
    assert [name for name in figures if not (FIGURE_DIR / name).exists()] == []
    assert [name for name in reports if not (REPORT_DIR / name).exists()] == []


def test_phase7_9_safe_json_exports_exist_and_are_nonempty() -> None:
    latent = json.loads((APP_DATA / "antigenlm_latent_atlas.safe.json").read_text())
    mapping = json.loads((APP_DATA / "structure_mapping.safe.json").read_text())
    assert latent["projection"]["n_exported_points"] > 0
    assert latent["projection"]["points"]
    assert mapping["mapping_qc"]
    assert mapping["signal_catalog"]


def test_phase7_9_public_tables_do_not_include_sequence_columns() -> None:
    forbidden = {"sequence", "raw_sequence", "refined_sequence", "ha_sequence", "na_sequence", "epi_isl", "strain_name"}
    for path in list(TABLE_DIR.glob("phase7_*.csv")) + list(TABLE_DIR.glob("phase8_*.csv")) + list(TABLE_DIR.glob("phase9_*.csv")):
        df = pd.read_csv(path, nrows=1, keep_default_na=False)
        assert forbidden.isdisjoint(set(df.columns))


def test_phase7_9_reports_tables_and_exports_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    paths = (
        list(REPORT_DIR.glob("phase7_*.md"))
        + list(REPORT_DIR.glob("phase8_*.md"))
        + list(REPORT_DIR.glob("phase9_*.md"))
        + list(TABLE_DIR.glob("phase7_*.csv"))
        + list(TABLE_DIR.glob("phase8_*.csv"))
        + list(TABLE_DIR.glob("phase9_*.csv"))
        + [APP_DATA / "antigenlm_latent_atlas.safe.json", APP_DATA / "structure_mapping.safe.json"]
    )
    leaks = [str(path.relative_to(PROJECT)) for path in paths if pattern.search(path.read_text(encoding="utf-8", errors="replace"))]
    assert leaks == []
