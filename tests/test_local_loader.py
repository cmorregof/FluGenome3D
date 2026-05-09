from __future__ import annotations

import re
from pathlib import Path

import pytest

from flugenome3d.local_loader import detect_relevant_files, load_local_file_set


PROJECT = Path(__file__).resolve().parents[1]


def _phase1_text_outputs() -> list[Path]:
    paths = list((PROJECT / "reports").glob("*.md"))
    paths.extend((PROJECT / "results" / "tables").glob("*.csv"))
    return paths


def test_local_paths_config_detects_phase1_sources_when_available() -> None:
    config = PROJECT / "config" / "local_paths.yml"
    if not config.exists():
        pytest.skip("local restricted config is not present")
    file_set = load_local_file_set(config)
    detected = detect_relevant_files(file_set)
    assert len(detected["paired_dataset_paths"]) >= 2
    assert detected["rich_metadata_path"] is not None
    assert detected["dedup_metadata_path"] is not None


def test_no_long_sequences_in_public_phase1_text_outputs() -> None:
    paths = _phase1_text_outputs()
    if not paths:
        pytest.skip("phase1 text outputs have not been generated")
    pattern = re.compile(r"[ACGTN]{80,}")
    leaks = []
    for path in paths:
        text = path.read_text(encoding="utf-8", errors="replace")
        if pattern.search(text):
            leaks.append(str(path.relative_to(PROJECT)))
    assert leaks == []
