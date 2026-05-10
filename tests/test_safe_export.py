from __future__ import annotations

import json
import re
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
APP_DATA = PROJECT / "app" / "data"


def test_safe_export_files_exist() -> None:
    expected = [
        "dataset_overview.safe.json",
        "representation_maps.safe.json",
        "metric_summaries.safe.json",
        "tokenization_summaries.safe.json",
        "stability_summaries.safe.json",
        "antigenlm_latent_atlas.safe.json",
        "structure_catalog.safe.json",
        "structure_mapping.safe.json",
        "claims_and_limits.safe.json",
        "data_governance.safe.json",
    ]
    assert [name for name in expected if not (APP_DATA / name).exists()] == []


def test_safe_exports_have_no_long_sequences() -> None:
    pattern = re.compile(r"[ACGTN]{80,}")
    leaks = []
    for path in APP_DATA.glob("*.safe.json"):
        if pattern.search(path.read_text(encoding="utf-8", errors="replace")):
            leaks.append(path.name)
    assert leaks == []


def test_safe_token_exports_are_short() -> None:
    payload = json.loads((APP_DATA / "tokenization_summaries.safe.json").read_text())
    tokens = [row["token"] for row in payload["top_tokens_by_group"] if "token" in row]
    assert tokens
    assert max(len(str(token)) for token in tokens) <= 6


def test_representation_points_use_safe_ids() -> None:
    payload = json.loads((APP_DATA / "representation_maps.safe.json").read_text())
    representation = payload["representations"][0]
    first_raw = representation["points"][0]
    if isinstance(first_raw, list):
        first = dict(zip(representation["point_schema"], first_raw, strict=True))
    else:
        first = first_raw
    assert str(first["id"]).startswith("pt_")
    assert "fg3d_" not in json.dumps(first)
    forbidden = {"sequence", "sequence_sha256", "accession", "isolate", "strain_name"}
    assert forbidden.isdisjoint(first.keys())


def test_antigenlm_points_use_safe_ids() -> None:
    payload = json.loads((APP_DATA / "antigenlm_latent_atlas.safe.json").read_text())
    projection = payload["projection"]
    first = dict(zip(projection["point_schema"], projection["points"][0], strict=True))
    assert str(first["id"]).startswith("lm_")
    assert "EPI_ISL" not in json.dumps(first)
    forbidden = {"sequence", "sequence_sha256", "accession", "isolate", "strain_name", "epi_isl"}
    assert forbidden.isdisjoint(first.keys())
