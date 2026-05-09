from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .qc import normalize_subtype
from .utils import load_yaml


@dataclass(frozen=True)
class LocalFileSet:
    parent_repo_root: Path
    local_data_roots: tuple[Path, ...]
    rich_metadata_path: Path | None
    dedup_metadata_path: Path | None
    paired_dataset_paths: tuple[Path, ...]


def _resolve_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def load_local_file_set(config_path: str | Path) -> LocalFileSet:
    config_path = Path(config_path)
    cfg = load_yaml(config_path)
    project_root = _resolve_path(cfg.get("project_root", "."), Path.cwd())
    parent_root = _resolve_path(cfg.get("parent_repo_root", ".."), project_root)

    roots = tuple(_resolve_path(p, project_root) for p in cfg.get("local_data_roots", []))
    metadata_paths = cfg.get("metadata_paths", {}) or {}
    rich_metadata = metadata_paths.get("rich_metadata")
    dedup_metadata = metadata_paths.get("dedup_metadata")
    paired_paths = cfg.get("paired_dataset_paths", []) or []

    return LocalFileSet(
        parent_repo_root=parent_root,
        local_data_roots=roots,
        rich_metadata_path=_resolve_path(rich_metadata, project_root) if rich_metadata else None,
        dedup_metadata_path=_resolve_path(dedup_metadata, project_root) if dedup_metadata else None,
        paired_dataset_paths=tuple(_resolve_path(p, project_root) for p in paired_paths),
    )


def detect_relevant_files(file_set: LocalFileSet) -> dict[str, list[Path] | Path | None]:
    fasta_files: list[Path] = []
    for root in file_set.local_data_roots:
        if root.is_file() and root.suffix.lower() in {".fa", ".fasta", ".fna"}:
            fasta_files.append(root)
        elif root.exists():
            fasta_files.extend(sorted(p for p in root.rglob("*") if p.suffix.lower() in {".fa", ".fasta", ".fna"}))

    paired = list(file_set.paired_dataset_paths)
    if not paired:
        paired = sorted(file_set.parent_repo_root.glob("data/processed_gisaid/dataset_*.json"))

    rich = file_set.rich_metadata_path
    if rich is None:
        candidates = sorted(file_set.parent_repo_root.glob("data/gisaid_metadata_private/*combined.csv"))
        rich = candidates[0] if candidates else None

    dedup = file_set.dedup_metadata_path
    if dedup is None:
        candidates = sorted(file_set.parent_repo_root.glob("data/gisaid_metadata_private/*joined_dedup_cache.csv"))
        dedup = candidates[0] if candidates else None

    return {
        "fasta_files": fasta_files,
        "paired_dataset_paths": paired,
        "rich_metadata_path": rich,
        "dedup_metadata_path": dedup,
    }


def validate_columns(df: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {', '.join(missing)}")


def _split_location(value: object) -> tuple[str | None, str | None]:
    text = "" if value is None else str(value).strip()
    if not text:
        return None, None
    parts = [part.strip() for part in text.split("/") if part.strip()]
    if len(parts) >= 2:
        return parts[0], parts[1]
    return parts[0], None


def major_clade_from(value: object) -> str | None:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() in {"nan", "none"}:
        return None
    if text.lower() == "unassigned":
        return "unassigned"
    parts = text.split(".")
    if text.startswith("6B.1A"):
        return "6B.1A"
    if text.startswith("3C.2a1b"):
        return "3C.2a1b"
    return ".".join(parts[:2]) if len(parts) > 1 else text


def load_rich_metadata(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    required = {"Isolate_Id", "Subtype", "Host", "Location", "Collection_Date", "Clade"}
    validate_columns(df, required, "rich metadata")
    location = df["Location"].map(_split_location)
    out = pd.DataFrame(
        {
            "epi_isl": df["Isolate_Id"].astype(str),
            "metadata_subtype": df["Subtype"].map(normalize_subtype),
            "host": df["Host"].astype(str),
            "host_normalized": df["Host"].astype(str).str.strip().str.lower(),
            "location_raw": df["Location"].astype(str),
            "region": [x[0] for x in location],
            "country": [x[1] for x in location],
            "metadata_collection_date": df["Collection_Date"].astype(str),
            "metadata_year": pd.to_numeric(df["Collection_Date"].astype(str).str.extract(r"^(\d{4})")[0], errors="coerce"),
            "clade_raw_rich": df["Clade"].astype(str),
            "major_clade_rich": df["Clade"].map(major_clade_from),
        }
    )
    out["host_is_human"] = out["host_normalized"].eq("human")
    return out.drop_duplicates(subset=["epi_isl"], keep="first")


def load_dedup_metadata(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"epi_isl", "subtype", "year", "month", "matched", "clade", "major_clade"}
    validate_columns(df, required, "dedup metadata")
    out = df.rename(
        columns={
            "subtype": "dedup_subtype",
            "year": "dedup_year",
            "month": "dedup_month",
            "matched": "dedup_matched",
            "clade_raw": "clade_raw_dedup",
            "clade": "clade_dedup",
            "major_clade": "major_clade_dedup",
        }
    ).copy()
    out["epi_isl"] = out["epi_isl"].astype(str)
    out["dedup_subtype"] = out["dedup_subtype"].map(normalize_subtype)
    return out.drop_duplicates(subset=["epi_isl"], keep="first")


def load_paired_datasets(paths: list[str | Path] | tuple[Path, ...]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    required = {"epi_isl", "subtype", "year", "month", "day", "ha_sequence", "na_sequence"}
    for path in paths:
        p = Path(path)
        with p.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        records = payload.get("paired_strains", [])
        if records:
            missing = required - set(records[0])
            if missing:
                raise ValueError(f"{p} paired_strains missing required keys: {', '.join(sorted(missing))}")
        for record in records:
            row = dict(record)
            row["source_file"] = p.name
            rows.append(row)
    if not rows:
        raise ValueError("No paired HA/NA records found in configured paired_dataset_paths.")
    df = pd.DataFrame(rows)
    validate_columns(df, required, "paired datasets")
    df["epi_isl"] = df["epi_isl"].astype(str)
    df["subtype"] = df["subtype"].map(normalize_subtype)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["month"] = pd.to_numeric(df["month"], errors="coerce").astype("Int64")
    df["day"] = pd.to_numeric(df["day"], errors="coerce").astype("Int64")
    return df
