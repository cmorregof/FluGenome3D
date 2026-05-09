from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO

SEQUENCE_EXTS = {".fa", ".fasta", ".fna", ".ffn", ".faa"}
TABLE_EXTS = {".csv", ".tsv", ".txt", ".parquet"}
JSON_EXTS = {".json", ".jsonl"}
COMPRESSED_EXTS = {".gz"}


@dataclass
class FileInventoryRecord:
    path: str
    suffix: str
    size_bytes: int
    kind: str
    n_records: int | None = None
    columns: str | None = None
    notes: str | None = None


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def classify_file(path: Path) -> str:
    suffixes = [s.lower() for s in path.suffixes]
    effective_suffix = suffixes[-2] if suffixes and suffixes[-1] in COMPRESSED_EXTS and len(suffixes) > 1 else path.suffix.lower()
    if effective_suffix in SEQUENCE_EXTS:
        return "sequence_fasta"
    if effective_suffix in TABLE_EXTS:
        return "table"
    if effective_suffix in JSON_EXTS:
        return "json"
    return "other"


def inspect_fasta(path: Path, max_records: int = 5000) -> tuple[int, str]:
    count = 0
    lengths: list[int] = []
    try:
        with _open_text(path) as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                count += 1
                lengths.append(len(rec.seq))
                if count >= max_records:
                    break
        note = ""
        if lengths:
            note = f"min_len={min(lengths)};max_len={max(lengths)};sampled={len(lengths)}"
        return count, note
    except Exception as exc:  # pragma: no cover - defensive inventory
        return 0, f"FASTA parse error: {exc}"


def inspect_table(path: Path) -> tuple[int | None, str | None, str]:
    try:
        if path.suffix.lower() == ".parquet":
            df = pd.read_parquet(path)
            return len(df), ",".join(map(str, df.columns[:80])), "parquet"
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        # sniff when possible
        with _open_text(path) as fh:
            sample = fh.read(4096)
            try:
                dialect = csv.Sniffer().sniff(sample)
                sep = dialect.delimiter
            except Exception:
                pass
        df = pd.read_csv(path, sep=sep, nrows=50)
        return None, ",".join(map(str, df.columns[:80])), f"sampled_rows=50;sep={repr(sep)}"
    except Exception as exc:  # pragma: no cover
        return None, None, f"table parse error: {exc}"


def inspect_json(path: Path) -> tuple[int | None, str | None, str]:
    try:
        with _open_text(path) as fh:
            first = fh.readline().strip()
        if not first:
            return 0, None, "empty json/jsonl"
        obj = json.loads(first)
        if isinstance(obj, dict):
            return None, ",".join(list(obj.keys())[:80]), "json/jsonl first object keys"
        return None, None, f"json first type={type(obj).__name__}"
    except Exception as exc:  # pragma: no cover
        return None, None, f"json parse error: {exc}"


def iter_candidate_files(roots: Iterable[str | Path], max_size_mb: float | None = None) -> Iterable[Path]:
    for root in roots:
        root_path = Path(root).expanduser()
        if not root_path.exists():
            continue
        if root_path.is_file():
            yield root_path
            continue
        for path in root_path.rglob("*"):
            if path.is_file():
                if max_size_mb is not None and path.stat().st_size > max_size_mb * 1024 * 1024:
                    # Still yield sequence/table files; note size in report.
                    pass
                yield path


def build_inventory(roots: Iterable[str | Path]) -> pd.DataFrame:
    records: list[FileInventoryRecord] = []
    for path in iter_candidate_files(roots):
        kind = classify_file(path)
        suffix = "".join(path.suffixes).lower()
        rec = FileInventoryRecord(
            path=str(path),
            suffix=suffix,
            size_bytes=path.stat().st_size,
            kind=kind,
        )
        if kind == "sequence_fasta":
            n, notes = inspect_fasta(path)
            rec.n_records = n
            rec.notes = notes
        elif kind == "table":
            n, cols, notes = inspect_table(path)
            rec.n_records = n
            rec.columns = cols
            rec.notes = notes
        elif kind == "json":
            n, cols, notes = inspect_json(path)
            rec.n_records = n
            rec.columns = cols
            rec.notes = notes
        records.append(rec)
    return pd.DataFrame([asdict(r) for r in records])
