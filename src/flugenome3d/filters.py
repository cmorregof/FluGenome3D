from __future__ import annotations

import re
from collections import defaultdict

import pandas as pd
from Bio.SeqRecord import SeqRecord

from .sequence_metrics import ambiguous_fraction


def infer_segment_from_text(text: str) -> str | None:
    t = text.lower()
    if "hemagglutinin" in t or re.search(r"\bha\b", t):
        return "HA"
    if "neuraminidase" in t or re.search(r"\bna\b", t):
        return "NA"
    return None


def infer_subtype_from_text(text: str) -> str | None:
    m = re.search(r"H\d+N\d+", text.upper())
    return m.group(0) if m else None


def filter_records(records: list[SeqRecord], filters: dict) -> list[SeqRecord]:
    max_ambig = float(filters.get("max_ambiguous_fraction", 0.02))
    min_length = filters.get("min_length", {}) or {}
    keep = []
    seen = set()
    for rec in records:
        seq = str(rec.seq).upper()
        text = rec.description
        seg = infer_segment_from_text(text) or "unknown"
        min_len = int(min_length.get(seg, 0))
        if min_len and len(seq) < min_len:
            continue
        if ambiguous_fraction(seq) > max_ambig:
            continue
        if filters.get("drop_exact_duplicates", True):
            if seq in seen:
                continue
            seen.add(seq)
        keep.append(rec)
    return keep


def summarize_records(records: list[SeqRecord]) -> pd.DataFrame:
    rows = []
    for rec in records:
        rows.append(
            {
                "seq_id": rec.id,
                "description": rec.description,
                "segment_inferred": infer_segment_from_text(rec.description),
                "subtype_inferred": infer_subtype_from_text(rec.description),
                "length": len(rec.seq),
                "ambiguous_fraction": ambiguous_fraction(str(rec.seq)),
            }
        )
    return pd.DataFrame(rows)
