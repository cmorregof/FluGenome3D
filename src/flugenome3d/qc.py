from __future__ import annotations

import hashlib
import re
from typing import Any


def normalize_subtype(value: object) -> str:
    text = "" if value is None else str(value)
    match = re.search(r"H(\d+)N(\d+)", text, re.IGNORECASE)
    if match:
        return f"H{match.group(1)}N{match.group(2)}"
    return text.strip().upper() if text.strip() else "UNKNOWN"


def normalize_segment(value: object) -> str:
    text = "" if value is None else str(value).strip().upper()
    if text in {"4", "HA", "SEG4", "4_HA", "HEMAGGLUTININ"}:
        return "HA"
    if text in {"6", "NA", "SEG6", "6_NA", "NEURAMINIDASE"}:
        return "NA"
    return text or "UNKNOWN"


def normalize_sequence(seq: object) -> str:
    return str(seq or "").upper().replace("U", "T")


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sequence_sha256(seq: object) -> str:
    return sha256_text(normalize_sequence(seq))


def stable_internal_id(*parts: object, prefix: str = "fg3d") -> str:
    payload = "|".join("" if part is None else str(part) for part in parts)
    return f"{prefix}_{sha256_text(payload)[:16]}"


def ambiguous_fraction(seq: object) -> float:
    normalized = normalize_sequence(seq)
    if not normalized:
        return 1.0
    return sum(base not in {"A", "C", "G", "T"} for base in normalized) / len(normalized)


def non_acgtn_count(seq: object) -> int:
    normalized = normalize_sequence(seq)
    return sum(base not in {"A", "C", "G", "T", "N"} for base in normalized)


def expected_length_bounds(filters: dict[str, Any], segment: str) -> tuple[int | None, int | None]:
    qc = filters.get("sequence_qc", {})
    length_cfg = qc.get("length_nt", {}).get(segment, {})
    min_len = length_cfg.get("min")
    max_len = length_cfg.get("max")
    return (int(min_len) if min_len is not None else None, int(max_len) if max_len is not None else None)


def max_ambiguous_fraction(filters: dict[str, Any]) -> float:
    return float(filters.get("sequence_qc", {}).get("max_ambiguous_fraction", 0.01))


def sequence_qc(seq: object, segment: str, filters: dict[str, Any]) -> dict[str, object]:
    normalized = normalize_sequence(seq)
    length = len(normalized)
    min_len, max_len = expected_length_bounds(filters, segment)
    ambig = ambiguous_fraction(normalized)
    non_acgtn = non_acgtn_count(normalized)
    pass_min = min_len is None or length >= min_len
    pass_max = max_len is None or length <= max_len
    pass_ambig = ambig <= max_ambiguous_fraction(filters)
    pass_alphabet = non_acgtn == 0
    return {
        "length": length,
        "ambiguous_fraction": ambig,
        "non_acgtn_count": non_acgtn,
        "passes_length_min": pass_min,
        "passes_length_max": pass_max,
        "passes_length": pass_min and pass_max,
        "passes_ambiguity": pass_ambig,
        "passes_alphabet": pass_alphabet,
        "passes_qc": pass_min and pass_max and pass_ambig and pass_alphabet,
    }
