from __future__ import annotations

from collections import Counter
from math import log2
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

DNA_ALPHABET = set("ACGTUacgtu")
AMBIGUOUS = set("NRYKMSWBDHVnrykmswbdhv-")
DINUCLEOTIDES = tuple(a + b for a in "ACGT" for b in "ACGT")


def clean_sequence(seq: object) -> str:
    """Normalize a nucleotide sequence for descriptive compositional metrics."""
    return "".join(str(seq or "").upper().split()).replace("U", "T")


def normalize_rna_to_dna(seq: object) -> str:
    return clean_sequence(seq)


def sequence_length(seq: object) -> int:
    return len(clean_sequence(seq))


def ambiguous_fraction(seq: object) -> float:
    s = clean_sequence(seq)
    if not s:
        return 1.0
    return sum(1 for ch in s if ch not in {"A", "C", "G", "T"}) / len(s)


def gc_content(seq: object) -> float:
    s = clean_sequence(seq)
    valid = [ch for ch in s if ch in "ACGT"]
    if not valid:
        return np.nan
    return (valid.count("G") + valid.count("C")) / len(valid)


def gc_fraction(seq: object) -> float:
    """Backward-compatible alias for GC content."""
    return gc_content(seq)


def dinucleotide_counts(seq: object) -> Counter[str]:
    s = clean_sequence(seq)
    return Counter(s[i : i + 2] for i in range(max(0, len(s) - 1)))


def dinucleotide_frequencies(seq: object) -> dict[str, float]:
    counts = dinucleotide_counts(seq)
    total = sum(counts.values())
    if total == 0:
        return {dinuc: np.nan for dinuc in DINUCLEOTIDES}
    return {dinuc: counts[dinuc] / total for dinuc in DINUCLEOTIDES}


def mono_counts(seq: object) -> Counter[str]:
    s = clean_sequence(seq)
    return Counter(ch for ch in s if ch in "ACGT")


def dinucleotide_odds_ratio(seq: object, dinuc: str) -> float:
    normalized = dinuc.upper().replace("U", "T")
    if len(normalized) != 2 or any(base not in "ACGT" for base in normalized):
        return np.nan
    return dinucleotide_odds_ratios(seq).get(normalized, np.nan)


def dinucleotide_odds_ratios(seq: object) -> dict[str, float]:
    s = clean_sequence(seq)
    mono = mono_counts(s)
    valid_base_total = sum(mono.values())
    if valid_base_total == 0:
        return {dinuc: np.nan for dinuc in DINUCLEOTIDES}
    counts = dinucleotide_counts(s)
    valid_dinuc_total = sum(counts[d] for d in DINUCLEOTIDES)
    if valid_dinuc_total == 0:
        return {dinuc: np.nan for dinuc in DINUCLEOTIDES}
    ratios = {}
    for dinuc in DINUCLEOTIDES:
        x, y = dinuc
        fx = mono[x] / valid_base_total
        fy = mono[y] / valid_base_total
        fxy = counts[dinuc] / valid_dinuc_total
        denominator = fx * fy
        ratios[dinuc] = fxy / denominator if denominator else np.nan
    return ratios


def cpg_observed_expected(seq: object) -> float:
    return dinucleotide_odds_ratio(seq, "CG")


def upa_observed_expected(seq: object) -> float:
    # DNA alphabet: TA is the DNA proxy for RNA UpA.
    return dinucleotide_odds_ratio(seq, "TA")


def kmer_counts(seq: str, k: int, step: int = 1) -> Counter[str]:
    s = clean_sequence(seq)
    return Counter(
        s[i : i + k]
        for i in range(0, len(s) - k + 1, step)
        if set(s[i : i + k]) <= set("ACGT")
    )


def shannon_entropy(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total == 0:
        return np.nan
    probs = [v / total for v in counts.values() if v > 0]
    return -sum(p * log2(p) for p in probs)


def sequence_record_metrics(records: Iterable[SeqRecord]) -> pd.DataFrame:
    rows = []
    for rec in records:
        seq = str(rec.seq)
        rows.append(
            {
                "seq_id": rec.id,
                "description": rec.description,
                "length": sequence_length(seq),
                "ambiguous_fraction": ambiguous_fraction(seq),
                "gc_fraction": gc_content(seq),
                "cpg_oe": cpg_observed_expected(seq),
                "upa_oe": upa_observed_expected(seq),
                "k3_entropy": shannon_entropy(kmer_counts(seq, 3)),
                "k6_entropy": shannon_entropy(kmer_counts(seq, 6)),
            }
        )
    return pd.DataFrame(rows)


def translate_qc(seq: str) -> dict[str, object]:
    s = normalize_rna_to_dna(seq)
    trimmed = s[: len(s) - (len(s) % 3)]
    if not trimmed:
        return {"aa_length": 0, "internal_stop_count": np.nan, "terminal_stop": False}
    aa = str(Seq(trimmed).translate(to_stop=False))
    internal = aa[:-1].count("*") if len(aa) > 1 else 0
    return {"aa_length": len(aa), "internal_stop_count": internal, "terminal_stop": aa.endswith("*")}
