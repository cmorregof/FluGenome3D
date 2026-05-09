from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq

from .sequence_metrics import ambiguous_fraction, normalize_rna_to_dna

STANDARD_TABLE = CodonTable.unambiguous_dna_by_id[1]
CODON_TO_AA = {**STANDARD_TABLE.forward_table, **{codon: "*" for codon in STANDARD_TABLE.stop_codons}}
AA_TO_CODONS: dict[str, list[str]] = defaultdict(list)
for codon, aa in STANDARD_TABLE.forward_table.items():
    AA_TO_CODONS[aa].append(codon)
STOP_CODONS = set(STANDARD_TABLE.stop_codons)
SENSE_CODONS = tuple(sorted(STANDARD_TABLE.forward_table))
ALL_CODONS = tuple(sorted(CODON_TO_AA))


def is_dna_sequence(seq: object) -> bool:
    s = normalize_rna_to_dna(seq)
    return bool(s) and set(s) <= set("ACGTN")


def trim_to_codon_frame(seq: object) -> str:
    s = normalize_rna_to_dna(seq)
    n = len(s) - (len(s) % 3)
    return s[:n]


def has_valid_codon_frame(seq: object) -> bool:
    s = normalize_rna_to_dna(seq)
    return bool(s) and len(s) % 3 == 0


def translate_sequence(seq: object) -> str:
    trimmed = trim_to_codon_frame(seq)
    if not trimmed:
        return ""
    return str(Seq(trimmed).translate(to_stop=False))


def internal_stop_count(seq: object) -> int:
    aa = translate_sequence(seq)
    if len(aa) <= 1:
        return 0
    return aa[:-1].count("*")


def codon_counts(seq: object) -> Counter[str]:
    s = trim_to_codon_frame(seq)
    return Counter(s[i : i + 3] for i in range(0, len(s), 3) if set(s[i : i + 3]) <= set("ACGT"))


def codon_frequencies(seq: object) -> dict[str, float]:
    counts = codon_counts(seq)
    return codon_frequencies_from_counts(counts)


def codon_frequencies_from_counts(counts: Counter[str]) -> dict[str, float]:
    total = sum(counts.values())
    return {codon: counts[codon] / total if total else np.nan for codon in ALL_CODONS}


def translation_qc(seq: object, max_ambiguous_fraction: float = 0.01) -> dict[str, object]:
    s = normalize_rna_to_dna(seq)
    frame_fail = not has_valid_codon_frame(s)
    ambiguous_fail = ambiguous_fraction(s) > max_ambiguous_fraction
    translation_fail = False
    internal_stops: int | None = None
    aa_length: int | None = None
    try:
        aa = translate_sequence(s)
        aa_length = len(aa)
        internal_stops = aa[:-1].count("*") if len(aa) > 1 else 0
    except Exception:
        translation_fail = True
    internal_stop_fail = bool(internal_stops and internal_stops > 0)
    return {
        "frame_fail": frame_fail,
        "ambiguous_fail": ambiguous_fail,
        "internal_stop_fail": internal_stop_fail,
        "translation_fail": translation_fail,
        "internal_stop_count": internal_stops if internal_stops is not None else np.nan,
        "aa_length": aa_length if aa_length is not None else np.nan,
        "codon_total": sum(codon_counts(s).values()),
        "codon_reliable": not (frame_fail or ambiguous_fail or internal_stop_fail or translation_fail),
    }


def rscu(seq: object) -> dict[str, float]:
    counts = codon_counts(seq)
    return rscu_from_counts(counts)


def rscu_from_counts(counts: Counter[str]) -> dict[str, float]:
    values: dict[str, float] = {}
    for aa, codons in AA_TO_CODONS.items():
        total = sum(counts[c] for c in codons)
        expected = total / len(codons) if codons else 0
        for codon in codons:
            values[codon] = counts[codon] / expected if expected else np.nan
    return values


def codon_usage_table(records) -> pd.DataFrame:
    rows = []
    for rec in records:
        counts = codon_counts(str(rec.seq))
        total = sum(counts.values())
        row = {"seq_id": rec.id, "codon_total": total}
        for codon in ALL_CODONS:
            row[f"codon_{codon}"] = counts[codon]
            row[f"codon_freq_{codon}"] = counts[codon] / total if total else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def rscu_table(records) -> pd.DataFrame:
    rows = []
    for rec in records:
        row = {"seq_id": rec.id}
        row.update({f"rscu_{codon}": val for codon, val in rscu(str(rec.seq)).items()})
        rows.append(row)
    return pd.DataFrame(rows)
