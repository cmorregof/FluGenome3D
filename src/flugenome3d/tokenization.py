from __future__ import annotations

from collections import Counter
from math import log2
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
try:
    from tokenizers import Tokenizer
    from tokenizers.models import BPE
    from tokenizers.pre_tokenizers import Split
    from tokenizers.trainers import BpeTrainer
except ImportError:  # optional dependency for BPE mode
    Tokenizer = None
    BPE = None
    Split = None
    BpeTrainer = None

from .sequence_metrics import normalize_rna_to_dna


def codon_tokenize_with_positions(seq: str) -> list[tuple[str, int, int]]:
    s = normalize_rna_to_dna(seq)
    n = len(s) - (len(s) % 3)
    return [(s[i : i + 3], i, i + 3) for i in range(0, n, 3)]


def codon_tokenize(seq: str) -> list[str]:
    return [token for token, _, _ in codon_tokenize_with_positions(seq)]


def overlapping_kmer_tokenize(seq: str, k: int) -> list[str]:
    return [token for token, _, _ in fixed_kmer_tokenize(seq, k=k, step=1)]


def non_overlapping_kmer_tokenize(seq: str, k: int, offset: int = 0) -> list[str]:
    return [token for token, _, _ in non_overlapping_kmer_tokenize_with_positions(seq, k=k, offset=offset)]


def non_overlapping_kmer_tokenize_with_positions(seq: str, k: int, offset: int = 0) -> list[tuple[str, int, int]]:
    s = normalize_rna_to_dna(seq)
    if k <= 0:
        return []
    start = max(0, int(offset))
    return [(s[i : i + k], i, i + k) for i in range(start, len(s) - k + 1, k)]


def frame_aware_kmer_tokenize(seq: str, k: int, frame: int = 0) -> list[str]:
    return [token for token, _, _ in frame_aware_kmer_tokenize_with_positions(seq, k=k, frame=frame)]


def frame_aware_kmer_tokenize_with_positions(seq: str, k: int, frame: int = 0) -> list[tuple[str, int, int]]:
    return non_overlapping_kmer_tokenize_with_positions(seq, k=k, offset=frame)


def fixed_kmer_tokenize(seq: str, k: int = 6, step: int | None = None) -> list[tuple[str, int, int]]:
    s = normalize_rna_to_dna(seq)
    if step is None:
        step = k
    return [(s[i : i + k], i, i + k) for i in range(0, len(s) - k + 1, step)]


def token_entropy(tokens: list[str]) -> float:
    counts = Counter(tokens)
    total = sum(counts.values())
    if not total:
        return np.nan
    return -sum((v / total) * log2(v / total) for v in counts.values())


def token_contains_cpg(token: str) -> bool:
    return "CG" in normalize_rna_to_dna(token)


def token_contains_upa(token: str) -> bool:
    return "TA" in normalize_rna_to_dna(token)


def token_crosses_codon_boundary(start: int, end: int, frame: int = 0) -> bool:
    # A token crosses an internal codon boundary if it spans positions with different floor(pos/3)
    return any((pos - frame) % 3 == 0 for pos in range(start + 1, end))


def crosses_codon_boundary(start: int, end: int) -> bool:
    return token_crosses_codon_boundary(start, end, frame=0)


def token_metrics(seq_id: str, seq: str, method: str, tokens_with_pos: list[tuple[str, int, int]]) -> dict[str, object]:
    tokens = [t for t, _, _ in tokens_with_pos]
    lengths = [len(t) for t in tokens]
    n_cross = sum(crosses_codon_boundary(i, j) for _, i, j in tokens_with_pos)
    n_cpg = sum("CG" in t for t in tokens)
    n_upa = sum(("TA" in t) or ("UA" in t) for t in tokens)
    return {
        "seq_id": seq_id,
        "method": method,
        "sequence_length": len(seq),
        "n_tokens": len(tokens),
        "tokens_per_kb": len(tokens) / (len(seq) / 1000) if seq else np.nan,
        "mean_token_length": float(np.mean(lengths)) if lengths else np.nan,
        "median_token_length": float(np.median(lengths)) if lengths else np.nan,
        "vocab_size": len(set(tokens)),
        "vocab_entropy": token_entropy(tokens),
        "cross_codon_boundary_fraction": n_cross / len(tokens) if tokens else np.nan,
        "cpg_token_fraction": n_cpg / len(tokens) if tokens else np.nan,
        "upa_token_fraction": n_upa / len(tokens) if tokens else np.nan,
    }


def audit_records(records: list[SeqRecord], k_values: tuple[int, ...] = (3, 4, 6)) -> pd.DataFrame:
    rows = []
    for rec in records:
        seq = str(rec.seq)
        rows.append(token_metrics(rec.id, seq, "codon", codon_tokenize_with_positions(seq)))
        for k in k_values:
            rows.append(token_metrics(rec.id, seq, f"kmer{k}_nonoverlap", fixed_kmer_tokenize(seq, k=k, step=k)))
            rows.append(token_metrics(rec.id, seq, f"kmer{k}_overlap", fixed_kmer_tokenize(seq, k=k, step=1)))
    return pd.DataFrame(rows)


def tokenize_dataset(panel: pd.DataFrame, tokenizer_config: dict[str, object]) -> pd.DataFrame:
    rows = []
    tokenizer = str(tokenizer_config.get("tokenizer"))
    k = int(tokenizer_config.get("k", 3))
    offset = int(tokenizer_config.get("offset", 0))
    seq_col = str(tokenizer_config.get("sequence_column", "sequence"))
    id_col = str(tokenizer_config.get("id_column", "internal_sequence_id"))
    for row in panel.itertuples(index=False):
        seq = getattr(row, seq_col)
        seq_id = getattr(row, id_col)
        if tokenizer == "codon":
            tokens = codon_tokenize(seq)
        elif tokenizer == "overlapping_kmer":
            tokens = overlapping_kmer_tokenize(seq, k)
        elif tokenizer == "non_overlapping_kmer":
            tokens = non_overlapping_kmer_tokenize(seq, k, offset=offset)
        elif tokenizer == "frame_aware_kmer":
            tokens = frame_aware_kmer_tokenize(seq, k, frame=offset)
        else:
            raise ValueError(f"Unknown tokenizer: {tokenizer}")
        rows.append({"internal_sequence_id": seq_id, "tokenizer": tokenizer, "k": k, "offset": offset, "n_tokens": len(tokens)})
    return pd.DataFrame(rows)


def train_basic_bpe_from_records(records: list[SeqRecord], vocab_size: int = 600, out_path: str | Path | None = None):
    if Tokenizer is None:
        raise ImportError("The optional dependency 'tokenizers' is required for BPE mode. Install with `pip install tokenizers`.")
    tokenizer = Tokenizer(BPE(unk_token="[UNK]"))
    tokenizer.pre_tokenizer = Split(pattern="", behavior="isolated")
    trainer = BpeTrainer(vocab_size=vocab_size, special_tokens=["[UNK]"])
    seqs = [normalize_rna_to_dna(str(rec.seq)) for rec in records]
    tokenizer.train_from_iterator(seqs, trainer=trainer)
    if out_path is not None:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        tokenizer.save(str(out_path))
    return tokenizer
