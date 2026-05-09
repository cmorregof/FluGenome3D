from __future__ import annotations

import numpy as np
import pandas as pd

from flugenome3d.codon_usage import ALL_CODONS, CODON_TO_AA
from flugenome3d.representations import (
    build_codon_frequency_matrix,
    build_kmer_count_matrix,
    build_kmer_tfidf_matrix,
    build_rscu_matrix,
)


def test_kmer_count_matrix_dimensions_and_feature_names() -> None:
    matrix, names = build_kmer_count_matrix(["ACGT"], k=2)
    assert matrix.shape == (1, 16)
    assert all(len(name) == 2 for name in names)
    assert matrix.sum() == 3


def test_kmer_tfidf_has_no_nan() -> None:
    matrix, _ = build_kmer_tfidf_matrix(["ACGTACGT", "AAAAACCC"], k=3)
    assert not np.isnan(matrix.data).any()


def test_codon_frequency_matrix_has_64_columns() -> None:
    row = {f"codon_freq_{codon}": 0.0 for codon in ALL_CODONS}
    row["codon_freq_ATG"] = 1.0
    matrix, names = build_codon_frequency_matrix(pd.DataFrame([row]))
    assert matrix.shape == (1, 64)
    assert len(names) == 64


def test_rscu_matrix_excludes_stops() -> None:
    sense_codons = [codon for codon in ALL_CODONS if CODON_TO_AA[codon] != "*"]
    row = {f"rscu_{codon}": 1.0 for codon in sense_codons}
    matrix, names = build_rscu_matrix(pd.DataFrame([row]))
    assert matrix.shape == (1, len(sense_codons))
    assert "TAA" not in names
    assert "TAG" not in names
    assert "TGA" not in names
