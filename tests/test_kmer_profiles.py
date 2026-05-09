from __future__ import annotations

from flugenome3d.kmer_profiles import kmer_counts, kmer_entropy, kmer_frequencies


def test_kmer_counts_sum_possible_windows_for_unambiguous_sequence() -> None:
    seq = "ACGTACGT"
    counts = kmer_counts(seq, 3)
    assert sum(counts.values()) == len(seq) - 3 + 1


def test_kmer_frequencies_sum_to_one() -> None:
    freqs = kmer_frequencies("ACGTACGT", 4)
    assert round(sum(freqs.values()), 12) == 1.0


def test_kmer_entropy_runs() -> None:
    assert kmer_entropy("ACGTACGT", 3) > 0
