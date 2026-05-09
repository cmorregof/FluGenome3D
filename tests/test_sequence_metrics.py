import numpy as np

from flugenome3d.sequence_metrics import (
    ambiguous_fraction,
    cpg_observed_expected,
    dinucleotide_counts,
    dinucleotide_odds_ratio,
    gc_content,
    gc_fraction,
    kmer_counts,
    upa_observed_expected,
)


def test_gc_fraction():
    assert gc_fraction("ACGT") == 0.5
    assert gc_content("ACGT") == 0.5


def test_ambiguous_fraction():
    assert ambiguous_fraction("ACGTNN") == 2 / 6


def test_kmer_counts():
    counts = kmer_counts("ACGTAC", 3)
    assert counts["ACG"] == 1
    assert counts["CGT"] == 1


def test_dinucleotide_odds_ratio_runs():
    val = dinucleotide_odds_ratio("ACGTCGTA", "CG")
    assert val >= 0


def test_cpg_upa_observed_expected_zero_denominator_does_not_explode():
    assert np.isnan(cpg_observed_expected("AAAA"))
    assert np.isnan(upa_observed_expected("CCCC"))


def test_dinucleotide_counts_sum_adjacent_pairs():
    counts = dinucleotide_counts("ACGTNN")
    assert sum(counts.values()) == 5
