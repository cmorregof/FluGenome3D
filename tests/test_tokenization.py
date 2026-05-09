from flugenome3d.tokenization import (
    codon_tokenize,
    codon_tokenize_with_positions,
    crosses_codon_boundary,
    fixed_kmer_tokenize,
    frame_aware_kmer_tokenize,
    non_overlapping_kmer_tokenize,
    overlapping_kmer_tokenize,
    token_contains_cpg,
    token_contains_upa,
    token_crosses_codon_boundary,
)


def test_codon_tokenize():
    toks = codon_tokenize("AAACCCGGG")
    assert toks == ["AAA", "CCC", "GGG"]


def test_codon_tokenize_with_positions():
    toks = codon_tokenize_with_positions("ATGAAA")
    assert toks == [("ATG", 0, 3), ("AAA", 3, 6)]


def test_fixed_kmer_tokenize():
    toks = fixed_kmer_tokenize("ATGAAA", k=3, step=1)
    assert len(toks) == 4


def test_overlapping_kmer_count():
    toks = overlapping_kmer_tokenize("ATGAAA", k=3)
    assert len(toks) == 4


def test_non_overlapping_kmer_count():
    toks = non_overlapping_kmer_tokenize("ATGAAAC", k=3)
    assert toks == ["ATG", "AAA"]


def test_frame_aware_kmer_tokenize():
    toks = frame_aware_kmer_tokenize("AATGAAAC", k=3, frame=1)
    assert toks == ["ATG", "AAA"]


def test_token_cpg_upa_detection():
    assert token_contains_cpg("ACGA") is True
    assert token_contains_cpg("ATTA") is False
    assert token_contains_upa("ATAC") is True
    assert token_contains_upa("ACGC") is False


def test_crosses_codon_boundary():
    assert crosses_codon_boundary(1, 4) is True
    assert crosses_codon_boundary(0, 3) is False


def test_token_crosses_codon_boundary_with_frame():
    assert token_crosses_codon_boundary(1, 4, frame=0) is True
    assert token_crosses_codon_boundary(1, 4, frame=1) is False
