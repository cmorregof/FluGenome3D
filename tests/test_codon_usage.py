from flugenome3d.codon_usage import codon_counts, rscu, trim_to_codon_frame


def test_codon_counts():
    counts = codon_counts("ATGATGTAA")
    assert counts["ATG"] == 2
    assert counts["TAA"] == 1
    assert sum(counts.values()) == len(trim_to_codon_frame("ATGATGTAA")) // 3


def test_rscu_contains_atg():
    vals = rscu("ATGATGTAA")
    assert "ATG" in vals


def test_rscu_excludes_stop_codons():
    vals = rscu("ATGATGTAA")
    assert "TAA" not in vals
    assert "TAG" not in vals
    assert "TGA" not in vals
