from __future__ import annotations

import pandas as pd

from flugenome3d.structure_mapping import align_consensus_to_pdb, parse_rcsb_fasta, residue_signal_summary


def test_parse_rcsb_fasta() -> None:
    text = ">1ABC_1|Chain A|Protein name|Organism\nACDE\nFG\n"
    entries = parse_rcsb_fasta(text)
    assert entries[0]["entity"] == "1ABC_1"
    assert entries[0]["chains"] == "Chain A"
    assert entries[0]["sequence"] == "ACDEFG"


def test_align_consensus_to_pdb_reports_identity() -> None:
    chunks, metrics = align_consensus_to_pdb("MKTAAACCC", "TAAA")
    assert chunks
    assert metrics["mapped_residues"] >= 4
    assert metrics["identity"] == 1


def test_residue_signal_summary_no_sequence_columns() -> None:
    panel = pd.DataFrame(
        {
            "subtype": ["H1N1", "H1N1"],
            "protein": ["HA", "HA"],
            "refined_sequence": ["ATGAAACCC", "ATGAAGCCA"],
        }
    )
    summary = residue_signal_summary(panel)
    assert not summary.empty
    assert "refined_sequence" not in summary.columns
    assert {"gc_fraction_codon", "cpg_codon_fraction", "upa_codon_fraction", "aa_entropy"}.issubset(summary.columns)
