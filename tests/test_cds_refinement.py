from __future__ import annotations

from flugenome3d.cds_refinement import normalize_nucleotide_sequence, rescue_cds, try_frame_offsets


def _test_config(max_trim_left: int = 2, max_trim_right: int = 2) -> dict:
    return {
        "remove_gaps": True,
        "max_trim_left": max_trim_left,
        "max_trim_right": max_trim_right,
        "max_ambiguous_fraction": 0.01,
        "max_internal_stops": 0,
        "max_ambiguous_amino_acids": 0,
        "expected_aa_length": {"HA": {"min": 1, "max": 20}, "NA": {"min": 1, "max": 20}},
    }


def test_normalize_nucleotide_sequence_removes_gap_and_converts_u() -> None:
    assert normalize_nucleotide_sequence("acg-u") == "ACGT"


def test_try_frame_offsets_finds_valid_frame_for_synthetic_sequence() -> None:
    candidates = try_frame_offsets("AATGAAACCC")
    assert any(candidate.frame == 1 and candidate.internal_stop_count == 0 for candidate in candidates)


def test_rescue_cds_does_not_accept_internal_stop_without_allowed_trim() -> None:
    result = rescue_cds("ATGTAAATG", protein="HA", subtype="H1N1", config=_test_config(max_trim_left=0, max_trim_right=0))
    assert result.status == "unrescued"
    assert result.internal_stops_before > 0


def test_rescue_cds_reports_method_and_frame() -> None:
    result = rescue_cds("ATGAAA", protein="HA", subtype="H1N1", config=_test_config())
    assert result.status == "strict_pass"
    assert result.rescue_method == "strict_frame0"
    assert result.chosen_frame == 0
