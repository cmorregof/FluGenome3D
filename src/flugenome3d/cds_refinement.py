from __future__ import annotations

from dataclasses import asdict, dataclass
from math import inf
from typing import Any

from Bio.Seq import Seq

from .sequence_metrics import ambiguous_fraction


@dataclass(frozen=True)
class TranslationCandidate:
    sequence: str
    protein_sequence: str
    frame: int
    trim_left: int
    trim_right: int
    length_nt: int
    translation_length: int
    internal_stop_count: int
    ambiguous_aa_count: int
    score: float
    method: str


@dataclass(frozen=True)
class RescueResult:
    status: str
    rescue_method: str
    chosen_frame: int | None
    trim_left: int
    trim_right: int
    length_before: int
    length_after: int
    internal_stops_before: int
    internal_stops_after: int | None
    translation_length: int | None
    ambiguous_fraction: float
    ambiguous_aa_count: int | None
    gaps_removed: bool
    gap_count_removed: int
    notes: str
    refined_sequence: str | None = None

    def public_dict(self) -> dict[str, object]:
        out = asdict(self)
        out.pop("refined_sequence", None)
        return out


def normalize_nucleotide_sequence(seq: object, remove_gaps: bool = True) -> str:
    normalized = "".join(str(seq or "").upper().split()).replace("U", "T")
    if remove_gaps:
        normalized = normalized.replace("-", "")
    return normalized


def _translate(seq: str) -> str:
    if not seq:
        return ""
    return str(Seq(seq).translate(to_stop=False))


def _internal_stop_count(protein_seq: str) -> int:
    if len(protein_seq) <= 1:
        return 0
    return protein_seq[:-1].count("*")


def score_translation_candidate(protein_seq: str) -> float:
    internal_stops = _internal_stop_count(protein_seq)
    ambiguous = protein_seq.count("X")
    terminal_stop_penalty = 0 if (not protein_seq or protein_seq.endswith("*")) else 1
    return internal_stops * 1000 + ambiguous * 100 + terminal_stop_penalty


def _expected_aa_bounds(protein: str, config: dict[str, Any]) -> tuple[int | None, int | None]:
    cfg = config.get("expected_aa_length", {}).get(protein, {})
    min_len = cfg.get("min")
    max_len = cfg.get("max")
    return (int(min_len) if min_len is not None else None, int(max_len) if max_len is not None else None)


def _length_in_bounds(aa_length: int, protein: str, config: dict[str, Any]) -> bool:
    min_len, max_len = _expected_aa_bounds(protein, config)
    if min_len is not None and aa_length < min_len:
        return False
    if max_len is not None and aa_length > max_len:
        return False
    return True


def _candidate_from_trim(seq: str, trim_left: int, trim_right: int, method: str) -> TranslationCandidate | None:
    end = len(seq) - trim_right if trim_right else len(seq)
    candidate = seq[trim_left:end]
    if not candidate or len(candidate) % 3 != 0:
        return None
    if set(candidate) - set("ACGTN"):
        return None
    protein_seq = _translate(candidate)
    return TranslationCandidate(
        sequence=candidate,
        protein_sequence=protein_seq,
        frame=trim_left % 3,
        trim_left=trim_left,
        trim_right=trim_right,
        length_nt=len(candidate),
        translation_length=len(protein_seq),
        internal_stop_count=_internal_stop_count(protein_seq),
        ambiguous_aa_count=protein_seq.count("X"),
        score=score_translation_candidate(protein_seq),
        method=method,
    )


def try_frame_offsets(seq: object) -> list[TranslationCandidate]:
    normalized = normalize_nucleotide_sequence(seq)
    candidates: list[TranslationCandidate] = []
    for frame in range(3):
        remaining = len(normalized) - frame
        if remaining <= 0:
            continue
        trim_right = remaining % 3
        candidate = _candidate_from_trim(normalized, frame, trim_right, f"frame_offset_{frame}")
        if candidate is not None:
            candidates.append(candidate)
    return candidates


def trim_to_valid_cds_candidate(seq: object, max_trim_left: int = 2, max_trim_right: int = 2) -> list[TranslationCandidate]:
    normalized = normalize_nucleotide_sequence(seq)
    candidates: list[TranslationCandidate] = []
    for trim_left in range(max_trim_left + 1):
        for trim_right in range(max_trim_right + 1):
            candidate = _candidate_from_trim(normalized, trim_left, trim_right, f"trim_l{trim_left}_r{trim_right}")
            if candidate is not None:
                candidates.append(candidate)
    return candidates


def _accepts_candidate(candidate: TranslationCandidate, protein: str, config: dict[str, Any]) -> bool:
    if candidate.internal_stop_count > int(config.get("max_internal_stops", 0)):
        return False
    if candidate.ambiguous_aa_count > int(config.get("max_ambiguous_amino_acids", 0)):
        return False
    if not _length_in_bounds(candidate.translation_length, protein, config):
        return False
    return True


def rescue_cds(seq: object, protein: str, subtype: str, config: dict[str, Any]) -> RescueResult:
    del subtype  # subtype-specific rules can be added later without changing the public signature.
    remove_gaps = bool(config.get("remove_gaps", True))
    original = "".join(str(seq or "").upper().split()).replace("U", "T")
    normalized = normalize_nucleotide_sequence(seq, remove_gaps=remove_gaps)
    gap_count = original.count("-") if remove_gaps else 0
    gaps_removed = gap_count > 0
    length_before = len(original)
    max_ambig = float(config.get("max_ambiguous_fraction", 0.01))
    ambig = ambiguous_fraction(normalized)
    before_candidate = _candidate_from_trim(normalized, 0, len(normalized) % 3, "naive_frame0")
    before_stops = before_candidate.internal_stop_count if before_candidate is not None else 0

    if not normalized:
        return RescueResult(
            status="unrescued",
            rescue_method="empty_sequence",
            chosen_frame=None,
            trim_left=0,
            trim_right=0,
            length_before=length_before,
            length_after=0,
            internal_stops_before=before_stops,
            internal_stops_after=None,
            translation_length=None,
            ambiguous_fraction=ambig,
            ambiguous_aa_count=None,
            gaps_removed=gaps_removed,
            gap_count_removed=gap_count,
            notes="No nucleotide sequence after normalization.",
        )

    if ambig > max_ambig:
        return RescueResult(
            status="unrescued",
            rescue_method="ambiguous_fraction_fail",
            chosen_frame=None,
            trim_left=0,
            trim_right=0,
            length_before=length_before,
            length_after=len(normalized),
            internal_stops_before=before_stops,
            internal_stops_after=None,
            translation_length=None,
            ambiguous_fraction=ambig,
            ambiguous_aa_count=None,
            gaps_removed=gaps_removed,
            gap_count_removed=gap_count,
            notes="Ambiguous nucleotide fraction exceeds configured threshold.",
        )

    strict = _candidate_from_trim(normalized, 0, 0, "strict_frame0")
    if strict is not None and _accepts_candidate(strict, protein, config):
        return RescueResult(
            status="strict_pass",
            rescue_method="strict_frame0",
            chosen_frame=0,
            trim_left=0,
            trim_right=0,
            length_before=length_before,
            length_after=strict.length_nt,
            internal_stops_before=before_stops,
            internal_stops_after=strict.internal_stop_count,
            translation_length=strict.translation_length,
            ambiguous_fraction=ambig,
            ambiguous_aa_count=strict.ambiguous_aa_count,
            gaps_removed=gaps_removed,
            gap_count_removed=gap_count,
            notes="Sequence passed strict CDS checks without rescue.",
            refined_sequence=strict.sequence,
        )

    max_trim_left = int(config.get("max_trim_left", 2))
    max_trim_right = int(config.get("max_trim_right", 2))
    candidates = trim_to_valid_cds_candidate(normalized, max_trim_left=max_trim_left, max_trim_right=max_trim_right)
    accepted = [candidate for candidate in candidates if _accepts_candidate(candidate, protein, config)]
    if accepted:
        best = min(accepted, key=lambda c: (c.score, c.trim_left + c.trim_right, c.trim_left, c.trim_right))
        method = "gap_normalized_" + best.method if gaps_removed else best.method
        return RescueResult(
            status="rescued",
            rescue_method=method,
            chosen_frame=best.frame,
            trim_left=best.trim_left,
            trim_right=best.trim_right,
            length_before=length_before,
            length_after=best.length_nt,
            internal_stops_before=before_stops,
            internal_stops_after=best.internal_stop_count,
            translation_length=best.translation_length,
            ambiguous_fraction=ambig,
            ambiguous_aa_count=best.ambiguous_aa_count,
            gaps_removed=gaps_removed,
            gap_count_removed=gap_count,
            notes="Sequence passed conservative trim/frame rescue criteria.",
            refined_sequence=best.sequence,
        )

    best_score = min((candidate.score for candidate in candidates), default=inf)
    return RescueResult(
        status="unrescued",
        rescue_method="no_candidate_passed",
        chosen_frame=None,
        trim_left=0,
        trim_right=0,
        length_before=length_before,
        length_after=len(normalized),
        internal_stops_before=before_stops,
        internal_stops_after=None,
        translation_length=None,
        ambiguous_fraction=ambig,
        ambiguous_aa_count=None,
        gaps_removed=gaps_removed,
        gap_count_removed=gap_count,
        notes=f"No transparent candidate passed rescue criteria; best_score={best_score}.",
    )
