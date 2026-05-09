# FluGenome3D Phase 3 CDS refinement report

This phase was necessary because Phase 2 found that naive codon-frame assumptions were not reliable for all HA/NA nucleotide records. Non-coding sequence-context metrics remain valid on raw nucleotide context, while codon usage and RSCU require a separate CDS-aware subset.

## Representations

- `raw_nucleotide_context`: original normalized nucleotide strings for GC, CpG/UpA, dinucleotide odds and k-mer entropy. No coding frame is required.
- `cds_strict`: sequences that pass conservative frame, ambiguity, translation, internal-stop and protein-length checks without trimming.
- `cds_rescued`: sequences that fail naive QC but pass transparent 0-2 nt trim/frame rescue rules.

## Why naive QC failed

Likely contributors include UTR or non-coding flanks, partial records, gaps or ambiguous characters, nonzero CDS frame offsets, aligned-vs-CDS nucleotide representations, non-multiple-of-3 lengths, and records not directly translatable from position 0. This phase audits those possibilities without assuming a single cause.

| subtype | protein | n_sequences | n_naive_qc_pass | n_any_qc_fail | n_frame_fail | n_internal_stop_fail | n_ambiguous_fail | n_translation_fail | mean_length_nt | median_length_nt | mean_internal_stop_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1N1 | HA | 5000 | 2354 | 2646 | 1154 | 2645 | 0 | 0 | 1726.44 | 1724 | 16.1648 |
| H1N1 | NA | 5000 | 2394 | 2606 | 1783 | 2499 | 0 | 0 | 1423.27 | 1418 | 14.3674 |
| H3N2 | HA | 5000 | 1922 | 3078 | 1540 | 3064 | 0 | 0 | 1721.38 | 1721 | 19.1378 |
| H3N2 | NA | 5000 | 2053 | 2947 | 1847 | 2936 | 0 | 0 | 1428.86 | 1431 | 14.0322 |

## Length modulo 3

| subtype | protein | length_mod3 | n_sequences |
| --- | --- | --- | --- |
| H1N1 | HA | 0 | 3846 |
| H1N1 | HA | 1 | 640 |
| H1N1 | HA | 2 | 514 |
| H1N1 | NA | 0 | 3217 |
| H1N1 | NA | 1 | 437 |
| H1N1 | NA | 2 | 1346 |
| H3N2 | HA | 0 | 3460 |
| H3N2 | HA | 1 | 861 |
| H3N2 | HA | 2 | 679 |
| H3N2 | NA | 0 | 3153 |
| H3N2 | NA | 1 | 1182 |
| H3N2 | NA | 2 | 665 |

## Rescue rules

- Normalize to uppercase DNA alphabet and convert U to T.
- Remove gaps only when configured, while recording `gaps_removed` and `gap_count_removed`.
- Try deterministic trims of 0-2 nt at the left and 0-2 nt at the right.
- Accept only candidates with no internal stops, ambiguity under threshold, no ambiguous amino acids, and expected HA/NA protein-length range.
- Leave all other sequences as `unrescued`.
- No reverse-complement rescue, alignment rescue, reference-guided repair, prediction, optimization or biological validation is performed.

## Rescue outcome

- Strict pass sequences: 8675
- Rescued sequences: 299
- Unrescued sequences: 11026
- Final refined CDS panel size: 8974
- Codon-reliable sequences after refinement: 8974

| subtype | protein | status | n_sequences |
| --- | --- | --- | --- |
| H1N1 | HA | rescued | 3 |
| H1N1 | HA | strict_pass | 2345 |
| H1N1 | HA | unrescued | 2652 |
| H1N1 | NA | rescued | 186 |
| H1N1 | NA | strict_pass | 2391 |
| H1N1 | NA | unrescued | 2423 |
| H3N2 | HA | rescued | 74 |
| H3N2 | HA | strict_pass | 1914 |
| H3N2 | HA | unrescued | 3012 |
| H3N2 | NA | rescued | 36 |
| H3N2 | NA | strict_pass | 2025 |
| H3N2 | NA | unrescued | 2939 |

## Refined codon/RSCU QC

Codon usage and RSCU are now reported on the strict/rescued CDS subset with explicit rescue status and QC flags. These are more reliable than Phase 2 naive codon summaries, but they are still not biological validation of CDS boundaries.

| subtype | protein | n_refined_sequences | n_codon_reliable_after_refinement | reliable_fraction_after_refinement | n_frame_fail_after_refinement | n_internal_stop_fail_after_refinement | n_translation_fail_after_refinement | n_strict_pass | n_rescued | n_unrescued |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1N1 | HA | 2348 | 2348 | 1 | 0 | 0 | 0 | 2345 | 3 | 2652 |
| H1N1 | NA | 2577 | 2577 | 1 | 0 | 0 | 0 | 2391 | 186 | 2423 |
| H3N2 | HA | 1988 | 1988 | 1 | 0 | 0 | 0 | 1914 | 74 | 3012 |
| H3N2 | NA | 2061 | 2061 | 1 | 0 | 0 | 0 | 2025 | 36 | 2939 |

## Outputs

Local gitignored panels:

- `data/processed/panels/mvp_cds_strict_panel.parquet`
- `data/processed/panels/mvp_cds_rescued_panel.parquet`
- `data/processed/panels/mvp_cds_refined_panel.parquet`
- `data/processed/metrics/mvp_cds_refined_codon_metrics.parquet`

Public aggregate tables:

- `results/tables/phase3_cds_qc_failure_breakdown.csv`
- `results/tables/phase3_length_mod3_by_group.csv`
- `results/tables/phase3_internal_stop_distribution.csv`
- `results/tables/phase3_rescue_status_by_group.csv`
- `results/tables/phase3_refined_codon_usage_summary.csv`
- `results/tables/phase3_refined_rscu_summary.csv`
- `results/tables/phase3_refined_translation_qc_summary.csv`

Figures:

- `results/figures/fig8_cds_qc_failure_breakdown.png`
- `results/figures/fig9_length_mod3_by_group.png`
- `results/figures/fig10_refined_rscu_heatmap.png`
- `results/figures/fig11_rescue_status_by_group.png`

## Permitted claims

- We identified limitations of naive codon-frame assumptions in HA/NA sequence records.
- We separated non-coding sequence-context metrics from CDS-dependent codon analyses.
- Codon usage and RSCU are reported only on strict/rescued CDS subsets with explicit QC.
- This phase improves reliability of downstream tokenization and codon-level analysis.

## Prohibited claims

- Rescued CDS are biologically validated.
- Codon usage explains antigenic drift.
- RSCU predicts fitness.
- Frame rescue identifies functional variants.
- Translation QC implies antigenic relevance.

## Recommendation for Phase 4

Use `raw_nucleotide_context` for tokenizer-independent composition and k-mer analyses, and use `mvp_cds_refined_panel` only for CDS-dependent codon/RSCU summaries. The next phase should focus on tokenizer-ready representations over the MVP/refined panels without GROVER/BPE unless explicitly approved.
