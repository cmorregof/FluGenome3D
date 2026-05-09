# FluGenome3D Phase 5 tokenization audit report

This phase audits how deterministic tokenization choices segment Influenza A HA/NA sequences. It is a descriptive tokenization audit before learned BPE or GROVER tokenizers.

## Representation vs. tokenization

Phase 4 treated sequences as vector representations. Phase 5 treats tokenization itself as the object of study: vocabulary size, token entropy, CpG/UpA-containing token fractions, codon-boundary crossing and Jensen-Shannon distances between HA/NA and subtype groups.

## Datasets used

- Raw nucleotide tokenizers use `mvp_panel`: 20,000 HA/NA sequences.
- CDS-aware tokenizers use `mvp_cds_refined_panel`: 8,974 refined CDS sequences.
- Codon and frame-aware claims are restricted to the CDS-refined panel.

## Tokenizers audited

| tokenizer | dataset | source | family | mode | k | offset | n_sequences | observed_vocab_size | total_tokens |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cds_codon | mvp_cds_refined_panel | cds_refined | codon | codon | 3 | 0 | 8974 | 68 | 4639220 |
| cds_frame_k3 | mvp_cds_refined_panel | cds_refined | frame_aware_kmer | frame_aware | 3 | 0 | 8974 | 68 | 4639220 |
| cds_frame_k6 | mvp_cds_refined_panel | cds_refined | frame_aware_kmer | frame_aware | 6 | 0 | 8974 | 3058 | 2317384 |
| cds_nonoverlap_k3_offset0 | mvp_cds_refined_panel | cds_refined | non_overlapping_kmer | non_overlapping | 3 | 0 | 8974 | 68 | 4639220 |
| cds_nonoverlap_k3_offset1 | mvp_cds_refined_panel | cds_refined | non_overlapping_kmer | non_overlapping | 3 | 1 | 8974 | 68 | 4630246 |
| cds_nonoverlap_k3_offset2 | mvp_cds_refined_panel | cds_refined | non_overlapping_kmer | non_overlapping | 3 | 2 | 8974 | 68 | 4630246 |
| raw_nonoverlap_k3 | mvp_panel | raw_nucleotide | non_overlapping_kmer | non_overlapping | 3 | 0 | 20000 | 114 | 10496741 |
| raw_nonoverlap_k4 | mvp_panel | raw_nucleotide | non_overlapping_kmer | non_overlapping | 4 | 0 | 20000 | 355 | 7868050 |
| raw_nonoverlap_k5 | mvp_panel | raw_nucleotide | non_overlapping_kmer | non_overlapping | 5 | 0 | 20000 | 1141 | 6294280 |
| raw_nonoverlap_k6 | mvp_panel | raw_nucleotide | non_overlapping_kmer | non_overlapping | 6 | 0 | 20000 | 4113 | 5243363 |
| raw_overlap_k3 | mvp_panel | raw_nucleotide | overlapping_kmer | overlapping | 3 | 0 | 20000 | 121 | 31459751 |
| raw_overlap_k4 | mvp_panel | raw_nucleotide | overlapping_kmer | overlapping | 4 | 0 | 20000 | 458 | 31439751 |

_Table truncated to 12 of 14 rows._

## Entropy and effective vocabulary

Token entropy and effective vocabulary vary by tokenizer family, k, protein and subtype. Larger k generally increases the observed vocabulary available to the audit, while non-overlapping tokenizers produce fewer tokens per sequence.

| tokenizer | mean_group_entropy_bits |
| --- | --- |
| raw_overlap_k6 | 10.0835 |
| raw_overlap_k5 | 9.16638 |
| raw_nonoverlap_k6 | 7.90346 |
| cds_frame_k6 | 7.89993 |
| raw_nonoverlap_k5 | 7.87911 |
| raw_overlap_k4 | 7.62741 |
| raw_nonoverlap_k4 | 7.25876 |
| raw_overlap_k3 | 5.79937 |
| cds_nonoverlap_k3_offset1 | 5.67102 |
| cds_nonoverlap_k3_offset2 | 5.65158 |
| raw_nonoverlap_k3 | 5.63631 |
| cds_codon | 5.61835 |

_Table truncated to 12 of 14 rows._

| tokenizer | protein_subtype | observed_vocab_size | group_effective_vocab_size |
| --- | --- | --- | --- |
| cds_codon | HA-H1N1 | 66 | 48.1941 |
| cds_codon | HA-H3N2 | 64 | 50.199 |
| cds_codon | NA-H1N1 | 65 | 49.0558 |
| cds_codon | NA-H3N2 | 63 | 50.8799 |
| cds_frame_k3 | HA-H1N1 | 66 | 48.1941 |
| cds_frame_k3 | HA-H3N2 | 64 | 50.199 |
| cds_frame_k3 | NA-H1N1 | 65 | 49.0558 |
| cds_frame_k3 | NA-H3N2 | 63 | 50.8799 |
| cds_frame_k6 | HA-H1N1 | 1886 | 431.37 |
| cds_frame_k6 | HA-H3N2 | 1702 | 375.52 |
| cds_frame_k6 | NA-H1N1 | 1772 | 373.161 |
| cds_frame_k6 | NA-H3N2 | 1442 | 289.269 |
| cds_nonoverlap_k3_offset0 | HA-H1N1 | 66 | 48.1941 |
| cds_nonoverlap_k3_offset0 | HA-H3N2 | 64 | 50.199 |
| cds_nonoverlap_k3_offset0 | NA-H1N1 | 65 | 49.0558 |
| cds_nonoverlap_k3_offset0 | NA-H3N2 | 63 | 50.8799 |

## CpG and UpA/TA token analysis

CpG-containing and UpA/TA-containing token fractions are summarized descriptively by tokenizer and group. Because records are stored as DNA alphabet, TA is used as the DNA proxy for UpA.

| tokenizer | protein_subtype | mean_cpg_token_fraction | mean_upa_token_fraction |
| --- | --- | --- | --- |
| cds_codon | HA-H1N1 | 0.0178382 | 0.127308 |
| cds_codon | HA-H3N2 | 0.0275682 | 0.0941805 |
| cds_codon | NA-H1N1 | 0.0201546 | 0.109489 |
| cds_codon | NA-H3N2 | 0.0248026 | 0.0940546 |
| cds_frame_k3 | HA-H1N1 | 0.0178382 | 0.127308 |
| cds_frame_k3 | HA-H3N2 | 0.0275682 | 0.0941805 |
| cds_frame_k3 | NA-H1N1 | 0.0201546 | 0.109489 |
| cds_frame_k3 | NA-H3N2 | 0.0248026 | 0.0940546 |
| cds_frame_k6 | HA-H1N1 | 0.0653312 | 0.288692 |
| cds_frame_k6 | HA-H3N2 | 0.0877964 | 0.239721 |
| cds_frame_k6 | NA-H1N1 | 0.0739212 | 0.24095 |
| cds_frame_k6 | NA-H3N2 | 0.0735961 | 0.223695 |
| cds_nonoverlap_k3_offset0 | HA-H1N1 | 0.0178382 | 0.127308 |
| cds_nonoverlap_k3_offset0 | HA-H3N2 | 0.0275682 | 0.0941805 |
| cds_nonoverlap_k3_offset0 | NA-H1N1 | 0.0201546 | 0.109489 |
| cds_nonoverlap_k3_offset0 | NA-H3N2 | 0.0248026 | 0.0940546 |

## Codon-boundary crossing

Codon-boundary crossing is interpreted as CDS-aware only for tokenizers built on `mvp_cds_refined_panel`. For raw nucleotide tokenizers, boundary fractions are reported as position-0 segmentation diagnostics and should not be treated as coding-frame evidence.

| tokenizer | boundary_context | mean_codon_boundary_crossing_fraction |
| --- | --- | --- |
| cds_codon | cds_refined_codon_aligned | 0 |
| cds_frame_k3 | cds_refined_frame0 | 0 |
| cds_frame_k6 | cds_refined_frame0 | 1 |
| cds_nonoverlap_k3_offset0 | cds_refined_offset_sensitivity | 0 |
| cds_nonoverlap_k3_offset1 | cds_refined_offset_sensitivity | 1 |
| cds_nonoverlap_k3_offset2 | cds_refined_offset_sensitivity | 1 |
| raw_nonoverlap_k3 | position0_proxy_not_cds | 0 |
| raw_nonoverlap_k4 | position0_proxy_not_cds | 1 |
| raw_nonoverlap_k5 | position0_proxy_not_cds | 1 |
| raw_nonoverlap_k6 | position0_proxy_not_cds | 1 |
| raw_overlap_k3 | position0_proxy_not_cds | 0.666342 |
| raw_overlap_k4 | position0_proxy_not_cds | 1 |
| raw_overlap_k5 | position0_proxy_not_cds | 1 |
| raw_overlap_k6 | position0_proxy_not_cds | 1 |

## Jensen-Shannon distances

Jensen-Shannon distances compare group-level token distributions across HA-H1N1, NA-H1N1, HA-H3N2 and NA-H3N2. Higher values indicate stronger descriptive separation of aggregate token distributions, not predictive performance.

| tokenizer | mean_pairwise_js_distance |
| --- | --- |
| cds_frame_k6 | 0.894248 |
| raw_nonoverlap_k6 | 0.78688 |
| raw_overlap_k6 | 0.684813 |
| raw_nonoverlap_k5 | 0.549605 |
| raw_overlap_k5 | 0.423933 |
| raw_nonoverlap_k4 | 0.320845 |
| raw_overlap_k4 | 0.224535 |
| cds_nonoverlap_k3_offset2 | 0.219128 |
| cds_codon | 0.212845 |
| cds_frame_k3 | 0.212845 |
| cds_nonoverlap_k3_offset0 | 0.212845 |
| cds_nonoverlap_k3_offset1 | 0.210403 |

_Table truncated to 12 of 14 rows._

## Top enriched tokens

Top token tables contain only short public-safe tokens: k-mers of length <= 6 or codons. These tokens are descriptive enrichment summaries, not antigenic, escape, fitness or vaccine markers.

| tokenizer | protein_subtype | rank | token | group_frequency | background_frequency | enrichment_ratio | contains_cpg | contains_upa |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cds_codon | HA-H1N1 | 1 | CAC | 0.013196 | 0.00588458 | 2.24247 | False | False |
| cds_codon | HA-H1N1 | 2 | CTA | 0.0215403 | 0.0101275 | 2.12692 | False | True |
| cds_codon | HA-H1N1 | 3 | GTA | 0.0281845 | 0.0133761 | 2.10708 | False | True |
| cds_codon | HA-H1N1 | 4 | GCA | 0.028881 | 0.0137618 | 2.09863 | False | False |
| cds_codon | HA-H1N1 | 5 | CTG | 0.0213765 | 0.0103835 | 2.0587 | False | False |
| cds_codon | HA-H1N1 | 6 | CTC | 0.00779514 | 0.0042994 | 1.81308 | False | False |
| cds_codon | HA-H1N1 | 7 | CCG | 0.00779364 | 0.00439371 | 1.77382 | True | False |
| cds_codon | HA-H1N1 | 8 | GCC | 0.0155798 | 0.00880194 | 1.77004 | False | False |
| cds_codon | HA-H1N1 | 9 | ACA | 0.0426171 | 0.024785 | 1.71947 | False | False |
| cds_codon | HA-H1N1 | 10 | TAC | 0.0230541 | 0.0136914 | 1.68384 | False | True |
| cds_codon | HA-H3N2 | 1 | TGA | 0.00173843 | 4.27135e-06 | 406.903 | False | False |
| cds_codon | HA-H3N2 | 2 | CTT | 0.0167209 | 0.00423035 | 3.9526 | False | False |
| cds_codon | HA-H3N2 | 3 | CGA | 0.0063231 | 0.00265507 | 2.38152 | True | False |
| cds_codon | HA-H3N2 | 4 | ACG | 0.00825755 | 0.00354124 | 2.33183 | True | False |
| cds_codon | HA-H3N2 | 5 | CAA | 0.0321619 | 0.0158467 | 2.02956 | False | False |
| cds_codon | HA-H3N2 | 6 | ATC | 0.0309867 | 0.0162744 | 1.90401 | False | False |

## Figures

- `results/figures/fig18_token_entropy_by_tokenizer.png`
- `results/figures/fig19_effective_vocab_by_group.png`
- `results/figures/fig20_cpg_upa_token_fraction.png`
- `results/figures/fig21_codon_boundary_crossing.png`
- `results/figures/fig22_token_js_distance_heatmap.png`
- `results/figures/fig23_top_tokens_by_group.png`

## What this phase can support

- We audit how deterministic tokenization choices segment HA/NA sequences.
- Token entropy and effective vocabulary vary by protein/subtype group.
- CpG/UpA-containing token fractions are summarized descriptively.
- Codon-boundary crossing is evaluated as CDS-aware only in CDS-refined sequences.
- This phase establishes transparent tokenization baselines before learned BPE/GROVER tokenizers.

## What this phase does not support

- It does not identify antigenic sites.
- It does not predict escape, viral evolution, vaccine relevance or fitness.
- It does not replace biological language models.
- CpG/UpA token fractions should not be interpreted as pathogenicity or antigenicity evidence.

## Limitations

- Raw nucleotide tokenizers do not require a validated coding frame and therefore cannot support codon-frame claims.
- CDS-aware tokenizers inherit the Phase 3 rescue criteria and should be interpreted together with CDS refinement QC.
- Token enrichment is sensitive to k, overlap mode and group composition.
- Figures summarize aggregate token distributions; they do not inspect individual biological variants.

## Recommendation for Phase 6

Use this tokenization baseline to decide whether Phase 6 should add a learned tokenizer audit. If approved, BPE should be introduced as an explicitly optional learned baseline, still separate from GROVER and still without prediction, antigenicity or vaccine claims.
