# FluGenome3D Phase 6 tokenization stability report

This phase evaluates the stability of deterministic tokenization metrics under stratified bootstrap and temporal grouping. It is a descriptive robustness audit before learned BPE or GROVER tokenizers.

## Why robustness before learned tokenizers

Phase 5 compared deterministic tokenization choices. Phase 6 checks whether those descriptive differences are stable under resampling, whether group distances are consistent through time windows, and whether top-token enrichment is robust to sampling variation. This avoids treating a single full-panel estimate as if it were automatically stable.

## Tokenizers evaluated

raw_overlap_k3, raw_overlap_k6, raw_nonoverlap_k6, cds_codon, cds_frame_k6, cds_nonoverlap_k3_offset0, cds_nonoverlap_k3_offset1, cds_nonoverlap_k3_offset2

## Bootstrap design

- Bootstrap replicates: 100
- Sequence metric bootstrap: stratified by HA-H1N1, NA-H1N1, HA-H3N2 and NA-H3N2.
- Maximum sampled rows per group per bootstrap: 1000
- JS and top-token bootstrap: aggregate group token distributions were resampled within each tokenizer/group because Phase 5 does not store ordered token lists per sequence.

## JS distance stability

Mean pairwise Jensen-Shannon distance is summarized with bootstrap confidence intervals. Higher values indicate stronger descriptive separation of aggregate token distributions, not predictive performance.

| tokenizer | mean | std | ci_lower | ci_upper | ci_width |
| --- | --- | --- | --- | --- | --- |
| cds_frame_k6 | 0.895221 | 0.000241411 | 0.894746 | 0.895671 | 0.000924654 |
| raw_nonoverlap_k6 | 0.789584 | 0.000404968 | 0.788726 | 0.790331 | 0.00160529 |
| raw_overlap_k6 | 0.688757 | 0.000417339 | 0.688027 | 0.689468 | 0.00144075 |
| cds_nonoverlap_k3_offset2 | 0.219342 | 0.000719688 | 0.218075 | 0.220684 | 0.00260925 |
| cds_codon | 0.213097 | 0.000789478 | 0.211668 | 0.214611 | 0.00294291 |
| cds_nonoverlap_k3_offset0 | 0.21305 | 0.000765663 | 0.211656 | 0.214433 | 0.0027771 |
| cds_nonoverlap_k3_offset1 | 0.210726 | 0.000855567 | 0.209118 | 0.212405 | 0.00328732 |
| raw_overlap_k3 | 0.111653 | 0.000770245 | 0.110257 | 0.113245 | 0.00298811 |

## Top-token stability

Top-token stability is measured as Jaccard overlap between bootstrap top-20 enriched tokens and the full-panel top-20 reference for the same tokenizer/group.

| tokenizer | mean_top_token_jaccard |
| --- | --- |
| cds_nonoverlap_k3_offset1 | 0.932619 |
| raw_overlap_k3 | 0.927398 |
| cds_codon | 0.906035 |
| cds_nonoverlap_k3_offset2 | 0.903914 |
| cds_nonoverlap_k3_offset0 | 0.903651 |
| cds_frame_k6 | 0.366073 |
| raw_nonoverlap_k6 | 0.23725 |
| raw_overlap_k6 | 0.203458 |

## Temporal stability

Temporal bins used: 2009-2014, 2015-2019, 2020+, pre-2009

Temporal metrics are reported only for tokenizer/window combinations where all four protein-subtype groups pass the minimum group count threshold.

| row_type | tokenizer | time_window | protein_subtype | comparison | n_sequences | mean_token_entropy_bits | mean_effective_vocab_size | mean_cpg_token_fraction | mean_upa_token_fraction | js_distance |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| group_metric | cds_codon | 2009-2014 | HA-H1N1 |  | 985 | 5.57565 | 47.6918 | 0.0191748 | 0.133043 |  |
| group_metric | cds_codon | 2009-2014 | HA-H3N2 |  | 676 | 5.63381 | 49.654 | 0.0273707 | 0.0933161 |  |
| group_metric | cds_codon | 2009-2014 | NA-H1N1 |  | 1053 | 5.60622 | 48.7132 | 0.0206244 | 0.110882 |  |
| group_metric | cds_codon | 2009-2014 | NA-H3N2 |  | 791 | 5.6652 | 50.7459 | 0.0249278 | 0.0941946 |  |
| group_metric | cds_codon | 2015-2019 | HA-H1N1 |  | 749 | 5.56623 | 47.3815 | 0.01819 | 0.124097 |  |
| group_metric | cds_codon | 2015-2019 | HA-H3N2 |  | 507 | 5.6401 | 49.8703 | 0.0265014 | 0.0984158 |  |
| group_metric | cds_codon | 2015-2019 | NA-H1N1 |  | 809 | 5.60292 | 48.6016 | 0.0209999 | 0.110334 |  |
| group_metric | cds_codon | 2015-2019 | NA-H3N2 |  | 597 | 5.65624 | 50.4316 | 0.0254528 | 0.0914191 |  |
| group_metric | cds_codon | 2020+ | HA-H1N1 |  | 379 | 5.57359 | 47.6237 | 0.0183064 | 0.120578 |  |
| group_metric | cds_codon | 2020+ | HA-H3N2 |  | 279 | 5.63788 | 49.7941 | 0.0285662 | 0.0941945 |  |
| group_metric | cds_codon | 2020+ | NA-H1N1 |  | 406 | 5.5959 | 48.3658 | 0.02034 | 0.105531 |  |
| group_metric | cds_codon | 2020+ | NA-H3N2 |  | 288 | 5.65632 | 50.4341 | 0.0254791 | 0.0906134 |  |

Temporal mean pairwise JS distances:

| tokenizer | time_window | comparison | js_distance |
| --- | --- | --- | --- |
| cds_codon | 2009-2014 | mean_pairwise_nonidentical | 0.223673 |
| cds_codon | 2015-2019 | mean_pairwise_nonidentical | 0.223634 |
| cds_codon | 2020+ | mean_pairwise_nonidentical | 0.218789 |
| cds_codon | pre-2009 | mean_pairwise_nonidentical | 0.213403 |
| cds_frame_k6 | 2009-2014 | mean_pairwise_nonidentical | 0.917409 |
| cds_frame_k6 | 2015-2019 | mean_pairwise_nonidentical | 0.919145 |
| cds_frame_k6 | 2020+ | mean_pairwise_nonidentical | 0.917334 |
| cds_frame_k6 | pre-2009 | mean_pairwise_nonidentical | 0.911726 |
| cds_nonoverlap_k3_offset0 | 2009-2014 | mean_pairwise_nonidentical | 0.223673 |
| cds_nonoverlap_k3_offset0 | 2015-2019 | mean_pairwise_nonidentical | 0.223634 |
| cds_nonoverlap_k3_offset0 | 2020+ | mean_pairwise_nonidentical | 0.218789 |
| cds_nonoverlap_k3_offset0 | pre-2009 | mean_pairwise_nonidentical | 0.213403 |
| cds_nonoverlap_k3_offset1 | 2009-2014 | mean_pairwise_nonidentical | 0.220758 |
| cds_nonoverlap_k3_offset1 | 2015-2019 | mean_pairwise_nonidentical | 0.220267 |
| cds_nonoverlap_k3_offset1 | 2020+ | mean_pairwise_nonidentical | 0.217159 |
| cds_nonoverlap_k3_offset1 | pre-2009 | mean_pairwise_nonidentical | 0.216823 |

## Robustness ranking

The ranking is a composite descriptive score combining higher JS distance, narrower JS confidence intervals, higher top-token Jaccard stability, narrower entropy intervals and coverage.

| rank | tokenizer | robustness_score | mean_js_distance | js_ci_width | mean_top_token_jaccard | n_sequences |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | raw_overlap_k3 | 0.623896 | 0.111653 | 0.00298811 | 0.927398 | 20000 |
| 2 | raw_overlap_k6 | 0.616149 | 0.688757 | 0.00144075 | 0.203458 | 20000 |
| 3 | cds_frame_k6 | 0.611779 | 0.895221 | 0.000924654 | 0.366073 | 8974 |
| 4 | cds_nonoverlap_k3_offset2 | 0.527566 | 0.219342 | 0.00260925 | 0.903914 | 8974 |
| 5 | raw_nonoverlap_k6 | 0.52469 | 0.789584 | 0.00160529 | 0.23725 | 20000 |
| 6 | cds_nonoverlap_k3_offset0 | 0.509747 | 0.21305 | 0.0027771 | 0.903651 | 8974 |
| 7 | cds_codon | 0.501144 | 0.213097 | 0.00294291 | 0.906035 | 8974 |
| 8 | cds_nonoverlap_k3_offset1 | 0.482838 | 0.210726 | 0.00328732 | 0.932619 | 8974 |

## Figures

- `results/figures/fig24_bootstrap_js_distance_ci.png`
- `results/figures/fig25_token_entropy_stability.png`
- `results/figures/fig26_top_token_jaccard_stability.png`
- `results/figures/fig27_temporal_token_entropy.png`
- `results/figures/fig28_tokenizer_robustness_ranking.png`

## What this phase can support

- We evaluate the stability of deterministic tokenization metrics under stratified bootstrap.
- k=6 tokenizers can be described as more or less stable based on bootstrap confidence intervals and top-token overlap.
- Top-token enrichment is assessed descriptively through Jaccard stability.
- Temporal stability is evaluated only where metadata supports it.
- This phase identifies robust deterministic baselines before learned tokenizers.

## What this phase does not support

- Stable tokens are not antigenic markers.
- Token stability does not predict evolution, escape, vaccine relevance or fitness.
- Temporal token changes do not imply selection.
- This phase does not validate GROVER or any learned tokenizer.

## Recommendation for Phase 7

Proceed to a learned-tokenizer audit only if the goal is methodological comparison. A conservative Phase 7 would introduce local BPE as an optional learned baseline and compare it against the robust deterministic baselines identified here, still without GROVER, prediction or biological efficacy claims.
