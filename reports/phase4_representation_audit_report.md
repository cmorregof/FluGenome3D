# FluGenome3D Phase 4 representation audit report

This phase compares simple nucleotide and codon-level representations of Influenza A HA/NA before any GROVER or BPE tokenization. The analysis is descriptive and representation-focused.

## Dataset separation

- Raw nucleotide representations use `mvp_panel`, because k-mer and compositional features do not require a validated CDS frame.
- CDS/codon representations use `mvp_cds_refined_panel` and `mvp_cds_refined_codon_metrics`, because codon-frequency and RSCU features require explicit CDS QC.

## Representations built

| representation | source | n_sequences | n_features | nnz | sparsity | pca_explained_variance_pc1 | pca_explained_variance_pc2 | pca_explained_variance_total_2pc | umap_method |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| kmer3_freq | raw_nucleotide | 20000 | 64 | 1279894 | 8.28125e-05 | 0.46206 | 0.212951 | 0.675011 |  |
| kmer4_freq | raw_nucleotide | 20000 | 256 | 4871551 | 0.0485252 | 0.35046 | 0.232325 | 0.582785 |  |
| kmer5_freq | raw_nucleotide | 20000 | 1024 | 13918094 | 0.320406 | 0.277943 | 0.256372 | 0.534315 |  |
| kmer3_tfidf | raw_nucleotide | 20000 | 64 | 1279894 | 8.28125e-05 | 0.456997 | 0.217556 | 0.674553 |  |
| kmer4_tfidf | raw_nucleotide | 20000 | 256 | 4871551 | 0.0485252 | 0.347838 | 0.242105 | 0.589943 | pca_fallback_for_umap |
| codon_freq | cds_refined | 8974 | 64 | 531774 | 0.0741064 | 0.477261 | 0.205298 | 0.682559 |  |
| rscu | cds_refined | 8974 | 61 | 522800 | 0.0449641 | 0.461943 | 0.244373 | 0.706316 |  |

## PCA and UMAP summary

PCA was run for all representations. UMAP was attempted for `kmer4_tfidf` on a stratified sample; if UMAP is unavailable, the script records and plots a deterministic PCA fallback under the UMAP figure filename.

## Silhouette scores

Silhouette scores are computed on 2D PCA embeddings for protein, subtype and protein-subtype labels. These scores are descriptive separation diagnostics, not predictive performance metrics.

| representation | label_type | silhouette | space |
| --- | --- | --- | --- |
| kmer3_freq | protein_subtype | 0.834207 | pca_2d |
| kmer3_tfidf | protein_subtype | 0.83179 | pca_2d |
| kmer4_freq | protein_subtype | 0.821966 | pca_2d |
| rscu | protein_subtype | 0.818114 | pca_2d |
| kmer4_tfidf | protein_subtype | 0.785164 | pca_2d |

## Centroid distances

Pairwise centroid distances are stored in `results/tables/phase4_group_centroid_distances.csv`. Figure 16 visualizes group distances for k-mer TF-IDF, codon frequency and RSCU representations.

| representation | group_a | group_b | distance |
| --- | --- | --- | --- |
| kmer3_freq | HA-H1N1 | HA-H1N1 | 0 |
| kmer3_freq | HA-H1N1 | HA-H3N2 | 0.0276674 |
| kmer3_freq | HA-H1N1 | NA-H1N1 | 0.0371195 |
| kmer3_freq | HA-H1N1 | NA-H3N2 | 0.0391377 |
| kmer3_freq | HA-H3N2 | HA-H1N1 | 0.0276674 |
| kmer3_freq | HA-H3N2 | HA-H3N2 | 0 |
| kmer3_freq | HA-H3N2 | NA-H1N1 | 0.0354595 |
| kmer3_freq | HA-H3N2 | NA-H3N2 | 0.03488 |
| kmer3_freq | NA-H1N1 | HA-H1N1 | 0.0371195 |
| kmer3_freq | NA-H1N1 | HA-H3N2 | 0.0354595 |
| kmer3_freq | NA-H1N1 | NA-H1N1 | 0 |
| kmer3_freq | NA-H1N1 | NA-H3N2 | 0.0288007 |

## Top features

Top k-mer, codon-frequency and RSCU features are aggregated by group. These are high-weight descriptive features, not mutation, fitness, antigenicity or escape markers.

| representation | protein | subtype | group | rank | feature | mean_value | k |
| --- | --- | --- | --- | --- | --- | --- | --- |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 1 | AAA | 0.0516294 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 2 | GAA | 0.0321504 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 3 | AAT | 0.0320418 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 4 | CAA | 0.0314341 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 5 | ACA | 0.0305524 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 6 | ATG | 0.029231 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 7 | AGA | 0.0286199 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 8 | AAG | 0.0272342 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 9 | TGG | 0.0254875 | 3 |
| kmer3_freq | HA | H1N1 | HA|H1N1 | 10 | ATT | 0.0252658 | 3 |
| kmer3_freq | HA | H3N2 | HA|H3N2 | 1 | AAA | 0.0472405 | 3 |
| kmer3_freq | HA | H3N2 | HA|H3N2 | 2 | CAA | 0.0366196 | 3 |

## Figures

- `results/figures/fig12_kmer_pca_by_group.png`
- `results/figures/fig13_kmer_umap_by_group.png`
- `results/figures/fig14_codon_pca_by_group.png`
- `results/figures/fig15_rscu_pca_by_group.png`
- `results/figures/fig16_group_centroid_distance_heatmap.png`
- `results/figures/fig17_representation_silhouette_comparison.png`

## Limitations

- UMAP is optional and may fall back to PCA if the dependency is not installed.
- CDS-dependent representations cover only the refined CDS subset.
- Feature matrices are local/gitignored and should not be treated as public sequence redistribution artifacts.
- This phase does not make claims about antigenic drift, vaccine relevance, escape, fitness, causality or prediction.

## Recommendation for Phase 5

Use these simple representations as baselines for tokenizer audits. The next phase can compare codon/k-mer tokenization behavior against these matrices before considering GROVER or BPE.
