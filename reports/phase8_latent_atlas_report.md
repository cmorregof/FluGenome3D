# Phase 8 Latent Atlas report

Phase 8 turns Phase 7 outputs into a visual atlas concept for the app: deterministic nucleotide/codon representations and the learned AntigenLM layer are shown as complementary views of the same HA/NA research landscape.

## Representation families

| family | representation | n_sequences | n_features | pca_2pc_variance | protein_subtype_silhouette | molecular_proxy | app_role | global_pca_n95 | global_participation_ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| deterministic_baseline | kmer3_freq | 20000 | 64 | 0.6750106811523438 | 0.8342072563188496 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | kmer4_freq | 20000 | 256 | 0.5827852636575699 | 0.8219655935456675 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | kmer5_freq | 20000 | 1024 | 0.5343146920204163 | 0.7388535435339182 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | kmer3_tfidf | 20000 | 64 | 0.6745526492595673 | 0.8317900525093036 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | kmer4_tfidf | 20000 | 256 | 0.5899432450532913 | 0.7851643202087798 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | codon_freq | 8974 | 64 | 0.6825591027736664 | 0.7848540089786886 | not_computed_here | baseline comparison | nan | nan |
| deterministic_baseline | rscu | 8974 | 61 | 0.7063157707452774 | 0.8181137504266337 | not_computed_here | baseline comparison | nan | nan |
| learned_latent | AntigenLM HA+NA embeddings | 111756 | 384 | nan | nan | HA+NA Hamming rho mean 0.761 | learned biological representation layer | nan | 1.910691053477228 |

## Design decision

The AntigenLM layer should be foregrounded as a learned representation, while k-mer, codon-frequency and RSCU spaces remain transparent baselines. This avoids treating learned embeddings as magic and gives interviewers a clear ladder from interpretable biology to learned geometry.

## App behavior

- The Representation Projector remains the baseline feature-space explorer.
- The new Latent Atlas view focuses on AntigenLM PCA points, molecular-geometry summaries, clade enrichment and temporal locality.
- All displayed point IDs are hash-based internal IDs with minimal metadata.

## Boundaries

This phase does not add GROVER, BPE, prediction, antigenicity, escape, vaccine, pathogenicity, fitness or sequence optimization claims.
