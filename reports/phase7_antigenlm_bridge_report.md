# Phase 7 AntigenLM bridge report

Phase 7 imports the parent-repository AntigenLM geometry as a derived representation layer. It does not export raw sequences, accessions, isolate names, sequence hashes, FASTA, checkpoints, or restricted Parquet records.

## Cache summary

- Records represented: 111,756
- Embedding dimension: 384
- Year range: 2000-2022
- Exported latent atlas points: 24,000
- Projection: pca_3d

## Why this matters

The child repo already compares nucleotide, codon and token-level baselines. Phase 7 adds the learned AntigenLM layer from the thesis repo so the app can show how simple biological features relate to a learned influenza representation.

## Molecular geometry summary

| metric | subtype | rho_mean | rho_sd | valid_pairs_mean | omitted_pairs_mean | n_runs |
| --- | --- | --- | --- | --- | --- | --- |
| hamming_ha | H1N1 | 0.8280156332049321 | 0.0010822611171559169 | 199914.33333333334 | 85.66666666666667 | 3 |
| hamming_ha | H3N2 | 0.6081682079220184 | 0.0006412968546132165 | 199988.66666666666 | 11.333333333333334 | 3 |
| hamming_ha_na | H1N1 | 0.8537804922836947 | 0.0007202431218106303 | 199903.0 | 97.0 | 3 |
| hamming_ha_na | H3N2 | 0.6677909555897313 | 0.001661629990917984 | 199537.66666666666 | 462.3333333333333 | 3 |
| temporal | H1N1 | 0.06117937914963284 | 0.001080895044009855 | 200000.0 | 0.0 | 3 |
| temporal | H3N2 | 0.15397261168018747 | 0.0009366835297727553 | 200000.0 | 0.0 | 3 |

## PCA effective dimension

| group | n | n80 | n90 | n95 | n99 | participation_ratio | top1_evr | top2_evr |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| global | None | None | None | None | None | 1.910691053477228 | None | None |
| H1N1 | None | None | None | None | None | 1.4305587487342268 | None | None |
| H3N2 | None | None | None | None | None | 1.5294304210698815 | None | None |

## Clade enrichment summary

| subtype | label | k | classes | n_labeled | mean_precision | random_baseline | enrichment_vs_random | permutation_p05 | permutation_p95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1N1 | clade | 5 | 17 | 36156 | 0.9148749089366567 | 0.2408784157539551 | 3.798077574003784 | 0.23864078162318683 | 0.24265750500272038 |
| H1N1 | clade | 10 | 17 | 36156 | 0.893256569416128 | 0.24028284839768413 | 3.717521143821838 | 0.2389599637242145 | 0.24211028554394798 |
| H1N1 | clade | 20 | 17 | 36156 | 0.8646618682189291 | 0.24042021610060113 | 3.5964607396289883 | 0.23923635647192779 | 0.24166657977196243 |
| H3N2 | clade | 5 | 39 | 40082 | 0.8609146899708697 | 0.06893701245779486 | 12.48842471231177 | 0.06799296712442779 | 0.06984579692051603 |
| H3N2 | clade | 10 | 39 | 40082 | 0.8324670570140091 | 0.06893950734327961 | 12.07532645785813 | 0.06828248258849577 | 0.06971586847925398 |
| H3N2 | clade | 20 | 39 | 40082 | 0.7927227208488915 | 0.06884761239459108 | 11.514164301087233 | 0.06832653970483632 | 0.06948789259575315 |
| H1N1 | major_clade | 5 | 5 | 36156 | 0.9735897861509946 | 0.491794815060663 | 1.9796666340024385 | 0.48877467931870794 | 0.4937940216338839 |
| H1N1 | major_clade | 10 | 5 | 36156 | 0.9658594942199324 | 0.4907899103883173 | 1.967969336320985 | 0.4893657088024206 | 0.49317243553947004 |
| H1N1 | major_clade | 20 | 5 | 36156 | 0.9559079307565453 | 0.4909337316074787 | 1.9471221250709092 | 0.4896987882822127 | 0.49252471429660594 |
| H3N2 | major_clade | 5 | 16 | 40082 | 0.9048988764044943 | 0.18176405036342166 | 4.978426012158216 | 0.17984885559717018 | 0.18331077819392425 |
| H3N2 | major_clade | 10 | 16 | 40082 | 0.8853547537031886 | 0.18135239425843686 | 4.881957899279294 | 0.1802050621140407 | 0.1827235622087352 |
| H3N2 | major_clade | 20 | 16 | 40082 | 0.8578913788927558 | 0.1812929328210502 | 4.732072924980255 | 0.18040786554401045 | 0.18241834901427964 |

## App boundary

- The app may show reduced coordinates with hash-based IDs and minimal metadata.
- The app may show aggregate Spearman, PCA, TwoNN, temporal-locality and clade-enrichment summaries.
- The app must not show raw sequences, source accessions, isolate names, sequence hashes or checkpoints.

## Claims allowed

- AntigenLM embeddings are summarized as a learned representation layer.
- Latent distances can be compared descriptively with molecular and temporal proxies.
- The learned layer is compared against deterministic k-mer/codon/RSCU baselines.

## Claims not allowed

- This phase does not predict antigenicity, vaccine relevance, escape, pathogenicity, fitness or evolution.
- This phase does not validate sequence generation or optimize sequences.
