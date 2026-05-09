# FluGenome3D Phase 2 sequence-context report

This phase provides a descriptive sequence-context audit over the local Phase 1 panels. It compares compositional and codon-level patterns across HA/NA and H1N1/H3N2. It does not implement GROVER, BPE tokenization, 3D structure mapping, prediction, antigenicity, vaccine, escape, fitness, or sequence optimization.

## Panel coverage

- Smoke panel local metrics: `data/processed/metrics/smoke_sequence_metrics.parquet` and `data/processed/metrics/smoke_codon_metrics.parquet`.
- MVP panel local metrics: `data/processed/metrics/mvp_sequence_metrics.parquet` and `data/processed/metrics/mvp_codon_metrics.parquet`.
- MVP sequences analyzed: 20000.
- MVP codon/translation-QC reliable sequences: 8723.

## Composition summary

GC fraction, CpG observed/expected and UpA observed/expected are summarized as descriptive features. UpA is measured as TA in DNA alphabet.

| panel | subtype | protein | metric | n | mean | median | std | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mvp | H1N1 | HA | sequence_length | 5000 | 1726.44 | 1724 | 26.9794 | 1699.9 | 1777 |
| mvp | H1N1 | NA | sequence_length | 5000 | 1423.27 | 1418 | 15.3434 | 1410 | 1458 |
| mvp | H3N2 | HA | sequence_length | 5000 | 1721.38 | 1721 | 19.5223 | 1701 | 1762 |
| mvp | H3N2 | NA | sequence_length | 5000 | 1428.86 | 1431 | 18.7013 | 1410 | 1466 |
| mvp | H1N1 | HA | gc_content | 5000 | 0.408105 | 0.407534 | 0.00329822 | 0.40388 | 0.414494 |
| mvp | H1N1 | NA | gc_content | 5000 | 0.41728 | 0.41773 | 0.00347295 | 0.41118 | 0.421986 |
| mvp | H3N2 | HA | gc_content | 5000 | 0.417756 | 0.417436 | 0.00380236 | 0.411629 | 0.424456 |
| mvp | H3N2 | NA | gc_content | 5000 | 0.427813 | 0.427875 | 0.00271821 | 0.423207 | 0.432056 |
| mvp | H1N1 | HA | cpg_oe | 5000 | 0.35979 | 0.365502 | 0.027403 | 0.29857 | 0.392281 |
| mvp | H1N1 | NA | cpg_oe | 5000 | 0.386817 | 0.37708 | 0.0297099 | 0.348115 | 0.445748 |
| mvp | H3N2 | HA | cpg_oe | 5000 | 0.415592 | 0.414496 | 0.0259926 | 0.374981 | 0.463658 |
| mvp | H3N2 | NA | cpg_oe | 5000 | 0.378997 | 0.379202 | 0.0213998 | 0.343195 | 0.41203 |
| mvp | H1N1 | HA | upa_oe | 5000 | 0.697976 | 0.700386 | 0.031788 | 0.63961 | 0.744517 |
| mvp | H1N1 | NA | upa_oe | 5000 | 0.673665 | 0.67444 | 0.012567 | 0.651822 | 0.693034 |
| mvp | H3N2 | HA | upa_oe | 5000 | 0.584712 | 0.585323 | 0.0233703 | 0.546405 | 0.620013 |
| mvp | H3N2 | NA | upa_oe | 5000 | 0.649213 | 0.646259 | 0.0215937 | 0.617676 | 0.688386 |

## Dinucleotide odds ratios

The table below shows the highest mean dinucleotide odds ratios in the MVP panel by group. Full 16-dinucleotide summaries are in `results/tables/phase2_dinucleotide_odds_summary.csv`.

| subtype | protein | dinucleotide | n | mean | median | q05 | q95 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| H3N2 | NA | TG | 5000 | 1.43208 | 1.43199 | 1.39184 | 1.47473 |
| H3N2 | HA | TG | 5000 | 1.4148 | 1.41595 | 1.36603 | 1.45646 |
| H3N2 | NA | CA | 5000 | 1.37593 | 1.37666 | 1.33991 | 1.41101 |
| H1N1 | HA | TG | 5000 | 1.37239 | 1.37498 | 1.30339 | 1.44277 |
| H1N1 | NA | CA | 5000 | 1.36704 | 1.36802 | 1.32399 | 1.40413 |
| H1N1 | NA | TG | 5000 | 1.35837 | 1.35822 | 1.33179 | 1.3883 |
| H1N1 | HA | CA | 5000 | 1.31979 | 1.3258 | 1.25942 | 1.35333 |
| H3N2 | HA | CA | 5000 | 1.31251 | 1.31193 | 1.28577 | 1.33989 |
| H3N2 | HA | GG | 5000 | 1.27353 | 1.27133 | 1.23882 | 1.31364 |
| H1N1 | NA | GG | 5000 | 1.26682 | 1.26987 | 1.22796 | 1.29292 |
| H1N1 | HA | GG | 5000 | 1.25375 | 1.25591 | 1.20128 | 1.29086 |
| H1N1 | NA | TC | 5000 | 1.15901 | 1.1637 | 1.09547 | 1.2033 |
| H3N2 | HA | TT | 5000 | 1.14589 | 1.1429 | 1.11758 | 1.19011 |
| H3N2 | NA | CC | 5000 | 1.13621 | 1.13551 | 1.08938 | 1.18794 |
| H3N2 | NA | GG | 5000 | 1.11837 | 1.11624 | 1.08914 | 1.15413 |
| H1N1 | HA | TC | 5000 | 1.10876 | 1.11633 | 1.05735 | 1.14973 |

_Table truncated to 16 of 64 rows._

## Translation and codon QC

Codon usage and RSCU are reported only for sequences that pass the naive frame check, ambiguity check, translation check and internal-stop check. These checks use the available nucleotide strings as provided; they do not infer CDS boundaries.

| panel | subtype | protein | n_sequences | n_frame_fail | n_ambiguous_fail | n_internal_stop_fail | n_translation_fail | n_codon_reliable | reliable_fraction | mean_internal_stop_count |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mvp | H1N1 | HA | 5000 | 1154 | 0 | 2645 | 0 | 2354 | 0.4708 | 16.1648 |
| mvp | H1N1 | NA | 5000 | 1783 | 0 | 2499 | 0 | 2394 | 0.4788 | 14.3674 |
| mvp | H3N2 | HA | 5000 | 1540 | 0 | 3064 | 0 | 1922 | 0.3844 | 19.1378 |
| mvp | H3N2 | NA | 5000 | 1847 | 0 | 2936 | 0 | 2053 | 0.4106 | 14.0322 |

## Codon usage and RSCU

Codon usage summaries use translation-QC reliable sequences only. Stop codons are included in codon usage counts but excluded from RSCU.

Top codons by mean frequency:

| subtype | protein | codon | amino_acid | n_reliable_sequences | mean_frequency | median_frequency |
| --- | --- | --- | --- | --- | --- | --- |
| H1N1 | HA | AAT | N | 2354 | 0.0487018 | 0.0493827 |
| H1N1 | NA | AAT | N | 2394 | 0.0480809 | 0.0489362 |
| H1N1 | HA | AAA | K | 2354 | 0.0465162 | 0.047619 |
| H3N2 | HA | AAA | K | 1922 | 0.0459924 | 0.0440917 |
| H1N1 | NA | ATA | I | 2394 | 0.0445369 | 0.0446809 |
| H3N2 | HA | AAT | N | 1922 | 0.044204 | 0.0440917 |
| H1N1 | HA | GAA | E | 2354 | 0.0427983 | 0.042328 |
| H1N1 | HA | ACA | T | 2354 | 0.0426256 | 0.042328 |
| H3N2 | NA | ATA | I | 2053 | 0.0385691 | 0.0382979 |
| H1N1 | NA | GGA | G | 2394 | 0.0369036 | 0.0361702 |
| H1N1 | NA | TGG | W | 2394 | 0.0340325 | 0.0340426 |
| H3N2 | NA | AAT | N | 2053 | 0.0337943 | 0.0340426 |
| H3N2 | HA | AAC | N | 1922 | 0.0335404 | 0.0335097 |
| H1N1 | NA | AAC | N | 2394 | 0.032665 | 0.0319149 |
| H3N2 | NA | AAA | K | 2053 | 0.0323766 | 0.0319149 |
| H3N2 | HA | CAA | Q | 1922 | 0.0321255 | 0.031746 |

_Table truncated to 16 of 244 rows._

Top codons by mean RSCU:

| subtype | protein | codon | amino_acid | n_reliable_sequences | mean_rscu | median_rscu |
| --- | --- | --- | --- | --- | --- | --- |
| H1N1 | HA | AGA | R | 2354 | 4.4972 | 4.66667 |
| H1N1 | NA | AGA | R | 2394 | 3.50502 | 3.52941 |
| H3N2 | HA | AGA | R | 1922 | 2.98108 | 3 |
| H3N2 | NA | AGA | R | 2053 | 2.50032 | 2.45455 |
| H1N1 | HA | ACA | T | 2354 | 2.49247 | 2.48649 |
| H1N1 | HA | TCA | S | 2354 | 2.26975 | 2.29787 |
| H3N2 | NA | CCT | P | 2053 | 2.1344 | 2.10526 |
| H1N1 | NA | GCT | A | 2394 | 2.04579 | 2 |
| H3N2 | NA | TTG | L | 2053 | 2.03852 | 2.07692 |
| H3N2 | NA | AGG | R | 2053 | 1.99592 | 1.90909 |
| H1N1 | HA | GCA | A | 2354 | 1.99471 | 2.06061 |
| H1N1 | HA | CCA | P | 2354 | 1.95425 | 1.89474 |
| H1N1 | NA | CCA | P | 2394 | 1.92988 | 2 |
| H3N2 | HA | AGC | S | 1922 | 1.87886 | 1.90244 |
| H3N2 | HA | AGG | R | 1922 | 1.87594 | 1.92857 |
| H3N2 | NA | CAT | H | 2053 | 1.79631 | 1.8 |

_Table truncated to 16 of 244 rows._

## K-mer entropy

K-mer profiles were summarized for k=3,4,5. Versionable outputs keep entropy and top-kmer aggregate summaries only, not row-level sequence strings.

| panel | subtype | protein | k | n_sequences | mean_entropy | median_entropy | mean_valid_kmers | top1_kmer | top1_frequency | top2_kmer | top2_frequency | top3_kmer | top3_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mvp | H1N1 | HA | 3 | 5000 | 5.77202 | 5.77278 | 1724.39 | AAA | 0.0516589 | GAA | 0.0321563 | AAT | 0.0320417 |
| mvp | H1N1 | HA | 4 | 5000 | 7.59622 | 7.59639 | 1723.38 | AAAA | 0.0199223 | GAAA | 0.0134196 | AAAT | 0.0130836 |
| mvp | H1N1 | HA | 5 | 5000 | 9.15873 | 9.15787 | 1722.36 | AAAAA | 0.00765995 | GGAAA | 0.00739833 | AAGAA | 0.00576418 |
| mvp | H1N1 | NA | 3 | 5000 | 5.80614 | 5.80609 | 1421.24 | AAT | 0.0380207 | TGG | 0.0329407 | AAA | 0.0315844 |
| mvp | H1N1 | NA | 4 | 5000 | 7.63003 | 7.63015 | 1420.23 | TGGA | 0.012409 | AATG | 0.0119304 | ACAA | 0.0110182 |
| mvp | H1N1 | NA | 5 | 5000 | 9.14971 | 9.14961 | 1419.23 | AGACA | 0.00509137 | AATGG | 0.0050467 | ACAAT | 0.00461012 |
| mvp | H3N2 | HA | 3 | 5000 | 5.78731 | 5.78664 | 1719.32 | AAA | 0.0472469 | CAA | 0.0366156 | AAT | 0.0363861 |
| mvp | H3N2 | HA | 4 | 5000 | 7.61633 | 7.61494 | 1718.3 | CAAA | 0.0161662 | AAAA | 0.0160881 | AAAT | 0.0155154 |
| mvp | H3N2 | HA | 5 | 5000 | 9.17549 | 9.17348 | 1717.29 | AAATG | 0.00679245 | CAAAA | 0.0062148 | AACAA | 0.00567465 |
| mvp | H3N2 | NA | 3 | 5000 | 5.83118 | 5.83106 | 1426.77 | AAA | 0.0342685 | CAA | 0.0301957 | ATG | 0.0277085 |
| mvp | H3N2 | NA | 4 | 5000 | 7.66631 | 7.66619 | 1425.75 | GAAA | 0.0108585 | AAAA | 0.0102671 | AACA | 0.0100328 |
| mvp | H3N2 | NA | 5 | 5000 | 9.18098 | 9.17933 | 1424.72 | GGAAA | 0.00521308 | ACAAT | 0.0042039 | TGATG | 0.00418762 |

## Figures

- `results/figures/fig4_gc_cpg_upa_by_group.png`: GC, CpG O/E and UpA O/E by subtype/protein group.
- `results/figures/fig5_dinucleotide_odds_heatmap.png`: 16-dinucleotide mean odds-ratio heatmap.
- `results/figures/fig6_codon_usage_rscu_heatmap.png`: codon/translation QC summary, because at least one MVP group was below the codon-QC threshold.
- `results/figures/fig7_kmer_entropy_by_group.png`: k-mer entropy by group for k=3,4,5.

## What this phase can say

- This phase provides a descriptive sequence-context audit.
- We compare compositional and codon-level patterns across HA/NA and H1N1/H3N2.
- CpG/UpA and dinucleotide odds ratios are summarized as descriptive features.
- Codon usage is reported only after translation/frame QC.

## What this phase cannot say

- It does not predict antigenic drift.
- It does not identify escape mutations.
- It does not explain pathogenicity.
- It does not predict vaccine candidates.
- It does not claim codon usage explains fitness.
- It does not claim CpG/UpA determines antigenicity.
- It does not make causal claims.

## Limitations

- The translation QC is a naive check on the available nucleotide strings and does not infer curated CDS start/end coordinates.
- Ambiguous bases reduce valid k-mer and codon counts.
- The full panel remains a later extension; Phase 2 figures are anchored on the MVP panel.
- Aggregated summaries are public-safe, while row-level metrics remain local/gitignored.
