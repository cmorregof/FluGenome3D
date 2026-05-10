# Phase 9 structure mapping QC report

Phase 9 creates the first explicit bridge from refined HA/NA CDS positions to public PDB polymer sequences. It is an alignment-QC and residue-signal layer, not a validated antigenic or functional residue map.

## What was built

- Public RCSB FASTA sequences were loaded for 3LZG, 3VUN, 3NSS and 6BR6.
- Refined CDS sequences were translated locally, without exporting sequences.
- Per-residue aggregate signals were computed: codon GC fraction, CpG codon fraction, UpA/TA codon fraction and amino-acid entropy.
- Local consensus amino-acid sequences were aligned to each public PDB polymer sequence.

## Mapping QC

| pdb_id | protein | subtype | pdb_entity | chains | pdb_sequence_length | local_consensus_length | identity | mapped_residues | coverage_pdb | coverage_local | local_start | local_end | mapping_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3LZG | HA | H1N1 | 3LZG_1 | Chains A, C, E, G, I, K | 329 | 573 | 0.9694189602446484 | 327 | 0.993920972644377 | 0.5706806282722513 | 18 | 344 | alignment_qc_available |
| 3LZG | HA | H1N1 | 3LZG_2 | Chains B, D, F, H, J, L | 177 | 573 | 0.9827586206896551 | 174 | 0.9830508474576272 | 0.3036649214659686 | 345 | 518 | alignment_qc_available |
| 3VUN | HA | H3N2 | 3VUN_1 | Chains A, C, E | 329 | 577 | 0.8054711246200608 | 329 | 1.0 | 0.5701906412478336 | 17 | 345 | alignment_qc_available |
| 3VUN | HA | H3N2 | 3VUN_2 | Chains B, D, F | 175 | 577 | 0.9542857142857143 | 175 | 1.0 | 0.30329289428076256 | 346 | 520 | alignment_qc_available |
| 3NSS | NA | H1N1 | 3NSS_1 | Chains A, B | 388 | 477 | 0.979381443298969 | 388 | 1.0 | 0.8134171907756813 | 82 | 469 | alignment_qc_available |
| 6BR6 | NA | H3N2 | 6BR6_1 | Chain A | 387 | 476 | 0.9844961240310077 | 387 | 1.0 | 0.8130252100840336 | 83 | 469 | alignment_qc_available |

## Residue-signal catalog

| subtype | protein | group | n_local_positions | n_mapped_positions | mean_gc_fraction_codon | mean_cpg_codon_fraction | mean_upa_codon_fraction | mean_aa_entropy | max_aa_entropy |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1N1 | HA | HA-H1N1 | 573 | 501 | 0.4059569629894852 | 0.017649717110994168 | 0.12997203927223194 | 0.44978960064921686 | 1.7794533569415123 |
| H1N1 | NA | NA-H1N1 | 477 | 388 | 0.41692828430134393 | 0.019874572897162622 | 0.1142137648811306 | 0.5123695991933372 | 1.8285363438280822 |
| H3N2 | HA | HA-H3N2 | 577 | 504 | 0.4167718183338029 | 0.027096722449079222 | 0.09256928050103044 | 0.3775318982803106 | 1.954434002924965 |
| H3N2 | NA | NA-H3N2 | 476 | 387 | 0.42608076880383317 | 0.02449146042442076 | 0.10035435599794695 | 0.21414215797956393 | 1.8020456244375826 |

## Interpretation

The mapping is now more than pending: alignment QC is available. However, residue coloring in the 3D viewer still requires an additional chain/residue-number validation layer before FluGenome3D metrics are painted onto atoms.

## Boundaries

- These mapped signals are descriptive sequence-context summaries.
- They are not antigenic sites, vaccine markers, escape sites, pathogenicity markers, fitness estimates or causal explanations.
- No raw sequences, FASTA records, accessions or isolate names are exported.
