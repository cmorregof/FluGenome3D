# FluGenome3D Fase 0: implementation plan

This plan starts only after explicit approval of Fase 0. No Fase 1 analysis code should be added before that approval.

## Current Fase 0 status

Completed in this audit:

- Identified local data roots in the parent repository.
- Counted FASTA records, processed paired records, metadata rows, missing fields, ambiguity, length ranges, and exact duplicates.
- Proposed smoke, MVP, and full local panels.
- Defined GitHub-safe vs local-only outputs.
- Chose MVP figures.
- Decided GROVER status.
- Chose the first structure-mapping signal and representative PDB IDs.

Not done in Fase 0:

- No analysis modules implemented.
- No raw data copied into FluGenome3D.
- No sequences generated, optimized, transformed for design, or redistributed.
- No predictive or antigenicity claims.

## Phase 1 gate: local configuration only

Goal: let FluGenome3D read parent data without copying it.

Tasks after approval:

1. Create a local, uncommitted `config/local_paths.yml`.
2. Point it to parent data roots:
   - `../data/gisaid`
   - `../data/gisaid_metadata_private`
   - `../data/processed_gisaid`
3. Keep `config/local_paths.yml` gitignored if it contains private absolute paths.
4. Add a read-only path check to fail before writing outside FluGenome3D.

Expected outputs:

- No public data outputs.
- Optional local log confirming accessible files and sizes.

## Phase 2: inventory and schema loader

Goal: reproducible inventory without exposing raw rows.

Tasks:

1. Parse FASTA headers into a local restricted manifest.
2. Parse combined metadata CSV and dedup cache.
3. Parse processed JSON top-level structure and record counts.
4. Produce an aggregate inventory table under `results/tables`.
5. Produce local-only row-level Parquet files under `data/interim` if needed.

Validation:

- Counts must match Fase 0 within expected tolerance.
- No sequence strings or accessions in public aggregate outputs.

## Phase 3: dataset builder

Goal: define smoke, MVP, full-valid, and full-dedup panels.

Tasks:

1. Build `smoke_local`: 24 paired strains, 48 sequences.
2. Build `mvp_balanced_local`: 10,000 paired strains, 20,000 sequences.
3. Build `full_valid_paired_local`: 111,756 paired strains.
4. Build `full_dedup_pair_local`: 82,306 exact HA+NA pairs.
5. Track selection rules in local manifests.

Validation:

- Balanced subtype counts for smoke and MVP.
- Year-bin coverage reported.
- Exact duplicate rates reported before and after dedup.
- Host/location/clade coverage reported separately.

## Phase 4: sequence QC and sequence-context metrics

Goal: descriptive sequence context only.

Tasks:

1. Compute length, GC, GC1/GC2/GC3, N fraction, non-ACGTN counts.
2. Compute valid CDS frame checks where applicable.
3. Compute aggregate summaries by subtype, segment, year bin, and clade where metadata coverage allows.
4. Generate Figure 1 and Figure 2.

Validation:

- Processed lengths should remain within audited filter ranges.
- No row-level restricted output committed.

## Phase 5: codon usage and dinucleotide bias

Goal: codon and dinucleotide descriptive audit.

Tasks:

1. Compute codon counts, codon frequencies, and RSCU.
2. Compute CpG, UpA, and all 16 dinucleotide observed/expected ratios.
3. Aggregate by subtype, segment, year bin, and optionally major clade.
4. Generate Figure 3 and Figure 4.

Validation:

- Exclude or flag sequences with frame problems.
- Report n sequences per aggregate group.
- Do not infer fitness, antigenicity, immune escape, or optimization.

## Phase 6: tokenizer-dependent audit

Goal: compare tokenization behavior without making language-model claims.

MVP status:

- Codon and fixed k-mer tokenizers are allowed in MVP.
- Local BPE is optional and exploratory.
- Actual GROVER tokenizer/model integration is Phase 2 unless model/tokenizer provenance and licensing are documented.

Tasks:

1. Compare codon tokenizer vs k-mer tokenizers.
2. Optionally train a local BPE tokenizer only on local restricted data and keep tokenizer artifacts gitignored.
3. Summarize tokens/kb, entropy, token length, codon-boundary crossing, CpG-token fraction, and UpA-token fraction.
4. Generate the tokenizer panel of Figure 5.

Validation:

- Do not train or fine-tune a predictive model.
- Do not claim that tokenizers identify antigenic, vaccine, or escape features.

## Phase 7: lightweight structure-aware visualization

Goal: map a descriptive sequence-context scalar onto representative public structures.

Primary signal:

- Windowed CpG/UpA log-odds or observed/expected ratio per codon/residue coordinate.

Representative PDB IDs:

- H1 HA: `3LZG`
- H3 HA: `3VUN`
- N1 NA: `3NSS`
- N2 NA: `6BR6`

Tasks:

1. Download structures from public PDB sources into gitignored `results/structures`.
2. Align consensus/reference protein sequence to PDB polymer sequence.
3. Map only unambiguous aligned residues.
4. Render static PNG and/or local HTML visualization.
5. Generate the structure panel of Figure 5.

Validation:

- Gray out unmapped residues.
- Report alignment coverage.
- Keep interpretation as spatial visualization of sequence-context metrics only.

## Phase 8: reproducibility and packaging

Goal: make the project usable as a satellite repo without redistributing restricted data.

Tasks:

1. Add a reproducible command sequence for local mode.
2. Add smoke tests using synthetic or public redistributable toy data only.
3. Add checks that public outputs contain no sequences, no EPI_ISL rows, and no isolate names.
4. Document the exact local data DOI/acknowledgement text in `data_manifest/dataset_doi.txt` only if redistribution terms allow that text.

Validation:

- Run tests.
- Run smoke pipeline.
- Inspect `git status` before any commit.
- Confirm that `data/raw`, `data/interim`, `data/processed`, private metadata, and structures remain gitignored.

## Stop conditions

Stop and ask for approval if any task would:

- Copy GISAID raw data into the child repository.
- Commit accession-level or sequence-level files.
- Require installing or using an external model/tokenizer with unclear license.
- Move from descriptive audit into prediction, sequence generation, codon optimization, vaccine, antigenicity, escape, or design claims.
