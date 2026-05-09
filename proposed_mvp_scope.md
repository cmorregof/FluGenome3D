# FluGenome3D Fase 0: proposed MVP scope

This scope is descriptive only. It makes no predictive, generative, antigenicity, vaccine, escape, transmissibility, pathogenicity, or optimization claims.

## Recommended data products

### 1. Smoke test panel

Purpose: run the full pipeline quickly and verify file handling, QC, plotting, and structure mapping.

Recommended local panel:

- 24 paired strains = 48 sequence records.
- Balanced by subtype: 12 H1N1 and 12 H3N2 paired strains.
- Each paired strain contributes one HA and one NA.
- Prefer exact HA+NA deduplicated records from the 82,306 dedup cache.
- Prefer records with rich metadata available: host, location, collection date, and clade.
- Stratify across broad time bins when available: 2000-2005, 2006-2011, 2012-2016, 2017-2022.
- Use a fixed random seed, e.g. 42.

Governance:

- Keep any smoke panel derived from GISAID local/gitignored.
- A public smoke test for GitHub should use redistributable non-GISAID examples or a script that tells users how to fetch their own public records.

### 2. MVP CV-ready panel

Purpose: produce robust descriptive figures suitable for a CV/project page without relying on the full dataset.

Recommended local panel:

- 10,000 paired strains = 20,000 sequence records.
- Balanced by subtype: 5,000 H1N1 and 5,000 H3N2 paired strains.
- Use exact HA+NA deduplicated records where possible.
- Prefer records in the rich metadata subset so host/location/date/clade are available.
- Stratify by subtype x year, with caps per year to avoid dominance by densely sampled years.
- For clade plots, use only records with nonmissing `major_clade`; report clade coverage explicitly.

Why 10,000 paired strains:

- Available local data supports thousands of paired records per subtype.
- It is small enough for iterative figures, tokenization diagnostics, and structure mapping.
- It avoids letting exact duplicate HA/NA pairs dominate representation summaries.

### 3. Full local analysis panel

Purpose: final local-only descriptive analysis after the MVP is approved.

Recommended local panels:

- `full_valid_paired`: all 111,756 valid paired strains from `../data/processed_gisaid/dataset_H1N1.json` and `dataset_H3N2.json`.
- `full_dedup_pair`: all 82,306 exact HA+NA deduplicated paired strains.
- `metadata_rich_subset`: 81,943 records with rich host/location/date/clade metadata.

Use cases:

- Use `full_valid_paired` for sampling, year coverage, and local thesis continuity.
- Use `full_dedup_pair` for representation, distance, and duplicate-sensitive descriptive summaries.
- Use `metadata_rich_subset` for host/location/clade figures.

## GitHub/governance boundary

### Can go to GitHub

- Source code, tests, documentation, and config examples.
- `data_manifest/dataset_doi.txt` with allowed acknowledgement/DOI text only.
- Aggregate inventory tables with no sequence strings, no isolate names, no accession rows, and no row-level GISAID metadata.
- Aggregate figures by subtype/segment/year/clade when bins are sufficiently coarse and not reconstructable.
- PDB identifiers and scripts/configs that download structures from public sources.
- Public demo instructions or tiny redistributable examples not derived from restricted GISAID data.

### Must stay local/gitignored

- Raw FASTA files from `../data/gisaid`.
- Original `.xls` and `.csv` metadata exports from `../data/gisaid_metadata_private`.
- EPI_ISL/accession-level tables, even if they contain only identifiers.
- Processed JSON files containing `ha_sequence` and `na_sequence`.
- Any Parquet/CSV/JSON table with raw sequences, accession-level rows, isolate names, or per-record restricted metadata.
- BLAST databases and local indexes if added later.
- Downloaded structure files if kept under `results/structures`, as currently gitignored; the public project can keep PDB IDs instead.

## Exact MVP figures

### Figure 1: `fig01_qc_inventory_by_subtype_segment`

What it shows:

- Raw FASTA records by subtype and segment.
- Near-complete records by subtype and segment.
- Valid paired records by subtype.
- Exact HA+NA deduplicated paired records.
- Rich metadata join coverage.

Rationale:

- Establishes the real data basis before any biological interpretation.
- Makes duplicate and metadata coverage visible.

### Figure 2: `fig02_length_and_ambiguity_qc`

What it shows:

- HA and NA length distributions by subtype.
- Expected local filter windows for HA and NA.
- Counts or fractions of ambiguous and non-ACGTN records.

Rationale:

- Supports sequence-context and QC decisions.
- Shows why the processed paired JSONs are safer than raw FASTA for MVP metrics.

### Figure 3: `fig03_codon_usage_rscu_by_subtype_segment_year`

What it shows:

- Codon usage or RSCU summaries by subtype, segment, and year bin.
- Optional second panel: PCA/UMAP of aggregate codon usage vectors, labeled only by subtype/segment/year bin.

Rationale:

- Directly addresses codon usage without optimizing or recommending sequences.
- Uses nucleotide CDS-like HA/NA data already present locally.

### Figure 4: `fig04_cpg_upa_dinucleotide_bias`

What it shows:

- CpG and UpA observed/expected ratios by subtype, segment, and year bin.
- Companion 16-dinucleotide heatmap aggregated by subtype and segment.

Rationale:

- Directly addresses CpG/UpA/dinucleotide bias.
- Keeps interpretation descriptive.

### Figure 5: `fig05_tokenizer_and_structure_context`

What it shows:

- Left panel: tokenizer-dependent diagnostics for codon, fixed k-mer, and optional BPE tokenizers: tokens/kb, entropy, fraction crossing codon boundaries, fraction containing CpG/UpA.
- Right panel or companion HTML: structure-aware mapping of one descriptive sequence-context signal onto representative HA/NA structures.

Rationale:

- Keeps tokenizer behavior in the MVP without making it the scientific bottleneck.
- Connects nucleotide context to protein coordinates only as a visualization, not as a functional claim.

## GROVER tokenizer decision

Decision: GROVER itself is Phase 2; tokenizer-dependent audit is optional but allowed in the MVP.

MVP core tokenizers:

- Codon tokenizer.
- Fixed k-mer tokenizers, e.g. k=3,4,6 with overlapping and non-overlapping modes.
- Local BPE tokenizer only as an exploratory diagnostic if it is trained locally on allowed data and kept gitignored.

Why GROVER is not core for MVP:

- No actual GROVER tokenizer/model artifact was found as a local dependency for this repo.
- The current FluGenome3D code treats GROVER dependencies as optional (`grover = ["transformers", "torch"]`).
- The central MVP claims can be completed with interpretable codon/k-mer tokenizers.
- Adding real GROVER should be done only after licensing, model provenance, and exact tokenizer behavior are documented.

## Structure mapping decision

Signal to map:

- Primary MVP signal: local windowed CpG/UpA log-odds or CpG/UpA observed/expected ratio aggregated per codon/residue position.
- Secondary signal if needed: GC3 fraction or dinucleotide-bias residual.

Mapping rule:

- Compute sequence-context metrics in nucleotide coordinates.
- Convert codon index to amino-acid residue index.
- Align the subtype/segment consensus or reference sequence to the PDB polymer sequence.
- Map only aligned positions with unambiguous residue correspondence.
- Leave gaps/unmapped residues gray.

Representative structures:

- H1 HA: PDB `3LZG`, 2009 H1N1 hemagglutinin, A/California/04/2009.
- H3 HA: PDB `3VUN`, H3N2 hemagglutinin, A/Aichi/2/1968.
- N1 NA: PDB `3NSS`, 2009 H1N1 neuraminidase, A/California/04/2009.
- N2 NA: PDB `6BR6`, H3N2 neuraminidase, A/Perth/16/2009.

RCSB references:

- [3LZG](https://www.rcsb.org/structure/3LZG)
- [3VUN](https://www.rcsb.org/structure/3VUN)
- [3NSS](https://www.rcsb.org/structure/3NSS)
- [6BR6](https://www.rcsb.org/structure/6BR6)

Interpretation boundary:

- The structure figure is a descriptive coordinate map of sequence-context metrics.
- It must not be described as antigenicity, escape, fitness, vaccine, binding, or functional prediction.
