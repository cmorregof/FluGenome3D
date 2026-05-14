# Methods

FluGenome3D is a geometric audit and research-visualization lab for thesis-derived Influenza A HA/NA artifacts. The methods are organized around conservative inspection layers rather than biological forecasting.

## 1. Data basis

The local analysis starts from paired HA/NA records and associated metadata available in the thesis workspace. Public app views use only derived summaries of:

- paired HA/NA records;
- subtype and segment/protein groups;
- temporal and geographic bins;
- QC windows for length, ambiguity, duplication and metadata completeness;
- CDS reliability where codon-dependent summaries are used.

Raw FASTA, accessions, isolate names and restricted per-record metadata remain local.

## 2. Sequence-context metrics

Sequence-context summaries are descriptive aggregate metrics. They include:

- GC fraction;
- CpG observed/expected;
- UpA observed/expected, using DNA TA as an RNA UpA proxy;
- dinucleotide summaries;
- codon and RSCU summaries where CDS/frame QC permits;
- entropy and aggregate distribution summaries.

These metrics describe composition and representation behavior. They do not establish antigenicity, immune relevance, pathogenicity or fitness.

## 3. Tokenization audit

The tokenization audit evaluates deterministic tokenizers before learned tokenizers:

- overlapping k-mers;
- non-overlapping k-mers;
- codon tokenization;
- frame-aware variants;
- bootstrap and stability summaries where available.

The goal is to document how segmentation choices alter token distributions, effective vocabulary, entropy and group-level distances.

## 4. Representation geometry

Representation views include interpretable feature spaces and AntigenLM-derived embeddings from the parent thesis repository. The deployed app shows safe reduced-coordinate views:

- PCA coordinates;
- t-SNE coordinates;
- hashed point identifiers;
- minimal coarse metadata.

These maps support visual inspection of neighborhoods, coarse organization and outliers. They are descriptive representation geometry, not biological validation.

## 5. Structure context

The structure layer uses public RCSB/PDB entries for HA/NA context:

- `3LZG`;
- `3VUN`;
- `3NSS`;
- `6BR6`.

The app reports alignment QC, chain/coverage summaries and mapping status. It does not claim residue-level biological effects unless those mappings are independently validated.

## 6. Safe export

The public app is built from a safe derived-data layer:

- aggregate summaries;
- binned temporal/geographic summaries;
- reduced PCA/t-SNE coordinates;
- short token summaries;
- hash-based visual identifiers;
- public PDB identifiers;
- alignment-QC summaries;
- explanatory guide chunks.

It excludes raw sequences, FASTA files, accessions, isolate names, source sequence hashes, unrestricted sample-level metadata, long tokens and any table that could reasonably reconstruct restricted records.
