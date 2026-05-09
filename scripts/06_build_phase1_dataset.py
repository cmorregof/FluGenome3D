#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from flugenome3d.dataset_builder import (
    build_pair_table,
    build_phase1_panels,
    write_aggregate_tables,
    write_panels,
    write_phase1_figures,
)
from flugenome3d.local_loader import (
    detect_relevant_files,
    load_dedup_metadata,
    load_local_file_set,
    load_paired_datasets,
    load_rich_metadata,
)
from flugenome3d.utils import load_yaml


def _markdown_table(df: pd.DataFrame, max_rows: int = 20) -> str:
    shown = df.head(max_rows).copy()
    if shown.empty:
        return "_No rows._"
    cols = list(shown.columns)
    lines = [
        "| " + " | ".join(cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for _, row in shown.iterrows():
        values = []
        for col in cols:
            value = row[col]
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    if len(df) > max_rows:
        lines.append(f"\n_Table truncated to {max_rows} of {len(df)} rows._")
    return "\n".join(lines)


def write_phase1_report(
    report_path: str | Path,
    detected: dict,
    table_paths: dict[str, Path],
    panel_paths: dict[str, Path],
    figure_paths: dict[str, Path],
    full_choice: str,
) -> None:
    dataset_summary = pd.read_csv(table_paths["dataset_summary"], keep_default_na=False)
    panel_summary = pd.read_csv(table_paths["panel_summary"], keep_default_na=False)
    missingness = pd.read_csv(table_paths["missingness_summary"], keep_default_na=False)
    duplicates = pd.read_csv(table_paths["duplicate_summary"], keep_default_na=False)
    length_qc = pd.read_csv(table_paths["length_qc_summary"], keep_default_na=False)
    year_counts = pd.read_csv(table_paths["year_subtype_counts"], keep_default_na=False)
    year_full = year_counts[year_counts["panel"] == "full_panel"]

    report = f"""# FluGenome3D Phase 1 dataset report

Generated from local restricted parent-repository data. This report contains only aggregate counts and no raw nucleotide sequences.

## Detected datasets

- FASTA files detected: {len(detected.get("fasta_files", []))}
- Paired HA/NA JSON datasets detected: {len(detected.get("paired_dataset_paths", []))}
- Rich metadata CSV detected: {detected.get("rich_metadata_path") is not None}
- Dedup metadata CSV detected: {detected.get("dedup_metadata_path") is not None}

## Columns used

- Paired JSON: `epi_isl`, `subtype`, `year`, `month`, `day`, `ha_sequence`, `na_sequence`.
- Rich metadata: `Isolate_Id`, `Subtype`, `Host`, `Location`, `Collection_Date`, `Clade`.
- Dedup metadata: `epi_isl`, `subtype`, `year`, `month`, `matched`, `clade`, `major_clade`.

## Filtering and labels

- Subtypes retained: H1N1 and H3N2.
- Segment labels: paired HA and paired NA.
- Date availability: year, month, and day required for `paired_pass_qc`.
- HA expected length: 1650-1800 nt.
- NA expected length: 1200-1600 nt.
- Maximum ambiguous fraction: 0.01.
- Alphabet filter: A/C/G/T/N only.
- Duplicate definition: exact SHA256 identity of the HA+NA sequence-hash pair.

## Pair counts by stage

{_markdown_table(dataset_summary)}

## Panel definitions

- `smoke_panel`: small deduplicated local panel for fast tests.
- `mvp_panel`: 10,000-pair deduplicated panel, balanced by subtype and spread across years where possible, preferring rich human metadata.
- `full_panel`: {full_choice}.

{_markdown_table(panel_summary)}

## Duplicate accounting

{_markdown_table(duplicates)}

## Metadata missingness

{_markdown_table(missingness)}

## Year distribution

The full panel spans {int(year_full["year"].min()) if not year_full.empty else 0}-{int(year_full["year"].max()) if not year_full.empty else 0}. Counts by year/subtype are stored in `results/tables/phase1_year_subtype_counts.csv`.

## Length and ambiguity QC

{_markdown_table(length_qc)}

## Local-only outputs

- `{panel_paths["smoke"]}`
- `{panel_paths["mvp"]}`
- `{panel_paths["full"]}`

These Parquet files may contain raw sequences and accession-level fields. They must remain gitignored.

## GitHub-safe aggregate outputs

- `results/tables/phase1_dataset_summary.csv`
- `results/tables/phase1_panel_summary.csv`
- `results/tables/phase1_missingness_summary.csv`
- `results/tables/phase1_year_subtype_counts.csv`
- `results/tables/phase1_length_qc_summary.csv`
- `results/tables/phase1_duplicate_summary.csv`
- `{figure_paths["overview"]}`
- `{figure_paths["year_subtype"]}`
- `{figure_paths["length_qc"]}`

## Limitations

- This phase constructs descriptive local panels only.
- Host is confirmed from the rich metadata subset; not every valid paired record has rich metadata coverage.
- Exact deduplication removes identical HA+NA pairs but does not collapse near-duplicates.
- No GROVER tokenizer, 3D structure mapping, predictive modeling, antigenicity, vaccine, escape, fitness, or optimization analysis is implemented in Phase 1.
"""
    path = Path(report_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build FluGenome3D Phase 1 local analytic panels and aggregate reports.")
    parser.add_argument("--config", default="config/local_paths.yml")
    parser.add_argument("--filters", default="config/filters.yml")
    parser.add_argument("--panel-outdir", default="data/processed/panels")
    parser.add_argument("--tables-outdir", default="results/tables")
    parser.add_argument("--figures-outdir", default="results/figures")
    parser.add_argument("--report", default="reports/phase1_dataset_report.md")
    args = parser.parse_args()

    filters = load_yaml(args.filters)
    file_set = load_local_file_set(args.config)
    detected = detect_relevant_files(file_set)

    paired_paths = detected["paired_dataset_paths"]
    rich_path = detected["rich_metadata_path"]
    dedup_path = detected["dedup_metadata_path"]
    if not paired_paths:
        raise SystemExit("No paired HA/NA JSON datasets were detected.")
    if rich_path is None:
        raise SystemExit("No rich metadata CSV was detected.")
    if dedup_path is None:
        raise SystemExit("No dedup metadata CSV was detected.")

    paired = load_paired_datasets(paired_paths)
    rich = load_rich_metadata(rich_path)
    dedup = load_dedup_metadata(dedup_path)
    pair_table = build_pair_table(paired, rich, dedup, filters)
    panels = build_phase1_panels(pair_table, filters)

    panel_paths = write_panels(panels, args.panel_outdir)
    table_paths = write_aggregate_tables(pair_table, panels, args.tables_outdir)
    figure_paths = write_phase1_figures(table_paths, panels, args.figures_outdir)

    write_phase1_report(
        args.report,
        detected,
        table_paths,
        panel_paths,
        figure_paths,
        full_choice="all 82,306 exact HA+NA deduplicated pairs passing Phase 1 QC",
    )

    summary = pd.read_csv(table_paths["dataset_summary"])
    print(summary.to_string(index=False))
    print(f"Wrote panels: {', '.join(str(p) for p in panel_paths.values())}")
    print(f"Wrote report: {args.report}")


if __name__ == "__main__":
    main()
