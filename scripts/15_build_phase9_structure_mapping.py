#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

from flugenome3d.structure_mapping import build_structure_mapping_tables


PROJECT = Path(__file__).resolve().parents[1]


STRUCTURES = [
    {"pdb_id": "3LZG", "label": "HA H1N1 reference structure", "protein": "HA", "subtype_context": "H1N1"},
    {"pdb_id": "3VUN", "label": "HA H3N2 reference structure", "protein": "HA", "subtype_context": "H3N2"},
    {"pdb_id": "3NSS", "label": "NA N1 reference structure", "protein": "NA", "subtype_context": "H1N1"},
    {"pdb_id": "6BR6", "label": "NA N2 reference structure", "protein": "NA", "subtype_context": "H3N2"},
]


def markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return ""
    columns = list(df.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in df.to_dict(orient="records"):
        lines.append("| " + " | ".join(str(row.get(column, "")) for column in columns) + " |")
    return "\n".join(lines)


def read_config(path: str | Path) -> dict:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def signal_catalog(signals: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (subtype, protein), group in signals.groupby(["subtype", "protein"]):
        mapped = mapping[(mapping["subtype"] == subtype) & (mapping["protein"] == protein)]
        rows.append(
            {
                "subtype": subtype,
                "protein": protein,
                "group": f"{protein}-{subtype}",
                "n_local_positions": group["local_residue_index"].nunique(),
                "n_mapped_positions": mapped["local_residue_index"].nunique() if not mapped.empty else 0,
                "mean_gc_fraction_codon": group["gc_fraction_codon"].mean(),
                "mean_cpg_codon_fraction": group["cpg_codon_fraction"].mean(),
                "mean_upa_codon_fraction": group["upa_codon_fraction"].mean(),
                "mean_aa_entropy": group["aa_entropy"].mean(),
                "max_aa_entropy": group["aa_entropy"].max(),
            }
        )
    return pd.DataFrame(rows)


def write_report(qc: pd.DataFrame, signals: pd.DataFrame, catalog: pd.DataFrame) -> None:
    out = PROJECT / "reports" / "phase9_structure_mapping_report.md"
    lines = [
        "# Phase 9 structure mapping QC report",
        "",
        "Phase 9 creates the first explicit bridge from refined HA/NA CDS positions to public PDB polymer sequences. It is an alignment-QC and residue-signal layer, not a validated antigenic or functional residue map.",
        "",
        "## What was built",
        "",
        "- Public RCSB FASTA sequences were loaded for 3LZG, 3VUN, 3NSS and 6BR6.",
        "- Refined CDS sequences were translated locally, without exporting sequences.",
        "- Per-residue aggregate signals were computed: codon GC fraction, CpG codon fraction, UpA/TA codon fraction and amino-acid entropy.",
        "- Local consensus amino-acid sequences were aligned to each public PDB polymer sequence.",
        "",
        "## Mapping QC",
        "",
        markdown_table(qc) if not qc.empty else "No mapping QC rows available.",
        "",
        "## Residue-signal catalog",
        "",
        markdown_table(catalog) if not catalog.empty else "No residue signal catalog available.",
        "",
        "## Interpretation",
        "",
        "The mapping is now more than pending: alignment QC is available. However, residue coloring in the 3D viewer still requires an additional chain/residue-number validation layer before FluGenome3D metrics are painted onto atoms.",
        "",
        "## Boundaries",
        "",
        "- These mapped signals are descriptive sequence-context summaries.",
        "- They are not antigenic sites, vaccine markers, escape sites, pathogenicity markers, fitness estimates or causal explanations.",
        "- No raw sequences, FASTA records, accessions or isolate names are exported.",
    ]
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_figure(qc: pd.DataFrame, catalog: pd.DataFrame) -> None:
    fig_dir = PROJECT / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    if not qc.empty:
        best = qc.sort_values(["pdb_id", "mapped_residues"], ascending=[True, False]).drop_duplicates("pdb_id")
        plt.figure(figsize=(8, 4.5), facecolor="#080807")
        ax = plt.gca()
        ax.set_facecolor("#01050a")
        ax.bar(best["pdb_id"], best["coverage_pdb"], color="#5cdce2")
        ax.set_ylim(0, 1.05)
        ax.set_title("Best public PDB polymer coverage", color="#edf7f4")
        ax.set_ylabel("PDB sequence coverage", color="#edf7f4")
        ax.tick_params(colors="#9fb3ae")
        plt.tight_layout()
        plt.savefig(fig_dir / "fig33_structure_mapping_qc.png", dpi=220)
        plt.close()

    if not catalog.empty:
        plt.figure(figsize=(8, 4.5), facecolor="#080807")
        ax = plt.gca()
        ax.set_facecolor("#01050a")
        x = catalog["group"]
        ax.plot(x, catalog["mean_cpg_codon_fraction"], marker="o", color="#79d99c", label="CpG codon fraction")
        ax.plot(x, catalog["mean_upa_codon_fraction"], marker="o", color="#5cdce2", label="UpA/TA codon fraction")
        ax.set_title("Residue signal catalog", color="#edf7f4")
        ax.set_ylabel("Mean fraction", color="#edf7f4")
        ax.tick_params(colors="#9fb3ae", axis="x", rotation=20)
        ax.legend(facecolor="#071017", edgecolor="#1b5561", labelcolor="#edf7f4")
        plt.tight_layout()
        plt.savefig(fig_dir / "fig34_structure_signal_catalog.png", dpi=220)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Build Phase 9 structure mapping QC artifacts.")
    parser.add_argument("--config", default=str(PROJECT / "config" / "phase7_9.yml"))
    args = parser.parse_args()
    cfg = read_config(args.config).get("phase9_structure_mapping", {})
    url_template = cfg.get("rcsb_fasta_url_template", "https://www.rcsb.org/fasta/entry/{pdb_id}/display")

    panel = pd.read_parquet(PROJECT / "data" / "processed" / "panels" / "mvp_cds_refined_panel.parquet")
    qc, signals, mapping = build_structure_mapping_tables(panel, STRUCTURES, url_template=url_template)

    local_dir = PROJECT / "data" / "processed" / "structure_mapping"
    local_dir.mkdir(parents=True, exist_ok=True)
    signals.to_parquet(local_dir / "phase9_residue_signal_summary.parquet", index=False)
    mapping.to_parquet(local_dir / "phase9_structure_mapping_table.parquet", index=False)

    tables = PROJECT / "results" / "tables"
    qc.to_csv(tables / "phase9_structure_mapping_qc.csv", index=False)
    catalog = signal_catalog(signals, mapping)
    catalog.to_csv(tables / "phase9_structure_signal_catalog.csv", index=False)

    write_report(qc, signals, catalog)
    make_figure(qc, catalog)
    print(f"Phase 9 complete: {len(qc)} alignment rows, {len(mapping):,} mapped residue rows.")


if __name__ == "__main__":
    main()
