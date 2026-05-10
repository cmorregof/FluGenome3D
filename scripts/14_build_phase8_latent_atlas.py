#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "results" / "tables"
REPORTS = PROJECT / "reports"
FIGURES = PROJECT / "results" / "figures"


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


def read_csv(name: str) -> pd.DataFrame:
    path = TABLES / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def build_comparison() -> pd.DataFrame:
    phase4 = read_csv("phase4_representation_summary.csv")
    phase4_sil = read_csv("phase4_silhouette_scores.csv")
    phase7_cache = read_csv("phase7_antigenlm_cache_summary.csv")
    phase7_pca = read_csv("phase7_antigenlm_pca_summary.csv")
    phase7_spearman = read_csv("phase7_antigenlm_spearman_summary.csv")
    rows = []
    if not phase4.empty:
        sil = phase4_sil[phase4_sil["label_type"] == "protein_subtype"] if not phase4_sil.empty else pd.DataFrame()
        for row in phase4.itertuples(index=False):
            match = sil[sil["representation"] == row.representation] if not sil.empty else pd.DataFrame()
            rows.append(
                {
                    "family": "deterministic_baseline",
                    "representation": row.representation,
                    "n_sequences": row.n_sequences,
                    "n_features": row.n_features,
                    "pca_2pc_variance": row.pca_explained_variance_total_2pc,
                    "protein_subtype_silhouette": match["silhouette"].iloc[0] if not match.empty else None,
                    "molecular_proxy": "not_computed_here",
                    "app_role": "baseline comparison",
                }
            )
    if not phase7_cache.empty:
        pca_global = phase7_pca[phase7_pca["group"] == "global"] if not phase7_pca.empty else pd.DataFrame()
        spearman = phase7_spearman[phase7_spearman["metric"] == "hamming_ha_na"] if not phase7_spearman.empty else pd.DataFrame()
        rows.append(
            {
                "family": "learned_latent",
                "representation": "AntigenLM HA+NA embeddings",
                "n_sequences": phase7_cache["n_records"].iloc[0],
                "n_features": phase7_cache["embedding_dim"].iloc[0],
                "pca_2pc_variance": None,
                "protein_subtype_silhouette": None,
                "molecular_proxy": f"HA+NA Hamming rho mean {spearman['rho_mean'].mean():.3f}" if not spearman.empty else "available in parent metrics",
                "app_role": "learned biological representation layer",
                "global_pca_n95": pca_global["n95"].iloc[0] if not pca_global.empty else None,
                "global_participation_ratio": pca_global["participation_ratio"].iloc[0] if not pca_global.empty else None,
            }
        )
    return pd.DataFrame(rows)


def write_report(comparison: pd.DataFrame) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Phase 8 Latent Atlas report",
        "",
        "Phase 8 turns Phase 7 outputs into a visual atlas concept for the app: deterministic nucleotide/codon representations and the learned AntigenLM layer are shown as complementary views of the same HA/NA research landscape.",
        "",
        "## Representation families",
        "",
        markdown_table(comparison) if not comparison.empty else "No comparison table available.",
        "",
        "## Design decision",
        "",
        "The AntigenLM layer should be foregrounded as a learned representation, while k-mer, codon-frequency and RSCU spaces remain transparent baselines. This avoids treating learned embeddings as magic and gives interviewers a clear ladder from interpretable biology to learned geometry.",
        "",
        "## App behavior",
        "",
        "- The Representation Projector remains the baseline feature-space explorer.",
        "- The new Latent Atlas view focuses on AntigenLM PCA points, molecular-geometry summaries, clade enrichment and temporal locality.",
        "- All displayed point IDs are hash-based internal IDs with minimal metadata.",
        "",
        "## Boundaries",
        "",
        "This phase does not add GROVER, BPE, prediction, antigenicity, escape, vaccine, pathogenicity, fitness or sequence optimization claims.",
    ]
    (REPORTS / "phase8_latent_atlas_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_figure(comparison: pd.DataFrame) -> None:
    if comparison.empty:
        return
    FIGURES.mkdir(parents=True, exist_ok=True)
    plot = comparison.copy()
    plot["n_sequences"] = pd.to_numeric(plot["n_sequences"], errors="coerce")
    plt.figure(figsize=(8, 4.8), facecolor="#080807")
    ax = plt.gca()
    ax.set_facecolor("#01050a")
    colors = ["#5cdce2" if family == "learned_latent" else "#79d99c" for family in plot["family"]]
    ax.barh(plot["representation"], plot["n_sequences"], color=colors)
    ax.set_title("Representation atlas coverage", color="#edf7f4")
    ax.set_xlabel("Sequences / records represented", color="#edf7f4")
    ax.tick_params(colors="#9fb3ae")
    for label in ax.get_yticklabels():
        label.set_color("#edf7f4")
    plt.tight_layout()
    plt.savefig(FIGURES / "fig32_representation_atlas_coverage.png", dpi=220)
    plt.close()


def main() -> None:
    comparison = build_comparison()
    comparison.to_csv(TABLES / "phase8_representation_family_comparison.csv", index=False)
    write_report(comparison)
    make_figure(comparison)
    print(f"Phase 8 complete: {len(comparison)} representation rows summarized.")


if __name__ == "__main__":
    main()
