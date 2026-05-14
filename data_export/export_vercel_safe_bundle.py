#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import svds


PROJECT = Path(__file__).resolve().parents[1]
APP_DATA = PROJECT / "app" / "data"
SAFE_SALT = "flugenome3d-derived-data-v1"
LONG_SEQUENCE_RE = re.compile(r"[ACGTN]{80,}")
GROUP_ORDER = ["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"]
ATLAS_PANELS = {
    "mvp_panel": "Balanced MVP",
    "full_panel": "Full deduplicated",
}
GUIDE_SOURCE_FILES = [
    "docs/methods.md",
    "docs/biological_background.md",
    "docs/data_governance.md",
    "docs/reproducibility.md",
    "docs/claims_and_limitations.md",
    "docs/visual_guide.md",
    "reports/phase1_dataset_report.md",
    "reports/phase2_sequence_context_report.md",
    "reports/phase3_cds_refinement_report.md",
    "reports/phase4_representation_audit_report.md",
    "reports/phase5_tokenization_audit_report.md",
    "reports/phase6_tokenization_stability_report.md",
    "reports/phase7_antigenlm_bridge_report.md",
    "reports/phase8_latent_atlas_report.md",
    "reports/phase9_structure_mapping_report.md",
]
GUIDE_FORMULA_CARDS = [
    {
        "id": "formula_gc_fraction",
        "name": "GC fraction",
        "formula": "(G + C) / sequence length",
        "plain_language": "The share of positions that are G or C. It is a compact way to compare nucleotide composition across HA/NA groups.",
        "used_in": ["Sequence/Token Inspector", "Project Guide"],
        "data_requirement": "raw nucleotide context; no coding frame required",
        "claim_boundary": "composition summary only",
    },
    {
        "id": "formula_cpg_oe",
        "name": "CpG observed/expected",
        "formula": "f(CG) / (f(C) x f(G))",
        "plain_language": "Compares how often CG appears with how often it would be expected from C and G abundance.",
        "used_in": ["Sequence/Token Inspector"],
        "data_requirement": "raw nucleotide context; no coding frame required",
        "claim_boundary": "descriptive dinucleotide feature, not a causal statement",
    },
    {
        "id": "formula_upa_oe",
        "name": "UpA observed/expected",
        "formula": "f(TA) / (f(T) x f(A))",
        "plain_language": "Uses DNA TA as the proxy for RNA UpA and compares observed TA with single-base expectation.",
        "used_in": ["Sequence/Token Inspector"],
        "data_requirement": "raw nucleotide context; no coding frame required",
        "claim_boundary": "descriptive dinucleotide feature, not pathogenicity or antigenicity inference",
    },
    {
        "id": "formula_token_entropy",
        "name": "Token entropy",
        "formula": "-sum p(token) log2 p(token)",
        "plain_language": "Summarizes how spread out a tokenizer vocabulary is. Higher entropy means token usage is less concentrated.",
        "used_in": ["Sequence/Token Inspector", "Ask FluGenome3D"],
        "data_requirement": "token counts by tokenizer and group",
        "claim_boundary": "representation diagnostic only",
    },
    {
        "id": "formula_js_distance",
        "name": "Jensen-Shannon distance",
        "formula": "sqrt(JS divergence between token distributions)",
        "plain_language": "A bounded distance between two group-level token distributions. Larger values mean the token distributions are more different.",
        "used_in": ["Sequence/Token Inspector", "Tokenization stability"],
        "data_requirement": "aggregate token frequency distributions",
        "claim_boundary": "group comparison only; not prediction",
    },
    {
        "id": "formula_rscu",
        "name": "RSCU",
        "formula": "codon count / mean synonymous-codon count",
        "plain_language": "Compares each codon with the average usage of codons that encode the same amino acid.",
        "used_in": ["Project Guide", "CDS/codon summaries"],
        "data_requirement": "refined CDS panel with translation/frame QC",
        "claim_boundary": "codon-usage summary only",
    },
    {
        "id": "formula_pca",
        "name": "PCA coordinates",
        "formula": "orthogonal axes capturing largest variance directions",
        "plain_language": "A projection that turns high-dimensional features into 2D or 3D coordinates for visual inspection.",
        "used_in": ["Representation Projector", "AntigenLM Latent Atlas"],
        "data_requirement": "feature matrix or learned embeddings",
        "claim_boundary": "visual geometry, not biological validation",
    },
    {
        "id": "formula_tsne",
        "name": "t-SNE coordinates",
        "formula": "nonlinear neighbor-preserving projection",
        "plain_language": "A map that emphasizes local neighborhoods in high-dimensional embeddings, useful for visual cluster inspection.",
        "used_in": ["AntigenLM Latent Atlas"],
        "data_requirement": "sampled learned embeddings; exported only as hash-based coordinates",
        "claim_boundary": "visual neighborhood audit, not biological proof",
    },
]
GUIDE_GLOSSARY = [
    {
        "term": "HA",
        "short_definition": "Hemagglutinin, an Influenza A surface glycoprotein and one half of the HA/NA pair audited here.",
        "view": "Dataset Atlas",
        "ask": "What is HA in the FluGenome3D dataset?",
    },
    {
        "term": "NA",
        "short_definition": "Neuraminidase, the second surface glycoprotein paired with HA in this project.",
        "view": "Dataset Atlas",
        "ask": "What is NA in the FluGenome3D dataset?",
    },
    {
        "term": "MVP panel",
        "short_definition": "The balanced analysis panel used for the main descriptive figures and representation views.",
        "view": "Dataset Atlas",
        "ask": "Why does FluGenome3D use an MVP panel?",
    },
    {
        "term": "Refined CDS",
        "short_definition": "A subset of sequences with explicit coding-frame rescue or strict translation QC for codon-level analyses.",
        "view": "Sequence/Token Inspector",
        "ask": "Why is codon usage separated from raw sequence context?",
    },
    {
        "term": "Tokenizer",
        "short_definition": "A rule for splitting sequence text into small units such as codons or k-mers.",
        "view": "Sequence/Token Inspector",
        "ask": "What is a tokenizer in FluGenome3D?",
    },
    {
        "term": "AntigenLM embedding",
        "short_definition": "A learned vector representation from the parent thesis repository, displayed here only as derived coordinates.",
        "view": "AntigenLM Latent Atlas",
        "ask": "How should I read the AntigenLM latent atlas?",
    },
    {
        "term": "t-SNE",
        "short_definition": "A nonlinear projection that helps inspect local neighborhoods in the AntigenLM latent space.",
        "view": "AntigenLM Latent Atlas",
        "ask": "How should I compare PCA and t-SNE in the AntigenLM Latent Atlas?",
    },
    {
        "term": "PDB structure",
        "short_definition": "A public molecular structure record loaded from RCSB for visual inspection.",
        "view": "3D Molecular Viewer",
        "ask": "What does the structure viewer load from RCSB?",
    },
    {
        "term": "Mapping pending",
        "short_definition": "Alignment QC exists, but residue-level metric coloring waits for chain and residue-number validation.",
        "view": "3D Molecular Viewer",
        "ask": "What does structure mapping still need before residue coloring?",
    },
]
GUIDE_VIEW_PROMPTS = {
    "atlas": [
        {"label": "Why country-level only?", "question": "Why does the Dataset Atlas show country-level aggregates instead of sample-level locations?"},
        {"label": "MVP vs full panel", "question": "How should I interpret the MVP panel versus the full deduplicated panel?"},
    ],
    "latent": [
        {"label": "Read latent geometry", "question": "How should I read the AntigenLM latent atlas?"},
        {"label": "PCA vs t-SNE", "question": "How should I compare PCA and t-SNE in the AntigenLM Latent Atlas?"},
        {"label": "What does rho mean?", "question": "What does the HA+NA molecular rho summary mean in the latent atlas?"},
    ],
    "projector": [
        {"label": "PCA axes", "question": "What do PCA axes mean in the Representation Projector?"},
        {"label": "Compare representations", "question": "How do raw k-mer, codon frequency and RSCU representations differ?"},
    ],
    "inspector": [
        {"label": "CpG / UpA", "question": "What do CpG O/E and UpA O/E mean in FluGenome3D?"},
        {"label": "Token entropy", "question": "How should I interpret token entropy and effective vocabulary?"},
    ],
    "structure": [
        {"label": "Mapping status", "question": "What does structure mapping still need before residue coloring?"},
        {"label": "RCSB structures", "question": "What does the structure viewer load from RCSB and what does it not claim?"},
    ],
    "bridge": [
        {"label": "Bridge logic", "question": "How does Bridge View connect sequence context, representation space and structure?"},
        {"label": "Claim boundary", "question": "What can FluGenome3D say and not say from the Bridge View?"},
    ],
}


def read_csv(path: str) -> pd.DataFrame:
    full = PROJECT / path
    if not full.exists():
        return pd.DataFrame()
    return pd.read_csv(full, keep_default_na=False)


def read_parquet(path: str) -> pd.DataFrame:
    full = PROJECT / path
    if not full.exists():
        return pd.DataFrame()
    return pd.read_parquet(full)


def read_json(path: str) -> dict[str, Any]:
    full = PROJECT / path
    if not full.exists():
        return {}
    return json.loads(full.read_text(encoding="utf-8"))


def safe_id(value: str) -> str:
    digest = hashlib.sha256(f"{SAFE_SALT}::{value}".encode("utf-8")).hexdigest()
    return f"pt_{digest[:16]}"


def year_bin(year: Any) -> str:
    try:
        value = int(year)
    except (TypeError, ValueError):
        return "unknown"
    if value <= 2008:
        return "pre-2009"
    if value <= 2014:
        return "2009-2014"
    if value <= 2019:
        return "2015-2019"
    return "2020+"


def finite_number(value: Any, digits: int = 6) -> float | int | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(number):
        return None
    if number.is_integer():
        return int(number)
    return round(number, digits)


def records(df: pd.DataFrame, max_rows: int | None = None) -> list[dict[str, Any]]:
    if df.empty:
        return []
    use = df.head(max_rows).copy() if max_rows else df.copy()
    out = []
    for row in use.to_dict(orient="records"):
        clean = {}
        for key, value in row.items():
            if pd.isna(value):
                clean[key] = None
            elif isinstance(value, (np.integer,)):
                clean[key] = int(value)
            elif isinstance(value, (np.floating, float)):
                clean[key] = finite_number(value)
            else:
                clean[key] = value
        out.append(clean)
    return out


def build_geographic_atlas() -> dict[str, Any]:
    panel_frames: list[pd.DataFrame] = []
    for panel_id, panel_label in ATLAS_PANELS.items():
        df = read_parquet(f"data/processed/panels/{panel_id}.parquet")
        if df.empty:
            continue
        required = ["country", "region", "subtype", "year"]
        missing = [column for column in required if column not in df.columns]
        if missing:
            continue
        use = df[required].copy()
        use["panel"] = panel_id
        use["panel_label"] = panel_label
        use["country"] = use["country"].fillna("").astype(str).str.strip()
        use["region"] = use["region"].fillna("unknown").astype(str).str.strip()
        missing_country = {"", "nan", "none", "unknown"}
        use = use[~use["country"].str.lower().isin(missing_country)]
        panel_frames.append(use)

    if not panel_frames:
        return {
            "policy": "aggregate country/region counts only; no accessions, isolate names, sequence strings, or sample-level rows",
            "panels": [],
            "country_subtype_counts": [],
            "country_totals": [],
            "region_summary": [],
            "year_bin_country_counts": [],
        }

    atlas = pd.concat(panel_frames, ignore_index=True)
    atlas["year_bin"] = atlas["year"].map(year_bin)

    country_region = (
        atlas.groupby(["panel", "panel_label", "country", "region"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["panel", "country", "n"], ascending=[True, True, False])
        .drop_duplicates(["panel", "panel_label", "country"])
        [["panel", "panel_label", "country", "region"]]
    )

    country_subtype = (
        atlas.groupby(["panel", "panel_label", "country", "subtype"], dropna=False)
        .size()
        .reset_index(name="n_pairs")
        .merge(country_region, on=["panel", "panel_label", "country"], how="left")
        .sort_values(["panel", "country", "subtype"])
    )
    country_totals = (
        atlas.groupby(["panel", "panel_label", "country"], dropna=False)
        .agg(n_pairs=("subtype", "size"), year_min=("year", "min"), year_max=("year", "max"), n_subtypes=("subtype", "nunique"))
        .reset_index()
        .merge(country_region, on=["panel", "panel_label", "country"], how="left")
        .sort_values(["panel", "n_pairs"], ascending=[True, False])
    )
    region_summary = (
        atlas.groupby(["panel", "panel_label", "region", "subtype"], dropna=False)
        .size()
        .reset_index(name="n_pairs")
        .sort_values(["panel", "region", "subtype"])
    )
    year_country = (
        atlas.groupby(["panel", "panel_label", "year_bin", "country", "subtype"], dropna=False)
        .size()
        .reset_index(name="n_pairs")
        .sort_values(["panel", "year_bin", "country", "subtype"])
    )

    return {
        "policy": "aggregate country/region counts only; no accessions, isolate names, sequence strings, exact collection locations, or sample-level rows",
        "panels": [{"panel": panel, "label": label} for panel, label in ATLAS_PANELS.items()],
        "country_subtype_counts": records(country_subtype),
        "country_totals": records(country_totals),
        "region_summary": records(region_summary),
        "year_bin_country_counts": records(year_country),
    }


def stratified_points(df: pd.DataFrame, max_points: int, seed: int = 42) -> pd.DataFrame:
    if len(df) <= max_points:
        return df.copy()
    rng = np.random.default_rng(seed)
    selected: list[int] = []
    groups = sorted(df["protein_subtype"].dropna().unique())
    per_group = max(1, max_points // max(1, len(groups)))
    for group in groups:
        idx = df.index[df["protein_subtype"] == group].to_numpy()
        take = min(per_group, len(idx))
        selected.extend(rng.choice(idx, size=take, replace=False).tolist())
    if len(selected) < max_points:
        remaining = np.setdiff1d(df.index.to_numpy(), np.asarray(selected), assume_unique=False)
        take = min(max_points - len(selected), len(remaining))
        if take:
            selected.extend(rng.choice(remaining, size=take, replace=False).tolist())
    return df.loc[sorted(selected)].copy()


def run_sparse_pca(matrix: sparse.spmatrix, n_components: int = 3) -> tuple[np.ndarray, list[float | None]]:
    x = matrix.tocsr().astype(np.float32)
    n_samples, n_features = x.shape
    if n_samples == 0:
        return np.empty((0, n_components), dtype=np.float32), [None] * n_components

    dense_mean = np.asarray(x.mean(axis=0)).ravel().astype(np.float32)
    centered = x.toarray().astype(np.float32, copy=False) - dense_mean
    max_components = max(1, min(n_components, n_samples - 1, n_features - 1 if n_features > 1 else 1))
    if max_components == 1:
        u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
        u = u[:, :1]
        singular_values = singular_values[:1]
    else:
        u, singular_values, _ = svds(centered, k=max_components)
        order = np.argsort(singular_values)[::-1]
        singular_values = singular_values[order]
        u = u[:, order]

    embedding = u * singular_values
    total_var = float(np.sum(centered * centered) / max(n_samples - 1, 1))
    explained = (singular_values**2) / max(n_samples - 1, 1)
    ratios = (explained / total_var).tolist() if total_var else [None] * len(singular_values)
    clean_ratios = [finite_number(value, digits=6) for value in ratios]
    if embedding.shape[1] < n_components:
        embedding = np.hstack([embedding, np.zeros((n_samples, n_components - embedding.shape[1]), dtype=np.float32)])
        clean_ratios.extend([None] * (n_components - len(clean_ratios)))
    return embedding[:, :n_components].astype(np.float32), clean_ratios[:n_components]


def build_dataset_overview() -> dict[str, Any]:
    dataset = read_csv("results/tables/phase1_dataset_summary.csv")
    panels = read_csv("results/tables/phase1_panel_summary.csv")
    year_counts = read_csv("results/tables/phase1_year_subtype_counts.csv")
    duplicates = read_csv("results/tables/phase1_duplicate_summary.csv")
    translation = read_csv("results/tables/phase2_translation_qc_summary.csv")
    refined = read_csv("results/tables/phase3_refined_translation_qc_summary.csv")
    rescue = read_csv("results/tables/phase3_rescue_status_by_group.csv")

    mvp_year = year_counts[year_counts["panel"].isin(["mvp_panel", "mvp"])] if not year_counts.empty else pd.DataFrame()
    if not mvp_year.empty:
        mvp_year["year_bin"] = mvp_year["year"].map(year_bin)
        temporal = (
            mvp_year.groupby(["year_bin", "subtype"], dropna=False)["n_pairs"]
            .sum()
            .reset_index()
            .sort_values(["year_bin", "subtype"])
        )
    else:
        temporal = pd.DataFrame()

    return {
        "schema_version": "safe-bundle-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "data_mode": "deployable_derived",
        "source": "real derived artifacts from FluGenome3D phases 0-6",
        "not_included": ["raw sequences", "FASTA", "restricted Parquet panels", "accessions", "isolate names", "long tokens"],
        "dataset_summary": records(dataset),
        "panel_summary": records(panels),
        "temporal_counts": records(temporal),
        "duplicate_summary": records(duplicates),
        "naive_translation_qc": records(translation),
        "cds_refined_qc": records(refined),
        "rescue_status": records(rescue),
        "geographic_atlas": build_geographic_atlas(),
    }


def build_representation_maps(max_points: int = 12000) -> dict[str, Any]:
    summary = read_csv("results/tables/phase4_representation_summary.csv")
    silhouettes = read_csv("results/tables/phase4_silhouette_scores.csv")
    representations = [
        {
            "id": "raw_kmer4_tfidf_pca",
            "matrix": "data/processed/representations/mvp_kmer4_tfidf.npz",
            "metadata": "data/processed/representations/mvp_raw_representation_metadata.parquet",
            "source_representation": "kmer4_tfidf",
            "label": "Raw k-mer TF-IDF PCA (k=4)",
            "source": "raw_nucleotide",
            "cds_status": "not_required",
            "description": "Raw nucleotide representation; CDS frame is not required.",
        },
        {
            "id": "raw_kmer5_freq_pca",
            "matrix": "data/processed/representations/mvp_kmer5_freq.npz",
            "metadata": "data/processed/representations/mvp_raw_representation_metadata.parquet",
            "source_representation": "kmer5_freq",
            "label": "Raw k-mer frequency PCA (k=5)",
            "source": "raw_nucleotide",
            "cds_status": "not_required",
            "description": "Raw nucleotide frequency baseline with a longer k-mer context.",
        },
        {
            "id": "cds_codon_freq_pca",
            "matrix": "data/processed/representations/mvp_codon_freq.npz",
            "metadata": "data/processed/representations/mvp_cds_representation_metadata.parquet",
            "source_representation": "codon_freq",
            "label": "CDS refined codon-frequency PCA",
            "source": "cds_refined",
            "cds_status": "refined_cds",
            "description": "Codon-level representation restricted to the refined CDS panel.",
        },
        {
            "id": "cds_rscu_pca",
            "matrix": "data/processed/representations/mvp_rscu.npz",
            "metadata": "data/processed/representations/mvp_cds_representation_metadata.parquet",
            "source_representation": "rscu",
            "label": "CDS refined RSCU PCA",
            "source": "cds_refined",
            "cds_status": "refined_cds",
            "description": "RSCU representation restricted to the refined CDS panel.",
        },
    ]
    reps = []
    for spec in representations:
        matrix_path = PROJECT / spec["matrix"]
        metadata = read_parquet(spec["metadata"])
        if not matrix_path.exists() or metadata.empty:
            continue
        matrix = sparse.load_npz(matrix_path)
        if matrix.shape[0] != len(metadata):
            continue
        embedding, explained = run_sparse_pca(matrix, n_components=3)
        df = metadata.copy()
        df["axis1"] = embedding[:, 0]
        df["axis2"] = embedding[:, 1]
        df["axis3"] = embedding[:, 2]
        df["representation"] = spec["source_representation"]
        sampled = stratified_points(df, max_points=max_points)
        points = []
        for row in sampled.itertuples(index=False):
            points.append(
                [
                    safe_id(str(row.internal_sequence_id)),
                    finite_number(row.axis1, digits=7),
                    finite_number(row.axis2, digits=7),
                    finite_number(row.axis3, digits=7),
                    row.protein,
                    row.subtype,
                    row.protein_subtype,
                    year_bin(row.year),
                    spec["cds_status"],
                ]
            )
        source_rep = str(spec["source_representation"])
        rep_summary = summary[summary["representation"] == source_rep].to_dict(orient="records")
        sil = silhouettes[silhouettes["representation"] == source_rep].to_dict(orient="records")
        reps.append(
            {
                "id": spec["id"],
                "label": spec["label"],
                "source_representation": source_rep,
                "source": spec["source"],
                "description": spec["description"],
                "projection": "pca_3d",
                "axis_labels": ["PC1", "PC2", "PC3"],
                "pca_explained_variance": explained,
                "point_schema": ["id", "x", "y", "z", "protein", "subtype", "group", "year_bin", "cds_status"],
                "n_source_points": int(len(df)),
                "n_exported_points": int(len(points)),
                "privacy": "hashed point IDs; no accessions, isolate names, sequences, sequence hashes, or exact locations",
                "summary": rep_summary[0] if rep_summary else {},
                "silhouette_scores": sil,
                "points": points,
            }
        )
    return {
        "schema_version": "safe-bundle-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_policy": "real PCA coordinates exported with hashed IDs and minimal metadata; no sequences, accessions, isolate names, sequence hashes, or exact locations",
        "max_points_per_representation": max_points,
        "representations": reps,
    }


def build_metric_summaries() -> dict[str, Any]:
    gc = read_csv("results/tables/phase2_gc_cpg_upa_summary.csv")
    gc_mvp = gc[gc["panel"] == "mvp"] if not gc.empty else pd.DataFrame()
    dinuc = read_csv("results/tables/phase2_dinucleotide_odds_summary.csv")
    dinuc_mvp = dinuc[dinuc["panel"] == "mvp"] if not dinuc.empty else pd.DataFrame()
    kmer_entropy = read_csv("results/tables/phase2_kmer_entropy_summary.csv")
    translation = read_csv("results/tables/phase3_refined_translation_qc_summary.csv")
    rscu = read_csv("results/tables/phase3_refined_rscu_summary.csv")
    codon = read_csv("results/tables/phase3_refined_codon_usage_summary.csv")
    return {
        "schema_version": "safe-bundle-v1",
        "metrics_policy": "aggregate descriptive metrics only; no per-sequence values and no sequence strings",
        "gc_cpg_upa_summary": records(gc_mvp),
        "dinucleotide_odds_summary": records(dinuc_mvp),
        "kmer_entropy_summary": records(kmer_entropy),
        "cds_translation_qc_summary": records(translation),
        "rscu_summary": records(rscu, max_rows=600),
        "codon_usage_summary": records(codon, max_rows=600),
    }


def build_tokenization_summaries() -> dict[str, Any]:
    top = read_csv("results/tables/phase5_top_tokens_by_group.csv")
    if not top.empty:
        top = top[top["token"].astype(str).str.len() <= 6].copy()
    return {
        "schema_version": "safe-bundle-v1",
        "token_policy": "only aggregate token metrics and tokens of length <= 6 are exported",
        "tokenizer_summary": records(read_csv("results/tables/phase5_tokenizer_summary.csv")),
        "entropy_by_group": records(read_csv("results/tables/phase5_token_entropy_by_group.csv")),
        "effective_vocab_by_group": records(read_csv("results/tables/phase5_effective_vocab_by_group.csv")),
        "cpg_upa_token_summary": records(read_csv("results/tables/phase5_cpg_upa_token_summary.csv")),
        "boundary_crossing_summary": records(read_csv("results/tables/phase5_codon_boundary_crossing_summary.csv")),
        "group_js_distances": records(read_csv("results/tables/phase5_group_js_distances.csv")),
        "top_tokens_by_group": records(top),
    }


def build_stability_summaries() -> dict[str, Any]:
    return {
        "schema_version": "safe-bundle-v1",
        "bootstrap_policy": "aggregate bootstrap summaries only; no token lists per sequence",
        "bootstrap_metric_summary": records(read_csv("results/tables/phase6_bootstrap_metric_summary.csv")),
        "js_distance_stability": records(read_csv("results/tables/phase6_js_distance_stability.csv")),
        "top_token_jaccard_stability": records(read_csv("results/tables/phase6_top_token_jaccard_stability.csv")),
        "temporal_token_summary": records(read_csv("results/tables/phase6_temporal_token_summary.csv")),
        "tokenizer_robustness_ranking": records(read_csv("results/tables/phase6_tokenizer_robustness_ranking.csv")),
    }


def build_antigenlm_latent_atlas() -> dict[str, Any]:
    local = read_json("data/processed/antigenlm/antigenlm_latent_atlas.local.json")
    spearman = read_csv("results/tables/phase7_antigenlm_spearman_summary.csv")
    pca = read_csv("results/tables/phase7_antigenlm_pca_summary.csv")
    twonn = read_csv("results/tables/phase7_antigenlm_twonn_summary.csv")
    temporal = read_csv("results/tables/phase7_antigenlm_temporal_locality_summary.csv")
    clade = read_csv("results/tables/phase7_antigenlm_clade_enrichment_summary.csv")
    random_baseline = read_csv("results/tables/phase7_random_embedding_baseline_summary.csv")
    comparison = read_csv("results/tables/phase8_representation_family_comparison.csv")
    pca_projection = {
        "id": "antigenlm_full_pca",
        "label": "AntigenLM PCA 3D",
        "description": "Linear PCA coordinates derived from the parent AntigenLM embedding cache.",
        "projection": "pca_3d",
        "axis_labels": ["PC1", "PC2", "PC3"],
        "point_schema": local.get("point_schema", ["id", "x", "y", "z", "subtype", "year_bin", "representation", "source"]),
        "pca_explained_variance": local.get("pca_explained_variance", []),
        "n_source_points": local.get("n_source_points", 0),
        "n_exported_points": local.get("n_exported_points", 0),
        "sampling": local.get("sampling", {}),
        "privacy": "hash-based point IDs and coarse metadata only; no sequences, accessions, isolate names, sequence hashes, or checkpoint weights",
        "points": local.get("points", []),
    }

    return {
        "schema_version": "safe-bundle-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "coordinate_policy": "real AntigenLM PCA and t-SNE coordinates exported with hash-based IDs and minimal metadata; no sequences, accessions, isolate names, sequence hashes, or checkpoint weights",
        "model_card": {
            "name": "AntigenLM-derived HA+NA embedding cache",
            "role": "learned influenza representation layer from the parent thesis repository",
            "display_claim": "descriptive latent geometry audit only",
            "not_included": ["raw sequences", "source identifiers", "isolate names", "checkpoint weights", "sequence hashes"],
        },
        "projection": pca_projection,
        "additional_projections": local.get("additional_projections", []),
        "cache_summary": local.get("cache", {}),
        "spearman_summary": records(spearman),
        "pca_summary": records(pca),
        "twonn_summary": records(twonn, max_rows=80),
        "temporal_locality_summary": records(temporal, max_rows=80),
        "clade_enrichment_summary": records(clade),
        "random_embedding_baseline": records(random_baseline),
        "representation_family_comparison": records(comparison),
    }


def build_structure_mapping_export() -> dict[str, Any]:
    qc = read_csv("results/tables/phase9_structure_mapping_qc.csv")
    catalog = read_csv("results/tables/phase9_structure_signal_catalog.csv")
    mapping = read_parquet("data/processed/structure_mapping/phase9_structure_mapping_table.parquet")
    if not mapping.empty:
        keep = [
            "pdb_id",
            "protein",
            "subtype",
            "group",
            "pdb_entity",
            "chains",
            "local_residue_index",
            "pdb_sequence_index",
            "gc_fraction_codon",
            "cpg_codon_fraction",
            "upa_codon_fraction",
            "aa_entropy",
            "dominant_aa_fraction",
        ]
        mapping = mapping[[column for column in keep if column in mapping.columns]].copy()
        mapping = mapping.groupby("pdb_id", group_keys=False).head(900)
    return {
        "schema_version": "safe-bundle-v1",
        "mapping_policy": "alignment QC and aggregate residue signal tracks only; no sequences, no FASTA, no accessions, no isolate names",
        "mapping_status": "alignment_qc_available_residue_coloring_pending_chain_number_validation",
        "mapping_qc": records(qc),
        "signal_catalog": records(catalog),
        "mapped_tracks": records(mapping),
        "limitations": [
            "PDB polymer sequence indices are not final residue-number selections for 3D coloring.",
            "HA entries may contain separate HA1 and HA2 polymer entities.",
            "Metric-to-residue coloring requires chain/residue-number validation before public display.",
        ],
    }


def build_structure_catalog() -> dict[str, Any]:
    structures = [
        {
            "pdb_id": "3LZG",
            "label": "HA H1N1 reference structure",
            "protein": "HA",
            "subtype_context": "H1N1",
            "rcsb_url": "https://www.rcsb.org/structure/3LZG",
            "pdb_download_url": "https://files.rcsb.org/download/3LZG.pdb",
            "mapping_status": "alignment_qc_available",
        },
        {
            "pdb_id": "3VUN",
            "label": "HA H3N2 reference structure",
            "protein": "HA",
            "subtype_context": "H3N2",
            "rcsb_url": "https://www.rcsb.org/structure/3VUN",
            "pdb_download_url": "https://files.rcsb.org/download/3VUN.pdb",
            "mapping_status": "alignment_qc_available",
        },
        {
            "pdb_id": "3NSS",
            "label": "NA N1 reference structure",
            "protein": "NA",
            "subtype_context": "H1N1",
            "rcsb_url": "https://www.rcsb.org/structure/3NSS",
            "pdb_download_url": "https://files.rcsb.org/download/3NSS.pdb",
            "mapping_status": "alignment_qc_available",
        },
        {
            "pdb_id": "6BR6",
            "label": "NA N2 reference structure",
            "protein": "NA",
            "subtype_context": "H3N2",
            "rcsb_url": "https://www.rcsb.org/structure/6BR6",
            "pdb_download_url": "https://files.rcsb.org/download/6BR6.pdb",
            "mapping_status": "alignment_qc_available",
        },
    ]
    return {
        "schema_version": "safe-bundle-v1",
        "viewer_policy": "structures are loaded from public RCSB PDB IDs; sequence-to-structure alignment QC is available, while residue coloring remains pending chain-number validation",
        "structures": structures,
    }


def build_claims_and_limits() -> dict[str, Any]:
    return {
        "schema_version": "safe-bundle-v1",
        "banner": "GISAID-safe derived-data layer for thesis-derived Influenza A HA/NA artifacts. Hash-based IDs, aggregate views, and no raw sequence release.",
        "data_statement": "Real thesis-derived artifacts exported as aggregate summaries, reduced coordinates, short tokens and hash-based visual IDs. Raw sequences stay local.",
        "allowed_claims": [
            "FluGenome3D provides descriptive exploration of derived HA/NA artifacts.",
            "It audits sequence-context summaries, deterministic tokenization behavior, and representation geometry.",
            "PCA and t-SNE maps are safe reduced-coordinate views for inspection.",
            "Structure views use public PDB entries and conservative alignment-QC status.",
            "CDS-dependent summaries are limited to records passing documented CDS/QC filters.",
            "The deployed app is designed to avoid raw sequence, accession, or isolate-name redistribution.",
        ],
        "prohibited_claims": [
            "FluGenome3D predicts antigenic drift.",
            "FluGenome3D identifies immune escape mutations.",
            "FluGenome3D recommends vaccine candidates.",
            "FluGenome3D predicts fitness, pathogenicity, transmissibility, or selection.",
            "FluGenome3D validates AntigenLM, GROVER, BPE, or any tokenizer as biologically causal.",
            "FluGenome3D provides viral design or sequence optimization guidance.",
        ],
    }


def build_data_governance() -> dict[str, Any]:
    return {
        "schema_version": "safe-bundle-v1",
        "safe_exports": [
            "dataset_overview.safe.json",
            "representation_maps.safe.json",
            "metric_summaries.safe.json",
            "tokenization_summaries.safe.json",
            "stability_summaries.safe.json",
            "antigenlm_latent_atlas.safe.json",
            "structure_catalog.safe.json",
            "structure_mapping.safe.json",
            "lab_guide.safe.json",
            "claims_and_limits.safe.json",
            "data_governance.safe.json",
        ],
        "published_data_classes": [
            "aggregate summaries",
            "binned temporal/geographic summaries",
            "reduced PCA/t-SNE coordinates",
            "short token summaries",
            "hash-based visual identifiers",
            "public PDB identifiers",
            "alignment-QC summaries",
            "explanatory guide chunks",
        ],
        "excluded_data_classes": [
            "raw HA/NA sequences",
            "FASTA files",
            "accessions",
            "source record identifiers",
            "isolate names",
            "source sequence hashes",
            "unrestricted sample-level metadata",
            "long tokens or sequence-containing fields",
            "tables that could reasonably reconstruct restricted records",
        ],
        "local_full_mode": {
            "path": "app/data-local/",
            "gitignored": True,
            "description": "Optional local-only artifacts may be placed here for private development. They are not required for Vercel and must not be committed.",
        },
        "validation": {
            "long_sequence_regex": "[ACGTN]{80,}",
            "max_public_token_length": 6,
        },
    }


def clean_guide_text(text: str) -> str:
    text = LONG_SEQUENCE_RE.sub("[redacted long nucleotide-like string]", text)
    text = re.sub(r"`{3}.*?`{3}", " ", text, flags=re.DOTALL)
    text = re.sub(r"\s+", " ", text.replace("|", " ")).strip()
    return text


def topic_tags(title: str, text: str, source: str) -> list[str]:
    haystack = f"{title} {text} {source}".lower()
    mapping = {
        "dataset": ["dataset", "panel", "metadata", "atlas", "country", "coverage"],
        "sequence_context": ["gc", "cpg", "upa", "dinucleotide", "k-mer", "sequence-context"],
        "cds": ["cds", "codon", "rscu", "translation", "frame", "rescue"],
        "tokenization": ["token", "tokenizer", "entropy", "vocab", "jaccard", "bootstrap"],
        "representation": ["representation", "pca", "umap", "projector", "centroid", "silhouette"],
        "antigenlm": ["antigenlm", "latent", "embedding", "learned"],
        "structure": ["structure", "pdb", "rcsb", "alignment", "residue", "3lzg", "3vun", "3nss", "6br6"],
        "governance": ["governance", "safe", "hash", "raw sequence", "fasta", "restricted"],
        "claims": ["claim", "predict", "vaccine", "escape", "fitness", "pathogenicity", "limitation"],
        "future": ["grover", "bpe", "phase", "future", "recommendation"],
    }
    tags = [tag for tag, needles in mapping.items() if any(needle in haystack for needle in needles)]
    return tags or ["project"]


def chunk_markdown_source(source: str, max_chunks_per_source: int = 8) -> list[dict[str, Any]]:
    path = PROJECT / source
    if not path.exists():
        return []
    raw = path.read_text(encoding="utf-8", errors="replace")
    sections: list[tuple[str, list[str]]] = []
    current_title = path.stem.replace("_", " ").replace("-", " ").title()
    current_lines: list[str] = []
    for line in raw.splitlines():
        if line.lstrip().startswith("#"):
            if current_lines:
                sections.append((current_title, current_lines))
            current_title = line.lstrip("#").strip() or current_title
            current_lines = []
        else:
            current_lines.append(line)
    if current_lines:
        sections.append((current_title, current_lines))

    chunks: list[dict[str, Any]] = []
    for index, (title, lines) in enumerate(sections[:max_chunks_per_source], start=1):
        text = clean_guide_text(" ".join(lines))
        if len(text) < 80:
            continue
        chunks.append(
            {
                "id": safe_id(f"guide::{source}::{index}").replace("pt_", "guide_"),
                "title": title[:120],
                "source": source,
                "section": title[:120],
                "topic_tags": topic_tags(title, text, source),
                "text": text[:950],
            }
        )
    return chunks


def manual_guide_cards() -> list[dict[str, Any]]:
    formula_cards = [
        {
            "id": formula["id"].replace("formula_", "guide_formula_"),
            "title": formula["name"],
            "source": "FluGenome3D guide card",
            "section": "Plain-language formulas",
            "topic_tags": topic_tags(formula["name"], f"{formula['plain_language']} {formula['claim_boundary']}", "formula"),
            "text": f"{formula['name']}: {formula['plain_language']} Formula: {formula['formula']}. Boundary: {formula['claim_boundary']}.",
        }
        for formula in GUIDE_FORMULA_CARDS
    ]
    glossary_cards = [
        {
            "id": safe_id(f"guide::glossary::{item['term']}").replace("pt_", "guide_"),
            "title": item["term"],
            "source": "FluGenome3D glossary card",
            "section": "Glossary",
            "topic_tags": topic_tags(item["term"], item["short_definition"], item["view"]),
            "text": f"{item['term']}: {item['short_definition']} Related view: {item['view']}.",
        }
        for item in GUIDE_GLOSSARY
    ]
    extra_cards = [
        {
            "id": "guide_future_grover_bpe",
            "title": "Where GROVER and BPE fit",
            "source": "FluGenome3D guide card",
            "section": "Future tokenizer work",
            "topic_tags": ["future", "tokenization"],
            "text": "GROVER and BPE are future learned-tokenizer comparisons. The current app deliberately builds deterministic baselines first, so learned tokenizers can be evaluated against transparent k-mer and codon behavior.",
        },
        {
            "id": "guide_data_governance_plain",
            "title": "What the safe derived-data layer means",
            "source": "FluGenome3D guide card",
            "section": "Data governance",
            "topic_tags": ["governance"],
            "text": "The deployable layer contains aggregate summaries, binned metadata summaries, reduced PCA/t-SNE coordinates, short token summaries and hash-based visual identifiers. It does not publish raw sequences, FASTA, accessions, source record identifiers, isolate names, source sequence hashes or unrestricted sample-level metadata.",
        },
        {
            "id": "guide_scientific_framing",
            "title": "Scientific framing",
            "source": "FluGenome3D guide card",
            "section": "Scientific framing",
            "topic_tags": ["project", "representation", "governance"],
            "text": "FluGenome3D is a geometric audit and research-visualization lab for thesis-derived Influenza A HA/NA artifacts. It converts local thesis analyses into safe visual layers for sequence context, deterministic tokenization, AntigenLM-derived latent maps and public structure-alignment QC. It is descriptive, exploratory and governance-aware, not a forecasting or biological-validation system.",
        },
    ]
    return formula_cards + glossary_cards + extra_cards


def build_lab_guide() -> dict[str, Any]:
    chunks = manual_guide_cards()
    for source in GUIDE_SOURCE_FILES:
        chunks.extend(chunk_markdown_source(source))
    return {
        "schema_version": "safe-bundle-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "guide_policy": "Grounded explanatory guide built only from public docs, reports, safe formulas and governance text. It contains no raw sequences and does not call external models.",
        "answer_boundary": {
            "can_explain": [
                "dataset and panel design",
                "sequence-context metrics",
                "CDS/codon QC",
                "deterministic tokenization",
                "AntigenLM representation audit",
                "public structure viewer and mapping status",
                "data governance and allowed claims",
            ],
            "cannot_claim": [
                "antigenicity prediction",
                "vaccine recommendation",
                "escape prediction",
                "fitness, pathogenicity or transmissibility inference",
                "causal biological validation",
                "raw sequence access",
            ],
        },
        "suggested_questions": [
            "What does CpG O/E mean in this app?",
            "Why is codon usage separated from raw sequence context?",
            "How should I read the AntigenLM latent atlas?",
            "What does structure mapping still need before residue coloring?",
            "What can FluGenome3D not claim?",
            "Where would GROVER or BPE fit next?",
        ],
        "formula_cards": GUIDE_FORMULA_CARDS,
        "glossary_terms": GUIDE_GLOSSARY,
        "view_prompts": GUIDE_VIEW_PROMPTS,
        "chunks": chunks,
    }


def walk_strings(value: Any) -> list[str]:
    strings: list[str] = []
    if isinstance(value, str):
        strings.append(value)
    elif isinstance(value, dict):
        for item in value.values():
            strings.extend(walk_strings(item))
    elif isinstance(value, list):
        for item in value:
            strings.extend(walk_strings(item))
    return strings


def validate_bundle(bundle: dict[str, Any]) -> None:
    problems: list[str] = []
    for name, payload in bundle.items():
        text = json.dumps(payload, ensure_ascii=False)
        if LONG_SEQUENCE_RE.search(text):
            problems.append(f"{name}: contains a long ACGTN-like string")
        if ".fasta" in text.lower() or ".fa\"" in text.lower() or ">HA" in text or ">NA" in text:
            problems.append(f"{name}: contains FASTA-like reference")
        if name in {"tokenization_summaries.safe.json", "stability_summaries.safe.json"}:
            for string in walk_strings(payload):
                if set(string.upper()).issubset(set("ACGTN")) and len(string) > 6:
                    problems.append(f"{name}: token-like string exceeds 6 nt")
                    break
    if problems:
        raise ValueError("Safe bundle validation failed:\n" + "\n".join(problems))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name in {"representation_maps.safe.json", "antigenlm_latent_atlas.safe.json", "structure_mapping.safe.json"}:
        text = json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
    else:
        text = json.dumps(payload, indent=2, ensure_ascii=False)
    path.write_text(text + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Export deployable FluGenome3D derived JSON bundle.")
    parser.add_argument("--out", default=str(APP_DATA), help="Output directory for *.safe.json files.")
    parser.add_argument("--max-points", type=int, default=12000, help="Max reduced-coordinate points per representation.")
    args = parser.parse_args()

    out_dir = Path(args.out)
    bundle = {
        "dataset_overview.safe.json": build_dataset_overview(),
        "representation_maps.safe.json": build_representation_maps(max_points=args.max_points),
        "metric_summaries.safe.json": build_metric_summaries(),
        "tokenization_summaries.safe.json": build_tokenization_summaries(),
        "stability_summaries.safe.json": build_stability_summaries(),
        "antigenlm_latent_atlas.safe.json": build_antigenlm_latent_atlas(),
        "structure_catalog.safe.json": build_structure_catalog(),
        "structure_mapping.safe.json": build_structure_mapping_export(),
        "lab_guide.safe.json": build_lab_guide(),
        "claims_and_limits.safe.json": build_claims_and_limits(),
        "data_governance.safe.json": build_data_governance(),
    }
    validate_bundle(bundle)
    for filename, payload in bundle.items():
        write_json(out_dir / filename, payload)
    print(f"Wrote {len(bundle)} safe JSON files to {out_dir}")


if __name__ == "__main__":
    main()
