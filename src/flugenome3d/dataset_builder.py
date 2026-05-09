from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

from .qc import sequence_qc, sequence_sha256, sha256_text, stable_internal_id
from .sampling import sample_by_subtype_year


LOCAL_SEQUENCE_COLUMNS = {"ha_sequence", "na_sequence"}


def build_pair_table(
    paired: pd.DataFrame,
    rich_metadata: pd.DataFrame | None,
    dedup_metadata: pd.DataFrame | None,
    filters: dict[str, Any],
) -> pd.DataFrame:
    df = paired.copy()
    df["internal_strain_id"] = [
        stable_internal_id(epi, subtype, year, month, day)
        for epi, subtype, year, month, day in zip(df["epi_isl"], df["subtype"], df["year"], df["month"], df["day"], strict=False)
    ]
    df["ha_sequence"] = df["ha_sequence"].astype(str).str.upper().str.replace("U", "T", regex=False)
    df["na_sequence"] = df["na_sequence"].astype(str).str.upper().str.replace("U", "T", regex=False)
    df["ha_sha256"] = df["ha_sequence"].map(sequence_sha256)
    df["na_sha256"] = df["na_sequence"].map(sequence_sha256)
    df["pair_sha256"] = [sha256_text(f"{h}|{n}") for h, n in zip(df["ha_sha256"], df["na_sha256"], strict=False)]

    for segment, col_prefix, seq_col in [("HA", "ha", "ha_sequence"), ("NA", "na", "na_sequence")]:
        metrics = df[seq_col].map(lambda seq: sequence_qc(seq, segment, filters))
        metric_df = pd.DataFrame(metrics.tolist(), index=df.index).add_prefix(f"{col_prefix}_")
        df = pd.concat([df, metric_df], axis=1)

    df["date_available"] = df[["year", "month", "day"]].notna().all(axis=1)
    allowed_subtypes = set(filters.get("subtypes", ["H1N1", "H3N2"]))
    df["subtype_allowed"] = df["subtype"].isin(allowed_subtypes)
    df["paired_pass_qc"] = df["ha_passes_qc"] & df["na_passes_qc"] & df["date_available"] & df["subtype_allowed"]

    duplicate_counts = df.groupby("pair_sha256")["pair_sha256"].transform("size")
    df["exact_pair_duplicate_count"] = duplicate_counts.astype(int)
    df["is_exact_pair_duplicate"] = df.duplicated("pair_sha256", keep="first")

    if rich_metadata is not None and not rich_metadata.empty:
        df = df.merge(rich_metadata, on="epi_isl", how="left", validate="many_to_one")
    else:
        df["host"] = pd.NA
        df["host_is_human"] = False
        df["major_clade_rich"] = pd.NA

    if dedup_metadata is not None and not dedup_metadata.empty:
        df = df.merge(dedup_metadata, on="epi_isl", how="left", validate="many_to_one")
    else:
        df["major_clade_dedup"] = pd.NA

    df["rich_metadata_available"] = df.get("host").notna()
    df["host_is_human"] = df.get("host_is_human", False).fillna(False).astype(bool)
    df["major_clade"] = df.get("major_clade_rich").combine_first(df.get("major_clade_dedup"))
    df["major_clade_available"] = df["major_clade"].notna() & ~df["major_clade"].astype(str).str.lower().isin({"", "nan", "missing"})
    return df


def deduplicate_exact_pairs(df: pd.DataFrame) -> pd.DataFrame:
    sort_cols = [
        "paired_pass_qc",
        "rich_metadata_available",
        "host_is_human",
        "major_clade_available",
        "subtype",
        "year",
        "month",
        "internal_strain_id",
    ]
    sorted_df = df.sort_values(sort_cols, ascending=[False, False, False, False, True, True, True, True])
    dedup = sorted_df.drop_duplicates("pair_sha256", keep="first").copy()
    dedup["is_exact_pair_duplicate"] = False
    return dedup.sort_values(["subtype", "year", "month", "internal_strain_id"]).reset_index(drop=True)


def build_phase1_panels(pair_table: pd.DataFrame, filters: dict[str, Any]) -> dict[str, pd.DataFrame]:
    panel_cfg = filters.get("panels", {})
    seed = int(filters.get("random_seed", 42))
    smoke_n = int(panel_cfg.get("smoke_pairs", 400))
    mvp_n = int(panel_cfg.get("mvp_pairs", 10000))

    eligible = pair_table[pair_table["paired_pass_qc"]].copy()
    full = deduplicate_exact_pairs(eligible)

    preferred = full
    if panel_cfg.get("prefer_rich_metadata", True):
        preferred = preferred[preferred["rich_metadata_available"]]
    if panel_cfg.get("require_human_host_for_mvp", True):
        preferred = preferred[preferred["host_is_human"]]
    if len(preferred) < mvp_n:
        preferred = full

    smoke = sample_by_subtype_year(preferred, smoke_n, random_seed=seed)
    mvp = sample_by_subtype_year(preferred, mvp_n, random_seed=seed + 1)
    return {"smoke": smoke, "mvp": mvp, "full": full}


def write_panels(panels: dict[str, pd.DataFrame], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "smoke": out / "smoke_panel.parquet",
        "mvp": out / "mvp_panel.parquet",
        "full": out / "full_panel.parquet",
    }
    for name, path in paths.items():
        panels[name].to_parquet(path, index=False)
    return paths


def _length_stats(values: pd.Series) -> dict[str, float | int]:
    clean = values.dropna()
    if clean.empty:
        return {"n": 0, "min": 0, "p05": 0, "median": 0, "mean": 0, "p95": 0, "max": 0}
    return {
        "n": int(clean.size),
        "min": float(clean.min()),
        "p05": float(clean.quantile(0.05)),
        "median": float(clean.median()),
        "mean": float(clean.mean()),
        "p95": float(clean.quantile(0.95)),
        "max": float(clean.max()),
    }


def dataset_summary(pair_table: pd.DataFrame, panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    stages = {
        "all_valid_paired": pair_table[pair_table["paired_pass_qc"]],
        "exact_pair_deduplicated": panels["full"],
        "metadata_rich_subset": pair_table[pair_table["rich_metadata_available"]],
        "human_metadata_subset": pair_table[pair_table["host_is_human"]],
        "smoke_panel": panels["smoke"],
        "mvp_panel": panels["mvp"],
        "full_panel": panels["full"],
    }
    for stage, df in stages.items():
        rows.append(
            {
                "stage": stage,
                "n_pairs": int(len(df)),
                "n_ha_sequences": int(len(df)),
                "n_na_sequences": int(len(df)),
                "n_sequence_records": int(len(df) * 2),
                "n_unique_pair_hashes": int(df["pair_sha256"].nunique()) if "pair_sha256" in df else 0,
                "n_duplicate_pair_records": int(len(df) - df["pair_sha256"].nunique()) if "pair_sha256" in df else 0,
                "n_h1n1": int((df["subtype"] == "H1N1").sum()) if "subtype" in df else 0,
                "n_h3n2": int((df["subtype"] == "H3N2").sum()) if "subtype" in df else 0,
                "n_rich_metadata": int(df["rich_metadata_available"].sum()) if "rich_metadata_available" in df else 0,
                "n_human_host": int(df["host_is_human"].sum()) if "host_is_human" in df else 0,
                "year_min": int(df["year"].min()) if len(df) else 0,
                "year_max": int(df["year"].max()) if len(df) else 0,
            }
        )
    return pd.DataFrame(rows)


def panel_summary(panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for panel, df in panels.items():
        grouped = df.groupby("subtype", dropna=False)
        for subtype, group in grouped:
            rows.append(
                {
                    "panel": f"{panel}_panel",
                    "subtype": subtype,
                    "n_pairs": int(len(group)),
                    "n_sequence_records": int(2 * len(group)),
                    "n_unique_pair_hashes": int(group["pair_sha256"].nunique()),
                    "n_duplicate_pair_records": int(len(group) - group["pair_sha256"].nunique()),
                    "n_rich_metadata": int(group["rich_metadata_available"].sum()),
                    "n_human_host": int(group["host_is_human"].sum()),
                    "year_min": int(group["year"].min()) if len(group) else 0,
                    "year_max": int(group["year"].max()) if len(group) else 0,
                }
            )
    return pd.DataFrame(rows)


def missingness_summary(pair_table: pd.DataFrame) -> pd.DataFrame:
    fields = [
        "subtype",
        "year",
        "month",
        "day",
        "host",
        "region",
        "country",
        "metadata_collection_date",
        "major_clade",
        "dedup_matched",
    ]
    total = len(pair_table)
    rows = []
    for field in fields:
        if field not in pair_table:
            continue
        missing = int(pair_table[field].isna().sum())
        rows.append({"field": field, "n_rows": total, "missing": missing, "missing_fraction": missing / total if total else 0})
    return pd.DataFrame(rows)


def year_subtype_counts(pair_table: pd.DataFrame, panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = {"all_valid": pair_table[pair_table["paired_pass_qc"]], "smoke_panel": panels["smoke"], "mvp_panel": panels["mvp"], "full_panel": panels["full"]}
    rows = []
    for panel, df in frames.items():
        grouped = df.groupby(["year", "subtype"], dropna=False).size().reset_index(name="n_pairs")
        grouped.insert(0, "panel", panel)
        rows.append(grouped)
    return pd.concat(rows, ignore_index=True)


def length_qc_summary(pair_table: pd.DataFrame, panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = {"all_valid": pair_table[pair_table["paired_pass_qc"]], "full_panel": panels["full"], "mvp_panel": panels["mvp"], "smoke_panel": panels["smoke"]}
    rows = []
    for panel, df in frames.items():
        for subtype, group in df.groupby("subtype", dropna=False):
            for segment, prefix in [("HA", "ha"), ("NA", "na")]:
                stats = _length_stats(group[f"{prefix}_length"])
                rows.append(
                    {
                        "panel": panel,
                        "subtype": subtype,
                        "segment": segment,
                        **stats,
                        "n_any_ambiguity": int((group[f"{prefix}_ambiguous_fraction"] > 0).sum()),
                        "mean_ambiguous_fraction": float(group[f"{prefix}_ambiguous_fraction"].mean()) if len(group) else 0,
                        "n_non_acgtn": int((group[f"{prefix}_non_acgtn_count"] > 0).sum()),
                        "n_pass_qc": int(group[f"{prefix}_passes_qc"].sum()),
                    }
                )
    return pd.DataFrame(rows)


def duplicate_summary(pair_table: pd.DataFrame, panels: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = {"all_valid": pair_table[pair_table["paired_pass_qc"]], "smoke_panel": panels["smoke"], "mvp_panel": panels["mvp"], "full_panel": panels["full"]}
    rows = []
    for panel, df in frames.items():
        rows.append(
            {
                "panel": panel,
                "n_pairs": int(len(df)),
                "n_unique_pair_hashes": int(df["pair_sha256"].nunique()),
                "n_duplicate_pair_records": int(len(df) - df["pair_sha256"].nunique()),
            }
        )
        for subtype, group in df.groupby("subtype", dropna=False):
            rows.append(
                {
                    "panel": f"{panel}:{subtype}",
                    "n_pairs": int(len(group)),
                    "n_unique_pair_hashes": int(group["pair_sha256"].nunique()),
                    "n_duplicate_pair_records": int(len(group) - group["pair_sha256"].nunique()),
                }
            )
    return pd.DataFrame(rows)


def write_aggregate_tables(pair_table: pd.DataFrame, panels: dict[str, pd.DataFrame], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    tables = {
        "dataset_summary": (dataset_summary(pair_table, panels), out / "phase1_dataset_summary.csv"),
        "panel_summary": (panel_summary(panels), out / "phase1_panel_summary.csv"),
        "missingness_summary": (missingness_summary(pair_table), out / "phase1_missingness_summary.csv"),
        "year_subtype_counts": (year_subtype_counts(pair_table, panels), out / "phase1_year_subtype_counts.csv"),
        "length_qc_summary": (length_qc_summary(pair_table, panels), out / "phase1_length_qc_summary.csv"),
        "duplicate_summary": (duplicate_summary(pair_table, panels), out / "phase1_duplicate_summary.csv"),
    }
    paths: dict[str, Path] = {}
    for name, (df, path) in tables.items():
        safe = df.drop(columns=[col for col in LOCAL_SEQUENCE_COLUMNS if col in df.columns], errors="ignore")
        safe.to_csv(path, index=False)
        paths[name] = path
    return paths


def _plot_panel_overview(dataset_df: pd.DataFrame, panel_df: pd.DataFrame, duplicate_df: pd.DataFrame, outpath: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    panel_rows = dataset_df[dataset_df["stage"].isin(["smoke_panel", "mvp_panel", "full_panel"])]
    axes[0, 0].bar(panel_rows["stage"], panel_rows["n_pairs"], color=["#4C78A8", "#F58518", "#54A24B"])
    axes[0, 0].set_title("Pairs by panel")
    axes[0, 0].set_ylabel("HA/NA pairs")
    axes[0, 0].tick_params(axis="x", rotation=20)

    subtype = panel_df.pivot_table(index="panel", columns="subtype", values="n_pairs", aggfunc="sum", fill_value=0)
    subtype.plot(kind="bar", stacked=True, ax=axes[0, 1], color=["#72B7B2", "#E45756"])
    axes[0, 1].set_title("Subtype composition")
    axes[0, 1].set_ylabel("pairs")
    axes[0, 1].tick_params(axis="x", rotation=20)

    axes[1, 0].bar(panel_rows["stage"], panel_rows["n_ha_sequences"], label="HA", color="#4C78A8")
    axes[1, 0].bar(panel_rows["stage"], panel_rows["n_na_sequences"], bottom=panel_rows["n_ha_sequences"], label="NA", color="#F58518")
    axes[1, 0].set_title("Sequence records by panel")
    axes[1, 0].set_ylabel("sequence records")
    axes[1, 0].tick_params(axis="x", rotation=20)
    axes[1, 0].legend()

    dup_rows = duplicate_df[duplicate_df["panel"].isin(["all_valid", "full_panel", "mvp_panel", "smoke_panel"])]
    axes[1, 1].bar(dup_rows["panel"], dup_rows["n_unique_pair_hashes"], label="unique pairs", color="#59A14F")
    axes[1, 1].bar(dup_rows["panel"], dup_rows["n_duplicate_pair_records"], bottom=dup_rows["n_unique_pair_hashes"], label="duplicate records", color="#E15759")
    axes[1, 1].set_title("Exact HA+NA duplicate accounting")
    axes[1, 1].set_ylabel("pairs")
    axes[1, 1].tick_params(axis="x", rotation=20)
    axes[1, 1].legend()

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def _plot_year_distribution(year_df: pd.DataFrame, outpath: Path) -> None:
    full = year_df[year_df["panel"] == "full_panel"]
    pivot = full.pivot_table(index="year", columns="subtype", values="n_pairs", aggfunc="sum", fill_value=0).sort_index()
    fig, ax = plt.subplots(figsize=(12, 5))
    pivot.plot(kind="bar", stacked=True, ax=ax, color=["#72B7B2", "#E45756"], width=0.9)
    ax.set_title("Full deduplicated panel by year and subtype")
    ax.set_xlabel("year")
    ax.set_ylabel("HA/NA pairs")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def _plot_length_qc(panels: dict[str, pd.DataFrame], outpath: Path) -> None:
    full = panels["full"]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for ax, segment, prefix in [(axes[0, 0], "HA", "ha"), (axes[0, 1], "NA", "na")]:
        for subtype, group in full.groupby("subtype"):
            ax.hist(group[f"{prefix}_length"], bins=40, alpha=0.55, label=subtype)
        ax.set_title(f"{segment} length distribution")
        ax.set_xlabel("nt")
        ax.set_ylabel("pairs")
        ax.legend()

    ambig_rows = []
    for subtype, group in full.groupby("subtype"):
        for segment, prefix in [("HA", "ha"), ("NA", "na")]:
            ambig_rows.append({"subtype_segment": f"{subtype} {segment}", "n_any_ambiguity": int((group[f"{prefix}_ambiguous_fraction"] > 0).sum())})
    ambig = pd.DataFrame(ambig_rows)
    axes[1, 0].bar(ambig["subtype_segment"], ambig["n_any_ambiguity"], color="#B279A2")
    axes[1, 0].set_title("Records with any ambiguity")
    axes[1, 0].set_ylabel("pairs")
    axes[1, 0].tick_params(axis="x", rotation=20)

    qc_rows = []
    for subtype, group in full.groupby("subtype"):
        for segment, prefix in [("HA", "ha"), ("NA", "na")]:
            qc_rows.append({"subtype_segment": f"{subtype} {segment}", "n_pass_qc": int(group[f"{prefix}_passes_qc"].sum()), "n_fail_qc": int((~group[f"{prefix}_passes_qc"]).sum())})
    qc = pd.DataFrame(qc_rows)
    axes[1, 1].bar(qc["subtype_segment"], qc["n_pass_qc"], label="pass", color="#59A14F")
    axes[1, 1].bar(qc["subtype_segment"], qc["n_fail_qc"], bottom=qc["n_pass_qc"], label="fail", color="#E15759")
    axes[1, 1].set_title("Length/ambiguity QC")
    axes[1, 1].set_ylabel("pairs")
    axes[1, 1].tick_params(axis="x", rotation=20)
    axes[1, 1].legend()

    fig.tight_layout()
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def write_phase1_figures(table_paths: dict[str, Path], panels: dict[str, pd.DataFrame], outdir: str | Path) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    dataset_df = pd.read_csv(table_paths["dataset_summary"])
    panel_df = pd.read_csv(table_paths["panel_summary"])
    duplicate_df = pd.read_csv(table_paths["duplicate_summary"])
    year_df = pd.read_csv(table_paths["year_subtype_counts"])

    paths = {
        "overview": out / "fig1_dataset_overview.png",
        "year_subtype": out / "fig2_year_subtype_distribution.png",
        "length_qc": out / "fig3_length_qc_distribution.png",
    }
    _plot_panel_overview(dataset_df, panel_df, duplicate_df, paths["overview"])
    _plot_year_distribution(year_df, paths["year_subtype"])
    _plot_length_qc(panels, paths["length_qc"])
    return paths
