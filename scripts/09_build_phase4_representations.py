#!/usr/bin/env python
from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy import sparse

from flugenome3d.representations import (
    build_codon_frequency_matrix,
    build_kmer_frequency_matrix,
    build_kmer_tfidf_matrix,
    build_rscu_matrix,
    save_feature_names,
)


REP_DIR = Path("data/processed/representations")
MVP_PANEL = Path("data/processed/panels/mvp_panel.parquet")
CDS_METRICS = Path("data/processed/metrics/mvp_cds_refined_codon_metrics.parquet")


def mvp_sequence_rows(panel: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for row in panel.itertuples(index=False):
        for protein, seq_col, hash_col in [("HA", "ha_sequence", "ha_sha256"), ("NA", "na_sequence", "na_sha256")]:
            seq = getattr(row, seq_col)
            rows.append(
                {
                    "internal_strain_id": row.internal_strain_id,
                    "internal_sequence_id": f"{row.internal_strain_id}_{protein}",
                    "sequence_sha256": getattr(row, hash_col),
                    "pair_sha256": row.pair_sha256,
                    "subtype": row.subtype,
                    "protein": protein,
                    "protein_subtype": f"{protein}-{row.subtype}",
                    "year": row.year,
                    "sequence": seq,
                }
            )
    return pd.DataFrame(rows)


def write_matrix(name: str, matrix, feature_names: list[str]) -> None:
    REP_DIR.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(REP_DIR / f"{name}.npz", matrix.tocsr())
    save_feature_names(feature_names, REP_DIR / f"feature_names_{name}.txt")
    print(f"Wrote {name}: shape={matrix.shape}, nnz={matrix.nnz}")


def main() -> None:
    REP_DIR.mkdir(parents=True, exist_ok=True)
    panel = pd.read_parquet(MVP_PANEL)
    seq_rows = mvp_sequence_rows(panel)
    seq_rows.drop(columns=["sequence"]).to_parquet(REP_DIR / "mvp_raw_representation_metadata.parquet", index=False)
    sequences = seq_rows["sequence"].tolist()

    for k in [3, 4, 5]:
        matrix, names = build_kmer_frequency_matrix(sequences, k=k, mode="overlapping")
        write_matrix(f"mvp_kmer{k}_freq", matrix, names)

    for k in [3, 4]:
        matrix, names = build_kmer_tfidf_matrix(sequences, k=k, mode="overlapping")
        write_matrix(f"mvp_kmer{k}_tfidf", matrix, names)

    cds = pd.read_parquet(CDS_METRICS)
    cds_meta_cols = ["internal_sequence_id", "sequence_sha256", "pair_sha256", "subtype", "protein", "year", "cds_status", "rescue_method"]
    cds_meta = cds[cds_meta_cols].copy()
    cds_meta["protein_subtype"] = cds_meta["protein"] + "-" + cds_meta["subtype"]
    cds_meta.to_parquet(REP_DIR / "mvp_cds_representation_metadata.parquet", index=False)

    codon_matrix, codon_names = build_codon_frequency_matrix(cds)
    write_matrix("mvp_codon_freq", codon_matrix, codon_names)

    rscu_matrix, rscu_names = build_rscu_matrix(cds)
    write_matrix("mvp_rscu", rscu_matrix, rscu_names)


if __name__ == "__main__":
    main()
