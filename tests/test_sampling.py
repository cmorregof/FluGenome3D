from __future__ import annotations

import pandas as pd

from flugenome3d.sampling import sample_by_subtype_year


def test_sample_by_subtype_year_balances_subtypes_and_respects_target() -> None:
    rows = []
    for subtype in ["H1N1", "H3N2"]:
        for year in range(2000, 2005):
            for idx in range(10):
                rows.append({"subtype": subtype, "year": year, "internal_strain_id": f"{subtype}_{year}_{idx}"})
    df = pd.DataFrame(rows)

    sampled = sample_by_subtype_year(df, target_total=20, random_seed=42)

    assert len(sampled) == 20
    assert sampled["subtype"].value_counts().to_dict() == {"H1N1": 10, "H3N2": 10}
    assert sampled.groupby(["subtype", "year"]).size().min() >= 1
