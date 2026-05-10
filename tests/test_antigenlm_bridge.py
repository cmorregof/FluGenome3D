from __future__ import annotations

import pandas as pd

from flugenome3d.antigenlm_bridge import safe_latent_id, stratified_sample_indices, year_bin


def test_safe_latent_id_is_hashed() -> None:
    value = safe_latent_id("EPI_ISL_example")
    assert value.startswith("lm_")
    assert "EPI_ISL" not in value
    assert value == safe_latent_id("EPI_ISL_example")


def test_year_bin() -> None:
    assert year_bin(2008) == "pre-2009"
    assert year_bin(2012) == "2009-2014"
    assert year_bin(2017) == "2015-2019"
    assert year_bin(2021) == "2020+"
    assert year_bin(None) == "unknown"


def test_stratified_sample_indices_keeps_subtypes() -> None:
    metadata = pd.DataFrame(
        {
            "subtype": ["H1N1"] * 10 + ["H3N2"] * 10,
            "year_bin": ["pre-2009"] * 5 + ["2020+"] * 5 + ["pre-2009"] * 5 + ["2020+"] * 5,
        }
    )
    idx = stratified_sample_indices(metadata, max_points=8, random_state=1)
    sampled = metadata.loc[idx]
    assert len(idx) == 8
    assert set(sampled["subtype"]) == {"H1N1", "H3N2"}
