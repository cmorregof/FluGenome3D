#!/usr/bin/env python
from __future__ import annotations

import argparse

from flugenome3d.structure_mapping import download_pdb_cif, make_py3dmol_html
from flugenome3d.utils import load_yaml


def main() -> None:
    parser = argparse.ArgumentParser(description="Download representative structures and make minimal py3Dmol demos.")
    parser.add_argument("--pdb-config", default="config/pdbs.yml")
    parser.add_argument("--metrics", default="results/tables/sequence_metrics.csv")
    parser.add_argument("--outdir", default="results/structures")
    args = parser.parse_args()

    cfg = load_yaml(args.pdb_config)
    for item in cfg.get("structures", []):
        pdb_id = item["id"]
        cif = download_pdb_cif(pdb_id, args.outdir)
        html = make_py3dmol_html(pdb_id, cif, f"{args.outdir}/{pdb_id}_demo.html")
        print(f"Wrote {html}")


if __name__ == "__main__":
    main()
