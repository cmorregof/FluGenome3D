#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from flugenome3d.inventory import build_inventory
from flugenome3d.utils import load_yaml


def main() -> None:
    parser = argparse.ArgumentParser(description="Discover local sequence/metadata files without modifying them.")
    parser.add_argument("--config", default="config/local_paths.yml")
    parser.add_argument("--out", default="results/tables/data_inventory.csv")
    args = parser.parse_args()

    cfg = load_yaml(args.config)
    roots = cfg.get("local_data_roots", [])
    if not roots:
        raise SystemExit("No local_data_roots configured. Copy config/local_paths.example.yml to config/local_paths.yml and edit it.")
    df = build_inventory(roots)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"Wrote inventory: {out} ({len(df)} files)")
    print(df.groupby("kind").size().to_string() if not df.empty else "No files found.")


if __name__ == "__main__":
    main()
