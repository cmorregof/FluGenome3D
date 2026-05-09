from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def read_fasta(path: str | Path) -> list[SeqRecord]:
    return list(SeqIO.parse(str(path), "fasta"))


def write_fasta(records: Iterable[SeqRecord], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(list(records), str(path), "fasta")


def read_metadata(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".parquet":
        return pd.read_parquet(p)
    sep = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    return pd.read_csv(p, sep=sep)
