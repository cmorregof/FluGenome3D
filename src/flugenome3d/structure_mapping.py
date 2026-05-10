from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import urlopen, urlretrieve

import numpy as np
import pandas as pd
from Bio import Align

from .codon_usage import translate_sequence


def download_pdb_cif(pdb_id: str, outdir: str | Path = "results/structures") -> Path:
    pdb_id = pdb_id.upper()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    dest = out / f"{pdb_id}.cif"
    if dest.exists():
        return dest
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    urlretrieve(url, dest)
    return dest


def make_py3dmol_html(pdb_id: str, cif_path: str | Path, out_html: str | Path) -> Path:
    # Minimal standalone py3Dmol HTML. Metrics can be added later by rewriting B-factors or adding styles.
    cif_text = Path(cif_path).read_text(encoding="utf-8", errors="replace")
    html = f"""
<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
  <title>FluGenome3D structure demo - {pdb_id}</title>
</head>
<body>
<h2>FluGenome3D structure demo: {pdb_id}</h2>
<p>Descriptive visualization only. No antigenicity, vaccine, fitness, escape or pathogenicity claims.</p>
<div id="viewer" style="width: 900px; height: 650px; position: relative;"></div>
<script>
let viewer = $3Dmol.createViewer("viewer", {{backgroundColor: "white"}});
let cif = `{cif_text.replace('`', '\\`')}`;
viewer.addModel(cif, "cif");
viewer.setStyle({{}}, {{cartoon: {{}}}});
viewer.zoomTo();
viewer.render();
</script>
</body>
</html>
"""
    out = Path(out_html)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out


FASTA_HEADER_RE = re.compile(r"^>(?P<entity>[^|]+)\|(?P<chains>[^|]+)\|(?P<name>[^|]*)\|(?P<organism>.*)$")


def parse_rcsb_fasta(text: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                entries.append(_fasta_entry(header, "".join(seq_parts)))
            header = line
            seq_parts = []
        else:
            seq_parts.append(line)
    if header:
        entries.append(_fasta_entry(header, "".join(seq_parts)))
    return entries


def _fasta_entry(header: str, sequence: str) -> dict[str, str]:
    match = FASTA_HEADER_RE.match(header)
    clean_sequence = "".join(ch for ch in sequence.upper() if ch.isalpha())
    if not match:
        return {"entity": header[1:].split("|")[0], "chains": "", "name": "", "organism": "", "sequence": clean_sequence}
    payload = match.groupdict()
    payload["sequence"] = clean_sequence
    return payload


def fetch_rcsb_fasta_entries(pdb_id: str, url_template: str = "https://www.rcsb.org/fasta/entry/{pdb_id}/display") -> list[dict[str, str]]:
    url = url_template.format(pdb_id=pdb_id.upper())
    with urlopen(url, timeout=30) as response:
        text = response.read().decode("utf-8")
    return parse_rcsb_fasta(text)


def shannon_entropy(values: Counter[str]) -> float:
    total = sum(values.values())
    if total == 0:
        return float("nan")
    entropy = 0.0
    for count in values.values():
        if count:
            p = count / total
            entropy -= p * math.log2(p)
    return entropy


def residue_signal_summary(panel: pd.DataFrame) -> pd.DataFrame:
    required = {"subtype", "protein", "refined_sequence"}
    missing = required - set(panel.columns)
    if missing:
        raise ValueError(f"CDS panel missing required columns: {sorted(missing)}")

    accum: dict[tuple[str, str, int], dict[str, object]] = {}
    for row in panel.itertuples(index=False):
        subtype = str(row.subtype)
        protein = str(row.protein)
        seq = str(row.refined_sequence).upper().replace("U", "T")
        aa = translate_sequence(seq)
        codon_count = min(len(seq) // 3, len(aa))
        for offset in range(codon_count):
            codon = seq[offset * 3 : offset * 3 + 3]
            residue = aa[offset] if offset < len(aa) else "X"
            key = (subtype, protein, offset + 1)
            if key not in accum:
                accum[key] = {
                    "subtype": subtype,
                    "protein": protein,
                    "group": f"{protein}-{subtype}",
                    "local_residue_index": offset + 1,
                    "n_codons_observed": 0,
                    "gc_sum": 0.0,
                    "cpg_sum": 0.0,
                    "upa_sum": 0.0,
                    "ambiguous_sum": 0.0,
                    "aa_counts": Counter(),
                }
            item = accum[key]
            item["n_codons_observed"] = int(item["n_codons_observed"]) + 1
            valid_bases = [base for base in codon if base in "ACGT"]
            item["gc_sum"] = float(item["gc_sum"]) + (sum(1 for base in valid_bases if base in "GC") / len(codon) if codon else 0)
            item["cpg_sum"] = float(item["cpg_sum"]) + (1.0 if "CG" in codon else 0.0)
            item["upa_sum"] = float(item["upa_sum"]) + (1.0 if "TA" in codon else 0.0)
            item["ambiguous_sum"] = float(item["ambiguous_sum"]) + (1.0 if set(codon) - set("ACGT") else 0.0)
            item["aa_counts"][residue] += 1

    rows = []
    for item in accum.values():
        n = int(item["n_codons_observed"])
        counts = item["aa_counts"]
        rows.append(
            {
                "subtype": item["subtype"],
                "protein": item["protein"],
                "group": item["group"],
                "local_residue_index": item["local_residue_index"],
                "n_codons_observed": n,
                "gc_fraction_codon": float(item["gc_sum"]) / n if n else np.nan,
                "cpg_codon_fraction": float(item["cpg_sum"]) / n if n else np.nan,
                "upa_codon_fraction": float(item["upa_sum"]) / n if n else np.nan,
                "ambiguous_codon_fraction": float(item["ambiguous_sum"]) / n if n else np.nan,
                "aa_entropy": shannon_entropy(counts),
                "dominant_aa_fraction": max(counts.values()) / n if n and counts else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["subtype", "protein", "local_residue_index"]).reset_index(drop=True)


def consensus_amino_acid_by_group(panel: pd.DataFrame) -> dict[tuple[str, str], str]:
    translated: dict[tuple[str, str], list[str]] = defaultdict(list)
    for row in panel.itertuples(index=False):
        key = (str(row.subtype), str(row.protein))
        translated[key].append(translate_sequence(str(row.refined_sequence)))
    consensus: dict[tuple[str, str], str] = {}
    for key, seqs in translated.items():
        max_len = max((len(seq) for seq in seqs), default=0)
        chars = []
        for pos in range(max_len):
            counts = Counter(seq[pos] for seq in seqs if pos < len(seq) and seq[pos] not in {"*", "X"})
            chars.append(counts.most_common(1)[0][0] if counts else "X")
        consensus[key] = "".join(chars)
    return consensus


def align_consensus_to_pdb(consensus_sequence: str, pdb_sequence: str) -> tuple[list[dict[str, int]], dict[str, float | int]]:
    if not consensus_sequence or not pdb_sequence:
        return [], {"identity": np.nan, "mapped_residues": 0, "coverage_pdb": 0, "coverage_local": 0}

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(consensus_sequence, pdb_sequence)[0]
    chunks = []
    matches = 0
    mapped = 0
    for (local_start, local_end), (pdb_start, pdb_end) in zip(alignment.aligned[0], alignment.aligned[1], strict=True):
        length = min(local_end - local_start, pdb_end - pdb_start)
        if length <= 0:
            continue
        for offset in range(length):
            local_res = local_start + offset + 1
            pdb_res = pdb_start + offset + 1
            chunks.append({"local_residue_index": int(local_res), "pdb_sequence_index": int(pdb_res)})
            mapped += 1
            if consensus_sequence[local_start + offset] == pdb_sequence[pdb_start + offset]:
                matches += 1
    metrics = {
        "identity": matches / mapped if mapped else np.nan,
        "mapped_residues": mapped,
        "coverage_pdb": mapped / len(pdb_sequence) if pdb_sequence else np.nan,
        "coverage_local": mapped / len(consensus_sequence) if consensus_sequence else np.nan,
        "local_start": min((chunk["local_residue_index"] for chunk in chunks), default=np.nan),
        "local_end": max((chunk["local_residue_index"] for chunk in chunks), default=np.nan),
    }
    return chunks, metrics


def build_structure_mapping_tables(
    panel: pd.DataFrame,
    structure_catalog: list[dict[str, str]],
    url_template: str = "https://www.rcsb.org/fasta/entry/{pdb_id}/display",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    signals = residue_signal_summary(panel)
    consensus = consensus_amino_acid_by_group(panel)
    qc_rows = []
    mapping_rows = []

    for structure in structure_catalog:
        pdb_id = structure["pdb_id"]
        subtype = structure["subtype_context"]
        protein = structure["protein"]
        local_consensus = consensus.get((subtype, protein), "")
        entries = fetch_rcsb_fasta_entries(pdb_id, url_template=url_template)
        for entry in entries:
            chunks, metrics = align_consensus_to_pdb(local_consensus, entry["sequence"])
            qc_rows.append(
                {
                    "pdb_id": pdb_id,
                    "protein": protein,
                    "subtype": subtype,
                    "pdb_entity": entry["entity"],
                    "chains": entry["chains"],
                    "pdb_sequence_length": len(entry["sequence"]),
                    "local_consensus_length": len(local_consensus),
                    **metrics,
                    "mapping_status": "alignment_qc_available" if metrics.get("mapped_residues", 0) else "unmapped",
                }
            )
            for chunk in chunks:
                mapping_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "protein": protein,
                        "subtype": subtype,
                        "group": f"{protein}-{subtype}",
                        "pdb_entity": entry["entity"],
                        "chains": entry["chains"],
                        **chunk,
                    }
                )

    qc = pd.DataFrame(qc_rows)
    mapping = pd.DataFrame(mapping_rows)
    if not mapping.empty and not signals.empty:
        mapping = mapping.merge(signals, on=["subtype", "protein", "group", "local_residue_index"], how="left")
    return qc, signals, mapping
