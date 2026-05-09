from __future__ import annotations

from pathlib import Path
from urllib.request import urlretrieve


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
