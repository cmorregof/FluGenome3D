import { promises as fs } from "fs";
import path from "path";
import { NextResponse } from "next/server";

const SAFE_FILES = new Set([
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
  "data_governance.safe.json"
]);

async function readJsonFile(filePath: string) {
  const text = await fs.readFile(filePath, "utf-8");
  return JSON.parse(text);
}

export async function GET(_request: Request, context: { params: Promise<{ file: string }> }) {
  const { file } = await context.params;
  if (!SAFE_FILES.has(file)) {
    return NextResponse.json({ error: "Unknown safe data file." }, { status: 404 });
  }

  const cwd = process.cwd();
  const localFile = file.replace(".safe.json", ".local.json");
  const localPath = path.join(cwd, "data-local", localFile);
  const safePath = path.join(cwd, "data", file);
  const useLocal = process.env.FLUGENOME3D_DATA_MODE === "local";

  if (useLocal) {
    try {
      return NextResponse.json(await readJsonFile(localPath), {
        headers: { "x-flugenome3d-data-mode": "local-full" }
      });
    } catch {
      // Local mode is opt-in and gitignored. Fall back to the deployable derived-data layer if local files are absent.
    }
  }

  try {
    return NextResponse.json(await readJsonFile(safePath), {
      headers: { "x-flugenome3d-data-mode": "derived-data" }
    });
  } catch {
    return NextResponse.json({ error: "Safe data file is missing. Run data_export/export_vercel_safe_bundle.py." }, { status: 500 });
  }
}
