export const SAFE_DATA_FILES = [
  "dataset_overview.safe.json",
  "representation_maps.safe.json",
  "metric_summaries.safe.json",
  "tokenization_summaries.safe.json",
  "stability_summaries.safe.json",
  "structure_catalog.safe.json",
  "claims_and_limits.safe.json",
  "data_governance.safe.json"
] as const;

export type SafeBundle = {
  dataset: Record<string, any>;
  representations: Record<string, any>;
  metrics: Record<string, any>;
  tokenization: Record<string, any>;
  stability: Record<string, any>;
  structures: Record<string, any>;
  claims: Record<string, any>;
  governance: Record<string, any>;
};

const bundleKeys: Array<keyof SafeBundle> = [
  "dataset",
  "representations",
  "metrics",
  "tokenization",
  "stability",
  "structures",
  "claims",
  "governance"
];

export async function loadSafeBundle(): Promise<{ bundle: SafeBundle; mode: string }> {
  const responses = await Promise.all(
    SAFE_DATA_FILES.map((file) => fetch(`/api/data/${file}`, { cache: "no-store" }))
  );

  const failed = responses.find((response) => !response.ok);
  if (failed) {
    throw new Error(`Failed to load safe data bundle: ${failed.status}`);
  }

  const payloads = await Promise.all(responses.map((response) => response.json()));
  const mode = responses[0].headers.get("x-flugenome3d-data-mode") ?? "vercel-safe";
  const bundle = Object.fromEntries(bundleKeys.map((key, index) => [key, payloads[index]])) as SafeBundle;
  return { bundle, mode };
}

export function formatNumber(value: unknown, digits = 2): string {
  if (value === null || value === undefined || value === "") return "NA";
  const number = Number(value);
  if (!Number.isFinite(number)) return String(value);
  if (Math.abs(number) >= 1000) return Math.round(number).toLocaleString();
  return number.toFixed(digits).replace(/\.?0+$/, "");
}

export function groupKey(row: Record<string, any>): string {
  return row.protein_subtype ?? row.group ?? `${row.protein ?? ""}-${row.subtype ?? ""}`;
}

export function uniqueValues<T>(values: T[]): T[] {
  return Array.from(new Set(values));
}
