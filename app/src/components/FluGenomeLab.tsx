"use client";

import dynamic from "next/dynamic";
import { useEffect, useMemo, useState } from "react";
import {
  Activity,
  Atom,
  BookOpen,
  Database,
  Dna,
  Home,
  MessageCircle,
  Network,
  ShieldCheck,
  SplitSquareVertical
} from "lucide-react";
import MolecularViewer from "@/components/MolecularViewer";
import { SafeBundle, formatNumber, groupKey, loadSafeBundle, uniqueValues } from "@/lib/safe-data";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false }) as any;

const palette = ["#5cdce2", "#79d99c", "#edf7f4", "#4d9fb6", "#7aa88d", "#b9785f"];
const groupPalette: Record<string, string> = {
  "HA-H1N1": "#5cdce2",
  "NA-H1N1": "#79d99c",
  "HA-H3N2": "#4d9fb6",
  "NA-H3N2": "#edf7f4"
};
const subtypePalette: Record<string, string> = {
  H1N1: "#79d99c",
  H3N2: "#5cdce2"
};
const orderedGroups = ["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"];

const sequenceMetricSpecs = [
  {
    id: "gc_content",
    label: "GC fraction",
    short: "GC",
    color: "#5cdce2",
    definition: "Share of nucleotide positions that are G or C. A value of 0.42 means roughly 42% GC."
  },
  {
    id: "cpg_oe",
    label: "CpG observed/expected",
    short: "CpG O/E",
    color: "#79d99c",
    definition: "Observed CpG frequency divided by the frequency expected from C and G abundance. Values below 1 indicate relative depletion."
  },
  {
    id: "upa_oe",
    label: "UpA observed/expected",
    short: "UpA O/E",
    color: "#b7d8d1",
    definition: "DNA TA is used as the proxy for RNA UpA. Values below 1 indicate relative depletion against single-base expectations."
  }
];

const views = [
  { id: "home", label: "Home / Overview", icon: Home },
  { id: "guide", label: "Project Guide", icon: BookOpen },
  { id: "ask", label: "Ask FluGenome3D", icon: MessageCircle },
  { id: "atlas", label: "Dataset Atlas", icon: Database },
  { id: "latent", label: "AntigenLM Latent Atlas", icon: Activity },
  { id: "projector", label: "Representation Projector", icon: SplitSquareVertical },
  { id: "inspector", label: "Sequence/Token Inspector", icon: Dna },
  { id: "structure", label: "3D Molecular Viewer", icon: Atom },
  { id: "bridge", label: "Bridge View", icon: Network }
] as const;

type ViewId = (typeof views)[number]["id"];

const plotLayout = {
  paper_bgcolor: "rgba(0,0,0,0)",
  plot_bgcolor: "rgba(3,9,16,0.88)",
  font: { color: "#edf7f4", family: "IBM Plex Mono, ui-monospace, monospace", size: 11 },
  margin: { l: 48, r: 24, t: 40, b: 48 },
  xaxis: {
    gridcolor: "rgba(237,247,244,0.08)",
    zerolinecolor: "rgba(92,218,226,0.20)",
    tickangle: -25,
    automargin: true
  },
  yaxis: { gridcolor: "rgba(237,247,244,0.08)", zerolinecolor: "rgba(92,218,226,0.20)", automargin: true },
  legend: { bgcolor: "rgba(1,5,10,0.62)", bordercolor: "rgba(92,218,226,0.16)", borderwidth: 1 }
};

function Card({ label, value, detail }: { label: string; value: string | number; detail?: string }) {
  return (
    <div className="rounded-lg border border-line bg-panel/90 p-4 shadow-[0_0_30px_rgba(0,0,0,0.22)]">
      <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-muted">{label}</div>
      <div className="mt-3 text-2xl font-semibold text-ivory">{value}</div>
      {detail ? <div className="mt-2 text-xs leading-5 text-muted">{detail}</div> : null}
    </div>
  );
}

function SectionTitle({ kicker, title, children }: { kicker: string; title: string; children?: React.ReactNode }) {
  return (
    <div className="mb-5 flex flex-col gap-2 md:flex-row md:items-end md:justify-between">
      <div>
        <div className="font-mono text-[11px] uppercase tracking-[0.25em] text-brass">{kicker}</div>
        <h2 className="mt-2 text-2xl font-semibold tracking-tight text-ivory">{title}</h2>
      </div>
      {children ? <div className="text-sm text-muted">{children}</div> : null}
    </div>
  );
}

function MiniTable({ rows, columns, limit = 8 }: { rows: Record<string, any>[]; columns: string[]; limit?: number }) {
  return (
    <div className="overflow-hidden rounded-lg border border-line">
      <table className="w-full border-collapse text-left text-xs">
        <thead className="bg-panel2 font-mono uppercase tracking-[0.18em] text-muted">
          <tr>
            {columns.map((column) => (
              <th key={column} className="px-3 py-3 font-medium">
                {column.replaceAll("_", " ")}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {rows.slice(0, limit).map((row, index) => (
            <tr key={index} className="border-t border-line bg-panel/55">
              {columns.map((column) => (
                <td key={column} className="px-3 py-2 text-muted">
                  {typeof row[column] === "number" ? formatNumber(row[column], 4) : String(row[column] ?? "")}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function buildRegionPlot(rows: Record<string, any>[]) {
  const regions = uniqueValues(rows.map((row) => String(row.region ?? "unknown"))).sort();
  return ["H1N1", "H3N2"].map((subtype) => ({
    type: "bar",
    name: subtype,
    x: regions,
    y: regions.map((region) => {
      const match = rows.find((row) => row.region === region && row.subtype === subtype);
      return Number(match?.n_pairs ?? 0);
    }),
    marker: { color: subtypePalette[subtype] }
  }));
}

function buildAtlasMapTraces(countryTotals: Record<string, any>[], subtypeCounts: Record<string, any>[]) {
  const maxSubtype = Math.max(...subtypeCounts.map((row) => Number(row.n_pairs ?? 0)), 1);
  const totalTrace = {
    type: "choropleth",
    name: "Total pairs",
    locations: countryTotals.map((row) => row.country),
    z: countryTotals.map((row) => Number(row.n_pairs ?? 0)),
    text: countryTotals.map((row) => `${row.country}<br>${formatNumber(row.n_pairs, 0)} deduplicated pairs`),
    locationmode: "country names",
    colorscale: [
      [0, "rgba(8, 33, 47, 0.15)"],
      [0.45, "rgba(32, 136, 151, 0.38)"],
      [1, "rgba(92, 220, 226, 0.78)"]
    ],
    marker: { line: { color: "rgba(237,247,244,0.18)", width: 0.35 } },
    colorbar: {
      title: "pairs",
      thickness: 9,
      len: 0.62,
      tickfont: { color: "#9fb3ae" },
      titlefont: { color: "#9fb3ae" }
    },
    hovertemplate: "%{text}<extra></extra>"
  };
  const subtypeTraces = ["H1N1", "H3N2"].map((subtype) => {
    const rows = subtypeCounts.filter((row) => row.subtype === subtype);
    return {
      type: "scattergeo",
      mode: "markers",
      name: subtype,
      locations: rows.map((row) => row.country),
      locationmode: "country names",
      text: rows.map((row) => `${row.country}<br>${subtype}: ${formatNumber(row.n_pairs, 0)} pairs`),
      marker: {
        color: subtypePalette[subtype],
        opacity: subtype === "H1N1" ? 0.64 : 0.54,
        size: rows.map((row) => 5 + 35 * Math.sqrt(Number(row.n_pairs ?? 0) / maxSubtype)),
        line: { color: "rgba(237,247,244,0.78)", width: 0.65 }
      },
      hovertemplate: "%{text}<extra></extra>"
    };
  });
  return [totalTrace, ...subtypeTraces];
}

function buildProjectorTraces(points: Record<string, any>[], colorBy: string, use3d: boolean) {
  const values = uniqueValues(points.map((point) => String(point[colorBy] ?? "unknown")));
  return values.map((value, index) => {
    const subset = points.filter((point) => String(point[colorBy] ?? "unknown") === value);
    return {
      type: use3d ? "scatter3d" : "scattergl",
      mode: "markers",
      name: value,
      x: subset.map((point) => point.x),
      y: subset.map((point) => point.y),
      ...(use3d ? { z: subset.map((point) => point.z) } : {}),
      text: subset.map((point) => `${point.group} · ${point.year_bin}`),
      customdata: subset,
      marker: {
        color: groupPalette[value] ?? palette[index % palette.length],
        size: use3d ? 3.2 : 5,
        opacity: 0.72,
        line: { color: "rgba(237,247,244,0.16)", width: use3d ? 0 : 0.5 }
      },
      hovertemplate: "%{text}<br>PC1=%{x:.3f}<br>PC2=%{y:.3f}" + (use3d ? "<br>PC3=%{z:.3f}" : "") + "<extra></extra>"
    };
  });
}

function buildLatentTraces(points: Record<string, any>[], colorBy: string, use3d: boolean) {
  const values = uniqueValues(points.map((point) => String(point[colorBy] ?? "unknown")));
  return values.map((value, index) => {
    const subset = points.filter((point) => String(point[colorBy] ?? "unknown") === value);
    return {
      type: use3d ? "scatter3d" : "scattergl",
      mode: "markers",
      name: value,
      x: subset.map((point) => point.x),
      y: subset.map((point) => point.y),
      ...(use3d ? { z: subset.map((point) => point.z) } : {}),
      text: subset.map((point) => `${point.subtype} · ${point.year_bin}`),
      customdata: subset,
      marker: {
        color: subtypePalette[value] ?? palette[index % palette.length],
        size: use3d ? 2.6 : 4,
        opacity: 0.58,
        line: { color: "rgba(237,247,244,0.10)", width: use3d ? 0 : 0.3 }
      },
      hovertemplate: "%{text}<br>PC1=%{x:.3f}<br>PC2=%{y:.3f}" + (use3d ? "<br>PC3=%{z:.3f}" : "") + "<extra></extra>"
    };
  });
}

function SilhouetteSummary({ rows }: { rows: Record<string, any>[] }) {
  return (
    <div className="space-y-2">
      <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-muted">Silhouette summary</div>
      {rows.map((row) => {
        const value = Number(row.silhouette ?? 0);
        return (
          <div key={`${row.label_type}-${row.space}`} className="rounded-md border border-line bg-bg/35 p-3">
            <div className="flex items-center justify-between gap-3">
              <span className="font-mono text-[11px] uppercase tracking-[0.16em] text-muted">{String(row.label_type).replace("_", " + ")}</span>
              <span className="font-mono text-sm text-ivory">{formatNumber(value, 3)}</span>
            </div>
            <div className="mt-2 h-1.5 overflow-hidden rounded-full bg-ink">
              <span className="block h-full rounded-full bg-teal" style={{ width: `${Math.max(0, Math.min(100, value * 100))}%` }} />
            </div>
            <div className="mt-1 text-[10px] uppercase tracking-[0.14em] text-muted">{row.space}</div>
          </div>
        );
      })}
    </div>
  );
}

function decodeRepresentationPoints(rep: Record<string, any> | undefined): Record<string, any>[] {
  const rawPoints = (rep?.points ?? []) as any[];
  const schema = rep?.point_schema as string[] | undefined;
  if (!schema) {
    return rawPoints as Record<string, any>[];
  }
  return rawPoints.map((point) => {
    if (!Array.isArray(point)) {
      return point as Record<string, any>;
    }
    return Object.fromEntries(schema.map((field, index) => [field, point[index]]));
  });
}

function friendlyTokenizer(tokenizer: unknown): string {
  const value = String(tokenizer ?? "NA");
  const labels: Record<string, string> = {
    raw_overlap_k3: "Raw overlapping k=3",
    raw_overlap_k6: "Raw overlapping k=6",
    raw_nonoverlap_k6: "Raw non-overlap k=6",
    cds_frame_k6: "CDS frame-aware k=6",
    cds_codon: "CDS codons",
    cds_frame_k3: "CDS frame-aware k=3",
    cds_nonoverlap_k3_offset0: "CDS non-overlap k=3 offset 0",
    cds_nonoverlap_k3_offset1: "CDS non-overlap k=3 offset 1",
    cds_nonoverlap_k3_offset2: "CDS non-overlap k=3 offset 2"
  };
  return labels[value] ?? value.replaceAll("_", " ");
}

function metricGroup(row: Record<string, any>): string {
  return row.protein_subtype ?? `${row.protein ?? ""}-${row.subtype ?? ""}`;
}

function metricRange(rows: Record<string, any>[], metric: string) {
  const values = rows.filter((row) => row.metric === metric).map((row) => ({ group: metricGroup(row), value: Number(row.mean) }));
  const finite = values.filter((row) => Number.isFinite(row.value));
  if (!finite.length) return { min: NaN, max: NaN, minGroup: "NA", maxGroup: "NA" };
  const min = finite.reduce((best, row) => (row.value < best.value ? row : best), finite[0]);
  const max = finite.reduce((best, row) => (row.value > best.value ? row : best), finite[0]);
  return { min: min.value, max: max.value, minGroup: min.group, maxGroup: max.group };
}

function dataModeLabel(mode: string) {
  return mode === "local-full" ? "Local full layer" : "Cryptographic data layer";
}

function useSafeData() {
  const [bundle, setBundle] = useState<SafeBundle | null>(null);
  const [mode, setMode] = useState("loading");
  const [error, setError] = useState("");

  useEffect(() => {
    loadSafeBundle()
      .then((result) => {
        setBundle(result.bundle);
        setMode(result.mode);
      })
      .catch((err) => setError(err instanceof Error ? err.message : "Unknown data loading error"));
  }, []);

  return { bundle, mode, error };
}

export default function FluGenomeLab() {
  const { bundle, mode, error } = useSafeData();
  const [active, setActive] = useState<ViewId>("home");

  if (error) {
    return <div className="p-8 text-rust">{error}</div>;
  }

  if (!bundle) {
    return (
      <main className="flex min-h-screen items-center justify-center bg-ink text-ivory">
        <div className="font-mono text-sm uppercase tracking-[0.3em] text-brass">Loading FluGenome3D data layer</div>
      </main>
    );
  }

  const banner = bundle.claims.banner;
  const dataStatement = bundle.claims.data_statement;

  return (
    <main className="lab-grid min-h-screen">
      <div className="mx-auto flex min-h-screen max-w-[1560px] gap-4 px-4 py-4">
        <aside className="hidden w-72 shrink-0 rounded-lg border border-line bg-ink/88 p-4 lg:block">
          <div className="mb-6">
            <div className="font-mono text-[10px] uppercase tracking-[0.28em] text-brass">SYSTEM STATUS</div>
            <h1 className="mt-3 text-2xl font-semibold text-ivory">FluGenome3D</h1>
            <p className="mt-3 text-xs leading-5 text-muted">Visual lab for real derived Influenza A HA/NA research artifacts.</p>
          </div>
          <nav className="space-y-2">
            {views.map((view) => {
              const Icon = view.icon;
              const selected = active === view.id;
              return (
                <button
                  key={view.id}
                  onClick={() => setActive(view.id)}
                  className={`flex w-full items-center gap-3 rounded-md border px-3 py-3 text-left text-sm transition ${
                    selected
                      ? "border-lineStrong bg-brassSoft/35 text-ivory"
                      : "border-line bg-panel/50 text-muted hover:border-teal hover:text-ivory"
                  }`}
                >
                  <Icon size={16} />
                  <span>{view.label}</span>
                </button>
              );
            })}
          </nav>
          <div className="mt-6 rounded-lg border border-line bg-panel/70 p-3">
            <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-muted">DATA MODE</div>
            <div className="mt-2 text-sm text-ivory">{dataModeLabel(mode)}</div>
          </div>
        </aside>

        <section className="min-w-0 flex-1">
          <div className="mb-4 rounded-lg border border-lineStrong bg-brassSoft/20 px-4 py-3 text-sm text-ivory">
            <div className="flex flex-col gap-2 md:flex-row md:items-center md:justify-between">
              <span>{banner}</span>
              <span className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">
                {mode === "local-full" ? "LOCAL FULL LAYER" : "HASHED DATA LAYER"}
              </span>
            </div>
            <div className="mt-1 text-xs text-muted">{dataStatement}</div>
          </div>

          <div className="mb-4 flex gap-2 overflow-x-auto rounded-lg border border-line bg-ink/80 p-2 lg:hidden">
            {views.map((view) => (
              <button
                key={view.id}
                onClick={() => setActive(view.id)}
                className={`whitespace-nowrap rounded-md px-3 py-2 text-xs ${
                  active === view.id ? "bg-brass text-bg" : "bg-panel text-muted"
                }`}
              >
                {view.label}
              </button>
            ))}
          </div>

          <div
            className={
              active === "home"
                ? "overflow-hidden rounded-xl border border-line bg-bg/80"
                : "rounded-lg border border-line bg-ink/82 p-4 md:p-6"
            }
          >
            {active === "home" ? <HomeOverview bundle={bundle} mode={mode} setActive={setActive} /> : null}
            {active === "guide" ? <ProjectGuide bundle={bundle} setActive={setActive} /> : null}
            {active === "ask" ? <AskFluGenomeGuide bundle={bundle} setActive={setActive} /> : null}
            {active === "atlas" ? <DatasetAtlas bundle={bundle} mode={mode} /> : null}
            {active === "latent" ? <LatentAtlas bundle={bundle} /> : null}
            {active === "projector" ? <RepresentationProjector bundle={bundle} /> : null}
            {active === "inspector" ? <SequenceTokenInspector bundle={bundle} /> : null}
            {active === "structure" ? <StructureView bundle={bundle} /> : null}
            {active === "bridge" ? <BridgeView bundle={bundle} /> : null}
          </div>
        </section>
      </div>
    </main>
  );
}

function HomeOverview({ bundle, mode, setActive }: { bundle: SafeBundle; mode: string; setActive: (view: ViewId) => void }) {
  const features = [
    ["Sequence context", "GC, CpG/UpA, dinucleotide and k-mer summaries from derived HA/NA analyses."],
    ["Learned representation", "AntigenLM latent geometry from the thesis repo, shown through hash-based reduced coordinates."],
    ["Molecular structure", "Public RCSB structures connected to alignment QC and aggregate residue signals."]
  ];

  return (
    <div className="bg-bg">
      <section className="flu-hero-bg relative min-h-[calc(100vh-7.5rem)] overflow-hidden px-6 py-10 md:px-10 md:py-14 xl:px-14">
        <div className="absolute inset-0 bg-[linear-gradient(90deg,rgba(1,5,10,0.96)_0%,rgba(1,5,10,0.82)_29%,rgba(1,5,10,0.32)_58%,rgba(1,5,10,0.12)_100%)]" />
        <div className="absolute inset-0 bg-[linear-gradient(0deg,rgba(1,5,10,0.88)_0%,rgba(1,5,10,0.22)_44%,rgba(1,5,10,0.36)_100%)]" />
        <div className="relative flex min-h-[calc(100vh-14rem)] max-w-3xl flex-col justify-center">
          <div className="mb-6 flex flex-wrap gap-2">
            <span className="rounded-full border border-lineStrong bg-brassSoft/20 px-3 py-1 font-mono text-[10px] uppercase tracking-[0.22em] text-brass">
              {dataModeLabel(mode)}
            </span>
            <span className="rounded-full border border-teal/30 bg-teal/10 px-3 py-1 font-mono text-[10px] uppercase tracking-[0.22em] text-teal">
              Real derived artifacts
            </span>
          </div>
          <h2 className="text-6xl font-semibold tracking-[-0.045em] text-ivory md:text-8xl">FluGenome3D</h2>
          <p className="mt-6 max-w-2xl text-xl leading-9 text-ivory/90 md:text-2xl">
            A reproducible explorer connecting Influenza A sequence context, tokenization, and structural visualization.
          </p>
          <p className="mt-5 max-w-2xl text-sm leading-7 text-muted md:text-base">
            FluGenome3D turns local HA/NA analyses into maps, projections, token summaries and public structure views. The app shares derived layers with hash-based identifiers, while raw sequences remain local.
          </p>
          <div className="mt-9 flex flex-wrap gap-3">
            <button
              onClick={() => setActive("atlas")}
              className="rounded-md border border-teal/45 bg-teal/18 px-5 py-3 text-sm font-medium text-ivory transition hover:border-teal hover:bg-teal/26"
            >
              Explore dataset
            </button>
            <button
              onClick={() => setActive("ask")}
              className="rounded-md border border-line bg-bg/46 px-5 py-3 text-sm text-muted backdrop-blur transition hover:border-teal hover:text-ivory"
            >
              Ask FluGenome3D
            </button>
            <button
              onClick={() => setActive("guide")}
              className="rounded-md border border-line bg-bg/46 px-5 py-3 text-sm text-muted backdrop-blur transition hover:border-teal hover:text-ivory"
            >
              Read project guide
            </button>
            <button
              onClick={() => setActive("structure")}
              className="rounded-md border border-line bg-bg/46 px-5 py-3 text-sm text-muted backdrop-blur transition hover:border-teal hover:text-ivory"
            >
              View molecular structures
            </button>
          </div>
          <div className="mt-8 max-w-2xl rounded-lg border border-line bg-bg/44 p-4 text-xs leading-6 text-muted backdrop-blur">
            <span className="font-mono uppercase tracking-[0.18em] text-brass">DATA LAYER</span>
            <span className="ml-2">{bundle.claims.data_statement}</span>
          </div>
        </div>
      </section>

      <div className="grid gap-3 border-t border-line bg-bg px-6 py-5 md:grid-cols-3 md:px-10 xl:px-14">
        {features.map(([title, text]) => (
          <div key={title} className="rounded-lg border border-line bg-panel/62 p-5 backdrop-blur">
            <div className="text-sm font-medium text-ivory">{title}</div>
            <div className="mt-2 text-xs leading-6 text-muted">{text}</div>
          </div>
        ))}
      </div>
    </div>
  );
}

function ProjectGuide({ bundle, setActive }: { bundle: SafeBundle; setActive: (view: ViewId) => void }) {
  const reps = bundle.representations.representations as Record<string, any>[];
  const tokenizers = bundle.tokenization.tokenizer_summary as Record<string, any>[];
  const formulas = [
    ["GC fraction", "(G + C) / sequence length", "A compact readout of base composition."],
    ["CpG O/E", "f(CG) / (f(C) x f(G))", "Compares CpG frequency against what C and G abundance would suggest."],
    ["UpA O/E", "f(TA) / (f(T) x f(A))", "Uses DNA TA as the proxy for RNA UpA."],
    ["Entropy", "-sum p(token) log2 p(token)", "Higher values mean token usage is more spread out."],
    ["RSCU", "codon count / synonymous-codon average", "A codon-usage summary used only on the refined CDS panel."],
    ["JS distance", "distance between token distributions", "A descriptive way to compare groups without making prediction claims."]
  ];
  const models = [
    ["Dataset Atlas", "Country-level aggregate coverage, panel sizes and CDS reliability."],
    ["Ask FluGenome3D", "A grounded guide that answers plain-language questions from safe docs, reports and exported summaries."],
    ["AntigenLM Latent Atlas", "Learned HA+NA embedding geometry from the parent thesis repository, exported as hash-based reduced coordinates."],
    ["PCA projector", `${formatNumber(reps.length, 0)} reduced-coordinate maps built from k-mer, codon and RSCU features.`],
    ["Token audit", `${formatNumber(tokenizers.length, 0)} deterministic tokenizers: codons, overlapping k-mers, non-overlapping k-mers and frame-aware k-mers.`],
    ["Bootstrap stability", "Stratified resampling checks whether token patterns are stable under repeated sampling."],
    ["Structure viewer", "Public RCSB structures are loaded for inspection; residue-level metric coloring remains a future validated mapping step."]
  ];

  return (
    <div>
      <SectionTitle kicker="PROJECT GUIDE" title="What FluGenome3D is trying to show">
        A readable map of the project logic, formulas and current models.
      </SectionTitle>

      <div className="grid gap-4 xl:grid-cols-[1.05fr_0.95fr]">
        <div className="rounded-lg border border-line bg-panel/75 p-5">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">PROJECT IDEA</div>
          <p className="mt-3 text-sm leading-7 text-muted">
            FluGenome3D is a visual research layer for Influenza A HA/NA. It starts from real local analyses, exports only derived data, and helps a viewer move from dataset coverage to sequence composition, tokenization, representation space and public molecular structures.
          </p>
          <p className="mt-3 text-sm leading-7 text-muted">
            The project is not trying to predict vaccine candidates, escape, pathogenicity or fitness. It is building a clear descriptive baseline: what is in the dataset, how sequences are represented, and which patterns are stable enough to inspect before learned tokenizers such as BPE or GROVER.
          </p>
          <div className="mt-5 flex flex-wrap gap-3">
            <button onClick={() => setActive("atlas")} className="rounded-md border border-teal/45 bg-teal/16 px-4 py-2 text-sm text-ivory hover:border-teal">
              Open atlas
            </button>
            <button onClick={() => setActive("ask")} className="rounded-md border border-teal/45 bg-teal/16 px-4 py-2 text-sm text-ivory hover:border-teal">
              Ask the guide
            </button>
            <button onClick={() => setActive("projector")} className="rounded-md border border-line bg-bg/45 px-4 py-2 text-sm text-muted hover:border-teal hover:text-ivory">
              Open projector
            </button>
            <button onClick={() => setActive("inspector")} className="rounded-md border border-line bg-bg/45 px-4 py-2 text-sm text-muted hover:border-teal hover:text-ivory">
              Open metrics
            </button>
          </div>
        </div>

        <div className="rounded-lg border border-line bg-panel/75 p-5">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">DATA LAYER</div>
          <div className="mt-3 grid gap-3 md:grid-cols-2">
            <div className="rounded-md border border-line bg-bg/30 p-3">
              <div className="text-sm font-semibold text-ivory">Shared app layer</div>
              <p className="mt-1 text-xs leading-5 text-muted">Aggregates, reduced coordinates, short tokens and hash-based point IDs.</p>
            </div>
            <div className="rounded-md border border-line bg-bg/30 p-3">
              <div className="text-sm font-semibold text-ivory">Local-only layer</div>
              <p className="mt-1 text-xs leading-5 text-muted">Raw sequences, restricted panels, detailed Parquet files and private metadata stay off the public app.</p>
            </div>
          </div>
          <p className="mt-4 text-xs leading-6 text-muted">
            "Cryptographic" here refers to hash-based internal IDs and a derived export boundary. It is a research-sharing layer, not a claim that raw data are published.
          </p>
        </div>
      </div>

      <div className="mt-4 grid gap-4 xl:grid-cols-[0.95fr_1.05fr]">
        <div className="rounded-lg border border-line bg-panel/75 p-5">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">FORMULAS, IN PLAIN LANGUAGE</div>
          <div className="mt-4 grid gap-3">
            {formulas.map(([name, formula, explanation]) => (
              <div key={name} className="rounded-md border border-line bg-bg/30 p-3">
                <div className="flex flex-col gap-1 md:flex-row md:items-baseline md:justify-between">
                  <div className="text-sm font-semibold text-ivory">{name}</div>
                  <code className="font-mono text-xs text-teal">{formula}</code>
                </div>
                <p className="mt-2 text-xs leading-5 text-muted">{explanation}</p>
              </div>
            ))}
          </div>
        </div>

        <div className="rounded-lg border border-line bg-panel/75 p-5">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">MODELS AND VIEWS USED NOW</div>
          <div className="mt-4 grid gap-3">
            {models.map(([name, explanation]) => (
              <div key={name} className="rounded-md border border-line bg-bg/30 p-3">
                <div className="text-sm font-semibold text-ivory">{name}</div>
                <p className="mt-1 text-xs leading-5 text-muted">{explanation}</p>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
}

type AskGuideResponse = {
  answer: string;
  citations: Array<{ marker: string; title: string; source: string; section: string; topic_tags?: string[] }>;
  matched_topics: string[];
  guardrails: string[];
};

function AskFluGenomeGuide({ bundle, setActive }: { bundle: SafeBundle; setActive: (view: ViewId) => void }) {
  const suggestions = (bundle.guide.suggested_questions ?? []) as string[];
  const [question, setQuestion] = useState("What does CpG O/E mean in this app?");
  const [response, setResponse] = useState<AskGuideResponse | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState("");

  async function askGuide(nextQuestion = question) {
    const trimmed = nextQuestion.trim();
    if (!trimmed) return;
    setQuestion(trimmed);
    setLoading(true);
    setError("");
    try {
      const result = await fetch("/api/ask", {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ question: trimmed }),
      });
      const payload = await result.json();
      if (!result.ok) throw new Error(payload.error ?? "Guide request failed");
      setResponse(payload as AskGuideResponse);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Guide request failed");
    } finally {
      setLoading(false);
    }
  }

  return (
    <div>
      <SectionTitle kicker="ASK FLUGENOME3D" title="A grounded guide for the visual lab">
        Plain-language answers from safe reports, formulas and exported summaries.
      </SectionTitle>

      <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_360px]">
        <div className="rounded-lg border border-line bg-panel/75 p-5">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">ASK A QUESTION</div>
          <p className="mt-3 max-w-3xl text-sm leading-7 text-muted">
            This guide is meant to make the app understandable without turning it into a black box. It answers from the project reports, guide cards and safe JSON manifests already shipped with FluGenome3D.
          </p>
          <div className="mt-5 flex flex-col gap-3 md:flex-row">
            <input
              value={question}
              onChange={(event) => setQuestion(event.target.value)}
              onKeyDown={(event) => {
                if (event.key === "Enter") askGuide();
              }}
              className="min-w-0 flex-1 rounded-md border border-line bg-bg/65 px-4 py-3 text-sm text-ivory outline-none transition placeholder:text-muted focus:border-teal"
              placeholder="Ask about CpG, RSCU, PCA, AntigenLM, structure mapping..."
            />
            <button
              onClick={() => askGuide()}
              disabled={loading}
              className="rounded-md border border-teal/45 bg-teal/18 px-5 py-3 text-sm font-medium text-ivory transition hover:border-teal hover:bg-teal/26 disabled:cursor-not-allowed disabled:opacity-60"
            >
              {loading ? "Reading..." : "Ask"}
            </button>
          </div>

          <div className="mt-4 flex flex-wrap gap-2">
            {suggestions.slice(0, 6).map((item) => (
              <button
                key={item}
                onClick={() => askGuide(item)}
                className="rounded-full border border-line bg-bg/40 px-3 py-1.5 text-xs text-muted transition hover:border-teal hover:text-ivory"
              >
                {item}
              </button>
            ))}
          </div>

          <div className="mt-5 min-h-[340px] rounded-lg border border-line bg-bg/42 p-5">
            {error ? <div className="text-sm text-rust">{error}</div> : null}
            {!response && !error ? (
              <div className="flex h-full min-h-[280px] flex-col justify-center text-sm leading-7 text-muted">
                <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-teal">GROUNDING LAYER READY</div>
                <p className="mt-3">
                  Try asking: what does a metric mean, why a panel exists, how to read a projection, what AntigenLM contributes, or what structure mapping still needs.
                </p>
              </div>
            ) : null}
            {response ? (
              <div>
                <div className="whitespace-pre-wrap text-sm leading-7 text-muted">{response.answer}</div>
                <div className="mt-5 grid gap-3 md:grid-cols-2">
                  {response.citations.map((citation) => (
                    <div key={`${citation.marker}-${citation.source}-${citation.section}`} className="rounded-md border border-line bg-panel/55 p-3">
                      <div className="font-mono text-[10px] uppercase tracking-[0.18em] text-teal">{citation.marker}</div>
                      <div className="mt-1 text-sm font-medium text-ivory">{citation.title}</div>
                      <div className="mt-1 text-xs leading-5 text-muted">{citation.source}</div>
                    </div>
                  ))}
                </div>
              </div>
            ) : null}
          </div>
        </div>

        <div className="space-y-4">
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">HOW IT WORKS</div>
            <div className="mt-3 space-y-3 text-sm leading-6 text-muted">
              <p>Retrieval is local to the app: no raw sequence access, no external LLM call, no hidden biological predictor.</p>
              <p>Answers cite safe guide chunks generated from docs, reports, formulas and governance manifests.</p>
              <p>When a question asks for prediction or causal meaning, the guide keeps the answer inside the descriptive scope.</p>
            </div>
          </div>
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">FAST PATHS</div>
            <div className="mt-3 grid gap-2">
              {[
                ["Project formulas", "guide"],
                ["Dataset atlas", "atlas"],
                ["AntigenLM layer", "latent"],
                ["Structure status", "structure"],
              ].map(([label, view]) => (
                <button
                  key={label}
                  onClick={() => setActive(view as ViewId)}
                  className="rounded-md border border-line bg-bg/35 px-3 py-2 text-left text-sm text-muted transition hover:border-teal hover:text-ivory"
                >
                  {label}
                </button>
              ))}
            </div>
          </div>
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">SOURCE CORPUS</div>
            <div className="mt-3 text-3xl font-semibold text-ivory">{formatNumber((bundle.guide.chunks ?? []).length, 0)}</div>
            <p className="mt-2 text-xs leading-5 text-muted">
              Safe explanation chunks. No FASTA, raw sequences, restricted panels, accessions or isolate names.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}

function DatasetAtlas({ bundle, mode }: { bundle: SafeBundle; mode: string }) {
  const datasetRows = bundle.dataset.dataset_summary as Record<string, any>[];
  const refined = bundle.dataset.cds_refined_qc as Record<string, any>[];
  const atlas = bundle.dataset.geographic_atlas ?? {};
  const panelOptions = (atlas.panels ?? []) as Record<string, any>[];
  const [atlasPanel, setAtlasPanel] = useState("full_panel");
  const mvp = datasetRows.find((row) => row.stage === "mvp_panel");
  const full = datasetRows.find((row) => row.stage === "full_panel");
  const countryTotals = ((atlas.country_totals ?? []) as Record<string, any>[]).filter((row) => row.panel === atlasPanel);
  const subtypeCounts = ((atlas.country_subtype_counts ?? []) as Record<string, any>[]).filter((row) => row.panel === atlasPanel);
  const regionRows = ((atlas.region_summary ?? []) as Record<string, any>[]).filter((row) => row.panel === atlasPanel);
  const topCountries = countryTotals.slice(0, 10).map((country) => {
    const h1n1 = subtypeCounts.find((row) => row.country === country.country && row.subtype === "H1N1")?.n_pairs ?? 0;
    const h3n2 = subtypeCounts.find((row) => row.country === country.country && row.subtype === "H3N2")?.n_pairs ?? 0;
    const total = Number(country.n_pairs ?? 0) || 1;
    return {
      country: country.country,
      region: country.region,
      n_pairs: country.n_pairs,
      h1n1_share: Number(h1n1) / total,
      h3n2_share: Number(h3n2) / total,
    };
  });
  const totalMappedPairs = countryTotals.reduce((sum, row) => sum + Number(row.n_pairs ?? 0), 0);
  const yearMin = Math.min(...countryTotals.map((row) => Number(row.year_min ?? Infinity)));
  const yearMax = Math.max(...countryTotals.map((row) => Number(row.year_max ?? -Infinity)));
  const mapTraces = useMemo(() => buildAtlasMapTraces(countryTotals, subtypeCounts), [countryTotals, subtypeCounts]);
  const regionTraces = useMemo(() => buildRegionPlot(regionRows), [regionRows]);

  return (
    <div>
      <SectionTitle kicker="DATASET ATLAS" title="Geographic coverage of HA/NA pairs">
        Country-level aggregates from real derived artifacts. No raw sequences, accessions or isolate names.
      </SectionTitle>
      <div className="mb-4 grid gap-3 md:grid-cols-[1.2fr_1fr_1fr_1fr]">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Atlas layer</span>
          <select
            value={atlasPanel}
            onChange={(event) => setAtlasPanel(event.target.value)}
            className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory"
          >
            {panelOptions.map((item) => (
              <option key={item.panel} value={item.panel}>
                {item.label}
              </option>
            ))}
          </select>
        </label>
        <Card label="Mapped countries" value={formatNumber(countryTotals.length, 0)} detail="Country aggregates only" />
        <Card label="Mapped pairs" value={formatNumber(totalMappedPairs, 0)} detail={atlasPanel === "full_panel" ? "Full deduplicated panel" : "Balanced MVP panel"} />
        <Card label="Temporal span" value={`${Number.isFinite(yearMin) ? yearMin : "NA"}-${Number.isFinite(yearMax) ? yearMax : "NA"}`} detail={dataModeLabel(mode)} />
      </div>

      <div className="grid gap-4 xl:grid-cols-[1fr_350px]">
        <div className="plot-shell overflow-hidden rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={mapTraces}
            layout={{
              ...plotLayout,
              title: atlasPanel === "full_panel" ? "Full deduplicated HA/NA atlas" : "Balanced MVP HA/NA atlas",
              height: 470,
              margin: { l: 8, r: 8, t: 44, b: 8 },
              geo: {
                domain: { x: [0, 1], y: [0.18, 1] },
                center: { lat: 8, lon: 5 },
                projection: { type: "natural earth", scale: 1.34 },
                bgcolor: "rgba(0,0,0,0)",
                showframe: false,
                showcoastlines: true,
                coastlinecolor: "rgba(237,247,244,0.20)",
                showcountries: true,
                countrycolor: "rgba(237,247,244,0.13)",
                showland: true,
                landcolor: "rgba(8, 27, 38, 0.92)",
                showocean: true,
                oceancolor: "rgba(1, 5, 10, 0.96)",
                lataxis: { showgrid: false },
                lonaxis: { showgrid: false }
              }
            }}
            config={{ responsive: true, displaylogo: false, scrollZoom: false, staticPlot: true }}
            style={{ width: "100%", height: "470px" }}
          />
        </div>
        <div className="space-y-4">
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">Top countries</div>
            <div className="mt-3 space-y-3">
              {topCountries.map((row) => (
                <div key={row.country}>
                  <div className="flex items-center justify-between gap-3 text-xs">
                    <span className="font-medium text-ivory">{row.country}</span>
                    <span className="font-mono text-muted">{formatNumber(row.n_pairs, 0)}</span>
                  </div>
                  <div className="mt-2 flex h-2 overflow-hidden rounded-full bg-bg/80">
                    <span className="h-full bg-sage" style={{ width: `${Math.round(row.h1n1_share * 100)}%` }} />
                    <span className="h-full bg-teal/75" style={{ width: `${Math.round(row.h3n2_share * 100)}%` }} />
                  </div>
                  <div className="mt-1 text-[10px] uppercase tracking-[0.14em] text-muted">{row.region}</div>
                </div>
              ))}
            </div>
            <div className="mt-4 flex gap-4 text-[11px] text-muted">
              <span><span className="mr-1 inline-block h-2 w-2 rounded-full bg-sage" />H1N1</span>
              <span><span className="mr-1 inline-block h-2 w-2 rounded-full bg-teal" />H3N2</span>
            </div>
          </div>
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">CDS RELIABILITY</div>
            <MiniTable
              rows={refined}
              columns={["subtype", "protein", "n_refined_sequences", "n_strict_pass", "n_rescued", "n_unrescued"]}
            />
          </div>
        </div>
      </div>

      <div className="mt-4 grid gap-4 xl:grid-cols-[0.95fr_1.05fr]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={regionTraces}
            layout={{ ...plotLayout, barmode: "stack", title: "Regional subtype composition", height: 320 }}
            config={{ responsive: true, displaylogo: false }}
            useResizeHandler
            style={{ width: "100%", height: "320px" }}
          />
        </div>
        <div className="rounded-lg border border-line bg-panel/75 p-4">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">Panel interpretation</div>
          <div className="mt-3 grid gap-3 md:grid-cols-2">
            <div className="rounded-lg border border-line bg-bg/36 p-4">
              <div className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">MVP panel</div>
              <div className="mt-2 text-2xl font-semibold text-ivory">{formatNumber(mvp?.n_sequence_records ?? 0, 0)}</div>
              <div className="mt-2 text-xs leading-5 text-muted">Balanced subset for CV-ready sequence-context figures.</div>
            </div>
            <div className="rounded-lg border border-line bg-bg/36 p-4">
              <div className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Full panel</div>
              <div className="mt-2 text-2xl font-semibold text-ivory">{formatNumber(full?.n_pairs ?? 0, 0)}</div>
              <div className="mt-2 text-xs leading-5 text-muted">Deduplicated HA+NA pairs for coverage context.</div>
            </div>
          </div>
          <p className="mt-4 text-sm leading-6 text-muted">
            Smoke-test rows are intentionally not foregrounded here. This atlas emphasizes the balanced MVP panel and the full deduplicated panel as the useful research views.
          </p>
          <p className="mt-3 text-xs leading-6 text-muted">
            Map locations are country-level aggregates from safe derived metadata. They are descriptive coverage summaries, not sample-level locations.
          </p>
        </div>
      </div>
    </div>
  );
}

function LatentAtlas({ bundle }: { bundle: SafeBundle }) {
  const atlas = bundle.antigenlm ?? {};
  const projection = atlas.projection ?? {};
  const points = useMemo(() => decodeRepresentationPoints(projection), [projection]);
  const [colorBy, setColorBy] = useState("subtype");
  const [projectionMode, setProjectionMode] = useState<"3d" | "2d">("3d");
  const [selected, setSelected] = useState<Record<string, any> | null>(null);
  const use3d = projectionMode === "3d";
  const traces = useMemo(() => buildLatentTraces(points, colorBy, use3d), [points, colorBy, use3d]);
  const cache = atlas.cache_summary ?? {};
  const spearman = (atlas.spearman_summary ?? []) as Record<string, any>[];
  const hammingHaNa = spearman.filter((row) => row.metric === "hamming_ha_na");
  const temporal = spearman.filter((row) => row.metric === "temporal");
  const clade = (atlas.clade_enrichment_summary ?? []) as Record<string, any>[];
  const pca = (atlas.pca_summary ?? []) as Record<string, any>[];
  const comparison = (atlas.representation_family_comparison ?? []) as Record<string, any>[];
  const explained = (projection.pca_explained_variance ?? []) as Array<number | null>;
  const globalPca = pca.find((row) => row.group === "global");
  const meanHamming = hammingHaNa.reduce((sum, row) => sum + Number(row.rho_mean ?? 0), 0) / Math.max(hammingHaNa.length, 1);
  const meanTemporal = temporal.reduce((sum, row) => sum + Number(row.rho_mean ?? 0), 0) / Math.max(temporal.length, 1);

  useEffect(() => setSelected(null), [colorBy, projectionMode]);

  return (
    <div>
      <SectionTitle kicker="ANTIGENLM LATENT ATLAS" title="Learned HA/NA representation layer">
        AntigenLM embeddings from the parent thesis repo, exported as hash-based PCA coordinates.
      </SectionTitle>

      <div className="mb-4 grid gap-3 md:grid-cols-2 xl:grid-cols-4">
        <Card label="Latent records" value={formatNumber(cache.n_records ?? 0, 0)} detail={`${formatNumber(cache.embedding_dim ?? 0, 0)}-dimensional AntigenLM embeddings`} />
        <Card label="Exported points" value={formatNumber(projection.n_exported_points ?? 0, 0)} detail="Hash IDs and coarse metadata only" />
        <Card label="HA+NA molecular rho" value={formatNumber(meanHamming, 3)} detail="Mean Spearman vs Hamming proxy" />
        <Card label="Global PCA n95" value={formatNumber(globalPca?.n95 ?? "NA", 0)} detail="Components for 95% variance in parent audit" />
      </div>

      <div className="mb-4 grid gap-3 xl:grid-cols-[1fr_1fr_1fr]">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Color by</span>
          <select value={colorBy} onChange={(event) => setColorBy(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {["subtype", "year_bin", "source", "representation"].map((value) => (
              <option key={value} value={value}>
                {value}
              </option>
            ))}
          </select>
        </label>
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Projection</span>
          <select value={projectionMode} onChange={(event) => setProjectionMode(event.target.value as "3d" | "2d")} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            <option value="3d">3D PCA</option>
            <option value="2d">2D PCA</option>
          </select>
        </label>
        <div className="rounded-lg border border-line bg-panel/70 p-4 text-sm leading-6 text-muted">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-brass">Reading this layer</span>
          <p className="mt-2">The learned space is compared with molecular proxies and temporal locality. It is not a prediction or sequence-generation panel.</p>
        </div>
      </div>

      <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_360px]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={traces}
            layout={
              use3d
                ? {
                    ...plotLayout,
                    title: projection.label ?? "AntigenLM latent PCA",
                    height: 620,
                    margin: { l: 8, r: 8, t: 42, b: 8 },
                    scene: {
                      xaxis: { title: "PC1", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", backgroundcolor: "rgba(1,5,10,0.72)" },
                      yaxis: { title: "PC2", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", backgroundcolor: "rgba(1,5,10,0.72)" },
                      zaxis: { title: "PC3", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", backgroundcolor: "rgba(1,5,10,0.72)" },
                      camera: { eye: { x: 1.55, y: 1.35, z: 0.95 } }
                    }
                  }
                : {
                    ...plotLayout,
                    title: projection.label ?? "AntigenLM latent PCA",
                    height: 620,
                    xaxis: { ...plotLayout.xaxis, title: "PC1", tickformat: ".2f", nticks: 7, tickangle: 0 },
                    yaxis: { ...plotLayout.yaxis, title: "PC2", tickformat: ".2f", nticks: 7 },
                  }
            }
            config={{ responsive: true, displaylogo: false, scrollZoom: false }}
            useResizeHandler
            style={{ width: "100%", height: "620px" }}
            onClick={(event: any) => setSelected(event.points?.[0]?.customdata ?? null)}
          />
        </div>

        <div className="space-y-4">
          <div className="rounded-lg border border-line bg-panel/80 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">LATENT POINT</div>
            {selected ? (
              <div className="mt-4 space-y-3 text-sm text-muted">
                {["id", "subtype", "year_bin", "representation", "source"].map((field) => (
                  <div key={field} className="flex justify-between gap-3 border-b border-line pb-2">
                    <span className="font-mono uppercase tracking-[0.14em]">{field}</span>
                    <span className="text-right text-ivory">{String(selected[field])}</span>
                  </div>
                ))}
              </div>
            ) : (
              <p className="mt-4 text-sm leading-6 text-muted">Select a point to inspect its hash ID and coarse metadata. No source identifiers are exposed.</p>
            )}
          </div>
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">PCA variance</div>
            <div className="mt-3 grid grid-cols-3 gap-2">
              {["PC1", "PC2", "PC3"].map((axis, index) => (
                <div key={axis} className="rounded-md border border-line bg-bg/35 p-2 text-center">
                  <div className="font-mono text-[10px] text-muted">{axis}</div>
                  <div className="text-sm font-semibold text-ivory">{explained[index] == null ? "NA" : `${formatNumber(Number(explained[index]) * 100, 1)}%`}</div>
                </div>
              ))}
            </div>
          </div>
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">Interpretation</div>
            <p className="mt-3 text-sm leading-6 text-muted">
              The parent audit found stronger latent correlation with HA+NA Hamming distance than with global time distance. In plain terms: this learned space is organized more by molecular similarity than by a simple calendar line.
            </p>
            <p className="mt-2 text-xs leading-5 text-muted">
              Mean HA+NA rho: <span className="text-ivory">{formatNumber(meanHamming, 3)}</span>. Mean temporal rho:{" "}
              <span className="text-ivory">{formatNumber(meanTemporal, 3)}</span>.
            </p>
          </div>
        </div>
      </div>

      <div className="mt-4 grid gap-4 xl:grid-cols-[1fr_1fr]">
        <div>
          <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">MOLECULAR GEOMETRY</div>
          <MiniTable rows={spearman.map((row) => ({ metric: row.metric, subtype: row.subtype, rho_mean: row.rho_mean, valid_pairs: row.valid_pairs_mean }))} columns={["metric", "subtype", "rho_mean", "valid_pairs"]} limit={8} />
        </div>
        <div>
          <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">CLADE / BASELINE CHECKS</div>
          <MiniTable rows={clade.map((row) => ({ subtype: row.subtype, label: row.label, k: row.k, precision: row.mean_precision, enrichment: row.enrichment_vs_random }))} columns={["subtype", "label", "k", "precision", "enrichment"]} limit={8} />
        </div>
      </div>

      <div className="mt-4 rounded-lg border border-line bg-panel/75 p-4">
        <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">REPRESENTATION LADDER</div>
        <p className="mt-2 text-sm leading-6 text-muted">
          FluGenome3D now shows a ladder from interpretable baselines to a learned biological representation: raw k-mers, codon/RSCU vectors, deterministic tokenizers, and AntigenLM embeddings.
        </p>
        <div className="mt-3">
          <MiniTable rows={comparison} columns={["family", "representation", "n_sequences", "n_features", "app_role"]} limit={8} />
        </div>
      </div>
    </div>
  );
}

function RepresentationProjector({ bundle }: { bundle: SafeBundle }) {
  const reps = ((bundle.representations.representations as Record<string, any>[]) ?? []).filter((item) => !String(item.id).includes("umap"));
  const [repId, setRepId] = useState(reps[0]?.id ?? "");
  const [colorBy, setColorBy] = useState("group");
  const [projectionMode, setProjectionMode] = useState<"3d" | "2d">("3d");
  const [selected, setSelected] = useState<Record<string, any> | null>(null);
  const rep = reps.find((item) => item.id === repId) ?? reps[0];
  const points = useMemo(() => decodeRepresentationPoints(rep), [rep]);
  const has3d = points.some((point) => Number.isFinite(Number(point.z)));
  const use3d = projectionMode === "3d" && has3d;
  const explained = (rep?.pca_explained_variance ?? []) as Array<number | null>;

  const traces = useMemo(() => buildProjectorTraces(points, colorBy, use3d), [points, colorBy, use3d]);
  useEffect(() => setSelected(null), [repId, colorBy, projectionMode]);
  const projectorLayout = useMemo(
    () =>
      use3d
        ? {
            ...plotLayout,
            title: rep?.label ?? "Representation",
            height: 620,
            margin: { l: 8, r: 8, t: 42, b: 8 },
            scene: {
              xaxis: { title: "PC1", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", zerolinecolor: "rgba(92,218,226,0.20)", backgroundcolor: "rgba(1,5,10,0.72)" },
              yaxis: { title: "PC2", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", zerolinecolor: "rgba(92,218,226,0.20)", backgroundcolor: "rgba(1,5,10,0.72)" },
              zaxis: { title: "PC3", tickformat: ".2f", nticks: 5, gridcolor: "rgba(237,247,244,0.08)", zerolinecolor: "rgba(92,218,226,0.20)", backgroundcolor: "rgba(1,5,10,0.72)" },
              camera: { eye: { x: 1.45, y: 1.45, z: 0.95 } }
            }
          }
        : {
            ...plotLayout,
            title: rep?.label ?? "Representation",
            height: 620,
            xaxis: { ...plotLayout.xaxis, title: "PC1", tickformat: ".2f", nticks: 7, tickangle: 0 },
            yaxis: { ...plotLayout.yaxis, title: "PC2", tickformat: ".2f", nticks: 7 },
          },
    [rep?.label, use3d]
  );

  return (
    <div>
      <SectionTitle kicker="REPRESENTATION PROJECTOR" title="Reduced-coordinate maps">
        PCA coordinates are real derived artifacts with hashed IDs.
      </SectionTitle>
      <div className="mb-4 grid gap-3 xl:grid-cols-[1.2fr_0.8fr_0.8fr_1fr]">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Representation</span>
          <select value={repId} onChange={(event) => setRepId(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {reps.map((item) => (
              <option key={item.id} value={item.id}>
                {item.label}
              </option>
            ))}
          </select>
          <div className="mt-2 text-xs leading-5 text-muted">{rep?.description ?? "Safe derived coordinate map."}</div>
        </label>
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Color by</span>
          <select value={colorBy} onChange={(event) => setColorBy(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {["group", "protein", "subtype", "year_bin", "cds_status"].map((value) => (
              <option key={value} value={value}>
                {value}
              </option>
            ))}
          </select>
        </label>
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Projection</span>
          <select value={projectionMode} onChange={(event) => setProjectionMode(event.target.value as "3d" | "2d")} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            <option value="3d">3D PCA</option>
            <option value="2d">2D PCA</option>
          </select>
          <div className="mt-2 text-xs leading-5 text-muted">Use 3D to inspect separation without overcrowded axes.</div>
        </label>
        <Card label="Exported points" value={formatNumber(rep?.n_exported_points ?? 0, 0)} detail={rep?.privacy ?? "safe point metadata"} />
      </div>
      <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_360px]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={traces}
            layout={projectorLayout}
            config={{ responsive: true, displaylogo: false, scrollZoom: false }}
            useResizeHandler
            style={{ width: "100%", height: "620px" }}
            onClick={(event: any) => setSelected(event.points?.[0]?.customdata ?? null)}
          />
        </div>
        <div className="rounded-lg border border-line bg-panel/80 p-4">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">POINT / MAP INSPECTOR</div>
          {selected ? (
            <div className="mt-4 space-y-3 text-sm text-muted">
              {["id", "group", "protein", "subtype", "year_bin", "cds_status"].map((field) => (
                <div key={field} className="flex justify-between gap-3 border-b border-line pb-2">
                  <span className="font-mono uppercase tracking-[0.14em]">{field}</span>
                  <span className="text-right text-ivory">{String(selected[field])}</span>
                </div>
              ))}
            </div>
          ) : (
            <p className="mt-4 text-sm leading-6 text-muted">Select a point to inspect its safe hash ID and minimal group metadata.</p>
          )}
          <div className="mt-6 rounded-lg border border-line bg-bg/25 p-3">
            <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-muted">PCA variance</div>
            <div className="mt-2 grid grid-cols-3 gap-2 text-center">
              {["PC1", "PC2", "PC3"].map((axis, index) => (
                <div key={axis} className="rounded-md border border-line bg-panel/55 px-2 py-2">
                  <div className="font-mono text-[10px] text-muted">{axis}</div>
                  <div className="text-sm font-semibold text-ivory">{explained[index] == null ? "NA" : `${formatNumber(Number(explained[index]) * 100, 1)}%`}</div>
                </div>
              ))}
            </div>
          </div>
          <div className="mt-4">
            <SilhouetteSummary rows={rep?.silhouette_scores ?? []} />
          </div>
        </div>
      </div>
    </div>
  );
}

function SequenceTokenInspector({ bundle }: { bundle: SafeBundle }) {
  const gcRows = (bundle.metrics.gc_cpg_upa_summary as Record<string, any>[]).filter((row) =>
    ["gc_content", "cpg_oe", "upa_oe"].includes(row.metric)
  );
  const entropy = bundle.tokenization.entropy_by_group as Record<string, any>[];
  const ranking = bundle.stability.tokenizer_robustness_ranking as Record<string, any>[];
  const topTokens = bundle.tokenization.top_tokens_by_group as Record<string, any>[];
  const cdsQc = bundle.metrics.cds_translation_qc_summary as Record<string, any>[];
  const mvpSequenceCount = Math.max(...sequenceMetricSpecs.map((spec) => gcRows.filter((row) => row.metric === spec.id).reduce((sum, row) => sum + Number(row.n ?? 0), 0)));
  const refinedCdsCount = cdsQc.reduce((sum, row) => sum + Number(row.n_refined_sequences ?? 0), 0);
  const gc = metricRange(gcRows, "gc_content");
  const cpg = metricRange(gcRows, "cpg_oe");
  const upa = metricRange(gcRows, "upa_oe");
  const allCpgBelowOne = gcRows.filter((row) => row.metric === "cpg_oe").every((row) => Number(row.mean) < 1);
  const allUpaBelowOne = gcRows.filter((row) => row.metric === "upa_oe").every((row) => Number(row.mean) < 1);
  const separationLeader = ranking.reduce((best, row) => (Number(row.mean_js_distance ?? 0) > Number(best?.mean_js_distance ?? -1) ? row : best), ranking[0]);
  const tokenStabilityLeader = ranking.reduce((best, row) => (Number(row.mean_top_token_jaccard ?? 0) > Number(best?.mean_top_token_jaccard ?? -1) ? row : best), ranking[0]);
  const entropyRows = entropy.slice(0, 8).map((row) => ({
    tokenizer: friendlyTokenizer(row.tokenizer),
    group: row.protein_subtype,
    entropy_bits: row.mean_entropy_bits,
    effective_vocab: row.mean_effective_vocab_size
  }));
  const rankingRows = ranking.slice(0, 6).map((row) => ({
    rank: row.rank,
    tokenizer: friendlyTokenizer(row.tokenizer),
    score: row.robustness_score,
    js_distance: row.mean_js_distance,
    top_token_overlap: row.mean_top_token_jaccard
  }));
  const gcPlot = sequenceMetricSpecs.map((spec) => ({
    type: "bar",
    name: spec.label,
    x: orderedGroups,
    y: orderedGroups.map((group) => {
      const row = gcRows.find((candidate) => candidate.metric === spec.id && metricGroup(candidate) === group);
      return Number(row?.mean ?? 0);
    }),
    marker: { color: spec.color },
    hovertemplate: `${spec.label}<br>%{x}<br>mean=%{y:.3f}<extra>${spec.definition}</extra>`
  }));

  return (
    <div>
      <SectionTitle kicker="SEQUENCE / TOKEN INSPECTOR" title="Sequence context, translated">
        Real aggregate metrics, no sequences. Descriptive only.
      </SectionTitle>
      <div className="grid gap-3 md:grid-cols-2 xl:grid-cols-4">
        <Card label="MVP sequences" value={formatNumber(mvpSequenceCount, 0)} detail="HA/NA records in the sequence-context audit" />
        <Card label="Refined CDS sequences" value={formatNumber(refinedCdsCount, 0)} detail="Used only for codon-aware token summaries" />
        <Card label="Robustness leader" value={friendlyTokenizer(ranking[0]?.tokenizer)} detail={`Composite score ${formatNumber(ranking[0]?.robustness_score, 3)}`} />
        <Card label="Safe top tokens" value={formatNumber(topTokens.length, 0)} detail="Only short tokens, length <= 6" />
      </div>

      <div className="mt-4 grid gap-3 xl:grid-cols-3">
        <div className="rounded-lg border border-line bg-panel/75 p-4">
          <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">GC window</div>
          <p className="mt-2 text-sm leading-6 text-muted">
            The four MVP groups sit in a narrow GC range: <span className="text-ivory">{formatNumber(gc.min * 100, 1)}%</span> to{" "}
            <span className="text-ivory">{formatNumber(gc.max * 100, 1)}%</span>.
          </p>
        </div>
        <div className="rounded-lg border border-line bg-panel/75 p-4">
          <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">CpG readout</div>
          <p className="mt-2 text-sm leading-6 text-muted">
            CpG O/E stays below 1 across groups, ranging from <span className="text-ivory">{formatNumber(cpg.min, 3)}</span> to{" "}
            <span className="text-ivory">{formatNumber(cpg.max, 3)}</span>.
          </p>
        </div>
        <div className="rounded-lg border border-line bg-panel/75 p-4">
          <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">UpA readout</div>
          <p className="mt-2 text-sm leading-6 text-muted">
            UpA O/E also stays below 1 across groups, ranging from <span className="text-ivory">{formatNumber(upa.min, 3)}</span> to{" "}
            <span className="text-ivory">{formatNumber(upa.max, 3)}</span>.
          </p>
        </div>
      </div>

      <div className="mt-5 grid gap-4 xl:grid-cols-[1.05fr_0.95fr]">
        <div className="space-y-4">
          <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
            <Plot
              data={gcPlot}
              layout={{
                ...plotLayout,
                barmode: "group",
                title: "MVP sequence-context metrics by HA/NA group",
                height: 390,
                yaxis: { ...plotLayout.yaxis, title: "Mean value", range: [0, 1] },
                xaxis: { ...plotLayout.xaxis, tickangle: 0 }
              }}
              config={{ responsive: true, displaylogo: false }}
              useResizeHandler
              style={{ width: "100%", height: "390px" }}
            />
          </div>

          <div className="grid gap-3 md:grid-cols-3">
            {sequenceMetricSpecs.map((spec) => (
              <div key={spec.id} className="rounded-lg border border-line bg-panel/70 p-4">
                <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">{spec.short}</div>
                <div className="mt-2 text-sm font-semibold text-ivory">{spec.label}</div>
                <p className="mt-2 text-xs leading-5 text-muted">{spec.definition}</p>
              </div>
            ))}
          </div>

          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">READING THE RESULT</div>
            <div className="mt-3 space-y-3 text-sm leading-6 text-muted">
              <p>
                GC fraction is tightly bounded across the four MVP groups, from{" "}
                <span className="text-ivory">{formatNumber(gc.min * 100, 1)}%</span> in <span className="text-ivory">{gc.minGroup}</span> to{" "}
                <span className="text-ivory">{formatNumber(gc.max * 100, 1)}%</span> in <span className="text-ivory">{gc.maxGroup}</span>.
              </p>
              <p>
                CpG O/E ranges from <span className="text-ivory">{formatNumber(cpg.min, 3)}</span> in <span className="text-ivory">{cpg.minGroup}</span> to{" "}
                <span className="text-ivory">{formatNumber(cpg.max, 3)}</span> in <span className="text-ivory">{cpg.maxGroup}</span>
                {allCpgBelowOne ? ", so CpG is below single-base expectation in each displayed group." : "."}
              </p>
              <p>
                UpA O/E ranges from <span className="text-ivory">{formatNumber(upa.min, 3)}</span> in <span className="text-ivory">{upa.minGroup}</span> to{" "}
                <span className="text-ivory">{formatNumber(upa.max, 3)}</span> in <span className="text-ivory">{upa.maxGroup}</span>
                {allUpaBelowOne ? ", so UpA is also below expectation in each displayed group." : "."}
              </p>
              <p className="text-xs">
                These are compositional summaries, not claims about antigenicity, pathogenicity, escape, vaccine relevance, or fitness.
              </p>
            </div>
          </div>
        </div>

        <div className="space-y-4">
          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">TOKEN METRICS, PLAINLY</div>
            <div className="mt-3 grid gap-3 md:grid-cols-2">
              <div className="rounded-md border border-line bg-bg/30 p-3">
                <div className="text-sm font-semibold text-ivory">Entropy</div>
                <p className="mt-1 text-xs leading-5 text-muted">Higher entropy means token usage is more spread across the vocabulary.</p>
              </div>
              <div className="rounded-md border border-line bg-bg/30 p-3">
                <div className="text-sm font-semibold text-ivory">Effective vocabulary</div>
                <p className="mt-1 text-xs leading-5 text-muted">The equivalent number of equally common tokens needed to produce that entropy.</p>
              </div>
              <div className="rounded-md border border-line bg-bg/30 p-3">
                <div className="text-sm font-semibold text-ivory">JS distance</div>
                <p className="mt-1 text-xs leading-5 text-muted">A descriptive distance between group-level token distributions.</p>
              </div>
              <div className="rounded-md border border-line bg-bg/30 p-3">
                <div className="text-sm font-semibold text-ivory">Top-token overlap</div>
                <p className="mt-1 text-xs leading-5 text-muted">Jaccard overlap of top tokens under bootstrap. Higher means more stable top-token lists.</p>
              </div>
            </div>
          </div>

          <div>
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">TOKEN ENTROPY SNAPSHOT</div>
            <MiniTable rows={entropyRows} columns={["tokenizer", "group", "entropy_bits", "effective_vocab"]} limit={8} />
          </div>

          <div className="rounded-lg border border-line bg-panel/75 p-4">
            <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">ROBUSTNESS INTERPRETATION</div>
            <p className="mt-3 text-sm leading-6 text-muted">
              The composite leader is <span className="text-ivory">{friendlyTokenizer(ranking[0]?.tokenizer)}</span>. It balances coverage, low bootstrap variance and stable top-token lists.
            </p>
            <p className="mt-2 text-sm leading-6 text-muted">
              The largest mean JS distance here is <span className="text-ivory">{friendlyTokenizer(separationLeader?.tokenizer)}</span>{" "}
              ({formatNumber(separationLeader?.mean_js_distance, 3)}), while the most stable top-token overlap is{" "}
              <span className="text-ivory">{friendlyTokenizer(tokenStabilityLeader?.tokenizer)}</span> ({formatNumber(tokenStabilityLeader?.mean_top_token_jaccard, 3)}).
            </p>
            <div className="mt-4">
              <MiniTable rows={rankingRows} columns={["rank", "tokenizer", "score", "js_distance", "top_token_overlap"]} limit={6} />
            </div>
          </div>
        </div>
      </div>
    </div>
  );
}

function StructureView({ bundle }: { bundle: SafeBundle }) {
  const structures = bundle.structures.structures as Record<string, any>[];
  const mapping = bundle.structureMapping ?? {};
  const mappingQc = (mapping.mapping_qc ?? []) as Record<string, any>[];
  const signalCatalog = (mapping.signal_catalog ?? []) as Record<string, any>[];
  const [pdbId, setPdbId] = useState(structures[0]?.pdb_id ?? "");
  const structure = structures.find((item) => item.pdb_id === pdbId) ?? structures[0];
  const qcRows = mappingQc.filter((row) => row.pdb_id === pdbId);
  const bestQc = qcRows.reduce((best, row) => (Number(row.mapped_residues ?? 0) > Number(best?.mapped_residues ?? -1) ? row : best), qcRows[0]);
  const catalogRow = signalCatalog.find((row) => row.protein === structure?.protein && row.subtype === structure?.subtype_context);

  return (
    <div>
      <SectionTitle kicker="STRUCTURE VIEW" title="Public structures with alignment QC">
        Public RCSB structures plus a derived residue-signal bridge. Metric coloring still waits for chain-number validation.
      </SectionTitle>
      <div className="mb-4 grid gap-3 xl:grid-cols-[0.8fr_1.1fr_1.1fr]">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Structure</span>
          <select value={pdbId} onChange={(event) => setPdbId(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {structures.map((item) => (
              <option key={item.pdb_id} value={item.pdb_id}>
                {item.pdb_id} · {item.label}
              </option>
            ))}
          </select>
          <div className="mt-3 text-xs leading-5 text-muted">
            {structure?.protein} · {structure?.subtype_context} context ·{" "}
            <a className="text-teal underline-offset-4 hover:underline" href={structure?.rcsb_url} target="_blank" rel="noreferrer">
              open RCSB
            </a>
          </div>
        </label>
        <div className="rounded-lg border border-line bg-panel/70 p-4 text-sm leading-6 text-muted">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-brass">WHAT IS LOADED</span>
          <p className="mt-2">
            {structure?.label}. The viewer loads the public coordinate file from RCSB and lets you rotate, zoom and switch styles. This is safe to show because it is not a restricted sequence artifact.
          </p>
        </div>
        <div className="rounded-lg border border-line bg-panel/70 p-4 text-sm leading-6 text-muted">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-brass">MAPPING STATUS</span>
          <p className="mt-2">
            Alignment QC is now available: local refined CDS positions were aligned to public PDB polymer sequences. Residue coloring remains pending until chain IDs, residue numbering and missing residues are validated.
          </p>
        </div>
      </div>

      <div className="mb-4 grid gap-3 md:grid-cols-2 xl:grid-cols-4">
        <Card label="Best mapped residues" value={formatNumber(bestQc?.mapped_residues ?? 0, 0)} detail={bestQc?.pdb_entity ? `${bestQc.pdb_entity} · ${bestQc.chains}` : "Public polymer sequence"} />
        <Card label="PDB coverage" value={bestQc?.coverage_pdb == null ? "NA" : `${formatNumber(Number(bestQc.coverage_pdb) * 100, 1)}%`} detail="Sequence alignment coverage, not final atom coloring" />
        <Card label="Mean CpG codon signal" value={formatNumber(catalogRow?.mean_cpg_codon_fraction ?? "NA", 3)} detail="Aggregate over refined CDS residues" />
        <Card label="Mean UpA codon signal" value={formatNumber(catalogRow?.mean_upa_codon_fraction ?? "NA", 3)} detail="DNA TA proxy for RNA UpA" />
      </div>

      <div className="mb-4 grid gap-3 xl:grid-cols-4">
        {[
          ["1", "Align sequence to structure", "Match each local HA/NA sequence position to the reference PDB sequence."],
          ["2", "Resolve residue numbering", "Handle PDB chain IDs, insertion codes, missing residues and subtype-specific numbering."],
          ["3", "Choose mapped metric", "Define whether a residue receives codon, token, entropy or group-level summaries."],
          ["4", "Validate before coloring", "Only then should the viewer paint FluGenome3D-derived values onto the structure."]
        ].map(([step, title, detail]) => (
          <div key={step} className="rounded-lg border border-line bg-panel/70 p-4">
            <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-teal">Mapping step {step}</div>
            <div className="mt-2 text-sm font-semibold text-ivory">{title}</div>
            <p className="mt-2 text-xs leading-5 text-muted">{detail}</p>
          </div>
        ))}
      </div>

      {structure ? <MolecularViewer structure={structure as any} /> : null}

      <div className="mt-4 grid gap-4 xl:grid-cols-[1.1fr_0.9fr]">
        <div>
          <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">ALIGNMENT QC FOR SELECTED STRUCTURE</div>
          <MiniTable rows={qcRows} columns={["pdb_entity", "chains", "pdb_sequence_length", "mapped_residues", "identity", "coverage_pdb", "local_start", "local_end"]} limit={6} />
        </div>
        <div>
          <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-brass">RESIDUE SIGNAL CATALOG</div>
          <MiniTable rows={signalCatalog} columns={["group", "n_local_positions", "n_mapped_positions", "mean_gc_fraction_codon", "mean_aa_entropy"]} limit={6} />
        </div>
      </div>
    </div>
  );
}

function BridgeView({ bundle }: { bundle: SafeBundle }) {
  const reps = bundle.representations.representations as Record<string, any>[];
  const structures = bundle.structures.structures as Record<string, any>[];
  const [group, setGroup] = useState("HA-H1N1");
  const [repId, setRepId] = useState(reps[0]?.id ?? "");
  const rep = reps.find((item) => item.id === repId) ?? reps[0];
  const points = decodeRepresentationPoints(rep).filter((point) => point.group === group);
  const [protein, subtype] = group.split("-");
  const structure = structures.find((item) => item.protein === protein && item.subtype_context === subtype) ?? structures[0];
  const metricRows = (bundle.metrics.gc_cpg_upa_summary as Record<string, any>[]).filter((row) => `${row.protein}-${row.subtype}` === group);
  const stabilityRows = (bundle.stability.tokenizer_robustness_ranking as Record<string, any>[]).slice(0, 4);
  const latentCache = bundle.antigenlm?.cache_summary ?? {};

  return (
    <div>
      <SectionTitle kicker="BRIDGE MODE" title="Sequence context to representation to structure">
        Integrated view of safe derived artifacts.
      </SectionTitle>
      <div className="mb-4 grid gap-3 md:grid-cols-2">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Group</span>
          <select value={group} onChange={(event) => setGroup(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {["HA-H1N1", "NA-H1N1", "HA-H3N2", "NA-H3N2"].map((value) => (
              <option key={value} value={value}>
                {value}
              </option>
            ))}
          </select>
        </label>
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Representation</span>
          <select value={repId} onChange={(event) => setRepId(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {reps.map((item) => (
              <option key={item.id} value={item.id}>
                {item.label}
              </option>
            ))}
          </select>
        </label>
      </div>
      <div className="grid gap-4 xl:grid-cols-[1fr_1fr]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={[
              {
                type: "scattergl",
                mode: "markers",
                name: group,
                x: points.map((point) => point.x),
                y: points.map((point) => point.y),
                marker: { color: groupPalette[group], size: 6, opacity: 0.78 },
                text: points.map((point) => `${point.id} · ${point.year_bin}`)
              }
            ]}
            layout={{ ...plotLayout, title: `${group} mini projector`, height: 360 }}
            config={{ responsive: true, displaylogo: false }}
            useResizeHandler
            style={{ width: "100%", height: "360px" }}
          />
          <div className="mt-3 grid gap-3 md:grid-cols-2">
            <MiniTable rows={metricRows} columns={["metric", "mean", "median", "q05", "q95"]} limit={5} />
            <MiniTable rows={stabilityRows} columns={["rank", "tokenizer", "robustness_score"]} limit={4} />
          </div>
        </div>
        <div className="rounded-lg border border-line bg-panel/75 p-4">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-brass">ASSOCIATED STRUCTURE</div>
          <div className="mt-3 text-xl font-semibold text-ivory">{structure?.pdb_id}</div>
          <p className="mt-2 text-sm leading-6 text-muted">
            {structure?.label}. Mapping status is {structure?.mapping_status}; this bridge shows alignment QC but does not claim residue-level functional interpretation.
          </p>
          <div className="mt-5 rounded-lg border border-line bg-ink p-4 text-xs leading-6 text-muted">
            <ShieldCheck className="mb-3 text-brass" size={18} />
            Data are real derived research artifacts. Raw sequences, FASTA, accessions, isolate names and restricted Parquet panels are not redistributed.
          </div>
          <div className="mt-4 rounded-lg border border-line bg-ink p-4 text-xs leading-6 text-muted">
            <span className="font-mono uppercase tracking-[0.18em] text-teal">AntigenLM layer</span>
            <p className="mt-2">
              Parent latent cache represented {formatNumber(latentCache.n_records ?? 0, 0)} HA+NA records. The learned layer is used as descriptive geometry, not prediction.
            </p>
          </div>
        </div>
      </div>
    </div>
  );
}
