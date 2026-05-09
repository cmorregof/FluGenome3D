"use client";

import dynamic from "next/dynamic";
import { useEffect, useMemo, useState } from "react";
import {
  Activity,
  Atom,
  Boxes,
  Database,
  Dna,
  Network,
  ShieldCheck,
  SplitSquareVertical
} from "lucide-react";
import MolecularViewer from "@/components/MolecularViewer";
import { SafeBundle, formatNumber, groupKey, loadSafeBundle, uniqueValues } from "@/lib/safe-data";

const Plot = dynamic(() => import("react-plotly.js"), { ssr: false }) as any;

const palette = ["#d7a84a", "#7fa69a", "#b87554", "#eee4cf", "#a58bd6", "#7aa0c4"];
const groupPalette: Record<string, string> = {
  "HA-H1N1": "#d7a84a",
  "NA-H1N1": "#7fa69a",
  "HA-H3N2": "#b87554",
  "NA-H3N2": "#eee4cf"
};

const views = [
  { id: "atlas", label: "Dataset Atlas", icon: Database },
  { id: "projector", label: "Representation Projector", icon: SplitSquareVertical },
  { id: "inspector", label: "Sequence/Token Inspector", icon: Dna },
  { id: "structure", label: "3D Molecular Viewer", icon: Atom },
  { id: "bridge", label: "Bridge View", icon: Network }
] as const;

type ViewId = (typeof views)[number]["id"];

const plotLayout = {
  paper_bgcolor: "rgba(0,0,0,0)",
  plot_bgcolor: "rgba(17,16,13,0.96)",
  font: { color: "#eee4cf", family: "IBM Plex Mono, ui-monospace, monospace", size: 11 },
  margin: { l: 48, r: 24, t: 40, b: 48 },
  xaxis: { gridcolor: "rgba(238,228,207,0.12)", zerolinecolor: "rgba(238,228,207,0.2)" },
  yaxis: { gridcolor: "rgba(238,228,207,0.12)", zerolinecolor: "rgba(238,228,207,0.2)" },
  legend: { bgcolor: "rgba(8,8,6,0.5)", bordercolor: "rgba(238,228,207,0.14)", borderwidth: 1 }
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
        <div className="font-mono text-[11px] uppercase tracking-[0.25em] text-amber">{kicker}</div>
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

function buildTemporalPlot(rows: Record<string, any>[]) {
  const windows = ["pre-2009", "2009-2014", "2015-2019", "2020+"];
  return ["H1N1", "H3N2"].map((subtype, index) => ({
    type: "bar",
    name: subtype,
    x: windows,
    y: windows.map((window) => {
      const match = rows.find((row) => row.year_bin === window && row.subtype === subtype);
      return Number(match?.n_pairs ?? 0);
    }),
    marker: { color: palette[index] }
  }));
}

function buildScatterTraces(points: Record<string, any>[], colorBy: string) {
  const values = uniqueValues(points.map((point) => String(point[colorBy] ?? "unknown")));
  return values.map((value, index) => {
    const subset = points.filter((point) => String(point[colorBy] ?? "unknown") === value);
    return {
      type: "scattergl",
      mode: "markers",
      name: value,
      x: subset.map((point) => point.x),
      y: subset.map((point) => point.y),
      text: subset.map((point) => `${point.group} · ${point.year_bin}`),
      customdata: subset,
      marker: {
        color: groupPalette[value] ?? palette[index % palette.length],
        size: 5,
        opacity: 0.72,
        line: { color: "rgba(238,228,207,0.16)", width: 0.5 }
      }
    };
  });
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
  const [active, setActive] = useState<ViewId>("atlas");

  if (error) {
    return <div className="p-8 text-rust">{error}</div>;
  }

  if (!bundle) {
    return (
      <main className="flex min-h-screen items-center justify-center bg-ink text-ivory">
        <div className="font-mono text-sm uppercase tracking-[0.3em] text-amber">Loading FluGenome3D safe bundle</div>
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
            <div className="font-mono text-[10px] uppercase tracking-[0.28em] text-amber">SYSTEM STATUS</div>
            <h1 className="mt-3 text-2xl font-semibold text-ivory">FluGenome3D</h1>
            <p className="mt-3 text-xs leading-5 text-muted">Visual lab for real derived Influenza A HA/NA artifacts.</p>
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
                      ? "border-amber bg-amber/12 text-ivory"
                      : "border-line bg-panel/50 text-muted hover:border-mint hover:text-ivory"
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
            <div className="mt-2 text-sm text-ivory">{mode === "local-full" ? "Local full" : "Vercel safe"}</div>
          </div>
        </aside>

        <section className="min-w-0 flex-1">
          <div className="mb-4 rounded-lg border border-amber/40 bg-amber/10 px-4 py-3 text-sm text-ivory">
            <div className="flex flex-col gap-2 md:flex-row md:items-center md:justify-between">
              <span>{banner}</span>
              <span className="font-mono text-[10px] uppercase tracking-[0.22em] text-amber">
                {mode === "local-full" ? "LOCAL FULL MODE" : "VERCEL SAFE MODE"}
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
                  active === view.id ? "bg-amber text-ink" : "bg-panel text-muted"
                }`}
              >
                {view.label}
              </button>
            ))}
          </div>

          <div className="rounded-lg border border-line bg-ink/82 p-4 md:p-6">
            {active === "atlas" ? <DatasetAtlas bundle={bundle} mode={mode} /> : null}
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

function DatasetAtlas({ bundle, mode }: { bundle: SafeBundle; mode: string }) {
  const datasetRows = bundle.dataset.dataset_summary as Record<string, any>[];
  const panelRows = bundle.dataset.panel_summary as Record<string, any>[];
  const temporal = bundle.dataset.temporal_counts as Record<string, any>[];
  const refined = bundle.dataset.cds_refined_qc as Record<string, any>[];
  const allValid = datasetRows.find((row) => row.stage === "all_valid_paired");
  const mvp = datasetRows.find((row) => row.stage === "mvp_panel");
  const full = datasetRows.find((row) => row.stage === "full_panel");

  return (
    <div>
      <SectionTitle kicker="DATASET ATLAS" title="Restricted data, safe surface">
        Real derived artifacts. No raw sequences in app data.
      </SectionTitle>
      <div className="grid gap-3 md:grid-cols-2 xl:grid-cols-4">
        <Card label="Valid paired strains" value={formatNumber(allValid?.n_pairs ?? 0, 0)} detail="Phase 1 local audit" />
        <Card label="MVP sequences" value={formatNumber(mvp?.n_sequence_records ?? 0, 0)} detail="Balanced HA/NA x subtype" />
        <Card label="Full dedup pairs" value={formatNumber(full?.n_pairs ?? 0, 0)} detail="Exact HA+NA deduplicated" />
        <Card label="App mode" value={mode === "local-full" ? "Local full" : "Vercel safe"} detail="Safe JSON boundary active" />
      </div>

      <div className="mt-5 grid gap-4 xl:grid-cols-[1.35fr_1fr]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={buildTemporalPlot(temporal)}
            layout={{ ...plotLayout, barmode: "stack", title: "MVP temporal distribution by subtype", height: 390 }}
            config={{ responsive: true, displaylogo: false }}
            className="h-full w-full"
            useResizeHandler
            style={{ width: "100%", height: "390px" }}
          />
        </div>
        <div className="space-y-4">
          <div>
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-amber">PANEL COUNTS</div>
            <MiniTable rows={panelRows} columns={["panel", "subtype", "n_pairs", "n_sequence_records", "year_min", "year_max"]} />
          </div>
          <div>
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-amber">CDS RELIABILITY</div>
            <MiniTable
              rows={refined}
              columns={["subtype", "protein", "n_refined_sequences", "n_strict_pass", "n_rescued", "n_unrescued"]}
            />
          </div>
        </div>
      </div>
    </div>
  );
}

function RepresentationProjector({ bundle }: { bundle: SafeBundle }) {
  const reps = bundle.representations.representations as Record<string, any>[];
  const [repId, setRepId] = useState(reps[0]?.id ?? "");
  const [colorBy, setColorBy] = useState("group");
  const [selected, setSelected] = useState<Record<string, any> | null>(null);
  const rep = reps.find((item) => item.id === repId) ?? reps[0];
  const points = (rep?.points ?? []) as Record<string, any>[];

  const traces = useMemo(() => buildScatterTraces(points, colorBy), [points, colorBy]);

  return (
    <div>
      <SectionTitle kicker="REPRESENTATION PROJECTOR" title="Reduced-coordinate maps">
        Coordinates are real derived artifacts with hashed IDs.
      </SectionTitle>
      <div className="mb-4 grid gap-3 md:grid-cols-3">
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
        <Card label="Exported points" value={formatNumber(rep?.n_exported_points ?? 0, 0)} detail={rep?.privacy ?? "safe point metadata"} />
      </div>
      <div className="grid gap-4 xl:grid-cols-[1fr_300px]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={traces}
            layout={{ ...plotLayout, title: rep?.label ?? "Representation", height: 610 }}
            config={{ responsive: true, displaylogo: false }}
            useResizeHandler
            style={{ width: "100%", height: "610px" }}
            onClick={(event: any) => setSelected(event.points?.[0]?.customdata ?? null)}
          />
        </div>
        <div className="rounded-lg border border-line bg-panel/80 p-4">
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-amber">POINT / MAP INSPECTOR</div>
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
          <div className="mt-6">
            <MiniTable rows={rep?.silhouette_scores ?? []} columns={["label_type", "silhouette", "space"]} limit={4} />
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
  const gcPlot = ["gc_content", "cpg_oe", "upa_oe"].map((metric, index) => ({
    type: "bar",
    name: metric,
    x: gcRows.filter((row) => row.metric === metric).map(groupKey),
    y: gcRows.filter((row) => row.metric === metric).map((row) => Number(row.mean)),
    marker: { color: palette[index] }
  }));

  return (
    <div>
      <SectionTitle kicker="SEQUENCE CONTEXT" title="Metrics without sequences">
        Aggregate GC/CpG/UpA, token entropy and robustness summaries.
      </SectionTitle>
      <div className="grid gap-3 md:grid-cols-3">
        <Card label="Safe top tokens" value={formatNumber(topTokens.length, 0)} detail="All token strings length <= 6" />
        <Card label="Robustness top" value={ranking[0]?.tokenizer ?? "NA"} detail={`Score ${formatNumber(ranking[0]?.robustness_score, 3)}`} />
        <Card label="Tokenizers audited" value={formatNumber((bundle.tokenization.tokenizer_summary ?? []).length, 0)} detail="Deterministic only" />
      </div>
      <div className="mt-5 grid gap-4 xl:grid-cols-[1.1fr_0.9fr]">
        <div className="plot-shell rounded-lg border border-line bg-panel/75 p-3">
          <Plot
            data={gcPlot}
            layout={{ ...plotLayout, barmode: "group", title: "MVP sequence-context metrics by group", height: 420 }}
            config={{ responsive: true, displaylogo: false }}
            useResizeHandler
            style={{ width: "100%", height: "420px" }}
          />
        </div>
        <div className="space-y-4">
          <div>
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-amber">TOKEN ENTROPY</div>
            <MiniTable rows={entropy} columns={["tokenizer", "protein_subtype", "mean_entropy_bits", "mean_effective_vocab_size"]} limit={8} />
          </div>
          <div>
            <div className="mb-2 font-mono text-[11px] uppercase tracking-[0.22em] text-amber">ROBUSTNESS RANKING</div>
            <MiniTable rows={ranking} columns={["rank", "tokenizer", "robustness_score", "mean_js_distance", "mean_top_token_jaccard"]} limit={8} />
          </div>
        </div>
      </div>
    </div>
  );
}

function StructureView({ bundle }: { bundle: SafeBundle }) {
  const structures = bundle.structures.structures as Record<string, any>[];
  const [pdbId, setPdbId] = useState(structures[0]?.pdb_id ?? "");
  const structure = structures.find((item) => item.pdb_id === pdbId) ?? structures[0];

  return (
    <div>
      <SectionTitle kicker="STRUCTURE VIEW" title="Public PDB viewer">
        Mapping is pending unless explicitly validated.
      </SectionTitle>
      <div className="mb-4 grid gap-3 md:grid-cols-[1fr_2fr]">
        <label className="rounded-lg border border-line bg-panel/70 p-3 text-sm">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-muted">Structure</span>
          <select value={pdbId} onChange={(event) => setPdbId(event.target.value)} className="mt-2 w-full rounded-md border border-line bg-ink px-3 py-2 text-ivory">
            {structures.map((item) => (
              <option key={item.pdb_id} value={item.pdb_id}>
                {item.pdb_id} · {item.label}
              </option>
            ))}
          </select>
        </label>
        <div className="rounded-lg border border-line bg-panel/70 p-4 text-sm leading-6 text-muted">
          <span className="font-mono text-[10px] uppercase tracking-[0.2em] text-amber">MAPPING STATUS</span>
          <p className="mt-2">
            {structure?.label}. The viewer loads public RCSB coordinates. Sequence-to-structure metric mapping is marked{" "}
            <span className="text-ivory">{structure?.mapping_status}</span>.
          </p>
        </div>
      </div>
      {structure ? <MolecularViewer structure={structure as any} /> : null}
    </div>
  );
}

function BridgeView({ bundle }: { bundle: SafeBundle }) {
  const reps = bundle.representations.representations as Record<string, any>[];
  const structures = bundle.structures.structures as Record<string, any>[];
  const [group, setGroup] = useState("HA-H1N1");
  const [repId, setRepId] = useState(reps[0]?.id ?? "");
  const rep = reps.find((item) => item.id === repId) ?? reps[0];
  const points = ((rep?.points ?? []) as Record<string, any>[]).filter((point) => point.group === group);
  const [protein, subtype] = group.split("-");
  const structure = structures.find((item) => item.protein === protein && item.subtype_context === subtype) ?? structures[0];
  const metricRows = (bundle.metrics.gc_cpg_upa_summary as Record<string, any>[]).filter((row) => `${row.protein}-${row.subtype}` === group);
  const stabilityRows = (bundle.stability.tokenizer_robustness_ranking as Record<string, any>[]).slice(0, 4);

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
          <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-amber">ASSOCIATED STRUCTURE</div>
          <div className="mt-3 text-xl font-semibold text-ivory">{structure?.pdb_id}</div>
          <p className="mt-2 text-sm leading-6 text-muted">
            {structure?.label}. Mapping status is {structure?.mapping_status}; this bridge does not claim residue-level correspondence.
          </p>
          <div className="mt-5 rounded-lg border border-line bg-ink p-4 text-xs leading-6 text-muted">
            <ShieldCheck className="mb-3 text-amber" size={18} />
            Data are real derived research artifacts. Raw sequences, FASTA, accessions, isolate names and restricted Parquet panels are not redistributed.
          </div>
        </div>
      </div>
    </div>
  );
}
