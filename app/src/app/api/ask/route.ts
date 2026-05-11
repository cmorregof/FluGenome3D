import { promises as fs } from "fs";
import path from "path";
import { NextResponse } from "next/server";

type GuideChunk = {
  id: string;
  title: string;
  source: string;
  section: string;
  topic_tags?: string[];
  text: string;
};

type GuidePayload = {
  guide_policy?: string;
  answer_boundary?: Record<string, string[]>;
  suggested_questions?: string[];
  chunks?: GuideChunk[];
};

const STOP_WORDS = new Set([
  "a",
  "al",
  "and",
  "are",
  "como",
  "con",
  "de",
  "del",
  "does",
  "el",
  "en",
  "es",
  "este",
  "for",
  "is",
  "it",
  "la",
  "las",
  "lo",
  "los",
  "me",
  "mean",
  "para",
  "por",
  "que",
  "se",
  "the",
  "this",
  "to",
  "un",
  "una",
  "what",
  "why",
  "y"
]);

const ALIASES: Record<string, string[]> = {
  ai: ["guide", "antigenlm", "learned", "model"],
  antigenlm: ["antigenlm", "latent", "embedding", "learned"],
  bpe: ["bpe", "tokenizer", "future"],
  cpg: ["cpg", "observed", "expected", "dinucleotide"],
  gc: ["gc", "composition", "sequence"],
  grover: ["grover", "tokenizer", "future"],
  js: ["jensen", "shannon", "distance", "token"],
  pca: ["pca", "projection", "representation", "coordinate"],
  rscu: ["rscu", "codon", "cds"],
  sne: ["tsne", "projection", "neighborhood", "latent"],
  structure: ["structure", "pdb", "rcsb", "mapping", "residue"],
  tsne: ["tsne", "t", "sne", "projection", "neighborhood", "latent"],
  "t-sne": ["tsne", "projection", "neighborhood", "latent"],
  upa: ["upa", "ta", "observed", "expected"],
};

const PROHIBITED_SCOPE_TERMS = [
  "antigenic",
  "antigenicity",
  "escape",
  "fitness",
  "pathogenicity",
  "pathogenic",
  "transmissibility",
  "vaccine",
  "vacuna",
  "virulence",
  "predict",
  "prediction",
  "predecir",
  "causal",
  "causality",
  "optimize",
  "optimizar",
];

function normalizeText(text: string) {
  return text
    .normalize("NFD")
    .replace(/[\u0300-\u036f]/g, "")
    .toLowerCase();
}

function tokenize(text: string): string[] {
  return normalizeText(text)
    .split(/[^a-z0-9]+/)
    .filter((token) => token.length > 1 && !STOP_WORDS.has(token));
}

function expandedQueryTokens(question: string): Set<string> {
  const base = new Set(tokenize(question));
  for (const token of Array.from(base)) {
    for (const alias of ALIASES[token] ?? []) base.add(alias);
  }
  return base;
}

function chunkTerms(chunk: GuideChunk) {
  return new Set(tokenize(`${chunk.title} ${chunk.section} ${chunk.topic_tags?.join(" ") ?? ""} ${chunk.text}`));
}

function scoreChunk(chunk: GuideChunk, query: Set<string>) {
  const terms = chunkTerms(chunk);
  const titleTerms = new Set(tokenize(`${chunk.title} ${chunk.section}`));
  const tagTerms = new Set(tokenize(chunk.topic_tags?.join(" ") ?? ""));
  let score = 0;
  for (const token of query) {
    if (terms.has(token)) score += 1;
    if (titleTerms.has(token)) score += 3;
    if (tagTerms.has(token)) score += 4;
  }
  return score;
}

function firstReadableSentence(text: string) {
  const clean = text.replace(/\s+/g, " ").trim();
  const match = clean.match(/^(.{80,360}?[.!?])\s/);
  return (match?.[1] ?? clean.slice(0, 360)).trim();
}

function scopeBoundary(question: string) {
  const normal = normalizeText(question);
  return PROHIBITED_SCOPE_TERMS.some((term) => normal.includes(term));
}

function sourceLabel(chunk: GuideChunk, index: number) {
  return `[${index + 1}] ${chunk.title} (${chunk.source})`;
}

function viewHint(chunks: GuideChunk[]) {
  const tags = new Set(chunks.flatMap((chunk) => chunk.topic_tags ?? []));
  if (tags.has("structure")) return "Open 3D Molecular Viewer or Bridge View for the structure-side context.";
  if (tags.has("antigenlm")) return "Open AntigenLM Latent Atlas for the learned representation layer.";
  if (tags.has("dataset")) return "Open Dataset Atlas for panel design, geography and coverage.";
  if (tags.has("tokenization") || tags.has("sequence_context")) return "Open Sequence/Token Inspector for metrics, entropy and token stability.";
  if (tags.has("representation")) return "Open Representation Projector to inspect reduced-coordinate maps.";
  return "Open Project Guide for the plain-language map of the project.";
}

function composeAnswer(question: string, chunks: GuideChunk[], guide: GuidePayload) {
  const boundary = scopeBoundary(question);
  const top = chunks.slice(0, 4);
  const snippets = top.map((chunk, index) => `- ${firstReadableSentence(chunk.text)} ${sourceLabel(chunk, index)}`);
  const opening = boundary
    ? "Short answer: FluGenome3D should stay in descriptive mode for that question. It can explain the representation, tokenization, sequence-context or structure-QC evidence, but it should not turn those summaries into antigenicity, vaccine, escape, fitness or pathogenicity claims."
    : `Short answer: ${firstReadableSentence(top[0]?.text ?? guide.guide_policy ?? "FluGenome3D is a descriptive visual lab built from safe derived research artifacts.")}`;

  return [
    opening,
    "",
    "How to read it:",
    snippets.join("\n"),
    "",
    "Where to look in the app:",
    `- ${viewHint(top)}`,
    "",
    "Boundary:",
    "- The guide is grounded in safe docs, reports and exported summaries. It does not access raw sequences, FASTA, restricted Parquet files or external model APIs."
  ].join("\n");
}

async function readGuide(): Promise<GuidePayload> {
  const filePath = path.join(process.cwd(), "data", "lab_guide.safe.json");
  const text = await fs.readFile(filePath, "utf-8");
  return JSON.parse(text) as GuidePayload;
}

export async function POST(request: Request) {
  let question = "";
  try {
    const payload = (await request.json()) as { question?: unknown };
    question = String(payload.question ?? "").trim();
  } catch {
    return NextResponse.json({ error: "Expected JSON payload with a question field." }, { status: 400 });
  }

  if (question.length < 4) {
    return NextResponse.json({ error: "Ask a slightly longer question." }, { status: 400 });
  }

  const guide = await readGuide();
  const chunks = guide.chunks ?? [];
  const query = expandedQueryTokens(question);
  const ranked = chunks
    .map((chunk) => ({ chunk, score: scoreChunk(chunk, query) }))
    .filter((item) => item.score > 0)
    .sort((a, b) => b.score - a.score)
    .map((item) => item.chunk);
  const selected = (ranked.length ? ranked : chunks).slice(0, 4);

  return NextResponse.json({
    answer: composeAnswer(question, selected, guide),
    citations: selected.map((chunk, index) => ({
      marker: `[${index + 1}]`,
      id: chunk.id,
      title: chunk.title,
      source: chunk.source,
      section: chunk.section,
      topic_tags: chunk.topic_tags ?? [],
    })),
    matched_topics: Array.from(new Set(selected.flatMap((chunk) => chunk.topic_tags ?? []))).slice(0, 8),
    guardrails: [
      "Grounded in safe derived artifacts only.",
      "No raw sequences or restricted panels are available to this route.",
      "Descriptive interpretation only; no antigenicity, vaccine, escape, fitness or pathogenicity prediction.",
    ],
  });
}

export async function GET() {
  const guide = await readGuide();
  return NextResponse.json({
    guide_policy: guide.guide_policy,
    suggested_questions: guide.suggested_questions ?? [],
    n_chunks: guide.chunks?.length ?? 0,
  });
}
