"use client";

import { RotateCcw } from "lucide-react";
import { useCallback, useEffect, useRef, useState } from "react";

type StructureRecord = {
  pdb_id: string;
  label: string;
  pdb_download_url: string;
  rcsb_url: string;
  mapping_status: string;
  protein: string;
  subtype_context: string;
};

type ViewerState = "idle" | "loading" | "ready" | "error";
type ViewerStyle = "cartoon" | "surface" | "stick";
type ColorMode = "default" | "chain" | "pending";

function controlClass(active: boolean) {
  return `rounded-md border px-3 py-2 text-xs transition ${
    active
      ? "border-lineStrong bg-brassSoft/45 text-ivory"
      : "border-line bg-panel2 text-muted hover:border-teal hover:text-ivory"
  }`;
}

export default function MolecularViewer({ structure }: { structure: StructureRecord }) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const viewerRef = useRef<any>(null);
  const threeDmolRef = useRef<any>(null);
  const [viewerState, setViewerState] = useState<ViewerState>("idle");
  const [message, setMessage] = useState("Ready to load public RCSB coordinates.");
  const [styleMode, setStyleMode] = useState<ViewerStyle>("cartoon");
  const [colorMode, setColorMode] = useState<ColorMode>("default");

  const cleanupViewer = useCallback(() => {
    const viewer = viewerRef.current;
    if (viewer) {
      try {
        viewer.spin(false);
        viewer.clear?.();
        viewer.removeAllModels?.();
        viewer.removeAllSurfaces?.();
      } catch {
        // 3Dmol cleanup APIs vary by version; clearing the container below is the hard boundary.
      }
    }
    viewerRef.current = null;
    if (containerRef.current) {
      containerRef.current.innerHTML = "";
    }
  }, []);

  const applyStyle = useCallback((viewer: any, threeDmol: any) => {
    viewer.setStyle({}, {});
    viewer.removeAllSurfaces?.();

    const chainStyle = colorMode === "chain" ? { colorscheme: "chainHetatm" } : {};
    const pendingColor = colorMode === "pending" ? "#87a99b" : undefined;

    if (styleMode === "cartoon") {
      viewer.setStyle(
        {},
        {
          cartoon: {
            opacity: 0.94,
            ...(pendingColor ? { color: pendingColor } : { color: "spectrum" }),
            ...chainStyle
          }
        }
      );
    }

    if (styleMode === "surface") {
      viewer.setStyle({}, { cartoon: { opacity: 0.24, color: "#a8a092", ...chainStyle } });
      viewer.addSurface(threeDmol.SurfaceType.VDW, {
        opacity: 0.58,
        color: pendingColor ?? "#87a99b"
      });
    }

    if (styleMode === "stick") {
      viewer.setStyle(
        {},
        {
          stick: {
            radius: 0.14,
            ...(pendingColor ? { color: pendingColor } : chainStyle)
          }
        }
      );
    }

    viewer.render();
  }, [colorMode, styleMode]);

  const resetCamera = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    viewer.zoomTo();
    viewer.render();
  }, []);

  useEffect(() => {
    let cancelled = false;

    async function loadStructure() {
      cleanupViewer();
      if (!containerRef.current) return;
      setViewerState("loading");
      setMessage("Loading RCSB structure...");

      try {
        const mod = await import("3dmol");
        const threeDmol = (mod as any).default ?? mod;
        threeDmolRef.current = threeDmol;

        const viewer = threeDmol.createViewer(containerRef.current, {
          backgroundColor: "#080807"
        });
        viewerRef.current = viewer;

        const response = await fetch(structure.pdb_download_url);
        if (!response.ok) throw new Error(`RCSB returned ${response.status}`);
        const pdb = await response.text();
        if (cancelled) return;

        viewer.addModel(pdb, "pdb");
        applyStyle(viewer, threeDmol);
        viewer.zoomTo();
        viewer.spin(true);
        viewer.render();
        setViewerState("ready");
        setMessage(`${structure.pdb_id} loaded from RCSB. Mapping status: ${structure.mapping_status}.`);
      } catch (error) {
        if (cancelled) return;
        cleanupViewer();
        setViewerState("error");
        setMessage(error instanceof Error ? error.message : "Could not load structure.");
      }
    }

    loadStructure();
    return () => {
      cancelled = true;
      cleanupViewer();
    };
  }, [cleanupViewer, structure]);

  useEffect(() => {
    const viewer = viewerRef.current;
    const threeDmol = threeDmolRef.current;
    if (!viewer || !threeDmol || viewerState !== "ready") return;
    applyStyle(viewer, threeDmol);
  }, [applyStyle, viewerState]);

  useEffect(() => {
    const node = containerRef.current;
    if (!node) return;
    const observer = new ResizeObserver(() => {
      const viewer = viewerRef.current;
      if (!viewer) return;
      try {
        viewer.resize();
        viewer.render();
      } catch {
        // Resize is an enhancement; failure should not break the safe fallback.
      }
    });
    observer.observe(node);
    return () => observer.disconnect();
  }, []);

  const fallback = viewerState === "error";

  return (
    <div className="overflow-hidden rounded-lg border border-line bg-ink shadow-[0_0_36px_rgba(0,0,0,0.28)]">
      <div className="flex flex-col gap-3 border-b border-line bg-panel2 p-3 xl:flex-row xl:items-center xl:justify-between">
        <div>
          <div className="font-mono text-[10px] uppercase tracking-[0.22em] text-brass">PUBLIC STRUCTURE VIEWER</div>
          <div className="mt-1 text-sm text-muted">
            Structure viewer loads public RCSB coordinates. Mapping status: pending unless explicitly validated.
          </div>
        </div>
        <div className="flex flex-wrap gap-2">
          {(["cartoon", "surface", "stick"] as ViewerStyle[]).map((mode) => (
            <button key={mode} onClick={() => setStyleMode(mode)} className={controlClass(styleMode === mode)}>
              {mode}
            </button>
          ))}
          {(["default", "chain", "pending"] as ColorMode[]).map((mode) => (
            <button key={mode} onClick={() => setColorMode(mode)} className={controlClass(colorMode === mode)}>
              {mode}
            </button>
          ))}
          <button onClick={resetCamera} className="inline-flex items-center gap-2 rounded-md border border-line bg-panel2 px-3 py-2 text-xs text-muted transition hover:border-teal hover:text-ivory">
            <RotateCcw size={13} />
            reset
          </button>
        </div>
      </div>

      <div className="relative h-[520px] w-full overflow-hidden bg-bg">
        <div ref={containerRef} className="molecular-viewer absolute inset-0 h-full w-full overflow-hidden" />

        {viewerState === "loading" ? (
          <div className="absolute inset-0 flex items-center justify-center bg-bg/84">
            <div className="rounded-lg border border-line bg-panel px-5 py-4 text-center">
              <div className="font-mono text-[11px] uppercase tracking-[0.24em] text-brass">Loading RCSB structure...</div>
              <div className="mt-2 text-sm text-muted">{structure.pdb_id} public coordinates</div>
            </div>
          </div>
        ) : null}

        {fallback ? (
          <div className="absolute inset-0 flex items-center justify-center bg-bg/92 p-6">
            <div className="max-w-lg rounded-lg border border-lineStrong bg-panel p-5">
              <div className="font-mono text-[11px] uppercase tracking-[0.22em] text-rust">Could not load structure. Try another PDB.</div>
              <h3 className="mt-3 text-2xl font-semibold text-ivory">{structure.pdb_id}</h3>
              <p className="mt-2 text-sm leading-6 text-muted">
                {structure.label}. Protein {structure.protein}, subtype context {structure.subtype_context}. The app keeps the structure panel contained and falls back safely if 3Dmol or RCSB loading fails.
              </p>
              <a className="mt-4 inline-block rounded-md border border-lineStrong px-3 py-2 text-sm text-ivory hover:border-teal" href={structure.rcsb_url} target="_blank" rel="noreferrer">
                Open RCSB entry
              </a>
            </div>
          </div>
        ) : null}
      </div>

      <div className="border-t border-line px-4 py-3 text-xs leading-5 text-muted">
        <span className="font-mono uppercase tracking-[0.18em] text-brass">{viewerState}</span>
        <span className="ml-3">{message}</span>
      </div>
    </div>
  );
}
