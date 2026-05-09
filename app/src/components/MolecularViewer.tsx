"use client";

import { useEffect, useRef, useState } from "react";

type StructureRecord = {
  pdb_id: string;
  label: string;
  pdb_download_url: string;
  mapping_status: string;
};

export default function MolecularViewer({ structure }: { structure: StructureRecord }) {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [status, setStatus] = useState("initializing viewer");

  useEffect(() => {
    let cancelled = false;
    async function loadStructure() {
      if (!containerRef.current) return;
      containerRef.current.innerHTML = "";
      setStatus(`loading ${structure.pdb_id} from RCSB`);
      try {
        const mod = await import("3dmol");
        const threeDmol = (mod as any).default ?? mod;
        const viewer = threeDmol.createViewer(containerRef.current, {
          backgroundColor: "#080806"
        });
        const response = await fetch(structure.pdb_download_url);
        if (!response.ok) throw new Error(`RCSB returned ${response.status}`);
        const pdb = await response.text();
        if (cancelled) return;
        viewer.addModel(pdb, "pdb");
        viewer.setStyle({}, { cartoon: { color: "spectrum", opacity: 0.92 } });
        viewer.addSurface(threeDmol.SurfaceType.VDW, {
          opacity: 0.12,
          color: "#d7a84a"
        });
        viewer.zoomTo();
        viewer.spin(true);
        viewer.render();
        setStatus(`${structure.pdb_id} loaded. Mapping status: ${structure.mapping_status}`);
      } catch (error) {
        setStatus(`viewer unavailable: ${error instanceof Error ? error.message : "unknown error"}`);
      }
    }

    loadStructure();
    return () => {
      cancelled = true;
    };
  }, [structure]);

  return (
    <div className="overflow-hidden rounded-lg border border-line bg-ink">
      <div ref={containerRef} className="h-[480px] w-full" />
      <div className="border-t border-line px-4 py-3 font-mono text-[11px] uppercase tracking-[0.2em] text-muted">
        {status}
      </div>
    </div>
  );
}
