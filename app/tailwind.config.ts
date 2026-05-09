import type { Config } from "tailwindcss";

const config: Config = {
  content: ["./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        bg: "#080807",
        ink: "#0b0b09",
        graphite: "#12110e",
        panel: "rgba(18, 17, 14, 0.82)",
        panel2: "rgba(24, 23, 19, 0.72)",
        ivory: "#e8e1d2",
        muted: "#a8a092",
        brass: "#87a99b",
        brassSoft: "#27473f",
        teal: "#87a99b",
        sage: "#9aab91",
        rust: "#b9785f",
        dangerSoft: "#8a4d45",
        line: "rgba(135, 169, 155, 0.16)",
        lineStrong: "rgba(135, 169, 155, 0.36)"
      },
      fontFamily: {
        sans: ["Inter", "ui-sans-serif", "system-ui", "sans-serif"],
        mono: ["IBM Plex Mono", "ui-monospace", "SFMono-Regular", "monospace"]
      }
    }
  },
  plugins: []
};

export default config;
