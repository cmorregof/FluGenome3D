import type { Config } from "tailwindcss";

const config: Config = {
  content: ["./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        bg: "#01050a",
        ink: "#030910",
        graphite: "#06111a",
        panel: "rgba(5, 17, 24, 0.76)",
        panel2: "rgba(7, 26, 36, 0.68)",
        ivory: "#edf7f4",
        muted: "#9fb3ae",
        brass: "#5cdce2",
        brassSoft: "#0d3c43",
        teal: "#5cdce2",
        sage: "#79d99c",
        rust: "#b9785f",
        dangerSoft: "#8a4d45",
        line: "rgba(92, 218, 226, 0.16)",
        lineStrong: "rgba(92, 218, 226, 0.36)"
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
