import type { Config } from "tailwindcss";

const config: Config = {
  content: ["./src/**/*.{ts,tsx}"],
  theme: {
    extend: {
      colors: {
        ink: "#080806",
        panel: "#11100d",
        panel2: "#171511",
        ivory: "#eee4cf",
        muted: "#9b9383",
        amber: "#d7a84a",
        mint: "#7fa69a",
        rust: "#b87554",
        line: "rgba(238, 228, 207, 0.14)"
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
