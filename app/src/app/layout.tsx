import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "FluGenome3D Visual Lab",
  description: "Visual lab for real derived Influenza A HA/NA sequence, tokenization and structure artifacts."
};

export default function RootLayout({ children }: Readonly<{ children: React.ReactNode }>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
