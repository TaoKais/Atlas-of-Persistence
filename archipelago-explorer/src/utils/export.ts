import Papa from "papaparse";
import type { CalculatedEntity, ControlMode, Neighbor, StableHandling } from "./analysis";

function download(name: string, content: string, type: string) {
  const link = document.createElement("a");
  link.href = URL.createObjectURL(new Blob([content], { type }));
  link.download = name; link.click(); URL.revokeObjectURL(link.href);
}
export const exportCsv = (rows: CalculatedEntity[]) => download("archipelago-calculated.csv", Papa.unparse(rows), "text/csv");
export const exportJson = (state: unknown) => download("archipelago-project.json", JSON.stringify(state, null, 2), "application/json");
export function exportReport(rows: CalculatedEntity[], base: number, stable: StableHandling, control: ControlMode, neighbors: Neighbor[]) {
  const report = `# Archipelago Explorer Report

## Scope
Exploratory visualization only. This report does not claim new physics, hidden symmetries, or physical validation.

## Configuration
- Entities plotted: ${rows.length}
- Logarithmic base: ${base}
- Stable handling: ${stable}
- Robustness control: ${control}

## Top persistent neighbors
${neighbors.slice(0, 10).map((row) => `- ${row.particle_a} / ${row.particle_b}: ${(100 * row.fraction).toFixed(1)}%`).join("\n")}

## Warning
If a structure also appears in randomized controls, it is likely a mathematical artifact.
`;
  download("archipelago-report.md", report, "text/markdown");
}
