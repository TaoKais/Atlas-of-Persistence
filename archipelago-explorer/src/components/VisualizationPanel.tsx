import Plot from "react-plotly.js";
import Plotly from "plotly.js-dist-min";
import type { CalculatedEntity, Gap, Neighbor, SpiralFit } from "../utils/analysis";

export type View = "frequency" | "persistence" | "circle" | "cylindrical" | "radial" | "gaps" | "neighbors" | "family_centroids" | "interaction_centroids";
const views: [View, string][] = [["frequency", "Frequency vs lifetime"], ["persistence", "Persistence index"], ["circle", "Unit circle phase"], ["cylindrical", "Cylindrical helicoid 3D"], ["radial", "Radial helicoid 3D"], ["gaps", "Gap centers"], ["neighbors", "Persistent neighbor network"], ["family_centroids", "Family centroid evolution"], ["interaction_centroids", "Interaction centroid evolution"]];
const familyColors: Record<string, string> = { lepton: "#4cc9f0", meson: "#f72585", baryon: "#fca311", boson: "#80ed99", quark: "#b5179e" };
const color = (row: CalculatedEntity) => familyColors[row.family] ?? "#adb5bd";
const config = { responsive: true, displaylogo: false };

function centroidTraces(rows: CalculatedEntity[], category: "family" | "dominant_interaction", minBase: number, maxBase: number) {
  const bases = Array.from({ length: 40 }, (_, index) => minBase + index * (maxBase - minBase) / 39);
  const groups = [...new Set(rows.map((row) => row[category]))];
  return groups.map((group) => ({ type: "scatter", mode: "lines", name: group, x: bases, y: bases.map((base) => {
    const selected = rows.filter((row) => row[category] === group);
    const vectors = selected.map((row) => { const theta = 2 * Math.PI * ((Math.log(row.N) / Math.log(base)) % 1 + 1) % 1; return [Math.cos(theta), Math.sin(theta)]; });
    return Math.atan2(vectors.reduce((s, v) => s + v[1], 0), vectors.reduce((s, v) => s + v[0], 0)) * 180 / Math.PI;
  }) }));
}

export function VisualizationPanel({ rows, gaps, neighbors, view, onView, minBase, maxBase, spiral, showSpiralOverlay }: { rows: CalculatedEntity[]; gaps: Gap[]; neighbors: Neighbor[]; view: View; onView: (view: View) => void; minBase: number; maxBase: number; spiral: SpiralFit; showSpiralOverlay: boolean }) {
  let data: any[] = [], layout: any = { paper_bgcolor: "#101722", plot_bgcolor: "#101722", font: { color: "#dbe8f4" }, margin: { t: 48, r: 20, b: 52, l: 62 } };
  if (view === "frequency") data = [{ type: "scatter", mode: "markers+text", x: rows.map((r) => Math.log10(r.compton_frequency_hz)), y: rows.map((r) => Math.log10(r.effective_lifetime_s)), text: rows.map((r) => r.name), textposition: "top center", marker: { color: rows.map(color), size: 9 } }];
  if (view === "persistence") data = [{ type: "bar", x: rows.map((r) => r.name), y: rows.map((r) => r.log10_N), marker: { color: rows.map(color) } }];
  if (view === "circle") data = [{ type: "scatter", mode: "markers+text", x: rows.map((r) => r.x), y: rows.map((r) => r.y), text: rows.map((r) => r.name), textposition: "top center", marker: { color: rows.map(color), size: 10 } }], layout.yaxis = { scaleanchor: "x" };
  if (view === "cylindrical") data = [{ type: "scatter3d", mode: "markers+text", x: rows.map((r) => r.x), y: rows.map((r) => r.y), z: rows.map((r) => r.log10_N), text: rows.map((r) => r.name), marker: { color: rows.map(color), size: 4 } }];
  if (view === "radial") data = [{ type: "scatter3d", mode: "markers+text", x: rows.map((r) => r.radial_x), y: rows.map((r) => r.radial_y), z: rows.map((r) => r.log10_N), text: rows.map((r) => r.name), marker: { color: rows.map(color), size: 4 } }];
  if (view === "radial" && showSpiralOverlay) data.push({ type: "scatter3d", mode: "lines", name: "exploratory golden spiral reference", x: spiral.overlay.x, y: spiral.overlay.y, z: spiral.overlay.z, line: { color: "#fca311", width: 5 } });
  if (view === "gaps" && gaps.length) { const centers = gaps.slice(0, 4).map((gap) => gap.center_deg * Math.PI / 180), x = centers.map(Math.cos), y = centers.map(Math.sin); data = [
    { type: "scatter", mode: "lines", x: [...x, x[0]], y: [...y, y[0]], line: { color: "#fca311", width: 3 }, name: "top-4 polygon" },
    { type: "scatter", mode: "markers+text", x, y, text: gaps.slice(0, 4).map((gap) => `${gap.from} → ${gap.to}`), textposition: "top center", marker: { color: "#f72585", size: 11 }, name: "gap centers" }
  ]; layout.yaxis = { scaleanchor: "x" }; }
  if (view === "neighbors") { const top = neighbors.slice(0, 14), names = [...new Set(top.flatMap((n) => [n.particle_a, n.particle_b]))], angles = names.map((_, i) => 2 * Math.PI * i / names.length), positions = Object.fromEntries(names.map((name, i) => [name, [Math.cos(angles[i]), Math.sin(angles[i])]])); data = [
    ...top.map((n) => ({ type: "scatter", mode: "lines", x: [positions[n.particle_a][0], positions[n.particle_b][0]], y: [positions[n.particle_a][1], positions[n.particle_b][1]], line: { width: 1 + 5 * n.fraction, color: "#4cc9f0" }, opacity: .25 + .75 * n.fraction, hoverinfo: "text", text: `${n.particle_a} - ${n.particle_b}: ${n.fraction.toFixed(3)}` })),
    { type: "scatter", mode: "markers+text", x: names.map((name) => positions[name][0]), y: names.map((name) => positions[name][1]), text: names, textposition: "top center", marker: { size: 11, color: "#f72585" } }]; layout.yaxis = { scaleanchor: "x" }; }
  if (view === "family_centroids") data = centroidTraces(rows, "family", minBase, maxBase);
  if (view === "interaction_centroids") data = centroidTraces(rows, "dominant_interaction", minBase, maxBase);
  layout.title = views.find(([id]) => id === view)?.[1];
  const download = (format: "png" | "svg") => Plotly.downloadImage("archipelago-plot", { format, filename: `archipelago-${view}` });
  return <section className="panel plot-panel"><div className="plot-toolbar"><h2>Visualization</h2><select value={view} onChange={(event) => onView(event.target.value as View)}>{views.map(([id, name]) => <option key={id} value={id}>{name}</option>)}</select><button onClick={() => download("png")}>PNG</button><button onClick={() => download("svg")}>SVG</button></div>
    <Plot divId="archipelago-plot" data={data} layout={{ ...layout, autosize: true }} config={config} useResizeHandler style={{ width: "100%", height: "620px" }} />
  </section>;
}
