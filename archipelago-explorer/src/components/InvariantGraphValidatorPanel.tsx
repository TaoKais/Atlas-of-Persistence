import { useEffect, useMemo, useState } from "react";
import Papa from "papaparse";

interface EdgeRow {
  particle_a: string;
  particle_b: string;
  tolerance_value: string;
  persistence_fraction: string;
}

interface SummaryRow {
  graph_type: string;
  tolerance_value: string;
  control_type: string;
  empirical_p_value: string;
  interpretation: string;
}

const graphTypes = ["angular_closeness", "3d_distance", "knn", "multi_representation", "multi_base_persistent"];
const representations = ["unit_circle", "cylindrical_helicoid", "radial_helicoid", "physical_plane", "physical_3d"];
const stableModes = ["finite_lifetime_only", "stable_excluded", "stable_truncated_1e25", "stable_truncated_1e30", "stable_truncated_1e35", "stable_as_infinity_marker"];
const controls = ["shuffle_lifetimes", "shuffle_masses", "shuffle_N_values", "random_logN_same_range", "random_phase_same_z", "family_label_shuffle", "interaction_label_shuffle"];

function parse<T>(text: string): T[] {
  return Papa.parse<T>(text, { header: true, skipEmptyLines: true }).data;
}

export function InvariantGraphValidatorPanel() {
  const [edges, setEdges] = useState<EdgeRow[]>([]);
  const [summary, setSummary] = useState<SummaryRow[]>([]);
  const [graphType, setGraphType] = useState("angular_closeness");
  const [tolerance, setTolerance] = useState("10");
  const [baseMin, setBaseMin] = useState(1.2);
  const [baseMax, setBaseMax] = useState(50);
  const [representation, setRepresentation] = useState("physical_3d");
  const [stableMode, setStableMode] = useState("finite_lifetime_only");
  const [control, setControl] = useState("shuffle_N_values");
  const [runs, setRuns] = useState(1000);

  useEffect(() => {
    fetch("/data/invariant_edges.csv").then((response) => response.ok ? response.text() : "").then((text) => text && setEdges(parse<EdgeRow>(text))).catch(() => setEdges([]));
    fetch("/data/graph_validation_summary.csv").then((response) => response.ok ? response.text() : "").then((text) => text && setSummary(parse<SummaryRow>(text))).catch(() => setSummary([]));
  }, []);

  const visibleEdges = useMemo(() => edges
    .filter((row) => row.tolerance_value === tolerance)
    .sort((a, b) => Number(b.persistence_fraction) - Number(a.persistence_fraction))
    .slice(0, 8), [edges, tolerance]);

  const selectedSummary = useMemo(() => summary.find((row) =>
    row.graph_type === graphType && row.tolerance_value === tolerance && row.control_type === control
  ) ?? summary.find((row) => row.graph_type === graphType && row.control_type === control), [summary, graphType, tolerance, control]);

  const pValue = selectedSummary ? Number(selectedSummary.empirical_p_value) : NaN;

  return <section className="panel invariant-validator">
    <h2>Invariant Graph Validator</h2>
    <div className="control-grid">
      <label>Graph type<select value={graphType} onChange={(event) => setGraphType(event.target.value)}>{graphTypes.map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
      <label>Tolerance<select value={tolerance} onChange={(event) => setTolerance(event.target.value)}>{["5", "10", "15", "20", "30"].map((value) => <option key={value} value={value}>{value} deg</option>)}</select></label>
      <label>Base min<input type="number" value={baseMin} min={1.2} max={baseMax} step={0.1} onChange={(event) => setBaseMin(Number(event.target.value))} /></label>
      <label>Base max<input type="number" value={baseMax} min={baseMin} max={50} step={0.1} onChange={(event) => setBaseMax(Number(event.target.value))} /></label>
      <label>Representation<select value={representation} onChange={(event) => setRepresentation(event.target.value)}>{representations.map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
      <label>Stable handling<select value={stableMode} onChange={(event) => setStableMode(event.target.value)}>{stableModes.map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
      <label>Random control<select value={control} onChange={(event) => setControl(event.target.value)}>{controls.map((value) => <option key={value} value={value}>{value}</option>)}</select></label>
      <label>Monte Carlo<input type="number" value={runs} min={1} max={10000} step={100} onChange={(event) => setRuns(Number(event.target.value))} /></label>
    </div>
    <div className="metrics"><div><strong>{Number.isFinite(pValue) ? pValue.toFixed(4) : "n/a"}</strong><span>empirical p-value</span></div><div><strong>{selectedSummary?.interpretation ?? "Insufficient data"}</strong><span>robustness label</span></div><div><strong>{tolerance} deg</strong><span>tolerance used</span></div><div><strong>{runs}</strong><span>Monte Carlo runs</span></div></div>
    <div className="mini-graph">{visibleEdges.map((row) => <div className="edge-row" key={`${row.particle_a}-${row.particle_b}-${row.tolerance_value}`}><span>{row.particle_a} / {row.particle_b}</span><meter min="0" max="1" value={Number(row.persistence_fraction)} /><strong>{Number(row.persistence_fraction).toFixed(2)}</strong></div>)}</div>
    <div className="button-grid"><a className="button-link" href="/data/invariant_edges.csv" download>Edges CSV</a><a className="button-link" href="/data/graph_validation_summary.csv" download>Summary CSV</a><a className="button-link" href="/data/robust_clusters.csv" download>Clusters CSV</a></div>
    <div className="warning">Persistent graph structure does not imply physical interaction. Every displayed result is tied to the selected tolerance and randomized-control p-value.</div>
  </section>;
}
