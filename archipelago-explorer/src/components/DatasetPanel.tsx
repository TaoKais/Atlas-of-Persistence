import type { Entity, StableHandling } from "../utils/analysis";

interface Props {
  entities: Entity[]; family: string; interaction: string; stableHandling: StableHandling; manualTau: number;
  onUpload: (file: File) => void; onFamily: (value: string) => void; onInteraction: (value: string) => void;
  onStableHandling: (value: StableHandling) => void; onManualTau: (value: number) => void;
}
export function DatasetPanel(props: Props) {
  const families = [...new Set(props.entities.map((entity) => entity.family))].sort();
  const interactions = [...new Set(props.entities.map((entity) => entity.dominant_interaction))].sort();
  return <section className="panel">
    <h2>Dataset</h2>
    <label>Upload CSV<input type="file" accept=".csv,text/csv" onChange={(event) => event.target.files?.[0] && props.onUpload(event.target.files[0])} /></label>
    <div className="inline">
      <label>Family<select value={props.family} onChange={(event) => props.onFamily(event.target.value)}><option value="">all</option>{families.map((item) => <option key={item}>{item}</option>)}</select></label>
      <label>Interaction<select value={props.interaction} onChange={(event) => props.onInteraction(event.target.value)}><option value="">all</option>{interactions.map((item) => <option key={item}>{item}</option>)}</select></label>
    </div>
    <label>Stable-particle handling<select value={props.stableHandling} onChange={(event) => props.onStableHandling(event.target.value as StableHandling)}>
      <option value="exclude">exclude</option><option value="lower_bound">lower bound or manual fallback</option>
      <option value="manual">manual tau</option><option value="infinity_marker">infinity marker (excluded from finite plots)</option>
    </select></label>
    <label>Manual stable tau (s)<input type="number" value={props.manualTau} onChange={(event) => props.onManualTau(Number(event.target.value))} /></label>
    <details><summary>Preview source rows ({props.entities.length})</summary>
      <div className="table-wrap"><table><thead><tr><th>name</th><th>family</th><th>mass MeV</th><th>lifetime s</th><th>stability</th></tr></thead>
      <tbody>{props.entities.slice(0, 14).map((row) => <tr key={row.name}><td>{row.name}</td><td>{row.family}</td><td>{row.mass_mev}</td><td>{row.lifetime_s ?? "missing"}</td><td>{row.stability}</td></tr>)}</tbody></table></div>
    </details>
  </section>;
}
