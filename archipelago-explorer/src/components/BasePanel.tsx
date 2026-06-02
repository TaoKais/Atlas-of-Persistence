import { baseOptions } from "../utils/analysis";

interface Props { baseName: string; base: number; minBase: number; maxBase: number; onBaseName: (name: string) => void; onBase: (base: number) => void; onMin: (base: number) => void; onMax: (base: number) => void; }
export function BasePanel(props: Props) {
  return <section className="panel">
    <h2>Logarithmic Base</h2>
    <label>Preset<select value={props.baseName} onChange={(event) => props.onBaseName(event.target.value)}>
      {baseOptions.map(([name]) => <option key={name}>{name}</option>)}<option value="custom">custom</option>
    </select></label>
    <label>Current base: {props.base.toFixed(4)}<input type="range" min="1.2" max="50" step="0.01" value={props.base} onChange={(event) => props.onBase(Number(event.target.value))} /></label>
    <div className="inline"><label>scan min<input type="number" min="1.01" step="0.1" value={props.minBase} onChange={(event) => props.onMin(Number(event.target.value))} /></label>
    <label>scan max<input type="number" step="0.1" value={props.maxBase} onChange={(event) => props.onMax(Number(event.target.value))} /></label></div>
  </section>;
}
