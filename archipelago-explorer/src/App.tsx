import { useEffect, useMemo, useState } from "react";
import { DatasetPanel } from "./components/DatasetPanel";
import { FormulaPanel } from "./components/FormulaPanel";
import { BasePanel } from "./components/BasePanel";
import { VisualizationPanel, type View } from "./components/VisualizationPanel";
import { AnalysisPanel } from "./components/AnalysisPanel";
import { applyControl, baseOptions, calculate, gaps, metrics, neighborScan, parseCsv, recalculatePhase, type ControlMode, type Entity, type StableHandling } from "./utils/analysis";
import { exportCsv, exportJson, exportReport } from "./utils/export";
import "./styles.css";

export default function App() {
  const [entities, setEntities] = useState<Entity[]>([]), [error, setError] = useState("");
  const [family, setFamily] = useState(""), [interaction, setInteraction] = useState("");
  const [stableHandling, setStableHandling] = useState<StableHandling>("exclude"), [manualTau, setManualTau] = useState(1e35);
  const [baseName, setBaseName] = useState("10"), [base, setBase] = useState(10), [minBase, setMinBase] = useState(1.2), [maxBase, setMaxBase] = useState(50);
  const [formula, setFormula] = useState("cylindrical"), [view, setView] = useState<View>("cylindrical"), [control, setControl] = useState<ControlMode>("observed");
  const [sweeping, setSweeping] = useState(false);
  useEffect(() => { fetch("/data/particles_reference.csv").then((response) => response.text()).then((text) => setEntities(parseCsv(text))).catch((reason) => setError(String(reason))); }, []);
  useEffect(() => { if (!sweeping) return; const id = setInterval(() => setBase((value) => value >= maxBase ? minBase : value + (maxBase - minBase) / 90), 110); return () => clearInterval(id); }, [sweeping, minBase, maxBase]);
  const selected = useMemo(() => entities.filter((row) => (!family || row.family === family) && (!interaction || row.dominant_interaction === interaction)), [entities, family, interaction]);
  const rows = useMemo(() => recalculatePhase(applyControl(calculate(selected, base, stableHandling, manualTau), control), base), [selected, base, stableHandling, manualTau, control]);
  const gapRows = useMemo(() => gaps(rows), [rows]), neighborRows = useMemo(() => neighborScan(rows, minBase, maxBase), [rows, minBase, maxBase]);
  const stats = useMemo(() => metrics(rows), [rows]);
  const chooseBase = (name: string) => { setBaseName(name); const preset = baseOptions.find(([label]) => label === name); if (preset) setBase(preset[1]); };
  const upload = (file: File) => file.text().then((text) => { setEntities(parseCsv(text)); setError(""); }).catch((reason) => setError(String(reason)));
  const state = { base, baseName, minBase, maxBase, stableHandling, manualTau, control, family, interaction, formula, view };
  return <><header><div><h1>Archipelago Explorer</h1><p>Frequency, Persistence and Regime Visualization Toolkit</p></div><span className="scope">Exploratory only · no claims of new physics</span></header>
    <main>{error && <div className="error">{error}</div>}<div className="sidebar">
      <DatasetPanel entities={entities} family={family} interaction={interaction} stableHandling={stableHandling} manualTau={manualTau} onUpload={upload} onFamily={setFamily} onInteraction={setInteraction} onStableHandling={setStableHandling} onManualTau={setManualTau} />
      <FormulaPanel formulaId={formula} onChange={setFormula} />
      <BasePanel baseName={baseName} base={base} minBase={minBase} maxBase={maxBase} onBaseName={chooseBase} onBase={(value) => { setBase(value); setBaseName("custom"); }} onMin={setMinBase} onMax={setMaxBase} />
      <section className="panel"><h2>Robustness Controls</h2><select value={control} onChange={(event) => setControl(event.target.value as ControlMode)}><option value="observed">observed dataset</option><option value="shuffle_lifetimes">shuffle lifetimes</option><option value="shuffle_masses">shuffle masses</option><option value="randomize_n">randomize N</option></select><button onClick={() => setSweeping(!sweeping)}>{sweeping ? "Stop" : "Animate"} base sweep</button></section>
      <section className="panel"><h2>Export</h2><div className="button-grid"><button onClick={() => exportCsv(rows)}>CSV</button><button onClick={() => exportJson(state)}>JSON state</button><button onClick={() => exportReport(rows, base, stableHandling, control, neighborRows)}>Markdown report</button></div></section>
    </div><div className="workspace"><VisualizationPanel rows={rows} gaps={gapRows} neighbors={neighborRows} view={view} onView={setView} minBase={minBase} maxBase={maxBase} /><AnalysisPanel stats={stats} gaps={gapRows} neighbors={neighborRows} /></div></main>
    <footer>Archipelago Explorer provides representations for inspection. Mathematical patterns require independent physical justification.</footer></>;
}
