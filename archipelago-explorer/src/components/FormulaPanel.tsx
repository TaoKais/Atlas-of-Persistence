import katex from "katex";
import "katex/dist/katex.min.css";
import { formulas } from "../formulas/catalog";

export function FormulaPanel({ formulaId, onChange }: { formulaId: string; onChange: (id: string) => void }) {
  const formula = formulas.find((item) => item.id === formulaId) ?? formulas[0];
  return <section className="panel">
    <h2>Formula Framework</h2>
    <select value={formula.id} onChange={(event) => onChange(event.target.value)}>
      {formulas.map((item) => <option key={item.id} value={item.id}>{item.name}</option>)}
    </select>
    <span className={`tag ${formula.kind}`}>{formula.kind}</span>
    <div className="formula" dangerouslySetInnerHTML={{ __html: katex.renderToString(formula.latex) }} />
    <p>{formula.units}</p>
  </section>;
}
