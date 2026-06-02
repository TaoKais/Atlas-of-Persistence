import type { ControlMode, SpiralFit, SpiralFitOrdering } from "../utils/analysis";

export function GoldenRatioPanel({ fit, ordering, overlay, onOrdering, onOverlay, onNeighborhood, onControl }: {
  fit: SpiralFit; ordering: SpiralFitOrdering; overlay: boolean;
  onOrdering: (ordering: SpiralFitOrdering) => void; onOverlay: (value: boolean) => void;
  onNeighborhood: () => void; onControl: (control: ControlMode) => void;
}) {
  return <section className="panel">
    <h2>Golden Ratio / Spiral Test</h2>
    <label>Spiral fit ordering<select value={ordering} onChange={(event) => onOrdering(event.target.value as SpiralFitOrdering)}>
      <option value="sorted_theta">sorted by theta</option><option value="sorted_log10_N">sorted by log10(N)</option>
      <option value="persistent_neighbor_path">persistent-neighbor path</option><option value="nearest_neighbor_3d">nearest-neighbor path in 3D</option>
    </select></label>
    <label><input type="checkbox" checked={overlay} onChange={(event) => onOverlay(event.target.checked)} /> show golden spiral overlay</label>
    <div className="button-grid"><button onClick={onNeighborhood}>Phi-neighborhood scan</button><button onClick={() => onControl("randomize_n")}>Compare randomized control</button></div>
    <div className="metrics"><div><strong>{fit.fittedK.toFixed(5)}</strong><span>fitted k</span></div><div><strong>{fit.goldenK.toFixed(5)}</strong><span>golden k</span></div><div><strong>{fit.relativeError.toFixed(4)}</strong><span>relative error</span></div><div><strong>{fit.R2.toFixed(4)}</strong><span>R^2</span></div></div>
    <div className="warning">Golden spiral fit is exploratory. A visual spiral does not imply physical significance.</div>
  </section>;
}
