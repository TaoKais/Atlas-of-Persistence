import type { Gap, Neighbor } from "../utils/analysis";

export function AnalysisPanel({ stats, gaps, neighbors }: { stats: { R: number; entropy: number; symmetry: number; largestGap: number; clusterCount: number }; gaps: Gap[]; neighbors: Neighbor[] }) {
  return <section className="panel">
    <h2>Analysis</h2>
    <div className="metrics"><div><strong>{stats.R.toFixed(4)}</strong><span>resultant R</span></div><div><strong>{stats.entropy.toFixed(4)}</strong><span>angular entropy</span></div><div><strong>{stats.symmetry.toFixed(4)}</strong><span>top-4 symmetry</span></div><div><strong>{stats.largestGap.toFixed(2)}°</strong><span>largest gap</span></div><div><strong>{stats.clusterCount}</strong><span>phase clusters</span></div></div>
    <h3>Largest gaps</h3><ol>{gaps.slice(0, 5).map((gap) => <li key={`${gap.from}-${gap.to}`}>{gap.from} → {gap.to}: {gap.gap_deg.toFixed(2)}°</li>)}</ol>
    <h3>Persistent neighbors</h3><ol>{neighbors.slice(0, 6).map((row) => <li key={`${row.particle_a}-${row.particle_b}`}>{row.particle_a} / {row.particle_b}: {(100 * row.fraction).toFixed(1)}%</li>)}</ol>
    <div className="warning">If a structure also appears in randomized controls, it is likely a mathematical artifact. Visual structure is not evidence of physical significance.</div>
  </section>;
}
