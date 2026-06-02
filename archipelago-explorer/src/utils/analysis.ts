import Papa from "papaparse";
import { constants, MEV_TO_J } from "../physics/constants";

export type StableHandling = "exclude" | "lower_bound" | "manual" | "infinity_marker";
export type ControlMode = "observed" | "shuffle_lifetimes" | "shuffle_masses" | "randomize_n";
export interface Entity {
  name: string; family: string; mass_mev: number; lifetime_s?: number; stability: string;
  dominant_interaction: string; spin?: string; charge?: string; decay_mode_count?: number;
  source?: string; notes?: string;
}
export interface CalculatedEntity extends Entity {
  effective_lifetime_s: number; compton_frequency_hz: number; N: number; log10_N: number;
  theta: number; theta_deg: number; x: number; y: number; radial_x: number; radial_y: number;
}
export interface Gap { from: string; to: string; gap_deg: number; center_deg: number; }
export interface Neighbor { particle_a: string; particle_b: string; fraction: number; mean_distance: number; }

export const baseOptions = [
  ["2", 2], ["e", Math.E], ["pi", Math.PI], ["phi", constants.phi.value], ["10", 10],
  ["2pi", 2 * Math.PI], ["sqrt(2pi)", Math.sqrt(2 * Math.PI)], ["pi^2", Math.PI ** 2],
  ["e^pi", Math.E ** Math.PI],
] as const;

const positive = (value: unknown) => Number.isFinite(Number(value)) && Number(value) > 0;
const frac = (value: number) => value - Math.floor(value);
export const angularDistance = (left: number, right: number) =>
  Math.abs((((left - right) * 180 / Math.PI + 180) % 360 + 360) % 360 - 180);

export function parseCsv(text: string): Entity[] {
  const parsed = Papa.parse<Record<string, string>>(text, { header: true, skipEmptyLines: true });
  if (parsed.errors.length) throw new Error(parsed.errors[0].message);
  const required = ["name", "family", "mass_mev", "lifetime_s", "stability", "dominant_interaction"];
  const missing = required.filter((key) => !parsed.meta.fields?.includes(key));
  if (missing.length) throw new Error(`Missing required fields: ${missing.join(", ")}`);
  return parsed.data.map((row) => ({
    name: row.name, family: row.family, mass_mev: Number(row.mass_mev),
    lifetime_s: positive(row.lifetime_s) ? Number(row.lifetime_s) : undefined,
    stability: row.stability, dominant_interaction: row.dominant_interaction,
    spin: row.spin, charge: row.charge,
    decay_mode_count: Number(row.decay_mode_count || row.representative_decay_mode_count || 0),
    source: row.source, notes: row.notes,
  })).filter((row) => row.name && positive(row.mass_mev));
}

export function calculate(entities: Entity[], base: number, handling: StableHandling, manualTau: number): CalculatedEntity[] {
  return entities.flatMap((entity) => {
    const stable = entity.stability.toLowerCase() === "stable";
    let tau = entity.lifetime_s;
    if (stable && handling === "exclude") return [];
    if (stable && handling === "manual") tau = manualTau;
    if (stable && handling === "lower_bound") tau = entity.lifetime_s ?? manualTau;
    if (stable && handling === "infinity_marker") return [];
    if (!positive(tau)) return [];
    const compton = entity.mass_mev * MEV_TO_J / constants.h.value;
    const N = compton * Number(tau);
    const log10_N = Math.log10(N);
    const theta = 2 * Math.PI * frac(Math.log(N) / Math.log(base));
    const radius = 0; // replaced after normalization
    return [{ ...entity, effective_lifetime_s: Number(tau), compton_frequency_hz: compton, N, log10_N,
      theta, theta_deg: theta * 180 / Math.PI, x: Math.cos(theta), y: Math.sin(theta),
      radial_x: radius, radial_y: radius }];
  }).map((entity, _, rows) => {
    const logs = rows.map((row) => row.log10_N);
    const radius = (entity.log10_N - Math.min(...logs)) / (Math.max(...logs) - Math.min(...logs) || 1);
    return { ...entity, radial_x: radius * Math.cos(entity.theta), radial_y: radius * Math.sin(entity.theta) };
  });
}

export function applyControl(rows: CalculatedEntity[], mode: ControlMode): CalculatedEntity[] {
  if (mode === "observed") return rows;
  const reversed = [...rows].reverse();
  return rows.map((row, index) => {
    const donor = reversed[index];
    let compton = row.compton_frequency_hz, tau = row.effective_lifetime_s, N = row.N;
    if (mode === "shuffle_lifetimes") tau = donor.effective_lifetime_s;
    if (mode === "shuffle_masses") compton = donor.compton_frequency_hz;
    if (mode === "randomize_n") N = 10 ** (Math.min(...rows.map((x) => x.log10_N)) + (index * 0.61803398875 % 1) * (Math.max(...rows.map((x) => x.log10_N)) - Math.min(...rows.map((x) => x.log10_N))));
    else N = compton * tau;
    return { ...row, compton_frequency_hz: compton, effective_lifetime_s: tau, N, log10_N: Math.log10(N) };
  });
}

export function recalculatePhase(rows: CalculatedEntity[], base: number): CalculatedEntity[] {
  const logs = rows.map((row) => row.log10_N), min = Math.min(...logs), span = Math.max(...logs) - min || 1;
  return rows.map((row) => {
    const theta = 2 * Math.PI * frac(Math.log(row.N) / Math.log(base));
    const radius = (row.log10_N - min) / span;
    return { ...row, theta, theta_deg: theta * 180 / Math.PI, x: Math.cos(theta), y: Math.sin(theta),
      radial_x: radius * Math.cos(theta), radial_y: radius * Math.sin(theta) };
  });
}

export function gaps(rows: CalculatedEntity[]): Gap[] {
  const sorted = [...rows].sort((a, b) => a.theta - b.theta);
  return sorted.map((row, index) => {
    const next = sorted[(index + 1) % sorted.length];
    const delta = ((next.theta - row.theta) + 2 * Math.PI) % (2 * Math.PI);
    return { from: row.name, to: next.name, gap_deg: delta * 180 / Math.PI,
      center_deg: ((row.theta + delta / 2) % (2 * Math.PI)) * 180 / Math.PI };
  }).sort((a, b) => b.gap_deg - a.gap_deg);
}

export function metrics(rows: CalculatedEntity[]) {
  if (!rows.length) return { R: 0, entropy: 0, symmetry: 0, largestGap: 0, clusterCount: 0 };
  const meanX = rows.reduce((sum, row) => sum + row.x, 0) / rows.length;
  const meanY = rows.reduce((sum, row) => sum + row.y, 0) / rows.length;
  const R = Math.hypot(meanX, meanY);
  const bins = Array(12).fill(0);
  rows.forEach((row) => bins[Math.min(11, Math.floor(row.theta / (2 * Math.PI) * 12))]++);
  const entropy = -bins.filter(Boolean).reduce((sum, count) => { const p = count / rows.length; return sum + p * Math.log2(p); }, 0);
  const top = gaps(rows);
  const clusterCount = Math.max(1, top.filter((gap) => gap.gap_deg > 30).length);
  const angles = top.slice(0, 4).map((gap) => gap.center_deg).sort((a, b) => a - b);
  const spacing = angles.map((angle, index) => ((angles[(index + 1) % angles.length] - angle) + 360) % 360);
  const symmetry = 1 / (1 + spacing.reduce((sum, value) => sum + Math.abs(value - 90), 0) / spacing.length);
  return { R, entropy, symmetry, largestGap: top[0]?.gap_deg ?? 0, clusterCount };
}

export function neighborScan(rows: CalculatedEntity[], minBase: number, maxBase: number, samples = 80): Neighbor[] {
  const bases = Array.from({ length: samples }, (_, index) => minBase + index * (maxBase - minBase) / (samples - 1));
  return Array.from({ length: rows.length }, (_, left) => left).flatMap((left) =>
    Array.from({ length: rows.length - left - 1 }, (_, offset) => left + offset + 1).map((right) => {
      const distances = bases.map((base) => {
        const phased = recalculatePhase([rows[left], rows[right]], base);
        return angularDistance(phased[0].theta, phased[1].theta);
      });
      return { particle_a: rows[left].name, particle_b: rows[right].name,
        fraction: distances.filter((distance) => distance < 20).length / distances.length,
        mean_distance: distances.reduce((sum, distance) => sum + distance, 0) / distances.length };
    })
  ).sort((a, b) => b.fraction - a.fraction || a.mean_distance - b.mean_distance);
}
