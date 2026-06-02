import { describe, expect, it } from "vitest";
import { calculate, gaps, neighborScan, parseCsv, spiralFit } from "../src/utils/analysis";

const csv = `name,family,mass_mev,lifetime_s,stability,dominant_interaction
muon,lepton,105.658,2.19e-6,unstable,weak
tau,lepton,1776.86,2.9e-13,unstable,weak
electron,lepton,0.511,,stable,stable`;
describe("analysis", () => {
  it("parses required CSV fields", () => expect(parseCsv(csv)).toHaveLength(3));
  it("excludes stable rows when requested", () => expect(calculate(parseCsv(csv), 10, "exclude", 1e35)).toHaveLength(2));
  it("includes stable rows with manual tau", () => expect(calculate(parseCsv(csv), Math.PI, "manual", 1e35)).toHaveLength(3));
  it("computes circular gaps and neighbors", () => { const rows = calculate(parseCsv(csv), 2, "manual", 1e35); expect(gaps(rows)).toHaveLength(3); expect(neighborScan(rows, 1.2, 5, 8)).toHaveLength(3); });
  it("computes the declared golden spiral slope", () => { const fit = spiralFit(calculate(parseCsv(csv), 2, "manual", 1e35), "sorted_log10_N"); expect(fit.goldenK).toBeCloseTo(2 * Math.log((1 + Math.sqrt(5)) / 2) / Math.PI); });
});
