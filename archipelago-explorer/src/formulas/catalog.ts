export type FormulaKind = "standard" | "exploratory";
export interface FormulaDefinition {
  id: string;
  name: string;
  latex: string;
  units: string;
  kind: FormulaKind;
}

export const formulas: FormulaDefinition[] = [
  { id: "planck", name: "Planck relation", latex: "E = h f", units: "E in J, f in Hz", kind: "standard" },
  { id: "compton", name: "Compton frequency", latex: "f_C = \\frac{m c^2}{h}", units: "Hz", kind: "standard" },
  { id: "omega", name: "Reduced Compton angular frequency", latex: "\\omega_C = \\frac{m c^2}{\\hbar}", units: "rad s^{-1}", kind: "standard" },
  { id: "width", name: "Lifetime-width relation", latex: "\\tau \\approx \\frac{\\hbar}{\\Gamma}", units: "s", kind: "standard" },
  { id: "quality", name: "Quality factor", latex: "Q = \\frac{f}{\\Delta f}", units: "dimensionless", kind: "standard" },
  { id: "heisenberg", name: "Heisenberg energy-time relation", latex: "\\Delta E\\,\\Delta t \\geq \\frac{\\hbar}{2}", units: "J s", kind: "standard" },
  { id: "debroglie", name: "de Broglie wavelength", latex: "\\lambda = \\frac{h}{p}", units: "m", kind: "standard" },
  { id: "compactness", name: "Schwarzschild compactness", latex: "\\Phi = \\frac{G M}{R c^2}", units: "dimensionless", kind: "standard" },
  { id: "radius", name: "Schwarzschild radius", latex: "R_s = \\frac{2 G M}{c^2}", units: "m", kind: "standard" },
  { id: "redshift", name: "Gravitational redshift approximation", latex: "f_\\infty = f_{local}\\sqrt{1 - 2\\Phi}", units: "Hz", kind: "standard" },
  { id: "cycles", name: "Persistence cycles", latex: "N = f_C\\tau", units: "dimensionless cycle count", kind: "exploratory" },
  { id: "persistence", name: "Persistence logarithm", latex: "P = \\log_{10}(N)", units: "dimensionless", kind: "exploratory" },
  { id: "phase", name: "Phase mapping", latex: "\\theta = 2\\pi\\,\\mathrm{frac}(\\log_b(N))", units: "rad", kind: "exploratory" },
  { id: "complex", name: "Complex phase", latex: "z = \\exp(i\\theta)", units: "dimensionless", kind: "exploratory" },
  { id: "cylindrical", name: "Cylindrical helicoid", latex: "x=\\cos\\theta,\\ y=\\sin\\theta,\\ z=\\log_{10}(N)", units: "exploratory coordinates", kind: "exploratory" },
  { id: "radial", name: "Radial helicoid", latex: "r=\\widehat{\\log_{10}(N)},\\ x=r\\cos\\theta,\\ y=r\\sin\\theta", units: "exploratory coordinates", kind: "exploratory" },
  { id: "gaps", name: "Gap-center geometry", latex: "\\Delta\\theta_i = \\theta_{i+1} - \\theta_i", units: "rad", kind: "exploratory" },
  { id: "neighbors", name: "Neighbor persistence", latex: "S_{ij}=\\frac{\\#\\{b:d_{ij}(b)<\\epsilon\\}}{\\#\\{b\\}}", units: "fraction", kind: "exploratory" },
];
