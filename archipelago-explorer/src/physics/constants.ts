/** CODATA/SI constants used by the client-side exploratory calculator. */
export const constants = {
  c: { value: 299_792_458, unit: "m s^-1", label: "speed of light" },
  h: { value: 6.626_070_15e-34, unit: "J s", label: "Planck constant" },
  hbar: { value: 1.054_571_817e-34, unit: "J s", label: "reduced Planck constant" },
  G: { value: 6.674_30e-11, unit: "m^3 kg^-1 s^-2", label: "gravitational constant" },
  k_B: { value: 1.380_649e-23, unit: "J K^-1", label: "Boltzmann constant" },
  e: { value: 1.602_176_634e-19, unit: "C", label: "elementary charge" },
  alpha: { value: 7.297_352_5643e-3, unit: "dimensionless", label: "fine-structure constant" },
  pi: { value: Math.PI, unit: "dimensionless", label: "pi" },
  phi: { value: (1 + Math.sqrt(5)) / 2, unit: "dimensionless", label: "golden ratio" },
  planckMass: { value: 2.176_434e-8, unit: "kg", label: "Planck mass" },
  planckTime: { value: 5.391_247e-44, unit: "s", label: "Planck time" },
  planckFrequency: { value: 1 / 5.391_247e-44, unit: "Hz", label: "Planck frequency" },
} as const;

export const MEV_TO_J = 1e6 * constants.e.value;
