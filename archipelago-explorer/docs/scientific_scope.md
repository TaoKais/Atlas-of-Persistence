# Scientific Scope

## Purpose

Archipelago Explorer is an exploratory visualization tool. It does not claim
to discover new particles, new laws, hidden symmetries, or physical
relationships. It provides an interface for checking whether visual patterns
survive changes of formula, logarithmic base, dataset selection, stable-particle
handling, and randomized controls.

## Standard Physics

The standard formula catalog includes established relations such as the Planck
relation, Compton frequency, reduced Compton angular frequency, lifetime-width
relation, de Broglie wavelength, Schwarzschild compactness, Schwarzschild
radius, and a gravitational-redshift approximation. Each entry is labeled
`standard` and includes a units note.

## Exploratory Representations

The following are explicitly exploratory descriptors and visual mappings:

- Persistence cycles: `N = f_C tau`.
- Persistence logarithm: `P = log10(N)`.
- Circular phase: `theta = 2 pi frac(log_b(N))`.
- Complex phase and unit-circle plots.
- Cylindrical and radial helicoids.
- Circular gaps, gap-center polygons, and persistent-neighbor scores.

These can be mathematically useful without being physically fundamental.

## Why Persistence Cycles Are Dimensionless

`f_C` has units of `s^-1` and `tau` has units of `s`. Their product `N = f_C
tau` is therefore dimensionless. It can be interpreted as a cycle count under
the chosen characteristic frequency.

## Projection Artifacts

Phase projection deliberately folds a logarithmic axis around a circle. A
helix appears when the wrapped phase is plotted against the original
logarithmic height because winding is encoded by the coordinate definition.
Patterns that also appear after shuffled-lifetime, shuffled-mass, or randomized
`N` controls are likely mathematical artifacts.

## Stable Particles

Stable particles do not have a finite measured total lifetime. The app requires
an explicit handling choice:

- Exclude stable rows.
- Use an available lower bound, with a manual fallback.
- Assign a manual exploratory lifetime.
- Mark infinity and exclude the row from finite plots.

Manual values are assumptions, not measurements.

## Interpretation Limit

Visual grouping is not evidence. Any claim of physical significance requires a
predeclared hypothesis, uncertainty treatment, appropriate controls,
independent datasets, and review by domain experts.
