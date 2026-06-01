# Persistence Framework

## Scope

This module investigates a descriptive question:

> Does `log10(f * tau)` reveal persistence structures that are not obvious when
> analyzing mass or lifetime independently?

It does not propose a new physical theory. The focused particle table uses
measured inputs and derived quantities. A visual gap or algorithmic cluster is
not evidence of a new physical law.

## Established Physics

### Energy and frequency

The Planck relation is

```text
E = h f
```

For rest energy `E = m c^2`, the corresponding Compton frequency is

```text
f_C = m c^2 / h
```

The implementation accepts mass-energy in `MeV` and converts it to joules
before dividing by the exact SI value of `h`.

### Particle lifetimes

Unstable particles have experimentally measured mean lifetimes or widths. A
width can be related to a lifetime through established quantum mechanics:

```text
tau = hbar / Gamma
```

Electron and proton decay have not been observed. Their rows use experimental
lower bounds, not invented finite lifetimes:

- Electron: total lifetime lower bound of `6.6e28 years`.
- Proton: partial-lifetime lower bound of `2.4e34 years` for
  `p -> e+ pi0`. This is not a measured total proton lifetime.

Consequently, electron and proton values of `N` and `log10(N)` are also lower
bounds.

### Quality-factor concepts

For a resonance, a quality factor compares a characteristic frequency with a
decay or damping scale. The exploratory persistence index is related in spirit
to cycle-count reasoning:

```text
N = f * tau
```

For the Compton-frequency particle case, `2 pi N` is related to the ratio of
rest energy to width. This does not make `N` a causal explanation of lifetime.

## Exploratory Concepts

### Persistence cycles

`N = f * tau` is dimensionless. It counts characteristic cycles completed
during a selected lifetime, coherence time, or persistence duration.

### Persistence landscapes

A persistence landscape plots `log10(N)` across a selected catalogue. The
logarithm supports inspection across many orders of magnitude.

### Stability islands

The analysis sorts `log10(N)` values and computes adjacent gaps:

```text
Delta_i = log10(N_(i+1)) - log10(N_i)
```

For visualization only, the generated clustering analysis splits the ordered
sample where a gap exceeds `mean(gaps) + population_stddev(gaps)`. This is a
transparent descriptive rule, not a validated physical classifier.

### Frequency-defined identities

The same calculation can be applied to measured frequencies from different
systems. `FrequencyType` records the selected definition:

- `COMPTON`
- `SCHUMANN`
- `ORBITAL`
- `PLASMA`
- `CUSTOM`

This supports particles, atoms, planets, stars, and resonant systems when a
measured characteristic frequency and measured duration are available. Values
from different frequency definitions must not be compared without stating the
definition and physical context.

## Numerical Results

Run:

```powershell
python -m src.persistence
```

The generator writes:

- `data/persistence_table.csv`
- `output/persistence/persistence_gaps.csv`
- `output/persistence/persistence_clusters.csv`
- `output/persistence/landscape_metrics.csv`
- `output/persistence/summary.md`
- `output/persistence/*.svg`

In the current eight-particle focused sample, the largest adjacent
`log10(N)` gap is between neutron and electron, with `Delta = 30.1105`. The
electron value is a lower bound, so the observed gap is not an exact distance.
The visualization is useful for forming questions, but the sample does not
support a discovery claim. `landscape_metrics.csv` compares log-scale ranges
and largest adjacent gaps for mass, lifetime, and `N`. This makes the research
question numerically inspectable without treating `N` as an independent
predictor: `N` contains the selected lifetime by construction.

## Data Sources

- [PDG particle listings, 2025 update](https://pdg.lbl.gov/2025/listings/contents_listings.html)
- [Electron lifetime bound, Borexino](https://doi.org/10.1103/PhysRevLett.115.231802)
- [Proton partial-lifetime bound, Super-Kamiokande](https://doi.org/10.1103/PhysRevD.102.112011)
- [NIST CODATA constants](https://physics.nist.gov/cuu/Constants/)

The source for each particle input is also recorded in
`data/persistence_particles.csv`.
