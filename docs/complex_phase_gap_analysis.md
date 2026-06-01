# Complex Phase Gap Landscape

## Scope

This exploratory analysis asks whether mapping particle persistence cycles into
complex phase space reveals descriptive clustering or recurring angular gaps:

```text
N = f_C * tau
log_b(N) = log(N) / log(b)
phase_fraction = log_b(N) mod 1
z = exp(i * 2 pi * phase_fraction)
```

The secondary question is whether particles with similar dominant interaction
or stability class occupy similar phase regions for any selected base. The
calculation does not claim discovery.

## Why Add Phase

`log10(N)` is useful because particle persistence spans many orders of
magnitude. A linear ordering of logarithmic distances can still hide the
fractional part of `log_b(N)`. Wrapping that fractional part onto a unit circle
tests whether entities appear near similar positions modulo powers of a chosen
base.

Complex phase may reveal candidate harmonic or resonant-looking structure
because values separated by integer powers of a base map to the same angle.
That is a mathematical property of the transform. It is not evidence that a
particle system physically uses base `2`, `e`, `pi`, `phi`, or `10`.

## Established Physics Check

For particle rest energy:

```text
f_C = E / h = m c^2 / h
N = f_C * tau
```

When an unstable state's lifetime is represented by `tau = hbar / Gamma`, the
cycle count obeys:

```text
N = (E / h) * (hbar / Gamma) = E / (2 pi Gamma)
2 pi N = E / Gamma
```

This is consistent with the usual energy-to-width ratio and quality-factor
reasoning. The circular phase transform is an added exploratory
representation; it is not an established particle-physics observable.

Stable particles are omitted unless the input catalogue supplies an explicit
finite lower-bound lifetime. The broad catalogue currently omits electron and
proton from this phase analysis because their lifetime cells are blank.

## Circular Metrics

For each base, the mean resultant length is:

```text
R = abs(mean(exp(i * theta)))
```

`R` near `0` indicates dispersed phases. `R` near `1` indicates alignment.
Circular variance is `1 - R`.

Adjacent circular gaps are computed after sorting phase angles. A large gap is
an under-populated angular interval for the selected catalogue and base. It
may change with catalogue membership, measurement revisions, or the arbitrary
choice of base.

## Interpretation Limit

Phase alignment does not imply a physical law, resonance mechanism, preferred
logarithmic base, or causal explanation of particle lifetime. The output is a
numerical screening tool for forming questions that would require independent
physical justification and statistical testing on a pre-specified sample.

## Current Numerical Observation

For the current 29-entity finite-lifetime catalogue, `R` ranges from about
`0.026` for base `10` to `0.268` for base `phi`. This does not show strong
global phase alignment. The largest circular gap varies by base, from about
`31.9` degrees for base `10` to `79.5` degrees for base `2`.

The generated plots allow interaction-group patterns to be inspected, but this
small curated sample is not sufficient to conclude that similar interactions
or stability classes occupy preferred phase regions. A next statistical step
would compare pre-specified metrics against explicit null catalogues while
accounting for the multiple selected bases.

## Reproduction

Run:

```powershell
python -m src.complex_phase
python -m unittest discover -s tests -v
```

Tables are written to `data/complex_phase_*.csv`. Figures are written to
`outputs/complex_phase/`.
