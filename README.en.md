# Atlas of Persistence

**English** | [Espanol](README.es.md)

Exploratory mathematical and computational framework for studying identity
persistence, characteristic frequency, stability, gravitational compactness,
and exergy. The project generates testable hypotheses and reproducible
calculations; it does not propose a demonstrated physical theory.

## Minimal Conceptual Result

The Compton frequency

```text
f_C = m c^2 / h
```

is proportional to mass. Therefore, sorting particles by `f_C` or by mass
produces exactly the same order and the same logarithmic gaps. `f_C` is a valid
physical coordinate, but it cannot classify persistence better than mass when
used in isolation.

The atlas explores parameters that do incorporate additional information:

| Parameter | Definition | Exploratory interpretation |
| --- | --- | --- |
| `N_C` | `f_C * tau` | Compton cycles during the mean lifetime |
| `Q_C` | `2 pi f_C tau` | inverse relative width `Gamma/(mc^2)` |
| `rho_E` | `E_total/(mc^2)` | total relativistic energy relative to rest energy |
| `Phi` | `GM/(Rc^2)` | gravitational compactness |
| `delta_X` | `X_destroyed/X_input` | fraction of destroyed exergy |

`N_C` and `Q_C` are not causal explanations of stability: they are
dimensionless descriptors derived from mass and mean lifetime. The decay
channel count included in the data is an exploratory variable, not a rigorous
phase-space volume.

## Execution

Requires Python 3.10 or later and installs no third-party dependencies:

```powershell
python -m src.atlas_persistence
python -m unittest discover -s tests -v
```

The first command regenerates the CSV files and summary in `output/`.

## Structure

```text
data/       traceable, editable input data
src/        calculations and CLI
tests/      basic mathematical and physical checks
output/     regenerable reports
```

## Data and Sources

The particle values are a pedagogical selection of central values published by
the Particle Data Group (PDG). Constants are based on CODATA 2022/NIST. Each CSV
includes a `source` column.

- [PDG particle listings, 2025 update](https://pdg.lbl.gov/2025/listings/contents_listings.html)
- [PDG summary tables, 2025 update](https://pdg.lbl.gov/2025/tables/contents_tables.html)
- [NIST CODATA values](https://physics.nist.gov/cuu/Constants/)

The exergy scenarios and some compact objects are reference examples for
validating formulas, not fitted observations.

## Falsifiable Questions

1. After controlling for mass, charge, spin, and conservation laws, determine
   whether an accessible-channel descriptor improves out-of-sample lifetime
   prediction.
2. Compare observed logarithmic gaps with simulated catalogues that account for
   experimental bias. A visual gap alone does not establish a privileged
   frequency.
3. Evaluate whether `delta_X` predicts loss of reconstruction capacity in
   concrete thermodynamic processes better than isolated energy efficiency.
4. Treat `Phi` as a regime coordinate and test where a Newtonian approximation
   exceeds a predefined error tolerance.

