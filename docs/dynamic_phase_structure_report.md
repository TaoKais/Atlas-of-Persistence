# Dynamic Phase Structure Scan

## Scope

This exploratory scan evaluates 31 particles across 5000 evenly spaced bases from `1.2` to `50.0`. Electron and proton use documented lifetime lower bounds. Their phases are therefore bound-based coordinates, not measured total-lifetime phases.

The projection is descriptive. It does not establish a preferred logarithmic base, resonance mechanism, physical significance, or new physics.

## Answers

- **Is pi genuinely special?** No statistically justified special status is established. Its local top-4 symmetry percentile is `0.282`.
- **Is phi genuinely special?** No statistically justified special status is established. Its local top-4 symmetry percentile is `0.063`.
- **Is base 2 genuinely special?** No statistically justified special status is established. Its local top-4 symmetry percentile is `0.942`.
- **Are observed structures robust?** The metrics vary repeatedly under continuous base changes. Any selected-base geometry must be interpreted against the full scan.
- **Are the largest gaps stable?** No. Their magnitudes and bounding particles change across the scan.
- **Do persistent clusters exist?** Recurring neighbor groups exist descriptively at the declared `20.0` degree threshold; they are dataset relationships, not physical resonances.
- **Does any base maximize organization beyond neighboring bases?** The strongest sampled top-4 symmetry occurs at `43.79139828` with score `0.794921`. This is a scan maximum, not evidence of significance; nearby-base and null-model testing would be required for a stronger claim.

## Scan Extremes

- Highest sampled top-4 symmetry: base `43.79139828`, score `0.794921`.
- Lowest sampled top-4 symmetry: base `4.76311262`, score `0.010818`.

## Most Persistent Neighbor Pairs

| Pair | Fraction of bases |
| --- | ---: |
| B0, Bs0 | 0.997200 |
| W, Z | 0.970800 |
| B_plus, Bs0 | 0.963200 |
| omega_782, top | 0.960600 |
| B_plus, B0 | 0.943000 |
| J_psi, Higgs | 0.920400 |
| pi0, eta | 0.348400 |
| Sigma_plus, pi_plus | 0.342200 |

## Bonus Projection Comparison

`theta = 2*pi*frac(ln(N))` and `theta = arg(exp(i*ln(N)))` are mathematically equivalent up to angle wrapping. They are also exactly the base-`e` member of `theta = 2*pi*frac(log_b(N))`. Circular geometry is therefore caused by the circular projection; its detailed layout changes with logarithmic base.

## Reproduction

```powershell
python -m src.dynamic_phase_structure
```
