# Gap-Pair Structure Analysis

## Scope

This follow-up examines the five endpoint pairs behind the strongest
gap-center configuration from the previous analysis: base `2`, top `5`.

The selected pairs are:

| from | to | circular gap (degrees) |
| --- | --- | ---: |
| eta | Higgs | 79.490362 |
| K_long | Z | 32.553305 |
| pi_plus | rho_770 | 29.978960 |
| B_plus | Lambda | 27.858678 |
| K_short | neutron | 26.779036 |

For each endpoint, the generated table records mass, Compton frequency,
lifetime, cycle count `N`, `log2(N)`, and `log10(N)`. It also records pairwise
ratios, signed differences, and geometric means.

## Shared-Scale Check

To test whether the five pairs share a multiplicative separation scale, each
ratio is folded to be at least one so that the comparison does not depend on
pair orientation. For each metric:

```text
multiplicative_spread = maximum_folded_ratio / minimum_folded_ratio
```

A conservative descriptive candidate flag requires the folded ratios to stay
within a factor of two. This threshold is a screening rule, not a statistical
test.

| metric | minimum ratio | maximum ratio | multiplicative spread | shared scale |
| --- | ---: | ---: | ---: | --- |
| mass | 1.888152 | 228.5247 | 121.0309 | no |
| Compton frequency | 1.888152 | 228.5247 | 121.0309 | no |
| lifetime | 159.7680 | 1.937879e17 | 1.212933e15 | no |
| N | 13.72937 | 1.062751e15 | 7.740713e13 | no |

Mass and Compton-frequency ratios match because Compton frequency is directly
proportional to mass. None of the four positive-valued metrics passes the
shared-scale screening rule.

## Interpretation

The five circular phase gaps do not correspond to one repeated physical
endpoint ratio. Their visual arrangement on the base-2 phase circle is a
modular-logarithmic pattern: it does not imply a common mass scale, lifetime
scale, or persistence-cycle scale.

This remains an exploratory catalogue-specific comparison. It does not claim
discovery, establish a physical law, or show that base `2` is physically
preferred.

## Reproduction

Run:

```powershell
python -m src.gap_pair_structure
python -m unittest discover -s tests -v
```

Detailed pair metrics are written to `data/gap_pair_structure.csv`. The
cross-pair scale summary is written to `data/gap_pair_structure_summary.csv`.
