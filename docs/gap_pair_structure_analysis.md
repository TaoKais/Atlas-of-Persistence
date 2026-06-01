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

## Frequency-Specific Check

Compton frequency is:

```text
f_C = m c^2 / h
```

It expresses the rest-energy scale as a frequency. It is not independent of
mass, so the Compton-frequency ratios exactly reproduce the mass ratios.

The frequency-specific output also measures each pair separation in octaves:

```text
octave_separation = log2(max(f_1 / f_2, f_2 / f_1))
```

| pair | folded frequency ratio | octave separation | nearest integer octaves | residual |
| --- | ---: | ---: | ---: | ---: |
| eta -> Higgs | 228.5247 | 7.836206 | 8 | 0.163794 |
| K_long -> Z | 183.2516 | 7.517682 | 8 | 0.482318 |
| pi_plus -> rho_770 | 5.554617 | 2.473687 | 2 | 0.473687 |
| B_plus -> Lambda | 4.731998 | 2.242450 | 2 | 0.242450 |
| K_short -> neutron | 1.888152 | 0.916975 | 1 | 0.083025 |

There is a coarse descriptive grouping near `8`, `2`, and `1` octaves, but
there is no single repeated interval. Two residuals are close to half an
octave, so this does not support a repeated octave-harmonic relationship.

The shared frequency property is therefore limited and expected: all five
pairs can be expressed through rest-energy Compton frequencies, and those
frequencies inherit the same ratios as mass. The current data do not show an
additional common frequency law linking the five gaps.

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
The frequency-specific comparison is written to
`data/gap_pair_frequency_structure.csv`.
