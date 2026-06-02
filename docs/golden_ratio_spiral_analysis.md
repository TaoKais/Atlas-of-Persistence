# Golden Ratio Spiral Robustness Analysis

## Scope

The golden ratio is not assumed physically relevant. This analysis is falsification-oriented:
the goal is to test whether phi is special, not to prove it. The score is heuristic only and
must not be presented as physical evidence.

## Questions

1. **Does phi produce unusually coherent geometry?** Phi ranks 4962 of 5017
   scanned bases by the declared heuristic score. This does not establish unusual coherence.
2. **Is the apparent spiral actually logarithmic?** Under the fixed `sorted_log10_N`
   ordering, the sparse radial projection has R^2 = 0.0382.
   The fit depends on ordering and should not be overinterpreted.
3. **Is it close to a golden spiral?** Phi's fitted k is -0.109693; golden k is
   0.306349; relative error is 1.3581. Phi ranks
   4967 of 5017 by golden-spiral error, so the scan directly tests whether
   other bases fit at least as well.
4. **Does the pattern survive randomized controls?** The control table records empirical
   score p-values: random_angular_phase=0.871, random_logN_distribution=0.743, shuffle_N_values=0.010, shuffle_lifetimes=0.713, shuffle_masses=0.762.
   Similar randomized scores support an artifact interpretation.
5. **Does the result depend on stable particle handling?** Yes. See
   `data/golden_ratio_stable_mode_metrics.csv`; truncations are numerical assumptions, while
   infinity markers are excluded from finite logarithmic fits.
6. **Which bases outperform phi?** The smallest golden-spiral errors include:

| base | golden_spiral_error |
| --- | --- |
| 9.33171 | 1.31983e-06 |
| 13.6758 | 0.000144793 |
| 9.65385 | 0.00023746 |
| 16.4677 | 0.000312289 |
| 10.8741 | 0.000317477 |
| 9.66361 | 0.000367999 |
| 13.6855 | 0.000445823 |
| 9.34147 | 0.000457924 |

7. **Geometry, artifact, or potentially interesting signal?** The conservative interpretation
   is visualization geometry whose appearance is sensitive to base, point ordering, and stable
   handling. This exploratory analysis does not claim discovery or physical significance.

## Ordering Sensitivity At Phi

| ordering | fitted_k | golden_spiral_relative_error | spiral_fit_R2 |
| --- | --- | --- | --- |
| sorted_theta | 0.138356 | 0.548371 | 0.042576 |
| sorted_log10_N | -0.109693 | 1.35806 | 0.0381602 |
| persistent_neighbor_path | 0.275883 | 0.0994478 | 0.392793 |
| nearest_neighbor_3d | -0.168139 | 1.54885 | 0.0798063 |

## Reproduce

```bash
python -m src.golden_ratio_spiral_analysis
```
