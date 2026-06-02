# Cylindrical Helicoid Persistence Across Bases

## Scope

This exploratory visualization uses the `finite_lifetime_only` dataset. Stable rows are excluded so truncation assumptions do not drive the comparison. Visual grouping does not imply a physical relationship, validation, or discovery.

## Answers

- **Which base produces the clearest visible grouping?** By the declared family-compactness metric, base `phi` is the smallest-distance case (`6.621317`). This is a descriptive ranking, not a privileged base.
- **Do groups survive across bases?** Some local relationships survive, but angular arrangement changes visibly. The common `z = log10(N)` height preserves layering across every base.
- **Are the same neighbors repeated?** Yes. The table below lists pairs close in all selected bases.
- **Does base pi remain visually interesting?** Base `pi` has family compactness `6.644403`, cluster count `5`, and average nearest-neighbor distance `1.158664`. It is one exploratory view among several, not evidence of special status.
- **Is base 10 only a control or does it show structure?** Base `10` shows the clearest direct wrapped-log10 helicoid because its angle and height use the same logarithm. Its family compactness is `6.633493`. That structure is primarily coordinate-induced.
- **Does the cylindrical helicoid preserve more useful information than the 2D circle?** It preserves `log10(N)` height, which distinguishes particles that overlap angularly. Whether that is useful depends on the analysis question; it is not independent information because height is derived from `N`.
- **Is visible structure dominated by log10(N) height?** Largely yes. Height is unchanged across bases and often dominates Euclidean 3D separation. Base changes rotate points around the cylinder.

## Repeated Close Neighbors

| Particle A | Particle B | Fraction of bases close | Mean 3D distance |
| --- | --- | ---: | ---: |
| B0 | Bs0 | 1.000000 | 0.090481 |
| W | Z | 1.000000 | 0.305426 |
| B_plus | Bs0 | 1.000000 | 0.347101 |
| omega_782 | top | 1.000000 | 0.359216 |
| B_plus | B0 | 1.000000 | 0.433763 |
| J_psi | Higgs | 1.000000 | 0.507507 |

## Clustering Method

DBSCAN uses Euclidean helicoid coordinates with `eps = 1.6` and `min_samples = 2`. Close-pair tracking uses distance `< 1.6`. These are exploratory display-scale choices.

## Warning

The helicoid representation is exploratory. Visual grouping does not imply physical relation.

## Reproduction

```powershell
python -m src.cylindrical_helicoid_bases
```
