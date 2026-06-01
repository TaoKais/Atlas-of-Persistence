# Gap-Center Polygon Geometry

## Scope

This exploratory analysis asks whether the centers of the largest circular
gaps in the complex phase landscape form recognizable descriptive shapes. For
each logarithmic base (`2`, `e`, `pi`, `phi`, and `10`), it selects the largest
`3`, `4`, `5`, and `6` gaps, computes their angular midpoints, and connects the
midpoints in angular order on the unit circle.

The resulting polygon is a visualization aid. It is not an established
physical observable and does not claim discovery.

## What Gap-Center Geometry Means

Each midpoint marks the center of an under-populated angular interval. The
connected polygon summarizes how the largest empty sectors are distributed
around the transformed phase circle. A near-regular polygon indicates that
the selected gaps are spaced relatively evenly. A displaced centroid or a
dominant angular separation indicates a more uneven distribution.

The circular symmetry score compares the angular gaps between selected
midpoints with the ideal spacing `360 / K`:

```text
symmetry_error = mean(abs(actual_gap - expected_gap))
symmetry_score = 1 / (1 + symmetry_error)
```

Higher values indicate spacing closer to a regular polygon. The score is
descriptive and depends on the selected catalogue, logarithmic base, and `K`.

## Why Examine Gaps

Entity positions show where transformed persistence phases are occupied.
Gaps instead show where the catalogue has no nearby transformed phases. For
an exploratory landscape, the empty sectors can expose angular organization
that is less obvious in a dense scatter plot.

Gaps may change when catalogue membership changes, when input measurements are
revised, or when a different logarithmic base is selected. They should not be
treated as intrinsic physical boundaries.

## Current Numerical Comparison

For the current finite-lifetime catalogue, base `2` gives the strongest
descriptive symmetry score: `0.077859` for the top `5` gaps. Its centroid
radius is `0.111966`, and the rule-based label is `balanced radial structure`.

Base `e` produces the next clearest balanced case for the top `3` gaps, with a
symmetry score of `0.067578` and centroid radius `0.119540`.

The best score for each base is:

| base | best top K | best symmetry score | qualitative shape |
| --- | ---: | ---: | --- |
| 2 | 5 | 0.077859 | balanced radial structure |
| e | 3 | 0.067578 | balanced radial structure |
| pi | 6 | 0.061768 | irregular |
| phi | 6 | 0.040293 | asymmetric sectoring |
| 10 | 6 | 0.030115 | asymmetric sectoring |

Under these simple metrics, base `2` produces the clearest structure among the
five tested bases. The scores remain modest, and several configurations are
irregular or asymmetric. This is a catalogue-specific visual comparison, not
evidence that base `2` is physically preferred.

## Interpretation Limit

Connecting gap centers is exploratory. Recognizable polygons, cross-like
arrangements, triangular sectoring, or balanced radial patterns can arise from
the transform, sample selection, and ranking procedure. Visual geometry does
not imply a physical law, a resonance mechanism, or a causal explanation of
particle lifetime.

## Reproduction

Run:

```powershell
python -m src.gap_geometry
python -m unittest discover -s tests -v
```

The summary table is written to `data/gap_geometry_summary.csv`. Figures are
written to `outputs/gap_geometry/`.
