# Stable Particle Sensitivity

## Scope

This sensitivity test compares finite-lifetime particles with controlled stable-particle truncations. A truncation is a numerical assumption, not a measured lifetime. No physical significance is claimed.

## Variants

- `finite_lifetime_only`: excludes rows marked stable and excludes missing lifetimes.
- `stable_truncated`: includes stable particles with `tau_stable` values of `1e25`, `1e30`, `1e35`, and `1e40` seconds. Plots use `1e35 s` unless stated otherwise.

## Answers

- **Was electron included in previous centroid plots?** Yes. The preceding dynamic base scan included electron and proton using documented lifetime lower bounds. This test replaces that handling with explicit finite-only and controlled-truncation variants.
- **Does including electron dominate the lepton centroid?** It materially changes the lepton centroid. The mean absolute circular shift across the five bases is `32.669` degrees. With only muon and tau in the finite-only variant, adding electron can substantially redirect the three-particle centroid.
- **Does proton dominate the baryon centroid?** It affects the baryon centroid, with a mean absolute circular shift of `28.016` degrees across the five bases. The term `dominate` is not consistently justified because the direction and magnitude vary by base.
- **Are previous structures robust without stable particles?** They are only partially robust. The continuous top-4 symmetry curves have Pearson correlation `0.365` and mean absolute score difference `0.007064`. Stable handling changes some peak details.
- **Which conclusions survive both variants?** The phase layout remains strongly base-dependent; selected-base geometry remains exploratory; and stable-particle treatment must be stated explicitly whenever family centroids or gap symmetry are interpreted.

## Interpretation Limit

This is a numerical sensitivity analysis. It does not establish new physics, a preferred base, or a stable-particle lifetime model.

## Reproduction

```powershell
python -m src.stable_particle_sensitivity
```
