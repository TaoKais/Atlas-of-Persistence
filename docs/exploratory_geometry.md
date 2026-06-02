# Exploratory 3D Persistence Geometry

## Scope

This study searches for unexpected geometry without validating a hypothesis. It uses the `finite_lifetime_only` dataset so stable-particle truncation does not drive the result. All structures are exploratory visual geometry only. No new physics is claimed.

## Why a Helix Appears

A cylindrical helix is induced by the coordinates. At base `10`, `z = log10(N)` and `theta = 2*pi*frac(log10(N))`, so increasing `z` necessarily winds around the unit circle. Other bases rescale the winding rate. Randomized controls retain a helix-like trace because the same wrapping rule is applied after shuffling. The helix itself is therefore primarily a visualization, logarithmic, and circular-projection artifact.

## What Geometries Appear?

- Cylindrical, radial, and conical views show winding layers because phase is a wrapped logarithm.
- Toroidal and spherical projections reveal occupancy bands and voids, but their locations move with the chosen projection and base.
- The density cloud shows frequently visited radial regions as base changes.
- Particle worldlines and family tubes orbit and cross through base space.
- Gap-center constellation symmetry changes with base rather than remaining fixed.

## What Disappears Under Controls?

Shuffling lifetimes or masses changes point ordering, local clusters, gaps, and neighbor identities. The broad helix-like winding survives because it is encoded by the transform. Dataset-specific occupancy patterns do not survive unchanged.

## Which Structures Are Robust?

The robust structure is the wrapped-logarithm winding itself. Local gaps, polygons, bands, and family trajectories are conditional on the dataset, base, and projection. Persistent neighbor pairs are the most reproducible dataset-specific relationships in this scan.

## Persistent Neighbor Pairs

| Pair | Fraction of bases |
| --- | ---: |
| B0, Bs0 | 0.995833 |
| W, Z | 0.970833 |
| B_plus, Bs0 | 0.962500 |
| omega_782, top | 0.958333 |
| B_plus, B0 | 0.941667 |
| J_psi, Higgs | 0.916667 |
| pi0, eta | 0.350000 |
| muon, K_plus | 0.337500 |

## Randomized Controls

| Variant | Particle count | Endpoint turn difference |
| --- | ---: | ---: |
| observed | 29 | -0.617267 |
| shuffled_lifetimes | 29 | -0.215888 |
| shuffled_masses | 29 | -0.530274 |

## Interpretation Limit

This is an exploratory geometry study. It does not establish a physical helix, preferred base, resonance, or new physics.

## Reproduction

```powershell
python -m src.exploratory_3d_geometry
```
