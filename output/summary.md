# Generated exploratory summary

## Dataset

- Particles: 31
- Stable identities represented without an assigned decay lifetime: 2
- Unstable identities: 29

## Immediate checks

- `f_C = mc^2/h` is exactly proportional to mass; its isolated ranking cannot
  improve on mass ranking.
- Largest sampled logarithmic gap: `electron` to
  `muon`, ratio `206.768`.
- Shortest sampled lifetime: `Z`, `2.64e-25 s`.
- Longest finite sampled lifetime: `neutron`, `878.4 s`.
- Schwarzschild horizon reference: `Phi = 0.5`;
  equivalently `2 Phi = 1`.

## Out-of-sample lifetime validation

- Validation: leave-one-out ridge regression over unstable identities only.
- Mass-only MAE: `4.884` log10 seconds.
- Mass plus representative modes MAE:
  `4.995` log10 seconds.
- Mass, modes, family, and dominant interaction MAE:
  `2.878` log10 seconds.
- In this sample, the curated mode-count proxy does not improve on mass alone.
- Dominant interaction improves prediction, but it may encode information close
  to the decay mechanism and is not evidence of a new persistence law.

## Interpretation limits

- Gaps depend on catalogue selection and observational bias.
- Stable particles need lower bounds, not invented finite lifetimes.
- `representative_decay_mode_count` is a curated exploratory annotation, not a
  complete channel count or phase-space volume.
- `N_C = f_C tau` contains the target lifetime. It is descriptive and must not
  be used as a lifetime predictor.
- The exergy scenarios are formula demonstrations, not fitted experiments.
