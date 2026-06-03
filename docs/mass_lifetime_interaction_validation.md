# Mass-Lifetime Interaction Validation

Tolerances used: coefficient tolerance +/-0.05, p-value threshold 0.05, bootstrap confidence 0.95, weak exponent reference -5.0 +/-2.0, minimum regression group n=3, permutation runs=10000.

## Section 1 - What was tested

This validation tests the mass-lifetime correlation for finite-lifetime unstable particles only, with n=29. Stable particles, missing lifetimes, infinite lifetimes, and manually truncated stable values are excluded by default.

The Compton frequency is `f_C = m c^2 / h`, and persistence cycles are `N = f_C tau`. Because `f_C` is proportional to mass, correlations involving `N` and mass are partly tautological. Lifetime itself is not tautological because `tau` is measured independently of the definition of `f_C`.

Interaction type matters because decay width follows `Gamma = hbar / tau` and, schematically, `Gamma` depends on matrix elements and phase space. Weak three-body beta-like decays can show an approximate `Q^5` or muon-like `m^-5` lifetime scaling, but that exponent is not universal.

## Section 2 - Results

- Global Pearson rho: -0.555603, p=0.00169983, n=29, 95% bootstrap CI=[-0.780806, -0.268095], tolerance +/-0.05.
- Global Spearman rho: -0.517551, p=0.00449955, n=29, 95% bootstrap CI=[-0.783567, -0.151214], tolerance +/-0.05.
- Global regression slope alpha: -4.62306, n=29.
- Mass-only adjusted R2: 0.28309, AIC=105.118, BIC=107.853.
- Mass plus interaction adjusted R2: 0.78532, AIC=72.734, BIC=79.5705.
- Interaction permutation p-value for Model 3 vs Model 1: 9.999e-05, n=29, permutation runs=10000.

Per-interaction slopes:

- electromagnetic: status=insufficient_data, n=2
- mixed: status=fit, n=3, alpha=-1.13892, 95% CI=[-5.33437, 0.486555]
- strong: status=fit, n=4, alpha=-0.873408, 95% CI=[-17.5358, 301.729]
- weak: status=fit, n=20, alpha=-6.57648, 95% CI=[-8.21111, -4.08754]

## Section 3 - Physical interpretation

Mass tends to open phase space and often increases decay width, reducing lifetime. But lifetime is not determined by mass alone. It also depends on coupling constants, available decay channels, phase space, conservation laws, selection rules, and interaction type.

The weak exponent comparison uses alpha approximately -5 only as a loose reference motivated by muon-like weak decay scaling. It should not be applied as a universal weak-decay law.

## Section 4 - Tautology audit

Tautological or partly tautological:

- `N` includes mass through `f_C`.
- `logN = log_tau + log_mass + log10(MeV_J / h)`.
- Correlations between mass and `N` reuse mass by construction.

Not tautological:

- Correlation between mass and measured lifetime.
- Interaction-dependent residual structure after fitting lifetime against mass.
- Improvement of models when interaction labels are included, subject to permutation controls.

In this run, mass explains R2=0.308694 of `log_lifetime`, while mass explains the reported `logN` variance partly through the definitional mass term. Interaction labels explain additional variance only if the permutation p-value is below 0.05.

## Section 5 - Valid conclusion

After removing the definitional contribution of Compton frequency to N, the measured lifetime still shows a statistically meaningful relationship with mass and interaction class. This supports the interpretation that lifetime is governed by energy scale, phase space, coupling strength, and decay topology, not by mass alone.

This conclusion is evaluated with n=29, coefficient tolerance +/-0.05, p-value threshold 0.05, and bootstrap confidence 0.95.

## Section 6 - Warnings

- Small dataset warning: per-interaction groups are small, and groups with n < 3 are marked insufficient.
- PDG/catalog selection bias warning: the sample is curated and not a complete particle catalogue.
- Stable particle exclusion warning: stable particles are excluded from default regression to avoid artificial truncation.
- Mixed interaction ambiguity warning: `mixed` labels combine multiple mechanisms and should not be overinterpreted.
- No new physics claim: these are descriptive validation checks on known measured quantities.

## Reproduction

```powershell
python -m src.mass_lifetime_interactions
```
