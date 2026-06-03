# License Migration Report

## Migration Target

License: MIT

## Previous License Detected

- Repository root: no existing `LICENSE`, `COPYING`, or `NOTICE` file was detected.
- First-party package metadata: `archipelago-explorer/package.json` did not declare a license.
- First-party README/docs/source scan: no GPL, LGPL, AGPL, Apache-2.0, BSD, MPL, custom license block, `SPDX-License-Identifier`, or copyright header was detected in first-party source files.
- Third-party dependency metadata in `archipelago-explorer/package-lock.json` contains dependency license declarations that are not repository license declarations and were preserved.

## Files Modified

- `LICENSE`
- `README.md`
- `README.en.md`
- `README.es.md`
- `archipelago-explorer/README.md`
- `archipelago-explorer/package.json`
- `archipelago-explorer/package-lock.json`
- `docs/LICENSE_MIGRATION_REPORT.md`
- `docs/LICENSE_AUDIT.md`

## Source Header Handling

No first-party source files contained license headers. No code functionality was modified.

## Files Requiring Manual Review

- `data/*.csv`: dataset rows contain external source references such as PDG and NIST/CODATA. These data-source references are attribution/provenance references, not repository license grants.
- `outputs/**` and `output/**`: generated figures, animations, and reports should be regenerated from MIT-licensed first-party code when release provenance matters.
- `archipelago-explorer/package-lock.json`: dependency license metadata includes permissive and attribution licenses that remain attached to third-party packages.
- Python requirements are unpinned (`matplotlib>=3.10`, `numpy>=2.4`, `pandas>=3.0`, `imageio-ffmpeg>=0.6`), so exact transitive dependency licenses can vary by resolver date and platform.

## Possible Third-Party Licensing Conflicts

- No GPL, LGPL, AGPL, or MPL dependency declaration was found in the checked manifests or npm lockfile.
- `archipelago-explorer/package-lock.json` includes permissive third-party licenses: MIT, ISC, BSD-2-Clause, BSD-3-Clause, Apache-2.0, 0BSD, Unlicense, Zlib, BlueOak-1.0.0.
- `archipelago-explorer/package-lock.json` includes `CC-BY-4.0` for `caniuse-lite`; this is an attribution license and should be retained in third-party notices for redistributed dependency bundles.
- `archipelago-explorer/package-lock.json` includes `SEE LICENSE IN LICENSE.txt` for `@plotly/mapbox-gl` and `mapbox-gl`; inspect those package license files before redistributing bundled artifacts.
- No first-party file was identified as legally blocked from MIT relicensing during this scan. Third-party packages and external datasets are not relicensed by this migration.

## Validation

Validation scan targets:

- `LICENSE*`, `COPYING*`, `NOTICE*`, `README*`
- package metadata and lockfiles
- `docs/**`
- `src/**`, `tests/**`, `archipelago-explorer/src/**`

Expected remaining non-MIT license references are limited to third-party dependency metadata in `archipelago-explorer/package-lock.json` and this audit documentation.
