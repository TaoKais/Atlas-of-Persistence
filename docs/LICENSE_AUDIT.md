# License Audit

## Repository License

License: MIT

The repository now contains a root `LICENSE` file with the standard MIT License text and first-party metadata/readme references identify the repository license as MIT.

## Dependency License Summary

Python direct requirements:

| Dependency | Requirement | License assessment |
| --- | --- | --- |
| `matplotlib` | `>=3.10` | Permissive-style project license; exact transitive set depends on resolver output. |
| `numpy` | `>=2.4` | BSD-style permissive license. |
| `pandas` | `>=3.0` | BSD-style permissive license. |
| `imageio-ffmpeg` | `>=0.6` | BSD-style permissive project license; bundled ffmpeg provenance should be checked for binary redistribution. |

Npm dependency lockfile summary from `archipelago-explorer/package-lock.json`:

| License | Package count |
| --- | ---: |
| MIT | 322 |
| ISC | 47 |
| BSD-3-Clause | 24 |
| BSD-2-Clause | 9 |
| Apache-2.0 | 4 |
| BlueOak-1.0.0 | 2 |
| SEE LICENSE IN LICENSE.txt | 2 |
| 0BSD | 1 |
| CC-BY-4.0 | 1 |
| Unlicense | 1 |
| Zlib | 1 |

No GPL, LGPL, AGPL, or MPL dependency declaration was found in the checked manifests or npm lockfile.

## Remaining Risks

- This is a repository license migration, not legal advice.
- Unpinned Python requirements may resolve to different transitive dependency sets over time.
- `@plotly/mapbox-gl` and `mapbox-gl` use `SEE LICENSE IN LICENSE.txt`; inspect bundled package license files before redistribution.
- `caniuse-lite` is listed as `CC-BY-4.0`; preserve attribution notices when redistributing bundles containing it.
- External scientific data references in CSV files and documentation remain source/provenance references and are not converted into MIT-licensed third-party data.
- Generated outputs should be regenerated for release builds if provenance or reproducibility is required.

## Migration Status

Complete for first-party repository files found in this scan.

Repository license = MIT.
