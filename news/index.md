# Changelog

## IRTsimrel 0.2.0

### Breaking Changes

- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  is now the primary function name for Stochastic Approximation
  Calibration (was
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)).
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
  remains available as a deprecated alias.
- [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  now supports `"info"`/`"tilde"` targets only. Direct MSEM-based
  targets (`"msem"`/`"bar"`) are rejected because the MSEM objective can
  be non-monotone under EQC root-finding; use
  [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  for direct MSEM targeting.
- [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  first parameter renamed from `eqc_result` to `result` and now accepts
  both `eqc_result` and `sac_result` objects.
- Default `item_source` for
  [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  and
  [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  changed from `"irw"` to `"parametric"` so core workflows do not
  require the optional external IRW integration.

### New Features

- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  now stores calibrated item parameters (`beta_vec`, `lambda_base`,
  `lambda_scaled`, `items_base`, `items_calib`, `theta_quad`) in the
  result object, enabling direct use with
  [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md).
- [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  replaces
  [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
  as the primary comparison function (old name remains as deprecated
  alias).
- [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  now records calibration and simulation provenance, including result
  class, metric, calibration status, item source/design, calibration
  call, seed, sample size, and latent settings.
- [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  now returns richer diagnostics, including achieved reliability
  differences, metric/model/test-length metadata, and EQC/SAC status
  fields.

### Bug Fixes

- Fixed [`library()`](https://rdrr.io/r/base/library.html) calls inside
  package functions
  ([`plot.item_params()`](https://joonho112.github.io/IRTsimrel/reference/plot.item_params.md),
  [`plot.sac_result()`](https://joonho112.github.io/IRTsimrel/reference/plot.sac_result.md))
  — now use proper namespace qualification.
- Fixed S3 method signatures:
  [`as.numeric.latent_G()`](https://joonho112.github.io/IRTsimrel/reference/as.numeric.latent_G.md)
  and
  [`as.data.frame.item_params()`](https://joonho112.github.io/IRTsimrel/reference/as.data.frame.item_params.md)
  now include required `...` parameter.
- Seed handling now saves and restores `.Random.seed` through an
  internal helper.
- Fixed [`par()`](https://rdrr.io/r/graphics/par.html) side effects in
  base R plot methods.
- Fixed non-ASCII characters (μ, σ) in `sim_latentG.R`.

### Documentation

- Completed a vignette overhaul, expanding from 6 vignettes (~3,400
  lines) to 12 source vignettes (~10,200 lines) organized as
  getting-started material, two applied/methodological tracks, and a
  reference article:
  - **Getting Started**: `introduction`, `quick-start`
  - **Applied Researchers Track**: `applied-guide`,
    `latent-distributions`, `item-parameters`, `simulation-design`,
    `case-studies`
  - **Methodological Researchers Track**: `theory-reliability`,
    `algorithm-eqc`, `algorithm-sac`, `validation`
  - **Reference**: `api-reference`
- Added explicit figure alt text to rendered README and vignette
  figures.
- Added comprehensive roxygen documentation for all S3 methods.
- Updated `README.Rmd` with correct function names, metrics table,
  direct pkgdown article links, and the 2026 arXiv v2 reference.
- Added hex sticker logo (`man/figures/logo.png`).
- Updated pkgdown site navigation and reference index, including a
  dedicated deprecated-aliases section for
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
  and
  [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md).

### Design Improvements

- Proper [`summary()`](https://rdrr.io/r/base/summary.html) objects for
  `eqc_result`, `sac_result`, `item_params`, and `latent_G` (each with
  dedicated `print.summary.*` method).
- Verbose output now uses
  [`message()`](https://rdrr.io/r/base/message.html) instead of
  [`cat()`](https://rdrr.io/r/base/cat.html) (suppressible with
  [`suppressMessages()`](https://rdrr.io/r/base/message.html)).
- [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  now validates matching `target_rho`, warns on model/n_items/metric
  differences.
- Convenience auto-wrapping for `latent_params` shape parameters.
- Full [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) diagnostics
  stored in `misc$uniroot_result`.

### New Utilities

- [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
  screens achievable reliability ranges before calibration.
- [`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
  computes and plots reliability as a function of scaling factor.
- [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)
  computes both reliability metrics in a single pass.
- [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) S3 methods for
  `eqc_result` and `sac_result`.
- [`as.double()`](https://rdrr.io/r/base/double.html) is now available
  as a synonym for [`as.numeric()`](https://rdrr.io/r/base/numeric.html)
  on `latent_G` objects.

### Known Limitations

- Saved objects from IRTsimrel 0.1.x whose class is literally
  `"spc_result"` may not dispatch every S3 helper method directly. New
  calls to
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
  return `"sac_result"` objects through the migration alias, and
  [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  plus comparison helpers still accept legacy `"spc_result"` objects.
  Recreate legacy calibration objects with
  [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  if direct [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`coef()`](https://rdrr.io/r/stats/coef.html), or
  [`predict()`](https://rdrr.io/r/stats/predict.html) dispatch is
  needed.

### Testing

- Expanded the testthat suite across core functions, S3 methods, edge
  cases, and error handling.

### Infrastructure

- Fixed `.Rbuildignore` to exclude `dev/`, `log/`, `pkg.lock`, and
  release-process helper files from source package builds.
- Fixed `DESCRIPTION`: removed `LazyData`, `Remotes`; added `graphics`,
  `grDevices`, `utils` to Imports; added `MASS` to Suggests.
- Removed unused imports
  ([`stats::ecdf`](https://rdrr.io/r/stats/ecdf.html),
  [`stats::setNames`](https://rdrr.io/r/stats/setNames.html)) from
  NAMESPACE.
- Updated `_pkgdown.yml` with new function references, S3 helper method
  sections, and deprecated-alias grouping.
- Added a release checklist with GitHub release note and CRAN-comment
  guidance.
- Version bumped to 0.2.0.

## IRTsimrel 0.1.0

- Initial release with EQC and SPC calibration algorithms.
- 12 built-in latent distribution shapes with pre-standardization.
- IRW integration for realistic item difficulty generation.
- Copula method for correlated difficulty-discrimination generation.
- TAM-based validation utilities.
