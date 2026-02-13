# Changelog

## IRTsimrel 0.2.0

### Breaking Changes

- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  is now the primary function name for Stochastic Approximation
  Calibration (was
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)).
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  remains available as a deprecated alias.
- [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  default reliability metric changed from `"msem"` to `"info"` (see Lee,
  2025, Section 4.3 on non-monotone MSEM objective for EQC).
- [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
  first parameter renamed from `eqc_result` to `result` and now accepts
  both `eqc_result` and `sac_result` objects.
- Default `item_source` for
  [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
  and
  [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  changed from `"irw"` to `"parametric"` (since irw is an optional
  Suggests dependency).

### New Features

- [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  now stores calibrated item parameters (`beta_vec`, `lambda_base`,
  `lambda_scaled`, `items_base`, `items_calib`, `theta_quad`) in the
  result object, enabling direct use with
  [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md).
- [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  replaces
  [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
  as the primary comparison function (old name remains as deprecated
  alias).

### Bug Fixes

- Fixed [`library()`](https://rdrr.io/r/base/library.html) calls inside
  package functions
  ([`plot.item_params()`](https://joonho112.github.io/IRTsimrel/reference/plot.item_params.md),
  `plot.spc_result()`) — now use proper namespace qualification.
- Fixed S3 method signatures:
  [`as.numeric.latent_G()`](https://joonho112.github.io/IRTsimrel/reference/as.numeric.latent_G.md)
  and
  [`as.data.frame.item_params()`](https://joonho112.github.io/IRTsimrel/reference/as.data.frame.item_params.md)
  now include required `...` parameter.
- Fixed [`set.seed()`](https://rdrr.io/r/base/Random.html) side effects
  — all 6 instances now properly save and restore `.Random.seed` via
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html).
- Fixed [`par()`](https://rdrr.io/r/graphics/par.html) side effects in
  base R plot methods.
- Fixed non-ASCII characters (μ, σ) in `sim_latentG.R`.

### Documentation

- Complete vignette overhaul: expanded from 6 vignettes (~3,400 lines)
  to 12 vignettes (~10,200 lines) organized in a two-track structure:
  - **Getting Started**: `introduction`, `quick-start`
  - **Applied Researchers Track**: `applied-guide`,
    `latent-distributions`, `item-parameters`, `simulation-design`,
    `case-studies`
  - **Methodological Researchers Track**: `theory-reliability`,
    `algorithm-eqc`, `algorithm-sac`, `validation`, `api-reference`
- All vignettes use `eval=TRUE` code chunks with working examples.
- Added comprehensive roxygen documentation for all S3 methods.
- Created `NEWS.md`.
- Updated `README.Rmd` with correct function names, metrics table, and
  arXiv reference.
- Added hex sticker logo (`man/figures/logo.png`).
- Updated pkgdown site with two-track article navigation and
  comprehensive reference index.

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

### Testing

- Added comprehensive test suite with testthat (519 tests covering all
  core functions, S3 methods, edge cases, and error handling).

### Infrastructure

- Fixed `.Rbuildignore` to exclude `dev/`, `pkg.lock`, `.claude/`.
- Fixed `DESCRIPTION`: removed `LazyData`, `Remotes`; added `graphics`,
  `grDevices`, `utils` to Imports; added `MASS` to Suggests.
- Removed unused imports
  ([`stats::ecdf`](https://rdrr.io/r/stats/ecdf.html),
  [`stats::setNames`](https://rdrr.io/r/stats/setNames.html)) from
  NAMESPACE.
- Updated `_pkgdown.yml` with new function references and S3 method
  sections.
- Updated `README.Rmd` with Phase 6 features, corrected metrics table,
  and arXiv URL.
- Version bumped to 0.2.0.
- R CMD check passes with 0 errors, 0 warnings, 0 notes.

## IRTsimrel 0.1.0

- Initial release with EQC and SPC calibration algorithms.
- 12 built-in latent distribution shapes with pre-standardization.
- IRW integration for realistic item difficulty generation.
- Copula method for correlated difficulty-discrimination generation.
- TAM-based validation utilities.
