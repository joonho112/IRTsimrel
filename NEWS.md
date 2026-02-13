# IRTsimrel 0.2.0

## Breaking Changes
* `sac_calibrate()` is now the primary function name for Stochastic Approximation
  Calibration (was `spc_calibrate()`). `spc_calibrate()` remains available as a
  deprecated alias.
* `eqc_calibrate()` default reliability metric changed from `"msem"` to `"info"`
  (see Lee, 2025, Section 4.3 on non-monotone MSEM objective for EQC).
* `simulate_response_data()` first parameter renamed from `eqc_result` to `result`
  and now accepts both `eqc_result` and `sac_result` objects.
* Default `item_source` for `eqc_calibrate()` and `sac_calibrate()` changed from
  `"irw"` to `"parametric"` (since irw is an optional Suggests dependency).

## New Features
* `sac_calibrate()` now stores calibrated item parameters (`beta_vec`, `lambda_base`,
  `lambda_scaled`, `items_base`, `items_calib`, `theta_quad`) in the result object,
  enabling direct use with `simulate_response_data()`.
* `compare_eqc_sac()` replaces `compare_eqc_spc()` as the primary comparison function
  (old name remains as deprecated alias).

## Bug Fixes
* Fixed `library()` calls inside package functions (`plot.item_params()`,
  `plot.spc_result()`) — now use proper namespace qualification.
* Fixed S3 method signatures: `as.numeric.latent_G()` and
  `as.data.frame.item_params()` now include required `...` parameter.
* Fixed `set.seed()` side effects — all 6 instances now properly save and restore
  `.Random.seed` via `on.exit()`.
* Fixed `par()` side effects in base R plot methods.
* Fixed non-ASCII characters (μ, σ) in `sim_latentG.R`.

## Documentation
* Complete vignette overhaul: expanded from 6 vignettes (~3,400 lines) to 12
  vignettes (~10,200 lines) organized in a two-track structure:
  - **Getting Started**: `introduction`, `quick-start`
  - **Applied Researchers Track**: `applied-guide`, `latent-distributions`,
    `item-parameters`, `simulation-design`, `case-studies`
  - **Methodological Researchers Track**: `theory-reliability`, `algorithm-eqc`,
    `algorithm-sac`, `validation`, `api-reference`
* All vignettes use `eval=TRUE` code chunks with working examples.
* Added comprehensive roxygen documentation for all S3 methods.
* Created `NEWS.md`.
* Updated `README.Rmd` with correct function names, metrics table, and arXiv reference.
* Added hex sticker logo (`man/figures/logo.png`).
* Updated pkgdown site with two-track article navigation and comprehensive
  reference index.

## Design Improvements
* Proper `summary()` objects for `eqc_result`, `sac_result`, `item_params`, and
  `latent_G` (each with dedicated `print.summary.*` method).
* Verbose output now uses `message()` instead of `cat()` (suppressible with
  `suppressMessages()`).
* `compare_eqc_sac()` now validates matching `target_rho`, warns on
  model/n_items/metric differences.
* Convenience auto-wrapping for `latent_params` shape parameters.
* Full `uniroot()` diagnostics stored in `misc$uniroot_result`.

## New Utilities
* `check_feasibility()` screens achievable reliability ranges before calibration.
* `rho_curve()` computes and plots reliability as a function of scaling factor.
* `compute_rho_both()` computes both reliability metrics in a single pass.
* `coef()` and `predict()` S3 methods for `eqc_result` and `sac_result`.

## Testing
* Added comprehensive test suite with testthat (519 tests covering all core functions,
  S3 methods, edge cases, and error handling).

## Infrastructure
* Fixed `.Rbuildignore` to exclude `dev/`, `pkg.lock`, `.claude/`.
* Fixed `DESCRIPTION`: removed `LazyData`, `Remotes`; added `graphics`, `grDevices`,
  `utils` to Imports; added `MASS` to Suggests.
* Removed unused imports (`stats::ecdf`, `stats::setNames`) from NAMESPACE.
* Updated `_pkgdown.yml` with new function references and S3 method sections.
* Updated `README.Rmd` with Phase 6 features, corrected metrics table, and arXiv URL.
* Version bumped to 0.2.0.
* R CMD check passes with 0 errors, 0 warnings, 0 notes.

# IRTsimrel 0.1.0

* Initial release with EQC and SPC calibration algorithms.
* 12 built-in latent distribution shapes with pre-standardization.
* IRW integration for realistic item difficulty generation.
* Copula method for correlated difficulty-discrimination generation.
* TAM-based validation utilities.
