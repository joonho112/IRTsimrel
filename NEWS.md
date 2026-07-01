# IRTsimrel 0.2.0

## Breaking Changes
* `sac_calibrate()` is now the primary function name for Stochastic Approximation
  Calibration (was `spc_calibrate()`). `spc_calibrate()` remains available as a
  deprecated alias.
* `eqc_calibrate()` now supports `"info"`/`"tilde"` targets only. Direct
  MSEM-based targets (`"msem"`/`"bar"`) are rejected because the MSEM objective
  can be non-monotone under EQC root-finding; use `sac_calibrate()` for direct
  MSEM targeting.
* `simulate_response_data()` first parameter renamed from `eqc_result` to `result`
  and now accepts both `eqc_result` and `sac_result` objects.
* Default `item_source` for `eqc_calibrate()` and `sac_calibrate()` changed from
  `"irw"` to `"parametric"` so core workflows do not require the optional
  external IRW integration.

## New Features
* `sac_calibrate()` now stores calibrated item parameters (`beta_vec`, `lambda_base`,
  `lambda_scaled`, `items_base`, `items_calib`, `theta_quad`) in the result object,
  enabling direct use with `simulate_response_data()`.
* `compare_eqc_sac()` replaces `compare_eqc_spc()` as the primary comparison function
  (old name remains as deprecated alias).
* `simulate_response_data()` now records calibration and simulation provenance,
  including result class, metric, calibration status, item source/design,
  calibration call, seed, sample size, and latent settings.
* `compare_eqc_sac()` now returns richer diagnostics, including achieved
  reliability differences, metric/model/test-length metadata, and EQC/SAC
  status fields.

## Bug Fixes
* Fixed `library()` calls inside package functions (`plot.item_params()`,
  `plot.sac_result()`) — now use proper namespace qualification.
* Fixed S3 method signatures: `as.numeric.latent_G()` and
  `as.data.frame.item_params()` now include required `...` parameter.
* Seed handling now saves and restores `.Random.seed` through an internal helper.
* Fixed `par()` side effects in base R plot methods.
* Fixed non-ASCII characters (μ, σ) in `sim_latentG.R`.

## Documentation
* Completed a vignette overhaul, expanding from 6 vignettes (~3,400 lines) to
  12 source vignettes (~10,200 lines) organized as getting-started material,
  two applied/methodological tracks, and a reference article:
  - **Getting Started**: `introduction`, `quick-start`
  - **Applied Researchers Track**: `applied-guide`, `latent-distributions`,
    `item-parameters`, `simulation-design`, `case-studies`
  - **Methodological Researchers Track**: `theory-reliability`, `algorithm-eqc`,
    `algorithm-sac`, `validation`
  - **Reference**: `api-reference`
* Added explicit figure alt text to rendered README and vignette figures.
* Added comprehensive roxygen documentation for all S3 methods.
* Updated `README.Rmd` with correct function names, metrics table, direct
  pkgdown article links, and the 2026 arXiv v2 reference.
* Added hex sticker logo (`man/figures/logo.png`).
* Updated pkgdown site navigation and reference index, including a dedicated
  deprecated-aliases section for `spc_calibrate()` and `compare_eqc_spc()`.

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
* `as.double()` is now available as a synonym for `as.numeric()` on
  `latent_G` objects.

## Known Limitations
* Saved objects from IRTsimrel 0.1.x whose class is literally `"spc_result"` may
  not dispatch every S3 helper method directly. New calls to `spc_calibrate()`
  return `"sac_result"` objects through the migration alias, and
  `simulate_response_data()` plus comparison helpers still accept legacy
  `"spc_result"` objects. Recreate legacy calibration objects with
  `sac_calibrate()` if direct `summary()`, `plot()`, `coef()`, or `predict()`
  dispatch is needed.

## Testing
* Expanded the testthat suite across core functions, S3 methods, edge cases,
  and error handling.

## Infrastructure
* Fixed `.Rbuildignore` to exclude `dev/`, `log/`, `pkg.lock`, and
  release-process helper files from source package builds.
* Fixed `DESCRIPTION`: removed `LazyData`, `Remotes`; added `graphics`, `grDevices`,
  `utils` to Imports; added `MASS` to Suggests.
* Removed unused imports (`stats::ecdf`, `stats::setNames`) from NAMESPACE.
* Updated `_pkgdown.yml` with new function references, S3 helper method sections,
  and deprecated-alias grouping.
* Added a release checklist with GitHub release note and CRAN-comment guidance.
* Version bumped to 0.2.0.

# IRTsimrel 0.1.0

* Initial release with EQC and SPC calibration algorithms.
* 12 built-in latent distribution shapes with pre-standardization.
* IRW integration for realistic item difficulty generation.
* Copula method for correlated difficulty-discrimination generation.
* TAM-based validation utilities.
