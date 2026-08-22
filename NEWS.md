# IRTsimrel 0.3.0

## 3PL support

* Added end-to-end support for the unidimensional three-parameter logistic
  model with logistic scaling constant `D = 1`. The shared model engine,
  reliability reducers, item generation, EQC, SAC, response simulation, and S3
  result methods now carry item lower asymptotes.
* Added fixed, beta, and uniform guessing generators through
  `guessing_params`; custom forms use `custom_params$guessing`. Values must
  satisfy `0 <= guessing < 1`.
* Clarified parameter names: item guessing is stored as `guessing` (mathematical
  `g_i`), while `c`/`c_star` is the global discrimination scale. Calibration
  applies `lambda_i(c) = c * lambda_i0` and never rescales guessing.
* The `guessing = 0` 3PL path is computationally identical to the 2PL path.
  Probability and Fisher-information calculations use stable log-domain
  kernels for extreme linear predictors.

## Calibration and estimand contracts

* EQC now scans the complete user-supplied interval on `log(c)`, records all
  resolved roots and local extrema, and applies an explicit root policy. EQC
  remains a fixed-form, average-information (`info`/`tilde`) calibrator.
* SAC now runs a topology preflight before Robbins--Monro updates, constrains
  the trajectory to the selected branch, and performs an independent final
  evaluation. New controls are appended after the historical argument prefix:
  `root_policy`, `preflight_controls`, and `evaluation_controls`.
* SAC explicitly distinguishes `resample_items = FALSE` (fixed-form estimand)
  from `resample_items = TRUE` (item-superpopulation estimand). Results store
  calibration/evaluation design metadata, achieved distributions and standard
  errors, branch diagnostics, and RNG provenance.
* EQC result warm starts are accepted only for matching metric/model/item count
  and the same fixed-form contract. Same-seed SAC replay and caller RNG
  restoration are enforced by independent random-number streams.
* `check_feasibility()` and `rho_curve()` use the shared topology scanner rather
  than assuming global monotonicity. Multiple roots remain possible for
  adversarial or poorly covered designs.

## Reliability, migration, and performance

* Reliability utilities now accept 3PL guessing vectors and optional quadrature
  weights, and return numerical/estimand diagnostics on request. The default
  scalar/list shapes remain backward compatible.
* New result schema fields include `guessing_vec`, `schema_version`,
  `item_scope`, `estimand_signature`, and `design_signature`. Legacy Rasch/2PL
  results without guessing are interpreted lazily as zero-guessing designs and
  are not mutated. Incomplete legacy `spc_result` objects now fail with an
  actionable instruction to rerun `sac_calibrate()`.
* `spc_calibrate()` and `compare_eqc_spc()` remain deprecated forwarding aliases.
  The historical 18-argument prefix of `sac_calibrate()` is unchanged; new
  controls were appended to preserve positional calls.
* Reliability reduction automatically chunks transient item-information
  matrices while preserving exact scalar results. Matrix-returning internal
  kernels retain their previous full-output contract.

## Validation boundaries and known limitations

* Analytic identities, finite-difference derivatives, extreme-predictor cases,
  independent holdouts, stochastic sensitivity runs, and external TAM
  probability/information calculations were used to validate the D = 1 3PL
  implementation.
* TAM 3PL validation is an EAP-only diagnostic in the validated environment.
  `TAM::tam.wle()` did not support the fitted `tam.mml.3pl` object, so IRTsimrel
  does not claim 3PL WLE support. The public `compute_reliability_tam()` helper
  continues to provide WLE and EAP output for its documented Rasch/2PL paths.
* Analytic information/MSEM reliability and fitted-score EAP reliability are
  different estimands; equality is neither assumed nor used as an acceptance
  criterion.
* Population MSEM reliability is non-integrable for the package's built-in
  Student-t heavy-tail distribution. SAC raises a classed error for that
  combination; use information reliability or explicitly define a finite,
  truncated empirical estimand.
* The 3PL release remains unidimensional. Multidimensional models are outside
  the 0.3.0 scope.

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
