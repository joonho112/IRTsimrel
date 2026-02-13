# Stochastic Approximation Calibration (Algorithm 2: SAC)

`sac_calibrate()` implements Algorithm 2 (Stochastic Approximation
Calibration, SAC) for reliability-targeted IRT simulation using the
Robbins-Monro stochastic approximation algorithm.

Given a target marginal reliability \\\rho^\*\\, a latent distribution
generator
[`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
(for \\G\\) and an item parameter generator
[`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
(for \\H\\), the function iteratively searches for a global
discrimination scale \\c^\* \> 0\\ such that the population reliability
\\\rho(c)\\ of the Rasch/2PL model is approximately equal to
\\\rho^\*\\.

SAC complements EQC (Algorithm 1) by:

1.  Providing an independent validation of EQC calibration results.

2.  Enabling calibration to the exact marginal reliability \\\bar{w}\\
    (not just the average-information approximation \\\tilde{\rho}\\).

3.  Handling complex data-generating processes where analytic
    information functions may be unavailable.

The algorithm uses the Robbins-Monro update rule: \$\$c\_{n+1} = c_n -
a_n \cdot (\hat{\rho}\_n - \rho^\*)\$\$

where \\a_n = a / (n + A)^\gamma\\ is a decreasing step size sequence
satisfying \\\sum a_n = \infty\\ and \\\sum a_n^2 \< \infty\\.

## Usage

``` r
sac_calibrate(
  target_rho,
  n_items,
  model = c("rasch", "2pl"),
  latent_shape = "normal",
  item_source = "parametric",
  latent_params = list(),
  item_params = list(),
  reliability_metric = c("msem", "info", "bar", "tilde"),
  c_init = NULL,
  M_per_iter = 500L,
  M_pre = 10000L,
  n_iter = 300L,
  burn_in = NULL,
  step_params = list(),
  c_bounds = c(0.01, 20),
  resample_items = TRUE,
  seed = NULL,
  verbose = FALSE
)

spc_calibrate(...)

# S3 method for class 'sac_result'
print(x, digits = 4, ...)

# S3 method for class 'sac_result'
summary(object, ...)
```

## Arguments

- target_rho:

  Numeric in (0, 1). Target marginal reliability \\\rho^\*\\.

- n_items:

  Integer. Number of items in the test form.

- model:

  Character. Measurement model: `"rasch"` or `"2pl"`. For `"rasch"`, all
  baseline discriminations are set to 1 before scaling.

- latent_shape:

  Character. Shape argument passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
  (e.g. `"normal"`, `"bimodal"`, `"heavy_tail"`, ...).

- item_source:

  Character. Source argument passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
  (e.g. `"parametric"`, `"irw"`, `"hierarchical"`, `"custom"`). Defaults
  to `"parametric"` since the irw package is an optional dependency
  (listed in Suggests). Use `"irw"` for empirically-grounded
  difficulties when the irw package is installed.

- latent_params:

  List. Additional arguments passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- item_params:

  List. Additional arguments passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- reliability_metric:

  Character. Reliability definition used inside SAC:

  `"msem"`

  :   MSEM-based marginal reliability (theoretically exact, targets
      \\\bar{w}\\).

  `"info"`

  :   Average-information reliability (faster, targets
      \\\tilde{\rho}\\).

  Synonyms: `"bar"` for `"msem"`, `"tilde"` for `"info"`.

- c_init:

  Numeric, `eqc_result` object, or NULL. Initial value for the scaling
  factor \\c_0\\.

  - If an `eqc_result` object is provided, its `c_star` is used (warm
    start).

  - If a numeric value is provided, it is used directly.

  - If NULL, initialized using Analytic Pre-Calibration (APC).

  Providing a warm start from EQC greatly accelerates convergence.

- M_per_iter:

  Integer. Number of Monte Carlo samples per iteration for estimating
  reliability. Default: 500. Larger values reduce variance but increase
  computation time.

- M_pre:

  Integer. Number of Monte Carlo samples for pre-calculating the latent
  variance \\\sigma^2\_\theta\\. Default: 10000. This variance is fixed
  throughout the iterations for stability. This is a CRITICAL parameter
  for numerical stability.

- n_iter:

  Integer. Total number of Robbins-Monro iterations. Default: 300.

- burn_in:

  Integer. Number of initial iterations to discard before Polyak-Ruppert
  averaging. Default: `floor(n_iter / 2)`.

- step_params:

  List. Parameters controlling the step size sequence:

  `a`

  :   Base step size (default: 1.0)

  `A`

  :   Stabilization constant (default: 50)

  `gamma`

  :   Decay exponent (default: 0.67, i.e., 2/3)

- c_bounds:

  Numeric length-2 vector. Projection bounds for \\c\\. Iterates are
  clipped to this interval after each update. Default: c(0.01, 20).

- resample_items:

  Logical. If TRUE (default), resample item parameters at each
  iteration. If FALSE, fix item parameters across all iterations
  (reduces variance but may introduce bias).

- seed:

  Optional integer for reproducibility.

- verbose:

  Logical or integer. If TRUE or \>= 1, print progress messages. If \>=
  2, print detailed iteration-level output.

- ...:

  Arguments passed to `sac_calibrate()`.

- x:

  An object of class `"sac_result"`.

- digits:

  Integer. Number of decimal places for printing.

- object:

  An object of class `"sac_result"`.

## Value

An object of class `"sac_result"` (a list) with elements:

- `c_star`:

  Calibrated discrimination scale (Polyak-Ruppert average).

- `c_final`:

  Final iterate \\c\_{n\_{iter}}\\.

- `target_rho`:

  Target reliability \\\rho^\*\\.

- `achieved_rho`:

  Estimated reliability at \\c^\*\\ (post-calibration).

- `theta_var`:

  Pre-calculated latent variance used throughout.

- `trajectory`:

  Numeric vector of all iterates.

- `rho_trajectory`:

  Numeric vector of reliability estimates.

- `init_method`:

  Character indicating initialization method.

- `metric`:

  Reliability metric used.

- `convergence`:

  List with convergence diagnostics.

- `beta_vec`:

  Item difficulties from the final post-calibration draw.

- `lambda_base`:

  Baseline (unscaled) item discriminations.

- `lambda_scaled`:

  Scaled item discriminations (`lambda_base * c_star`).

- `items_base`:

  Baseline item_params object (scale = 1).

- `items_calib`:

  Calibrated item_params object (discriminations scaled by c_star).

- `theta_quad`:

  Theta sample used for post-calibration reliability estimate.

The input object, invisibly.

An object of class `"summary.sac_result"` containing key calibration
results.

## References

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics, 22*(3), 400–407.

Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
approximation by averaging. *SIAM Journal on Control and Optimization,
30*(4), 838–855.

## See also

[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
for the faster deterministic Algorithm 1,
[`compute_rho_bar`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
and
[`compute_rho_tilde`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
for reliability computation utilities,
[`compare_eqc_sac`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
for comparing EQC and SAC results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Basic SAC calibration
sac_result <- sac_calibrate(
  target_rho = 0.75,
  n_items = 20,
  model = "rasch",
  n_iter = 200,
  seed = 12345,
  verbose = TRUE
)
print(sac_result)
plot(sac_result)

# Example 2: Warm start from EQC (RECOMMENDED)
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "2pl",
  M = 10000,
  seed = 42
)

sac_result <- sac_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "2pl",
  c_init = eqc_result,  # Direct EQC object passing!
  n_iter = 100,
  seed = 42
)

# Compare EQC and SAC results
cat(sprintf("EQC c* = %.4f, SAC c* = %.4f\n",
            eqc_result$c_star, sac_result$c_star))
} # }
```
