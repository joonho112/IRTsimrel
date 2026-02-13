# Complete API Reference with Examples

## Overview

This vignette provides a comprehensive reference for every exported
function in **IRTsimrel** v0.2.0. Functions are organized into five
categories:

1.  **Calibration Functions** – The core algorithms for
    reliability-targeted simulation (`eqc_calibrate`, `sac_calibrate`,
    `compare_eqc_sac`).
2.  **Simulation Functions** – Generators for latent abilities, item
    parameters, and response data (`sim_latentG`, `sim_item_params`,
    `simulate_response_data`, `compare_shapes`).
3.  **Reliability Functions** – Low-level utilities for computing and
    exploring reliability (`compute_rho_bar`, `compute_rho_tilde`,
    `compute_rho_both`, `compute_apc_init`, `check_feasibility`,
    `rho_curve`, `compute_reliability_tam`).
4.  **S3 Methods** – Print, summary, plot, coef, predict, and coercion
    methods for all object classes.
5.  **Deprecated Functions** – Legacy aliases retained for backward
    compatibility (`spc_calibrate`, `compare_eqc_spc`).

For conceptual introductions, see
[`vignette("introduction")`](https://joonho112.github.io/IRTsimrel/articles/introduction.md).
For applied workflow guidance, see
[`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md).
For the mathematical theory, see
[`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md).

**Citation:** Lee, J. (2025). Reliability-targeted simulation of item
response data: Solving the inverse design problem. *arXiv preprint
arXiv:2512.16012*.

------------------------------------------------------------------------

## 1. Calibration Functions

### 1.1 `eqc_calibrate()` – Empirical Quadrature Calibration

Implements Algorithm 1 (EQC/SQC) from Lee (2025). Given a target
marginal reliability, EQC draws a large fixed quadrature sample of
abilities and item parameters, then solves the scalar root-finding
problem via Brent’s method
([`uniroot()`](https://rdrr.io/r/stats/uniroot.html)). This is the
recommended starting point for most simulation studies.

#### Signature

``` r
eqc_calibrate(
  target_rho,
  n_items,
  model             = c("rasch", "2pl"),
  latent_shape      = "normal",
  item_source       = "parametric",
  latent_params     = list(),
  item_params       = list(),
  reliability_metric = c("info", "tilde", "msem", "bar"),
  M                 = 10000L,
  c_bounds          = c(0.3, 3),
  tol               = 1e-4,
  seed              = NULL,
  verbose           = FALSE
)
```

#### Parameters

| Parameter            | Type            | Default                                      | Description                                                                                                                                 |
|:---------------------|:----------------|:---------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------|
| `target_rho`         | numeric         | *required*                                   | Target marginal reliability in (0, 1).                                                                                                      |
| `n_items`            | integer         | *required*                                   | Number of test items.                                                                                                                       |
| `model`              | character       | `"rasch"`                                    | Measurement model: `"rasch"` or `"2pl"`.                                                                                                    |
| `latent_shape`       | character       | `"normal"`                                   | Shape passed to [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).                                          |
| `item_source`        | character       | `"parametric"`                               | Source passed to [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md) (e.g., `"parametric"`, `"irw"`). |
| `latent_params`      | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).                                 |
| `item_params`        | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).                         |
| `reliability_metric` | character       | `"info"`                                     | Metric: `"info"`/`"tilde"` (recommended) or `"msem"`/`"bar"`.                                                                               |
| `M`                  | integer         | `10000L`                                     | Quadrature sample size.                                                                                                                     |
| `c_bounds`           | numeric(2)      | `c(0.3, 3)`                                  | Search bounds for scaling factor c.                                                                                                         |
| `tol`                | numeric         | `1e-4`                                       | Tolerance for [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).                                                                          |
| `seed`               | integer or NULL | `NULL`                                       | Random seed for reproducibility.                                                                                                            |
| `verbose`            | logical         | `FALSE`                                      | Print progress messages.                                                                                                                    |

#### Return Value

An object of class `"eqc_result"` (a list) with elements:

| Element         | Description                                       |
|:----------------|:--------------------------------------------------|
| `c_star`        | Calibrated discrimination scale factor.           |
| `target_rho`    | The target reliability supplied.                  |
| `achieved_rho`  | Empirical quadrature estimate at `c_star`.        |
| `metric`        | Internal metric label (`"info"` or `"msem"`).     |
| `model`         | Model used (`"rasch"` or `"2pl"`).                |
| `n_items`       | Number of items.                                  |
| `M`             | Quadrature sample size.                           |
| `theta_quad`    | Length-M vector of quadrature abilities.          |
| `theta_var`     | Sample variance of `theta_quad`.                  |
| `beta_vec`      | Item difficulties.                                |
| `lambda_base`   | Baseline (unscaled) discriminations.              |
| `lambda_scaled` | Discriminations scaled by `c_star`.               |
| `items_base`    | `item_params` object at scale = 1.                |
| `items_calib`   | `item_params` object with scaled discriminations. |
| `call`          | The matched call.                                 |
| `misc`          | List of convergence diagnostics.                  |

#### Examples

``` r
# Rasch model with 25 items targeting rho = 0.80
eqc_res <- eqc_calibrate(
  target_rho = 0.80,
  n_items    = 25,
  model      = "rasch",
  M          = 5000L,
  seed       = 42
)
eqc_res
#> 
#> =======================================================
#>   Empirical Quadrature Calibration (EQC) Results
#> =======================================================
#> 
#> Calibration Summary:
#>   Model                        : RASCH
#>   Target reliability (rho*)    : 0.8000
#>   Achieved reliability         : 0.8000
#>   Absolute error               : 1.19e-07
#>   Scaling factor (c*)          : 0.8995
#> 
#> Design Parameters:
#>   Number of items (I)          : 25
#>   Quadrature points (M)        : 5000
#>   Reliability metric           : Average-information (tilde)
#>   Latent variance              : 1.0099
#> 
#> Convergence:
#>   Root status                  : uniroot_success
#>   Search bracket               : [0.300, 3.000]
#>   Bracket reliabilities        : [0.3539, 0.9550]
#> 
#> Parameter Summaries:
#>   theta:        mean = -0.014, sd = 1.005
#>   beta:         mean = 0.000, sd = 0.861, range = [-2.17, 1.45]
#>   lambda_base:  mean = 1.000, sd = 0.000
#>   lambda_scaled: mean = 0.899, sd = 0.000
```

``` r
# 2PL model with bimodal latent distribution
eqc_2pl <- eqc_calibrate(
  target_rho   = 0.85,
  n_items      = 30,
  model        = "2pl",
  latent_shape = "bimodal",
  M            = 5000L,
  seed         = 42
)
cat(sprintf("c* = %.4f, achieved rho = %.4f\n",
            eqc_2pl$c_star, eqc_2pl$achieved_rho))
#> c* = 1.0072, achieved rho = 0.8500
```

See
[`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)
for a detailed walk-through.

------------------------------------------------------------------------

### 1.2 `sac_calibrate()` – Stochastic Approximation Calibration

Implements Algorithm 2 (SAC) from Lee (2025). Uses the Robbins–Monro
stochastic approximation framework with Polyak–Ruppert averaging. SAC
provides an independent validation of EQC and can target the exact
MSEM-based marginal reliability.

#### Signature

``` r
sac_calibrate(
  target_rho,
  n_items,
  model              = c("rasch", "2pl"),
  latent_shape       = "normal",
  item_source        = "parametric",
  latent_params      = list(),
  item_params        = list(),
  reliability_metric = c("msem", "info", "bar", "tilde"),
  c_init             = NULL,
  M_per_iter         = 500L,
  M_pre              = 10000L,
  n_iter             = 300L,
  burn_in            = NULL,
  step_params        = list(),
  c_bounds           = c(0.01, 20),
  resample_items     = TRUE,
  seed               = NULL,
  verbose            = FALSE
)
```

#### Parameters

| Parameter            | Type                           | Default                                      | Description                                                                                                         |
|:---------------------|:-------------------------------|:---------------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| `target_rho`         | numeric                        | *required*                                   | Target marginal reliability in (0, 1).                                                                              |
| `n_items`            | integer                        | *required*                                   | Number of test items.                                                                                               |
| `model`              | character                      | `"rasch"`                                    | Measurement model: `"rasch"` or `"2pl"`.                                                                            |
| `latent_shape`       | character                      | `"normal"`                                   | Shape passed to [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).                  |
| `item_source`        | character                      | `"parametric"`                               | Source passed to [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).         |
| `latent_params`      | list                           | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).         |
| `item_params`        | list                           | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md). |
| `reliability_metric` | character                      | `"msem"`                                     | Metric: `"msem"`/`"bar"` or `"info"`/`"tilde"`.                                                                     |
| `c_init`             | numeric, `eqc_result`, or NULL | `NULL`                                       | Initial scaling factor. If `NULL`, uses APC. If an `eqc_result` object, uses its `c_star` (warm start).             |
| `M_per_iter`         | integer                        | `500L`                                       | Monte Carlo samples per iteration.                                                                                  |
| `M_pre`              | integer                        | `10000L`                                     | Samples for pre-calculating latent variance.                                                                        |
| `n_iter`             | integer                        | `300L`                                       | Total Robbins–Monro iterations.                                                                                     |
| `burn_in`            | integer or NULL                | `NULL`                                       | Iterations to discard before averaging. Default: `floor(n_iter / 2)`.                                               |
| `step_params`        | list                           | [`list()`](https://rdrr.io/r/base/list.html) | Step size parameters: `a` (base, default 1.0), `A` (stabilization, default 50), `gamma` (decay, default 0.67).      |
| `c_bounds`           | numeric(2)                     | `c(0.01, 20)`                                | Projection bounds for iterates.                                                                                     |
| `resample_items`     | logical                        | `TRUE`                                       | Re-draw item parameters each iteration.                                                                             |
| `seed`               | integer or NULL                | `NULL`                                       | Random seed.                                                                                                        |
| `verbose`            | logical or integer             | `FALSE`                                      | Print progress (2 for iteration-level detail).                                                                      |

#### Return Value

An object of class `"sac_result"` (a list) with elements:

| Element          | Description                                 |
|:-----------------|:--------------------------------------------|
| `c_star`         | Polyak–Ruppert averaged scaling factor.     |
| `c_final`        | Final iterate value.                        |
| `target_rho`     | Target reliability.                         |
| `achieved_rho`   | Post-calibration reliability estimate.      |
| `theta_var`      | Pre-calculated latent variance.             |
| `trajectory`     | Numeric vector of all iterates.             |
| `rho_trajectory` | Reliability estimates at each iteration.    |
| `metric`         | Internal metric label.                      |
| `model`          | Model used.                                 |
| `n_items`        | Number of items.                            |
| `n_iter`         | Total iterations run.                       |
| `burn_in`        | Burn-in used.                               |
| `M_per_iter`     | Samples per iteration.                      |
| `M_pre`          | Pre-calculation sample size.                |
| `step_params`    | Step size parameters used.                  |
| `c_bounds`       | Projection bounds.                          |
| `c_init`         | Initial value used.                         |
| `init_method`    | Initialization method label.                |
| `convergence`    | List of convergence diagnostics.            |
| `beta_vec`       | Item difficulties (post-calibration draw).  |
| `lambda_base`    | Baseline discriminations.                   |
| `lambda_scaled`  | Scaled discriminations.                     |
| `items_base`     | Baseline `item_params` object.              |
| `items_calib`    | Calibrated `item_params` object.            |
| `theta_quad`     | Theta sample for post-calibration estimate. |
| `call`           | The matched call.                           |

#### Examples

``` r
# SAC with APC initialization (default)
sac_res <- sac_calibrate(
  target_rho = 0.80,
  n_items    = 25,
  model      = "rasch",
  n_iter     = 200L,
  M_per_iter = 500L,
  M_pre      = 5000L,
  seed       = 42,
  verbose    = FALSE
)
cat(sprintf("SAC c* = %.4f, achieved rho = %.4f\n",
            sac_res$c_star, sac_res$achieved_rho))
#> SAC c* = 1.0504, achieved rho = 0.8292
```

``` r
# Warm start from EQC result (recommended workflow)
sac_warm <- sac_calibrate(
  target_rho = 0.80,
  n_items    = 25,
  model      = "rasch",
  c_init     = eqc_res,
  n_iter     = 100L,
  M_per_iter = 500L,
  M_pre      = 5000L,
  seed       = 42,
  verbose    = FALSE
)
cat(sprintf("EQC c* = %.4f, SAC c* = %.4f\n",
            eqc_res$c_star, sac_warm$c_star))
#> EQC c* = 0.8995, SAC c* = 0.9232
```

See
[`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md)
for convergence tuning details.

------------------------------------------------------------------------

### 1.3 `compare_eqc_sac()` – Compare Calibration Results

Computes agreement diagnostics between EQC and SAC calibration results
targeting the same reliability.

#### Signature

``` r
compare_eqc_sac(eqc_result, sac_result, verbose = TRUE)
```

#### Parameters

| Parameter    | Type         | Default    | Description                                                                                        |
|:-------------|:-------------|:-----------|:---------------------------------------------------------------------------------------------------|
| `eqc_result` | `eqc_result` | *required* | Output from [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md). |
| `sac_result` | `sac_result` | *required* | Output from [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md). |
| `verbose`    | logical      | `TRUE`     | Print comparison summary.                                                                          |

#### Return Value

A list (returned invisibly) with components:

| Element      | Description                                      |
|:-------------|:-------------------------------------------------|
| `c_eqc`      | Calibrated `c*` from EQC.                        |
| `c_sac`      | Calibrated `c*` from SAC.                        |
| `diff_abs`   | Absolute difference between the two `c*` values. |
| `diff_pct`   | Percent difference relative to EQC.              |
| `agreement`  | Logical. `TRUE` if percent difference \< 5%.     |
| `target_rho` | Shared target reliability.                       |

#### Examples

``` r
# Compare results from the two algorithms
comp <- compare_eqc_sac(eqc_res, sac_warm, verbose = TRUE)
#> Warning in compare_eqc_sac(eqc_res, sac_warm, verbose = TRUE): Reliability
#> metric differs between EQC ('info') and SAC ('msem').
#> 
#> =======================================================
#>   EQC vs SAC Comparison
#> =======================================================
#> 
#>   Target reliability  : 0.8000
#>   EQC c*              : 0.899499
#>   SAC c*              : 0.923198
#>   Absolute difference : 0.023698
#>   Percent difference  : 2.63%
#>   Agreement (< 5%)    : YES
#> 
cat(sprintf("Agreement: %s (%.2f%% difference)\n",
            ifelse(comp$agreement, "YES", "NO"), comp$diff_pct))
#> Agreement: YES (2.63% difference)
```

------------------------------------------------------------------------

## 2. Simulation Functions

### 2.1 `sim_latentG()` – Simulate Latent Ability Distributions

Generates latent abilities from a flexible family of pre-standardized
distributions. Each built-in shape is mathematically constructed to have
mean 0 and variance 1, so changes in distributional shape do not alter
the scale.

#### Signature

``` r
sim_latentG(
  n,
  shape              = c("normal", "bimodal", "trimodal", "multimodal",
                          "skew_pos", "skew_neg", "heavy_tail",
                          "light_tail", "uniform", "floor", "ceiling",
                          "custom"),
  sigma              = 1,
  mu                 = 0,
  xcov               = NULL,
  beta               = NULL,
  shape_params       = list(),
  mixture_spec       = NULL,
  standardize_custom = TRUE,
  seed               = NULL,
  return_z           = TRUE
)
```

#### Parameters

| Parameter            | Type            | Default                                      | Description                                        |
|:---------------------|:----------------|:---------------------------------------------|:---------------------------------------------------|
| `n`                  | integer         | *required*                                   | Number of persons.                                 |
| `shape`              | character       | `"normal"`                                   | Distributional shape (see table below).            |
| `sigma`              | numeric         | `1`                                          | Standard deviation of the residual latent trait.   |
| `mu`                 | numeric         | `0`                                          | Grand mean of latent abilities.                    |
| `xcov`               | matrix or NULL  | `NULL`                                       | Optional covariate matrix (`n` rows).              |
| `beta`               | numeric or NULL | `NULL`                                       | Regression coefficients for `xcov`.                |
| `shape_params`       | list            | [`list()`](https://rdrr.io/r/base/list.html) | Shape-specific parameters (see below).             |
| `mixture_spec`       | list or NULL    | `NULL`                                       | For `shape = "custom"`: `weights`, `means`, `sds`. |
| `standardize_custom` | logical         | `TRUE`                                       | Standardize custom mixtures to mean 0, variance 1. |
| `seed`               | integer or NULL | `NULL`                                       | Random seed.                                       |
| `return_z`           | logical         | `TRUE`                                       | Include standardized draws in output.              |

**Available shapes:**

| Shape          | Description                           | Key Parameters                                                    |
|:---------------|:--------------------------------------|:------------------------------------------------------------------|
| `"normal"`     | Standard normal                       | –                                                                 |
| `"bimodal"`    | Symmetric two-component mixture       | `delta` (mode separation, default 0.8)                            |
| `"trimodal"`   | Symmetric three-component mixture     | `w0` (central weight, default 1/3), `m` (side means, default 1.2) |
| `"multimodal"` | Symmetric four-component mixture      | `m1`, `m2`, `w_inner`                                             |
| `"skew_pos"`   | Right-skewed (standardized Gamma)     | `k` (shape, default 4)                                            |
| `"skew_neg"`   | Left-skewed (negated Gamma)           | `k` (shape, default 4)                                            |
| `"heavy_tail"` | Heavy-tailed (standardized Student-t) | `df` (degrees of freedom, default 5)                              |
| `"light_tail"` | Light-tailed (platykurtic mixture)    | –                                                                 |
| `"uniform"`    | Uniform on \[-sqrt(3), sqrt(3)\]      | –                                                                 |
| `"floor"`      | Floor effect                          | `w_floor` (default 0.3), `m_floor` (default -1.5)                 |
| `"ceiling"`    | Ceiling effect                        | `w_ceil` (default 0.3), `m_ceil` (default 1.5)                    |
| `"custom"`     | User-specified mixture                | See `mixture_spec`                                                |

#### Return Value

An object of class `"latent_G"` (a list) with elements:

| Element          | Description                                      |
|:-----------------|:-------------------------------------------------|
| `theta`          | Numeric vector of simulated latent abilities.    |
| `z`              | Standardized draws (if `return_z = TRUE`).       |
| `eta_cov`        | Covariate linear predictor (0 if no covariates). |
| `mu`             | Grand mean.                                      |
| `sigma`          | Scale parameter.                                 |
| `shape`          | Shape label.                                     |
| `shape_params`   | Shape parameters used.                           |
| `n`              | Sample size.                                     |
| `sample_moments` | List: `mean`, `sd`, `skewness`, `kurtosis`.      |

#### Examples

``` r
# Standard normal abilities
g_norm <- sim_latentG(n = 2000, shape = "normal", seed = 42)
cat(sprintf("Mean = %.3f, SD = %.3f, Skew = %.3f, Kurt = %.3f\n",
            g_norm$sample_moments$mean,
            g_norm$sample_moments$sd,
            g_norm$sample_moments$skewness,
            g_norm$sample_moments$kurtosis))
#> Mean = -0.016, SD = 0.994, Skew = 0.013, Kurt = 0.056
```

``` r
# Bimodal distribution with strong separation
g_bimod <- sim_latentG(
  n            = 2000,
  shape        = "bimodal",
  shape_params = list(delta = 0.9),
  seed         = 42
)
cat(sprintf("Bimodal: Mean = %.3f, SD = %.3f, Kurt = %.3f\n",
            g_bimod$sample_moments$mean,
            g_bimod$sample_moments$sd,
            g_bimod$sample_moments$kurtosis))
#> Bimodal: Mean = -0.017, SD = 0.989, Kurt = -1.298
```

``` r
# Positively skewed distribution
g_skew <- sim_latentG(
  n            = 2000,
  shape        = "skew_pos",
  shape_params = list(k = 3),
  seed         = 42
)
cat(sprintf("Skew_pos: Skewness = %.3f\n",
            g_skew$sample_moments$skewness))
#> Skew_pos: Skewness = 1.152
```

``` r
# Custom three-component mixture
g_custom <- sim_latentG(
  n            = 2000,
  shape        = "custom",
  mixture_spec = list(
    weights = c(0.3, 0.5, 0.2),
    means   = c(-1.5, 0, 2),
    sds     = c(0.5, 0.7, 0.5)
  ),
  seed = 42
)
cat(sprintf("Custom: Mean = %.3f, SD = %.3f\n",
            g_custom$sample_moments$mean,
            g_custom$sample_moments$sd))
#> Custom: Mean = 0.011, SD = 0.996
```

See
[`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md)
for detailed shape comparisons.

------------------------------------------------------------------------

### 2.2 `sim_item_params()` – Simulate Item Parameters

Generates item difficulty and discrimination parameters for IRT models.
Supports multiple sources (parametric, IRW, hierarchical, custom) and
methods for inducing difficulty–discrimination correlation.

#### Signature

``` r
sim_item_params(
  n_items,
  model                 = c("rasch", "2pl"),
  source                = c("irw", "parametric", "hierarchical", "custom"),
  method                = c("copula", "conditional", "independent"),
  n_forms               = 1L,
  difficulty_params     = list(),
  discrimination_params = list(),
  hierarchical_params   = list(),
  custom_params         = list(),
  scale                 = 1,
  center_difficulties   = TRUE,
  seed                  = NULL
)
```

#### Parameters

| Parameter               | Type            | Default                                      | Description                                                                                 |
|:------------------------|:----------------|:---------------------------------------------|:--------------------------------------------------------------------------------------------|
| `n_items`               | integer         | *required*                                   | Number of items per form.                                                                   |
| `model`                 | character       | `"rasch"`                                    | `"rasch"` or `"2pl"`.                                                                       |
| `source`                | character       | `"irw"`                                      | `"irw"`, `"parametric"`, `"hierarchical"`, or `"custom"`.                                   |
| `method`                | character       | `"copula"`                                   | For 2PL: `"copula"` (recommended), `"conditional"`, `"independent"`.                        |
| `n_forms`               | integer         | `1L`                                         | Number of parallel test forms.                                                              |
| `difficulty_params`     | list            | [`list()`](https://rdrr.io/r/base/list.html) | For parametric: `mu` (default 0), `sigma` (default 1), `distribution` (default `"normal"`). |
| `discrimination_params` | list            | [`list()`](https://rdrr.io/r/base/list.html) | For 2PL: `mu_log` (default 0), `sigma_log` (default 0.3), `rho` (default -0.3).             |
| `hierarchical_params`   | list            | [`list()`](https://rdrr.io/r/base/list.html) | For hierarchical: `mu` (2-vector), `tau` (2-vector), `rho`.                                 |
| `custom_params`         | list            | [`list()`](https://rdrr.io/r/base/list.html) | For custom: `beta` (vector or function), `lambda` (vector or function).                     |
| `scale`                 | numeric         | `1`                                          | Global discrimination scaling factor.                                                       |
| `center_difficulties`   | logical         | `TRUE`                                       | Center difficulties to sum to zero.                                                         |
| `seed`                  | integer or NULL | `NULL`                                       | Random seed.                                                                                |

#### Return Value

An object of class `"item_params"` (a list) with elements:

| Element    | Description                                                            |
|:-----------|:-----------------------------------------------------------------------|
| `data`     | Data frame: `form_id`, `item_id`, `beta`, `lambda`, `lambda_unscaled`. |
| `model`    | Model type.                                                            |
| `source`   | Source used.                                                           |
| `method`   | Method for discriminations (NA for Rasch or hierarchical).             |
| `n_items`  | Items per form.                                                        |
| `n_forms`  | Forms generated.                                                       |
| `scale`    | Scale factor applied.                                                  |
| `centered` | Whether difficulties were centered.                                    |
| `params`   | Parameters used for generation.                                        |
| `achieved` | Achieved statistics (correlations, moments).                           |

#### Examples

``` r
# Rasch items with parametric difficulties
items_rasch <- sim_item_params(
  n_items = 25,
  model   = "rasch",
  source  = "parametric",
  seed    = 42
)
items_rasch
#> Item Parameters Object
#> ======================
#>   Model          : RASCH
#>   Source         : parametric
#>   Items per form : 25
#>   Number of forms: 1
#>   Scale factor   : 1.000
#>   Centered       : Yes
#> 
#> Difficulty (beta):
#>   Mean: 0.0000, SD: 1.3064, Range: [-2.844, 2.099]
```

``` r
# 2PL items with copula-induced correlation
items_2pl <- sim_item_params(
  n_items               = 30,
  model                 = "2pl",
  source                = "parametric",
  method                = "copula",
  discrimination_params = list(rho = -0.3, mu_log = 0, sigma_log = 0.3),
  seed                  = 42
)
cat(sprintf("Achieved Spearman r(beta, log-lambda) = %.3f\n",
            items_2pl$achieved$overall$cor_spearman_pooled))
#> Achieved Spearman r(beta, log-lambda) = -0.378
```

``` r
# Hierarchical 2PL (Glas & van der Linden style)
items_hier <- sim_item_params(
  n_items              = 25,
  model                = "2pl",
  source               = "hierarchical",
  hierarchical_params  = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3),
  seed                 = 42
)
head(as.data.frame(items_hier))
#>   form_id item_id       beta    lambda lambda_unscaled
#> 1       1       1  1.1864454 0.9930935       0.9930935
#> 2       1       2 -0.7521311 1.1116693       1.1116693
#> 3       1       3  0.2039055 1.4755926       1.4755926
#> 4       1       4  0.4317000 0.8526922       0.8526922
#> 5       1       5  0.2238748 1.1270078       1.1270078
#> 6       1       6 -0.3070683 0.9052717       0.9052717
```

``` r
# Custom user-supplied parameters
items_custom <- sim_item_params(
  n_items       = 10,
  model         = "2pl",
  source        = "custom",
  custom_params = list(
    beta   = seq(-2, 2, length.out = 10),
    lambda = rep(1.2, 10)
  ),
  seed = 42
)
#> Warning in cor(df$beta, log(df$lambda_unscaled)): the standard deviation is
#> zero
#> Warning in cor(df$beta, log(df$lambda_unscaled), method = "spearman"): the
#> standard deviation is zero
#> Warning in cor(data$beta, log(data$lambda_unscaled)): the standard deviation is
#> zero
#> Warning in cor(data$beta, log(data$lambda_unscaled), method = "spearman"): the
#> standard deviation is zero
head(as.data.frame(items_custom))
#>   form_id item_id       beta lambda lambda_unscaled
#> 1       1       1 -2.0000000    1.2             1.2
#> 2       1       2 -1.5555556    1.2             1.2
#> 3       1       3 -1.1111111    1.2             1.2
#> 4       1       4 -0.6666667    1.2             1.2
#> 5       1       5 -0.2222222    1.2             1.2
#> 6       1       6  0.2222222    1.2             1.2
```

See
[`vignette("item-parameters")`](https://joonho112.github.io/IRTsimrel/articles/item-parameters.md)
for a complete discussion of sources and methods.

------------------------------------------------------------------------

### 2.3 `simulate_response_data()` – Generate Response Matrices

Simulates binary (0/1) item response data using calibrated parameters
from
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
or
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).
This is the final step before external validation (e.g., with TAM).

#### Signature

``` r
simulate_response_data(
  result,
  n_persons,
  latent_shape  = "normal",
  latent_params = list(),
  seed          = NULL
)
```

#### Parameters

| Parameter       | Type                         | Default                                      | Description                                                                                                 |
|:----------------|:-----------------------------|:---------------------------------------------|:------------------------------------------------------------------------------------------------------------|
| `result`        | `eqc_result` or `sac_result` | *required*                                   | Calibration result object.                                                                                  |
| `n_persons`     | integer                      | *required*                                   | Number of simulated examinees.                                                                              |
| `latent_shape`  | character                    | `"normal"`                                   | Shape for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).                |
| `latent_params` | list                         | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md). |
| `seed`          | integer or NULL              | `NULL`                                       | Random seed.                                                                                                |

#### Return Value

A list with components:

| Element           | Description                                                           |
|:------------------|:----------------------------------------------------------------------|
| `response_matrix` | N x I matrix of binary responses (column names: `item1`, …, `itemI`). |
| `theta`           | True latent abilities (length N).                                     |
| `beta`            | Item difficulties (length I).                                         |
| `lambda`          | Scaled item discriminations (length I).                               |

#### Examples

``` r
# Generate 500 examinees from EQC calibration
sim_data <- simulate_response_data(
  result    = eqc_res,
  n_persons = 500,
  seed      = 123
)
cat(sprintf("Response matrix: %d persons x %d items\n",
            nrow(sim_data$response_matrix),
            ncol(sim_data$response_matrix)))
#> Response matrix: 500 persons x 25 items
cat(sprintf("Mean proportion correct: %.3f\n",
            mean(sim_data$response_matrix)))
#> Mean proportion correct: 0.502
```

``` r
# Generate from SAC result with skewed abilities
sim_data2 <- simulate_response_data(
  result       = sac_warm,
  n_persons    = 500,
  latent_shape = "skew_pos",
  seed         = 123
)
cat(sprintf("Theta skewness: %.3f\n",
            mean(((sim_data2$theta - mean(sim_data2$theta)) /
                    sd(sim_data2$theta))^3)))
#> Theta skewness: 0.906
```

------------------------------------------------------------------------

### 2.4 `compare_shapes()` – Compare Distribution Shapes

Generates and compares multiple latent distributions side-by-side using
faceted density plots. Requires the **ggplot2** package.

#### Signature

``` r
compare_shapes(
  n      = 2000,
  shapes = c("normal", "bimodal", "trimodal",
             "skew_pos", "skew_neg", "heavy_tail",
             "uniform"),
  sigma  = 1,
  seed   = NULL
)
```

#### Parameters

| Parameter | Type             | Default     | Description                   |
|:----------|:-----------------|:------------|:------------------------------|
| `n`       | integer          | `2000`      | Sample size per distribution. |
| `shapes`  | character vector | (see above) | Shapes to compare.            |
| `sigma`   | numeric          | `1`         | Common scale parameter.       |
| `seed`    | integer or NULL  | `NULL`      | Random seed.                  |

#### Return Value

A `ggplot` object with faceted density plots, one panel per shape. A
dashed red line shows the N(0, 1) reference density.

#### Example

``` r
p <- compare_shapes(
  n      = 2000,
  shapes = c("normal", "bimodal", "skew_pos", "heavy_tail"),
  seed   = 42
)
print(p)
```

![Comparison of latent distribution
shapes.](api-reference_files/figure-html/compare-shapes-1.png)

Comparison of latent distribution shapes.

See
[`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md)
for extended shape comparisons.

------------------------------------------------------------------------

## 3. Reliability Functions

### 3.1 `compute_rho_bar()` – MSEM-Based Marginal Reliability

Computes the MSEM-based marginal reliability using the harmonic mean of
test information.

#### Signature

``` r
compute_rho_bar(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)
```

#### Parameters

| Parameter     | Type            | Default    | Description                                                          |
|:--------------|:----------------|:-----------|:---------------------------------------------------------------------|
| `c`           | numeric         | *required* | Global discrimination scaling factor.                                |
| `theta_vec`   | numeric vector  | *required* | Abilities.                                                           |
| `beta_vec`    | numeric vector  | *required* | Item difficulties.                                                   |
| `lambda_base` | numeric vector  | *required* | Baseline discriminations (before scaling).                           |
| `theta_var`   | numeric or NULL | `NULL`     | Pre-calculated theta variance. If `NULL`, computed from `theta_vec`. |

#### Return Value

A numeric scalar: the MSEM-based reliability.

#### Example

``` r
set.seed(1)
theta <- rnorm(2000)
beta  <- rnorm(20)
lambda0 <- rep(1, 20)

rho_bar_val <- compute_rho_bar(c = 1.0, theta, beta, lambda0)
cat(sprintf("rho_bar at c=1: %.4f\n", rho_bar_val))
#> rho_bar at c=1: 0.7856
```

------------------------------------------------------------------------

### 3.2 `compute_rho_tilde()` – Average-Information Reliability

Computes the average-information reliability using the arithmetic mean
of test information. Guaranteed monotone in c, making it the recommended
metric for EQC.

#### Signature

``` r
compute_rho_tilde(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)
```

Parameters and return value are identical to
[`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
above; only the internal formula differs.

#### Example

``` r
rho_tilde_val <- compute_rho_tilde(c = 1.0, theta, beta, lambda0)
cat(sprintf("rho_tilde at c=1: %.4f\n", rho_tilde_val))
#> rho_tilde at c=1: 0.7931

# Jensen's inequality: rho_tilde >= rho_bar
cat(sprintf("rho_tilde >= rho_bar: %s\n",
            ifelse(rho_tilde_val >= rho_bar_val, "TRUE", "FALSE")))
#> rho_tilde >= rho_bar: TRUE
```

------------------------------------------------------------------------

### 3.3 `compute_rho_both()` – Both Metrics in a Single Pass

Computes both reliability metrics from a single set of matrix
computations, avoiding redundant M x I operations.

#### Signature

``` r
compute_rho_both(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)
```

Parameters are identical to
[`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md).
The return value is a named list.

#### Return Value

| Element     | Description                      |
|:------------|:---------------------------------|
| `rho_tilde` | Average-information reliability. |
| `rho_bar`   | MSEM-based reliability.          |

#### Example

``` r
both <- compute_rho_both(c = 1.0, theta, beta, lambda0)
cat(sprintf("rho_tilde = %.4f, rho_bar = %.4f, gap = %.4f\n",
            both$rho_tilde, both$rho_bar,
            both$rho_tilde - both$rho_bar))
#> rho_tilde = 0.7931, rho_bar = 0.7856, gap = 0.0075
```

------------------------------------------------------------------------

### 3.4 `compute_apc_init()` – Analytic Pre-Calibration

Computes an initial scaling factor using a closed-form approximation
under Gaussian Rasch assumptions. Used internally by
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
when `c_init = NULL`.

#### Signature

``` r
compute_apc_init(target_rho, n_items, sigma_beta = 1.0)
```

#### Parameters

| Parameter    | Type    | Default    | Description              |
|:-------------|:--------|:-----------|:-------------------------|
| `target_rho` | numeric | *required* | Target reliability.      |
| `n_items`    | integer | *required* | Number of items.         |
| `sigma_beta` | numeric | `1.0`      | SD of item difficulties. |

#### Return Value

A numeric scalar: the initial scaling factor (bounded to \[0.1, 10\]).

#### Example

``` r
# Compare APC estimates across test lengths
for (I in c(10, 20, 30, 50)) {
  c0 <- compute_apc_init(target_rho = 0.80, n_items = I)
  cat(sprintf("  I = %2d: c_init = %.4f\n", I, c0))
}
#>   I = 10: c_init = 2.0988
#>   I = 20: c_init = 1.4841
#>   I = 30: c_init = 1.2117
#>   I = 50: c_init = 0.9386
```

------------------------------------------------------------------------

### 3.5 `check_feasibility()` – Feasibility Screening

Screens whether a target reliability is achievable for a given test
design by computing the range of reliabilities across a range of scaling
factors. Run this before calibration to avoid wasting time on infeasible
targets.

#### Signature

``` r
check_feasibility(
  n_items,
  model         = c("rasch", "2pl"),
  latent_shape  = "normal",
  item_source   = "parametric",
  c_bounds      = c(0.1, 10),
  M             = 10000L,
  seed          = NULL,
  latent_params = list(),
  item_params   = list(),
  verbose       = TRUE
)
```

#### Parameters

| Parameter       | Type            | Default                                      | Description                                                                                                         |
|:----------------|:----------------|:---------------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| `n_items`       | integer         | *required*                                   | Number of items.                                                                                                    |
| `model`         | character       | `"rasch"`                                    | `"rasch"` or `"2pl"`.                                                                                               |
| `latent_shape`  | character       | `"normal"`                                   | Latent distribution shape.                                                                                          |
| `item_source`   | character       | `"parametric"`                               | Item parameter source.                                                                                              |
| `c_bounds`      | numeric(2)      | `c(0.1, 10)`                                 | Scaling factor range to evaluate.                                                                                   |
| `M`             | integer         | `10000L`                                     | Monte Carlo sample size.                                                                                            |
| `seed`          | integer or NULL | `NULL`                                       | Random seed.                                                                                                        |
| `latent_params` | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).         |
| `item_params`   | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md). |
| `verbose`       | logical         | `TRUE`                                       | Print results.                                                                                                      |

#### Return Value

An object of class `"feasibility_check"` (returned invisibly) with:

| Element          | Description                    |
|:-----------------|:-------------------------------|
| `rho_range_info` | Achievable range of rho_tilde. |
| `rho_range_msem` | Achievable range of rho_bar.   |
| `n_items`        | Number of items.               |
| `model`          | Model.                         |
| `latent_shape`   | Latent shape.                  |
| `c_bounds`       | Evaluated range.               |
| `M`              | Sample size.                   |
| `theta_var`      | Estimated latent variance.     |

#### Example

``` r
feas <- check_feasibility(
  n_items = 25,
  model   = "rasch",
  M       = 5000L,
  seed    = 42,
  verbose = FALSE
)
cat(sprintf("rho_tilde range: [%.4f, %.4f]\n",
            feas$rho_range_info[1], feas$rho_range_info[2]))
#> rho_tilde range: [0.0591, 0.9872]
cat(sprintf("rho_bar range:   [%.4f, %.4f]\n",
            feas$rho_range_msem[1], feas$rho_range_msem[2]))
#> rho_bar range:   [0.0002, 0.9146]

# Is rho = 0.85 achievable with rho_tilde?
cat(sprintf("rho=0.85 feasible (info): %s\n",
            0.85 >= feas$rho_range_info[1] && 0.85 <= feas$rho_range_info[2]))
#> rho=0.85 feasible (info): TRUE
```

------------------------------------------------------------------------

### 3.6 `rho_curve()` – Reliability as a Function of Scaling Factor

Computes and optionally plots the reliability curve across a grid of
scaling factor values. Helps visualize the relationship between
discrimination strength and measurement precision.

#### Signature

``` r
rho_curve(
  c_values      = seq(0.1, 5, length.out = 50),
  n_items,
  model         = c("rasch", "2pl"),
  latent_shape  = "normal",
  item_source   = "parametric",
  metric        = c("both", "info", "msem"),
  M             = 5000L,
  seed          = NULL,
  latent_params = list(),
  item_params   = list(),
  plot          = TRUE
)
```

#### Parameters

| Parameter       | Type            | Default                                      | Description                                                                                                         |
|:----------------|:----------------|:---------------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| `c_values`      | numeric vector  | `seq(0.1, 5, length.out = 50)`               | Grid of scaling factor values.                                                                                      |
| `n_items`       | integer         | *required*                                   | Number of items.                                                                                                    |
| `model`         | character       | `"rasch"`                                    | `"rasch"` or `"2pl"`.                                                                                               |
| `latent_shape`  | character       | `"normal"`                                   | Latent distribution shape.                                                                                          |
| `item_source`   | character       | `"parametric"`                               | Item parameter source.                                                                                              |
| `metric`        | character       | `"both"`                                     | Which metric(s): `"both"`, `"info"`, or `"msem"`.                                                                   |
| `M`             | integer         | `5000L`                                      | Monte Carlo sample size.                                                                                            |
| `seed`          | integer or NULL | `NULL`                                       | Random seed.                                                                                                        |
| `latent_params` | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).         |
| `item_params`   | list            | [`list()`](https://rdrr.io/r/base/list.html) | Additional arguments for [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md). |
| `plot`          | logical         | `TRUE`                                       | Create a plot.                                                                                                      |

#### Return Value

A data frame of class `"rho_curve"` (returned invisibly when
`plot = TRUE`) with columns `c`, `rho_tilde`, and/or `rho_bar` depending
on the `metric` argument.

#### Example

``` r
curve_data <- rho_curve(
  c_values = seq(0.2, 4, length.out = 40),
  n_items  = 25,
  model    = "rasch",
  metric   = "both",
  M        = 5000L,
  seed     = 42,
  plot     = TRUE
)
```

![Reliability curve for a 25-item Rasch
test.](api-reference_files/figure-html/rho-curve-plot-1.png)

Reliability curve for a 25-item Rasch test.

``` r
# Inspect the data
head(curve_data)
#> Reliability Curve
#> =================
#>   Items: 25 | Model: RASCH | Metric: both
#>   c range: [0.20, 0.69] (6 points)
#>   rho_tilde range: [0.1989, 0.7174]
#>   rho_bar range  : [0.1989, 0.7139]
#> 
#>           c rho_tilde   rho_bar
#> 1 0.2000000 0.1988965 0.1988657
#> 2 0.2974359 0.3500759 0.3498742
#> 3 0.3948718 0.4807284 0.4800963
#> 4 0.4923077 0.5827614 0.5814095
#> 5 0.5897436 0.6596796 0.6573608
#> 6 0.6871795 0.7174113 0.7139425
```

------------------------------------------------------------------------

### 3.7 `compute_reliability_tam()` – TAM-Based Validation

Fits a Rasch or 2PL model using the **TAM** package and returns WLE and
EAP reliability estimates. This function requires the **TAM** package to
be installed.

#### Signature

``` r
compute_reliability_tam(resp, model = c("rasch", "2pl"), verbose = FALSE, ...)
```

#### Parameters

| Parameter | Type                 | Default    | Description                                           |
|:----------|:---------------------|:-----------|:------------------------------------------------------|
| `resp`    | matrix or data.frame | *required* | Binary response matrix (0/1).                         |
| `model`   | character            | `"rasch"`  | `"rasch"` or `"2pl"`.                                 |
| `verbose` | logical              | `FALSE`    | Print TAM fitting messages.                           |
| `...`     |                      |            | Additional arguments passed to TAM fitting functions. |

#### Return Value

A list with components:

| Element   | Description                                                               |
|:----------|:--------------------------------------------------------------------------|
| `rel_wle` | WLE reliability.                                                          |
| `rel_eap` | EAP reliability.                                                          |
| `mod`     | Fitted TAM model object.                                                  |
| `wle`     | Output from [`TAM::tam.wle()`](https://rdrr.io/pkg/TAM/man/tam.wle.html). |

#### Example

``` r
# Requires TAM package
tam_rel <- compute_reliability_tam(
  resp    = sim_data$response_matrix,
  model   = "rasch",
  verbose = FALSE
)
cat(sprintf("WLE reliability: %.4f\n", tam_rel$rel_wle))
cat(sprintf("EAP reliability: %.4f\n", tam_rel$rel_eap))
```

See
[`vignette("validation")`](https://joonho112.github.io/IRTsimrel/articles/validation.md)
for a complete validation workflow.

------------------------------------------------------------------------

## 4. S3 Methods

IRTsimrel defines S3 methods for six object classes. This section
documents each method with its signature and a brief example.

### 4.1 Methods for `eqc_result` Objects

Objects returned by
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md).

#### `print.eqc_result()`

``` r
print.eqc_result(x, digits = 4, ...)
```

Displays calibration summary, design parameters, convergence
diagnostics, and parameter summaries. Returns `x` invisibly.

``` r
print(eqc_res)
#> 
#> =======================================================
#>   Empirical Quadrature Calibration (EQC) Results
#> =======================================================
#> 
#> Calibration Summary:
#>   Model                        : RASCH
#>   Target reliability (rho*)    : 0.8000
#>   Achieved reliability         : 0.8000
#>   Absolute error               : 1.19e-07
#>   Scaling factor (c*)          : 0.8995
#> 
#> Design Parameters:
#>   Number of items (I)          : 25
#>   Quadrature points (M)        : 5000
#>   Reliability metric           : Average-information (tilde)
#>   Latent variance              : 1.0099
#> 
#> Convergence:
#>   Root status                  : uniroot_success
#>   Search bracket               : [0.300, 3.000]
#>   Bracket reliabilities        : [0.3539, 0.9550]
#> 
#> Parameter Summaries:
#>   theta:        mean = -0.014, sd = 1.005
#>   beta:         mean = 0.000, sd = 0.861, range = [-2.17, 1.45]
#>   lambda_base:  mean = 1.000, sd = 0.000
#>   lambda_scaled: mean = 0.899, sd = 0.000
```

#### `summary.eqc_result()`

``` r
summary.eqc_result(object, ...)
```

Returns an object of class `"summary.eqc_result"` containing a compact
subset of key results.

``` r
s <- summary(eqc_res)
s
#> Summary: Empirical Quadrature Calibration (EQC)
#> ================================================
#>   Model            : RASCH
#>   Metric           : Average-information (tilde)
#>   Number of items  : 25
#>   Quadrature (M)   : 5000
#>   Latent variance  : 1.0099
#> 
#> Calibration Results:
#>   Target rho*      : 0.8000
#>   Achieved rho     : 0.8000
#>   Absolute error   : 1.19e-07
#>   Scaling factor c*: 0.8995
#>   Root status      : uniroot_success
```

#### `coef.eqc_result()`

``` r
coef.eqc_result(object, ...)
```

Returns a data frame with columns: `item_id`, `beta`, `lambda_base`,
`lambda_scaled`, and `c_star`.

``` r
item_df <- coef(eqc_res)
head(item_df)
#>   item_id         beta lambda_base lambda_scaled    c_star
#> 1       1  0.197732269           1     0.8994993 0.8994993
#> 2       2  1.096799859           1     0.8994993 0.8994993
#> 3       3  0.436545084           1     0.8994993 0.8994993
#> 4       4 -0.013038730           1     0.8994993 0.8994993
#> 5       5 -0.199801302           1     0.8994993 0.8994993
#> 6       6  0.007700326           1     0.8994993 0.8994993
```

#### `predict.eqc_result()`

``` r
predict.eqc_result(object, newdata = NULL, ...)
```

If `newdata` is `NULL`, returns `achieved_rho`. If `newdata` is a
numeric vector of scaling factor values, computes reliability at each
value using the stored quadrature sample.

``` r
# Achieved reliability
predict(eqc_res)
#> [1] 0.8000001

# Reliability at several scaling factors
predict(eqc_res, newdata = c(0.5, 1.0, 1.5, 2.0))
#>     c=0.5     c=1.0     c=1.5     c=2.0 
#> 0.5896711 0.8259511 0.8972175 0.9280310
```

------------------------------------------------------------------------

### 4.2 Methods for `sac_result` Objects

Objects returned by
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).

#### `print.sac_result()`

``` r
print.sac_result(x, digits = 4, ...)
```

Displays calibration summary, algorithm settings, and convergence
diagnostics.

``` r
print(sac_res)
#> 
#> =======================================================
#>   Stochastic Approximation Calibration (SAC) Results
#> =======================================================
#> 
#> Calibration Summary:
#>   Model                        : RASCH
#>   Target reliability (rho*)    : 0.8000
#>   Achieved reliability         : 0.8292
#>   Absolute error               : 2.92e-02
#>   Scaling factor (c*)          : 1.0504
#> 
#> Algorithm Settings:
#>   Number of items (I)          : 25
#>   M per iteration              : 500
#>   M for variance pre-calc      : 5000
#>   Total iterations             : 200
#>   Burn-in                      : 100
#>   Reliability metric           : MSEM-based (bar/w)
#>   Step params: a=1.00, A=50, gamma=0.67
#> 
#> Convergence Diagnostics:
#>   Initialization method        : apc_warm_start
#>   Initial c_0                  : 1.3274
#>   Final iterate c_n            : 1.0195
#>   Polyak-Ruppert c*            : 1.0504
#>   Pre-calculated theta_var     : 1.0099
#>   Converged                    : Yes
#>   Post-burn-in SD              : 0.0210
#>   Final gradient (rho - rho*)  : +0.0194
```

#### `summary.sac_result()`

``` r
summary.sac_result(object, ...)
```

Returns a `"summary.sac_result"` object with compact results.

``` r
s_sac <- summary(sac_res)
s_sac
#> Summary: Stochastic Approximation Calibration (SAC)
#> ====================================================
#>   Model            : RASCH
#>   Metric           : MSEM-based (bar/w)
#>   Number of items  : 25
#>   Iterations       : 200
#>   Burn-in          : 100
#>   M per iteration  : 500
#>   M pre-calc       : 5000
#>   Init method      : apc_warm_start
#> 
#> Calibration Results:
#>   Target rho*      : 0.8000
#>   Achieved rho     : 0.8292
#>   Absolute error   : 2.92e-02
#>   Scaling factor c*: 1.0504
#>   Converged        : Yes
#>   Post-burn-in SD  : 0.0210
```

#### `coef.sac_result()`

``` r
coef.sac_result(object, ...)
```

Returns a data frame with columns: `item_id`, `beta`, `lambda_base`,
`lambda_scaled`, and `c_star`.

``` r
head(coef(sac_res))
#>   item_id        beta lambda_base lambda_scaled  c_star
#> 1       1  0.03774168           1       1.05045 1.05045
#> 2       2 -1.18785344           1       1.05045 1.05045
#> 3       3 -0.21413874           1       1.05045 1.05045
#> 4       4 -0.81661157           1       1.05045 1.05045
#> 5       5  1.52992507           1       1.05045 1.05045
#> 6       6 -0.38999393           1       1.05045 1.05045
```

#### `predict.sac_result()`

``` r
predict.sac_result(object, newdata = NULL, theta_vec = NULL, ...)
```

Like
[`predict.eqc_result()`](https://joonho112.github.io/IRTsimrel/reference/predict.eqc_result.md),
but accepts an optional `theta_vec` argument for using a different
ability sample.

``` r
predict(sac_res)
#> [1] 0.8291582
predict(sac_res, newdata = c(0.5, 1.0, 1.5))
#>     c=0.5     c=1.0     c=1.5 
#> 0.5886405 0.8189501 0.8833365
```

#### `plot.sac_result()`

``` r
plot.sac_result(x, type = c("both", "trajectory", "c", "rho"), ...)
```

Plots the Robbins–Monro convergence trajectory. Types `"trajectory"` and
`"c"` are synonyms showing the scaling factor path. Type `"rho"` shows
reliability estimates across iterations. Type `"both"` combines them.
Uses **ggplot2** if available, falling back to base R graphics.

``` r
plot(sac_res, type = "c")
```

![SAC convergence trajectory (scaling
factor).](api-reference_files/figure-html/plot-sac-traj-1.png)

SAC convergence trajectory (scaling factor).

``` r
plot(sac_res, type = "rho")
```

![SAC reliability estimates across
iterations.](api-reference_files/figure-html/plot-sac-rho-1.png)

SAC reliability estimates across iterations.

------------------------------------------------------------------------

### 4.3 Methods for `latent_G` Objects

Objects returned by
[`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

#### `print.latent_G()`

``` r
print.latent_G(x, ...)
```

Prints shape, sample size, target parameters, and sample moments.

``` r
print(g_bimod)
#> Latent Ability Distribution (G-family)
#> =======================================
#>   Shape     : bimodal
#>   n         : 2000
#>   Target mu : 0.000
#>   Target sigma: 1.000
#> 
#> Sample Moments:
#>   Mean      : -0.0172
#>   SD        : 0.9885
#>   Skewness  : 0.0126
#>   Kurtosis  : -1.2981 (excess)
```

#### `summary.latent_G()`

``` r
summary.latent_G(object, ...)
```

Returns a `"summary.latent_G"` object with detailed statistics including
quantiles.

``` r
summary(g_bimod)
#> Summary: Latent Ability Distribution
#> ====================================
#>   Shape      : bimodal
#>   n          : 2000
#>   Target     : mu = 0.00, sigma = 1.00
#>   Covariates : No
#> 
#> Sample Statistics:
#>   Mean       : -0.0172
#>   SD         : 0.9885
#>   Median     : -0.0657
#>   Skewness   : 0.0126
#>   Kurtosis   : -1.2981 (excess)
#>   Range      : [-2.2867, 2.4625]
#> 
#> Quantiles:
#>    2.5%      5%     25%     50%     75%     95%   97.5% 
#> -1.6338 -1.4650 -0.9025 -0.0657  0.8813  1.4319  1.5588
```

#### `plot.latent_G()`

``` r
plot.latent_G(x, type = c("both", "histogram", "density"),
              show_normal = TRUE, bins = 50, ...)
```

Plots the latent distribution as histogram, density, or both. Uses
**ggplot2** if available.

| Parameter     | Default  | Description                                         |
|:--------------|:---------|:----------------------------------------------------|
| `type`        | `"both"` | Plot type: `"both"`, `"histogram"`, or `"density"`. |
| `show_normal` | `TRUE`   | Overlay a normal reference density (dashed red).    |
| `bins`        | `50`     | Number of histogram bins.                           |

``` r
plot(g_bimod, type = "both", bins = 40)
```

![Bimodal latent ability
distribution.](api-reference_files/figure-html/plot-latentg-1.png)

Bimodal latent ability distribution.

#### `as.numeric.latent_G()`

``` r
as.numeric.latent_G(x, ...)
```

Extracts the `theta` vector for use with other functions.

``` r
# Extract theta vector directly (equivalent to as.numeric dispatch)
theta_vec <- g_norm$theta
cat(sprintf("Length: %d, Mean: %.3f\n", length(theta_vec), mean(theta_vec)))
#> Length: 2000, Mean: -0.016
```

------------------------------------------------------------------------

### 4.4 Methods for `item_params` Objects

Objects returned by
[`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

#### `print.item_params()`

``` r
print.item_params(x, ...)
```

Prints model, source, method, summary statistics of difficulties and
discriminations, and achieved correlations.

``` r
print(items_2pl)
#> Item Parameters Object
#> ======================
#>   Model          : 2PL
#>   Source         : parametric
#>   Method         : copula
#>   Items per form : 30
#>   Number of forms: 1
#>   Scale factor   : 1.000
#>   Centered       : Yes
#> 
#> Difficulty (beta):
#>   Mean: 0.0000, SD: 1.2550, Range: [-2.725, 2.218]
#> 
#> Discrimination (lambda, scaled):
#>   Mean: 1.0867, SD: 0.2885, Range: [0.440, 1.731]
#> 
#> Correlation (beta, log-lambda):
#>   Target (rho): -0.300
#>   Achieved Pearson : -0.361
#>   Achieved Spearman: -0.378
```

#### `summary.item_params()`

``` r
summary.item_params(object, ...)
```

Returns a `"summary.item_params"` object with detailed parameter
summaries.

``` r
summary(items_2pl)
#> Summary: Item Parameters
#> ========================
#>   Model          : 2PL
#>   Source         : parametric
#>   Method         : copula
#>   Items per form : 30
#>   Number of forms: 1
#>   Scale factor   : 1.000
#>   Centered       : Yes
#> 
#> Difficulty (beta):
#>   Mean     : 0.0000
#>   SD       : 1.2550
#>   Min      : -2.7250
#>   Max      : 2.2181
#>   Quantiles: Q25=-0.4681, Q50=-0.1690, Q75=1.0014
#> 
#> Discrimination (lambda):
#>   Before scaling: Mean=1.0867, SD=0.2885
#>   After scaling (c=1.000): Mean=1.0867, SD=0.2885
#>   Range [0.4396, 1.7309]
#> 
#> Correlation (beta, log-lambda):
#>   Target (rho)     : -0.3000
#>   Achieved Pearson : -0.3611
#>   Achieved Spearman: -0.3784
```

#### `plot.item_params()`

``` r
plot.item_params(x, type = c("scatter", "density", "both"), ...)
```

For 2PL models, creates difficulty vs. discrimination scatter plots
and/or density plots. Uses **ggplot2** if available.

| Parameter | Default     | Description                            |
|:----------|:------------|:---------------------------------------|
| `type`    | `"scatter"` | `"scatter"`, `"density"`, or `"both"`. |

``` r
plot(items_2pl, type = "scatter")
#> `geom_smooth()` using formula = 'y ~ x'
```

![Difficulty vs. discrimination scatter
plot.](api-reference_files/figure-html/plot-items-scatter-1.png)

Difficulty vs. discrimination scatter plot.

#### `as.data.frame.item_params()`

``` r
as.data.frame.item_params(x, row.names = NULL, optional = FALSE, ...)
```

Extracts the item parameter data frame.

``` r
df <- as.data.frame(items_rasch)
head(df)
#>   form_id item_id       beta lambda lambda_unscaled
#> 1       1       1  1.1834223      1               1
#> 2       1       2 -0.7522343      1               1
#> 3       1       3  0.1755922      1               1
#> 4       1       4  0.4453264      1               1
#> 5       1       5  0.2167322      1               1
#> 6       1       6 -0.2936607      1               1
```

------------------------------------------------------------------------

### 4.5 Methods for `feasibility_check` Objects

Objects returned by
[`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md).

#### `print.feasibility_check()`

``` r
print.feasibility_check(x, digits = 4, ...)
```

Prints achievable reliability ranges and design information.

``` r
print(feas)
#> 
#> =======================================================
#>   Feasibility Check: Achievable Reliability Range
#> =======================================================
#> 
#>   Number of items  : 25
#>   Model            : RASCH
#>   Latent shape     : normal
#>   Latent variance  : 1.0099
#>   c range          : [0.10, 10.00]
#>   Monte Carlo M    : 5000
#> 
#> Achievable Reliability Ranges:
#>   rho_tilde (info) : [0.0591, 0.9872]
#>   rho_bar   (msem) : [0.0002, 0.9146]
#> 
#> Note: rho_tilde >= rho_bar always (Jensen's inequality).
#>   Use rho_tilde range for EQC targets.
#>   Use rho_bar range for SAC targets.
```

------------------------------------------------------------------------

### 4.6 Methods for `rho_curve` Objects

Objects returned by
[`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md).

#### `print.rho_curve()`

``` r
print.rho_curve(x, ...)
```

Prints a compact summary of the reliability curve data and the first few
rows.

``` r
print(curve_data)
#> Reliability Curve
#> =================
#>   Items: 25 | Model: RASCH | Metric: both
#>   c range: [0.20, 4.00] (40 points)
#>   rho_tilde range: [0.1989, 0.9672]
#>   rho_bar range  : [0.1989, 0.9146]
#> 
#>           c rho_tilde   rho_bar
#> 1 0.2000000 0.1988965 0.1988657
#> 2 0.2974359 0.3500759 0.3498742
#> 3 0.3948718 0.4807284 0.4800963
#> 4 0.4923077 0.5827614 0.5814095
#> 5 0.5897436 0.6596796 0.6573608
#> 6 0.6871795 0.7174113 0.7139425
#>   ... (34 more rows)
```

------------------------------------------------------------------------

### 4.7 Summary Print Methods

Each summary class has its own `print` method that displays formatted
output.

#### `print.summary.eqc_result()`

``` r
print.summary.eqc_result(x, digits = 4, ...)
```

``` r
print(summary(eqc_res))
#> Summary: Empirical Quadrature Calibration (EQC)
#> ================================================
#>   Model            : RASCH
#>   Metric           : Average-information (tilde)
#>   Number of items  : 25
#>   Quadrature (M)   : 5000
#>   Latent variance  : 1.0099
#> 
#> Calibration Results:
#>   Target rho*      : 0.8000
#>   Achieved rho     : 0.8000
#>   Absolute error   : 1.19e-07
#>   Scaling factor c*: 0.8995
#>   Root status      : uniroot_success
```

#### `print.summary.sac_result()`

``` r
print.summary.sac_result(x, digits = 4, ...)
```

``` r
print(summary(sac_res))
#> Summary: Stochastic Approximation Calibration (SAC)
#> ====================================================
#>   Model            : RASCH
#>   Metric           : MSEM-based (bar/w)
#>   Number of items  : 25
#>   Iterations       : 200
#>   Burn-in          : 100
#>   M per iteration  : 500
#>   M pre-calc       : 5000
#>   Init method      : apc_warm_start
#> 
#> Calibration Results:
#>   Target rho*      : 0.8000
#>   Achieved rho     : 0.8292
#>   Absolute error   : 2.92e-02
#>   Scaling factor c*: 1.0504
#>   Converged        : Yes
#>   Post-burn-in SD  : 0.0210
```

#### `print.summary.item_params()`

``` r
print.summary.item_params(x, digits = 4, ...)
```

``` r
print(summary(items_2pl))
#> Summary: Item Parameters
#> ========================
#>   Model          : 2PL
#>   Source         : parametric
#>   Method         : copula
#>   Items per form : 30
#>   Number of forms: 1
#>   Scale factor   : 1.000
#>   Centered       : Yes
#> 
#> Difficulty (beta):
#>   Mean     : 0.0000
#>   SD       : 1.2550
#>   Min      : -2.7250
#>   Max      : 2.2181
#>   Quantiles: Q25=-0.4681, Q50=-0.1690, Q75=1.0014
#> 
#> Discrimination (lambda):
#>   Before scaling: Mean=1.0867, SD=0.2885
#>   After scaling (c=1.000): Mean=1.0867, SD=0.2885
#>   Range [0.4396, 1.7309]
#> 
#> Correlation (beta, log-lambda):
#>   Target (rho)     : -0.3000
#>   Achieved Pearson : -0.3611
#>   Achieved Spearman: -0.3784
```

#### `print.summary.latent_G()`

``` r
print.summary.latent_G(x, digits = 4, ...)
```

``` r
print(summary(g_norm))
#> Summary: Latent Ability Distribution
#> ====================================
#>   Shape      : normal
#>   n          : 2000
#>   Target     : mu = 0.00, sigma = 1.00
#>   Covariates : No
#> 
#> Sample Statistics:
#>   Mean       : -0.0156
#>   SD         : 0.9941
#>   Median     : -0.0131
#>   Skewness   : 0.0128
#>   Kurtosis   : 0.0560 (excess)
#>   Range      : [-3.3717, 3.5847]
#> 
#> Quantiles:
#>    2.5%      5%     25%     50%     75%     95%   97.5% 
#> -1.9835 -1.6446 -0.6691 -0.0131  0.6608  1.5768  1.8871
```

------------------------------------------------------------------------

## 5. Deprecated Functions

These functions are retained for backward compatibility and will be
removed in a future release. They issue a deprecation warning when
called.

### 5.1 `spc_calibrate()`

``` r
spc_calibrate(...)
```

Deprecated alias for
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).
All arguments are passed through unchanged. Use
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
for new code.

``` r
# Deprecated: use sac_calibrate() instead
result <- spc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  n_iter = 200, seed = 42
)
```

### 5.2 `compare_eqc_spc()`

``` r
compare_eqc_spc(...)
```

Deprecated alias for
[`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md).
All arguments are passed through unchanged. Use
[`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
for new code.

``` r
# Deprecated: use compare_eqc_sac() instead
compare_eqc_spc(eqc_result, sac_result)
```

------------------------------------------------------------------------

## 6. Cross-Reference Table

The table below maps each function to the vignette that covers it in
detail.

| Function                                                                                                  | Category    | Primary Vignette                                                                                             |
|:----------------------------------------------------------------------------------------------------------|:------------|:-------------------------------------------------------------------------------------------------------------|
| [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)                     | Calibration | [`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)               |
| [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)                     | Calibration | [`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md)               |
| [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)                 | Calibration | [`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)               |
| [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)                         | Simulation  | [`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md) |
| [`compare_shapes()`](https://joonho112.github.io/IRTsimrel/reference/compare_shapes.md)                   | Simulation  | [`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md) |
| [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)                 | Simulation  | [`vignette("item-parameters")`](https://joonho112.github.io/IRTsimrel/articles/item-parameters.md)           |
| [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)   | Simulation  | [`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)               |
| [`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)                 | Reliability | [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)     |
| [`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)               | Reliability | [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)     |
| [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)               | Reliability | [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)     |
| [`compute_apc_init()`](https://joonho112.github.io/IRTsimrel/reference/compute_apc_init.md)               | Reliability | [`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md)               |
| [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)             | Reliability | [`vignette("simulation-design")`](https://joonho112.github.io/IRTsimrel/articles/simulation-design.md)       |
| [`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)                             | Reliability | [`vignette("simulation-design")`](https://joonho112.github.io/IRTsimrel/articles/simulation-design.md)       |
| [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md) | Validation  | [`vignette("validation")`](https://joonho112.github.io/IRTsimrel/articles/validation.md)                     |
| [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)                     | Deprecated  | –                                                                                                            |
| [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)                 | Deprecated  | –                                                                                                            |

------------------------------------------------------------------------

## 7. Complete Workflow Example

This section demonstrates a full simulation study using the IRTsimrel
API.

### Step 1: Check Feasibility

``` r
feas_check <- check_feasibility(
  n_items = 30,
  model   = "rasch",
  M       = 5000L,
  seed    = 42,
  verbose = FALSE
)
cat(sprintf("For 30 Rasch items, achievable rho_tilde: [%.3f, %.3f]\n",
            feas_check$rho_range_info[1], feas_check$rho_range_info[2]))
#> For 30 Rasch items, achievable rho_tilde: [0.070, 0.989]
```

### Step 2: Calibrate with EQC

``` r
eqc_wf <- eqc_calibrate(
  target_rho = 0.85,
  n_items    = 30,
  model      = "rasch",
  M          = 5000L,
  seed       = 42
)
cat(sprintf("EQC: c* = %.4f, rho = %.4f\n",
            eqc_wf$c_star, eqc_wf$achieved_rho))
#> EQC: c* = 0.9943, rho = 0.8500
```

### Step 3: Validate with SAC

``` r
sac_wf <- sac_calibrate(
  target_rho = 0.85,
  n_items    = 30,
  model      = "rasch",
  c_init     = eqc_wf,
  n_iter     = 150L,
  M_per_iter = 500L,
  M_pre      = 5000L,
  seed       = 42,
  verbose    = FALSE
)
cat(sprintf("SAC: c* = %.4f, rho = %.4f\n",
            sac_wf$c_star, sac_wf$achieved_rho))
#> SAC: c* = 1.0298, rho = 0.8417
```

### Step 4: Compare Results

``` r
comp_wf <- compare_eqc_sac(eqc_wf, sac_wf, verbose = TRUE)
#> Warning in compare_eqc_sac(eqc_wf, sac_wf, verbose = TRUE): Reliability metric
#> differs between EQC ('info') and SAC ('msem').
#> 
#> =======================================================
#>   EQC vs SAC Comparison
#> =======================================================
#> 
#>   Target reliability  : 0.8500
#>   EQC c*              : 0.994331
#>   SAC c*              : 1.029780
#>   Absolute difference : 0.035449
#>   Percent difference  : 3.57%
#>   Agreement (< 5%)    : YES
#> 
```

### Step 5: Generate Response Data

``` r
resp_wf <- simulate_response_data(
  result    = eqc_wf,
  n_persons = 1000,
  seed      = 123
)
cat(sprintf("Generated %d x %d response matrix\n",
            nrow(resp_wf$response_matrix),
            ncol(resp_wf$response_matrix)))
#> Generated 1000 x 30 response matrix
cat(sprintf("Mean score: %.2f / %d items\n",
            mean(rowSums(resp_wf$response_matrix)),
            ncol(resp_wf$response_matrix)))
#> Mean score: 15.00 / 30 items
```

### Step 6: Extract Calibrated Parameters

``` r
params_wf <- coef(eqc_wf)
head(params_wf)
#>   item_id        beta lambda_base lambda_scaled    c_star
#> 1       1  0.17744288           1     0.9943313 0.9943313
#> 2       2  1.07651047           1     0.9943313 0.9943313
#> 3       3  0.41625569           1     0.9943313 0.9943313
#> 4       4 -0.03332812           1     0.9943313 0.9943313
#> 5       5 -0.22009069           1     0.9943313 0.9943313
#> 6       6 -0.01258907           1     0.9943313 0.9943313
cat(sprintf("\nScaled discrimination: mean = %.3f, sd = %.3f\n",
            mean(params_wf$lambda_scaled),
            sd(params_wf$lambda_scaled)))
#> 
#> Scaled discrimination: mean = 0.994, sd = 0.000
```

### Step 7: Explore Reliability Curve

``` r
rc_wf <- rho_curve(
  c_values = seq(0.2, 3, length.out = 30),
  n_items  = 30,
  model    = "rasch",
  M        = 5000L,
  seed     = 42,
  plot     = TRUE
)
```

![Reliability curve for the workflow
example.](api-reference_files/figure-html/workflow-rho-curve-1.png)

Reliability curve for the workflow example.

### Step 8: Validate with TAM (optional)

``` r
# Requires TAM package
tam_wf <- compute_reliability_tam(
  resp  = resp_wf$response_matrix,
  model = "rasch"
)
cat(sprintf("WLE reliability: %.4f\n", tam_wf$rel_wle))
cat(sprintf("EAP reliability: %.4f\n", tam_wf$rel_eap))
cat(sprintf("Target reliability: %.4f\n", eqc_wf$target_rho))
```

------------------------------------------------------------------------

## 8. Complete Export Inventory

The table below lists every symbol exported from IRTsimrel (from the
NAMESPACE file), organized alphabetically.

| Export                      | Type                  |
|:----------------------------|:----------------------|
| `as.data.frame.item_params` | S3 method             |
| `as.numeric.latent_G`       | S3 method             |
| `check_feasibility`         | Function              |
| `coef.eqc_result`           | S3 method             |
| `coef.sac_result`           | S3 method             |
| `compare_eqc_sac`           | Function              |
| `compare_eqc_spc`           | Function (deprecated) |
| `compare_shapes`            | Function              |
| `compute_apc_init`          | Function              |
| `compute_reliability_tam`   | Function              |
| `compute_rho_bar`           | Function              |
| `compute_rho_both`          | Function              |
| `compute_rho_tilde`         | Function              |
| `eqc_calibrate`             | Function              |
| `plot.item_params`          | S3 method             |
| `plot.latent_G`             | S3 method             |
| `plot.sac_result`           | S3 method             |
| `predict.eqc_result`        | S3 method             |
| `predict.sac_result`        | S3 method             |
| `print.eqc_result`          | S3 method             |
| `print.feasibility_check`   | S3 method             |
| `print.item_params`         | S3 method             |
| `print.latent_G`            | S3 method             |
| `print.rho_curve`           | S3 method             |
| `print.sac_result`          | S3 method             |
| `print.summary.eqc_result`  | S3 method             |
| `print.summary.item_params` | S3 method             |
| `print.summary.latent_G`    | S3 method             |
| `print.summary.sac_result`  | S3 method             |
| `rho_curve`                 | Function              |
| `sac_calibrate`             | Function              |
| `sim_item_params`           | Function              |
| `sim_latentG`               | Function              |
| `simulate_response_data`    | Function              |
| `spc_calibrate`             | Function (deprecated) |
| `summary.eqc_result`        | S3 method             |
| `summary.item_params`       | S3 method             |
| `summary.latent_G`          | S3 method             |
| `summary.sac_result`        | S3 method             |

------------------------------------------------------------------------

## Session Information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] IRTsimrel_0.2.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.7.1        nlme_3.1-168       cli_3.6.5          knitr_1.51        
#>  [5] rlang_1.1.7        xfun_0.56          S7_0.2.1           textshaping_1.0.4 
#>  [9] jsonlite_2.0.0     labeling_0.4.3     glue_1.8.0         htmltools_0.5.9   
#> [13] ragg_1.5.0         sass_0.4.10        scales_1.4.0       rmarkdown_2.30    
#> [17] grid_4.5.2         evaluate_1.0.5     jquerylib_0.1.4    MASS_7.3-65       
#> [21] fastmap_1.2.0      yaml_2.3.12        lifecycle_1.0.5    compiler_4.5.2    
#> [25] RColorBrewer_1.1-3 fs_1.6.6           mgcv_1.9-3         lattice_0.22-7    
#> [29] farver_2.1.2       systemfonts_1.3.1  digest_0.6.39      R6_2.6.1          
#> [33] splines_4.5.2      Matrix_1.7-4       bslib_0.10.0       withr_3.0.2       
#> [37] tools_4.5.2        gtable_0.3.6       pkgdown_2.2.0      ggplot2_4.0.2     
#> [41] cachem_1.1.0       desc_1.4.3
```
