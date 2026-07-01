# Check Feasibility of Target Reliability

Screens whether a given target reliability is achievable for a
particular test design (number of items, model, latent distribution,
item source) by computing the range of achievable reliabilities across a
range of scaling factors.

This function is useful for determining whether a planned simulation
study is feasible before running the (potentially expensive) calibration
algorithms.

## Usage

``` r
check_feasibility(
  n_items,
  model = c("rasch", "2pl"),
  latent_shape = "normal",
  item_source = "parametric",
  c_bounds = c(0.1, 10),
  M = 10000L,
  seed = NULL,
  latent_params = list(),
  item_params = list(),
  target_rho = NULL,
  verbose = TRUE
)

# S3 method for class 'feasibility_check'
print(x, digits = 4, ...)
```

## Arguments

- n_items:

  Integer. Number of items in the test form.

- model:

  Character. Measurement model: `"rasch"` or `"2pl"`.

- latent_shape:

  Character. Shape argument passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- item_source:

  Character. Source argument passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- c_bounds:

  Numeric length-2 vector. Range of scaling factors to evaluate.
  Default: `c(0.1, 10)`.

- M:

  Integer. Monte Carlo sample size for theta. Default: 10000.

- seed:

  Optional integer for reproducibility.

- latent_params:

  List. Additional arguments passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- item_params:

  List. Additional arguments passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- target_rho:

  Optional numeric scalar in (0, 1). If supplied, the result includes
  status fields indicating whether the target is below, inside, or above
  the achievable range for each metric.

- verbose:

  Logical. If TRUE, print results.

- x:

  An object of class `"feasibility_check"`.

- digits:

  Integer. Number of decimal places for printing.

- ...:

  Additional arguments (ignored).

## Value

An object of class `"feasibility_check"` (a list) with:

- `rho_range_info`:

  Numeric length-2 vector: achievable range of average-information
  reliability (\\\tilde{\rho}\\).

- `rho_range_msem`:

  Numeric length-2 vector: achievable range of MSEM-based reliability
  (\\\bar{w}\\).

- `rho_bounds_info`:

  Endpoint reliabilities for average-information reliability.

- `rho_bounds_msem`:

  Endpoint reliabilities for MSEM-based reliability.

- `rho_info_max_c`:

  Scaling factor at the maximum sampled info reliability within
  `c_bounds`.

- `rho_msem_max_c`:

  Scaling factor at the maximum sampled MSEM reliability within
  `c_bounds`.

- `target_status_info`:

  If `target_rho` is supplied, one of `"below_lower"`,
  `"boundary_lower"`, `"feasible"`, `"boundary_upper"`, or
  `"above_upper"` for \\\tilde{\rho}\\.

- `target_status_msem`:

  Analogous status for \\\bar{w}\\.

- `n_items`:

  Number of items.

- `model`:

  Model used.

- `latent_shape`:

  Latent distribution shape.

- `c_bounds`:

  Scaling factor bounds evaluated.

- `M`:

  Monte Carlo sample size.

- `theta_var`:

  Estimated latent variance.

The input object, invisibly.

## Details

For the average-information metric (\\\tilde{\rho}\\), practical EQC
search intervals are usually monotone. With very wide finite empirical
quadrature bounds, however, logistic probabilities can saturate
numerically. The function therefore reports endpoint reliabilities
separately and uses an interior maximization pass for the feasible upper
range.

For the MSEM-based metric (\\\bar{w}\\), the reliability can be
non-monotone at extreme scaling factors. The function uses
[`optimize`](https://rdrr.io/r/stats/optimize.html) to find the maximum
within the bounds, and the range is \\\[\min(\bar{w}(c\_{min}),
\bar{w}(c\_{max})), \bar{w}\_{max}\]\\.

## See also

[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
[`rho_curve`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)

## Examples

``` r
# Check feasibility for 25-item Rasch test
feas <- check_feasibility(n_items = 25, model = "rasch",
                          target_rho = 0.90, seed = 42,
                          M = 5000, verbose = FALSE)
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
#> Target rho*        : 0.9000
#>   info status      : feasible
#>   msem status      : feasible
#> 
#> Note: rho_tilde >= rho_bar always (Jensen's inequality).
#>   Use rho_tilde range for EQC targets.
#>   Use rho_bar range for SAC targets.
#> 

# Metric-specific target status
feas$target_status_info
#> [1] "feasible"
feas$target_status_msem
#> [1] "feasible"
```
