# Compute Reliability as a Function of Scaling Factor

Computes and optionally plots the reliability curve \\\rho(c)\\ across a
grid of scaling factor values. This visualization helps understand how
reliability varies with the discrimination scaling factor and aids in
selecting appropriate target reliability values.

## Usage

``` r
rho_curve(
  c_values = seq(0.1, 5, length.out = 50),
  n_items,
  model = c("rasch", "2pl"),
  latent_shape = "normal",
  item_source = "parametric",
  metric = c("both", "info", "msem"),
  M = 5000L,
  seed = NULL,
  latent_params = list(),
  item_params = list(),
  plot = TRUE
)

# S3 method for class 'rho_curve'
print(x, ...)
```

## Arguments

- c_values:

  Numeric vector. Grid of scaling factor values to evaluate. Default:
  `seq(0.1, 5, length.out = 50)`.

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

- metric:

  Character. Which reliability metric(s) to compute: `"both"`, `"info"`,
  or `"msem"`.

- M:

  Integer. Monte Carlo sample size. Default: 5000.

- seed:

  Optional integer for reproducibility.

- latent_params:

  List. Additional arguments passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- item_params:

  List. Additional arguments passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- plot:

  Logical. If TRUE (default), create a plot of the curve.

- x:

  An object of class `"rho_curve"`.

- ...:

  Additional arguments (ignored).

## Value

A data frame of class `"rho_curve"` with columns:

- `c`:

  Scaling factor values.

- `rho_tilde`:

  Average-information reliability (if metric includes "info").

- `rho_bar`:

  MSEM-based reliability (if metric includes "msem").

The input object, invisibly.

## Details

The function generates a single set of theta and item parameters, then
evaluates the reliability at each value of `c_values`. When
`metric = "both"`, it uses
[`compute_rho_both`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)
for efficiency.

## See also

[`check_feasibility`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md),
[`compute_rho_both`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)

## Examples

``` r
# Basic usage: plot reliability curve for 25-item Rasch test
curve_data <- rho_curve(n_items = 25, model = "rasch", seed = 42,
                        M = 3000, plot = FALSE)
head(curve_data)
#> Reliability Curve
#> =================
#>   Items: 25 | Model: RASCH | Metric: both
#>   c range: [0.10, 0.60] (6 points)
#>   rho_tilde range: [0.0592, 0.6645]
#>   rho_bar range  : [0.0592, 0.6622]
#> 
#>     c  rho_tilde    rho_bar
#> 1 0.1 0.05922061 0.05921989
#> 2 0.2 0.19896249 0.19893104
#> 3 0.3 0.35354114 0.35333047
#> 4 0.4 0.48588780 0.48523096
#> 5 0.5 0.58813434 0.58675218
#> 6 0.6 0.66450137 0.66217929
```
