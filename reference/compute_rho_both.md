# Compute Both Reliability Metrics in a Single Pass

Computes both the average-information reliability (\\\tilde{\rho}\\) and
the MSEM-based marginal reliability (\\\bar{w}\\) from a single set of
test information values, avoiding redundant matrix computation.

This is a performance optimization over calling
[`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
and
[`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
separately, since both share the same \\M \times I\\ matrix
computations.

## Usage

``` r
compute_rho_both(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)
```

## Arguments

- c:

  Numeric scalar. Global discrimination scaling factor.

- theta_vec:

  Numeric vector of abilities \\\theta_m\\.

- beta_vec:

  Numeric vector of item difficulties \\\beta_i\\.

- lambda_base:

  Numeric vector of baseline discriminations \\\lambda\_{i,0}\\ (before
  scaling by `c`).

- theta_var:

  Optional numeric. Pre-calculated variance of theta. If NULL, computed
  from `theta_vec`.

## Value

A named list with components:

- `rho_tilde`:

  Average-information reliability.

- `rho_bar`:

  MSEM-based marginal reliability.

## Details

By Jensen's inequality, \\\tilde{\rho} \geq \bar{w}\\ always holds. See
[`compute_rho_bar`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
for details on each metric.

## Examples

``` r
set.seed(1)
theta <- rnorm(1000)
beta  <- rnorm(20)
lambda0 <- rep(1, 20)

both <- compute_rho_both(1, theta, beta, lambda0)
both$rho_tilde
#> [1] 0.7870079
both$rho_bar
#> [1] 0.7795393

# Verify: rho_tilde >= rho_bar (Jensen's inequality)
both$rho_tilde >= both$rho_bar
#> [1] TRUE
```
