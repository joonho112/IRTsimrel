# Compute Marginal Reliability from Simulated Test Information

Low-level utilities for mapping a discrimination scale \\c\\ and
simulated person/item parameters to marginal reliability.

These functions implement the same reliability definitions used inside
both
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
and
[`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md):

- `compute_rho_bar()`: MSEM-based marginal reliability \\\bar{w}(c) =
  \sigma\_\theta^2 / (\sigma\_\theta^2 + E\[1/\mathcal{J}(\theta;c)\])\\

- `compute_rho_tilde()`: Average-information reliability
  \\\tilde{\rho}(c) = \sigma\_\theta^2 \bar{\mathcal{J}}(c) /
  (\sigma\_\theta^2 \bar{\mathcal{J}}(c) + 1)\\

## Usage

``` r
compute_rho_bar(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)

compute_rho_tilde(c, theta_vec, beta_vec, lambda_base, theta_var = NULL)
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

Numeric scalar reliability value.

## Details

The computation proceeds as:

1.  Form scaled discriminations \\\lambda_i(c) = c \cdot
    \lambda\_{i,0}\\.

2.  Compute item response probabilities \\p\_{mi} =
    \text{logit}^{-1}\\\lambda_i(c)(\theta_m - \beta_i)\\\\.

3.  Item information: \\\mathcal{J}\_{mi} = \lambda_i(c)^2
    p\_{mi}(1-p\_{mi})\\.

4.  Test information at each \\\theta_m\\: \\\mathcal{J}\_m = \sum_i
    \mathcal{J}\_{mi}\\.

5.  Reliability:

    - `compute_rho_bar()`: harmonic-mean-based MSEM \\\text{MSEM} =
      E\[1/\mathcal{J}\_m\]\\, \\\bar{w}(c) = \sigma\_\theta^2 /
      (\sigma\_\theta^2 + \text{MSEM})\\.

    - `compute_rho_tilde()`: arithmetic-mean-based information
      \\\bar{\mathcal{J}} = E\[\mathcal{J}\_m\]\\, \\\tilde{\rho}(c) =
      \sigma\_\theta^2 \bar{\mathcal{J}} / (\sigma\_\theta^2
      \bar{\mathcal{J}} + 1)\\.

A small floor (`1e-10`) is applied to test information to avoid
numerical problems when taking reciprocals.

## Examples

``` r
# Simple toy example
set.seed(1)
theta <- rnorm(1000)
beta  <- rnorm(20)
lambda0 <- rep(1, 20)

compute_rho_bar(1, theta, beta, lambda0)
#> [1] 0.7795393
compute_rho_tilde(1, theta, beta, lambda0)
#> [1] 0.7870079

# With pre-calculated theta variance (recommended for SPC)
theta_var_fixed <- var(rnorm(10000))  # Pre-calculate from large sample
compute_rho_bar(1, theta, beta, lambda0, theta_var = theta_var_fixed)
#> [1] 0.7684109
```
