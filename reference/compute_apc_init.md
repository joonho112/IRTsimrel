# Analytic Pre-Calibration (APC) Initialization

Computes an initial value for the scaling factor using the closed-form
approximation under Gaussian Rasch assumptions.

## Usage

``` r
compute_apc_init(target_rho, n_items, sigma_beta = 1)
```

## Arguments

- target_rho:

  Numeric. Target reliability.

- n_items:

  Integer. Number of items.

- sigma_beta:

  Numeric. SD of item difficulties (default: 1.0).

## Value

Numeric. Initial scaling factor c_init.

## Details

Under the Gaussian Rasch setting with \\\theta \sim N(0,1)\\ and \\\beta
\sim N(0, \sigma\_\beta^2)\\, the expected item information involves the
logistic-normal convolution: \$\$\kappa(\sigma^2) = \int
\frac{e^z}{(1+e^z)^2} \phi(z; 0, \sigma^2) dz\$\$

Approximating \\\kappa \approx 0.25 / \sqrt{1 + \sigma^2 \pi^2/3}\\, the
closed-form pre-calibration is: \$\$c\_{init} = \sqrt{\frac{\rho^\*}{I
\cdot \kappa \cdot (1 - \rho^\*)}}\$\$

## Examples

``` r
# Compute initial c for target reliability of 0.80 with 25 items
compute_apc_init(target_rho = 0.80, n_items = 25)
#> [1] 1.327405

# With different difficulty spread
compute_apc_init(target_rho = 0.75, n_items = 20, sigma_beta = 1.5)
#> [1] 1.432348
```
