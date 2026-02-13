# Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)

`eqc_calibrate()` implements Algorithm 1 (Empirical / Stochastic
Quadrature Calibration, EQC/SQC) for reliability-targeted IRT
simulation.

Given a target marginal reliability \\\rho^\*\\, a latent distribution
generator
[`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
(for \\G\\) and an item parameter generator
[`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
(for \\H\\), the function searches for a global discrimination scale
\\c^\* \> 0\\ such that the population reliability \\\rho(c)\\ of the
Rasch/2PL model is approximately equal to \\\rho^\*\\.

The key idea is to:

1.  Draw a large fixed "quadrature" sample \\\\\theta_m\\\_{m=1}^M \sim
    G\\ and item parameters \\\\(\beta_i, \lambda\_{i,0})\\\_{i=1}^I
    \sim H\\ once.

2.  For any scale \\c\\, form \\\lambda_i(c) = c \cdot \lambda\_{i,0}\\
    and compute the empirical approximation to population reliability
    \\\hat\rho_M(c)\\ from the test information function.

3.  Solve the scalar equation \\\hat\rho_M(c^\*) = \rho^\*\\ using
    deterministic root-finding (Brent's method via
    [`uniroot()`](https://rdrr.io/r/stats/uniroot.html)).

## Usage

``` r
eqc_calibrate(
  target_rho,
  n_items,
  model = c("rasch", "2pl"),
  latent_shape = "normal",
  item_source = "parametric",
  latent_params = list(),
  item_params = list(),
  reliability_metric = c("info", "tilde", "msem", "bar"),
  M = 10000L,
  c_bounds = c(0.3, 3),
  tol = 1e-04,
  seed = NULL,
  verbose = FALSE
)

# S3 method for class 'eqc_result'
print(x, digits = 4, ...)

# S3 method for class 'eqc_result'
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

  Character. Reliability definition used inside EQC:

  `"info"`

  :   Average-information reliability (default, recommended for EQC).
      Targets \\\tilde{\rho}\\, which is guaranteed to be monotone in
      \\c\\.

  `"msem"`

  :   MSEM-based marginal reliability (theoretically exact, but may have
      a non-monotone objective for EQC; see Details).

  Synonyms: `"tilde"` for `"info"`, `"bar"` for `"msem"`.

- M:

  Integer. Size of the empirical quadrature sample (default: 10000).

- c_bounds:

  Numeric length-2 vector. Search bounds for \\c\\. Default: c(0.3, 3).

- tol:

  Numeric. Tolerance for
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html). Default: 1e-4.

- seed:

  Optional integer for reproducibility.

- verbose:

  Logical. If TRUE, print progress messages.

- x:

  An object of class `"eqc_result"`.

- digits:

  Integer. Number of decimal places for printing.

- ...:

  Additional arguments passed to or from other methods.

- object:

  An object of class `"eqc_result"`.

## Value

An object of class `"eqc_result"` (a list) with elements:

- `c_star`:

  Calibrated discrimination scale \\c^\*\\.

- `target_rho`:

  Target reliability \\\rho^\*\\.

- `achieved_rho`:

  Empirical quadrature estimate \\\hat\rho_M(c^\*)\\.

- `metric`:

  Reliability metric used.

- `theta_quad`:

  Length-M vector of quadrature abilities.

- `theta_var`:

  Sample variance of theta_quad.

- `items_base`:

  item_params object with scale = 1 (baseline).

- `items_calib`:

  item_params object with discriminations scaled by c_star.

The input object, invisibly.

An object of class `"summary.eqc_result"` containing key calibration
results.

## Details

### Reliability Metrics

The function supports two reliability definitions:

- **Average-information** (`"info"`/`"tilde"`, **default**): Uses the
  arithmetic mean, \\\tilde{\rho}(c) = \sigma^2\_\theta
  \bar{\mathcal{J}}(c) / (\sigma^2\_\theta \bar{\mathcal{J}}(c) + 1)\\.
  By Jensen's inequality, \\\tilde{\rho} \geq \bar{w}\\, so this metric
  typically yields higher reliability values. **This is the recommended
  default for EQC** because the objective function \\\tilde{\rho}(c) -
  \rho^\*\\ is guaranteed to be monotone in \\c\\, ensuring that
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) will find the
  unique root.

- **MSEM-based** (`"msem"`/`"bar"`): Uses the harmonic mean of test
  information, \\\bar{w}(c) = \sigma^2\_\theta / (\sigma^2\_\theta +
  E\[1/\mathcal{J}(\theta;c)\])\\. This is theoretically exact but the
  objective \\\bar{w}(c) - \rho^\*\\ may be non-monotone (see Lee, 2025,
  Section 4.3), which can cause
  [`uniroot()`](https://rdrr.io/r/stats/uniroot.html) to fail or find an
  incorrect root. Use
  [`sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  if you need to target \\\bar{w}\\ directly.

### WLE vs EAP Reliability Interpretation

When validating with TAM, note that EAP reliability is systematically
higher than WLE reliability. This is not a bug but a mathematical
property of TAM's definitions. EAP reliability more directly corresponds
to the MSEM-based population reliability targeted by EQC. For
conservative inference, treat WLE as a lower bound and EAP as an upper
bound for true measurement precision.

## See also

[`sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
for the stochastic approximation alternative,
[`compute_rho_bar`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
and
[`compute_rho_tilde`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
for reliability computation utilities,
[`compute_reliability_tam`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)
for TAM validation.

## Examples

``` r
# Basic EQC calibration with parametric items (fast)
# \donttest{
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "parametric",
  M = 5000L,
  seed = 42
)
print(eqc_result)
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
#> 
# }

if (FALSE) { # \dontrun{
# EQC with IRW difficulties (requires irw package)
eqc_result2 <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  item_source = "irw",
  seed = 42,
  verbose = TRUE
)
} # }
```
