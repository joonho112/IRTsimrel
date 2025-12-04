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
  item_source = "irw",
  latent_params = list(),
  item_params = list(),
  reliability_metric = c("msem", "info", "bar", "tilde"),
  M = 10000L,
  c_bounds = c(0.3, 3),
  tol = 1e-04,
  seed = NULL,
  verbose = FALSE
)
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
  (e.g. `"irw"`, `"parametric"`, `"hierarchical"`, `"custom"`).

- latent_params:

  List. Additional arguments passed to
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- item_params:

  List. Additional arguments passed to
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- reliability_metric:

  Character. Reliability definition used inside EQC:

  `"msem"`

  :   MSEM-based marginal reliability (default, theoretically exact).

  `"info"`

  :   Average-information reliability (faster, more stable).

  Synonyms: `"bar"` for `"msem"`, `"tilde"` for `"info"`.

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

## Details

### Reliability Metrics

The function supports two reliability definitions:

- **MSEM-based** (`"msem"`/`"bar"`): Uses the harmonic mean of test
  information, \\\bar{w}(c) = \sigma^2\_\theta / (\sigma^2\_\theta +
  E\[1/\mathcal{J}(\theta;c)\])\\. This is theoretically exact but may
  have a lower ceiling for high reliability.

- **Average-information** (`"info"`/`"tilde"`): Uses the arithmetic
  mean, \\\tilde{\rho}(c) = \sigma^2\_\theta \bar{\mathcal{J}}(c) /
  (\sigma^2\_\theta \bar{\mathcal{J}}(c) + 1)\\. By Jensen's inequality,
  \\\tilde{\rho} \geq \bar{w}\\, so this metric typically yields higher
  reliability values.

### WLE vs EAP Reliability Interpretation

When validating with TAM, note that EAP reliability is systematically
higher than WLE reliability. This is not a bug but a mathematical
property of TAM's definitions. EAP reliability more directly corresponds
to the MSEM-based population reliability targeted by EQC. For
conservative inference, treat WLE as a lower bound and EAP as an upper
bound for true measurement precision.

## See also

[`spc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
for the stochastic approximation alternative,
[`compute_rho_bar`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
and
[`compute_rho_tilde`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
for reliability computation utilities,
[`compute_reliability_tam`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)
for TAM validation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic EQC calibration
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  seed = 42,
  verbose = TRUE
)
print(eqc_result)
} # }
```
