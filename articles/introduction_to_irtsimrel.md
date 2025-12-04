# Introduction to IRTsimrel

``` r
library(IRTsimrel)
```

## Overview

The **IRTsimrel** package provides a principled framework for
*reliability-targeted simulation* of Item Response Theory (IRT) data.
Instead of treating reliability as an implicit outcome of simulation
design choices, IRTsimrel allows researchers to specify a **target
reliability level** as an explicit input parameter.

### Why Reliability-Targeted Simulation?

In Monte Carlo simulation studies for IRT, researchers routinely vary:

- Sample size ($`N`$)
- Test length ($`I`$, number of items)
- Item parameter distributions ($`\boldsymbol{\lambda}`$,
  $`\boldsymbol{\beta}`$)
- Latent trait distribution shape ($`G`$)

However, **marginal reliability**—the fundamental measure of data
informativeness—is almost never directly controlled or reported. This
creates several problems:

1.  **Ecological validity threat**: Real-world assessments often have
    reliabilities of 0.5–0.7, but simulation studies may implicitly
    assume higher values.

2.  **Confounded model comparisons**: Conclusions about model
    superiority may only hold within certain reliability regimes.

3.  **Limited replicability**: Without knowing the implied reliability,
    exact replication is impossible.

### The Multilevel Modeling Analogy

In multilevel/hierarchical linear modeling, the **intraclass correlation
(ICC)** is always a primary design factor in simulation studies. ICC
determines the signal-to-noise ratio—the proportion of variance at the
cluster level versus residual level.

**Marginal reliability in IRT serves the same conceptual role as ICC in
multilevel models.** Both quantify the proportion of variance
attributable to the latent construct versus measurement error.

IRTsimrel brings IRT simulation methodology into alignment with this
best practice.

## The Key Idea: Separation of Structure and Scale

IRTsimrel implements a fundamental principle: **separate structure from
scale**.

- **Structure**: Realistic item characteristics (difficulty
  distributions, discrimination patterns, difficulty-discrimination
  correlations) and flexible latent trait distributions come from
  empirically-grounded generators.

- **Scale**: A global discrimination scaling factor $`c`$ is calibrated
  to achieve the target reliability.

The calibrated item discriminations are:

``` math
\lambda_i^* = c^* \cdot \lambda_{i,0}
```

where $`\lambda_{i,0}`$ are the baseline discriminations from the
generator and $`c^*`$ is the calibrated scaling factor.

## Two Calibration Algorithms

IRTsimrel provides two complementary algorithms:

### Algorithm 1: Empirical Quadrature Calibration (EQC)

[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
is the **recommended method for routine use**. It:

1.  Draws a large fixed “quadrature” sample from the latent and item
    parameter distributions
2.  Computes an empirical approximation to population reliability as a
    function of $`c`$
3.  Solves for $`c^*`$ using deterministic root-finding (Brent’s method)

**Advantages**: Fast, deterministic, highly accurate (typically within
±0.01 of target).

### Algorithm 2: Stochastic Approximation Calibration (SPC)

[`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)
uses the Robbins-Monro stochastic approximation algorithm. It’s useful
for:

- Independent validation of EQC results
- Complex data-generating processes
- Exact MSEM-based reliability targeting

**Advantages**: More general, targets exact marginal reliability,
applicable to complex DGPs.

## Quick Start Example

Let’s walk through a basic example. Suppose we want to generate Rasch
model data with:

- 25 items
- Target reliability of 0.80
- Normal latent trait distribution
- Realistic item difficulties from the Item Response Warehouse (IRW)

### Step 1: Calibrate Using EQC

``` r
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
```

    #> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

``` r
# View the results
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
#>   Absolute error               : 2.27e-06
#>   Scaling factor (c*)          : 0.9082
#> 
#> Design Parameters:
#>   Number of items (I)          : 25
#>   Quadrature points (M)        : 10000
#>   Reliability metric           : MSEM-based (bar/w)
#>   Latent variance              : 1.0123
#> 
#> Convergence:
#>   Root status                  : uniroot_success
#>   Search bracket               : [0.300, 3.000]
#>   Bracket reliabilities        : [0.3555, 0.8240]
#> 
#> Parameter Summaries:
#>   theta:        mean = -0.011, sd = 1.006
#>   beta:         mean = 0.000, sd = 0.670, range = [-1.48, 1.06]
#>   lambda_base:  mean = 1.000, sd = 0.000
#>   lambda_scaled: mean = 0.908, sd = 0.000
```

The output shows:

- **c\***: The calibrated scaling factor
- **Achieved reliability**: Very close to our target of 0.80
- **Absolute error**: Typically \< 0.001

### Step 2: Generate Response Data

Now use the calibrated parameters to generate item response data:

``` r
# Generate response data
sim_data <- simulate_response_data(
  eqc_result = eqc_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Examine the response matrix
dim(sim_data$response_matrix)
#> [1] 1000   25
head(sim_data$response_matrix[, 1:6])
#>      item1 item2 item3 item4 item5 item6
#> [1,]     0     1     0     1     0     1
#> [2,]     0     0     1     0     0     0
#> [3,]     1     1     1     1     1     1
#> [4,]     1     1     1     1     0     1
#> [5,]     0     1     0     1     0     0
#> [6,]     1     1     1     1     1     0
```

### Step 3: Validate with TAM (Optional)

For rigorous validation, you can estimate reliability using the TAM
package:

``` r
# Compute WLE and EAP reliability
tam_rel <- compute_reliability_tam(
  resp = sim_data$response_matrix,
  model = "rasch"
)

cat(sprintf("Target reliability: %.3f\n", eqc_result$target_rho))
cat(sprintf("WLE reliability:    %.3f\n", tam_rel$rel_wle))
cat(sprintf("EAP reliability:    %.3f\n", tam_rel$rel_eap))
```

**Note on WLE vs EAP Reliability**: TAM’s EAP reliability is typically
higher than WLE reliability. This is not a bug but reflects different
mathematical definitions. For conservative inference, treat WLE as a
lower bound and EAP as an upper bound for true measurement precision.

## Understanding Reliability Metrics

IRTsimrel supports two reliability definitions: \### MSEM-Based
Reliability ($`\bar{w}`$)

The theoretically exact marginal reliability based on the harmonic mean
of test information:

``` math
\bar{w}(c) = \frac{\sigma^2_\theta}{\sigma^2_\theta + \mathbb{E}[1/\mathcal{J}(\theta; c)]}
```

where $`\mathcal{J}(\theta; c)`$ is the test information function at
ability $`\theta`$.

### Average-Information Reliability ($`\tilde{\rho}`$)

Based on the arithmetic mean of test information:

``` math
\tilde{\rho}(c) = \frac{\sigma^2_\theta \bar{\mathcal{J}}(c)}{\sigma^2_\theta \bar{\mathcal{J}}(c) + 1}
```

where $`\bar{\mathcal{J}}(c) = \mathbb{E}[\mathcal{J}(\theta; c)]`$.

**Key relationship**: By Jensen’s inequality,
$`\tilde{\rho} \geq \bar{w}`$ always holds.

You can choose between these metrics:

``` r
# MSEM-based (default, theoretically exact)
eqc_msem <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  reliability_metric = "msem",  # or "bar"
  seed = 42
)

# Average-information (faster, higher values)
eqc_info <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  reliability_metric = "info",  # or "tilde"
  seed = 42
)
```

## Working with Different Latent Distributions

The
[`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
function generates latent abilities with various distributional shapes,
all pre-standardized to have mean 0 and variance 1:

``` r
# Compare different shapes
compare_shapes(
  n = 3000,
  shapes = c("normal", "bimodal", "trimodal", 
             "skew_pos", "heavy_tail", "uniform"),
  seed = 42
)
```

### Using Non-Normal Latent Distributions

You can specify any shape for calibration:

``` r
# Calibrate for a bimodal population
eqc_bimodal <- eqc_calibrate(
  target_rho = 0.75,
  n_items = 20,
  model = "rasch",
  latent_shape = "bimodal",
  latent_params = list(delta = 0.8),  # Mode separation
  seed = 42
)
```

Available shapes include:

| Shape        | Description                              |
|--------------|------------------------------------------|
| `normal`     | Standard normal N(0,1)                   |
| `bimodal`    | Symmetric two-component Gaussian mixture |
| `trimodal`   | Symmetric three-component mixture        |
| `skew_pos`   | Right-skewed (standardized Gamma)        |
| `skew_neg`   | Left-skewed (negated Gamma)              |
| `heavy_tail` | Heavy-tailed (standardized Student-t)    |
| `uniform`    | Uniform distribution                     |
| `floor`      | Floor effect                             |
| `ceiling`    | Ceiling effect                           |
| `custom`     | User-specified mixture                   |

## Working with the 2PL Model

For the two-parameter logistic (2PL) model, discriminations vary across
items:

``` r
# 2PL calibration with IRW-based difficulties and copula method
eqc_2pl <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 30,
  model = "2pl",
  latent_shape = "normal",
  item_source = "irw",
  item_params = list(
    discrimination_params = list(
      mu_log = 0,        # Mean of log-discrimination
      sigma_log = 0.3,   # SD of log-discrimination
      rho = -0.3         # Correlation between difficulty and log-discrimination
    )
  ),
  seed = 42
)
```

The copula method for generating item parameters preserves:

- Realistic difficulty distributions from the IRW
- Log-normal marginal for discriminations
- The empirically-observed negative correlation between difficulty and
  discrimination

## Validating with the SPC Algorithm

For rigorous validation, use SPC to independently verify EQC results:

``` r
# First run EQC
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42
)

# Validate with SPC (warm start from EQC)
spc_result <- spc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  c_init = eqc_result,  # Use EQC result for warm start
  n_iter = 200,
  seed = 42
)

# Compare results
compare_eqc_spc(eqc_result, spc_result)
```

When results agree (typically within 5%), you can be confident in your
calibration.

## Typical Workflow Summary

A complete reliability-targeted simulation workflow:

``` r
# 1. Define target reliability
target_rho <- 0.80

# 2. Calibrate using EQC
eqc_result <- eqc_calibrate(
  target_rho = target_rho,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  seed = 42
)

# 3. (Optional) Validate with SPC
spc_result <- spc_calibrate(
  target_rho = target_rho,
  n_items = 25,
  model = "rasch",
  c_init = eqc_result,
  n_iter = 200,
  seed = 42
)

# 4. Generate response data
sim_data <- simulate_response_data(
  eqc_result = eqc_result,
  n_persons = 1000,
  seed = 123
)

# 5. (Optional) Validate achieved reliability with TAM
tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")

# 6. Use the data for your simulation study
# ... your analysis code here ...
```

## Practical Recommendations

### Choosing the Right Algorithm

| Scenario                     | Recommended Algorithm                       |
|------------------------------|---------------------------------------------|
| Routine simulation work      | EQC                                         |
| Independent validation       | SPC with EQC warm start                     |
| Complex custom DGPs          | SPC                                         |
| Very high target reliability | Consider increasing `c_bounds` or `n_items` |

### Setting Target Reliability

- **Low reliability (0.50–0.65)**: Typical for short formative
  assessments
- **Moderate reliability (0.70–0.80)**: Common in educational research
- **High reliability (0.85+)**: Standardized tests; may require many
  items or high discriminations

### When Target Reliability Cannot Be Achieved

If EQC reports hitting the upper bound, consider:

1.  Using `reliability_metric = "info"` (yields higher values)
2.  Increasing `c_bounds[2]` beyond the default of 3
3.  Increasing `n_items`
4.  Accepting the achievable maximum for your design

## Further Resources

For more detailed information, see:

- `vignette("latent-distributions")`: Working with different latent
  trait distributions
- `vignette("item-parameters")`: Generating realistic item parameters
- `vignette("spc-algorithm")`: Advanced usage of the SPC algorithm
- [`vignette("validation")`](https://joonho112.github.io/IRTsimrel/articles/validation.md):
  Comprehensive validation procedures

## References

Lee, J.-H. (2025). Reliability-targeted simulation for item response
data. *arXiv Preprint*.

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics, 22*(3), 400–407.

Zhang, L., Fellinghauer, C., Geerlings, H., & Sijtsma, K. (2025).
Realistic simulation of item difficulties using the Item Response
Warehouse. *PsyArXiv*. <https://doi.org/10.31234/osf.io/r5mxv>
