# Algorithm 1: Empirical Quadrature Calibration (EQC)

``` r
library(IRTsimrel)
```

## Overview

Empirical Quadrature Calibration (EQC) is the **primary algorithm** in
IRTsimrel for reliability-targeted simulation. Given a target marginal
reliability $`\rho^*`$, EQC finds a global discrimination scaling factor
$`c^*`$ such that the population reliability equals the target.

### Why EQC is Recommended

| Feature | EQC Advantage |
|----|----|
| **Speed** | Deterministic root-finding with $`O(\log(1/\varepsilon))`$ convergence |
| **Accuracy** | Typically achieves target within ±0.01 |
| **Simplicity** | No step-size tuning or learning rates |
| **Reliability** | Guaranteed convergence via Brent’s method |

For most simulation studies, EQC is all you need. Use SPC (Algorithm 2)
only for validation or specialized applications.

## Mathematical Foundation

### The Inverse Problem

**Forward problem**: Given item parameters and a scaling factor $`c`$,
compute reliability $`\rho(c)`$.

**Inverse problem**: Given target reliability $`\rho^*`$, find $`c^*`$
such that $`\rho(c^*) = \rho^*`$.

### Discrimination Scaling

EQC solves the inverse problem by scaling all item discriminations by a
common factor:

``` math
\lambda_i(c) = c \cdot \lambda_i^{(0)}
```

where $`\lambda_i^{(0)}`$ are baseline discriminations from the item
generator.

**Why this works**:

1.  **Separates informativeness from structure**: Changing $`c`$
    modifies reliability without altering item difficulty distributions
    or latent trait shape
2.  **Monotonic relationship**: Reliability is strictly increasing in
    $`c`$
3.  **Reduces dimensionality**: Multi-dimensional design space becomes a
    1D root-finding problem

### Monotonicity Guarantee

**Key Property**: Under mild regularity conditions, $`\rho(c)`$ is
strictly monotonically increasing in $`c`$:

``` math
\frac{\partial \rho(c)}{\partial c} > 0 \quad \text{for all } c > 0
```

This guarantees:

- As $`c \to 0`$: $`\rho(c) \to 0`$
- As $`c \to \infty`$: $`\rho(c) \to 1`$
- **Existence and uniqueness**: For any $`\rho^* \in (0, 1)`$, there
  exists a unique $`c^*`$

## Algorithm Steps

EQC operates in three steps:

### Step 1: Generate Quadrature Samples

Draw large fixed samples from the latent and item parameter
distributions:

- $`\{\theta_m\}_{m=1}^M \sim G`$ (latent abilities)
- $`\{(\beta_i, \lambda_i^{(0)})\}_{i=1}^I \sim H`$ (item parameters)

These samples serve as an **empirical quadrature rule** for
approximating population expectations.

### Step 2: Define the Reliability Function

For any scale $`c`$, compute the empirical reliability:

``` math
\hat{\rho}_M(c) = \frac{1}{M} \sum_{m=1}^{M} \rho(c; \theta_m, \boldsymbol{\beta}, \boldsymbol{\lambda}(c))
```

This approximates the population reliability $`\mathbb{E}[\rho(c)]`$ by
the Law of Large Numbers.

### Step 3: Root-Finding

Solve the scalar equation:

``` math
g(c) = \hat{\rho}_M(c) - \rho^* = 0
```

using Brent’s method
([`uniroot()`](https://rdrr.io/r/stats/uniroot.html) in R), which
combines bisection, secant, and inverse quadratic interpolation for
robust convergence.

## Basic Usage

``` r
# Basic EQC calibration
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  seed = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

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

## Understanding the Output

The `eqc_result` object contains:

| Component       | Description                                            |
|-----------------|--------------------------------------------------------|
| `c_star`        | Calibrated scaling factor $`c^*`$                      |
| `target_rho`    | Target reliability $`\rho^*`$                          |
| `achieved_rho`  | Achieved reliability $`\hat{\rho}(c^*)`$               |
| `metric`        | Reliability metric used (`msem` or `info`)             |
| `theta_quad`    | Quadrature sample of abilities                         |
| `beta_vec`      | Item difficulties                                      |
| `lambda_base`   | Baseline discriminations                               |
| `lambda_scaled` | Calibrated discriminations $`c^* \cdot \lambda^{(0)}`$ |
| `items_base`    | Item parameters before calibration                     |
| `items_calib`   | Item parameters after calibration                      |

### Accessing Calibrated Parameters

``` r
# The calibrated scaling factor
cat(sprintf("c* = %.4f\n", eqc_result$c_star))
#> c* = 0.9082

# Item difficulties
head(eqc_result$beta_vec)
#> [1]  0.1312995 -1.2945477  0.3836801 -0.6358214  0.5071690  0.4314725

# Calibrated discriminations
head(eqc_result$lambda_scaled)
#> [1] 0.9082186 0.9082186 0.9082186 0.9082186 0.9082186 0.9082186

# Full calibrated item parameters
head(eqc_result$items_calib$data)
#>   form_id item_id       beta    lambda lambda_unscaled
#> 1       1       1  0.1312995 0.9082186               1
#> 2       1       2 -1.2945477 0.9082186               1
#> 3       1       3  0.3836801 0.9082186               1
#> 4       1       4 -0.6358214 0.9082186               1
#> 5       1       5  0.5071690 0.9082186               1
#> 6       1       6  0.4314725 0.9082186               1
```

### Convergence Information

``` r
# Root-finding status
cat(sprintf("Root status: %s\n", eqc_result$misc$root_status))
#> Root status: uniroot_success

# Achievable reliability range
cat(sprintf("Reliability at c_lower (%.2f): %.4f\n", 
            eqc_result$misc$c_bounds[1], 
            eqc_result$misc$rho_bounds["rho_L"]))
#> Reliability at c_lower (0.30): 0.3555
cat(sprintf("Reliability at c_upper (%.2f): %.4f\n", 
            eqc_result$misc$c_bounds[2], 
            eqc_result$misc$rho_bounds["rho_U"]))
#> Reliability at c_upper (3.00): 0.8240
```

## Key Parameters

### Target Reliability (`target_rho`)

The desired marginal reliability, typically between 0.50 and 0.90:

``` r
# Low reliability (formative assessment)
eqc_low <- eqc_calibrate(target_rho = 0.60, n_items = 20, seed = 1)

# Moderate reliability (typical research)
eqc_mod <- eqc_calibrate(target_rho = 0.75, n_items = 20, seed = 1)

# High reliability (standardized test)
eqc_high <- eqc_calibrate(target_rho = 0.85, n_items = 20, seed = 1)
```

### Quadrature Size (`M`)

The number of latent ability samples for empirical integration. Larger
$`M`$ reduces Monte Carlo error:

``` r
# Default: M = 10,000
eqc_default <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, M = 10000, seed = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

# Higher precision: M = 50,000
eqc_precise <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, M = 50000, seed = 42
)

cat(sprintf("M = 10,000: c* = %.6f\n", eqc_default$c_star))
#> M = 10,000: c* = 0.908219
cat(sprintf("M = 50,000: c* = %.6f\n", eqc_precise$c_star))
#> M = 50,000: c* = 1.156505
```

**Recommendation**: $`M \geq 10,000`$ for routine use; $`M \geq 50,000`$
for high-precision work.

### Search Bounds (`c_bounds`)

The interval $`[c_{\min}, c_{\max}]`$ for root-finding:

``` r
# Default bounds
eqc_default <- eqc_calibrate(
  target_rho = 0.80, n_items = 25,
  c_bounds = c(0.3, 3),  # Default
  seed = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

# Extended bounds for high reliability
eqc_extended <- eqc_calibrate(
  target_rho = 0.85, n_items = 25,
  c_bounds = c(0.1, 5),
  seed = 42
)
#> Warning in eqc_calibrate(target_rho = 0.85, n_items = 25, c_bounds = c(0.1, : Target rho* = 0.850 exceeds the maximum achievable reliability for this configuration.
#>   - Items: 25
#>   - At c = 5.00 (upper bound): rho = 0.0440
#>   - Gap: 0.8060
#> 
#> Suggestions:
#>   1. Use reliability_metric = 'info' (typically yields higher values)
#>   2. Increase c_bounds[2] beyond 5.0
#>   3. Increase n_items for higher achievable reliability
#>   4. Accept achieved rho = 0.0440 as maximum for this design
#> 
#> Returning c_star = 5.00 (upper bound).
```

**When to extend bounds**:

- Target reliability near achievable maximum → increase `c_bounds[2]`
- Very low target reliability → decrease `c_bounds[1]`

## Reliability Metrics

EQC supports two reliability definitions:

### MSEM-Based Reliability ($`\bar{w}`$)

The theoretically exact marginal reliability based on the **harmonic
mean** of test information:

``` math
\bar{w}(c) = \frac{\sigma^2_\theta}{\sigma^2_\theta + \mathbb{E}[1/\mathcal{J}(\theta; c)]}
```

``` r
eqc_msem <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  reliability_metric = "msem",  # Default, or use "bar"
  seed = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.
```

### Average-Information Reliability ($`\tilde{\rho}`$)

Based on the **arithmetic mean** of test information:

``` math
\tilde{\rho}(c) = \frac{\sigma^2_\theta \bar{\mathcal{J}}(c)}{\sigma^2_\theta \bar{\mathcal{J}}(c) + 1}
```

``` r
eqc_info <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  reliability_metric = "info",  # Or use "tilde"
  seed = 42
)
```

### Key Relationship

By Jensen’s inequality, $`\tilde{\rho} \geq \bar{w}`$ always holds:

``` r
cat(sprintf("MSEM metric c*: %.4f\n", eqc_msem$c_star))
#> MSEM metric c*: 0.9082
cat(sprintf("Info metric c*: %.4f\n", eqc_info$c_star))
#> Info metric c*: 0.8830
```

The `info` metric typically requires a **smaller** $`c^*`$ to achieve
the same target, because it yields higher reliability values for a given
$`c`$.

**Recommendation**:

- Use `"msem"` (default) for theoretical correctness
- Use `"info"` when targeting very high reliability that `"msem"` cannot
  achieve

## Working with Different Models

### Rasch Model

``` r
eqc_rasch <- eqc_calibrate(
  target_rho = 0.75,
  n_items = 30,
  model = "rasch",
  seed = 42
)

# All baseline discriminations are 1
unique(eqc_rasch$lambda_base)
#> [1] 1
```

### 2PL Model

``` r
eqc_2pl <- eqc_calibrate(
  target_rho = 0.75,
  n_items = 30,
  model = "2pl",
  item_source = "irw",
  item_params = list(
    discrimination_params = list(
      mu_log = 0,
      sigma_log = 0.3,
      rho = -0.3
    )
  ),
  seed = 42
)

# Discriminations vary across items
summary(eqc_2pl$lambda_base)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.6887  0.7896  0.9104  1.0097  1.2525  1.6613
```

## Working with Different Latent Distributions

``` r
# Normal distribution
eqc_normal <- eqc_calibrate(
  target_rho = 0.80, n_items = 25,
  latent_shape = "normal",
  seed = 42
)

# Bimodal distribution
eqc_bimodal <- eqc_calibrate(
  target_rho = 0.80, n_items = 25,
  latent_shape = "bimodal",
  latent_params = list(delta = 0.8),
  seed = 42
)

# Heavy-tailed distribution
eqc_heavy <- eqc_calibrate(
  target_rho = 0.80, n_items = 25,
  latent_shape = "heavy_tail",
  latent_params = list(df = 5),
  seed = 42
)
```

Different latent shapes require different $`c^*`$ values to achieve the
same target reliability:

``` r
cat(sprintf("Normal:     c* = %.4f\n", eqc_normal$c_star))
cat(sprintf("Bimodal:    c* = %.4f\n", eqc_bimodal$c_star))
cat(sprintf("Heavy-tail: c* = %.4f\n", eqc_heavy$c_star))
```

## Generating Response Data

After calibration, generate item response data:

``` r
# Generate response data
sim_data <- simulate_response_data(
  eqc_result = eqc_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Response matrix dimensions
dim(sim_data$response_matrix)
#> [1] 1000   25

# First few responses
head(sim_data$response_matrix[, 1:6])
#>      item1 item2 item3 item4 item5 item6
#> [1,]     0     1     0     1     0     1
#> [2,]     0     0     1     0     0     0
#> [3,]     1     1     1     1     1     1
#> [4,]     1     1     1     1     0     1
#> [5,]     0     1     0     1     0     0
#> [6,]     1     1     1     1     1     0

# True abilities
head(sim_data$theta)
#> [1] -0.56047565 -0.23017749  1.55870831  0.07050839  0.12928774  1.71506499
```

## Validation with TAM

Validate achieved reliability using the TAM package:

``` r
# Compute WLE and EAP reliability
tam_rel <- compute_reliability_tam(
  resp = sim_data$response_matrix,
  model = "rasch"
)

cat(sprintf("Target reliability:  %.4f\n", eqc_result$target_rho))
cat(sprintf("EQC achieved rho:    %.4f\n", eqc_result$achieved_rho))
cat(sprintf("TAM WLE reliability: %.4f\n", tam_rel$rel_wle))
cat(sprintf("TAM EAP reliability: %.4f\n", tam_rel$rel_eap))
```

### Understanding WLE vs EAP Reliability

An important validation finding: TAM’s EAP reliability is systematically
**higher** than WLE reliability. This is not a bug but reflects
different mathematical definitions:

| Metric | Definition | Typical Pattern |
|----|----|----|
| **WLE** | Design-effect based: $`1 - \bar{s}^2 / V_{WLE}`$ | Conservative (lower bound) |
| **EAP** | Posterior variance based: $`V_{EAP} / (V_{EAP} + \bar{\sigma}^2)`$ | Liberal (upper bound) |

**Practical interpretation**:

- WLE reliability ≈ lower bound for true measurement precision
- EAP reliability ≈ upper bound for true measurement precision
- EQC’s theoretical reliability falls between these bounds

## Handling Boundary Cases

### Target Too High

When the target reliability exceeds what’s achievable with the given
configuration:

``` r
# Try to achieve very high reliability with few items
eqc_high <- eqc_calibrate(
  target_rho = 0.95,
  n_items = 15,
  model = "rasch",
  c_bounds = c(0.3, 3),
  seed = 42
)
#> Warning in eqc_calibrate(target_rho = 0.95, n_items = 15, model = "rasch", : Target rho* = 0.950 exceeds the maximum achievable reliability for this configuration.
#>   - Items: 15
#>   - At c = 3.00 (upper bound): rho = 0.7054
#>   - Gap: 0.2446
#> 
#> Suggestions:
#>   1. Use reliability_metric = 'info' (typically yields higher values)
#>   2. Increase c_bounds[2] beyond 3.0
#>   3. Increase n_items for higher achievable reliability
#>   4. Accept achieved rho = 0.7054 as maximum for this design
#> 
#> Returning c_star = 3.00 (upper bound).
```

The function returns `c_star = c_bounds[2]` (upper bound) with a
warning.

**Solutions**:

1.  Use `reliability_metric = "info"` (yields higher values)
2.  Increase `c_bounds[2]`
3.  Increase `n_items`
4.  Accept the achievable maximum

### Target Too Low

When the target is below the minimum achievable:

``` r
# Very low target with high discrimination items
eqc_low <- eqc_calibrate(
  target_rho = 0.30,
  n_items = 50,
  model = "2pl",
  c_bounds = c(0.3, 3),
  seed = 42
)
```

The function returns `c_star = c_bounds[1]` (lower bound) with a
warning.

**Solution**: Decrease `c_bounds[1]`.

## Verbose Mode

Enable verbose output to see algorithm progress:

``` r
eqc_verbose <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42,
  verbose = TRUE
)
#> Step 1: Generating quadrature samples...
#>   M (quad persons) = 10000
#>   I (items)        = 25
#>   theta: mean = -0.011, sd = 1.006, var = 1.012
#>   beta:  mean = 0.000, sd = 0.670
#>   lambda_base: mean = 1.000, sd = 0.000
#>   metric = msem
#> Step 2: Running root-finding algorithm...
#>   At c = 0.300: rho = 0.3555, g = -0.4445
#>   At c = 3.000: rho = 0.8240, g = 0.0240
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.
#>   c* = 0.908219
#>   Target rho    = 0.8000
#>   Achieved rho  = 0.8000
#>   Root status   = uniroot_success
```

## Validation with SPC

For rigorous validation, compare EQC results with the SPC algorithm:

``` r
# Run EQC
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
  c_init = eqc_result,  # Use EQC for warm start
  n_iter = 200,
  seed = 42
)

# Compare
compare_eqc_spc(eqc_result, spc_result)
```

When EQC and SPC agree (typically within 5%), you can be confident in
your calibration.

## Complete Workflow Example

``` r
# ==== 1. Define simulation parameters ====
target_rho <- 0.80
n_items <- 25
n_persons <- 1000
n_replications <- 100

# ==== 2. Calibrate once using EQC ====
eqc_result <- eqc_calibrate(
  target_rho = target_rho,
  n_items = n_items,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  M = 20000,
  seed = 42,
  verbose = TRUE
)

cat(sprintf("\nCalibrated c* = %.4f for target rho* = %.2f\n\n",
            eqc_result$c_star, target_rho))

# ==== 3. Run simulation replications ====
results <- vector("list", n_replications)

for (r in 1:n_replications) {
  # Generate data
  sim_data <- simulate_response_data(
    eqc_result = eqc_result,
    n_persons = n_persons,
    latent_shape = "normal",
    seed = r
  )
  
  # Your analysis here...
  # results[[r]] <- your_analysis(sim_data)
  
  if (r %% 10 == 0) cat(sprintf("Completed replication %d\n", r))
}

# ==== 4. Summarize results ====
# summary_stats <- summarize_results(results)
```

## Troubleshooting Guide

| Problem | Possible Cause | Solution |
|----|----|----|
| Warning: upper bound hit | Target too high | Use `metric = "info"`, increase `c_bounds[2]`, or increase `n_items` |
| Warning: lower bound hit | Target too low | Decrease `c_bounds[1]` |
| High Monte Carlo error | Small `M` | Increase `M` to 20,000+ |
| WLE reliability too low | Normal TAM behavior | EAP provides upper bound; WLE is conservative |
| Different results with same seed | Item parameter resampling | Fix item parameters or use same seed consistently |

## Theoretical Properties

### Convergence Rate

EQC uses Brent’s method, which achieves superlinear convergence:

``` math
|c_n - c^*| = O\left(\frac{1}{2^{2^n}}\right)
```

In practice, convergence to machine precision typically requires \< 50
function evaluations.

### Monte Carlo Error

The empirical reliability estimate has error:

``` math
|\hat{\rho}_M(c) - \rho(c)| = O\left(\frac{1}{\sqrt{M}}\right)
```

With $`M = 10,000`$, the Monte Carlo standard error is approximately
0.01.

### Consistency

As $`M \to \infty`$, the EQC solution converges to the true population
solution:

``` math
\hat{c}^*_M \xrightarrow{\text{a.s.}} c^* \quad \text{as } M \to \infty
```

## Summary

EQC is a fast, accurate, and reliable algorithm for reliability-targeted
IRT simulation:

- **When to use**: Routine simulation studies, initial calibration
- **Key parameters**: `target_rho`, `n_items`, `M`, `reliability_metric`
- **Validation**: Compare with TAM reliability and/or SPC algorithm
- **Typical accuracy**: Within ±0.01 of target reliability

For specialized applications or independent validation, see the SPC
algorithm vignette.

## References

Lee, J.-H. (2025). Reliability-targeted simulation for item response
data. *Manuscript in preparation*.

Brent, R. P. (1973). *Algorithms for Minimization Without Derivatives*.
Prentice-Hall.
