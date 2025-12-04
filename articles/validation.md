# Validating Calibration Results

``` r
library(IRTsimrel)
```

## Overview

After calibrating item parameters using EQC or SPC, it is essential to
**validate** that the achieved reliability matches your target.
IRTsimrel provides several validation tools:

| Function | Purpose |
|----|----|
| [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md) | Generate item response data from calibrated parameters |
| [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md) | Compute WLE/EAP reliability using the TAM package |
| [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md) | Compare EQC and SPC calibration results |

This vignette covers:

1.  Generating response data from calibrated parameters
2.  External validation using the TAM package
3.  Understanding WLE vs EAP reliability
4.  Comparing EQC and SPC results
5.  Monte Carlo validation workflow
6.  Troubleshooting common validation issues

## Generating Response Data

The
[`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
function generates item response matrices using calibrated parameters
from EQC.

### Basic Usage

``` r
# First, run EQC calibration
eqc_result <- eqc_calibrate(
target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  seed = 42
)
#> Note: Target rho* = 0.800 is near the achievable maximum (0.824) for this configuration.

# Generate response data
sim_data <- simulate_response_data(
  eqc_result = eqc_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Examine the output
names(sim_data)
#> [1] "response_matrix" "theta"           "beta"            "lambda"
```

### Output Structure

The function returns a list containing:

``` r
# Response matrix: N persons × I items
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

# True abilities (used for generating responses)
head(sim_data$theta)
#> [1] -0.56047565 -0.23017749  1.55870831  0.07050839  0.12928774  1.71506499
cat(sprintf("Theta mean: %.3f, SD: %.3f\n", 
            mean(sim_data$theta), sd(sim_data$theta)))
#> Theta mean: 0.016, SD: 0.992

# Item parameters used
cat(sprintf("Number of items: %d\n", length(sim_data$beta)))
#> Number of items: 25
cat(sprintf("Beta range: [%.2f, %.2f]\n", 
            min(sim_data$beta), max(sim_data$beta)))
#> Beta range: [-1.48, 1.06]
```

### Matching Latent Distribution

**Important**: Use the same latent distribution for simulation as you
used for calibration:

``` r
# Calibrate for bimodal population
eqc_bimodal <- eqc_calibrate(
  target_rho = 0.75,
  n_items = 20,
  latent_shape = "bimodal",
  latent_params = list(delta = 0.8),
  seed = 42
)

# Generate data with SAME distribution
sim_data <- simulate_response_data(
  eqc_result = eqc_bimodal,
  n_persons = 1000,
  latent_shape = "bimodal",           # Same shape!
  latent_params = list(delta = 0.8),  # Same params!
  seed = 123
)
```

If distributions don’t match, achieved reliability will differ from
target.

## External Validation with TAM

The
[`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)
function provides independent validation using the TAM package’s
official reliability functions.

### Requirements

``` r
# Install TAM if needed
install.packages("TAM")
```

### Computing Reliability

``` r
# Compute WLE and EAP reliability
tam_rel <- compute_reliability_tam(
  resp = sim_data$response_matrix,
  model = "rasch",
  verbose = FALSE
)

# View results
cat(sprintf("Target reliability:  %.4f\n", eqc_result$target_rho))
cat(sprintf("EQC achieved rho:    %.4f\n", eqc_result$achieved_rho))
cat(sprintf("TAM WLE reliability: %.4f\n", tam_rel$rel_wle))
cat(sprintf("TAM EAP reliability: %.4f\n", tam_rel$rel_eap))
```

### Output Components

``` r
names(tam_rel)

# rel_wle: WLE reliability coefficient
# rel_eap: EAP reliability coefficient
# mod: Fitted TAM model object
# wle: Output from TAM::tam.wle()
```

### Using with 2PL Model

``` r
# For 2PL data
eqc_2pl <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 30,
  model = "2pl",
  seed = 42
)

sim_2pl <- simulate_response_data(eqc_2pl, n_persons = 1000, seed = 123)

tam_2pl <- compute_reliability_tam(
  resp = sim_2pl$response_matrix,
  model = "2pl"  # Specify 2PL
)
```

## Understanding WLE vs EAP Reliability

A critical validation finding: **TAM’s EAP reliability is systematically
higher than WLE reliability**. This is not a bug—it reflects
fundamentally different mathematical definitions.

### The Two Definitions

**WLE Reliability** (design-effect based):
``` math
\rho_{\text{WLE}} = 1 - \frac{\bar{s}^2}{V_{\text{WLE}}}
```

where $`\bar{s}^2`$ is the average squared standard error and
$`V_{\text{WLE}}`$ is the variance of WLE estimates.

**EAP Reliability** (posterior variance based):
``` math
\rho_{\text{EAP}} = \frac{V_{\text{EAP}}}{V_{\text{EAP}} + \bar{\sigma}^2}
```

where $`V_{\text{EAP}}`$ is the variance of EAP estimates and
$`\bar{\sigma}^2`$ is the average posterior variance.

### Mathematical Relationship

Under TAM’s definitions,
$`\rho_{\text{EAP}} \geq \rho_{\text{WLE}}`$**always** holds. The
inequality is strict except in trivial cases.

### Practical Interpretation

| Metric          | Interpretation                     | Use Case            |
|-----------------|------------------------------------|---------------------|
| WLE reliability | Conservative lower bound           | Worst-case scenario |
| EAP reliability | Liberal upper bound                | Best-case scenario  |
| EQC/SPC target  | Theoretical population reliability | Design target       |

**Recommendation**: For conservative inference, focus on WLE
reliability. The EQC/SPC target typically falls between WLE and EAP.

### Expected Patterns

``` r
# Typical validation results
# Target: 0.80

# WLE reliability: 0.78 - 0.82  (slight underestimation common)
# EAP reliability: 0.82 - 0.86  (slight overestimation common)

# Both within ±0.03 of target is excellent
# WLE within ±0.05 of target is acceptable
```

## Comparing EQC and SPC Results

The
[`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
function provides a formal comparison between the two calibration
algorithms.

### Basic Comparison

``` r
# Run EQC
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42
)

# Run SPC with EQC warm start
spc_result <- spc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  c_init = eqc_result,
  n_iter = 200,
  seed = 42
)

# Compare
comparison <- compare_eqc_spc(eqc_result, spc_result)
```

### Interpreting the Comparison

``` r
# Output includes:
comparison$c_eqc      # EQC's c*
comparison$c_spc      # SPC's c*
comparison$diff_abs   # Absolute difference
comparison$diff_pct   # Percent difference
comparison$agreement  # TRUE if difference < 5%
```

### Expected Differences

EQC and SPC may produce slightly different results due to:

1.  **Different reliability formulas**: EQC uses fixed quadrature; SPC
    uses fresh samples
2.  **Stochastic noise**: SPC has Monte Carlo variance
3.  **Reliability metric**: If using different metrics, expect
    systematic differences

**Agreement criterion**: Differences \< 5% indicate good agreement.

### When They Disagree

If EQC and SPC differ by more than 5%:

``` r
# 1. Check if using same reliability metric
cat(sprintf("EQC metric: %s\n", eqc_result$metric))
cat(sprintf("SPC metric: %s\n", spc_result$metric))

# 2. Increase SPC iterations for more stable estimate
spc_longer <- spc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  c_init = eqc_result,
  n_iter = 500,  # More iterations
  burn_in = 250,
  seed = 42
)

# 3. Check SPC convergence
cat(sprintf("Converged: %s\n", spc_longer$convergence$converged))
cat(sprintf("Post-burn-in SD: %.4f\n", spc_longer$convergence$sd_post_burn))
```

## Monte Carlo Validation Workflow

For rigorous validation, use multiple replications to account for
sampling variability.

### Complete Validation Example

``` r
# ==== Configuration ====
target_rho <- 0.80
n_items <- 25
n_persons <- 1000
n_replications <- 50

# ==== Step 1: Calibrate ====
eqc_result <- eqc_calibrate(
  target_rho = target_rho,
  n_items = n_items,
  model = "rasch",
  latent_shape = "normal",
  item_source = "irw",
  M = 20000,  # Large quadrature for precision
  seed = 42
)

cat(sprintf("Calibrated c* = %.4f\n", eqc_result$c_star))
cat(sprintf("EQC achieved rho = %.4f\n\n", eqc_result$achieved_rho))

# ==== Step 2: Monte Carlo Validation ====
wle_rels <- numeric(n_replications)
eap_rels <- numeric(n_replications)

cat("Running Monte Carlo validation...\n")

for (r in 1:n_replications) {
  # Generate data
  sim_data <- simulate_response_data(
    eqc_result = eqc_result,
    n_persons = n_persons,
    latent_shape = "normal",
    seed = r
  )
  
  # Compute TAM reliability
  tam_rel <- compute_reliability_tam(
    resp = sim_data$response_matrix,
    model = "rasch",
    verbose = FALSE
  )
  
  wle_rels[r] <- tam_rel$rel_wle
  eap_rels[r] <- tam_rel$rel_eap
  
  if (r %% 10 == 0) cat(sprintf("  Completed %d/%d\n", r, n_replications))
}

# ==== Step 3: Summarize Results ====
cat("\n")
cat("=======================================================\n")
cat("  Monte Carlo Validation Results\n")
cat("=======================================================\n\n")

cat(sprintf("Target reliability: %.3f\n\n", target_rho))

cat("WLE Reliability:\n")
cat(sprintf("  Mean:   %.4f\n", mean(wle_rels)))
cat(sprintf("  SD:     %.4f\n", sd(wle_rels)))
cat(sprintf("  Range:  [%.4f, %.4f]\n", min(wle_rels), max(wle_rels)))
cat(sprintf("  MAE:    %.4f\n\n", mean(abs(wle_rels - target_rho))))

cat("EAP Reliability:\n")
cat(sprintf("  Mean:   %.4f\n", mean(eap_rels)))
cat(sprintf("  SD:     %.4f\n", sd(eap_rels)))
cat(sprintf("  Range:  [%.4f, %.4f]\n", min(eap_rels), max(eap_rels)))
cat(sprintf("  MAE:    %.4f\n", mean(abs(eap_rels - target_rho))))
```

### Interpreting Monte Carlo Results

**Success criteria**:

| Metric    | Good    | Acceptable | Investigate |
|-----------|---------|------------|-------------|
| MAE (WLE) | \< 0.02 | \< 0.03    | \> 0.05     |
| MAE (EAP) | \< 0.02 | \< 0.03    | \> 0.05     |
| SD (WLE)  | \< 0.02 | \< 0.03    | \> 0.04     |

### Visualizing Results

``` r
# Create validation plot
par(mfrow = c(1, 2))

# WLE distribution
hist(wle_rels, breaks = 15, col = "lightblue",
     main = "WLE Reliability Distribution",
     xlab = "Reliability", xlim = c(target_rho - 0.1, target_rho + 0.1))
abline(v = target_rho, col = "red", lwd = 2, lty = 2)
abline(v = mean(wle_rels), col = "blue", lwd = 2)
legend("topright", c("Target", "Mean"), 
       col = c("red", "blue"), lty = c(2, 1), lwd = 2)

# EAP distribution
hist(eap_rels, breaks = 15, col = "lightgreen",
     main = "EAP Reliability Distribution",
     xlab = "Reliability", xlim = c(target_rho - 0.1, target_rho + 0.1))
abline(v = target_rho, col = "red", lwd = 2, lty = 2)
abline(v = mean(eap_rels), col = "darkgreen", lwd = 2)
legend("topright", c("Target", "Mean"), 
       col = c("red", "darkgreen"), lty = c(2, 1), lwd = 2)

par(mfrow = c(1, 1))
```

## Validation Across Multiple Targets

Test calibration accuracy across a range of target reliabilities:

``` r
# Test multiple target levels
targets <- c(0.60, 0.70, 0.80, 0.85)
results <- data.frame(
  target = targets,
  c_star = NA,
  achieved_rho = NA,
  mean_wle = NA,
  mean_eap = NA,
  mae_wle = NA,
  mae_eap = NA
)

for (i in seq_along(targets)) {
  # Calibrate
  eqc <- eqc_calibrate(
    target_rho = targets[i],
    n_items = 25,
    model = "rasch",
    seed = 42
  )
  
  results$c_star[i] <- eqc$c_star
  results$achieved_rho[i] <- eqc$achieved_rho
  
  # Validate with 20 replications
  wle <- eap <- numeric(20)
  for (r in 1:20) {
    sim <- simulate_response_data(eqc, n_persons = 1000, seed = r)
    tam <- compute_reliability_tam(sim$response_matrix, "rasch", verbose = FALSE)
    wle[r] <- tam$rel_wle
    eap[r] <- tam$rel_eap
  }
  
  results$mean_wle[i] <- mean(wle)
  results$mean_eap[i] <- mean(eap)
  results$mae_wle[i] <- mean(abs(wle - targets[i]))
  results$mae_eap[i] <- mean(abs(eap - targets[i]))
  
  cat(sprintf("Target %.2f: c* = %.3f, MAE(WLE) = %.4f\n", 
              targets[i], eqc$c_star, results$mae_wle[i]))
}

# Summary table
print(results)
```

## Troubleshooting Validation Issues

### Issue 1: Large Discrepancy Between Target and Achieved

**Symptoms**: MAE \> 0.05, systematic bias

**Possible causes and solutions**:

``` r
# 1. Latent distribution mismatch
# Solution: Ensure same distribution for calibration and simulation
eqc <- eqc_calibrate(..., latent_shape = "bimodal", latent_params = list(delta = 0.8))
sim <- simulate_response_data(eqc, ..., 
                               latent_shape = "bimodal",  # MUST match
                               latent_params = list(delta = 0.8))

# 2. Small sample size
# Solution: Increase n_persons
sim <- simulate_response_data(eqc, n_persons = 2000, ...)  # Larger N

# 3. Insufficient quadrature
# Solution: Increase M in EQC
eqc <- eqc_calibrate(..., M = 50000)  # Larger quadrature
```

### Issue 2: WLE Reliability Much Lower Than Target

**Symptoms**: WLE consistently 0.03-0.05 below target

**Explanation**: This is often expected behavior due to the conservative
nature of WLE reliability.

``` r
# Check if EAP is closer to target
cat(sprintf("Target: %.3f\n", target_rho))
cat(sprintf("WLE:    %.3f (diff: %+.3f)\n", tam_rel$rel_wle, tam_rel$rel_wle - target_rho))
cat(sprintf("EAP:    %.3f (diff: %+.3f)\n", tam_rel$rel_eap, tam_rel$rel_eap - target_rho))

# If EAP is close but WLE is low, this is normal
# The true reliability is between WLE and EAP
```

### Issue 3: High Variability Across Replications

**Symptoms**: SD of reliability estimates \> 0.03

**Solutions**:

``` r
# 1. Increase sample size per replication
sim <- simulate_response_data(eqc, n_persons = 2000, ...)

# 2. Use more items
eqc <- eqc_calibrate(..., n_items = 40)

# 3. For heavy-tailed latent distributions, use even larger N
```

### Issue 4: EQC and SPC Disagree by \> 10%

**Symptoms**:
[`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
shows large percentage difference

**Solutions**:

``` r
# 1. Ensure both use same reliability metric
eqc <- eqc_calibrate(..., reliability_metric = "msem")
spc <- spc_calibrate(..., reliability_metric = "msem")

# 2. Increase SPC iterations
spc <- spc_calibrate(..., n_iter = 500, burn_in = 250)

# 3. Check SPC convergence
if (!spc$convergence$converged) {
  warning("SPC did not converge - increase n_iter")
}

# 4. Some difference is theoretically expected
# See SPC vignette for Jensen's inequality explanation
```

### Issue 5: Target Reliability Not Achievable

**Symptoms**: EQC hits upper bound, warning message

**Solutions**:

``` r
# 1. Use "info" metric (yields higher reliability for same c)
eqc <- eqc_calibrate(..., reliability_metric = "info")

# 2. Extend search bounds
eqc <- eqc_calibrate(..., c_bounds = c(0.1, 5))

# 3. Increase number of items
eqc <- eqc_calibrate(..., n_items = 40)

# 4. Accept maximum achievable reliability
# Check: eqc_result$misc$rho_bounds["rho_U"]
```

## Complete Validation Template

Use this template for your simulation studies:

``` r
# ==============================================================
# IRTsimrel Validation Template
# ==============================================================

library(IRTsimrel)
library(TAM)

# ---- Configuration ----
TARGET_RHO <- 0.80
N_ITEMS <- 25
N_PERSONS <- 1000
N_REPS <- 50
MODEL <- "rasch"
LATENT_SHAPE <- "normal"
SEED_CALIB <- 42

# ---- Step 1: Calibration ----
cat("Step 1: Running EQC calibration...\n")

eqc_result <- eqc_calibrate(
  target_rho = TARGET_RHO,
  n_items = N_ITEMS,
  model = MODEL,
  latent_shape = LATENT_SHAPE,
  item_source = "irw",
  M = 20000,
  seed = SEED_CALIB,
  verbose = TRUE
)

# ---- Step 2: SPC Validation (Optional) ----
cat("\nStep 2: Running SPC validation...\n")

spc_result <- spc_calibrate(
  target_rho = TARGET_RHO,
  n_items = N_ITEMS,
  model = MODEL,
  latent_shape = LATENT_SHAPE,
  item_source = "irw",
  c_init = eqc_result,
  n_iter = 200,
  seed = SEED_CALIB,
  verbose = TRUE
)

compare_eqc_spc(eqc_result, spc_result)

# ---- Step 3: Monte Carlo TAM Validation ----
cat("\nStep 3: Running Monte Carlo validation...\n")

wle_rels <- eap_rels <- numeric(N_REPS)

for (r in 1:N_REPS) {
  sim_data <- simulate_response_data(
    eqc_result = eqc_result,
    n_persons = N_PERSONS,
    latent_shape = LATENT_SHAPE,
    seed = r
  )
  
  tam_rel <- compute_reliability_tam(
    resp = sim_data$response_matrix,
    model = MODEL,
    verbose = FALSE
  )
  
  wle_rels[r] <- tam_rel$rel_wle
  eap_rels[r] <- tam_rel$rel_eap
  
  if (r %% 10 == 0) cat(sprintf("  %d/%d complete\n", r, N_REPS))
}

# ---- Step 4: Summary ----
cat("\n")
cat("=======================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("=======================================================\n\n")

cat(sprintf("Configuration:\n"))
cat(sprintf("  Target reliability: %.3f\n", TARGET_RHO))
cat(sprintf("  Items: %d, Persons: %d, Replications: %d\n\n", N_ITEMS, N_PERSONS, N_REPS))

cat(sprintf("Calibration:\n"))
cat(sprintf("  EQC c*: %.4f (achieved rho = %.4f)\n", eqc_result$c_star, eqc_result$achieved_rho))
cat(sprintf("  SPC c*: %.4f (agreement: %s)\n\n", 
            spc_result$c_star, 
            ifelse(abs(eqc_result$c_star - spc_result$c_star)/eqc_result$c_star < 0.05, "YES", "NO")))

cat(sprintf("TAM Validation:\n"))
cat(sprintf("  WLE: Mean = %.4f, SD = %.4f, MAE = %.4f\n", 
            mean(wle_rels), sd(wle_rels), mean(abs(wle_rels - TARGET_RHO))))
cat(sprintf("  EAP: Mean = %.4f, SD = %.4f, MAE = %.4f\n", 
            mean(eap_rels), sd(eap_rels), mean(abs(eap_rels - TARGET_RHO))))

cat(sprintf("\nValidation Status: "))
if (mean(abs(wle_rels - TARGET_RHO)) < 0.03) {
  cat("PASSED\n")
} else {
  cat("INVESTIGATE\n")
}
```

## Summary

Effective validation of reliability-targeted simulation requires:

1.  **Generate representative data** using
    [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
    with matching latent distributions
2.  **External validation** via TAM’s
    [`WLErel()`](https://rdrr.io/pkg/TAM/man/WLErel.html) and
    [`EAPrel()`](https://rdrr.io/pkg/TAM/man/WLErel.html) functions
3.  **Understand reliability metrics**: WLE is conservative, EAP is
    liberal, target falls between
4.  **Cross-validate algorithms** using
    [`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
    when possible
5.  **Use Monte Carlo replications** to account for sampling variability
6.  **Set realistic expectations**: MAE \< 0.03 is excellent, \< 0.05 is
    acceptable

Following these practices ensures your simulation studies achieve the
intended reliability levels with confidence.

## References

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427-450.

Robitzsch, A., Kiefer, T., & Wu, M. (2022). TAM: Test Analysis Modules.
R package version 4.1-4.
