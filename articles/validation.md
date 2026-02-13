# Validating Calibration Results

``` r
library(IRTsimrel)
set.seed(42)
```

## Overview

Validation is the process of confirming that a calibrated test design
actually achieves the intended reliability when data are generated from
it. In the IRTsimrel framework, validation operates at three levels:

| Level                  | Method                                                                                                                                           | What It Checks                                                  |
|:-----------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------|
| **1. Internal**        | [`predict()`](https://rdrr.io/r/stats/predict.html), [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md) | Achieved reliability matches target; Jensen’s gap is understood |
| **2. Cross-algorithm** | [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)                                                        | EQC and SAC agree on $c^{*}$                                    |
| **3. External**        | [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)                                        | Fitted-model reliability (WLE/EAP) matches theoretical target   |

**Reading time**: approximately 20 minutes.

**Prerequisites**: Basic familiarity with
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
and
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).
For mathematical foundations, see
[`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md).
For algorithm details, see
[`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)
and
[`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md).

## Level 1: Internal Validation

Internal validation checks whether the calibration output is
self-consistent, without fitting any external model to simulated data.

### The `predict()` Method

After running
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
the [`predict()`](https://rdrr.io/r/stats/predict.html) method evaluates
the reliability function at arbitrary scaling factors using the stored
quadrature samples.

``` r
# Calibrate
eqc_result <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  M = 5000L, seed = 42, verbose = FALSE
)

# Internal consistency check
cat(sprintf("Target reliability:  %.4f\n", eqc_result$target_rho))
#> Target reliability:  0.8000
cat(sprintf("Achieved reliability: %.4f\n", eqc_result$achieved_rho))
#> Achieved reliability: 0.8000
cat(sprintf("Absolute error:      %.6f\n",
            abs(eqc_result$achieved_rho - eqc_result$target_rho)))
#> Absolute error:      0.000000
```

The achieved reliability should match the target within the EQC
tolerance (default: $10^{- 4}$).

### Jensen’s Gap via `compute_rho_both()`

A critical internal diagnostic is measuring the gap between the two
reliability metrics at the calibrated $c^{*}$. This reveals how much the
choice of metric matters for the specific test configuration.

``` r
# Compute both metrics at c*
theta_val <- sim_latentG(5000, shape = "normal")$theta
items_val <- sim_item_params(25, model = "rasch", source = "parametric")
beta_val  <- items_val$data$beta
lambda_val <- rep(1, 25)

both <- compute_rho_both(eqc_result$c_star, theta_val, beta_val, lambda_val)

cat(sprintf("At calibrated c* = %.4f:\n", eqc_result$c_star))
#> At calibrated c* = 0.8995:
cat(sprintf("  rho_tilde (info): %.4f\n", both$rho_tilde))
#>   rho_tilde (info): 0.8000
cat(sprintf("  rho_bar   (msem): %.4f\n", both$rho_bar))
#>   rho_bar   (msem): 0.7937
cat(sprintf("  Jensen's gap:     %.4f\n", both$rho_tilde - both$rho_bar))
#>   Jensen's gap:     0.0063
```

**Interpretation guidelines**:

| Jensen’s gap | Interpretation                                           |
|:-------------|:---------------------------------------------------------|
| \< 0.01      | Negligible; metric choice does not matter                |
| 0.01–0.03    | Small; typical for normal latent distributions           |
| 0.03–0.05    | Moderate; consider which metric is more appropriate      |
| \> 0.05      | Large; investigate latent distribution or item structure |

### Jensen’s Gap Across Latent Distributions

The gap depends on how much test information varies across the ability
continuum. Non-standard latent distributions typically produce larger
gaps.

``` r
shapes <- c("normal", "bimodal", "heavy_tail", "skew_pos")
shape_pars <- list(
  normal     = list(),
  bimodal    = list(delta = 0.9),
  heavy_tail = list(df = 5),
  skew_pos   = list(k = 4)
)

cat("Jensen's gap by latent shape (c* from info metric, target = 0.80):\n")
#> Jensen's gap by latent shape (c* from info metric, target = 0.80):
cat(sprintf("  %-12s %-8s %-10s %-10s %-8s\n",
            "Shape", "c*", "rho_tilde", "rho_bar", "Gap"))
#>   Shape        c*       rho_tilde  rho_bar    Gap

for (sh in shapes) {
  eqc_sh <- eqc_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    latent_shape = sh, latent_params = shape_pars[[sh]],
    M = 5000L, seed = 42, verbose = FALSE
  )

  theta_sh <- sim_latentG(5000, shape = sh, shape_params = shape_pars[[sh]])$theta
  both_sh <- compute_rho_both(eqc_sh$c_star, theta_sh, beta_val, lambda_val)

  cat(sprintf("  %-12s %-8.4f %-10.4f %-10.4f %-8.4f\n",
              sh, eqc_sh$c_star, both_sh$rho_tilde, both_sh$rho_bar,
              both_sh$rho_tilde - both_sh$rho_bar))
}
#>   normal       0.8995   0.8006     0.7944     0.0062
#> Auto-wrapping shape parameter(s) {delta} into latent_params$shape_params.
#>   bimodal      0.9605   0.8133     0.8106     0.0027
#> Auto-wrapping shape parameter(s) {df} into latent_params$shape_params.
#>   heavy_tail   0.9365   0.8138     0.7952     0.0186
#> Auto-wrapping shape parameter(s) {k} into latent_params$shape_params.
#>   skew_pos     0.9132   0.8008     0.7885     0.0123
```

## Level 2: Cross-Algorithm Validation

Cross-algorithm validation uses two independent calibration algorithms
to solve the same problem. Agreement between them provides strong
evidence that the calibration is correct.

### The `compare_eqc_sac()` Function

``` r
# Run EQC
eqc_cross <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  M = 5000L, seed = 42, verbose = FALSE
)

# Run SAC with EQC warm start
sac_cross <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_cross, n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

# Formal comparison
comparison <- compare_eqc_sac(eqc_cross, sac_cross)
#> 
#> =======================================================
#>   EQC vs SAC Comparison
#> =======================================================
#> 
#>   Target reliability  : 0.8000
#>   EQC c*              : 0.899499
#>   SAC c*              : 0.912155
#>   Absolute difference : 0.012656
#>   Percent difference  : 1.41%
#>   Agreement (< 5%)    : YES
#> 
```

### Interpreting the Comparison

``` r
cat("\nComparison summary:\n")
#> 
#> Comparison summary:
cat(sprintf("  EQC c*:          %.4f\n", comparison$c_eqc))
#>   EQC c*:          0.8995
cat(sprintf("  SAC c*:          %.4f\n", comparison$c_sac))
#>   SAC c*:          0.9122
cat(sprintf("  Absolute diff:   %.4f\n", comparison$diff_abs))
#>   Absolute diff:   0.0127
cat(sprintf("  Percent diff:    %.2f%%\n", comparison$diff_pct))
#>   Percent diff:    1.41%
cat(sprintf("  Agreement (<5%%): %s\n", comparison$agreement))
#>   Agreement (<5%): TRUE
```

### Agreement Criteria

| Percent difference | Interpretation      | Action                                 |
|:-------------------|:--------------------|:---------------------------------------|
| \< 2%              | Excellent agreement | Proceed with confidence                |
| 2%–5%              | Good agreement      | Acceptable for most applications       |
| 5%–10%             | Marginal agreement  | Investigate; increase SAC iterations   |
| \> 10%             | Poor agreement      | Check metrics, convergence, and bounds |

### Warm Start Benefits

Using EQC to initialize SAC dramatically improves convergence and
agreement.

``` r
# SAC without warm start (APC initialization)
sac_cold <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = NULL,  # APC start
  n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

# SAC with EQC warm start
sac_warm <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_cross,  # EQC warm start
  n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

cat("Effect of warm start on SAC-EQC agreement:\n")
#> Effect of warm start on SAC-EQC agreement:
cat(sprintf("  SAC (APC start): c* = %.4f, diff from EQC = %.4f\n",
            sac_cold$c_star, abs(sac_cold$c_star - eqc_cross$c_star)))
#>   SAC (APC start): c* = 0.9793, diff from EQC = 0.0798
cat(sprintf("  SAC (EQC start): c* = %.4f, diff from EQC = %.4f\n",
            sac_warm$c_star, abs(sac_warm$c_star - eqc_cross$c_star)))
#>   SAC (EQC start): c* = 0.9122, diff from EQC = 0.0127
```

## Level 3: External Validation with TAM

External validation fits an IRT model to simulated response data using
an independent software package (TAM) and compares the recovered
reliability to the design target.

### Requirements

``` r
# Install TAM if needed
install.packages("TAM")
```

### Generating Response Data

The
[`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
function generates item response matrices from the calibrated
parameters.

``` r
sim_data <- simulate_response_data(
  result = eqc_cross,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

cat(sprintf("Response matrix: %d persons x %d items\n",
            nrow(sim_data$response_matrix),
            ncol(sim_data$response_matrix)))
#> Response matrix: 1000 persons x 25 items
cat(sprintf("Theta: mean = %.3f, SD = %.3f\n",
            mean(sim_data$theta), sd(sim_data$theta)))
#> Theta: mean = 0.016, SD = 0.992
cat(sprintf("Item difficulty range: [%.2f, %.2f]\n",
            min(sim_data$beta), max(sim_data$beta)))
#> Item difficulty range: [-2.17, 1.45]
```

### Computing TAM Reliability

``` r
# Compute WLE and EAP reliability
tam_rel <- compute_reliability_tam(
  resp = sim_data$response_matrix,
  model = "rasch",
  verbose = FALSE
)

cat(sprintf("External validation (TAM):\n"))
cat(sprintf("  Target reliability:  %.4f\n", eqc_cross$target_rho))
cat(sprintf("  EQC achieved rho:    %.4f\n", eqc_cross$achieved_rho))
cat(sprintf("  TAM WLE reliability: %.4f\n", tam_rel$rel_wle))
cat(sprintf("  TAM EAP reliability: %.4f\n", tam_rel$rel_eap))
```

## Understanding WLE vs EAP Reliability

A critical aspect of external validation is understanding why TAM
reports two different reliability coefficients that bracket the target.

### Mathematical Definitions

**WLE Reliability** (design-effect based):

$$\rho_{\text{WLE}} = 1 - \frac{{\bar{s}}^{2}}{V_{\text{WLE}}}$$

where
${\bar{s}}^{2} = N^{- 1}\sum_{k = 1}^{N}\text{SE}^{2}\left( {\widehat{\theta}}_{k}^{\text{WLE}} \right)$
is the average squared standard error and $V_{\text{WLE}}$ is the sample
variance of WLE estimates
$\{{\widehat{\theta}}_{k}^{\text{WLE}}\}_{k = 1}^{N}$.

**EAP Reliability** (posterior variance based):

$$\rho_{\text{EAP}} = \frac{V_{\text{EAP}}}{V_{\text{EAP}} + {\bar{\sigma}}^{2}}$$

where $V_{\text{EAP}}$ is the sample variance of EAP estimates and
${\bar{\sigma}}^{2} = N^{- 1}\sum_{k = 1}^{N}\sigma_{k}^{2}$ is the
average posterior variance.

### Why EAP $\geq$ WLE Always

The inequality $\rho_{\text{EAP}} \geq \rho_{\text{WLE}}$ is a
consequence of how posterior shrinkage affects variance estimation.

**EAP shrinkage**: EAP estimates are shrunk toward the prior mean,
which:

- Reduces the variance of estimates: $V_{\text{EAP}} < V_{\text{WLE}}$
- Reduces the average error variance:
  ${\bar{\sigma}}^{2} < {\bar{s}}^{2}$
- The error reduction is proportionally larger than the variance
  reduction, resulting in a higher reliability ratio.

**Connection to Jensen’s inequality**: The EAP reliability relates to
the MSEM-based population reliability $\bar{w}$, while the WLE
reliability relates more closely to the average-information reliability
$\widetilde{\rho}$. The ordering
$\rho_{\text{EAP}} \geq \rho_{\text{WLE}}$ parallels the theoretical
$\widetilde{\rho} \geq \bar{w}$ inequality, though the finite-sample
estimators introduce additional variability.

### Practical Interpretation

| TAM Metric      | Interpretation              | Typical Relationship to $\rho^{*}$    |
|:----------------|:----------------------------|:--------------------------------------|
| WLE reliability | Conservative lower bound    | $\rho_{\text{WLE}} \lesssim \rho^{*}$ |
| EAP reliability | Liberal upper bound         | $\rho_{\text{EAP}} \gtrsim \rho^{*}$  |
| Midpoint        | Reasonable central estimate | $\approx \rho^{*}$                    |

**Expected patterns** (for $\rho^{*} = 0.80$):

- WLE reliability: approximately 0.77–0.82
- EAP reliability: approximately 0.82–0.86
- Both within $\pm 0.03$ of target is excellent
- WLE within $\pm 0.05$ of target is acceptable

## Monte Carlo Validation Workflow

For rigorous validation, use multiple replications to account for
sampling variability in both data generation and model fitting.

``` r
# ==== Configuration ====
target_rho <- 0.80
n_items <- 25
n_persons <- 1000
n_reps <- 50

# ==== Step 1: Calibrate ====
eqc_result <- eqc_calibrate(
  target_rho = target_rho,
  n_items = n_items,
  model = "rasch",
  latent_shape = "normal",
  item_source = "parametric",
  M = 20000L,
  seed = 42,
  verbose = TRUE
)

cat(sprintf("Calibrated c* = %.4f\n\n", eqc_result$c_star))

# ==== Step 2: Monte Carlo Validation ====
wle_rels <- eap_rels <- numeric(n_reps)

for (r in 1:n_reps) {
  sim_data <- simulate_response_data(
    result = eqc_result,
    n_persons = n_persons,
    latent_shape = "normal",
    seed = r
  )

  tam_rel <- compute_reliability_tam(
    resp = sim_data$response_matrix,
    model = "rasch",
    verbose = FALSE
  )

  wle_rels[r] <- tam_rel$rel_wle
  eap_rels[r] <- tam_rel$rel_eap

  if (r %% 10 == 0) cat(sprintf("  Completed %d/%d\n", r, n_reps))
}

# ==== Step 3: Summarize ====
cat("\nMonte Carlo Validation Results\n")
cat("==============================\n")
cat(sprintf("Target: %.3f\n\n", target_rho))
cat(sprintf("WLE: Mean = %.4f, SD = %.4f, MAE = %.4f\n",
            mean(wle_rels), sd(wle_rels),
            mean(abs(wle_rels - target_rho))))
cat(sprintf("EAP: Mean = %.4f, SD = %.4f, MAE = %.4f\n",
            mean(eap_rels), sd(eap_rels),
            mean(abs(eap_rels - target_rho))))
```

### Success Criteria

| Metric    | Good    | Acceptable | Investigate |
|:----------|:--------|:-----------|:------------|
| MAE (WLE) | \< 0.02 | \< 0.03    | \> 0.05     |
| MAE (EAP) | \< 0.02 | \< 0.03    | \> 0.05     |
| SD (WLE)  | \< 0.02 | \< 0.03    | \> 0.04     |
| SD (EAP)  | \< 0.01 | \< 0.02    | \> 0.03     |

## The 960-Condition Validation Study

The paper (Lee, 2025) reports a comprehensive validation study that
provides empirical evidence for the package’s accuracy across a wide
range of conditions.

### Study Design

The validation study crossed the following factors:

| Factor             | Levels | Values                                  |
|:-------------------|:------:|:----------------------------------------|
| IRT model          |   3    | Rasch, 2PL-copula, 2PL-conditional      |
| Latent shape       |   4    | Normal, bimodal, heavy-tailed, skewed   |
| Test length        |   4    | $I \in \{ 10,20,30,50\}$                |
| Target reliability |   4    | $\rho^{*} \in \{ 0.60,0.70,0.80,0.85\}$ |
| Replications       |   5    | Independent seeds per condition         |

Total: $3 \times 4 \times 4 \times 4 \times 5 = 960$ conditions.

For each condition, the study:

1.  Ran EQC calibration ($M = 20,000$).
2.  Generated response data ($N = 1,000$).
3.  Fitted the model in TAM.
4.  Computed WLE and EAP reliability.

### Key Findings

The validation study established the following:

1.  **Overall MAE**: The mean absolute error (MAE) between the EQC
    target and TAM WLE reliability was approximately 0.02 across all 960
    conditions, confirming the package achieves its design targets.

2.  **Metric ordering**: EAP reliability was consistently higher than
    WLE reliability, as predicted by the theoretical analysis in
    [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md).

3.  **EQC-SAC agreement**: When both algorithms were run on the same
    condition, percent differences were typically below 5%, confirming
    the theoretical consistency results.

4.  **Latent shape sensitivity**: Non-normal latent distributions
    (bimodal, heavy-tailed) showed slightly larger MAE values
    (approximately 0.03) compared to normal distributions (approximately
    0.015), reflecting the increased Jensen’s gap.

5.  **Test length effect**: Longer tests ($I = 50$) produced more stable
    reliability estimates (smaller SD across replications) compared to
    shorter tests ($I = 10$).

These results provide confidence that the calibration algorithms are
accurate across a wide range of realistic simulation scenarios.

## Validation Under Misspecification

An important practical question is: what happens when the latent
distribution used for data generation does not match the one used for
calibration?

### Demonstration: Calibration-Simulation Mismatch

``` r
# Calibrate under normal assumption
eqc_normal <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric",
  latent_shape = "normal",
  M = 5000L, seed = 42, verbose = FALSE
)

# Evaluate reliability under different simulation distributions
sim_shapes <- c("normal", "bimodal", "heavy_tail", "skew_pos")
sim_pars <- list(
  normal     = list(),
  bimodal    = list(delta = 0.9),
  heavy_tail = list(df = 5),
  skew_pos   = list(k = 4)
)

cat("Misspecification analysis:\n")
#> Misspecification analysis:
cat("  Calibration: normal, target = 0.80\n\n")
#>   Calibration: normal, target = 0.80
cat(sprintf("  %-14s %-10s %-10s %-10s\n",
            "Sim shape", "rho_tilde", "rho_bar", "Deviation"))
#>   Sim shape      rho_tilde  rho_bar    Deviation

for (sh in sim_shapes) {
  theta_mis <- sim_latentG(5000, shape = sh, shape_params = sim_pars[[sh]])$theta
  both_mis <- compute_rho_both(
    eqc_normal$c_star, theta_mis,
    eqc_normal$beta_vec, eqc_normal$lambda_base
  )
  dev <- both_mis$rho_tilde - eqc_normal$target_rho
  cat(sprintf("  %-14s %-10.4f %-10.4f %+10.4f\n",
              sh, both_mis$rho_tilde, both_mis$rho_bar, dev))
}
#>   normal         0.7987     0.7923        -0.0013
#>   bimodal        0.7950     0.7928        -0.0050
#>   heavy_tail     0.7949     0.7761        -0.0051
#>   skew_pos       0.8006     0.7897        +0.0006
```

### Implications

When the simulation distribution differs from the calibration
distribution:

- **Same distribution**: Achieved reliability matches target closely.
- **Different distribution**: Reliability deviates from target,
  sometimes substantially.

**Recommendation**: Always use the same latent distribution for
calibration and data generation. If you need to study misspecification
effects, do so deliberately and document the mismatch.

### Sensitivity to Difficulty-Ability Alignment

``` r
# Calibrate with default difficulty distribution (centered at 0)
eqc_aligned <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric",
  M = 5000L, seed = 42, verbose = FALSE
)

# What if abilities are shifted?
theta_shifted <- sim_latentG(5000, shape = "normal")$theta + 1.5  # shifted right

both_aligned <- compute_rho_both(
  eqc_aligned$c_star,
  sim_latentG(5000, shape = "normal")$theta,
  eqc_aligned$beta_vec, eqc_aligned$lambda_base
)

both_shifted <- compute_rho_both(
  eqc_aligned$c_star,
  theta_shifted,
  eqc_aligned$beta_vec, eqc_aligned$lambda_base
)

cat("Effect of ability-difficulty misalignment:\n")
#> Effect of ability-difficulty misalignment:
cat(sprintf("  Aligned (theta ~ N(0,1)):   rho_tilde = %.4f\n",
            both_aligned$rho_tilde))
#>   Aligned (theta ~ N(0,1)):   rho_tilde = 0.8002
cat(sprintf("  Shifted (theta ~ N(1.5,1)): rho_tilde = %.4f\n",
            both_shifted$rho_tilde))
#>   Shifted (theta ~ N(1.5,1)): rho_tilde = 0.7651
cat(sprintf("  Reliability change:         %+.4f\n",
            both_shifted$rho_tilde - both_aligned$rho_tilde))
#>   Reliability change:         -0.0351
```

## Troubleshooting Guide

### Issue 1: Large Discrepancy Between Target and Achieved Reliability

**Symptoms**: `|achieved_rho - target_rho| > 0.01`

| Possible Cause        | Diagnostic                            | Solution                       |
|:----------------------|:--------------------------------------|:-------------------------------|
| Small $M$             | Check `length(eqc_result$theta_quad)` | Increase `M` to 10,000+        |
| Target near boundary  | Check `eqc_result$misc$rho_bounds`    | Extend `c_bounds`              |
| Upper/lower bound hit | Check for warning messages            | Adjust `c_bounds` or `n_items` |

### Issue 2: WLE Reliability Much Lower Than Target

**Symptoms**: $\rho_{\text{WLE}} < \rho^{*} - 0.05$ consistently

This is often *expected behavior* rather than a calibration failure. WLE
reliability is a conservative measure (see Section 5 above).

**Diagnostic**: Check if EAP reliability is closer to the target. If the
midpoint $\left( \rho_{\text{WLE}} + \rho_{\text{EAP}} \right)/2$ is
near $\rho^{*}$, the calibration is working correctly.

### Issue 3: High Variability Across Replications

**Symptoms**: SD of TAM reliability estimates $> 0.03$

| Possible Cause      | Solution                           |
|:--------------------|:-----------------------------------|
| Small sample size   | Increase `n_persons` to 2,000+     |
| Few items           | Increase `n_items`                 |
| Heavy-tailed latent | Use larger `n_persons`             |
| Non-standard model  | Verify item parameter distribution |

### Issue 4: EQC and SAC Disagree by $> 10\%$

**Symptoms**:
[`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
shows large difference

**Step-by-step diagnosis**:

1.  Verify both use the same `reliability_metric`.
2.  Check SAC convergence: `sac_result$convergence$converged`.
3.  Increase SAC iterations: `n_iter = 500`.
4.  Check if SAC hit bounds: `sac_result$convergence$hit_lower_bound`.
5.  Run SAC with EQC warm start if not already doing so.

### Issue 5: Target Reliability Not Achievable

**Symptoms**: EQC returns boundary $c^{*}$ with warning.

``` r
# Use check_feasibility() to diagnose
feas <- check_feasibility(
  n_items = 15, model = "rasch",
  item_source = "parametric",
  c_bounds = c(0.1, 10),
  M = 5000L, seed = 42
)
#> 
#> =======================================================
#>   Feasibility Check: Achievable Reliability Range
#> =======================================================
#> 
#>   Number of items  : 15
#>   Model            : RASCH
#>   Latent shape     : normal
#>   Latent variance  : 1.0099
#>   c range          : [0.10, 10.00]
#>   Monte Carlo M    : 5000
#> 
#> Achievable Reliability Ranges:
#>   rho_tilde (info) : [0.0363, 0.9797]
#>   rho_bar   (msem) : [0.0000, 0.8513]
#> 
#> Note: rho_tilde >= rho_bar always (Jensen's inequality).
#>   Use rho_tilde range for EQC targets.
#>   Use rho_bar range for SAC targets.
```

**Solutions**:

1.  Use `reliability_metric = "info"` (yields higher values than
    `"msem"`).
2.  Extend `c_bounds` (e.g., `c(0.01, 20)`).
3.  Increase `n_items`.
4.  Accept the maximum achievable reliability from
    [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md).

### Issue 6: Reproducibility Failures

**Symptoms**: Different results with the same seed.

| Possible Cause                | Solution                                                             |
|:------------------------------|:---------------------------------------------------------------------|
| Different R version           | Document [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) |
| Different package version     | Pin package version                                                  |
| Parallel/multithreaded RNG    | Use sequential execution                                             |
| Seed not passed to both steps | Ensure `seed` is set in both calibration and simulation              |

## Complete Validation Template

This is a self-contained template that can be copied and adapted for any
simulation study.

``` r
# ==============================================================
# IRTsimrel Validation Template
# ==============================================================

library(IRTsimrel)

# ---- Configuration ----
TARGET_RHO     <- 0.80
N_ITEMS        <- 25
N_PERSONS      <- 1000
N_REPS         <- 50
MODEL          <- "rasch"
LATENT_SHAPE   <- "normal"
LATENT_PARAMS  <- list()
ITEM_SOURCE    <- "parametric"
ITEM_PARAMS    <- list()
SEED_CALIB     <- 42

# ---- Step 1: Feasibility Check ----
cat("Step 1: Feasibility check...\n")
feas <- check_feasibility(
  n_items = N_ITEMS, model = MODEL,
  latent_shape = LATENT_SHAPE,
  item_source = ITEM_SOURCE,
  c_bounds = c(0.1, 10),
  M = 10000L, seed = SEED_CALIB,
  latent_params = LATENT_PARAMS,
  item_params = ITEM_PARAMS
)

stopifnot(TARGET_RHO >= feas$rho_range_info[1],
          TARGET_RHO <= feas$rho_range_info[2])

# ---- Step 2: EQC Calibration ----
cat("\nStep 2: EQC calibration...\n")
eqc_result <- eqc_calibrate(
  target_rho = TARGET_RHO,
  n_items = N_ITEMS,
  model = MODEL,
  latent_shape = LATENT_SHAPE,
  item_source = ITEM_SOURCE,
  latent_params = LATENT_PARAMS,
  item_params = ITEM_PARAMS,
  reliability_metric = "info",
  M = 20000L,
  seed = SEED_CALIB,
  verbose = TRUE
)

# ---- Step 3: SAC Cross-Validation ----
cat("\nStep 3: SAC cross-validation...\n")
sac_result <- sac_calibrate(
  target_rho = TARGET_RHO,
  n_items = N_ITEMS,
  model = MODEL,
  latent_shape = LATENT_SHAPE,
  item_source = ITEM_SOURCE,
  latent_params = LATENT_PARAMS,
  item_params = ITEM_PARAMS,
  reliability_metric = "info",
  c_init = eqc_result,
  n_iter = 300L,
  M_per_iter = 1000L,
  seed = SEED_CALIB,
  verbose = TRUE
)

comparison <- compare_eqc_sac(eqc_result, sac_result)
stopifnot(comparison$agreement)

# ---- Step 4: Jensen's Gap Check ----
cat("\nStep 4: Jensen's gap check...\n")
theta_check <- sim_latentG(10000, shape = LATENT_SHAPE,
                            params = LATENT_PARAMS)$theta
both_check <- compute_rho_both(
  eqc_result$c_star, theta_check,
  eqc_result$beta_vec, eqc_result$lambda_base
)
cat(sprintf("  Jensen's gap: %.4f\n",
            both_check$rho_tilde - both_check$rho_bar))

# ---- Step 5: Monte Carlo TAM Validation ----
cat("\nStep 5: Monte Carlo TAM validation...\n")

library(TAM)

wle_rels <- eap_rels <- numeric(N_REPS)

for (r in 1:N_REPS) {
  sim_data <- simulate_response_data(
    result = eqc_result,
    n_persons = N_PERSONS,
    latent_shape = LATENT_SHAPE,
    latent_params = LATENT_PARAMS,
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

# ---- Step 6: Summary Report ----
cat("\n")
cat("=======================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("=======================================================\n\n")

cat(sprintf("Configuration:\n"))
cat(sprintf("  Target: %.3f | Items: %d | Model: %s\n",
            TARGET_RHO, N_ITEMS, MODEL))
cat(sprintf("  Shape: %s | N: %d | Reps: %d\n\n",
            LATENT_SHAPE, N_PERSONS, N_REPS))

cat(sprintf("Calibration:\n"))
cat(sprintf("  EQC c*: %.4f (achieved = %.4f)\n",
            eqc_result$c_star, eqc_result$achieved_rho))
cat(sprintf("  SAC c*: %.4f (agreement: %s)\n\n",
            sac_result$c_star,
            ifelse(comparison$agreement, "YES", "NO")))

cat(sprintf("TAM Validation:\n"))
cat(sprintf("  WLE: Mean=%.4f SD=%.4f MAE=%.4f\n",
            mean(wle_rels), sd(wle_rels),
            mean(abs(wle_rels - TARGET_RHO))))
cat(sprintf("  EAP: Mean=%.4f SD=%.4f MAE=%.4f\n",
            mean(eap_rels), sd(eap_rels),
            mean(abs(eap_rels - TARGET_RHO))))

verdict <- ifelse(mean(abs(wle_rels - TARGET_RHO)) < 0.03,
                  "PASSED", "INVESTIGATE")
cat(sprintf("\nVerdict: %s\n", verdict))
```

## Summary

Effective validation of reliability-targeted simulation requires
multiple complementary approaches:

| Level                | Tool                                                                                                                                             | Key Question                             |
|:---------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------|
| **Internal**         | [`predict()`](https://rdrr.io/r/stats/predict.html), [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md) | Does the calibration achieve its target? |
| **Cross-algorithm**  | [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)                                                        | Do independent algorithms agree?         |
| **External**         | [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)                                        | Does fitted-model reliability match?     |
| **Misspecification** | Manual analysis                                                                                                                                  | How robust is the calibration?           |

**Best practices**:

1.  Always check the Jensen’s gap with
    [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)
    to understand the impact of metric choice.
2.  Use EQC warm start for SAC cross-validation.
3.  Report both WLE and EAP reliability from TAM validation.
4.  Use $N_{\text{persons}} \geq 1,000$ for stable TAM estimates.
5.  Match the latent distribution between calibration and simulation.
6.  Document the seed and package version for reproducibility.

## References

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. *arXiv preprint arXiv:2512.16012*.

Warm, T. A. (1989). Weighted likelihood estimation of ability in item
response theory. *Psychometrika, 54*(3), 427–450.

Robitzsch, A., Kiefer, T., & Wu, M. (2022). TAM: Test Analysis Modules.
R package version 4.1-4.

Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
in a microcomputer environment. *Applied Psychological Measurement,
6*(4), 431–444.
