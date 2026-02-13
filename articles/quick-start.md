# Quick Start: Your First Calibration in 5 Minutes

## 1. Overview

This vignette walks you through a complete IRTsimrel workflow in about
five minutes. By the end you will know how to:

1.  **Calibrate** item parameters to hit a target marginal reliability.
2.  **Check feasibility** and visualize the reliability curve.
3.  **Generate** a simulated binary response dataset.
4.  **Extract and inspect** calibrated parameters using S3 methods.

**Estimated time:** 5 minutes.

**Prerequisites:** Only the IRTsimrel package is needed. All code chunks
in this vignette use `eval = TRUE` and will run on any system with
IRTsimrel installed.

## 2. Calibrate

The core function is
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md).
Give it a target reliability, the number of items, and a measurement
model, and it returns a calibrated scaling factor $c^{*}$ that makes the
population reliability match your target.

``` r
result <- eqc_calibrate(
  target_rho = 0.80,
  n_items    = 20,
  model      = "rasch",
  seed       = 42,
  M          = 5000L
)
```

Print the result to see the key quantities:

``` r
result
#> 
#> =======================================================
#>   Empirical Quadrature Calibration (EQC) Results
#> =======================================================
#> 
#> Calibration Summary:
#>   Model                        : RASCH
#>   Target reliability (rho*)    : 0.8000
#>   Achieved reliability         : 0.8000
#>   Absolute error               : 7.94e-07
#>   Scaling factor (c*)          : 1.0183
#> 
#> Design Parameters:
#>   Number of items (I)          : 20
#>   Quadrature points (M)        : 5000
#>   Reliability metric           : Average-information (tilde)
#>   Latent variance              : 1.0099
#> 
#> Convergence:
#>   Root status                  : uniroot_success
#>   Search bracket               : [0.300, 3.000]
#>   Bracket reliabilities        : [0.3054, 0.9471]
#> 
#> Parameter Summaries:
#>   theta:        mean = -0.014, sd = 1.005
#>   beta:         mean = -0.000, sd = 0.758, range = [-2.15, 1.12]
#>   lambda_base:  mean = 1.000, sd = 0.000
#>   lambda_scaled: mean = 1.018, sd = 0.000
```

There are two numbers to focus on in the output:

- **`c*`** (scaling factor): All baseline discriminations are multiplied
  by this value. In a Rasch model the baseline discriminations are all
  1, so the calibrated discriminations equal $c^{*}$ directly. A larger
  $c^{*}$ means the test needs more discriminating items to reach the
  target.

- **`achieved_rho`**: The empirical reliability at the calibrated
  $c^{*}$, computed over the Monte Carlo quadrature sample. This should
  be very close to the target of 0.80 — typically the absolute error is
  less than $10^{- 4}$.

The `model` argument accepts `"rasch"` (all discriminations equal) or
`"2pl"` (log-normal discrimination distribution). For this quick start
we use the Rasch model; see
[`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)
for 2PL examples.

## 3. Check Feasibility

Before committing to a particular simulation design, it is good practice
to verify that your target reliability is actually achievable. Not all
combinations of test length, model, and latent distribution can produce
every reliability level.

[`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
reports the range of achievable reliabilities for a given configuration:

``` r
feas <- check_feasibility(
  n_items = 20,
  model   = "rasch",
  seed    = 42,
  M       = 5000L
)
#> 
#> =======================================================
#>   Feasibility Check: Achievable Reliability Range
#> =======================================================
#> 
#>   Number of items  : 20
#>   Model            : RASCH
#>   Latent shape     : normal
#>   Latent variance  : 1.0099
#>   c range          : [0.10, 10.00]
#>   Monte Carlo M    : 5000
#> 
#> Achievable Reliability Ranges:
#>   rho_tilde (info) : [0.0479, 0.9850]
#>   rho_bar   (msem) : [0.0000, 0.8818]
#> 
#> Note: rho_tilde >= rho_bar always (Jensen's inequality).
#>   Use rho_tilde range for EQC targets.
#>   Use rho_bar range for SAC targets.
```

The output shows two ranges — one for each reliability metric. If your
target falls within the `rho_tilde (info)` range, EQC can calibrate for
it. If it falls within the `rho_bar (msem)` range, SAC can calibrate for
it.

You can also visualize how reliability varies continuously with the
scaling factor by using
[`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md).
This plot shows both the average-information metric ($\widetilde{\rho}$,
blue) and the MSEM-based metric ($\bar{w}$, red):

``` r
curve_data <- rho_curve(
  n_items = 20,
  model   = "rasch",
  metric  = "both",
  M       = 5000L,
  seed    = 42,
  plot    = TRUE
)
```

![Reliability curve for a 20-item Rasch test. The average-information
metric (blue) always lies at or above the MSEM-based metric (red) due to
Jensen's inequality.](quick-start_files/figure-html/rho-curve-1.png)

Reliability curve for a 20-item Rasch test. The average-information
metric (blue) always lies at or above the MSEM-based metric (red) due to
Jensen’s inequality.

The two curves are close together for this configuration. The gap widens
for shorter tests or non-normal latent distributions.

## 4. Generate Data

Once you have a calibration result,
[`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
generates a binary response matrix using the calibrated item parameters.
The function draws fresh latent abilities from the specified
distribution and produces responses according to the 2PL (or Rasch)
model with the calibrated discriminations:

``` r
sim <- simulate_response_data(
  result    = result,
  n_persons = 500,
  seed      = 123
)
```

The result is a list with four components:

- `response_matrix`: an $N \times I$ matrix of binary (0/1) responses
- `theta`: the true latent abilities for each person
- `beta`: the item difficulties
- `lambda`: the scaled item discriminations

``` r
# Dimensions: 500 persons x 20 items
dim(sim$response_matrix)
#> [1] 500  20
```

``` r
# First 6 persons, first 8 items
sim$response_matrix[1:6, 1:8]
#>      item1 item2 item3 item4 item5 item6 item7 item8
#> [1,]     0     1     0     1     0     0     0     0
#> [2,]     0     1     0     1     1     0     1     0
#> [3,]     1     1     1     1     1     1     1     1
#> [4,]     1     0     0     1     0     1     1     1
#> [5,]     1     0     0     0     1     1     0     0
#> [6,]     1     1     1     1     1     1     1     1
```

You can verify that the item parameters match what you expect:

``` r
# All scaled discriminations should equal c*
summary(sim$lambda)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.018   1.018   1.018   1.018   1.018   1.018

# Difficulty distribution
summary(sim$beta)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -2.14978 -0.33975  0.04838  0.00000  0.37403  1.12028
```

A quick diagnostic: plot the proportion of correct responses per item.
Items with extreme difficulties should have very high or very low
proportions:

``` r
p_correct <- colMeans(sim$response_matrix)
item_order <- order(sim$beta)

barplot(
  p_correct[item_order],
  names.arg = item_order,
  col  = "steelblue",
  xlab = "Item (ordered by difficulty)",
  ylab = "Proportion Correct",
  main = "Proportion Correct per Item",
  ylim = c(0, 1),
  las  = 2,
  cex.names = 0.7
)
abline(h = 0.5, lty = 2, col = "gray50")
```

![Proportion correct per item, ordered by difficulty. Items near the
center of the difficulty distribution have proportions near 0.50, while
extreme items show floor or ceiling
effects.](quick-start_files/figure-html/item-difficulty-plot-1.png)

Proportion correct per item, ordered by difficulty. Items near the
center of the difficulty distribution have proportions near 0.50, while
extreme items show floor or ceiling effects.

## 5. Extract and Inspect

IRTsimrel result objects support standard S3 generics, so you can
interact with them the same way you would with `lm` or `glm` objects.

### 5.1 Summary

[`summary()`](https://rdrr.io/r/base/summary.html) returns a compact
overview of the calibration:

``` r
summary(result)
#> Summary: Empirical Quadrature Calibration (EQC)
#> ================================================
#>   Model            : RASCH
#>   Metric           : Average-information (tilde)
#>   Number of items  : 20
#>   Quadrature (M)   : 5000
#>   Latent variance  : 1.0099
#> 
#> Calibration Results:
#>   Target rho*      : 0.8000
#>   Achieved rho     : 0.8000
#>   Absolute error   : 7.94e-07
#>   Scaling factor c*: 1.0183
#>   Root status      : uniroot_success
```

### 5.2 Extract Coefficients

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a tidy data frame
of calibrated item parameters — one row per item:

``` r
item_pars <- coef(result)
head(item_pars)
#>   item_id         beta lambda_base lambda_scaled   c_star
#> 1       1  0.217684472           1      1.018304 1.018304
#> 2       2  1.116752062           1      1.018304 1.018304
#> 3       3  0.456497288           1      1.018304 1.018304
#> 4       4  0.006913473           1      1.018304 1.018304
#> 5       5 -0.179849098           1      1.018304 1.018304
#> 6       6  0.027652529           1      1.018304 1.018304
```

For a Rasch model, `lambda_base` is always 1 and `lambda_scaled` equals
$c^{*}$ for every item. For a 2PL model, `lambda_base` would vary across
items and `lambda_scaled` would be `lambda_base * c*`.

### 5.3 Predict Reliability at New Scaling Factors

[`predict()`](https://rdrr.io/r/stats/predict.html) evaluates the
reliability function at new scaling factor values. With no arguments, it
returns the achieved reliability at the calibrated $c^{*}$:

``` r
# Achieved reliability at c*
predict(result)
#> [1] 0.8000008
```

With `newdata`, it computes reliability at each specified value of $c$:

``` r
# Reliability at several scaling factors
predict(result, newdata = c(0.5, 1.0, 1.5, 2.0))
#>     c=0.5     c=1.0     c=1.5     c=2.0 
#> 0.5369213 0.7952739 0.8784820 0.9150086
```

This is useful for understanding how sensitive reliability is to the
discrimination level. For instance, at $c = 0.5$ the test is
substantially less reliable, while at $c = 2.0$ it is more reliable than
needed.

## 6. Validate (Optional)

For added confidence in your calibration, you can run SAC (Stochastic
Approximation Calibration) as an independent check. SAC uses a
completely different algorithm (Robbins-Monro stochastic approximation)
to solve the same calibration problem, so agreement between EQC and SAC
is strong evidence that both found the correct $c^{*}$.

Passing the EQC result as `c_init` provides a warm start that makes SAC
converge in very few iterations:

``` r
sac_result <- sac_calibrate(
  target_rho = 0.80,
  n_items    = 20,
  model      = "rasch",
  c_init     = result,
  n_iter     = 100L,
  M_per_iter = 500L,
  seed       = 42
)
```

Compare the two algorithms side by side:

``` r
comp <- compare_eqc_sac(result, sac_result)
#> Warning in compare_eqc_sac(result, sac_result): Reliability metric differs
#> between EQC ('info') and SAC ('msem').
#> 
#> =======================================================
#>   EQC vs SAC Comparison
#> =======================================================
#> 
#>   Target reliability  : 0.8000
#>   EQC c*              : 1.018304
#>   SAC c*              : 1.060003
#>   Absolute difference : 0.041699
#>   Percent difference  : 4.09%
#>   Agreement (< 5%)    : YES
#> 
```

The comparison reports the absolute and percent difference between the
two $c^{*}$ values. If the percent difference is below 5%, the two
algorithms are in agreement (and typically the difference is below 1%).

You can also visualize the SAC convergence trajectory to verify that the
iterations stabilized near the EQC warm start:

``` r
plot(sac_result, type = "both")
```

![SAC convergence trajectory. The top panel shows the scaling factor c
across iterations, and the bottom panel shows the per-iteration
reliability estimates. The warm start from EQC ensures rapid
convergence.](quick-start_files/figure-html/sac-plot-1.png)

SAC convergence trajectory. The top panel shows the scaling factor c
across iterations, and the bottom panel shows the per-iteration
reliability estimates. The warm start from EQC ensures rapid
convergence.

The top panel shows the scaling factor trajectory, which should
stabilize quickly when initialized from EQC. The bottom panel shows the
noisy per-iteration reliability estimates oscillating around the target
value (dashed red line). The Polyak-Ruppert average (reported as
$c^{*}$) smooths out the iteration-to-iteration noise.

## 7. Putting It All Together

Here is the complete workflow in a single code block, from calibration
to validated data generation:

``` r
library(IRTsimrel)

# Step 1: Check feasibility
feas <- check_feasibility(n_items = 20, model = "rasch", seed = 42)

# Step 2: Calibrate with EQC
eqc_res <- eqc_calibrate(
  target_rho = 0.80, n_items = 20, model = "rasch",
  seed = 42, M = 5000L
)

# Step 3: Validate with SAC (optional but recommended)
sac_res <- sac_calibrate(
  target_rho = 0.80, n_items = 20, model = "rasch",
  c_init = eqc_res, n_iter = 100L, seed = 42
)
compare_eqc_sac(eqc_res, sac_res)

# Step 4: Generate response data
sim_data <- simulate_response_data(
  result = eqc_res, n_persons = 1000, seed = 123
)

# Step 5: Use the data in your analysis
dim(sim_data$response_matrix)  # 1000 x 20
```

## 8. What’s Next?

You have now completed a full calibrate-generate-validate cycle. Here
are some directions for deeper exploration:

- **[`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)**:
  Comprehensive applied tutorial covering 2PL models, non-normal latent
  distributions, IRW-based item sources, and factorial simulation
  designs with multiple reliability levels.

- **[`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md)**:
  Explore all 12 latent distribution shapes available in
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
  and learn when to use each one for different research scenarios.

- **[`vignette("item-parameters")`](https://joonho112.github.io/IRTsimrel/articles/item-parameters.md)**:
  Parametric, IRW, hierarchical, and custom item generation methods,
  including correlated difficulty-discrimination parameters.

- **[`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)**:
  Mathematical foundations of the two reliability metrics
  ($\widetilde{\rho}$ and $\bar{w}$), Jensen’s inequality, and the
  theoretical justification for the calibration approach.

- **[`vignette("api-reference")`](https://joonho112.github.io/IRTsimrel/articles/api-reference.md)**:
  Full function reference with complete signatures, all arguments,
  return values, and runnable examples for every exported function.

## References

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. *arXiv preprint*, arXiv:2512.16012.
