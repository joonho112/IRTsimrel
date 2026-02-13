# IRTsimrel: Reliability-Targeted IRT Simulation

## 1. Introduction

In Monte Carlo simulation studies for Item Response Theory (IRT),
researchers carefully manipulate sample size, test length, item
parameter distributions, and latent trait shapes. Yet **marginal
reliability** — the single number that best summarises how informative a
test is — is almost never directly controlled. Instead, it emerges as an
implicit, unreported consequence of these design choices. The
**IRTsimrel** package closes that gap by letting you specify a target
reliability as an explicit input parameter, then calibrating the
data-generating process so the simulated test hits that target.

IRTsimrel implements the framework described in Lee (2025). It provides
two calibration algorithms (EQC and SAC), a rich set of
latent-distribution and item-parameter generators, and diagnostic tools
for feasibility screening and reliability visualization. Whether you are
running a large-scale Monte Carlo study or generating a single dataset
for classroom use, the package gives you principled control over the
signal-to-noise ratio in your simulated data.

This introductory vignette is the gateway to the rest of the package
documentation. It provides a high-level overview of the core ideas, the
two calibration algorithms, the two reliability metrics, and the full
set of exported functions. At the end, a **reading roadmap** helps you
find the vignette most relevant to your goals, whether you are an
applied researcher who just needs to simulate data, or a methodological
researcher who wants to understand the mathematical foundations.

## 2. The Core Challenge: Why Reliability Matters

### 2.1 The ICC Analogy

Multilevel modelling (MLM) researchers would never run a simulation
study without specifying the **intraclass correlation (ICC)**. The ICC
controls the signal-to-noise ratio — the proportion of total variance
attributable to the cluster-level random effect versus residual noise.
Every credible MLM simulation paper reports and systematically varies
ICC as a primary design factor.

**Marginal reliability in IRT serves the same conceptual role as ICC in
MLM.** Both quantities capture the ratio of signal variance to total
variance. In IRT the “signal” is the latent trait $\theta$ and the
“noise” is measurement error. Despite this parallel, the vast majority
of IRT simulation studies never report or control the implied
reliability of their simulated data.

Consider how ICC appears in a typical MLM simulation paper: the authors
might specify ICC values of 0.05, 0.15, and 0.30, and systematically
evaluate how model performance varies across these conditions. By
contrast, an IRT simulation paper might specify item discriminations of
$\lambda = 1.0$ and test length of $I = 30$, but never compute or report
the implied reliability of the resulting test. As Lee (2025) argues,
this is a serious methodological gap that undermines the
interpretability and generalizability of simulation results.

IRTsimrel brings IRT simulation methodology into alignment with this
long-standing best practice in multilevel modelling.

### 2.2 Three Consequences of the “Reliability Omission”

When reliability is left uncontrolled in IRT simulation studies, three
interrelated problems follow:

1.  **Ecological validity threat.** Operational assessments commonly
    have reliabilities between 0.50 and 0.70 (especially short screening
    instruments, formative assessments, and subscales). However,
    simulation studies frequently use item parameters that imply
    reliabilities above 0.90. Conclusions drawn under such optimistic
    conditions may not generalize to the lower-reliability regime that
    characterizes many real-world applications. When a researcher
    reports that a new estimator performs well, the reader has no way of
    knowing whether the same would hold at $\rho = 0.60$.

2.  **Confounded comparisons.** Claims that one model or estimator
    outperforms another may hold only within a narrow reliability
    regime. Two studies that ostensibly compare the same pair of methods
    may reach opposite conclusions simply because they used different
    (and unreported) reliability levels. Without systematically varying
    reliability as a factor in the simulation design, such conclusions
    are fragile and potentially misleading.

3.  **Limited replicability.** If the implied reliability is not
    documented, no other researcher can construct a comparable
    data-generating process. Even if all item parameter values are
    reported, subtle choices about the latent distribution shape or the
    correlation between difficulty and discrimination can change the
    implied reliability. This undermines the replicability of simulation
    studies — precisely the domain where replicability should be easiest
    to achieve.

### 2.3 A Formal Definition

To fix ideas, consider the **average-information marginal reliability**,
which IRTsimrel targets by default:

$$\widetilde{\rho}\; = \;\frac{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}}{\sigma_{\theta}^{2}\,\bar{\mathcal{J}} + 1}$$

where $\sigma_{\theta}^{2}$ is the latent trait variance and
$\bar{\mathcal{J}} = {\mathbb{E}}_{\theta}\left\lbrack \mathcal{J}(\theta) \right\rbrack$
is the marginal (average) test information. The test information at a
given ability level is the sum of item information functions:

$$\mathcal{J}(\theta)\; = \;\sum\limits_{i = 1}^{I}\lambda_{i}^{2}\, p_{i}(\theta)\,\left\lbrack 1 - p_{i}(\theta) \right\rbrack$$

where
$p_{i}(\theta) = \text{logit}^{- 1}\left\lbrack \lambda_{i}\left( \theta - \beta_{i} \right) \right\rbrack$
is the probability of a correct response to item $i$. This formula makes
explicit that reliability depends on both the **item parameters**
(through $\bar{\mathcal{J}}$) and the **latent distribution** (through
$\sigma_{\theta}^{2}$ and the expectation over $\theta$). Changing
either one changes the implied reliability, which is precisely why it
must be actively controlled.

The reliability formula has an intuitive interpretation: it is a ratio
of signal to total variance, just like ICC. When average test
information $\bar{\mathcal{J}}$ is large (high discrimination, many
items), most of the variance in observed scores is due to true
differences in $\theta$. When $\bar{\mathcal{J}}$ is small, measurement
error dominates.

## 3. How IRTsimrel Works

### 3.1 The Key Idea: Separation of Structure and Scale

IRTsimrel decouples two aspects of the data-generating process:

| Aspect        | What it controls                                                                                      | Source                                                                                                                                                                                         |
|---------------|-------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Structure** | Difficulty distribution, discrimination patterns, difficulty–discrimination correlation, latent shape | [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md), [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)                   |
| **Scale**     | Global discrimination level $\Rightarrow$ overall informativeness $\Rightarrow$ reliability           | Calibration via [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md) or [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md) |

This separation is the central innovation. Structure captures all the
qualitative features of the data-generating process — how difficulties
are spread, whether discriminations vary across items, whether the
latent distribution is skewed or multimodal. These features come from
realistic generators (parametric distributions, the Item Response
Warehouse, or user-supplied custom parameters) and are preserved exactly
throughout calibration.

Scale controls only the **overall level of informativeness**, and thus
reliability. It operates through a single **scaling factor** $c^{*}$.
Given baseline (unscaled) discriminations $\lambda_{i,0}$, the
calibrated discriminations are:

$$\lambda_{i}^{*}\; = \; c^{*} \cdot \lambda_{i,0}$$

Increasing $c^{*}$ amplifies all discriminations proportionally, raising
test information and therefore reliability, while preserving the
relative structure among items (e.g., the ratio
$\lambda_{2}/\lambda_{1}$ remains constant). The calibration algorithm
finds the value $c^{*}$ that makes the population reliability equal to
the user’s target $\rho^{*}$.

### 3.2 How It Works in Practice

The workflow has three steps:

1.  **Generate structure.** Call
    [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)
    to draw latent abilities from a specified distribution, and
    [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)
    to draw item parameters from a specified source. These provide the
    qualitative “shape” of the data-generating process.

2.  **Calibrate scale.** Call
    [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
    (or
    [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md))
    with your target reliability. The function uses the generated
    structure to find the scaling factor $c^{*}$ that achieves
    $\rho^{*} = \rho_{\text{target}}$.

3.  **Generate data.** Call
    [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
    with the calibration result to produce a binary response matrix
    whose implied population reliability matches your target.

### 3.3 A Minimal Example

Here is a complete EQC calibration in five lines of code:

``` r
# Calibrate a 20-item Rasch test for target reliability 0.80
result <- eqc_calibrate(
  target_rho = 0.80,
  n_items    = 20,
  model      = "rasch",
  seed       = 42,
  M          = 5000L
)
```

``` r
# Inspect the calibration result
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

The output reports the calibrated scaling factor `c*`, the achieved
reliability, and convergence diagnostics. The absolute error between
target and achieved reliability is typically on the order of $10^{- 5}$
or smaller, confirming that the calibration solved the root-finding
problem to high precision.

### 3.4 From Calibration to Data

With the calibration result in hand, generating response data takes one
more call:

``` r
# Generate 500 response vectors from the calibrated DGP
sim <- simulate_response_data(
  result    = result,
  n_persons = 500,
  seed      = 123
)

# The result is a 500 x 20 binary matrix
dim(sim$response_matrix)
#> [1] 500  20
```

The resulting dataset has exactly the reliability properties you
specified. You can use it directly in any downstream analysis: model
fitting, parameter recovery studies, comparison of estimators, or power
analysis.

## 4. Two Algorithms

IRTsimrel provides two calibration algorithms. They solve the same
problem — find $c^{*}$ such that $\rho\left( c^{*} \right) = \rho^{*}$ —
but differ in their approach and strengths.

### 4.1 Algorithm 1: EQC (Empirical Quadrature Calibration)

EQC draws a large, fixed Monte Carlo sample (“quadrature”) of $M$
abilities and $I$ item parameters, then uses Brent’s root-finding method
([`uniroot()`](https://rdrr.io/r/stats/uniroot.html)) to solve
${\widehat{\rho}}_{M}(c) = \rho^{*}$ deterministically.

**How it works:**

1.  Draw $\{\theta_{m}\}_{m = 1}^{M} \sim G$ and
    $\{\left( \beta_{i},\lambda_{i,0} \right)\}_{i = 1}^{I} \sim H$
    once.
2.  For any candidate $c$, form $\lambda_{i}(c) = c \cdot \lambda_{i,0}$
    and compute the empirical reliability ${\widehat{\rho}}_{M}(c)$.
3.  Solve ${\widehat{\rho}}_{M}\left( c^{*} \right) = \rho^{*}$ via
    [`uniroot()`](https://rdrr.io/r/stats/uniroot.html).

**Strengths:** Fast (typically under one second), deterministic, exact
up to Monte Carlo noise, reproducible with a fixed seed. The recommended
default for routine simulation studies.

**When to use EQC:**

- Routine simulation studies where you need fast, reliable calibration.
- Rasch or 2PL models with standard distributions.
- As initialization (warm start) for SAC.

### 4.2 Algorithm 2: SAC (Stochastic Approximation Calibration)

SAC uses the Robbins-Monro stochastic approximation algorithm. At each
iteration it draws a fresh mini-batch of abilities and items, computes a
noisy reliability estimate, and takes a gradient-like step toward the
target. After a burn-in phase, Polyak-Ruppert averaging yields the final
estimate $c^{*}$.

**The update rule:**

$$c_{n + 1}\; = \; c_{n}\; - \; a_{n} \cdot \left( {\widehat{\rho}}_{n} - \rho^{*} \right)$$

where $a_{n} = a/(n + A)^{\gamma}$ is a decreasing step size satisfying
the Robbins-Monro conditions $\sum a_{n} = \infty$ and
$\sum a_{n}^{2} < \infty$.

**Strengths:** Handles non-monotone objectives (needed for $\bar{w}$
targeting), provides convergence diagnostics, works with arbitrary DGPs,
and can independently validate EQC results.

**When to use SAC:**

- Independent validation of EQC results.
- Targeting the exact MSEM-based reliability $\bar{w}$, where the
  objective may be non-monotone and root-finding can fail.
- Complex data-generating processes where the test information function
  is not easily evaluated in closed form.

**Warm start from EQC:** For the best of both worlds, run EQC first,
then pass the result to SAC as a warm start. This combines the speed of
EQC with the robustness of SAC:

``` r
eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20, ...)
sac_res <- sac_calibrate(target_rho = 0.80, n_items = 20,
                          c_init = eqc_res, ...)
```

### 4.3 Decision Guide

| Scenario                                           | Recommended              | Rationale                                        |
|----------------------------------------------------|--------------------------|--------------------------------------------------|
| Standard Rasch/2PL, need quick answer              | EQC                      | Deterministic, fast, precise                     |
| Targeting $\widetilde{\rho}$ (average-information) | EQC                      | Monotone objective, guaranteed convergence       |
| Targeting $\bar{w}$ (MSEM-based)                   | SAC                      | Handles non-monotone objective correctly         |
| Validating an EQC result                           | SAC (warm start)         | Independent confirmation via different algorithm |
| Complex custom DGP                                 | SAC                      | Only requires pointwise reliability estimates    |
| Large-scale Monte Carlo study                      | EQC + optional SAC check | Speed for production, SAC for spot checks        |

For most users in most situations, **EQC is the right choice**. Use SAC
when you need extra assurance or when the mathematical properties of
your specific problem require it.

## 5. Two Reliability Metrics

IRTsimrel supports two population-level reliability metrics from Lee
(2025). Understanding the distinction is important for choosing the
right algorithm and interpreting results.

### 5.1 Average-Information Reliability ($\widetilde{\rho}$)

$$\widetilde{\rho}(c)\; = \;\frac{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c)}{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c) + 1}\qquad\text{where}\quad\bar{\mathcal{J}}(c) = {\mathbb{E}}_{\theta}\left\lbrack \mathcal{J}(\theta;\, c) \right\rbrack$$

This metric uses the **arithmetic mean** of the test information
function across the latent distribution. Because the arithmetic mean is
a linear functional of $c^{2}$ (via the item information formula),
$\widetilde{\rho}(c)$ is **monotonically increasing** in $c$. This
monotonicity guarantees that Brent’s method will find the unique root,
making $\widetilde{\rho}$ the natural target for EQC.

The
[`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
function computes this metric directly from a scaling factor, theta
sample, and item parameters. Inside
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
this is the default metric (`reliability_metric = "info"`).

### 5.2 MSEM-Based Reliability ($\bar{w}$)

$$\bar{w}(c)\; = \;\frac{\sigma_{\theta}^{2}}{\sigma_{\theta}^{2} + \text{MSEM}(c)}\qquad\text{where}\quad\text{MSEM}(c) = {\mathbb{E}}_{\theta}\!\left\lbrack \frac{1}{\mathcal{J}(\theta;\, c)} \right\rbrack$$

This metric uses the **harmonic mean** of test information (via the mean
squared error of measurement). It is the theoretically exact marginal
reliability defined in terms of the expected reciprocal of test
information. However, its objective $\bar{w}(c) - \rho^{*}$ can be
non-monotone at extreme scaling factors (see Lee, 2025, Section 4.3),
which is why SAC is preferred for targeting $\bar{w}$ directly.

The
[`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)
function computes this metric, and
[`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)
computes both metrics in a single pass for efficiency.

### 5.3 Jensen’s Inequality: $\widetilde{\rho} \geq \bar{w}$

By Jensen’s inequality applied to the convex function $f(x) = 1/x$, the
harmonic mean of test information is always less than or equal to the
arithmetic mean. This implies:

$$\widetilde{\rho}(c)\; \geq \;\bar{w}(c)\quad{\text{for all}\mspace{6mu}}c > 0$$

The gap between the two metrics is small when test information is
roughly constant across the latent distribution (e.g., a well-targeted
Rasch test with normally distributed abilities) and larger when
information varies substantially (e.g., short tests with skewed latent
distributions or extreme difficulty values).

### 5.4 Visualizing the Reliability Curve

The
[`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
function visualizes both metrics simultaneously as a function of the
scaling factor $c$:

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

![Reliability as a function of the scaling factor c for a 20-item Rasch
test. The gap between the two curves reflects Jensen's inequality: the
average-information metric always lies above or on the MSEM-based
metric.](introduction_files/figure-html/rho-curve-1.png)

Reliability as a function of the scaling factor c for a 20-item Rasch
test. The gap between the two curves reflects Jensen’s inequality: the
average-information metric always lies above or on the MSEM-based
metric.

This plot is useful for several purposes: understanding the feasible
range of reliabilities, seeing how sensitive reliability is to the
discrimination level, and verifying that the gap between
$\widetilde{\rho}$ and $\bar{w}$ is small for your particular
configuration.

## 6. What the Package Offers

### 6.1 Exported Functions

IRTsimrel exports 16 user-facing functions organized into five
categories:

**Simulation generators:**

| Function                                                                                                | Purpose                                                                   |
|---------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------|
| [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)                       | Generate latent abilities from a flexible distribution family (12 shapes) |
| [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)               | Generate item parameters (parametric, IRW, hierarchical, custom sources)  |
| [`simulate_response_data()`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md) | Generate a binary response matrix from a calibration result               |

**Calibration algorithms:**

| Function                                                                              | Purpose                                                                |
|---------------------------------------------------------------------------------------|------------------------------------------------------------------------|
| [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md) | Algorithm 1: deterministic calibration via Brent’s root-finding method |
| [`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md) | Algorithm 2: stochastic calibration via Robbins-Monro approximation    |

**Reliability computation:**

| Function                                                                                    | Purpose                                                                 |
|---------------------------------------------------------------------------------------------|-------------------------------------------------------------------------|
| [`compute_rho_tilde()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md) | Compute the average-information reliability at a given $c$              |
| [`compute_rho_bar()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_bar.md)   | Compute the MSEM-based marginal reliability at a given $c$              |
| [`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md) | Compute both metrics in a single matrix pass (performance optimization) |
| [`compute_apc_init()`](https://joonho112.github.io/IRTsimrel/reference/compute_apc_init.md) | Analytic pre-calibration: closed-form initial value for $c$             |

**Diagnostics and screening:**

| Function                                                                                                  | Purpose                                                              |
|-----------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------|
| [`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)             | Screen whether a target reliability is achievable for a given design |
| [`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)                             | Compute and plot $\rho(c)$ over a grid of scaling factor values      |
| [`compute_reliability_tam()`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md) | Validate calibration using TAM (WLE and EAP reliability)             |

**Comparison and visualization:**

| Function                                                                                  | Purpose                                              |
|-------------------------------------------------------------------------------------------|------------------------------------------------------|
| [`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md) | Compare EQC and SAC calibration results side by side |
| [`compare_shapes()`](https://joonho112.github.io/IRTsimrel/reference/compare_shapes.md)   | Visualize and compare latent distribution shapes     |

Two deprecated aliases
([`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
and
[`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md))
are retained for backward compatibility with code written for v0.1.0.

### 6.2 S3 Methods

Every core object class has a rich set of S3 methods for inspection and
extraction. These make it easy to work with IRTsimrel results using
standard R idioms.

**`latent_G` objects** (from
[`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md)):

- [`print()`](https://rdrr.io/r/base/print.html) — concise one-line
  summary
- [`summary()`](https://rdrr.io/r/base/summary.html) — distributional
  statistics (mean, SD, skewness, kurtosis)
- `print.summary()` — formatted display of summary
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) — density
  plot of the latent distribution
- [`as.numeric()`](https://rdrr.io/r/base/numeric.html) — extract raw
  theta vector

**`item_params` objects** (from
[`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md)):

- [`print()`](https://rdrr.io/r/base/print.html) — concise parameter
  summary
- [`summary()`](https://rdrr.io/r/base/summary.html) — statistics for
  each parameter type
- `print.summary()` — formatted display
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) — scatterplot
  of difficulty vs. discrimination
- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) —
  convert to a tidy data frame

**`eqc_result` objects** (from
[`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)):

- [`print()`](https://rdrr.io/r/base/print.html) — full calibration
  report
- [`summary()`](https://rdrr.io/r/base/summary.html) — key results in
  compact form
- `print.summary()` — formatted display
- [`coef()`](https://rdrr.io/r/stats/coef.html) — extract calibrated
  item parameters as a data frame
- [`predict()`](https://rdrr.io/r/stats/predict.html) — evaluate
  reliability at new scaling factor values

**`sac_result` objects** (from
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)):

- [`print()`](https://rdrr.io/r/base/print.html) — full calibration
  report with convergence diagnostics
- [`summary()`](https://rdrr.io/r/base/summary.html) — key results in
  compact form
- `print.summary()` — formatted display
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) — convergence
  trajectory (scaling factor and reliability)
- [`coef()`](https://rdrr.io/r/stats/coef.html) — extract calibrated
  item parameters as a data frame
- [`predict()`](https://rdrr.io/r/stats/predict.html) — evaluate
  reliability at new scaling factor values

**Other classes:**

- `feasibility_check`: [`print()`](https://rdrr.io/r/base/print.html)
  displays achievable reliability ranges
- `rho_curve`: [`print()`](https://rdrr.io/r/base/print.html) displays
  curve summary and first few rows

In total, the package provides **25 S3 methods** for interactive
exploration of results.

### 6.3 New in v0.2.0

Version 0.2.0 introduces several important new features:

- **[`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)**:
  Screen whether a target reliability is achievable before investing
  computation time in calibration. Reports the achievable range for both
  metrics.

- **[`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)**:
  Visualize how reliability varies with the scaling factor over a fine
  grid. Supports both metrics, customizable grid, and optional plotting.

- **[`compute_rho_both()`](https://joonho112.github.io/IRTsimrel/reference/compute_rho_both.md)**:
  Compute both the average-information and MSEM-based reliabilities in a
  single matrix pass, avoiding redundant computation.

- **[`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) methods**:
  Standard S3 generics for extracting item parameters and evaluating the
  reliability function at new points.

- **SAC warm start from EQC**: Pass an `eqc_result` object directly to
  `sac_calibrate(c_init = eqc_result)` for fast, well-initialized
  stochastic validation.

- **Renamed from SPC to SAC**: The stochastic algorithm is now called
  SAC (Stochastic Approximation Calibration) to match the paper. The old
  name
  [`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
  is retained as a deprecated alias.

Here is
[`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
in action:

``` r
# Screen what reliabilities are achievable for this design
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

## 7. Road Map to the Vignettes

The IRTsimrel documentation is organized into two parallel tracks so you
can find the right level of detail for your needs. All vignettes are
accessible from the package website or via `vignette("name")` in your R
console.

### 7.1 Applied Researchers Track

This track focuses on **how to use** the package in practice, with
emphasis on complete worked examples, actionable guidance, and minimal
mathematical detail. If you primarily want to generate data for your
simulation study, start here.

| Vignette                                                                                                     | Title                | Est. Time | What You Will Learn                                                          |
|--------------------------------------------------------------------------------------------------------------|----------------------|-----------|------------------------------------------------------------------------------|
| [`vignette("quick-start")`](https://joonho112.github.io/IRTsimrel/articles/quick-start.md)                   | Quick Start          | 5 min     | Calibrate, generate data, and extract results in a single pass               |
| [`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)               | Applied Guide        | 30 min    | Complete applied tutorial: 2PL models, non-normal distributions, IRW items   |
| [`vignette("latent-distributions")`](https://joonho112.github.io/IRTsimrel/articles/latent-distributions.md) | Latent Distributions | 20 min    | All 12 latent shapes, when to use each, custom mixtures, visualization       |
| [`vignette("item-parameters")`](https://joonho112.github.io/IRTsimrel/articles/item-parameters.md)           | Item Parameters      | 20 min    | Parametric, IRW, hierarchical, and custom item generation, correlated params |
| [`vignette("simulation-design")`](https://joonho112.github.io/IRTsimrel/articles/simulation-design.md)       | Simulation Design    | 25 min    | Factorial designs with crossed reliability levels, looping over conditions   |
| [`vignette("case-studies")`](https://joonho112.github.io/IRTsimrel/articles/case-studies.md)                 | Case Studies         | 30 min    | Complete worked examples: DIF detection, CAT simulation, equating, MIRT      |

### 7.2 Methodological Researchers Track

This track covers the **mathematical and algorithmic foundations** for
readers who want to understand or extend what happens under the hood. If
you are developing new calibration methods or writing a methods paper,
start here.

| Vignette                                                                                                 | Title              | Est. Time | What You Will Learn                                                              |
|----------------------------------------------------------------------------------------------------------|--------------------|-----------|----------------------------------------------------------------------------------|
| [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md) | Reliability Theory | 45 min    | Population-level reliability definitions, Jensen’s inequality, metric comparison |
| [`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)           | Algorithm 1: EQC   | 25 min    | Derivation of EQC, Brent’s method, Monte Carlo error bounds, diagnostics         |
| [`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md)           | Algorithm 2: SAC   | 25 min    | Robbins-Monro theory, step-size tuning, Polyak-Ruppert averaging, convergence    |
| [`vignette("validation")`](https://joonho112.github.io/IRTsimrel/articles/validation.md)                 | Validation         | 20 min    | TAM validation, EQC-vs-SAC agreement, large-scale Monte Carlo replication        |
| [`vignette("api-reference")`](https://joonho112.github.io/IRTsimrel/articles/api-reference.md)           | API Reference      | Reference | Every exported function with full signature, arguments, return values, examples  |

### 7.3 Recommended Reading Paths

Depending on your goal, here are four suggested paths through the
documentation:

**Path 1: “I need to simulate data right now.”**

> [`vignette("quick-start")`](https://joonho112.github.io/IRTsimrel/articles/quick-start.md)
> $\rightarrow$[`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)

Start with the 5-minute quick start for instant gratification, then move
to the full applied guide for more realistic scenarios including 2PL
models, non-normal distributions, and IRW item sources.

**Path 2: “I want to understand the mathematical foundations.”**

> [`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)
> $\rightarrow$[`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)
> $\rightarrow$[`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md)

Begin with the theoretical development of population-level reliability
metrics and Jensen’s inequality, then study the derivation and
convergence properties of each calibration algorithm.

**Path 3: “I am writing a methods paper and need simulation
conditions.”**

> [`vignette("applied-guide")`](https://joonho112.github.io/IRTsimrel/articles/applied-guide.md)
> $\rightarrow$[`vignette("simulation-design")`](https://joonho112.github.io/IRTsimrel/articles/simulation-design.md)
> $\rightarrow$[`vignette("case-studies")`](https://joonho112.github.io/IRTsimrel/articles/case-studies.md)

Learn the applied workflow first, then see how to set up factorial
simulation designs with crossed reliability levels, and finally study
complete worked examples spanning DIF detection, CAT simulation,
equating, and more.

**Path 4: “I just need to look up a specific function.”**

> [`vignette("api-reference")`](https://joonho112.github.io/IRTsimrel/articles/api-reference.md)

The API reference lists every exported function with its full signature,
all arguments, return value description, and runnable code examples. Use
it as a quick-lookup reference during your workflow.

### 7.4 Cross-References Between Tracks

The two tracks are not isolated. Many applied vignettes link to their
methodological counterparts for readers who want more depth, and
methodological vignettes link back to applied examples for concreteness.
For example:

- The **Applied Guide** links to **Reliability Theory** for the
  mathematical justification of why controlling reliability matters.
- The **Algorithm EQC** vignette links to the **Quick Start** for
  minimal working examples of the function calls.
- The **Case Studies** vignette links to **Simulation Design** for the
  general framework underlying each specific case.

## 8. Installation

Install the released version from CRAN (when available):

``` r
install.packages("IRTsimrel")
```

Or install the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("joonho112/IRTsimrel")
```

For full functionality, you may also want to install the optional
dependencies:

``` r
# For TAM-based validation of calibration results
install.packages("TAM")

# For empirically-grounded item difficulties from the Item Response Warehouse
# install.packages("remotes")
remotes::install_github("itemresponsewarehouse/irw")

# For enhanced plotting (ggplot2 themes) and multi-panel layouts
install.packages(c("ggplot2", "patchwork"))
```

After installation, verify that the package loads correctly:

``` r
library(IRTsimrel)
packageVersion("IRTsimrel")
```

## 9. Getting Help

- **Package documentation:** Every exported function has a help page.
  Try
  [`?eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
  [`?sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md),
  [`?sim_latentG`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md),
  etc.
- **Vignettes:** Browse `browseVignettes("IRTsimrel")` or visit the
  [package website](https://joonho112.github.io/IRTsimrel/).
- **GitHub issues:**
  [github.com/joonho112/IRTsimrel/issues](https://github.com/joonho112/IRTsimrel/issues)
  for bug reports and feature requests.
- **Paper:** Lee (2025) provides the full mathematical development,
  proofs, and extensive simulation evidence supporting the framework.

## 10. References

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. *arXiv preprint*, arXiv:2512.16012.

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics, 22*(3), 400–407.

Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
approximation by averaging. *SIAM Journal on Control and Optimization,
30*(4), 838–855.
