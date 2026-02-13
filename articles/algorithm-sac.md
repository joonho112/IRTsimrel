# Algorithm 2: Stochastic Approximation Calibration (SAC)

``` r
library(IRTsimrel)
set.seed(42)
```

## Overview

Stochastic Approximation Calibration (SAC) is the secondary algorithm in
IRTsimrel, designed to complement EQC (Algorithm 1). SAC uses the
Robbins-Monro stochastic approximation framework to find the
discrimination scaling factor $c^{*}$ that achieves a target
reliability.

**Reading time**: approximately 25 minutes.

**Positioning**: SAC serves primarily as a validation companion to EQC.
When both algorithms agree on $c^{*}$ (typically within 5%), the
calibration can be trusted with high confidence. SAC also provides
trajectory-based convergence diagnostics that EQC lacks.

**Prerequisites**: Familiarity with
[`vignette("theory-reliability")`](https://joonho112.github.io/IRTsimrel/articles/theory-reliability.md)
for the mathematical foundations. For applied workflows, see
[`vignette("introduction")`](https://joonho112.github.io/IRTsimrel/articles/introduction.md).

**Key references**: Lee (2025, arXiv:2512.16012), Section 2.5; Robbins
and Monro (1951); Polyak and Juditsky (1992); Kushner and Yin (2003).

### When to Use SAC

| Use Case                             | Recommendation                |
|:-------------------------------------|:------------------------------|
| Routine simulation studies           | Use EQC (faster, no tuning)   |
| Independent validation of EQC        | Use SAC with EQC warm start   |
| Convergence diagnostics              | SAC provides trajectory plots |
| Research on stochastic approximation | Use SAC                       |

## The Robbins-Monro Framework

### The Stochastic Root-Finding Problem

SAC targets the same root-finding problem as EQC, but uses a
fundamentally different approach. The goal is to solve:

$$g(c) = {\mathbb{E}}\!\left\lbrack \widehat{\rho}(c) \right\rbrack - \rho^{*} = 0$$

where $g(c) = \rho(c) - \rho^{*}$ is the reliability gap function. The
key distinction from EQC is that SAC never constructs a deterministic
approximation of $g$. Instead, at each iteration $n$, it observes a
noisy evaluation:

$$Y_{n} = {\widehat{\rho}}_{n}\left( c_{n} \right) = \rho\left( c_{n} \right) + \varepsilon_{n}$$

where ${\widehat{\rho}}_{n}$ is computed from fresh Monte Carlo samples
at each step, and $\varepsilon_{n}$ is zero-mean noise with bounded
variance:
${\mathbb{E}}\left\lbrack \varepsilon_{n} \mid c_{n} \right\rbrack = 0$
and
${\mathbb{E}}\left\lbrack \varepsilon_{n}^{2} \mid c_{n} \right\rbrack \leq \sigma_{\varepsilon}^{2}$.

### The Robbins-Monro Update Rule

The iterative update proceeds as:

$$c_{n + 1} = \Pi_{\lbrack c_{L},c_{U}\rbrack}\!\left\lbrack c_{n} - a_{n}\,\left( {\widehat{\rho}}_{n}\left( c_{n} \right) - \rho^{*} \right) \right\rbrack$$

where $\Pi_{\lbrack c_{L},c_{U}\rbrack}\lbrack \cdot \rbrack$ denotes
projection onto the feasible interval
$\left\lbrack c_{L},c_{U} \right\rbrack$ and $a_{n}$ is a decreasing
step size sequence.

**Intuition**: When
${\widehat{\rho}}_{n}\left( c_{n} \right) > \rho^{*}$ (current
reliability is too high), the update decreases $c_{n}$ (reducing
discrimination, hence reliability). When
${\widehat{\rho}}_{n}\left( c_{n} \right) < \rho^{*}$, the update
increases $c_{n}$.

### Robbins-Monro Conditions

The step size sequence must satisfy:

$$\sum\limits_{n = 1}^{\infty}a_{n} = \infty\quad\text{and}\quad\sum\limits_{n = 1}^{\infty}a_{n}^{2} < \infty$$

**First condition** ($\Sigma a_{n} = \infty$): Ensures the algorithm can
reach any point in the parameter space, regardless of the initial value.
Without this, the total adjustment
$\sum a_{n} \cdot \left| {\widehat{\rho}}_{n} - \rho^{*} \right|$ might
be bounded, preventing convergence from a distant starting point.

**Second condition** ($\Sigma a_{n}^{2} < \infty$): Controls the
accumulated noise. The total noise variance is proportional to
$\sum a_{n}^{2}\sigma_{\varepsilon}^{2}$, which must be finite for
convergence.

## Step Size Analysis

### The Step Size Sequence

SAC uses the parametric family:

$$a_{n} = \frac{a}{(n + A)^{\gamma}}$$

with three parameters controlling different aspects of the learning
dynamics.

### Role of Each Parameter

**$a > 0$ (scale constant)**: Controls the overall magnitude of updates.
Larger $a$ means more aggressive steps, leading to faster initial
progress but potentially more oscillation.

In the package: `step_params$a` (default: 1).

**$A \geq 0$ (stabilization constant)**: Dampens early iterations where
step sizes would otherwise be very large ($a_{1} = a/(1 + A)^{\gamma}$).
When $A > 0$, the effective step size starts smaller, providing
stability during the initial transient phase. This is particularly
important when `c_init` is far from $c^{*}$.

The package sets $A$ automatically based on the burn-in fraction.

**$\gamma \in (1/2,1\rbrack$ (decay exponent)**: Controls how quickly
step sizes decrease. In the package: `step_params$gamma` (default:
0.67).

- $\gamma = 0.51$: Slow decay, close to the boundary of the
  Robbins-Monro conditions. Maintains larger steps for longer, which can
  be beneficial when initial estimates are poor.
- $\gamma = 2/3$: A common compromise between convergence speed and
  stability.
- $\gamma = 1$: Fastest allowable decay. Reduces variance quickly but
  may under-correct if the initial value is far from the solution.

### Robbins-Monro Conditions Verification

For the parametric family $a_{n} = a/(n + A)^{\gamma}$:

$$\left. \sum\limits_{n = 1}^{\infty}\frac{a}{(n + A)^{\gamma}} = \infty\quad\Leftrightarrow\quad\gamma \leq 1 \right.$$

$$\left. \sum\limits_{n = 1}^{\infty}\frac{a^{2}}{(n + A)^{2\gamma}} < \infty\quad\Leftrightarrow\quad 2\gamma > 1\;\Leftrightarrow\;\gamma > 1/2 \right.$$

Therefore, any $\gamma \in (1/2,1\rbrack$ satisfies both conditions.

### Step Size Sensitivity

The effect of different step size configurations can be observed through
the convergence trajectory.

``` r
# Common setup
eqc_ref <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", M = 5000L, seed = 42, verbose = FALSE
)

# Run SAC with different step decay values
decays <- c(0.51, 0.67, 0.85, 1.0)
colors <- c("steelblue", "coral", "seagreen", "purple")

oldpar <- par(mar = c(4.5, 4.5, 3, 1))
on.exit(par(oldpar))

plot(NULL, xlim = c(1, 300), ylim = c(0.5, 2.5),
     xlab = "Iteration", ylab = "c_n",
     main = "SAC Trajectories: Step Decay Sensitivity")
abline(h = eqc_ref$c_star, lty = 2, col = "gray40", lwd = 1.5)

for (k in seq_along(decays)) {
  sac_k <- sac_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    c_init = 0.5,  # deliberately poor start
    n_iter = 300L, M_per_iter = 1000L,
    step_params = list(gamma = decays[k]),
    seed = 42, verbose = FALSE
  )
  lines(seq_along(sac_k$trajectory), sac_k$trajectory,
        col = colors[k], lwd = 1.5)
}

legend("topright",
       legend = c(sprintf("gamma = %.2f", decays), "EQC c*"),
       col = c(colors, "gray40"), lty = c(rep(1, 4), 2),
       lwd = c(rep(1.5, 4), 1.5), cex = 0.85)
```

![Figure 1: Effect of step decay parameter on SAC convergence
trajectories.](algorithm-sac_files/figure-html/step-size-sensitivity-1.png)

Figure 1: Effect of step decay parameter on SAC convergence
trajectories.

**Observations**:

- Slower decay ($\gamma = 0.51$) maintains larger oscillations longer
  but eventually converges.
- Faster decay ($\gamma = 1.0$) converges quickly when starting near the
  solution, but risks under-correcting from a distant start.
- The default $\gamma = 0.67$ provides a good compromise between speed
  and stability.

## Polyak-Ruppert Averaging

### Why Averaging?

The raw SAC iterates $c_{n}$ converge at rate
$O\left( n^{- \gamma} \right)$, which is suboptimal. Polyak and Juditsky
(1992) showed that a simple average of the iterates achieves the optimal
rate $O\left( n^{- 1/2} \right)$, regardless of $\gamma$ (as long as the
Robbins-Monro conditions hold).

### The Averaging Formula

After discarding a burn-in period of $B$ iterations, the Polyak-Ruppert
average is:

$${\bar{c}}_{N} = \frac{1}{N - B}\sum\limits_{n = B + 1}^{N}c_{n}$$

where $N$ is the total number of iterations and $B$ is the burn-in
count.

### Convergence Rate Comparison

| Estimator                      | Rate                           | Notes                            |
|:-------------------------------|:-------------------------------|:---------------------------------|
| Final iterate $c_{N}$          | $O\left( N^{- \gamma} \right)$ | Depends on step decay            |
| Polyak-Ruppert ${\bar{c}}_{N}$ | $O\left( N^{- 1/2} \right)$    | Optimal, independent of $\gamma$ |
| Truncated average              | $O\left( N^{- 1/2} \right)$    | Requires $B = o(N)$              |

The practical benefit is substantial. With $\gamma = 0.67$ (the package
default), the final iterate converges at rate
$O\left( n^{- 0.67} \right)$, while the Polyak-Ruppert average converges
at rate $O\left( n^{- 0.5} \right)$ — nearly the same exponent but with
considerably smaller constants due to variance reduction.

### Burn-In Selection

The burn-in period $B$ should be large enough to exclude the initial
transient (before the iterates have approached the solution) but small
enough to include a sufficient number of post-burn-in samples.

In the package: `burn_in` parameter, specified as an integer count
(default: `floor(n_iter / 2)`, i.e., 50% of iterations).

**Guidelines**:

- With EQC warm start: the default `burn_in = floor(n_iter / 2)` is
  sufficient (starting near the solution, transient is short).
- With cold start or APC: the default 50% burn-in is appropriate; for
  very long runs ($N > 1000$), a smaller fraction (e.g.,
  `burn_in = floor(n_iter / 5)`) can be used to include more
  post-burn-in samples.

## Initialization Strategies

The initial value $c_{0}$ significantly affects SAC’s convergence speed.
The package supports three initialization methods.

### Strategy 1: EQC Warm Start (Recommended)

Pass an `eqc_result` object directly to `c_init`:

``` r
eqc_for_init <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", M = 5000L, seed = 42, verbose = FALSE
)

sac_warm <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_for_init,  # EQC warm start
  n_iter = 200L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

cat(sprintf("EQC c*: %.4f\n", eqc_for_init$c_star))
#> EQC c*: 0.8995
cat(sprintf("SAC c* (warm start): %.4f\n", sac_warm$c_star))
#> SAC c* (warm start): 0.9113
cat(sprintf("Init method: %s\n", sac_warm$init_method))
#> Init method: eqc_warm_start
```

**Advantages**: Starting near the solution means the transient phase is
short, burn-in can be minimal, and fewer iterations are needed.

### Strategy 2: Analytic Pre-Calibration (APC)

When `c_init = NULL`, SAC uses a closed-form approximation derived under
Gaussian Rasch assumptions.

**Definition (APC Formula).** Under the assumptions $\theta \sim N(0,1)$
and $\beta_{i} \sim N\left( 0,\sigma_{\beta}^{2} \right)$ with
$\lambda_{i}^{(0)} = 1$ (Rasch), the approximate initial value is:

$$c_{\text{APC}} = \sqrt{\frac{\rho^{*}}{I \cdot \kappa\left( \sigma^{2} \right) \cdot \left( 1 - \rho^{*} \right)}}$$

where $\sigma^{2} = 1 + \sigma_{\beta}^{2}$ and the logistic-normal
convolution constant is approximated by:

$$\kappa\left( \sigma^{2} \right) \approx \frac{0.25}{\sqrt{1 + \sigma^{2}\pi^{2}/3}}$$

**Derivation sketch**: Under the Gaussian Rasch model, the expected item
information at a single item is
${\mathbb{E}}\left\lbrack \lambda^{2}p(1 - p) \right\rbrack = c^{2}\kappa\left( \sigma^{2} \right)$.
With $I$ items, the average test information is approximately
$Ic^{2}\kappa$, so
$\widetilde{\rho} \approx \sigma_{\theta}^{2}Ic^{2}\kappa/\left( \sigma_{\theta}^{2}Ic^{2}\kappa + 1 \right)$.
Setting this equal to $\rho^{*}$ and solving for $c$ yields Eq. (7).

``` r
# APC initialization (c_init = NULL triggers APC)
sac_apc <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = NULL,  # APC warm start
  n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

cat(sprintf("SAC c* (APC start): %.4f\n", sac_apc$c_star))
#> SAC c* (APC start): 0.9793
cat(sprintf("Init method: %s\n", sac_apc$init_method))
#> Init method: apc_warm_start

# Compare APC init quality
c_apc <- compute_apc_init(target_rho = 0.80, n_items = 25)
cat(sprintf("APC initial value: %.4f (vs true c* ~ %.4f)\n",
            c_apc, eqc_for_init$c_star))
#> APC initial value: 1.3274 (vs true c* ~ 0.8995)
```

### Strategy 3: User-Specified Value

A numeric value can be passed directly:

``` r
sac_user <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = 1.0,  # user-specified
  n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

cat(sprintf("SAC c* (user init = 1.0): %.4f\n", sac_user$c_star))
#> SAC c* (user init = 1.0): 0.9233
cat(sprintf("Init method: %s\n", sac_user$init_method))
#> Init method: user_specified
```

### Initialization Comparison

``` r
oldpar <- par(mar = c(4.5, 4.5, 3, 1))
on.exit(par(oldpar))

plot(NULL, xlim = c(1, 300), ylim = c(0, 3),
     xlab = "Iteration", ylab = "c_n",
     main = "Effect of Initialization on Convergence")
abline(h = eqc_for_init$c_star, lty = 2, col = "gray40", lwd = 1.5)

lines(seq_along(sac_warm$trajectory), sac_warm$trajectory,
      col = "steelblue", lwd = 1.5)
lines(seq_along(sac_apc$trajectory), sac_apc$trajectory,
      col = "coral", lwd = 1.5)
lines(seq_along(sac_user$trajectory), sac_user$trajectory,
      col = "seagreen", lwd = 1.5)

legend("topright",
       legend = c("EQC warm start", "APC", "User (c=1.0)", "EQC c*"),
       col = c("steelblue", "coral", "seagreen", "gray40"),
       lty = c(1, 1, 1, 2), lwd = 1.5, cex = 0.9)
```

![Figure 2: Convergence speed depends strongly on initialization
quality.](algorithm-sac_files/figure-html/init-comparison-1.png)

Figure 2: Convergence speed depends strongly on initialization quality.

## Convergence Theory

### Almost Sure Convergence (Theorem A.3)

**Theorem (SAC Convergence; cf. Theorem A.3, Lee 2025).** Suppose:

1.  The step sizes satisfy Eq. (4): $\sum a_{n} = \infty$,
    $\sum a_{n}^{2} < \infty$.
2.  The reliability function $\rho(c)$ is Lipschitz continuous on
    $\left\lbrack c_{L},c_{U} \right\rbrack$.
3.  $g(c) = \rho(c) - \rho^{*}$ has a unique root $c^{*}$ in
    $\left( c_{L},c_{U} \right)$.
4.  $g\prime\left( c^{*} \right) > 0$ (strict monotonicity at the root).
5.  The noise $\varepsilon_{n}$ satisfies
    ${\mathbb{E}}\left\lbrack \varepsilon_{n} \mid c_{n} \right\rbrack = 0$
    and
    ${\mathbb{E}}\left\lbrack \varepsilon_{n}^{2} \mid c_{n} \right\rbrack \leq \sigma_{\varepsilon}^{2} < \infty$.

Then the projected SAC iterates converge almost surely:

$$\left. c_{n}\overset{\text{a.s.}}{\rightarrow}c^{*}\quad{\text{as}\mspace{6mu}}n\rightarrow\infty \right.$$

**Proof sketch.** The proof follows the classical Robbins-Monro theory,
extended to the projected (constrained) setting by Kushner and Yin
(2003, Chapter 5). The key steps are:

1.  **Lyapunov function**: Define $V(c) = \left( c - c^{*} \right)^{2}$.
    Compute the conditional expectation of $V\left( c_{n + 1} \right)$
    given $c_{n}$: $$\begin{aligned}
    {{\mathbb{E}}\left\lbrack V\left( c_{n + 1} \right) \mid c_{n} \right\rbrack} & {= V\left( c_{n} \right) - 2a_{n}\left( c_{n} - c^{*} \right)g\left( c_{n} \right) + a_{n}^{2}{\mathbb{E}}\left\lbrack \left( {\widehat{\rho}}_{n} - \rho^{*} \right)^{2} \mid c_{n} \right\rbrack}
    \end{aligned}$$

2.  **Drift condition**: Since $g(c)\left( c - c^{*} \right) > 0$ for
    $c \neq c^{*}$ (the reliability function is increasing, so $g$ has
    the same sign as $c - c^{*}$), the drift term
    $- 2a_{n}\left( c_{n} - c^{*} \right)g\left( c_{n} \right) < 0$
    provides contraction toward $c^{*}$.

3.  **Noise control**: The accumulated noise
    $\sum a_{n}^{2}\sigma_{\varepsilon}^{2} < \infty$ (by the step size
    conditions), so the noise does not prevent convergence.

4.  **Almost sure convergence**: By the Robbins-Siegmund theorem (a
    supermartingale convergence result),
    $\left. V\left( c_{n} \right)\rightarrow V^{*} \right.$ a.s. for
    some random variable $V^{*}$. The drift condition then forces
    $V^{*} = 0$, giving $\left. c_{n}\rightarrow c^{*} \right.$ a.s.

### Asymptotic Normality of Polyak-Ruppert Average

**Theorem (Polyak-Ruppert CLT).** Under the conditions of Theorem A.3
and additional smoothness assumptions, the Polyak-Ruppert average
satisfies:

$$\sqrt{N - B}\,\left( {\bar{c}}_{N} - c^{*} \right)\overset{d}{\rightarrow}N\!\left( 0,\,\frac{\sigma_{\varepsilon}^{2}}{g\prime\left( c^{*} \right)^{2}} \right)$$

This achieves the Cramer-Rao lower bound for the stochastic root-finding
problem, meaning Polyak-Ruppert averaging is *asymptotically efficient*.

**Rate comparison**: For $\gamma = 0.67$ (the package default):

- Final iterate: MSE $\sim O\left( n^{- 1.34} \right)$
- Polyak-Ruppert average: MSE $\sim O\left( n^{- 1} \right)$

The final iterate converges faster in rate than the PR average in this
case, but the PR average achieves smaller constants and optimal
asymptotic efficiency. For $\gamma$ closer to $1/2$, the rate
improvement from averaging is more dramatic.

### Convergence Diagnostics

SAC provides automatic convergence assessment through several
statistics.

``` r
# Run SAC with enough iterations
sac_diag <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_for_init, n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

conv <- sac_diag$convergence
cat("Convergence diagnostics:\n")
#> Convergence diagnostics:
cat(sprintf("  Converged:          %s\n", conv$converged))
#>   Converged:          TRUE
cat(sprintf("  Mean (first half):  %.4f\n", conv$mean_first_half))
#>   Mean (first half):  0.9117
cat(sprintf("  Mean (second half): %.4f\n", conv$mean_second_half))
#>   Mean (second half): 0.9126
cat(sprintf("  SD (post-burn-in):  %.4f\n", conv$sd_post_burn))
#>   SD (post-burn-in):  0.0007
cat(sprintf("  Hit lower bound:    %s\n", conv$hit_lower_bound))
#>   Hit lower bound:    FALSE
cat(sprintf("  Hit upper bound:    %s\n", conv$hit_upper_bound))
#>   Hit upper bound:    FALSE
```

**Convergence criterion**: The algorithm is considered converged when:

$$\left| {\bar{c}}_{\text{first half}} - {\bar{c}}_{\text{second half}} \right| < 0.05$$

This mean-split diagnostic detects systematic drift, which would
indicate that the iterates have not yet settled around the solution.

**Post-burn-in SD**: Should be small relative to $c^{*}$. Values above
0.1 suggest more iterations or better initialization.

## Comparative Analysis: EQC vs SAC

### When They Agree

When both algorithms use the same reliability metric and the SAC run has
converged, agreement within 5% is typical and expected.

``` r
# Both targeting rho* = 0.80 with info metric
eqc_comp <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  M = 5000L, seed = 42, verbose = FALSE
)

sac_comp <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_comp, n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

comparison <- compare_eqc_sac(eqc_comp, sac_comp, verbose = FALSE)
cat(sprintf("EQC c*:   %.4f\n", comparison$c_eqc))
#> EQC c*:   0.8995
cat(sprintf("SAC c*:   %.4f\n", comparison$c_sac))
#> SAC c*:   0.9122
cat(sprintf("Diff:     %.4f (%.2f%%)\n", comparison$diff_abs, comparison$diff_pct))
#> Diff:     0.0127 (1.41%)
cat(sprintf("Agree:    %s\n", comparison$agreement))
#> Agree:    TRUE
```

### Sources of Disagreement

When EQC and SAC disagree by more than 5%, investigate the following:

1.  **Different reliability metrics**: Ensure both use `"info"` or both
    use `"msem"`. Different metrics target different quantities, so
    disagreement is expected.

2.  **SAC non-convergence**: Check `sac_result$convergence$converged`.
    If `FALSE`, increase `n_iter` or use EQC warm start.

3.  **Small M in EQC**: Low quadrature size introduces Monte Carlo
    error. Increase `M` in EQC.

4.  **Extreme target**: Very high ($> 0.95$) or very low ($< 0.50$)
    targets have flatter reliability curves, making both algorithms less
    precise.

### What to Do When They Disagree

``` r
# Strategy 1: Increase SAC iterations
sac_longer <- sac_calibrate(
  ..., n_iter = 500L, seed = 42
)

# Strategy 2: Increase EQC quadrature
eqc_precise <- eqc_calibrate(
  ..., M = 50000L, seed = 42
)

# Strategy 3: Verify same metric
stopifnot(eqc_result$metric == sac_result$metric)

# Strategy 4: Run multiple SAC replications
sac_reps <- replicate(5, {
  sac_calibrate(
    ..., seed = sample.int(10000, 1), verbose = FALSE
  )$c_star
})
cat(sprintf("SAC c* across 5 reps: mean = %.4f, sd = %.4f\n",
            mean(sac_reps), sd(sac_reps)))
```

## Diagnostic Examples

### Trajectory Analysis

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
provides visual diagnostics for SAC convergence.

``` r
# Run SAC from APC start for a more interesting trajectory
sac_traj <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = NULL,  # APC start
  n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

plot(sac_traj, type = "c")
```

![Figure 3: SAC convergence trajectory showing the scaling factor
iterations.](algorithm-sac_files/figure-html/trajectory-plot-1.png)

Figure 3: SAC convergence trajectory showing the scaling factor
iterations.

``` r
plot(sac_traj, type = "rho")
```

![Figure 4: SAC reliability trajectory showing convergence to the
target.](algorithm-sac_files/figure-html/rho-trajectory-plot-1.png)

Figure 4: SAC reliability trajectory showing convergence to the target.

### Replication Variability

SAC is inherently stochastic. Running multiple replications reveals the
variability of the estimate.

``` r
n_reps <- 8
sac_reps <- numeric(n_reps)

for (r in seq_len(n_reps)) {
  sac_r <- sac_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    c_init = eqc_for_init,
    n_iter = 200L, M_per_iter = 1000L,
    seed = r * 100, verbose = FALSE
  )
  sac_reps[r] <- sac_r$c_star
}

cat("SAC replication study (8 runs with EQC warm start):\n")
#> SAC replication study (8 runs with EQC warm start):
cat(sprintf("  Mean c*:  %.4f\n", mean(sac_reps)))
#>   Mean c*:  0.9160
cat(sprintf("  SD c*:    %.4f\n", sd(sac_reps)))
#>   SD c*:    0.0075
cat(sprintf("  Range:    [%.4f, %.4f]\n", min(sac_reps), max(sac_reps)))
#>   Range:    [0.9056, 0.9247]
cat(sprintf("  EQC c*:   %.4f\n", eqc_for_init$c_star))
#>   EQC c*:   0.8995
cat(sprintf("  |Mean - EQC|: %.4f\n", abs(mean(sac_reps) - eqc_for_init$c_star)))
#>   |Mean - EQC|: 0.0165
```

### Comparison Across Different Conditions

``` r
# Test across different target reliabilities
targets <- c(0.60, 0.70, 0.80, 0.90)
cat("EQC vs SAC across target reliabilities (25 Rasch items):\n")
#> EQC vs SAC across target reliabilities (25 Rasch items):
cat(sprintf("  %-8s %-10s %-10s %-10s %-8s\n",
            "Target", "EQC c*", "SAC c*", "|Diff|", "Agree"))
#>   Target   EQC c*     SAC c*     |Diff|     Agree

for (rho_t in targets) {
  eqc_t <- eqc_calibrate(
    target_rho = rho_t, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    M = 5000L, seed = 42, verbose = FALSE
  )
  sac_t <- sac_calibrate(
    target_rho = rho_t, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    c_init = eqc_t, n_iter = 200L, M_per_iter = 1000L,
    seed = 42, verbose = FALSE
  )
  diff_abs <- abs(eqc_t$c_star - sac_t$c_star)
  diff_pct <- 100 * diff_abs / eqc_t$c_star
  agree <- diff_pct < 5
  cat(sprintf("  %-8.2f %-10.4f %-10.4f %-10.4f %-8s\n",
              rho_t, eqc_t$c_star, sac_t$c_star, diff_abs,
              ifelse(agree, "YES", "NO")))
}
#>   0.60     0.5118     0.5145     0.0026     YES     
#>   0.70     0.6548     0.6603     0.0055     YES     
#>   0.80     0.8995     0.9113     0.0118     YES     
#>   0.90     1.5328     1.5537     0.0209     YES
```

### Effect of M_per_iter

The per-iteration Monte Carlo sample size controls the noise level in
each SAC step.

``` r
m_values <- c(200, 500, 1000, 2000)
cat("SAC sensitivity to M_per_iter (target = 0.80, 300 iter):\n")
#> SAC sensitivity to M_per_iter (target = 0.80, 300 iter):
cat(sprintf("  %-12s %-10s %-12s %-10s\n",
            "M_per_iter", "c*", "Post-burn SD", "Converged"))
#>   M_per_iter   c*         Post-burn SD Converged

for (m_val in m_values) {
  sac_m <- sac_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", reliability_metric = "info",
    c_init = eqc_for_init, n_iter = 300L,
    M_per_iter = as.integer(m_val),
    seed = 42, verbose = FALSE
  )
  cat(sprintf("  %-12d %-10.4f %-12.4f %-10s\n",
              m_val, sac_m$c_star,
              sac_m$convergence$sd_post_burn,
              sac_m$convergence$converged))
}
#>   200          0.9115     0.0005       TRUE      
#>   500          0.9115     0.0004       TRUE      
#>   1000         0.9122     0.0007       TRUE      
#>   2000         0.9127     0.0005       TRUE
```

## Deprecated Alias

The function was previously named
[`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).
For backward compatibility,
[`spc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
still works but is deprecated and will emit a warning. Use
[`sac_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)
in all new code.

Similarly,
[`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md)
is deprecated in favor of
[`compare_eqc_sac()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_sac.md).

## Summary

SAC is a theoretically grounded validation companion to EQC:

| Aspect                | Summary                                                  |
|:----------------------|:---------------------------------------------------------|
| **Method**            | Robbins-Monro stochastic approximation                   |
| **Convergence**       | Almost sure (Theorem A.3)                                |
| **Optimal rate**      | $O\left( n^{- 1/2} \right)$ via Polyak-Ruppert averaging |
| **Key advantage**     | Independent validation of EQC                            |
| **Key disadvantage**  | Requires tuning (step size, iterations, burn-in)         |
| **Default metric**    | Average-information (`"info"`)                           |
| **Best practice**     | Use EQC warm start, $n_{\text{iter}} \geq 200$           |
| **Convergence check** | `sac_result$convergence$converged`                       |

## References

Lee, J. (2025). Reliability-targeted simulation of item response data:
Solving the inverse design problem. *arXiv preprint arXiv:2512.16012*.

Robbins, H., & Monro, S. (1951). A stochastic approximation method. *The
Annals of Mathematical Statistics, 22*(3), 400–407.

Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
approximation by averaging. *SIAM Journal on Control and Optimization,
30*(4), 838–855.

Kushner, H. J., & Yin, G. G. (2003). *Stochastic Approximation and
Recursive Algorithms and Applications* (2nd ed.). Springer.

Spall, J. C. (2003). *Introduction to Stochastic Search and
Optimization*. Wiley.
