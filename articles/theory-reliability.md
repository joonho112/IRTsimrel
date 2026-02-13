# Mathematical Foundations: Reliability Theory and the Inverse Design Problem

``` r
library(IRTsimrel)
set.seed(42)
```

## Overview

This vignette provides the mathematical foundations for the
reliability-targeted IRT simulation framework implemented in the
IRTsimrel package. It develops the theory from first principles,
connecting classical IRT measurement theory to the inverse calibration
problem solved by the EQC and SAC algorithms.

**Reading time**: approximately 45–60 minutes.

**Prerequisites**: Familiarity with basic probability theory
(expectation, variance, convergence), linear algebra, and the logistic
function. Some exposure to item response theory is helpful but not
strictly required.

**Notation conventions**: Random variables are denoted by Greek or
capital Roman letters; vectors are boldfaced; hats denote estimators;
bars and tildes denote specific reliability definitions (defined below).

**Cross-references**: For algorithm-specific details, see
[`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md)
and
[`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md).
For applied examples and workflows, see
[`vignette("introduction")`](https://joonho112.github.io/IRTsimrel/articles/introduction.md).

## The Two-Parameter Logistic IRT Model

We work within the two-parameter logistic (2PL) item response theory
model. A test consists of $I$ items administered to examinees drawn from
a latent trait distribution $G$.

**Definition 1 (2PL Model).** Let $\theta \sim G$ denote a latent
ability parameter and let $\left( \beta_{i},\lambda_{i} \right)$ denote
the difficulty and discrimination parameters for item $i$. The
probability of a correct response is given by the item characteristic
curve (ICC):

$$P\left( X_{i} = 1 \mid \theta \right) = \psi(\lambda_{i}\left( \theta - \beta_{i} \right))$$

where $\psi(z) = \left( 1 + e^{- z} \right)^{- 1}$ is the standard
logistic function.

The Rasch model is the special case where $\lambda_{i} = 1$ for all
items $i$.

### Item and Test Information

The Fisher information provided by item $i$ at ability level $\theta$
is:

**Definition 2 (Item Information).** Under the 2PL model,

$$\mathcal{J}_{i}(\theta) = \lambda_{i}^{2}\, p_{i}(\theta)\,(1 - p_{i}(\theta))$$

where $p_{i}(\theta) = P\left( X_{i} = 1 \mid \theta \right)$ is given
by Eq. (1).

Since $p(1 - p)$ achieves its maximum of $1/4$ when $p = 1/2$ (i.e.,
when $\theta = \beta_{i}$), the maximum information from item $i$ is
$\lambda_{i}^{2}/4$.

**Definition 3 (Test Information).** The test information function is
the sum of item information functions:

$$\mathcal{J}(\theta) = \sum\limits_{i = 1}^{I}\mathcal{J}_{i}(\theta) = \sum\limits_{i = 1}^{I}\lambda_{i}^{2}\, p_{i}(\theta)\,(1 - p_{i}(\theta))$$

Test information measures the precision of ability estimation at a given
$\theta$. Higher information means smaller standard errors and more
precise measurement.

### The Conditional Standard Error

The asymptotic standard error of the maximum likelihood estimator
$\widehat{\theta}$ at $\theta$ is:

$$\text{SE}\left( \widehat{\theta} \mid \theta \right) = \frac{1}{\sqrt{\mathcal{J}(\theta)}}$$

This conditional error varies across the ability continuum, motivating
the need for *marginal* (population-averaged) reliability measures.

## Two Definitions of Marginal Reliability

Classical test theory defines reliability as the ratio of true-score
variance to observed-score variance. In IRT, multiple marginalizations
are possible, leading to distinct reliability definitions.

Let $\sigma_{\theta}^{2} = \text{Var}_{G}(\theta)$ denote the latent
trait variance under the population distribution $G$.

### Average-Information Reliability ($\widetilde{\rho}$)

**Definition 4 (Average-Information Reliability).** The
average-information marginal reliability is defined as:

$$\widetilde{\rho}(c) = \frac{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c)}{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c) + 1}$$

where
$\bar{\mathcal{J}}(c) = {\mathbb{E}}_{G}\left\lbrack \mathcal{J}(\theta;c) \right\rbrack$
is the average test information under $G$.

This definition first averages the test information, then converts to a
reliability-like coefficient. It can be interpreted as the reliability
one would obtain if the *average* measurement precision applied
uniformly across all examinees.

In the package, this metric is selected with
`reliability_metric = "info"` (or the synonym `"tilde"`).

### MSEM-Based Reliability ($\bar{w}$)

**Definition 5 (MSEM-Based Reliability).** The
marginal-standard-error-of- measurement (MSEM) based reliability is
defined as:

$$\bar{w}(c) = \frac{\sigma_{\theta}^{2}}{\sigma_{\theta}^{2} + {\mathbb{E}}_{G}\lbrack 1/\mathcal{J}(\theta;c)\rbrack}$$

This definition first takes the reciprocal of test information at each
$\theta$ (obtaining the conditional error variance), averages these
error variances, and then forms the reliability ratio. It is the
theoretically exact marginal reliability under the IRT model.

In the package, this metric is selected with
`reliability_metric = "msem"` (or the synonym `"bar"`).

### Jensen’s Inequality: $\widetilde{\rho} \geq \bar{w}$

The relationship between the two metrics is governed by Jensen’s
inequality.

**Theorem 1 (Jensen’s Inequality for Reliability Metrics).** For any
scaling factor $c > 0$ and any latent distribution $G$ with
$\sigma_{\theta}^{2} > 0$,

$$\widetilde{\rho}(c) \geq \bar{w}(c)$$

with equality if and only if $\mathcal{J}(\theta;c)$ is constant
$G$-almost surely.

**Proof.** The function $f(x) = 1/x$ is strictly convex on $(0,\infty)$.
By Jensen’s inequality applied to the random variable
$\mathcal{J}(\theta;c) > 0$:

$${\mathbb{E}}_{G}\!\left\lbrack \frac{1}{\mathcal{J}(\theta;c)} \right\rbrack \geq \frac{1}{{\mathbb{E}}_{G}\lbrack\mathcal{J}(\theta;c)\rbrack} = \frac{1}{\bar{\mathcal{J}}(c)}$$

Now compare the denominators of $\widetilde{\rho}$ and $\bar{w}$. The
denominator of $\widetilde{\rho}$ in Eq. (5) can be rewritten as:

$$\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c) + 1 = \sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c) + 1$$

while the denominator of $\bar{w}$ in Eq. (6) is:

$$\sigma_{\theta}^{2} + {\mathbb{E}}_{G}\!\left\lbrack \frac{1}{\mathcal{J}(\theta;c)} \right\rbrack$$

From Eq. (8):

$$\begin{aligned}
{\bar{w}(c)} & {= \frac{\sigma_{\theta}^{2}}{\sigma_{\theta}^{2} + {\mathbb{E}}\left\lbrack 1/\mathcal{J}(\theta;c) \right\rbrack} \leq \frac{\sigma_{\theta}^{2}}{\sigma_{\theta}^{2} + 1/\bar{\mathcal{J}}(c)}} \\
 & {= \frac{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c)}{\sigma_{\theta}^{2}\,\bar{\mathcal{J}}(c) + 1} = \widetilde{\rho}(c)}
\end{aligned}$$

where the inequality follows from
${\mathbb{E}}\lbrack 1/\mathcal{J}\rbrack \geq 1/\bar{\mathcal{J}}$
(Jensen) and the monotone decreasing nature of
$\left. x\mapsto\sigma_{\theta}^{2}/\left( \sigma_{\theta}^{2} + x \right) \right.$.
Equality holds iff $\mathcal{J}(\theta;c)$ is a.s. constant, i.e., when
information does not vary across $\theta$. $▫$

### Numerical Verification of Jensen’s Inequality

We can verify this inequality computationally using the package’s
low-level reliability functions.

``` r
# Generate a latent trait sample and item parameters
theta <- sim_latentG(5000, shape = "normal")$theta
items <- sim_item_params(25, model = "rasch", source = "parametric")
beta  <- items$data$beta
lambda_base <- rep(1, 25)

# Compute both metrics across a range of c values
c_grid <- seq(0.3, 3.0, by = 0.1)
results <- data.frame(
  c = c_grid,
  rho_tilde = numeric(length(c_grid)),
  rho_bar   = numeric(length(c_grid)),
  gap       = numeric(length(c_grid))
)

for (j in seq_along(c_grid)) {
  both <- compute_rho_both(c_grid[j], theta, beta, lambda_base)
  results$rho_tilde[j] <- both$rho_tilde
  results$rho_bar[j]   <- both$rho_bar
  results$gap[j]       <- both$rho_tilde - both$rho_bar
}

# Verify: gap is always non-negative
cat(sprintf("Jensen's gap range: [%.6f, %.6f]\n",
            min(results$gap), max(results$gap)))
#> Jensen's gap range: [0.000209, 0.042706]
cat(sprintf("All gaps >= 0: %s\n", all(results$gap >= -1e-12)))
#> All gaps >= 0: TRUE
```

``` r
oldpar <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
on.exit(par(oldpar))

# Panel 1: Both metrics
plot(results$c, results$rho_tilde, type = "l", col = "steelblue", lwd = 2,
     xlab = "Scaling factor c", ylab = "Reliability",
     main = "Two Reliability Metrics", ylim = c(0, 1))
lines(results$c, results$rho_bar, col = "coral", lwd = 2, lty = 2)
legend("bottomright",
       legend = c(expression(tilde(rho) ~ "(info)"),
                  expression(bar(w) ~ "(msem)")),
       col = c("steelblue", "coral"), lty = c(1, 2), lwd = 2, cex = 0.9)
```

![Figure 1: Jensen's gap between the two reliability metrics across
scaling
factors.](theory-reliability_files/figure-html/jensen-plot-1.png)

Figure 1: Jensen’s gap between the two reliability metrics across
scaling factors.

``` r

# Panel 2: The gap
plot(results$c, results$gap, type = "l", col = "darkred", lwd = 2,
     xlab = "Scaling factor c", ylab = expression(tilde(rho) - bar(w)),
     main = "Jensen's Gap")
abline(h = 0, lty = 3, col = "gray50")
```

![Figure 1: Jensen's gap between the two reliability metrics across
scaling
factors.](theory-reliability_files/figure-html/jensen-plot-2.png)

Figure 1: Jensen’s gap between the two reliability metrics across
scaling factors.

The gap is largest at intermediate $c$ values where information varies
most across the ability continuum, and shrinks toward zero for very
small $c$ (where both metrics approach zero) and very large $c$ (where
information becomes uniformly high).

### Scale Ordering Corollary

**Corollary 1 (Scale Ordering).** Let $c_{\widetilde{\rho}}^{*}$ and
$c_{\bar{w}}^{*}$ denote the calibrated scaling factors for the same
target $\rho^{*}$ under the average-information and MSEM metrics,
respectively. Then:

$$c_{\bar{w}}^{*} \geq c_{\widetilde{\rho}}^{*}$$

**Proof sketch.** Since $\widetilde{\rho}(c) \geq \bar{w}(c)$ for all
$c$, the MSEM metric reaches any given target level at a larger scaling
factor than the information metric does. Formally, if
$\widetilde{\rho}\left( c_{\widetilde{\rho}}^{*} \right) = \rho^{*}$,
then
$\bar{w}\left( c_{\widetilde{\rho}}^{*} \right) \leq \widetilde{\rho}\left( c_{\widetilde{\rho}}^{*} \right) = \rho^{*}$,
so $\bar{w}$ must increase further (requiring larger $c$) to reach
$\rho^{*}$. $▫$

This ordering is demonstrated numerically below.

``` r
# Calibrate with both metrics
eqc_info <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  M = 5000L, seed = 42
)

eqc_msem <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "msem",
  M = 5000L, seed = 42
)
#> Warning in eqc_calibrate(target_rho = 0.8, n_items = 25, model = "rasch", :
#> MSEM-based reliability (w-bar) may have a non-monotone objective for EQC (see
#> Lee, 2025, Section 4.3). Consider using 'info' (rho-tilde) for EQC, or use
#> sac_calibrate() which handles w-bar targeting correctly.

cat(sprintf("Target rho* = 0.80\n"))
#> Target rho* = 0.80
cat(sprintf("  c* (info metric): %.4f\n", eqc_info$c_star))
#>   c* (info metric): 0.8995
cat(sprintf("  c* (msem metric): %.4f\n", eqc_msem$c_star))
#>   c* (msem metric): 0.9232
cat(sprintf("  c*_msem >= c*_info: %s\n",
            eqc_msem$c_star >= eqc_info$c_star - 1e-6))
#>   c*_msem >= c*_info: TRUE
```

## Global Discrimination Scaling

The central design variable in IRTsimrel is the global discrimination
scaling factor $c > 0$.

### The Scaling Model

**Definition 6 (Global Discrimination Scaling).** Given baseline item
discriminations
${\mathbf{λ}}^{(0)} = \left( \lambda_{1}^{(0)},\ldots,\lambda_{I}^{(0)} \right)$
generated from an item parameter distribution $H$, the scaled
discriminations are:

$$\lambda_{i}(c) = c \cdot \lambda_{i}^{(0)},\quad i = 1,\ldots,I$$

This scaling has several desirable properties:

1.  **Separation of structure and informativeness**: Changing $c$
    modifies test information without altering the relative difficulty
    or discrimination structure among items.

2.  **Dimensionality reduction**: The multi-dimensional item parameter
    space is collapsed to a single scalar control variable $c$.

3.  **Smooth dependence**: Both $\widetilde{\rho}(c)$ and $\bar{w}(c)$
    are smooth (infinitely differentiable) functions of $c$ for $c > 0$.

### Monotonicity of $\widetilde{\rho}(c)$

**Proposition 1 (Monotonicity; cf. Proposition A.1, Lee 2025).** The
average-information reliability $\widetilde{\rho}(c)$ is strictly
monotonically increasing in $c$ for $c > 0$.

**Proof sketch.** We show $\partial\widetilde{\rho}/\partial c > 0$.
From Eq. (5), by the chain rule and the fact that $h(x) = x/(x + 1)$ is
strictly increasing for $x > 0$, it suffices to show that
$\bar{\mathcal{J}}(c) = {\mathbb{E}}_{G}\left\lbrack \mathcal{J}(\theta;c) \right\rbrack$
is strictly increasing in $c$.

Differentiating the test information:

$$\frac{\partial\mathcal{J}_{i}(\theta;c)}{\partial c} = 2c\,\left( \lambda_{i}^{(0)} \right)^{2}\, p_{i}\left( 1 - p_{i} \right) + c^{2}\,\left( \lambda_{i}^{(0)} \right)^{2}\, p_{i}\left( 1 - p_{i} \right)\left( 1 - 2p_{i} \right) \cdot \lambda_{i}^{(0)}\left( \theta - \beta_{i} \right)$$

After taking the expectation over $G$ and summing over items, the
dominant first-order term
$2c\,\left( \lambda_{i}^{(0)} \right)^{2}\,{\mathbb{E}}\left\lbrack p_{i}\left( 1 - p_{i} \right) \right\rbrack$
is strictly positive for each item $i$ (since
$p_{i}\left( 1 - p_{i} \right) > 0$ with positive probability), ensuring
$\partial\bar{\mathcal{J}}(c)/\partial c > 0$.

A more rigorous argument uses the representation
$\bar{\mathcal{J}}(c) = c^{2}\sum_{i}\left( \lambda_{i}^{(0)} \right)^{2}{\mathbb{E}}\left\lbrack p_{i}(\theta;c)\left( 1 - p_{i}(\theta;c) \right) \right\rbrack$
and the fact that
${\mathbb{E}}\left\lbrack p_{i}\left( 1 - p_{i} \right) \right\rbrack$
is bounded away from zero under regularity conditions. $▫$

### Existence and Uniqueness

**Corollary 2 (Existence and Uniqueness; cf. Corollary 2.1, Lee 2025).**
Under mild regularity conditions, for any target
$\rho^{*} \in \left( \rho_{\min},\rho_{\max} \right)$ where

$$\rho_{\min} = \lim\limits_{c\rightarrow 0^{+}}\widetilde{\rho}(c) = 0\quad\text{and}\quad\rho_{\max} = \lim\limits_{c\rightarrow\infty}\widetilde{\rho}(c) = 1,$$

there exists a unique $c^{*} > 0$ satisfying
$\widetilde{\rho}\left( c^{*} \right) = \rho^{*}$.

**Proof.** Since $\widetilde{\rho}$ is continuous (composition of
continuous functions), strictly increasing (Proposition 1), and maps
$(0,\infty)$ onto $(0,1)$, existence and uniqueness follow from the
intermediate value theorem and strict monotonicity. $▫$

### Non-Monotonicity of $\bar{w}(c)$

**Proposition 2 (Non-Monotonicity of MSEM; cf. Proposition A.5, Lee
2025).** The MSEM-based reliability $\bar{w}(c)$ is *not* necessarily
monotone in $c$ over the entire positive real line. In particular, for
sufficiently large $c$, $\bar{w}(c)$ can decrease.

**Intuition.** As $\left. c\rightarrow\infty \right.$, for most $\theta$
values far from any $\beta_{i}$, the item response probabilities
$p_{i}(\theta)$ approach 0 or 1, causing item information
$\left. \mathcal{J}_{i}(\theta) = \lambda_{i}^{2}p_{i}\left( 1 - p_{i} \right)\rightarrow 0 \right.$.
Only $\theta$ values very close to some $\beta_{i}$ retain high
information. The harmonic mean (used in $\bar{w}$) is dominated by the
low-information regions, causing
${\mathbb{E}}\lbrack 1/\mathcal{J}\rbrack$ to grow and $\bar{w}$ to
decrease.

The arithmetic mean (used in $\widetilde{\rho}$) is not as sensitive to
these low- information tails, which is why $\widetilde{\rho}$ remains
monotone.

This is why the `c_bounds` parameter matters in practice, especially
when using the `"msem"` metric.

``` r
# Demonstrate non-monotonicity with a heavy-tailed distribution.
# Non-monotonicity is most prominent with heavy tails because extreme
# theta values far from any beta_i have near-zero information, and
# the harmonic mean (used in rho_bar) is dominated by these regions.
theta_nm <- sim_latentG(5000, shape = "heavy_tail", shape_params = list(df = 3))$theta

# Use 10 items with concentrated difficulties for a more visible effect
items_nm <- sim_item_params(10, model = "rasch", source = "parametric")
beta_nm  <- items_nm$data$beta
lambda_nm <- rep(1, 10)

c_wide <- seq(0.1, 15, length.out = 100)
rho_t <- rho_b <- numeric(length(c_wide))

for (j in seq_along(c_wide)) {
  both <- compute_rho_both(c_wide[j], theta_nm, beta_nm, lambda_nm)
  rho_t[j] <- both$rho_tilde
  rho_b[j] <- both$rho_bar
}

oldpar <- par(mar = c(4.5, 4.5, 3, 1))
on.exit(par(oldpar))

plot(c_wide, rho_t, type = "l", col = "steelblue", lwd = 2,
     xlab = "Scaling factor c", ylab = "Reliability",
     main = "Monotonicity Comparison (Heavy-Tailed Latent, df = 3)",
     ylim = c(0, 1))
lines(c_wide, rho_b, col = "coral", lwd = 2, lty = 2)

# Mark the maximum of rho_bar
idx_max <- which.max(rho_b)
points(c_wide[idx_max], rho_b[idx_max], pch = 19, col = "coral", cex = 1.5)

legend("right",
       legend = c(expression(tilde(rho) ~ "(always monotone)"),
                  expression(bar(w) ~ "(non-monotone at extreme c)")),
       col = c("steelblue", "coral"), lty = c(1, 2), lwd = 2, cex = 0.85)
```

![Figure 2: Non-monotonicity of the MSEM metric at extreme scaling
factors (heavy-tailed latent distribution, df =
3).](theory-reliability_files/figure-html/non-monotone-msem-1.png)

Figure 2: Non-monotonicity of the MSEM metric at extreme scaling factors
(heavy-tailed latent distribution, df = 3).

``` r

cat(sprintf("rho_bar maximum at c = %.2f with value %.4f\n",
            c_wide[idx_max], rho_b[idx_max]))
#> rho_bar maximum at c = 0.25 with value 0.1285
cat(sprintf("rho_bar at c = %.1f: %.4f (decreasing)\n",
            max(c_wide), rho_b[length(c_wide)]))
#> rho_bar at c = 15.0: 0.0000 (decreasing)
```

## Achievable Reliability Bounds

Not every target reliability is achievable for a given test
configuration. The feasible range depends on the number of items, the
item parameter distribution, and the latent trait distribution.

### Upper Bound from the Logistic Kernel

**Proposition 3 (Information Upper Bound).** Under the 2PL model, the
maximum item information is:

$$\sup\limits_{\theta}\,\mathcal{J}_{i}(\theta) = \frac{\lambda_{i}^{2}}{4}$$

achieved when $\theta = \beta_{i}$ (i.e., $p_{i} = 1/2$). Therefore, the
maximum test information is:

$$\sup\limits_{\theta}\,\mathcal{J}(\theta;c) = \frac{c^{2}}{4}\sum\limits_{i = 1}^{I}\left( \lambda_{i}^{(0)} \right)^{2}$$

For the Rasch model with $\lambda_{i}^{(0)} = 1$, this simplifies to
$c^{2}I/4$.

**Proof.** The function $p(1 - p)$ on $\lbrack 0,1\rbrack$ has a unique
maximum of $1/4$ at $p = 1/2$. Under the 2PL model,
$p_{i}(\theta) = 1/2$ when $\theta = \beta_{i}$. $▫$

### Feasibility Conditions

The average test information is bounded above by the maximum:

$$\bar{\mathcal{J}}(c) = {\mathbb{E}}_{G}\!\left\lbrack \sum\limits_{i = 1}^{I}\lambda_{i}(c)^{2}\, p_{i}\left( 1 - p_{i} \right) \right\rbrack \leq \frac{c^{2}}{4}\sum\limits_{i = 1}^{I}\left( \lambda_{i}^{(0)} \right)^{2}$$

Therefore the maximum achievable average-information reliability for a
given $c$ is:

$${\widetilde{\rho}}_{\max}(c) = \frac{\sigma_{\theta}^{2} \cdot \left( c^{2}/4 \right)\sum\limits_{i}\left( \lambda_{i}^{(0)} \right)^{2}}{\sigma_{\theta}^{2} \cdot \left( c^{2}/4 \right)\sum\limits_{i}\left( \lambda_{i}^{(0)} \right)^{2} + 1}$$

As $\left. c\rightarrow\infty \right.$, this approaches 1, but in
practice we work within
$c \in \left\lbrack c_{\min},c_{\max} \right\rbrack$ and the achievable
range is finite.

### Connection to `check_feasibility()`

The
[`check_feasibility()`](https://joonho112.github.io/IRTsimrel/reference/check_feasibility.md)
function computes the achievable reliability range via Monte Carlo
evaluation at the boundary scaling factors.

``` r
# Check feasibility for a 25-item Rasch test
feas <- check_feasibility(
  n_items = 25, model = "rasch",
  latent_shape = "normal",
  item_source = "parametric",
  c_bounds = c(0.1, 10),
  M = 5000L, seed = 42
)
#> 
#> =======================================================
#>   Feasibility Check: Achievable Reliability Range
#> =======================================================
#> 
#>   Number of items  : 25
#>   Model            : RASCH
#>   Latent shape     : normal
#>   Latent variance  : 1.0099
#>   c range          : [0.10, 10.00]
#>   Monte Carlo M    : 5000
#> 
#> Achievable Reliability Ranges:
#>   rho_tilde (info) : [0.0591, 0.9872]
#>   rho_bar   (msem) : [0.0002, 0.9146]
#> 
#> Note: rho_tilde >= rho_bar always (Jensen's inequality).
#>   Use rho_tilde range for EQC targets.
#>   Use rho_bar range for SAC targets.
```

``` r
# Test which targets are feasible
targets <- c(0.50, 0.70, 0.80, 0.90, 0.95)
for (rho in targets) {
  feasible_info <- rho >= feas$rho_range_info[1] & rho <= feas$rho_range_info[2]
  feasible_msem <- rho >= feas$rho_range_msem[1] & rho <= feas$rho_range_msem[2]
  cat(sprintf("  rho* = %.2f: info %s, msem %s\n",
              rho,
              ifelse(feasible_info, "YES", "NO "),
              ifelse(feasible_msem, "YES", "NO ")))
}
#>   rho* = 0.50: info YES, msem YES
#>   rho* = 0.70: info YES, msem YES
#>   rho* = 0.80: info YES, msem YES
#>   rho* = 0.90: info YES, msem YES
#>   rho* = 0.95: info YES, msem NO
```

### Reliability Curves

The
[`rho_curve()`](https://joonho112.github.io/IRTsimrel/reference/rho_curve.md)
function visualizes how reliability changes with the scaling factor.

``` r
curve_data <- rho_curve(
  c_values = seq(0.1, 5, length.out = 60),
  n_items = 25, model = "rasch",
  latent_shape = "normal",
  item_source = "parametric",
  metric = "both",
  M = 5000L, seed = 42,
  plot = TRUE
)
```

![Figure 3: Reliability curve for a 25-item Rasch test under a normal
latent
distribution.](theory-reliability_files/figure-html/rho-curve-1.png)

Figure 3: Reliability curve for a 25-item Rasch test under a normal
latent distribution.

## EQC Theory: Deterministic Root-Finding

The Empirical Quadrature Calibration (EQC) algorithm solves the inverse
reliability problem using deterministic root-finding on an empirical
reliability function. Full algorithmic details are in
[`vignette("algorithm-eqc")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-eqc.md);
here we state the main theoretical results.

### Algorithm Sketch

EQC proceeds in three steps:

1.  **Quadrature sampling**: Draw $\{\theta_{m}\}_{m = 1}^{M} \sim G$
    and item parameters
    $\{\left( \beta_{i},\lambda_{i}^{(0)} \right)\}_{i = 1}^{I} \sim H$
    once.

2.  **Empirical reliability function**: For any $c > 0$, define:
    $${\widehat{\rho}}_{M}(c) = \rho(c;\,\{\theta_{m}\},\,{\mathbf{β}},\,{\mathbf{λ}}(c))$$

3.  **Root-finding**: Solve
    $g_{M}(c) = {\widehat{\rho}}_{M}(c) - \rho^{*} = 0$ using Brent’s
    method on $\left\lbrack c_{\min},c_{\max} \right\rbrack$.

### Consistency

**Theorem 2 (Consistency of EQC; cf. Theorem A.1, Lee 2025).** Let
$c^{*}$ be the true population solution and ${\widehat{c}}_{M}^{*}$ be
the EQC solution with quadrature size $M$. Then:

$$\left. {\widehat{c}}_{M}^{*}\overset{\text{a.s.}}{\rightarrow}c^{*}\quad{\text{as}\mspace{6mu}}M\rightarrow\infty \right.$$

**Proof sketch.** The empirical reliability ${\widehat{\rho}}_{M}(c)$
converges uniformly to $\rho(c)$ over the compact interval
$\left\lbrack c_{\min},c_{\max} \right\rbrack$ by the Uniform Law of
Large Numbers (ULLN). Since $\rho(c)$ is continuous and strictly
monotone, and the root of the limit function is unique (Corollary 2),
uniform convergence of the objective implies convergence of the root (by
the continuous mapping theorem for M-estimators). $▫$

### Asymptotic Normality

**Theorem 3 (Asymptotic Normality of EQC; cf. Theorem A.2, Lee 2025).**
Under regularity conditions, the EQC estimator satisfies:

$$\sqrt{M}\,({\widehat{c}}_{M}^{*} - c^{*})\overset{d}{\rightarrow}N\!\left( 0,\,\frac{\sigma_{g}^{2}}{g\prime\left( c^{*} \right)^{2}} \right)$$

where $g(c) = \rho(c) - \rho^{*}$, $g\prime\left( c^{*} \right)$ is the
derivative of the reliability function at the solution, and
$\sigma_{g}^{2}$ is the asymptotic variance of the empirical gradient.

**Proof sketch.** Apply the delta method to the implicit function
${\widehat{c}}_{M}^{*} = g_{M}^{- 1}\left( \rho^{*} \right)$. The
Central Limit Theorem gives
$\left. \sqrt{M}\left( {\widehat{\rho}}_{M}(c) - \rho(c) \right)\rightarrow{}_{d}N\left( 0,\sigma_{g}^{2} \right) \right.$,
and the implicit function theorem (applicable since
$g\prime\left( c^{*} \right) \neq 0$ by strict monotonicity) yields the
result. $▫$

### Monte Carlo Error Bounds

The practical implication of Theorem 3 is:

$$\left| {\widehat{c}}_{M}^{*} - c^{*} \right| = O_{p}\!\left( \frac{1}{\sqrt{M}} \right)$$

The reliability estimate itself has error:

$$\left| {\widehat{\rho}}_{M}(c) - \rho(c) \right| = O_{p}\!\left( \frac{1}{\sqrt{M}} \right)$$

| Quadrature size $M$ | Approx. MC std. error |
|--------------------:|:---------------------:|
|               1,000 |         ~0.03         |
|               5,000 |        ~0.014         |
|              10,000 |         ~0.01         |
|              50,000 |        ~0.004         |
|             100,000 |        ~0.003         |

## SAC Theory: Stochastic Approximation

The Stochastic Approximation Calibration (SAC) algorithm provides an
independent solution path using the Robbins-Monro framework. Full
algorithmic details are in
[`vignette("algorithm-sac")`](https://joonho112.github.io/IRTsimrel/articles/algorithm-sac.md).

### The Stochastic Root-Finding Problem

SAC targets the equation:

$$g(c) = {\mathbb{E}}\!\left\lbrack \widehat{\rho}(c) \right\rbrack - \rho^{*} = 0$$

At each iteration $n$, only a noisy observation
${\widehat{\rho}}_{n}\left( c_{n} \right) = \rho^{*} + g\left( c_{n} \right) + \varepsilon_{n}$
is available, where $\varepsilon_{n}$ is zero-mean Monte Carlo noise.

### Robbins-Monro Update

The iterative update rule is:

$$c_{n + 1} = \Pi_{\lbrack c_{L},c_{U}\rbrack}\!\left\lbrack c_{n} - a_{n}\,({\widehat{\rho}}_{n}\left( c_{n} \right) - \rho^{*}) \right\rbrack$$

where $\Pi_{\lbrack c_{L},c_{U}\rbrack}$ denotes projection onto the
constraint set and the step size sequence satisfies the Robbins-Monro
conditions:

$$a_{n} = \frac{a}{(n + A)^{\gamma}},\quad\sum\limits_{n = 1}^{\infty}a_{n} = \infty,\quad\sum\limits_{n = 1}^{\infty}a_{n}^{2} < \infty$$

These conditions are satisfied when $\gamma \in (1/2,1\rbrack$.

### Polyak-Ruppert Averaging

Rather than using the final iterate $c_{N}$, SAC computes the
Polyak-Ruppert average over post-burn-in iterates:

$${\bar{c}}_{N} = \frac{1}{N - B}\sum\limits_{n = B + 1}^{N}c_{n}$$

where $B$ is the burn-in period. This averaging achieves the optimal
convergence rate.

### Convergence

**Theorem 4 (SAC Convergence; cf. Theorem A.3, Lee 2025).** Under the
Robbins-Monro conditions (Eq. 19) and standard regularity assumptions
(boundedness of $g$, Lipschitz continuity, noise variance bounded), the
SAC iterates converge almost surely:

$$\left. c_{n}\overset{\text{a.s.}}{\rightarrow}c^{*}\quad{\text{as}\mspace{6mu}}n\rightarrow\infty \right.$$

Furthermore, the Polyak-Ruppert average satisfies:

$$\sqrt{N - B}\,\left( {\bar{c}}_{N} - c^{*} \right)\overset{d}{\rightarrow}N\!\left( 0,\,\frac{\sigma_{\varepsilon}^{2}}{g\prime\left( c^{*} \right)^{2}} \right)$$

achieving the optimal $O\left( n^{- 1/2} \right)$ rate, compared to
$O\left( n^{- \gamma} \right)$ for the final iterate.

**Proof sketch.** Almost sure convergence follows from the classical
Robbins-Monro theorem (Robbins and Monro, 1951), extended to the
projected case by Kushner and Yin (2003). The key conditions are: (i)
$g\left( c^{*} \right) = 0$, (ii) $g\prime\left( c^{*} \right) > 0$
(strict monotonicity from Proposition 1), (iii) the step size conditions
in Eq. (19), and (iv) bounded conditional variance
$\text{Var}\left( \varepsilon_{n} \mid c_{n} \right) \leq \sigma_{\varepsilon}^{2}$.
Asymptotic normality of the Polyak-Ruppert average follows from Polyak
and Juditsky (1992). $▫$

## Numerical Demonstrations

### Verifying Theoretical Predictions

We use the low-level reliability functions to verify key theoretical
properties.

``` r
# Common setup for demonstrations
M_demo <- 5000L
theta_demo <- sim_latentG(M_demo, shape = "normal")$theta
items_demo <- sim_item_params(30, model = "rasch", source = "parametric")
beta_demo  <- items_demo$data$beta
lambda_demo <- rep(1, 30)
theta_var_demo <- var(theta_demo)
```

#### Verification 1: Reliability increases with c (info metric)

``` r
c_test <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
rho_info <- sapply(c_test, function(cc) {
  compute_rho_tilde(cc, theta_demo, beta_demo, lambda_demo,
                    theta_var = theta_var_demo)
})

cat("Monotonicity check (info metric):\n")
#> Monotonicity check (info metric):
for (j in seq_along(c_test)) {
  cat(sprintf("  c = %.1f: rho_tilde = %.4f", c_test[j], rho_info[j]))
  if (j > 1) {
    cat(sprintf("  (diff = %+.4f)", rho_info[j] - rho_info[j-1]))
  }
  cat("\n")
}
#>   c = 0.5: rho_tilde = 0.6286
#>   c = 1.0: rho_tilde = 0.8456  (diff = +0.2171)
#>   c = 1.5: rho_tilde = 0.9082  (diff = +0.0626)
#>   c = 2.0: rho_tilde = 0.9351  (diff = +0.0269)
#>   c = 2.5: rho_tilde = 0.9498  (diff = +0.0147)
#>   c = 3.0: rho_tilde = 0.9590  (diff = +0.0092)
cat(sprintf("All differences positive: %s\n", all(diff(rho_info) > 0)))
#> All differences positive: TRUE
```

#### Verification 2: Jensen’s gap depends on information heterogeneity

Non-standard latent distributions produce larger Jensen’s gaps because
information varies more across the ability continuum.

``` r
shapes <- c("normal", "bimodal", "heavy_tail", "skew_pos")
shape_pars <- list(
  normal     = list(),
  bimodal    = list(delta = 0.9),
  heavy_tail = list(df = 5),
  skew_pos   = list(k = 4)
)

c_val <- 1.5
cat(sprintf("Jensen's gap at c = %.1f for different latent shapes:\n", c_val))
#> Jensen's gap at c = 1.5 for different latent shapes:

for (sh in shapes) {
  theta_sh <- sim_latentG(M_demo, shape = sh, shape_params = shape_pars[[sh]])$theta
  both_sh <- compute_rho_both(c_val, theta_sh, beta_demo, lambda_demo)
  gap <- both_sh$rho_tilde - both_sh$rho_bar
  cat(sprintf("  %-12s: rho_tilde = %.4f, rho_bar = %.4f, gap = %.4f\n",
              sh, both_sh$rho_tilde, both_sh$rho_bar, gap))
}
#>   normal      : rho_tilde = 0.9113, rho_bar = 0.9011, gap = 0.0102
#>   bimodal     : rho_tilde = 0.9070, rho_bar = 0.9040, gap = 0.0031
#>   heavy_tail  : rho_tilde = 0.9117, rho_bar = 0.8293, gap = 0.0824
#>   skew_pos    : rho_tilde = 0.9122, rho_bar = 0.8862, gap = 0.0260
```

#### Verification 3: APC initialization quality

The Analytic Pre-Calibration (APC) formula provides a closed-form
starting point. We compare it against the EQC solution.

``` r
targets <- c(0.60, 0.70, 0.80, 0.85, 0.90)
cat("APC initialization vs EQC solution (Rasch, 25 items, normal):\n")
#> APC initialization vs EQC solution (Rasch, 25 items, normal):
cat(sprintf("  %-10s %-10s %-10s %-10s\n", "Target", "c_APC", "c_EQC", "Ratio"))
#>   Target     c_APC      c_EQC      Ratio

for (rho_t in targets) {
  c_apc <- compute_apc_init(target_rho = rho_t, n_items = 25)
  eqc_t <- eqc_calibrate(
    target_rho = rho_t, n_items = 25, model = "rasch",
    item_source = "parametric", M = 5000L, seed = 42, verbose = FALSE
  )
  cat(sprintf("  %-10.2f %-10.4f %-10.4f %-10.2f\n",
              rho_t, c_apc, eqc_t$c_star, c_apc / eqc_t$c_star))
}
#>   0.60       0.8129     0.5118     1.59      
#>   0.70       1.0138     0.6548     1.55      
#>   0.80       1.3274     0.8995     1.48      
#>   0.85       1.5799     1.1199     1.41      
#>   0.90       1.9911     1.5328     1.30
```

#### Verification 4: Effect of test length on achievable reliability

``` r
n_items_grid <- c(10, 15, 20, 25, 30, 40, 50)
c_fixed <- 1.0

cat(sprintf("Reliability at c = %.1f for different test lengths:\n", c_fixed))
#> Reliability at c = 1.0 for different test lengths:
for (ni in n_items_grid) {
  items_ni <- sim_item_params(ni, model = "rasch", source = "parametric")
  beta_ni  <- items_ni$data$beta
  lambda_ni <- rep(1, ni)

  both_ni <- compute_rho_both(c_fixed, theta_demo, beta_ni, lambda_ni,
                               theta_var = theta_var_demo)
  cat(sprintf("  I = %2d: rho_tilde = %.4f, rho_bar = %.4f\n",
              ni, both_ni$rho_tilde, both_ni$rho_bar))
}
#>   I = 10: rho_tilde = 0.6514, rho_bar = 0.6406
#>   I = 15: rho_tilde = 0.7213, rho_bar = 0.7158
#>   I = 20: rho_tilde = 0.7797, rho_bar = 0.7738
#>   I = 25: rho_tilde = 0.8165, rho_bar = 0.8112
#>   I = 30: rho_tilde = 0.8510, rho_bar = 0.8442
#>   I = 40: rho_tilde = 0.8768, rho_bar = 0.8730
#>   I = 50: rho_tilde = 0.8995, rho_bar = 0.8961
```

### EQC vs SAC Agreement

``` r
eqc_demo <- eqc_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  M = 5000L, seed = 42, verbose = FALSE
)

sac_demo <- sac_calibrate(
  target_rho = 0.80, n_items = 25, model = "rasch",
  item_source = "parametric", reliability_metric = "info",
  c_init = eqc_demo, n_iter = 300L, M_per_iter = 1000L,
  seed = 42, verbose = FALSE
)

cat("EQC vs SAC comparison (target = 0.80, info metric):\n")
#> EQC vs SAC comparison (target = 0.80, info metric):
cat(sprintf("  EQC c*:        %.4f\n", eqc_demo$c_star))
#>   EQC c*:        0.8995
cat(sprintf("  SAC c*:        %.4f\n", sac_demo$c_star))
#>   SAC c*:        0.9122
cat(sprintf("  Abs diff:      %.4f\n", abs(eqc_demo$c_star - sac_demo$c_star)))
#>   Abs diff:      0.0127
cat(sprintf("  Pct diff:      %.2f%%\n",
            100 * abs(eqc_demo$c_star - sac_demo$c_star) / eqc_demo$c_star))
#>   Pct diff:      1.41%
```

### Verification 5: Reliability function under the 2PL model

The 2PL model introduces additional variability through heterogeneous
baseline discriminations.

``` r
items_2pl <- sim_item_params(25, model = "2pl", source = "parametric")
beta_2pl  <- items_2pl$data$beta
lambda_2pl <- if ("lambda_unscaled" %in% names(items_2pl$data)) {
  items_2pl$data$lambda_unscaled
} else {
  items_2pl$data$lambda
}

cat("2PL baseline discriminations:\n")
#> 2PL baseline discriminations:
cat(sprintf("  Range: [%.3f, %.3f]\n", min(lambda_2pl), max(lambda_2pl)))
#>   Range: [0.411, 2.144]
cat(sprintf("  Mean:  %.3f\n", mean(lambda_2pl)))
#>   Mean:  1.165
cat(sprintf("  SD:    %.3f\n", sd(lambda_2pl)))
#>   SD:    0.383

# Compare Rasch and 2PL reliability at same c
c_compare <- 1.0
both_rasch <- compute_rho_both(c_compare, theta_demo,
                                beta_demo[1:25], rep(1, 25),
                                theta_var = theta_var_demo)
both_2pl <- compute_rho_both(c_compare, theta_demo,
                              beta_2pl, lambda_2pl,
                              theta_var = theta_var_demo)

cat(sprintf("\nReliability comparison at c = %.1f:\n", c_compare))
#> 
#> Reliability comparison at c = 1.0:
cat(sprintf("  Rasch: rho_tilde = %.4f, rho_bar = %.4f, gap = %.4f\n",
            both_rasch$rho_tilde, both_rasch$rho_bar,
            both_rasch$rho_tilde - both_rasch$rho_bar))
#>   Rasch: rho_tilde = 0.8213, rho_bar = 0.8152, gap = 0.0062
cat(sprintf("  2PL:   rho_tilde = %.4f, rho_bar = %.4f, gap = %.4f\n",
            both_2pl$rho_tilde, both_2pl$rho_bar,
            both_2pl$rho_tilde - both_2pl$rho_bar))
#>   2PL:   rho_tilde = 0.8556, rho_bar = 0.8444, gap = 0.0112
```

Heterogeneous discriminations typically produce a larger Jensen’s gap
because the test information varies more across the ability continuum.

## The Inverse Design Problem

The theoretical framework above describes the *forward* problem: given a
scaling factor $c$, compute the reliability $\rho(c)$. The central
contribution of the IRTsimrel package is solving the *inverse* problem.

### Problem Statement

**Definition 7 (Inverse Reliability Problem).** Given a target marginal
reliability $\rho^{*} \in (0,1)$, find a global discrimination scaling
factor $c^{*} > 0$ such that:

$$\rho\left( c^{*} \right) = \rho^{*}$$

where $\rho$ denotes either $\widetilde{\rho}$ (average-information) or
$\bar{w}$ (MSEM-based), depending on the chosen metric.

### Why the Inverse Problem is Difficult

Several factors make this inverse problem non-trivial:

1.  **No closed-form solution**: The reliability function involves a
    logistic nonlinearity inside an expectation over $G$. Even in the
    simplest Rasch case with normal $G$, no analytic inversion is
    possible.

2.  **Distribution dependence**: The mapping
    $\left. c\mapsto\rho(c) \right.$ depends on the full distributions
    $G$ (latent traits) and $H$ (item parameters), not just their
    moments.

3.  **Metric ambiguity**: The two reliability definitions yield
    different $c^{*}$ values for the same target (Corollary 1).

4.  **Practical constraints**: The scaling factor must lie within bounds
    that produce realistic item parameters.

### The APC Approximation

The Analytic Pre-Calibration (APC) formula provides a rough initial
estimate under Gaussian Rasch assumptions. With $\theta \sim N(0,1)$ and
$\beta \sim N\left( 0,\sigma_{\beta}^{2} \right)$, the expected item
information involves the logistic-normal convolution constant:

$$\kappa\left( \sigma^{2} \right) = \int_{- \infty}^{\infty}\frac{e^{z}}{\left( 1 + e^{z} \right)^{2}}\,\phi\left( z;0,\sigma^{2} \right)\, dz \approx \frac{0.25}{\sqrt{1 + \sigma^{2}\pi^{2}/3}}$$

where $\sigma^{2} = 1 + \sigma_{\beta}^{2}$. The APC initial value is:

$$c_{\text{APC}} = \sqrt{\frac{\rho^{*}}{I \cdot \kappa \cdot \left( 1 - \rho^{*} \right)}}$$

This approximation is typically within a factor of 2 of the true $c^{*}$
and serves as a starting point for the stochastic approximation
algorithm (SAC).

## Complete Notation Table

The following table summarizes all mathematical notation used in the
package and its documentation.

|          Symbol           | Description                                                                                 |       Domain        |
|:-------------------------:|:--------------------------------------------------------------------------------------------|:-------------------:|
|         $\theta$          | Latent ability parameter                                                                    |    $\mathbb{R}$     |
|            $G$            | Latent trait distribution                                                                   |          –          |
|   $\sigma_{\theta}^{2}$   | Variance of $\theta$ under $G$                                                              |    $(0,\infty)$     |
|        $\beta_{i}$        | Difficulty parameter for item $i$                                                           |    $\mathbb{R}$     |
|       $\lambda_{i}$       | Discrimination parameter for item $i$                                                       |    $(0,\infty)$     |
|    $\lambda_{i}^{(0)}$    | Baseline (unscaled) discrimination                                                          |    $(0,\infty)$     |
|            $c$            | Global discrimination scaling factor                                                        |    $(0,\infty)$     |
|          $c^{*}$          | Calibrated scaling factor                                                                   |    $(0,\infty)$     |
|            $I$            | Number of items                                                                             |    $\mathbb{N}$     |
|            $M$            | Quadrature (Monte Carlo) sample size                                                        |    $\mathbb{N}$     |
|         $\psi(z)$         | Logistic function $1/\left( 1 + e^{- z} \right)$                                            |       $(0,1)$       |
|      $p_{i}(\theta)$      | $P\left( X_{i} = 1 \mid \theta \right)$; ICC                                                |       $(0,1)$       |
| $\mathcal{J}_{i}(\theta)$ | Item information function                                                                   |    $(0,\infty)$     |
|   $\mathcal{J}(\theta)$   | Test information function                                                                   |    $(0,\infty)$     |
|  $\bar{\mathcal{J}}(c)$   | Average test information ${\mathbb{E}}_{G}\left\lbrack \mathcal{J}(\theta;c) \right\rbrack$ |    $(0,\infty)$     |
|   $\widetilde{\rho}(c)$   | Average-information reliability                                                             |       $(0,1)$       |
|       $\bar{w}(c)$        | MSEM-based reliability                                                                      |       $(0,1)$       |
|        $\rho^{*}$         | Target reliability                                                                          |       $(0,1)$       |
| ${\widehat{\rho}}_{M}(c)$ | Empirical reliability (EQC)                                                                 |       $(0,1)$       |
| ${\widehat{\rho}}_{n}(c)$ | Noisy reliability estimate (SAC)                                                            |       $(0,1)$       |
|          $a_{n}$          | Step size at iteration $n$                                                                  |    $(0,\infty)$     |
|            $a$            | Step size constant (SAC)                                                                    |    $(0,\infty)$     |
|            $A$            | Stabilization constant (SAC)                                                                | $\lbrack 0,\infty)$ |
|         $\gamma$          | Step decay exponent (SAC)                                                                   |   $(1/2,1\rbrack$   |
|            $B$            | Burn-in period (SAC)                                                                        | ${\mathbb{N}}_{0}$  |
|            $N$            | Total number of SAC iterations                                                              |    $\mathbb{N}$     |
|      ${\bar{c}}_{N}$      | Polyak-Ruppert average                                                                      |    $(0,\infty)$     |
|         $\kappa$          | Logistic-normal convolution constant (APC)                                                  |      $(0,1/4)$      |
|            $H$            | Item parameter distribution                                                                 |          –          |
|       $\text{MSEM}$       | ${\mathbb{E}}\left\lbrack 1/\mathcal{J}(\theta;c) \right\rbrack$                            |    $(0,\infty)$     |

## Summary of Theoretical Results

| Result            | Statement                                                                     | Implication                                         |
|:------------------|:------------------------------------------------------------------------------|:----------------------------------------------------|
| **Theorem 1**     | $\widetilde{\rho}(c) \geq \bar{w}(c)$ (Jensen)                                | Info metric always yields higher reliability values |
| **Corollary 1**   | $c_{\bar{w}}^{*} \geq c_{\widetilde{\rho}}^{*}$                               | MSEM metric requires larger scaling factor          |
| **Proposition 1** | $\widetilde{\rho}(c)$ is strictly increasing                                  | Existence and uniqueness for info metric            |
| **Corollary 2**   | Unique $c^{*}$ exists for $\rho^{*} \in (0,1)$                                | EQC root-finding is well-defined                    |
| **Proposition 2** | $\bar{w}(c)$ is not always monotone                                           | `c_bounds` critical for MSEM metric                 |
| **Proposition 3** | $\max\mathcal{J}_{i} = \lambda_{i}^{2}/4$                                     | Upper bound on achievable reliability               |
| **Theorem 2**     | EQC consistent ($\left. {\widehat{c}}_{M}^{*}\rightarrow c^{*} \right.$ a.s.) | Arbitrarily precise via larger $M$                  |
| **Theorem 3**     | EQC asymptotically normal ($\sqrt{M}$ rate)                                   | Quantifiable uncertainty                            |
| **Theorem 4**     | SAC converges a.s., optimal rate with PR avg                                  | Independent validation pathway                      |

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

Brent, R. P. (1973). *Algorithms for Minimization Without Derivatives*.
Prentice-Hall.

Lord, F. M. (1980). *Applications of Item Response Theory to Practical
Testing Problems*. Erlbaum.

van der Linden, W. J. (Ed.). (2016). *Handbook of Item Response Theory*
(Vols. 1–3). CRC Press.

Birnbaum, A. (1968). Some latent trait models and their use in inferring
an examinee’s ability. In F. M. Lord & M. R. Novick (Eds.), *Statistical
Theories of Mental Test Scores* (pp. 395–479). Addison-Wesley.
