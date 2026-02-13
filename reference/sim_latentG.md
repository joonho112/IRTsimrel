# Simulate Latent Ability Distribution for IRT Studies (G-family)

`sim_latentG()` generates latent abilities (person parameters
\\\theta\\) for Item Response Theory (IRT) simulation studies. It
implements the population model \\\theta_p \sim G\\ where \\G\\ is a
flexible distribution family.

The function is designed with two key principles:

1.  **Pre-standardization:** Each distribution shape is mathematically
    constructed to have mean 0 and variance 1, ensuring that changing
    the shape does not inadvertently change the scale.

2.  **Separation of Structure and Scale:** The `sigma` parameter
    directly controls the standard deviation of the latent trait,
    independent of the distributional shape.

The generated abilities follow: \$\$\theta_p = \mu + X_p^\top \beta +
\sigma \cdot z_p\$\$ where \\z_p \sim G_0\\ with \\E\[z\]=0\\ and
\\Var\[z\]=1\\.

## Usage

``` r
sim_latentG(
  n,
  shape = c("normal", "bimodal", "trimodal", "multimodal", "skew_pos", "skew_neg",
    "heavy_tail", "light_tail", "uniform", "floor", "ceiling", "custom"),
  sigma = 1,
  mu = 0,
  xcov = NULL,
  beta = NULL,
  shape_params = list(),
  mixture_spec = NULL,
  standardize_custom = TRUE,
  seed = NULL,
  return_z = TRUE
)
```

## Arguments

- n:

  Integer. Number of persons (latent abilities) to generate.

- shape:

  Character. The distributional shape of the standardized component. One
  of:

  `"normal"`

  :   Standard normal \\N(0,1)\\

  `"bimodal"`

  :   Symmetric two-component Gaussian mixture with analytically
      standardized parameters

  `"trimodal"`

  :   Symmetric three-component Gaussian mixture

  `"multimodal"`

  :   Four-component Gaussian mixture

  `"skew_pos"`

  :   Right-skewed distribution via standardized Gamma

  `"skew_neg"`

  :   Left-skewed distribution (negated Gamma)

  `"heavy_tail"`

  :   Heavy-tailed distribution via standardized Student-t

  `"light_tail"`

  :   Light-tailed (platykurtic) mixture distribution

  `"uniform"`

  :   Uniform distribution on \\\[-\sqrt{3}, \sqrt{3}\]\\

  `"floor"`

  :   Distribution with floor effect (left-truncated feel)

  `"ceiling"`

  :   Distribution with ceiling effect (right-truncated feel)

  `"custom"`

  :   User-specified mixture distribution

- sigma:

  Numeric. Scale (standard deviation) of the residual latent trait.
  Since the standardized component has variance 1, `sigma` directly
  equals the marginal SD of the residual term. Default is 1.

- mu:

  Numeric. Grand mean of the latent ability distribution. In Rasch
  models this is often fixed to 0 for identification. Default is 0.

- xcov:

  Matrix or data.frame. Optional covariate matrix with `n` rows. If
  supplied, person-specific covariate effects are added as \\\eta =
  X\beta\\.

- beta:

  Numeric vector. Regression coefficients for `xcov`. Must have length
  equal to `ncol(xcov)`. Ignored if `xcov` is NULL.

- shape_params:

  List. Additional parameters controlling the shape. See Details for
  shape-specific parameters.

- mixture_spec:

  List. For `shape = "custom"`, specifies the mixture:

  `weights`

  :   Numeric vector of mixing proportions (must sum to 1)

  `means`

  :   Numeric vector of component means

  `sds`

  :   Numeric vector of component standard deviations

  The custom mixture is automatically standardized to have mean 0 and
  variance 1.

- standardize_custom:

  Logical. If TRUE (default), custom mixtures are post-standardized to
  ensure mean 0 and variance 1. If FALSE, the raw mixture is used (user
  must ensure proper standardization).

- seed:

  Integer. Random seed for reproducibility. If NULL (default), the
  current RNG state is used.

- return_z:

  Logical. If TRUE, include the standardized draws `z` in the output.
  Default is TRUE.

## Value

An object of class `"latent_G"` (a list) containing:

- `theta`:

  Numeric vector of length `n`, the simulated latent abilities

- `z`:

  Standardized draws (if `return_z = TRUE`)

- `eta_cov`:

  Covariate linear predictor (0 if no covariates)

- `mu`:

  Grand mean used

- `sigma`:

  Scale parameter used

- `shape`:

  Shape label

- `shape_params`:

  Shape parameters used

- `n`:

  Sample size

- `sample_moments`:

  List with sample mean, sd, skewness, kurtosis

## Details

### Pre-standardization Mathematics

Each built-in shape is constructed to have exactly mean 0 and variance
1:

**Bimodal:** Two-component mixture with modes at \\\pm\delta\\: \$\$z =
s \cdot \delta + \epsilon, \quad s \sim \text{Rademacher}, \quad
\epsilon \sim N(0, 1-\delta^2)\$\$ where the component variance
\\1-\delta^2\\ ensures \\Var\[z\] = \delta^2 + (1-\delta^2) = 1\\.

**Trimodal:** Three-component mixture with weights \\(w_L, w_0, w_R)\\
and means \\(-m, 0, m)\\. Component variance is \\\sigma_c^2 = 1 -
(1-w_0)m^2\\ to ensure unit total variance.

**Skewed:** Standardized Gamma distribution: \$\$z =
\frac{\Gamma(k, 1) - k}{\sqrt{k}}\$\$ which has \\E\[z\]=0\\ and
\\Var\[z\]=1\\ for any \\k \> 0\\.

**Heavy-tailed:** Standardized Student-t: \$\$z =
\frac{t\_\nu}{\sqrt{\nu/(\nu-2)}}\$\$ which has \\Var\[z\]=1\\ for \\\nu
\> 2\\.

### Shape-Specific Parameters

- `delta`:

  For "bimodal": mode separation, must satisfy \\0 \< \delta \< 1\\.
  Default: 0.8

- `w0`:

  For "trimodal": weight of central component, must satisfy \\0 \< w_0
  \< 1\\. Default: 1/3

- `m`:

  For "trimodal": magnitude of side component means. Must satisfy
  \\(1-w_0)m^2 \< 1\\. Default: 1.2

- `k`:

  For "skew_pos"/"skew_neg": Gamma shape parameter, controls skewness
  magnitude. Default: 4

- `df`:

  For "heavy_tail": degrees of freedom, must be \> 2. Default: 5

### Connection to IRT Framework

In the Rasch/2PL model, the latent distribution \\G\\ affects:

- Marginal reliability: \\\bar{w} = \sigma\_\theta^2 /
  (\sigma\_\theta^2 + \text{MSEM})\\

- Expected test information: \\\bar{\mathcal{J}} =
  E_G\[\mathcal{J}(\theta)\]\\

- Identifiability (see Appendix F of the manuscript)

This function serves as the generator for \\G\\ in reliability-targeted
simulation studies, allowing researchers to examine how distributional
shape affects model performance while holding scale constant.

## References

Baker, F. B., & Kim, S.-H. (2004). *Item Response Theory: Parameter
Estimation Techniques* (2nd ed.). Marcel Dekker.

Paganin, S., et al. (2022). Computational strategies and estimation
performance with Bayesian semiparametric item response theory models.
*Journal of Educational and Behavioral Statistics, 48*(2), 147-188.

## See also

[`summary.latent_G`](https://joonho112.github.io/IRTsimrel/reference/summary.latent_G.md)
for summary statistics,
[`plot.latent_G`](https://joonho112.github.io/IRTsimrel/reference/plot.latent_G.md)
for visualization,
[`compare_shapes`](https://joonho112.github.io/IRTsimrel/reference/compare_shapes.md)
for comparing multiple shapes.

## Examples

``` r
# Basic usage: standard normal abilities
sim1 <- sim_latentG(n = 1000, shape = "normal")
mean(sim1$theta)  # approximately 0
#> [1] -0.01441477
sd(sim1$theta)    # approximately 1
#> [1] 1.040669

# Bimodal distribution for heterogeneous population
sim2 <- sim_latentG(n = 1000, shape = "bimodal",
                    shape_params = list(delta = 0.9))

# Skewed distribution with larger scale
sim3 <- sim_latentG(n = 1000, shape = "skew_pos", sigma = 1.5)

# With covariate effects (e.g., group differences)
group <- rbinom(1000, 1, 0.5)
sim4 <- sim_latentG(n = 1000, shape = "normal",
                    xcov = data.frame(group = group),
                    beta = 0.5)

# Custom mixture distribution
sim5 <- sim_latentG(n = 1000, shape = "custom",
                    mixture_spec = list(
                      weights = c(0.3, 0.5, 0.2),
                      means = c(-1.5, 0, 2),
                      sds = c(0.5, 0.7, 0.5)
                    ))
```
