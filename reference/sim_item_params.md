# Simulate Item Parameters for IRT Studies

`sim_item_params()` generates item parameters (difficulty \\\beta\\ and
discrimination \\\lambda\\) for Item Response Theory (IRT) simulation
studies. It wraps the IRW `irw_simu_diff()` function for realistic
difficulty distributions and provides multiple methods for generating
correlated discriminations.

The function is designed with four key principles:

1.  **Realistic difficulties:** Integration with Item Response Warehouse
    (IRW) for empirically-grounded difficulty distributions.

2.  **Correlated parameters:** Support for the empirically observed
    negative correlation between difficulty and discrimination (Sweeney
    et al., 2022).

3.  **Marginal preservation:** Copula method preserves exact marginal
    distributions while achieving target correlation.

4.  **Reliability targeting:** Scale factor for subsequent calibration.

## Usage

``` r
sim_item_params(
  n_items,
  model = c("rasch", "2pl"),
  source = c("irw", "parametric", "hierarchical", "custom"),
  method = c("copula", "conditional", "independent"),
  n_forms = 1L,
  difficulty_params = list(),
  discrimination_params = list(),
  hierarchical_params = list(),
  custom_params = list(),
  scale = 1,
  center_difficulties = TRUE,
  seed = NULL
)
```

## Arguments

- n_items:

  Integer. Number of items to generate per form.

- model:

  Character. The data-generating model: "rasch" or "2pl".

- source:

  Character. Source for generating difficulties:

  `"irw"`

  :   Use IRW difficulty pool (realistic, empirical)

  `"parametric"`

  :   Generate from parametric distribution

  `"hierarchical"`

  :   Joint MVN for both parameters (Glas & van der Linden)

  `"custom"`

  :   User-supplied parameters or function

- method:

  Character. Method for generating discriminations (when model = "2pl"):

  `"copula"`

  :   Gaussian copula - preserves marginals exactly (RECOMMENDED)

  `"conditional"`

  :   Conditional normal regression on difficulty

  `"independent"`

  :   Independent generation (no correlation)

- n_forms:

  Integer. Number of test forms to generate. Default is 1. When \> 1,
  returns a data frame with form_id column.

- difficulty_params:

  List. Parameters for difficulty generation:

  For `source = "irw"`:

  :   `pool` - difficulty pool data frame

  For `source = "parametric"`:

  :   `mu`, `sigma`, `distribution`

- discrimination_params:

  List. Parameters for discrimination generation:

  `mu_log`

  :   Mean of log-discrimination (default: 0)

  `sigma_log`

  :   SD of log-discrimination (default: 0.3)

  `rho`

  :   Target correlation between \\\beta\\ and \\\log(\lambda)\\
      (default: -0.3)

- hierarchical_params:

  List. For source = "hierarchical":

  `mu`

  :   2-vector: means of \\(\log\lambda, \beta)\\

  `tau`

  :   2-vector: SDs

  `rho`

  :   Correlation

- custom_params:

  List. For source = "custom":

  `beta`

  :   Vector or function returning difficulties

  `lambda`

  :   Vector or function returning discriminations

- scale:

  Numeric. Global discrimination scaling factor for reliability
  targeting. Final discriminations are \\\lambda_i^\* = c \cdot
  \lambda_i\\. Default is 1.

- center_difficulties:

  Logical. If TRUE, center difficulties to sum to zero for
  identification. Default is TRUE.

- seed:

  Integer. Random seed for reproducibility.

## Value

An object of class `"item_params"` containing:

- `data`:

  Data frame with columns: form_id, item_id, beta, lambda,
  lambda_unscaled

- `model`:

  Model type used

- `source`:

  Source used for generation

- `method`:

  Method used for discrimination generation

- `n_items`:

  Number of items per form

- `n_forms`:

  Number of forms generated

- `scale`:

  Scale factor applied

- `centered`:

  Whether difficulties were centered

- `params`:

  Parameters used for generation

- `achieved`:

  Achieved statistics (correlations, moments)

## Details

### Why Copula Method is Recommended

When difficulties come from the IRW pool (which has realistic, often
non-normal marginal distributions), the conditional normal method can
distort the achieved correlation because it assumes linearity. The
Gaussian copula method:

1.  Transforms difficulties to uniform scale via empirical CDF

2.  Generates correlated uniforms through Gaussian copula

3.  Transforms back to desired marginals (log-normal for discrimination)

This guarantees:

- Exact preservation of difficulty marginal (whatever IRW provides)

- Exact log-normal marginal for discriminations

- Spearman correlation \\\approx \rho\\ (rank-based, robust to
  non-normality

### Connection to Reliability-Targeted Framework

The `scale` parameter implements "separation of structure and scale":

- **Structure**: Realistic item characteristics from IRW + correlation

- **Scale**: Global factor \\c\\ calibrated for target reliability

## References

Glas, C. A. W., & van der Linden, W. J. (2003). Computerized adaptive
testing with item cloning. *Applied Psychological Measurement, 27*(4),
247-261.

Sweeney, S. M., et al. (2022). An investigation of the nature and
consequence of the relationship between IRT difficulty and
discrimination. *EM:IP, 41*(4), 50-67.

Zhang, L., et al. (2025). Realistic simulation of item difficulties.
*PsyArXiv*.

## Examples

``` r
# Example 1: Rasch with IRW difficulties
items1 <- sim_item_params(n_items = 25, model = "rasch", source = "irw")

# Example 2: 2PL with copula method (recommended)
items2 <- sim_item_params(
  n_items = 30, model = "2pl", source = "irw",
  method = "copula",
  discrimination_params = list(rho = -0.3)
)

# Example 3: Multiple forms
items3 <- sim_item_params(
  n_items = 20, model = "2pl", n_forms = 5,
  source = "irw", method = "copula"
)
#> Warning: collapsing to unique 'x' values

# Example 4: Hierarchical 2PL
items4 <- sim_item_params(
  n_items = 25, model = "2pl", source = "hierarchical",
  hierarchical_params = list(mu = c(0, 0), tau = c(0.25, 1), rho = -0.3)
)
```
