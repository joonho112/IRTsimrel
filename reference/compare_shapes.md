# Compare Multiple Distribution Shapes

Generates and compares multiple latent ability distributions
side-by-side.

## Usage

``` r
compare_shapes(
  n = 2000,
  shapes = c("normal", "bimodal", "trimodal", "skew_pos", "skew_neg", "heavy_tail",
    "uniform"),
  sigma = 1,
  seed = NULL
)
```

## Arguments

- n:

  Integer. Sample size for each distribution.

- shapes:

  Character vector. Shapes to compare.

- sigma:

  Numeric. Common scale parameter.

- seed:

  Integer. Random seed.

## Value

A ggplot2 object with faceted density plots.

## Examples

``` r
if (FALSE) { # \dontrun{
compare_shapes(n = 3000)
} # }
```
