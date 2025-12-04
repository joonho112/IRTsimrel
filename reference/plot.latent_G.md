# Plot Method for latent_G Objects

Plot Method for latent_G Objects

## Usage

``` r
# S3 method for class 'latent_G'
plot(
  x,
  type = c("both", "histogram", "density"),
  show_normal = TRUE,
  bins = 50,
  ...
)
```

## Arguments

- x:

  A `latent_G` object from
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- type:

  Character. Type of plot: "histogram", "density", or "both".

- show_normal:

  Logical. Overlay a normal reference density?

- bins:

  Integer. Number of histogram bins.

- ...:

  Additional arguments passed to plotting functions.

## Value

A ggplot2 object (if available) or base R plot.
