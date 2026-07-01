# Plot method for item_params objects

Plot method for item_params objects

## Usage

``` r
# S3 method for class 'item_params'
plot(x, type = c("scatter", "density", "both"), ...)
```

## Arguments

- x:

  An `item_params` object from
  [`sim_item_params()`](https://joonho112.github.io/IRTsimrel/reference/sim_item_params.md).

- type:

  Character. Type of plot: `"scatter"`, `"density"`, or `"both"`.

- ...:

  Additional arguments (currently ignored).

## Value

A `ggplot` object if ggplot2 is available, or `NULL` invisibly when
using base R graphics fallback.
