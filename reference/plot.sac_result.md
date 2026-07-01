# Plot SAC Convergence Trajectory

Visualizes the Robbins-Monro iteration trajectory.

## Usage

``` r
# S3 method for class 'sac_result'
plot(x, type = c("both", "trajectory", "c", "rho"), ...)
```

## Arguments

- x:

  An `sac_result` object from
  [`sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md).

- type:

  Character. Plot type: `"trajectory"`, `"rho"`, `"c"`, or `"both"`.

- ...:

  Additional arguments (currently unused).

## Value

A `ggplot` object if ggplot2 is available for `type = "c"` or
`type = "rho"`. For `type = "both"`, returns a patchwork object when
patchwork is available; otherwise prints both ggplots and returns them
invisibly as a list. Returns `NULL` invisibly when using base R graphics
fallback.
