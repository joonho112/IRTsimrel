# Plot SPC Convergence Trajectory

Visualizes the Robbins-Monro iteration trajectory.

## Usage

``` r
# S3 method for class 'spc_result'
plot(x, type = c("both", "trajectory", "c", "rho"), ...)
```

## Arguments

- x:

  An `spc_result` object.

- type:

  Character. Plot type: `"trajectory"`, `"rho"`, `"c"`, or `"both"`.

- ...:

  Additional arguments (unused).

## Value

A ggplot object (if ggplot2 is available) or NULL (base R fallback).
