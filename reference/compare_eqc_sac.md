# Compare EQC and SAC Calibration Results

Compares the calibration results from EQC and SAC algorithms.

## Usage

``` r
compare_eqc_sac(eqc_result, sac_result, verbose = TRUE)

compare_eqc_spc(...)
```

## Arguments

- eqc_result:

  An object of class `"eqc_result"`.

- sac_result:

  An object of class `"sac_result"` (or `"spc_result"` for backward
  compatibility).

- verbose:

  Logical. If TRUE, print comparison summary.

- ...:

  Arguments passed to `compare_eqc_sac()`.

## Value

A list with comparison statistics (invisibly):

- `c_eqc`:

  Calibrated c\* from EQC.

- `c_sac`:

  Calibrated c\* from SAC.

- `diff_abs`:

  Absolute difference between c\* values.

- `diff_pct`:

  Percent difference relative to EQC.

- `agreement`:

  Logical. TRUE if difference is \< 5%.

- `target_rho`:

  Target reliability.

## See also

[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
[`sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Run both algorithms
eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
                            c_init = eqc_result, seed = 42)

# Compare results
compare_eqc_sac(eqc_result, sac_result)
} # }
```
