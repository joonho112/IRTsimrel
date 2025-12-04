# Compare EQC and SPC Calibration Results

Compares the calibration results from EQC and SPC algorithms.

## Usage

``` r
compare_eqc_spc(eqc_result, spc_result, verbose = TRUE)
```

## Arguments

- eqc_result:

  An object of class `"eqc_result"`.

- spc_result:

  An object of class `"spc_result"`.

- verbose:

  Logical. If TRUE, print comparison summary.

## Value

A list with comparison statistics (invisibly):

- `c_eqc`:

  Calibrated c\* from EQC.

- `c_spc`:

  Calibrated c\* from SPC.

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
[`spc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/spc_calibrate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Run both algorithms
eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
spc_result <- spc_calibrate(target_rho = 0.80, n_items = 25,
                            c_init = eqc_result, seed = 42)

# Compare results
compare_eqc_spc(eqc_result, spc_result)
} # }
```
