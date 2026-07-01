# Compare EQC and SAC Calibration Results

Compares the calibration results from EQC and SAC algorithms. The target
reliability must match exactly. Differences in model, test length, or
stored reliability metric are reported as warnings because they make
`c*` agreement a weaker diagnostic.

[`compare_eqc_spc()`](https://joonho112.github.io/IRTsimrel/reference/compare_eqc_spc.md)
is a deprecated backward-compatible alias for `compare_eqc_sac()`.

## Usage

``` r
compare_eqc_sac(eqc_result, sac_result, verbose = TRUE)
```

## Arguments

- eqc_result:

  An object of class `"eqc_result"`.

- sac_result:

  An object of class `"sac_result"` (or `"spc_result"` for backward
  compatibility).

- verbose:

  Logical. If TRUE, print comparison summary.

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

- `achieved_eqc`:

  Achieved reliability from EQC.

- `achieved_sac`:

  Achieved reliability from SAC.

- `achieved_diff_abs`:

  Absolute achieved-reliability difference.

- `metric_eqc`, `metric_sac`:

  Stored canonical reliability metrics.

- `model_eqc`, `model_sac`:

  Measurement models.

- `n_items_eqc`, `n_items_sac`:

  Test lengths.

- `eqc_status`:

  EQC root/calibration status when available.

- `sac_status`:

  SAC canonical calibration status when available.

- `sac_status_flags`:

  SAC multi-condition status flags when available.

## See also

[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
[`sac_calibrate`](https://joonho112.github.io/IRTsimrel/reference/sac_calibrate.md)

## Examples

``` r
# \donttest{
# Run both algorithms
eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
                            reliability_metric = "info",
                            c_init = eqc_result, seed = 42)

# Compare results
compare_eqc_sac(eqc_result, sac_result)
#> 
#> =======================================================
#>   EQC vs SAC Comparison
#> =======================================================
#> 
#>   Target reliability  : 0.8000
#>   EQC c*              : 0.905885
#>   SAC c*              : 0.912169
#>   Absolute difference : 0.006284
#>   Percent difference  : 0.69%
#>   EQC achieved rho    : 0.8000
#>   SAC achieved rho    : 0.7994
#>   Agreement (< 5%)    : YES
#>   EQC status          : uniroot_success
#>   SAC status          : ok
#>   SAC status flags    : ok
#> 
# }
```
