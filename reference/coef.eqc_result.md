# Extract Calibrated Item Parameters from EQC Results

Returns the calibrated item parameters as a data frame.

## Usage

``` r
# S3 method for class 'eqc_result'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"eqc_result"`.

- ...:

  Additional arguments (ignored).

## Value

A data frame with columns:

- `item_id`:

  Item identifier (1 to I).

- `beta`:

  Item difficulty.

- `lambda_base`:

  Baseline (unscaled) discrimination.

- `lambda_scaled`:

  Scaled discrimination (`lambda_base * c*`).

- `c_star`:

  Calibrated scaling factor (same for all items).

## Examples

``` r
# \donttest{
eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 25,
                          model = "rasch", seed = 42, M = 5000)
coef(eqc_res)
#>    item_id         beta lambda_base lambda_scaled    c_star
#> 1        1  0.197732269           1     0.8994993 0.8994993
#> 2        2  1.096799859           1     0.8994993 0.8994993
#> 3        3  0.436545084           1     0.8994993 0.8994993
#> 4        4 -0.013038730           1     0.8994993 0.8994993
#> 5        5 -0.199801302           1     0.8994993 0.8994993
#> 6        6  0.007700326           1     0.8994993 0.8994993
#> 7        7  1.020068726           1     0.8994993 0.8994993
#> 8        8  0.337624343           1     0.8994993 0.8994993
#> 9        9 -0.362269961           1     0.8994993 0.8994993
#> 10      10 -0.093862093           1     0.8994993 0.8994993
#> 11      11 -2.169735948           1     0.8994993 0.8994993
#> 12      12  0.403422109           1     0.8994993 0.8994993
#> 13      13  1.100326100           1     0.8994993 0.8994993
#> 14      14 -1.010084428           1     0.8994993 0.8994993
#> 15      15  0.222591967           1     0.8994993 0.8994993
#> 16      16 -0.358839904           1     0.8994993 0.8994993
#> 17      17 -0.661889194           1     0.8994993 0.8994993
#> 18      18  0.214349687           1     0.8994993 0.8994993
#> 19      19 -0.615836104           1     0.8994993 0.8994993
#> 20      20  0.049153126           1     0.8994993 0.8994993
#> 21      21  1.449782647           1     0.8994993 0.8994993
#> 22      22 -1.574952086           1     0.8994993 0.8994993
#> 23      23 -0.546968292           1     0.8994993 0.8994993
#> 24      24 -0.303196880           1     0.8994993 0.8994993
#> 25      25  1.374378677           1     0.8994993 0.8994993
# }
```
