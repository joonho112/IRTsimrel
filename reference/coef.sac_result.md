# Extract Calibrated Item Parameters from SAC Results

Returns the calibrated item parameters as a data frame.

## Usage

``` r
# S3 method for class 'sac_result'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"sac_result"`.

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
sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
                          model = "rasch", n_iter = 200, seed = 42)
coef(sac_res)
#>    item_id         beta lambda_base lambda_scaled   c_star
#> 1        1  0.499253218           1      0.973512 0.973512
#> 2        2  0.041044404           1      0.973512 0.973512
#> 3        3  0.596083945           1      0.973512 0.973512
#> 4        4  0.320660189           1      0.973512 0.973512
#> 5        5 -0.448996474           1      0.973512 0.973512
#> 6        6  0.635537612           1      0.973512 0.973512
#> 7        7 -0.040125538           1      0.973512 0.973512
#> 8        8  1.009535379           1      0.973512 0.973512
#> 9        9  1.397349921           1      0.973512 0.973512
#> 10      10 -0.718288593           1      0.973512 0.973512
#> 11      11 -0.724854292           1      0.973512 0.973512
#> 12      12 -1.168821333           1      0.973512 0.973512
#> 13      13  0.338043232           1      0.973512 0.973512
#> 14      14 -0.268043661           1      0.973512 0.973512
#> 15      15 -0.258680387           1      0.973512 0.973512
#> 16      16 -0.973434075           1      0.973512 0.973512
#> 17      17 -1.638272892           1      0.973512 0.973512
#> 18      18 -0.111653878           1      0.973512 0.973512
#> 19      19  1.514997537           1      0.973512 0.973512
#> 20      20 -0.001334315           1      0.973512 0.973512
# }
```
