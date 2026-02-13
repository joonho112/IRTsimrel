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
if (FALSE) { # \dontrun{
sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
                          model = "rasch", n_iter = 200, seed = 42)
coef(sac_res)
} # }
```
