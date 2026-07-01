# Predict Reliability at New Scaling Factor Values (EQC)

If `newdata` is NULL, returns the achieved reliability from calibration.
If `newdata` is a numeric vector of scaling factor values, computes the
reliability at each value using the stored quadrature sample and item
parameters.

## Usage

``` r
# S3 method for class 'eqc_result'
predict(object, newdata = NULL, ...)
```

## Arguments

- object:

  An object of class `"eqc_result"`.

- newdata:

  Optional numeric vector of scaling factor values. If NULL, returns
  `object$achieved_rho`.

- ...:

  Additional arguments (ignored).

## Value

If `newdata` is NULL, a single numeric value (achieved reliability). If
`newdata` is a numeric vector, a named numeric vector of reliability
values at each scaling factor.

## Examples

``` r
# \donttest{
eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 25,
                          model = "rasch", seed = 42, M = 5000)
predict(eqc_res)
#> [1] 0.8000001
predict(eqc_res, newdata = c(0.5, 1.0, 1.5, 2.0))
#>     c=0.5     c=1.0     c=1.5     c=2.0 
#> 0.5896711 0.8259511 0.8972175 0.9280310 
# }
```
