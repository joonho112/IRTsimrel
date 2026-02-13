# Predict Reliability at New Scaling Factor Values (SAC)

If `newdata` is NULL, returns the achieved reliability from calibration.
If `newdata` is a numeric vector of scaling factor values, computes the
reliability at each value using the stored theta sample and item
parameters.

## Usage

``` r
# S3 method for class 'sac_result'
predict(object, newdata = NULL, theta_vec = NULL, ...)
```

## Arguments

- object:

  An object of class `"sac_result"`.

- newdata:

  Optional numeric vector of scaling factor values. If NULL, returns
  `object$achieved_rho`.

- theta_vec:

  Optional numeric vector of abilities. If not provided, the stored
  `theta_quad` from the SAC result is used.

- ...:

  Additional arguments (ignored).

## Value

If `newdata` is NULL, a single numeric value (achieved reliability). If
`newdata` is a numeric vector, a named numeric vector of reliability
values at each scaling factor.

## Examples

``` r
if (FALSE) { # \dontrun{
sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
                          model = "rasch", n_iter = 200, seed = 42)
predict(sac_res)
predict(sac_res, newdata = c(0.5, 1.0, 1.5, 2.0))
} # }
```
