# Predict Reliability at New Scaling Factor Values (SAC)

If `newdata` is NULL, returns the achieved reliability from calibration.
If `newdata` is a numeric vector of scaling factor values, computes the
reliability at each value using the stored item parameters and either
the stored theta sample or a user-supplied theta sample. By default,
prediction uses the stored SAC theta sample with the pre-calculated
latent variance from calibration. If `theta_vec` is supplied,
reliability is recomputed on that ability sample and `theta_var` is set
to `stats::var(theta_vec)`.

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
  `theta_quad` from the SAC result is used with `object$theta_var`. If
  provided, `theta_vec` defines both the ability sample and the
  latent-variance basis: reliability is recomputed with
  `theta_var = stats::var(theta_vec)`.

- ...:

  Additional arguments (ignored).

## Value

If `newdata` is NULL, a single numeric value (achieved reliability). If
`newdata` is a numeric vector, a named numeric vector of reliability
values at each scaling factor.

## Examples

``` r
# \donttest{
sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
                          model = "rasch", n_iter = 200, seed = 42)
predict(sac_res)
#> [1] 0.7780128
predict(sac_res, newdata = c(0.5, 1.0, 1.5, 2.0))
#>     c=0.5     c=1.0     c=1.5     c=2.0 
#> 0.5349754 0.7847162 0.8588747 0.8857523 
# }
```
