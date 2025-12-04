# Compute WLE and EAP Reliability Using TAM

Fits a Rasch or 2PL model using TAM and computes WLE and EAP reliability
using the official `WLErel()` and `EAPrel()` functions.

## Usage

``` r
compute_reliability_tam(resp, model = c("rasch", "2pl"), verbose = FALSE, ...)
```

## Arguments

- resp:

  Matrix or data.frame of item responses (0/1).

- model:

  Character. `"rasch"` or `"2pl"`.

- verbose:

  Logical. If TRUE, print fitting messages.

- ...:

  Additional arguments passed to TAM fitting functions.

## Value

A list with components:

- `rel_wle`:

  WLE reliability.

- `rel_eap`:

  EAP reliability.

- `mod`:

  Fitted TAM model object.

- `wle`:

  Output from
  [`TAM::tam.wle()`](https://rdrr.io/pkg/TAM/man/tam.wle.html).

## Details

### WLE vs EAP Reliability

TAM defines these reliability coefficients differently:

- **WLE reliability**: \\1 - \bar{s}^2 / V\_{WLE}\\, based on design
  effect

- **EAP reliability**: \\V\_{EAP} / (V\_{EAP} + \bar{\sigma}^2)\\, based
  on posterior variance

Mathematically, \\\rho\_{EAP} \geq \rho\_{WLE}\\ always holds under
TAM's definitions. EAP reliability more closely corresponds to
MSEM-based population reliability. For conservative inference, treat WLE
as a lower bound.

## See also

[`simulate_response_data`](https://joonho112.github.io/IRTsimrel/reference/simulate_response_data.md)
for generating test data,
[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md)
for calibration.

## Examples

``` r
if (FALSE) { # \dontrun{
# Simulate response data from EQC results
sim_data <- simulate_response_data(eqc_result, n_persons = 500)

# Compute TAM reliability
tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
cat(sprintf("WLE reliability: %.4f\n", tam_rel$rel_wle))
cat(sprintf("EAP reliability: %.4f\n", tam_rel$rel_eap))
} # }
```
