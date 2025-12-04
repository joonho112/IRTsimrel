# Simulate Item Response Data from EQC Results

Generates item response data using the calibrated parameters from EQC.

## Usage

``` r
simulate_response_data(
  eqc_result,
  n_persons,
  latent_shape = "normal",
  latent_params = list(),
  seed = NULL
)
```

## Arguments

- eqc_result:

  An `eqc_result` object from
  [`eqc_calibrate()`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md).

- n_persons:

  Integer. Number of persons to simulate.

- latent_shape:

  Character. Shape argument for
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- latent_params:

  List. Additional arguments for
  [`sim_latentG()`](https://joonho112.github.io/IRTsimrel/reference/sim_latentG.md).

- seed:

  Optional integer for reproducibility.

## Value

A list containing:

- `response_matrix`:

  N x I matrix of binary responses

- `theta`:

  True abilities (N x 1)

- `beta`:

  Item difficulties (I x 1)

- `lambda`:

  Item discriminations (I x 1)

## See also

[`eqc_calibrate`](https://joonho112.github.io/IRTsimrel/reference/eqc_calibrate.md),
[`compute_reliability_tam`](https://joonho112.github.io/IRTsimrel/reference/compute_reliability_tam.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# First, run EQC calibration
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42
)

# Generate response data
sim_data <- simulate_response_data(
  eqc_result = eqc_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Use with TAM for validation
tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
} # }
```
