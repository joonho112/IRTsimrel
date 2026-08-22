% Please edit documentation in R/validation_utils.R
\name{simulate_response_data}
\alias{simulate_response_data}
\title{Simulate Item Response Data from Calibration Results}
\usage{
simulate_response_data(
  result,
  n_persons,
  latent_shape = "normal",
  latent_params = list(),
  seed = NULL
)
}
\arguments{
\item{result}{A calibration result object of class \code{"eqc_result"},
\code{"sac_result"}, or \code{"spc_result"} (for backward compatibility),
as returned by \code{eqc_calibrate()} or \code{sac_calibrate()}.}

\item{n_persons}{Integer. Number of persons to simulate.}

\item{latent_shape}{Character. Shape argument for \code{sim_latentG()}.}

\item{latent_params}{List. Additional arguments for \code{sim_latentG()}.}

\item{seed}{Optional integer for reproducibility.}
}
\value{
A list containing:
\describe{
\item{\code{response_matrix}}{N x I matrix of binary responses}
\item{\code{theta}}{True abilities (N x 1)}
\item{\code{beta}}{Item difficulties (I x 1)}
\item{\code{lambda}}{Item discriminations (I x 1)}
\item{\code{provenance}}{List describing the calibration result and
simulation settings used to generate the response data, including
result class, target and achieved reliability, metric, model, test
length, scalar \code{calibration_status}, character-vector
\code{status_flags}, item source/design, calibration call, simulation
seed, sample size, latent shape, and latent parameters.}
\item{\code{guessing}}{Item lower asymptotes (I x 1). This additive
component is a zero vector for Rasch and 2PL designs.}
}
}
\description{
Generates item response data using the calibrated parameters from
\code{eqc_calibrate()} or \code{sac_calibrate()} (formerly
\code{spc_calibrate()}).

The function uses the normalized item design stored in the calibration
result. For a fixed-form result this is the calibrated fixed form. For a
SAC result with \code{item_scope = "item_superpopulation"}, it is one stored
representative evaluation form; this function does not redraw a new item
form or simulate the aggregate item-superpopulation estimand.
}
\examples{
\dontrun{
# Example 1: Using EQC calibration result
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42
)

sim_data <- simulate_response_data(
  result = eqc_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Example 2: Using SAC calibration result
sac_result <- sac_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  reliability_metric = "info",
  c_init = eqc_result,
  resample_items = FALSE,
  n_iter = 200,
  seed = 42
)

sim_data2 <- simulate_response_data(
  result = sac_result,
  n_persons = 1000,
  latent_shape = "normal",
  seed = 123
)

# Use with TAM for validation when TAM is installed
if (requireNamespace("TAM", quietly = TRUE)) {
  tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
}
}

}
\seealso{
\code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}},
\code{\link{compute_reliability_tam}}
}
