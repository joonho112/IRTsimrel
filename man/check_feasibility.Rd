% Please edit documentation in R/reliability_utils.R
\name{check_feasibility}
\alias{check_feasibility}
\alias{print.feasibility_check}
\title{Check Feasibility of Target Reliability}
\usage{
check_feasibility(
  n_items,
  model = c("rasch", "2pl", "3pl"),
  latent_shape = "normal",
  item_source = "parametric",
  c_bounds = c(0.1, 10),
  M = 10000L,
  seed = NULL,
  latent_params = list(),
  item_params = list(),
  target_rho = NULL,
  verbose = TRUE
)

\method{print}{feasibility_check}(x, digits = 4, ...)
}
\arguments{
\item{n_items}{Integer. Number of items in the test form.}

\item{model}{Character. Measurement model: \code{"rasch"}, \code{"2pl"},
or \code{"3pl"}. For 3PL, supply lower-asymptote generation controls in
\code{item_params}; discrimination alone is varied under \code{D = 1}.}

\item{latent_shape}{Character. Shape argument passed to \code{sim_latentG()}.}

\item{item_source}{Character. Source argument passed to \code{sim_item_params()}.}

\item{c_bounds}{Numeric length-2 vector. Range of scaling factors to evaluate.
Default: \code{c(0.1, 10)}.}

\item{M}{Integer. Monte Carlo sample size for theta. Default: 10000.}

\item{seed}{Optional integer for reproducibility.}

\item{latent_params}{List. Additional arguments passed to \code{sim_latentG()}.}

\item{item_params}{List. Additional arguments passed to \code{sim_item_params()}.}

\item{target_rho}{Optional numeric scalar in (0, 1). If supplied, the result
includes status fields indicating whether the target is below, inside, or
above the achievable range for each metric.}

\item{verbose}{Logical. If TRUE, print results.}

\item{x}{An object of class \code{"feasibility_check"}.}

\item{digits}{Integer. Number of decimal places for printing.}

\item{...}{Additional arguments (ignored).}
}
\value{
An object of class \code{"feasibility_check"} (a list) with:
\describe{
\item{\code{rho_range_info}}{Numeric length-2 vector: empirical range of
average-information reliability (\eqn{\tilde{\rho}}) computed by the
configured finite scan; provisional when topology is unresolved.}
\item{\code{rho_range_msem}}{Numeric length-2 vector: empirical range of
MSEM-based reliability (\eqn{\bar{w}}) computed by the configured finite
scan; provisional when topology is unresolved.}
\item{\code{rho_bounds_info}}{Endpoint reliabilities for
average-information reliability.}
\item{\code{rho_bounds_msem}}{Endpoint reliabilities for MSEM-based
reliability.}
\item{\code{rho_info_max_c}}{Scaling factor at the maximum information
reliability detected by the configured scan within \code{c_bounds};
provisional when topology is unresolved.}
\item{\code{rho_msem_max_c}}{Scaling factor at the maximum MSEM
reliability detected by the configured scan within \code{c_bounds};
provisional when topology is unresolved.}
\item{\code{target_status_info}}{If \code{target_rho} is supplied, one of
\code{"below_lower"}, \code{"boundary_lower"}, \code{"feasible"},
\code{"boundary_upper"}, or \code{"above_upper"} for \eqn{\tilde{\rho}}.}
\item{\code{target_status_msem}}{Analogous status for \eqn{\bar{w}}.}
\item{\code{topology_info}, \code{topology_msem}}{Complete configured
log-grid scans, polished extrema, detected roots, monotone branches,
evaluation counts, and resolution status for each metric.}
\item{\code{target_status_info_canonical},
\code{target_status_msem_canonical}}{Canonical feasibility status:
\code{"feasible"}, \code{"infeasible_below_range"},
\code{"infeasible_above_range"}, or \code{"uncertain"}.}
\item{\code{root_count_info}, \code{root_count_msem},
\code{admissible_root_count_info},
\code{admissible_root_count_msem}}{Detected roots and roots on
increasing crossing/boundary branches.}
\item{\code{best_achievable_info},
\code{best_achievable_msem}}{Closest detected point to
\code{target_rho}, including scale, reliability, residual, absolute
error, and location type.}
\item{\code{n_items}}{Number of items.}
\item{\code{model}}{Model used.}
\item{\code{latent_shape}}{Latent distribution shape.}
\item{\code{c_bounds}}{Scaling factor bounds evaluated.}
\item{\code{M}}{Monte Carlo sample size.}
\item{\code{theta_var}}{Estimated latent variance.}
}

The input object, invisibly.
}
\description{
Screens whether a given target reliability is attainable in the generated
finite Monte Carlo design for a particular test design (number of items,
model, latent distribution, item source) by computing the range of
empirical reliabilities across a range of scaling factors.

This function is useful for determining whether a planned simulation
study is feasible before running the (potentially expensive) calibration
algorithms.
}
\details{
Either reliability metric can be non-monotone over a sufficiently broad
finite empirical integration interval. Both metrics therefore use the same
adaptive log-scale topology scan as the calibrators. Endpoint values remain
available in \code{rho_bounds_*}, while \code{rho_range_*} includes all
resolved interior extrema within the user-supplied bounds.

These ranges and target classifications are conditional on the generated
finite theta sample. In particular, for the built-in
\code{latent_shape = "heavy_tail"}, the population MSEM functional can be
non-integrable even though a finite Monte Carlo sample yields a numeric
\code{rho_range_msem}. Treat that MSEM range as finite-grid sensitivity
evidence, not as a population feasibility guarantee.
}
\examples{
# Check feasibility for 25-item Rasch test
feas <- check_feasibility(n_items = 25, model = "rasch",
                          target_rho = 0.90, seed = 42,
                          M = 5000, verbose = FALSE)
print(feas)

# Metric-specific target status
feas$target_status_info
feas$target_status_msem

}
\seealso{
\code{\link{eqc_calibrate}}, \code{\link{rho_curve}}
}
