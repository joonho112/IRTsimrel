% Please edit documentation in R/reliability_utils.R
\name{check_feasibility}
\alias{check_feasibility}
\alias{print.feasibility_check}
\title{Check Feasibility of Target Reliability}
\usage{
check_feasibility(
  n_items,
  model = c("rasch", "2pl"),
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

\item{model}{Character. Measurement model: \code{"rasch"} or \code{"2pl"}.}

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
\item{\code{rho_range_info}}{Numeric length-2 vector: achievable range of
average-information reliability (\eqn{\tilde{\rho}}).}
\item{\code{rho_range_msem}}{Numeric length-2 vector: achievable range of
MSEM-based reliability (\eqn{\bar{w}}).}
\item{\code{rho_bounds_info}}{Endpoint reliabilities for
average-information reliability.}
\item{\code{rho_bounds_msem}}{Endpoint reliabilities for MSEM-based
reliability.}
\item{\code{rho_info_max_c}}{Scaling factor at the maximum sampled info
reliability within \code{c_bounds}.}
\item{\code{rho_msem_max_c}}{Scaling factor at the maximum sampled MSEM
reliability within \code{c_bounds}.}
\item{\code{target_status_info}}{If \code{target_rho} is supplied, one of
\code{"below_lower"}, \code{"boundary_lower"}, \code{"feasible"},
\code{"boundary_upper"}, or \code{"above_upper"} for \eqn{\tilde{\rho}}.}
\item{\code{target_status_msem}}{Analogous status for \eqn{\bar{w}}.}
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
Screens whether a given target reliability is achievable for a particular
test design (number of items, model, latent distribution, item source)
by computing the range of achievable reliabilities across a range of
scaling factors.

This function is useful for determining whether a planned simulation
study is feasible before running the (potentially expensive) calibration
algorithms.
}
\details{
For the average-information metric (\eqn{\tilde{\rho}}), practical EQC
search intervals are usually monotone. With very wide finite empirical
quadrature bounds, however, logistic probabilities can saturate numerically.
The function therefore reports endpoint reliabilities separately and uses an
interior maximization pass for the feasible upper range.

For the MSEM-based metric (\eqn{\bar{w}}), the reliability can be
non-monotone at extreme scaling factors. The function uses
\code{\link[stats]{optimize}} to find the maximum within the bounds, and
the range is \eqn{[\min(\bar{w}(c_{min}), \bar{w}(c_{max})), \bar{w}_{max}]}.
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
