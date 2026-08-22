% Please edit documentation in R/reliability_utils.R
\name{compute_rho_both}
\alias{compute_rho_both}
\title{Compute Both Reliability Metrics in a Single Pass}
\usage{
compute_rho_both(
  c,
  theta_vec,
  beta_vec,
  lambda_base,
  theta_var = NULL,
  guessing = NULL,
  weights = NULL,
  return_diagnostics = FALSE
)
}
\arguments{
\item{c}{Numeric scalar. Global discrimination scaling factor.}

\item{theta_vec}{Numeric vector of abilities \eqn{\theta_m}.}

\item{beta_vec}{Numeric vector of item difficulties \eqn{\beta_i}.}

\item{lambda_base}{Numeric vector of baseline discriminations
\eqn{\lambda_{i,0}} (before scaling by \code{c}).}

\item{theta_var}{Optional numeric. Pre-calculated variance of theta.
If NULL, computed from \code{theta_vec}.}

\item{guessing}{Optional numeric scalar or item-length vector of 3PL lower
asymptotes in \verb{[0, 1)}. A scalar is recycled; \code{NULL} gives the Rasch/2PL
special case. The kernel uses \code{D = 1}.}

\item{weights}{Optional non-negative integration weights, one per theta
value. They are normalized internally.}

\item{return_diagnostics}{Logical. Return log-domain numerical diagnostics
in addition to the two legacy reliability fields.}
}
\value{
A named list with components:
\describe{
\item{\code{rho_tilde}}{Average-information reliability.}
\item{\code{rho_bar}}{MSEM-based marginal reliability.}
}
}
\description{
Computes both the average-information reliability (\eqn{\tilde{\rho}}) and
the MSEM-based marginal reliability (\eqn{\bar{w}}) from a single set of
test-information reductions, avoiding redundant kernel evaluation.

This is a performance optimization over calling \code{compute_rho_tilde()}
and \code{compute_rho_bar()} separately, since both share the same stable,
automatically chunked item-information reduction.
}
\details{
For the same information values and latent-variance basis, Jensen's
inequality implies \eqn{\tilde{\rho} \geq \bar{w}}.
See \code{\link{compute_rho_bar}} for details on each metric.
}
\examples{
set.seed(1)
theta <- rnorm(1000)
beta  <- rnorm(20)
lambda0 <- rep(1, 20)

both <- compute_rho_both(1, theta, beta, lambda0)
both$rho_tilde
both$rho_bar

# Verify: rho_tilde >= rho_bar (Jensen's inequality)
both$rho_tilde >= both$rho_bar

# The same call supports 3PL information.
compute_rho_both(1, theta, beta, lambda0, guessing = 0.20)

}
