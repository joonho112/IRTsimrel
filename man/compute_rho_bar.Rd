% Please edit documentation in R/reliability_utils.R
\name{compute_rho_bar}
\alias{compute_rho_bar}
\alias{compute_rho_tilde}
\title{Compute Marginal Reliability from Simulated Test Information}
\usage{
compute_rho_bar(
  c,
  theta_vec,
  beta_vec,
  lambda_base,
  theta_var = NULL,
  guessing = NULL,
  weights = NULL,
  return_diagnostics = FALSE
)

compute_rho_tilde(
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

\item{guessing}{Optional numeric scalar or item-length vector of 3PL
lower-asymptote parameters in \verb{[0, 1)}. A scalar is recycled. \code{NULL} is
equivalent to zero guessing and gives the Rasch/2PL special case.}

\item{weights}{Optional non-negative integration weights, one per value in
\code{theta_vec}. Weights are normalized internally.}

\item{return_diagnostics}{Logical. If \code{FALSE} (the default), return the
legacy scalar result. If \code{TRUE}, return the reliability together
with numerical diagnostics and estimand metadata.}
}
\value{
By default, a numeric scalar reliability value. With
\code{return_diagnostics = TRUE}, a named list containing the value,
logit-scale value, numerical summaries, and estimand metadata.
}
\description{
Low-level utilities for mapping a discrimination scale \eqn{c} and
simulated person/item parameters to marginal reliability.

These functions implement the same reliability definitions used
inside both \code{eqc_calibrate()} and \code{sac_calibrate()}:

\itemize{
\item \code{compute_rho_bar()}: MSEM-based marginal reliability
\eqn{\bar{w}(c) = \sigma_\theta^2 / (\sigma_\theta^2 + E[1/\mathcal{J}(\theta;c)])}
\item \code{compute_rho_tilde()}: Average-information reliability
\eqn{\tilde{\rho}(c) = \sigma_\theta^2 \bar{\mathcal{J}}(c) /
      (\sigma_\theta^2 \bar{\mathcal{J}}(c) + 1)}
}
}
\details{
The computation proceeds as:
\enumerate{
\item Form scaled discriminations \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}.
\item With lower asymptote \eqn{g_i}, compute
\eqn{p_{mi}=g_i+(1-g_i)\operatorname{logit}^{-1}
      \{\lambda_i(c)(\theta_m-\beta_i)\}} under the fixed \code{D = 1}
convention. Setting every \eqn{g_i=0} gives Rasch/2PL exactly.
\item Compute 3PL item information
\eqn{\mathcal{J}_{mi}=\lambda_i(c)^2
      (p_{mi}-g_i)^2(1-p_{mi})/\{(1-g_i)^2p_{mi}\}}. For \eqn{g_i=0},
this reduces to \eqn{\lambda_i(c)^2p_{mi}(1-p_{mi})}.
\item Test information at each \eqn{\theta_m}:
\eqn{\mathcal{J}_m = \sum_i \mathcal{J}_{mi}}.
\item Reliability:
\itemize{
\item \code{compute_rho_bar()}: harmonic-mean-based MSEM
\eqn{\text{MSEM} = E[1/\mathcal{J}_m]},
\eqn{\bar{w}(c) = \sigma_\theta^2 / (\sigma_\theta^2 + \text{MSEM})}.
\item \code{compute_rho_tilde()}: arithmetic-mean-based information
\eqn{\bar{\mathcal{J}} = E[\mathcal{J}_m]},
\eqn{\tilde{\rho}(c) = \sigma_\theta^2 \bar{\mathcal{J}} /
          (\sigma_\theta^2 \bar{\mathcal{J}} + 1)}.
}
}

Information and its expectations are evaluated in the log domain. No
artificial information floor is applied. Consequently, zero information
at a positively weighted node correctly implies infinite MSEM and zero
MSEM-based reliability.
}
\examples{
# Simple toy example
set.seed(1)
theta <- rnorm(1000)
beta  <- rnorm(20)
lambda0 <- rep(1, 20)

compute_rho_bar(1, theta, beta, lambda0)
compute_rho_tilde(1, theta, beta, lambda0)

# With pre-calculated theta variance (recommended for SAC)
theta_var_fixed <- var(rnorm(10000))  # Pre-calculate from large sample
compute_rho_bar(1, theta, beta, lambda0, theta_var = theta_var_fixed)

# 3PL lower asymptotes use the same public kernel.
g <- rep(0.20, length(beta))
compute_rho_tilde(1, theta, beta, lambda0, guessing = g)

}
