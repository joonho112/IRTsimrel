% Please edit documentation in R/spc_calibrate.R
\name{predict.sac_result}
\alias{predict.sac_result}
\title{Predict Reliability at New Scaling Factor Values (SAC)}
\usage{
\method{predict}{sac_result}(object, newdata = NULL, theta_vec = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{"sac_result"}.}

\item{newdata}{Optional numeric vector of scaling factor values.
If NULL, returns \code{object$achieved_rho}.}

\item{theta_vec}{Optional numeric vector of abilities. If not provided,
the stored \code{theta_quad} from the SAC result is used with
\code{object$theta_var}. If provided, \code{theta_vec} defines both the
ability sample and the latent-variance basis: reliability is recomputed
with \code{theta_var = stats::var(theta_vec)}.}

\item{...}{Additional arguments (ignored).}
}
\value{
If \code{newdata} is NULL, a single numeric value (achieved reliability).
If \code{newdata} is a numeric vector, a named numeric vector of reliability
values at each scaling factor.
}
\description{
If \code{newdata} is NULL, returns the achieved reliability from calibration.
If \code{newdata} is a numeric vector of scaling factor values, computes
the reliability at each value using the stored item parameters and either
the stored theta sample or a user-supplied theta sample.
By default, prediction uses the stored SAC theta sample with the
pre-calculated latent variance from calibration. If \code{theta_vec} is
supplied, reliability is recomputed on that ability sample and
\code{theta_var} is set to \code{stats::var(theta_vec)}.
}
\examples{
\donttest{
sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
                          model = "rasch", n_iter = 200, seed = 42)
predict(sac_res)
predict(sac_res, newdata = c(0.5, 1.0, 1.5, 2.0))
}

}
