% Please edit documentation in R/validation_utils.R
\name{compare_eqc_sac}
\alias{compare_eqc_sac}
\title{Compare EQC and SAC Calibration Results}
\usage{
compare_eqc_sac(eqc_result, sac_result, verbose = TRUE)
}
\arguments{
\item{eqc_result}{An object of class \code{"eqc_result"}.}

\item{sac_result}{An object of class \code{"sac_result"} (or \code{"spc_result"}
for backward compatibility).}

\item{verbose}{Logical. If TRUE, print comparison summary.}
}
\value{
A list with comparison statistics (invisibly):
\describe{
\item{\code{c_eqc}}{Calibrated c* from EQC.}
\item{\code{c_sac}}{Calibrated c* from SAC.}
\item{\code{diff_abs}}{Absolute difference between c* values.}
\item{\code{diff_pct}}{Percent difference relative to EQC.}
\item{\code{agreement}}{Logical. TRUE if difference is < 5\%.}
\item{\code{target_rho}}{Target reliability.}
\item{\code{achieved_eqc}}{Achieved reliability from EQC.}
\item{\code{achieved_sac}}{Achieved reliability from SAC.}
\item{\code{achieved_diff_abs}}{Absolute achieved-reliability difference.}
\item{\code{metric_eqc}, \code{metric_sac}}{Stored canonical reliability
metrics.}
\item{\code{model_eqc}, \code{model_sac}}{Measurement models.}
\item{\code{n_items_eqc}, \code{n_items_sac}}{Test lengths.}
\item{\code{eqc_status}}{EQC root/calibration status when available.}
\item{\code{sac_status}}{SAC canonical calibration status when available.}
\item{\code{sac_status_flags}}{SAC multi-condition status flags when
available.}
}
}
\description{
Compares the calibration results from EQC and SAC algorithms.
The target reliability must match exactly. Differences in model, test
length, or stored reliability metric are reported as warnings because they
make \code{c*} agreement a weaker diagnostic.

\code{compare_eqc_spc()} is a deprecated backward-compatible alias for
\code{compare_eqc_sac()}.
}
\examples{
\donttest{
# Run both algorithms
eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
                            reliability_metric = "info",
                            c_init = eqc_result, seed = 42)

# Compare results
compare_eqc_sac(eqc_result, sac_result)
}

}
\seealso{
\code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}}
}
