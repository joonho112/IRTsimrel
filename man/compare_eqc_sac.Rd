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
\item{\code{agreement}}{Logical. \code{TRUE} if the objects are comparable,
both calibrations are valid for agreement, and the difference is below
5 percent; \code{FALSE} if evaluated and above threshold; \code{NA} when
agreement is withheld.}
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
\item{\code{comparable}}{Logical indicating whether both objects target
the same estimand and design.}
\item{\code{comparability_reasons}}{Stable reason codes explaining a
non-comparable result.}
\item{\code{agreement_status}, \code{agreement_reasons}}{Whether agreement
was evaluated, and why it was withheld when it was not.}
\item{\code{eqc_contract}, \code{sac_contract}}{Normalized estimand and
design contracts used for the comparability decision.}
}
}
\description{
Compares the calibration results from EQC and SAC algorithms. The target
reliability must match exactly. Numeric estimates are always shown side by
side, but agreement is evaluated only when the stored estimand and design
contracts match: model, reliability metric, theta population, latent
specification, item scope, test length, scale convention, and either the
complete fixed item bank or normalized superpopulation generator.

\code{compare_eqc_spc()} is a deprecated backward-compatible alias for
\code{compare_eqc_sac()}.
}
\examples{
\donttest{
# Run both algorithms
eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
                            reliability_metric = "info",
                            c_init = eqc_result,
                            resample_items = FALSE,
                            seed = 42)

# Compare results
compare_eqc_sac(eqc_result, sac_result)
}

}
\seealso{
\code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}}
}
