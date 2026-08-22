% Please edit documentation in R/validation_utils.R
\name{compute_reliability_tam}
\alias{compute_reliability_tam}
\title{Compute WLE and EAP Reliability Using TAM}
\usage{
compute_reliability_tam(resp, model = c("rasch", "2pl"), verbose = FALSE, ...)
}
\arguments{
\item{resp}{Matrix or data.frame of item responses (0/1).}

\item{model}{Character. \code{"rasch"} or \code{"2pl"}. A 3PL value is deliberately
unsupported because the package's validated TAM contract requires both
WLE and EAP reliability, while the tested TAM 3PL path does not provide
the required WLE workflow.}

\item{verbose}{Logical. If TRUE, print fitting messages.}

\item{...}{Additional arguments passed to TAM fitting functions.}
}
\value{
A list with components:
\describe{
\item{\code{rel_wle}}{WLE reliability.}
\item{\code{rel_eap}}{EAP reliability.}
\item{\code{mod}}{Fitted TAM model object.}
\item{\code{wle}}{Output from \code{TAM::tam.wle()}.}
}
}
\description{
Fits a Rasch or 2PL model using \pkg{TAM} and computes WLE and EAP
reliability using \code{TAM::WLErel()} and \code{TAM::EAPrel()}. This helper does not
fit or score a 3PL model.
}
\details{
\subsection{WLE vs EAP Reliability}{

TAM defines these reliability coefficients differently:
\itemize{
\item \strong{WLE reliability}: \eqn{1 - \bar{s}^2 / V_{WLE}}, based on design effect
\item \strong{EAP reliability}: \eqn{V_{EAP} / (V_{EAP} + \bar{\sigma}^2)}, based on posterior variance
}

WLE and EAP use different estimators and variance bases, so neither is a
universal upper or lower bound for the other. EAP is the closer external
analogue to the MSEM-based population estimand, but it is not identical to
the information-based EQC estimand.
}

\subsection{3PL limitation}{

\code{compute_reliability_tam()} intentionally rejects \code{model = "3pl"}. The
package's Phase 7 external 3PL oracle used a direct TAM EAP-only fit; TAM 3PL
WLE was unavailable in that validated path. That oracle is evidence about
the probability/information implementation, not a public 3PL WLE API.
}
}
\examples{
\dontrun{
# Simulate response data from calibration results
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  seed = 42
)
sim_data <- simulate_response_data(result = eqc_result, n_persons = 500)

# Compute TAM reliability if TAM is installed
if (requireNamespace("TAM", quietly = TRUE)) {
  tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
  cat(sprintf("WLE reliability: \%.4f\n", tam_rel$rel_wle))
  cat(sprintf("EAP reliability: \%.4f\n", tam_rel$rel_eap))
}
}

}
\seealso{
\code{\link{simulate_response_data}} for generating test data,
\code{\link{eqc_calibrate}} for calibration.
}
