% Please edit documentation in R/eqc_calibrate.R
\name{eqc_calibrate}
\alias{eqc_calibrate}
\alias{print.eqc_result}
\alias{summary.eqc_result}
\title{Empirical Quadrature Calibration (Algorithm 1: EQC)}
\usage{
eqc_calibrate(
  target_rho,
  n_items,
  model = c("rasch", "2pl", "3pl"),
  latent_shape = "normal",
  item_source = "parametric",
  latent_params = list(),
  item_params = list(),
  reliability_metric = "info",
  M = 10000L,
  c_bounds = c(0.3, 3),
  tol = 1e-04,
  seed = NULL,
  verbose = FALSE,
  root_policy = "lowest_increasing",
  root_controls = list(),
  allow_tangent = FALSE,
  allow_plateau = FALSE
)

\method{print}{eqc_result}(x, digits = 4, ...)

\method{summary}{eqc_result}(object, ...)
}
\arguments{
\item{target_rho}{Numeric in (0, 1). Target marginal reliability \eqn{\rho^*}.}

\item{n_items}{Integer. Number of items in the test form.}

\item{model}{Character. Measurement model: \code{"rasch"}, \code{"2pl"},
or \code{"3pl"}.
For \code{"rasch"}, all baseline discriminations are set to 1 before scaling.
For \code{"3pl"}, \code{item_params$guessing_params} controls the item
lower asymptotes; calibration scales discrimination only and uses the
fixed logistic convention \code{D = 1}.}

\item{latent_shape}{Character. Shape argument passed to \code{sim_latentG()}
(e.g. \code{"normal"}, \code{"bimodal"}, \code{"heavy_tail"}, ...).}

\item{item_source}{Character. Source argument passed to \code{sim_item_params()}
(e.g. \code{"parametric"}, \code{"irw"}, \code{"hierarchical"}, \code{"custom"}).
Defaults to \code{"parametric"} so core workflows do not require external
item-pool data. Use \code{"irw"} for empirically-grounded difficulties
when the Item Response Warehouse package is installed.}

\item{latent_params}{List. Additional arguments passed to \code{sim_latentG()}.}

\item{item_params}{List. Additional arguments passed to \code{sim_item_params()}.}

\item{reliability_metric}{Character. Reliability definition used inside EQC:
\describe{
\item{\code{"info"}}{Average-information reliability (default, recommended for EQC).
Targets \eqn{\tilde{\rho}}, which is well behaved on practical EQC
search intervals and is bracket-checked before root-finding.}
\item{\code{"tilde"}}{Synonym for \code{"info"}.}
\item{\code{"msem"}, \code{"bar"}}{Rejected by \code{eqc_calibrate()} because
direct MSEM targeting can be non-monotone under EQC root-finding. Use
\code{\link{sac_calibrate}} for direct MSEM-based calibration.}
}}

\item{M}{Integer. Size of the empirical quadrature sample (default: 10000).}

\item{c_bounds}{Numeric length-2 vector. Search bounds for \eqn{c}. Default: c(0.3, 3).}

\item{tol}{Numeric. Root-polishing tolerance. Default: 1e-4.}

\item{seed}{Optional integer for reproducibility.}

\item{verbose}{Logical. If TRUE, print progress messages.}

\item{root_policy}{Character. Root-selection policy used when the fixed-form
reliability curve has multiple roots. The default,
\code{"lowest_increasing"}, selects the lowest-scale root with positive
local slope. \code{"nearest_increasing"} uses the unscaled baseline
\code{c = 1} as its reference. \code{"lowest_any"} and
\code{"highest_any"} select the lowest- or highest-scale admissible root
regardless of slope direction. Tangent and plateau roots remain
inadmissible under the \code{"*_any"} policies unless
\code{allow_tangent} or \code{allow_plateau}, respectively, is enabled.}

\item{root_controls}{Named list of advanced topology-scanner controls. The
supplied \code{c_bounds} are never expanded.}

\item{allow_tangent}{Logical. If \code{TRUE}, a detected tangent root may be
selected by an \code{"*_any"} policy. The default is \code{FALSE}.}

\item{allow_plateau}{Logical. If \code{TRUE}, the geometric midpoint of a
detected target plateau may be selected by an \code{"*_any"} policy. The
complete plateau interval remains recorded in the root table. The default
is \code{FALSE}.}

\item{x}{An object of class \code{"eqc_result"}.}

\item{digits}{Integer. Number of decimal places for printing.}

\item{...}{Additional arguments passed to or from other methods.}

\item{object}{An object of class \code{"eqc_result"}.}
}
\value{
An object of class \code{"eqc_result"} (a list) with elements:
\describe{
\item{\code{c_star}}{Calibrated discrimination scale \eqn{c^*}.}
\item{\code{target_rho}}{Target reliability \eqn{\rho^*}.}
\item{\code{achieved_rho}}{Empirical quadrature estimate \eqn{\hat\rho_M(c^*)}.}
\item{\code{metric}}{Reliability metric used.}
\item{\code{model}}{Measurement model.}
\item{\code{n_items}}{Number of items.}
\item{\code{M}}{Empirical quadrature sample size.}
\item{\code{theta_quad}}{Length-M vector of quadrature abilities.}
\item{\code{theta_var}}{Sample variance of theta_quad.}
\item{\code{beta_vec}}{Item difficulties from the baseline item design.}
\item{\code{lambda_base}}{Baseline (unscaled) item discriminations.}
\item{\code{lambda_scaled}}{Scaled item discriminations
(\code{lambda_base * c_star}).}
\item{\code{guessing_vec}}{Item lower asymptotes. New Rasch/2PL objects
store an all-zero vector; 3PL objects store the generated guessing
parameters.}
\item{\code{items_base}}{item_params object with scale = 1 (baseline).}
\item{\code{items_calib}}{item_params object with discriminations scaled by c_star.}
\item{\code{schema_version}, \code{item_scope}}{Additive result-schema and
fixed-form estimand metadata.}
\item{\code{estimand_signature}, \code{design_signature}}{Canonical
estimand and realized-design provenance used by comparison utilities.}
\item{\code{call}}{Matched function call.}
\item{\code{misc}}{List of root-finding diagnostics, including bounds,
endpoint reliabilities, global achievable range, root topology, branch
selection, and a legacy-compatible \code{uniroot_result} for an interior
crossing.}
}

The input object, invisibly.

An object of class \code{"summary.eqc_result"} containing key
calibration results.
}
\description{
\code{eqc_calibrate()} implements Algorithm 1 (Empirical Quadrature
Calibration, EQC) for reliability-targeted IRT simulation.

Given a target marginal reliability \eqn{\rho^*}, a latent distribution
generator \code{sim_latentG()} (for \eqn{G}) and an item parameter generator
\code{sim_item_params()} (for \eqn{H}), the function searches for a global
discrimination scale \eqn{c^* > 0} such that a selected empirical
quadrature approximation \eqn{\hat\rho_M(c)} for the Rasch, 2PL, or 3PL
model matches \eqn{\rho^*}. Population accuracy is a Monte Carlo
approximation and should be assessed on independent draws.

The key idea is to:
\enumerate{
\item Draw a large fixed "quadrature" sample
\eqn{\{\theta_m\}_{m=1}^M \sim G} and item parameters
\eqn{\{(\beta_i, \lambda_{i,0})\}_{i=1}^I \sim H} once.
\item For any scale \eqn{c}, form \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}
and compute the empirical approximation to population reliability
\eqn{\hat\rho_M(c)} from the test information function.
\item Scan the explicit scale interval on a logarithmic grid, enumerate
all detected roots of \eqn{\hat\rho_M(c^*) = \rho^*}, and apply the
requested branch-selection policy.
}
}
\details{
\subsection{Reliability Metrics}{

The function calibrates the average-information reliability definition:
\itemize{
\item \strong{Average-information} (\code{"info"}/\code{"tilde"}, \strong{default}): Uses the arithmetic mean,
\eqn{\tilde{\rho}(c) = \sigma^2_\theta \bar{\mathcal{J}}(c) / (\sigma^2_\theta \bar{\mathcal{J}}(c) + 1)}.
By Jensen's inequality, \eqn{\tilde{\rho} \geq \bar{w}}, so this metric typically
yields higher reliability values. \strong{This is the recommended default for EQC} because
the objective function \eqn{\tilde{\rho}(c) - \rho^*} is stable on practical
search intervals. The implementation nevertheless scans the complete
user-supplied interval on \code{log(c)}, polishes detected local extrema
and successfully bracketed sign-changing crossings, records any
near-target fallback separately, and applies an explicit root policy. It
does not assume that the empirical curve is globally monotone or unimodal.
}

\strong{MSEM-based} (\code{"msem"}/\code{"bar"}) reliability uses the harmonic mean of test information,
\eqn{\bar{w}(c) = \sigma^2_\theta / (\sigma^2_\theta + E[1/\mathcal{J}(\theta;c)])}.
This metric is rejected by \code{eqc_calibrate()} because the objective
\eqn{\bar{w}(c) - \rho^*} may be non-monotone (see the reliability-targeted
simulation manuscript, Section 4.3). Although the shared scanner can
diagnose such topology, MSEM remains outside EQC's calibration contract.
Use \code{\link{sac_calibrate}} if you need to target \eqn{\bar{w}} directly.
}

\subsection{WLE vs EAP Reliability Interpretation}{

WLE and EAP reliability use different estimators and variance bases, so they
are complementary external diagnostics rather than ordered versions of one
estimand. EAP is more closely related to MSEM-based population reliability;
neither coefficient is identical to EQC's information-based target. The
public \code{\link[=compute_reliability_tam]{compute_reliability_tam()}} helper supports Rasch/2PL only. Phase 7's
direct TAM 3PL check was EAP-only because a validated TAM 3PL WLE path was
unavailable.
}
}
\examples{
# Basic EQC calibration with parametric items (fast)
\donttest{
eqc_result <- eqc_calibrate(
  target_rho = 0.80,
  n_items = 25,
  model = "rasch",
  latent_shape = "normal",
  item_source = "parametric",
  M = 5000L,
  seed = 42
)
print(eqc_result)

# 3PL: guessing is generated once with the fixed item form and is not scaled.
eqc_3pl <- eqc_calibrate(
  target_rho = 0.70,
  n_items = 12L,
  model = "3pl",
  item_params = list(
    guessing_params = list(distribution = "fixed", value = 0.20)
  ),
  M = 1500L,
  c_bounds = c(0.1, 5),
  seed = 5102L
)
stopifnot(all(eqc_3pl$guessing_vec == 0.20))
}

\dontrun{
# EQC with IRW difficulties (requires irw package)
if (requireNamespace("irw", quietly = TRUE)) {
  eqc_result2 <- eqc_calibrate(
    target_rho = 0.80,
    n_items = 25,
    model = "rasch",
    item_source = "irw",
    seed = 42,
    verbose = TRUE
  )
}
}

}
\seealso{
\code{\link{sac_calibrate}} for the stochastic approximation alternative,
\code{\link{compute_rho_bar}} and \code{\link{compute_rho_tilde}} for
reliability computation utilities,
\code{\link{compute_reliability_tam}} for TAM validation.
}
