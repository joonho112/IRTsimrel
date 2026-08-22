% Please edit documentation in R/stochastic_calibration.R,
%   R/stochastic_methods.R
\name{sac_calibrate}
\alias{sac_calibrate}
\alias{print.sac_result}
\alias{summary.sac_result}
\alias{print.summary.sac_result}
\alias{coef.sac_result}
\alias{predict.sac_result}
\title{Stochastic Approximation Calibration (Algorithm 2: SAC)}
\usage{
sac_calibrate(
  target_rho,
  n_items,
  model = c("rasch", "2pl", "3pl"),
  latent_shape = "normal",
  item_source = "parametric",
  latent_params = list(),
  item_params = list(),
  reliability_metric = c("msem", "info", "bar", "tilde"),
  c_init = NULL,
  M_per_iter = 500L,
  M_pre = 10000L,
  n_iter = 300L,
  burn_in = NULL,
  step_params = list(),
  c_bounds = c(0.01, 20),
  resample_items = TRUE,
  seed = NULL,
  verbose = FALSE,
  root_policy = "lowest_increasing",
  preflight_controls = list(),
  evaluation_controls = list()
)

\method{print}{sac_result}(x, digits = 4, ...)

\method{summary}{sac_result}(object, ...)

\method{print}{summary.sac_result}(x, digits = 4, ...)

\method{coef}{sac_result}(object, ...)

\method{predict}{sac_result}(object, newdata = NULL, theta_vec = NULL, ...)
}
\arguments{
\item{target_rho}{Finite numeric scalar in \verb{(0, 1)}. Target marginal
reliability.}

\item{n_items}{Positive integer. Number of items per generated form.}

\item{model}{Character. One of \code{"rasch"}, \code{"2pl"}, or \code{"3pl"}.
Rasch fixes baseline discriminations at one. In the 3PL model,
\code{item_params$guessing_params} (or \code{custom_params$guessing} for a custom
source) supplies the lower asymptotes. SAC scales discrimination only and
uses the fixed logistic convention \code{D = 1}.}

\item{latent_shape}{Character shape passed to \code{\link[=sim_latentG]{sim_latentG()}}.}

\item{item_source}{Character source passed to \code{\link[=sim_item_params]{sim_item_params()}}. The
default, \code{"parametric"}, has no external data dependency; \code{"irw"}
requires the optional Item Response Warehouse package.}

\item{latent_params}{Named list of additional arguments passed to
\code{\link[=sim_latentG]{sim_latentG()}}.}

\item{item_params}{Named list of additional arguments passed to
\code{\link[=sim_item_params]{sim_item_params()}}.}

\item{reliability_metric}{Character reliability estimand. \code{"msem"} (the
default) and its synonym \code{"bar"} target reciprocal-information MSEM
reliability, \eqn{\bar w}. \code{"info"} and its synonym \code{"tilde"} target
average-information reliability, \eqn{\tilde\rho}. The stored \code{metric}
is canonicalized to \code{"msem"} or \code{"info"}. Population MSEM is rejected
for the built-in untruncated \code{"heavy_tail"} distribution because its
reciprocal-information expectation is not finite.}

\item{c_init}{\code{NULL}, a positive numeric scalar, or an \code{eqc_result}. With
\code{NULL}, 3PL uses the selected preflight root; Rasch/2PL uses APC when APC
lies inside the selected branch and otherwise uses the preflight root. A
numeric initializer is projected into the effective branch. An EQC warm
start is accepted only when target, model, item count, canonical metric,
fixed-form scope, realized item form, and increasing branch agree.}

\item{M_per_iter}{Positive integer. Ability draws per Robbins--Monro
iteration. Default \code{500L}.}

\item{M_pre}{Integer at least two. Draws used once to estimate the population
latent variance held fixed during calibration. Default \code{10000L}.}

\item{n_iter}{Positive integer. Robbins--Monro iterations. Default \code{300L}.}

\item{burn_in}{Non-negative integer less than \code{n_iter}. Iterations through
\code{burn_in} are excluded from Polyak--Ruppert averaging. The default is
\code{floor(n_iter / 2)}.}

\item{step_params}{Named list overriding \code{a = 1}, \code{A = 50}, or
\code{gamma = 0.67} in \eqn{a_n=a/(n+A)^\gamma}. \code{a} must be positive, \code{A}
non-negative, and \code{gamma} in \verb{(0.5, 1]}.}

\item{c_bounds}{Positive increasing length-two numeric vector. The preflight
scanner never expands these requested bounds. Iterates are projected to
the selected increasing branch after applying its safety margin.}

\item{resample_items}{Logical. \code{FALSE} calibrates one realized fixed form and
sets \code{item_scope = "fixed_form"}. \code{TRUE} (the default) draws a form at each
iteration, targets the item-superpopulation mean, and sets
\code{item_scope = "item_superpopulation"}.}

\item{seed}{Optional integer seed. When supplied, the caller's RNG state is
restored on exit; isolated deterministic substreams are used for
preflight and independent post-calibration evaluation.}

\item{verbose}{Logical or non-negative numeric verbosity level. Level one
reports checkpoints; level two reports every iteration.}

\item{root_policy}{Character root-selection policy: \code{"lowest_increasing"},
\code{"nearest_increasing"}, \code{"lowest_any"}, or \code{"highest_any"}. The default
selects the lowest admissible increasing root. Tangent and plateau roots
are not made admissible by SAC's public controls.}

\item{preflight_controls}{Named list of branch-preflight overrides:
\code{n_forms} (default \code{8L} for item-superpopulation calibration), \code{M}
(default \code{min(M_pre, max(500L, M_per_iter))}), \code{split_check} (default
\code{resample_items}), \code{split_log_tolerance} (default \code{log(1.5)}),
\code{root_controls} (a named list passed to the topology scanner),
\code{branch_margin} (default \code{0.01}), \code{max_consecutive_hits} (default \code{5L}),
and \code{max_post_burn_hit_rate} (default \code{0.10}). Unknown names are rejected.}

\item{evaluation_controls}{Named list controlling the independent final
evaluation: \code{n_forms} (default \code{20L}), \code{M} per form (default
\code{max(1000L, M_per_iter)}), \code{conf_level} (default \code{0.95}), and \code{probs}
(default \code{c(0.025, 0.5, 0.975)}). Unknown names are rejected. For a fixed
form, theta blocks are concatenated before computing \code{achieved_rho}; for
item-superpopulation calibration, form reliabilities are averaged.}

\item{x}{A \code{sac_result} or \code{summary.sac_result}, as appropriate.}

\item{digits}{Number of decimal places used by print methods.}

\item{...}{Additional arguments passed to or from methods. Current SAC
methods do not otherwise use them.}

\item{object}{A \code{sac_result}.}

\item{newdata}{Optional numeric vector of positive discrimination scales for
\code{predict.sac_result()}. \code{NULL} returns the stored \code{achieved_rho}.}

\item{theta_vec}{Optional finite ability vector for prediction. When
supplied, its sample variance defines the prediction variance basis.}
}
\value{
\code{sac_calibrate()} returns a \code{sac_result} list. Its historical fields
remain present and v0.3 metadata are appended:
\describe{
\item{Calibration}{\code{c_star} is the Polyak--Ruppert scale; \code{c_final} is the
last iterate; \code{target_rho}, \code{achieved_rho}, \code{achieved_se},
\code{achieved_distribution}, \code{metric}, and \code{calibration_status} describe the
independently evaluated result. \code{M_final} is the total number of final
evaluation theta draws.}
\item{Trajectories}{\code{trajectory} stores updated scales and
\code{rho_trajectory} stores reliabilities recomputed at those same scales.
\code{evaluation_trajectory}, \code{rho_update_trajectory},
\code{raw_trajectory}, \code{step_size_trajectory}, \code{gradient_trajectory},
\code{projected}, \code{projection_side}, \code{projection_source}, and
\code{iteration_trace} preserve the full update audit trail.}
\item{Algorithm settings}{\code{theta_var}, \code{model}, \code{n_items}, \code{n_iter},
\code{burn_in}, \code{M_per_iter}, \code{M_pre}, \code{step_params}, \code{c_bounds},
\code{requested_c_bounds}, \code{c_init}, \code{init_method}, \code{preflight_controls}, and
\code{evaluation_controls} record the effective run settings.}
\item{Item design}{\code{beta_vec}, \code{lambda_base}, \code{lambda_scaled},
\code{guessing_vec}, \code{items_base}, \code{items_calib}, \code{item_design}, and
\code{item_scope} describe the representative stored form. Rasch/2PL results
store zero guessing; 3PL results preserve the generated lower
asymptotes. For item-superpopulation runs the representative form is not
itself the estimand.}
\item{Contracts and diagnostics}{\code{schema_version}, \code{estimand_signature},
\code{design_signature}, \code{calibration_design}, \code{evaluation_design},
\code{convergence}, \code{preflight}, \code{branch}, \code{rng_provenance}, \code{warm_start},
\code{representative_achieved_rho}, \code{theta_quad}, and \code{call} preserve
estimand, topology, evaluation, RNG, and provenance diagnostics.}
}

\code{print()} returns its input invisibly; \code{summary()} returns a
\code{summary.sac_result}; \code{coef()} returns calibrated item parameters (including
\code{guessing} for 3PL); and \code{predict()} returns stored or recomputed reliability.
}
\description{
\code{sac_calibrate()} uses branch-guarded Robbins--Monro stochastic
approximation to choose a positive global discrimination multiplier that
targets a marginal reliability estimand for a unidimensional Rasch, 2PL, or
3PL item generator.
}
\details{
The first 18 arguments preserve the v0.2 positional API. The v0.3 branch
and evaluation controls are appended after \code{verbose}.
\subsection{Estimand and branch contract}{

SAC first scans the requested log-scale interval and must resolve an
admissible increasing root branch. A classed preflight error is raised when
the target is infeasible, a stable branch cannot be resolved, or only an
inadmissible boundary/decreasing branch is available. The effective branch
bounds are then used as projection guards. Repeated branch-bound contact or
a non-positive final branch slope produces \code{calibration_status = "branch_lost"} rather than silently treating the run as converged.

\code{reliability_metric} and \code{item_scope} are both part of the estimand. Thus an
EQC result (\code{metric = "info"}, fixed form) is directly comparable to SAC
only when SAC also uses \code{reliability_metric = "info"} and
\code{resample_items = FALSE} on that same form. Passing the compatible EQC object
as \code{c_init} without overriding item-generation arguments reuses its form.
}

\subsection{3PL and external validation}{

The 3PL probability and Fisher-information kernels use \code{D = 1}; the global
scale multiplies discrimination and never guessing. The package's public
\code{\link[=compute_reliability_tam]{compute_reliability_tam()}} helper intentionally supports Rasch and 2PL only
because the validated TAM workflow requires both WLE and EAP outputs. Phase
7's 3PL external oracle is an EAP-only validation path; TAM 3PL WLE is not a
supported public-package workflow. Do not interpret WLE and EAP as ordered
versions of the same estimand.

\code{spc_calibrate()} remains a deprecated forwarding alias.
}
}
\examples{
\donttest{
# Small fixed-form 3PL calibration with explicit lower asymptotes.
form <- list(
  custom_params = list(
    beta = c(-1.1, -0.55, -0.1, 0.35, 0.8, 1.25),
    lambda = c(0.75, 0.9, 1.05, 1.2, 1.35, 1.5),
    guessing = c(0.08, 0.12, 0.16, 0.20, 0.24, 0.28)
  ),
  center_difficulties = FALSE
)
fit <- sac_calibrate(
  target_rho = 0.40,
  n_items = 6L,
  model = "3pl",
  item_source = "custom",
  item_params = form,
  reliability_metric = "info",
  M_per_iter = 128L,
  M_pre = 512L,
  n_iter = 8L,
  burn_in = 4L,
  c_bounds = c(0.1, 3),
  resample_items = FALSE,
  seed = 7601L,
  preflight_controls = list(
    M = 256L,
    split_check = FALSE,
    root_controls = list(min_grid = 65L, max_evals = 257L)
  ),
  evaluation_controls = list(n_forms = 4L, M = 512L)
)
print(fit)
coef(fit)
predict(fit, newdata = c(0.5, 1))
}

}
\references{
Robbins, H., & Monro, S. (1951). A stochastic approximation method.
\emph{The Annals of Mathematical Statistics, 22}(3), 400--407.

Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
approximation by averaging. \emph{SIAM Journal on Control and Optimization,
30}(4), 838--855.
}
\seealso{
\code{\link[=eqc_calibrate]{eqc_calibrate()}}, \code{\link[=compute_rho_bar]{compute_rho_bar()}}, \code{\link[=compute_rho_tilde]{compute_rho_tilde()}},
\code{\link[=compare_eqc_sac]{compare_eqc_sac()}}, \code{\link[=sim_item_params]{sim_item_params()}}
}
