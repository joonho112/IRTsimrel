# =============================================================================
# reliability_utils.R
# =============================================================================
# Exported Low-Level Reliability Computation Functions
# These functions can be shared by EQC and SAC for consistency.
#
# Contents:
#   - compute_rho_bar(): MSEM-based marginal reliability
#   - compute_rho_tilde(): Average-information reliability
#   - .compute_rho_generic(): Internal engine (not exported)
#   - compute_apc_init(): Analytic Pre-Calibration initialization
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# =============================================================================


#' Compute Marginal Reliability from Simulated Test Information
#'
#' @description
#' Low-level utilities for mapping a discrimination scale \eqn{c} and
#' simulated person/item parameters to marginal reliability.
#'
#' These functions implement the same reliability definitions used
#' inside both `eqc_calibrate()` and `sac_calibrate()`:
#'
#' \itemize{
#'   \item \code{compute_rho_bar()}: MSEM-based marginal reliability
#'     \eqn{\bar{w}(c) = \sigma_\theta^2 / (\sigma_\theta^2 + E[1/\mathcal{J}(\theta;c)])}
#'   \item \code{compute_rho_tilde()}: Average-information reliability
#'     \eqn{\tilde{\rho}(c) = \sigma_\theta^2 \bar{\mathcal{J}}(c) /
#'       (\sigma_\theta^2 \bar{\mathcal{J}}(c) + 1)}
#' }
#'
#' @param c Numeric scalar. Global discrimination scaling factor.
#' @param theta_vec Numeric vector of abilities \eqn{\theta_m}.
#' @param beta_vec Numeric vector of item difficulties \eqn{\beta_i}.
#' @param lambda_base Numeric vector of baseline discriminations
#'   \eqn{\lambda_{i,0}} (before scaling by \code{c}).
#' @param theta_var Optional numeric. Pre-calculated variance of theta.
#'   If NULL, computed from \code{theta_vec}.
#' @param guessing Optional numeric scalar or item-length vector of 3PL
#'   lower-asymptote parameters in `[0, 1)`. A scalar is recycled. `NULL` is
#'   equivalent to zero guessing and gives the Rasch/2PL special case.
#' @param weights Optional non-negative integration weights, one per value in
#'   \code{theta_vec}. Weights are normalized internally.
#' @param return_diagnostics Logical. If \code{FALSE} (the default), return the
#'   legacy scalar result. If \code{TRUE}, return the reliability together
#'   with numerical diagnostics and estimand metadata.
#'
#' @return By default, a numeric scalar reliability value. With
#'   \code{return_diagnostics = TRUE}, a named list containing the value,
#'   logit-scale value, numerical summaries, and estimand metadata.
#'
#' @details
#' The computation proceeds as:
#' \enumerate{
#'   \item Form scaled discriminations \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}.
#'   \item With lower asymptote \eqn{g_i}, compute
#'     \eqn{p_{mi}=g_i+(1-g_i)\operatorname{logit}^{-1}
#'       \{\lambda_i(c)(\theta_m-\beta_i)\}} under the fixed `D = 1`
#'     convention. Setting every \eqn{g_i=0} gives Rasch/2PL exactly.
#'   \item Compute 3PL item information
#'     \eqn{\mathcal{J}_{mi}=\lambda_i(c)^2
#'       (p_{mi}-g_i)^2(1-p_{mi})/\{(1-g_i)^2p_{mi}\}}. For \eqn{g_i=0},
#'     this reduces to \eqn{\lambda_i(c)^2p_{mi}(1-p_{mi})}.
#'   \item Test information at each \eqn{\theta_m}:
#'     \eqn{\mathcal{J}_m = \sum_i \mathcal{J}_{mi}}.
#'   \item Reliability:
#'     \itemize{
#'       \item \code{compute_rho_bar()}: harmonic-mean-based MSEM
#'         \eqn{\text{MSEM} = E[1/\mathcal{J}_m]},
#'         \eqn{\bar{w}(c) = \sigma_\theta^2 / (\sigma_\theta^2 + \text{MSEM})}.
#'       \item \code{compute_rho_tilde()}: arithmetic-mean-based information
#'         \eqn{\bar{\mathcal{J}} = E[\mathcal{J}_m]},
#'         \eqn{\tilde{\rho}(c) = \sigma_\theta^2 \bar{\mathcal{J}} /
#'           (\sigma_\theta^2 \bar{\mathcal{J}} + 1)}.
#'     }
#' }
#'
#' Information and its expectations are evaluated in the log domain. No
#' artificial information floor is applied. Consequently, zero information
#' at a positively weighted node correctly implies infinite MSEM and zero
#' MSEM-based reliability.
#'
#' @examples
#' # Simple toy example
#' set.seed(1)
#' theta <- rnorm(1000)
#' beta  <- rnorm(20)
#' lambda0 <- rep(1, 20)
#'
#' compute_rho_bar(1, theta, beta, lambda0)
#' compute_rho_tilde(1, theta, beta, lambda0)
#'
#' # With pre-calculated theta variance (recommended for SAC)
#' theta_var_fixed <- var(rnorm(10000))  # Pre-calculate from large sample
#' compute_rho_bar(1, theta, beta, lambda0, theta_var = theta_var_fixed)
#'
#' # 3PL lower asymptotes use the same public kernel.
#' g <- rep(0.20, length(beta))
#' compute_rho_tilde(1, theta, beta, lambda0, guessing = g)
#'
#' @export
compute_rho_bar <- function(c,
                            theta_vec,
                            beta_vec,
                            lambda_base,
                            theta_var = NULL,
                            guessing = NULL,
                            weights = NULL,
                            return_diagnostics = FALSE) {
  .compute_rho_generic(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    theta_var = theta_var,
    metric_internal = "msem",
    guessing = guessing,
    weights = weights,
    return_diagnostics = return_diagnostics
  )
}


#' @rdname compute_rho_bar
#' @export
compute_rho_tilde <- function(c,
                              theta_vec,
                              beta_vec,
                              lambda_base,
                              theta_var = NULL,
                              guessing = NULL,
                              weights = NULL,
                              return_diagnostics = FALSE) {
  .compute_rho_generic(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    theta_var = theta_var,
    metric_internal = "info",
    guessing = guessing,
    weights = weights,
    return_diagnostics = return_diagnostics
  )
}


#' Internal: Normalize Reliability Weights
#'
#' @noRd
.irtsimrel_reducer_weights <- function(weights, n_nodes) {
  if (is.null(weights)) {
    normalized <- rep.int(1 / n_nodes, n_nodes)
    attr(normalized, "log_weights") <- rep.int(-log(n_nodes), n_nodes)
    return(normalized)
  }
  if (!is.numeric(weights) || length(weights) != n_nodes ||
      anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be a finite, non-negative numeric vector with length(theta_vec).")
  }
  positive <- weights > 0
  if (!any(positive)) {
    stop("weights must have a positive finite sum.")
  }

  # Normalize in the log domain. Dividing first in the natural domain can
  # erase a positive node whenever the raw weight ratio is below the smallest
  # representable double, even though that node can still contribute after it
  # is combined with an extreme log-information value.
  log_raw <- rep.int(-Inf, n_nodes)
  log_raw[positive] <- log(weights[positive])
  pivot <- max(log_raw[positive])
  log_total <- pivot + log(sum(exp(log_raw[positive] - pivot)))
  log_normalized <- log_raw - log_total
  normalized <- exp(log_normalized)
  attr(normalized, "log_weights") <- log_normalized
  normalized
}


#' Internal: Recover Normalized Log Weights
#'
#' @noRd
.irtsimrel_reducer_log_weights <- function(weights) {
  log_weights <- attr(weights, "log_weights", exact = TRUE)
  if (!is.null(log_weights) && length(log_weights) == length(weights) &&
      !anyNA(log_weights)) {
    return(as.numeric(log_weights))
  }

  log_weights <- rep.int(-Inf, length(weights))
  positive <- weights > 0
  log_weights[positive] <- log(weights[positive])
  log_weights
}


#' Internal: Stable Signed Weighted Sum
#'
#' @noRd
.irtsimrel_reducer_signed_log_sum <- function(values, log_weights) {
  active <- is.finite(log_weights) & values != 0
  if (!any(active)) return(0)

  log_terms <- log_weights[active] + log(abs(values[active]))
  pivot <- max(log_terms)
  scaled <- sum(sign(values[active]) * exp(log_terms - pivot))
  if (scaled == 0) return(0)

  sign(scaled) * exp(pivot + log(abs(scaled)))
}


#' Internal: Stable Log Absolute Difference
#'
#' @noRd
.irtsimrel_reducer_log_abs_diff <- function(x, center) {
  vapply(x, function(value) {
    if (identical(value, center) || value == center) return(-Inf)
    abs_value <- abs(value)
    abs_center <- abs(center)
    high <- max(abs_value, abs_center)
    low <- min(abs_value, abs_center)
    if (high == 0) return(-Inf)
    adjustment <- if (sign(value) == sign(center)) {
      log1p(-low / high)
    } else {
      log1p(low / high)
    }
    log(high) + adjustment
  }, numeric(1))
}


#' Internal: Resolve Theta Variance
#'
#' @noRd
.resolve_theta_var <- function(theta_vec, theta_var = NULL, weights = NULL) {
  source <- "supplied"
  if (is.null(theta_var)) {
    if (is.null(weights)) {
      theta_var <- stats::var(theta_vec)
      source <- "sample_variance"
    } else {
      weights_norm <- if (!is.null(attr(weights, "log_weights", exact = TRUE))) {
        weights
      } else {
        .irtsimrel_reducer_weights(weights, length(theta_vec))
      }
      log_weights <- .irtsimrel_reducer_log_weights(weights_norm)
      theta_mean <- .irtsimrel_reducer_signed_log_sum(theta_vec, log_weights)
      log_deviation <- .irtsimrel_reducer_log_abs_diff(theta_vec, theta_mean)
      log_theta_var <- .irtsimrel_reducer_log_sum_exp(
        log_weights + 2 * log_deviation
      )
      theta_var <- exp(log_theta_var)
      source <- "weighted_population_variance"
    }
  }

  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 0) {
    stop("theta_var must be a finite positive numeric scalar.")
  }

  structure(as.numeric(theta_var), source = source)
}


#' Internal: Compute Test Information
#'
#' @noRd
.compute_test_information <- function(c,
                                      theta_vec,
                                      beta_vec,
                                      lambda_base,
                                      eta_base = NULL,
                                      guessing = NULL,
                                      chunk_size = NULL) {
  if (!is.null(eta_base)) {
    eta_base <- as.matrix(eta_base)
    if (!identical(dim(eta_base), c(length(theta_vec), length(beta_vec)))) {
      stop("eta_base must have dimensions length(theta_vec) by length(beta_vec).")
    }
  }

  bank <- .irtsimrel_item_bank(
    beta = beta_vec,
    lambda_base = lambda_base,
    guessing = guessing
  )
  .irtsimrel_test_information(
    theta = theta_vec,
    bank = bank,
    scale = c,
    chunk_size = chunk_size
  )
}


#' Internal: Validate Reliability Inputs
#'
#' @noRd
.validate_reliability_inputs <- function(c,
                                         theta_vec,
                                         beta_vec,
                                         lambda_base,
                                         guessing = NULL,
                                         weights = NULL) {
  if (!is.numeric(c) || length(c) != 1L || is.na(c) || !is.finite(c)) {
    stop("c must be a finite numeric scalar.")
  }
  if (c < 0) {
    stop("c must be non-negative.")
  }
  if (!is.numeric(theta_vec) || length(theta_vec) < 2L || anyNA(theta_vec) ||
      any(!is.finite(theta_vec))) {
    stop("theta_vec must be a finite numeric vector with length >= 2.")
  }
  if (!is.numeric(beta_vec) || length(beta_vec) == 0L || anyNA(beta_vec) ||
      any(!is.finite(beta_vec))) {
    stop("beta_vec must be a non-empty finite numeric vector.")
  }
  if (!is.numeric(lambda_base) || length(lambda_base) == 0L ||
      anyNA(lambda_base) || any(!is.finite(lambda_base))) {
    stop("lambda_base must be a non-empty finite numeric vector.")
  }
  if (length(beta_vec) != length(lambda_base)) {
    stop("beta_vec and lambda_base must have the same length.")
  }
  if (any(lambda_base <= 0)) {
    stop("lambda_base must contain only positive values.")
  }

  weights_norm <- .irtsimrel_reducer_weights(weights, length(theta_vec))
  bank <- .irtsimrel_item_bank(
    beta = as.numeric(beta_vec),
    lambda_base = as.numeric(lambda_base),
    guessing = guessing
  )

  list(
    c = as.numeric(c),
    theta_vec = as.numeric(theta_vec),
    beta_vec = bank$beta,
    lambda_base = bank$lambda_base,
    guessing = bank$guessing,
    bank = bank,
    weights = weights_norm,
    weights_supplied = !is.null(weights)
  )
}


#' Internal: Stable Log-Sum-Exp
#'
#' @noRd
.irtsimrel_reducer_log_sum_exp <- function(x) {
  if (anyNA(x)) {
    stop("log-domain reduction received a missing value.")
  }
  if (any(is.infinite(x) & x > 0)) {
    return(Inf)
  }
  finite_x <- x[is.finite(x)]
  if (length(finite_x) == 0L) {
    return(-Inf)
  }
  pivot <- max(finite_x)
  pivot + log(sum(exp(x - pivot)))
}


#' Internal: Stable Weighted Log Mean
#'
#' @noRd
.irtsimrel_reducer_log_mean_exp <- function(log_values, weights) {
  log_weights <- .irtsimrel_reducer_log_weights(weights)
  positive <- is.finite(log_weights)
  .irtsimrel_reducer_log_sum_exp(
    log_values[positive] + log_weights[positive]
  )
}


#' Internal: Reliability Details from Log Test Information
#'
#' @noRd
.irtsimrel_reliability_details <- function(log_test_info,
                                           theta_var,
                                           weights,
                                           theta_var_source,
                                           weights_supplied,
                                           status = NULL) {
  log_mean_information <- .irtsimrel_reducer_log_mean_exp(log_test_info, weights)
  log_msem <- .irtsimrel_reducer_log_mean_exp(-log_test_info, weights)
  log_theta_var <- log(theta_var)
  rho_logit <- c(
    rho_tilde = log_theta_var + log_mean_information,
    rho_bar = log_theta_var - log_msem
  )
  rho <- stats::plogis(rho_logit)

  if (is.null(status)) {
    positive_nodes <- is.finite(.irtsimrel_reducer_log_weights(weights))
    relevant <- log_test_info[positive_nodes]
    if (all(is.finite(relevant))) {
      status <- "ok"
    } else if (any(is.infinite(relevant) & relevant < 0)) {
      status <- "zero_information"
    } else {
      status <- "nonfinite_information"
    }
  }

  below_floor <- log_test_info < log(1e-10)
  list(
    rho = rho,
    rho_logit = rho_logit,
    theta_var = unname(theta_var),
    theta_var_source = theta_var_source,
    log_mean_information = log_mean_information,
    log_msem = log_msem,
    log_info_range = unname(range(log_test_info)),
    below_legacy_floor_count = sum(below_floor),
    below_legacy_floor_rate = mean(below_floor),
    status = status,
    theta_measure = "empirical",
    approximation = if (weights_supplied) "weighted_nodes" else "finite_sample"
  )
}


#' Internal: Convert Test Information to Reliability
#'
#' @noRd
.compute_reliability_from_test_info <- function(test_info,
                                                theta_var,
                                                metric_internal = c("msem", "info"),
                                                weights = NULL) {
  metric_internal <- match.arg(metric_internal)
  if (!is.numeric(test_info) || length(test_info) == 0L || anyNA(test_info) ||
      any(test_info < 0)) {
    stop("test_info must be a non-negative numeric vector without missing values.")
  }
  weights_norm <- .irtsimrel_reducer_weights(weights, length(test_info))
  log_test_info <- log(test_info)
  theta_var <- .resolve_theta_var(c(-1, 1), theta_var)
  details <- .irtsimrel_reliability_details(
    log_test_info = log_test_info,
    theta_var = theta_var,
    weights = weights_norm,
    theta_var_source = attr(theta_var, "source"),
    weights_supplied = !is.null(weights)
  )
  if (metric_internal == "info") details$rho[["rho_tilde"]] else details$rho[["rho_bar"]]
}


#' Internal Engine for Reliability Computation
#'
#' @description
#' Internal function used by both `compute_rho_bar()` and `compute_rho_tilde()`.
#' Not exported.
#'
#' @param c Numeric scalar. Global discrimination scaling factor.
#' @param theta_vec Numeric vector of abilities.
#' @param beta_vec Numeric vector of item difficulties.
#' @param lambda_base Numeric vector of baseline discriminations.
#' @param theta_var Optional numeric. Pre-calculated variance of theta.
#' @param metric_internal Character. Either "msem" or "info".
#'
#' @return Numeric scalar reliability value.
#'
#' @noRd
.compute_rho_generic <- function(c,
                                 theta_vec,
                                 beta_vec,
                                 lambda_base,
                                 theta_var = NULL,
                                 metric_internal = c("msem", "info"),
                                 eta_base = NULL,
                                 guessing = NULL,
                                 weights = NULL,
                                 return_diagnostics = FALSE,
                                 chunk_size = NULL) {

  metric_internal <- match.arg(metric_internal)

  validated <- .validate_reliability_inputs(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    guessing = guessing,
    weights = weights
  )
  c <- validated$c
  theta_vec <- validated$theta_vec
  beta_vec <- validated$beta_vec
  lambda_base <- validated$lambda_base
  guessing <- validated$guessing
  weights_norm <- validated$weights

  if (!is.logical(return_diagnostics) || length(return_diagnostics) != 1L ||
      is.na(return_diagnostics)) {
    stop("return_diagnostics must be TRUE or FALSE.")
  }
  if (!is.null(eta_base)) {
    eta_base <- as.matrix(eta_base)
    if (!identical(dim(eta_base), c(length(theta_vec), length(beta_vec)))) {
      stop("eta_base must have dimensions length(theta_vec) by length(beta_vec).")
    }
  }

  theta_var <- .resolve_theta_var(
    theta_vec = theta_vec,
    theta_var = theta_var,
    weights = if (validated$weights_supplied) weights_norm else NULL
  )
  theta_var_source <- attr(theta_var, "source")

  if (c == 0) {
    log_test_info <- rep.int(-Inf, length(theta_vec))
    status <- "zero_scale"
  } else {
    log_test_info <- .irtsimrel_log_test_information(
      theta = theta_vec,
      bank = validated$bank,
      scale = c,
      chunk_size = chunk_size
    )
    status <- NULL
  }

  details <- .irtsimrel_reliability_details(
    log_test_info = log_test_info,
    theta_var = theta_var,
    weights = weights_norm,
    theta_var_source = theta_var_source,
    weights_supplied = validated$weights_supplied,
    status = status
  )
  metric_name <- if (metric_internal == "info") "rho_tilde" else "rho_bar"

  if (!return_diagnostics) {
    return(unname(details$rho[[metric_name]]))
  }

  details$metric <- metric_internal
  details$rho <- unname(details$rho[[metric_name]])
  details$rho_logit <- unname(details$rho_logit[[metric_name]])
  details
}


#' Compute Both Reliability Metrics in a Single Pass
#'
#' @description
#' Computes both the average-information reliability (\eqn{\tilde{\rho}}) and
#' the MSEM-based marginal reliability (\eqn{\bar{w}}) from a single set of
#' test-information reductions, avoiding redundant kernel evaluation.
#'
#' This is a performance optimization over calling \code{compute_rho_tilde()}
#' and \code{compute_rho_bar()} separately, since both share the same stable,
#' automatically chunked item-information reduction.
#'
#' @param c Numeric scalar. Global discrimination scaling factor.
#' @param theta_vec Numeric vector of abilities \eqn{\theta_m}.
#' @param beta_vec Numeric vector of item difficulties \eqn{\beta_i}.
#' @param lambda_base Numeric vector of baseline discriminations
#'   \eqn{\lambda_{i,0}} (before scaling by \code{c}).
#' @param theta_var Optional numeric. Pre-calculated variance of theta.
#'   If NULL, computed from \code{theta_vec}.
#' @param guessing Optional numeric scalar or item-length vector of 3PL lower
#'   asymptotes in `[0, 1)`. A scalar is recycled; `NULL` gives the Rasch/2PL
#'   special case. The kernel uses `D = 1`.
#' @param weights Optional non-negative integration weights, one per theta
#'   value. They are normalized internally.
#' @param return_diagnostics Logical. Return log-domain numerical diagnostics
#'   in addition to the two legacy reliability fields.
#'
#' @return A named list with components:
#' \describe{
#'   \item{\code{rho_tilde}}{Average-information reliability.}
#'   \item{\code{rho_bar}}{MSEM-based marginal reliability.}
#' }
#'
#' @details
#' For the same information values and latent-variance basis, Jensen's
#' inequality implies \eqn{\tilde{\rho} \geq \bar{w}}.
#' See \code{\link{compute_rho_bar}} for details on each metric.
#'
#' @examples
#' set.seed(1)
#' theta <- rnorm(1000)
#' beta  <- rnorm(20)
#' lambda0 <- rep(1, 20)
#'
#' both <- compute_rho_both(1, theta, beta, lambda0)
#' both$rho_tilde
#' both$rho_bar
#'
#' # Verify: rho_tilde >= rho_bar (Jensen's inequality)
#' both$rho_tilde >= both$rho_bar
#'
#' # The same call supports 3PL information.
#' compute_rho_both(1, theta, beta, lambda0, guessing = 0.20)
#'
#' @export
compute_rho_both <- function(c,
                             theta_vec,
                             beta_vec,
                             lambda_base,
                             theta_var = NULL,
                             guessing = NULL,
                             weights = NULL,
                             return_diagnostics = FALSE) {

  validated <- .validate_reliability_inputs(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    guessing = guessing,
    weights = weights
  )
  c <- validated$c
  theta_vec <- validated$theta_vec
  weights_norm <- validated$weights
  if (!is.logical(return_diagnostics) || length(return_diagnostics) != 1L ||
      is.na(return_diagnostics)) {
    stop("return_diagnostics must be TRUE or FALSE.")
  }

  theta_var <- .resolve_theta_var(
    theta_vec = theta_vec,
    theta_var = theta_var,
    weights = if (validated$weights_supplied) weights_norm else NULL
  )
  theta_var_source <- attr(theta_var, "source")

  if (c == 0) {
    log_test_info <- rep.int(-Inf, length(theta_vec))
    status <- "zero_scale"
  } else {
    log_test_info <- .irtsimrel_log_test_information(
      theta = theta_vec,
      bank = validated$bank,
      scale = c
    )
    status <- NULL
  }

  details <- .irtsimrel_reliability_details(
    log_test_info = log_test_info,
    theta_var = theta_var,
    weights = weights_norm,
    theta_var_source = theta_var_source,
    weights_supplied = validated$weights_supplied,
    status = status
  )
  legacy <- list(
    rho_tilde = unname(details$rho[["rho_tilde"]]),
    rho_bar = unname(details$rho[["rho_bar"]])
  )

  if (!return_diagnostics) {
    return(legacy)
  }

  details$rho <- unname(details$rho)
  details$rho_logit <- unname(details$rho_logit)
  c(legacy, details)
}


# =============================================================================
# Analytic Pre-Calibration (APC) Initialization
# =============================================================================

#' Analytic Pre-Calibration (APC) Initialization
#'
#' @description
#' Computes an initial value for the scaling factor using the closed-form
#' approximation under Gaussian Rasch assumptions.
#'
#' @param target_rho Numeric. Target reliability.
#' @param n_items Integer. Number of items.
#' @param sigma_beta Numeric. SD of item difficulties (default: 1.0).
#'
#' @return Numeric. Initial scaling factor c_init.
#'
#' @details
#' Under the Gaussian Rasch setting with \eqn{\theta \sim N(0,1)} and
#' \eqn{\beta \sim N(0, \sigma_\beta^2)}, the expected item information
#' involves the logistic-normal convolution:
#' \deqn{\kappa(\sigma^2) = \int \frac{e^z}{(1+e^z)^2} \phi(z; 0, \sigma^2) dz}
#'
#' Approximating \eqn{\kappa \approx 0.25 / \sqrt{1 + \sigma^2 \pi^2/3}},
#' the closed-form pre-calibration is:
#' \deqn{c_{init} = \sqrt{\frac{\rho^*}{I \cdot \kappa \cdot (1 - \rho^*)}}}
#'
#' @examples
#' # Compute initial c for target reliability of 0.80 with 25 items
#' compute_apc_init(target_rho = 0.80, n_items = 25)
#'
#' # With different difficulty spread
#' compute_apc_init(target_rho = 0.75, n_items = 20, sigma_beta = 1.5)
#'
#' @export
compute_apc_init <- function(target_rho, n_items, sigma_beta = 1.0) {

  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      is.na(target_rho) || !is.finite(target_rho) ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a scalar in (0, 1).")
  }

  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")

  if (!is.numeric(sigma_beta) || length(sigma_beta) != 1L ||
      is.na(sigma_beta) || !is.finite(sigma_beta) || sigma_beta < 0) {
    stop("`sigma_beta` must be a non-negative finite scalar.")
  }

  # Logistic-normal convolution approximation
  # kappa(sigma^2) for sigma^2 = 1 + sigma_beta^2
  sigma_sq <- 1 + sigma_beta^2
  kappa <- 0.25 / sqrt(1 + sigma_sq * pi^2 / 3)

  # Closed-form inversion
  c_init <- sqrt(target_rho / (n_items * kappa * (1 - target_rho)))

  # Bound to reasonable range
  c_init <- max(0.1, min(10, c_init))

  c_init
}


# =============================================================================
# Feasibility Screening
# =============================================================================

.classify_reliability_target <- function(target_rho,
                                         rho_range,
                                         boundary_tol = sqrt(.Machine$double.eps)) {
  if (is.null(target_rho)) {
    return(NULL)
  }

  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a scalar in (0, 1).")
  }

  if (!is.numeric(rho_range) || length(rho_range) < 2L) {
    stop("`rho_range` must be a numeric vector with lower and upper bounds.")
  }

  if (abs(target_rho - rho_range[1]) <= boundary_tol) {
    "boundary_lower"
  } else if (abs(target_rho - rho_range[2]) <= boundary_tol) {
    "boundary_upper"
  } else if (target_rho < rho_range[1]) {
    "below_lower"
  } else if (target_rho > rho_range[2]) {
    "above_upper"
  } else {
    "feasible"
  }
}

# Select the scale attached to the global maximum reported by a topology scan.
# The scanner's range is authoritative; re-evaluating its candidate locations
# through the caller's memoized objective keeps this adapter independent of the
# scanner's internal extrema-table representation.
.irtsimrel_topology_max_point <- function(topology, fn) {
  candidate_c <- topology$c_bounds

  extract_scale <- function(x) {
    if (is.null(x) || !(is.data.frame(x) || is.list(x))) {
      return(numeric())
    }
    nms <- names(x)
    scale_name <- intersect(
      c("c", "scale", "c_value", "location", "root"), nms
    )
    if (length(scale_name) == 0L) {
      return(numeric())
    }
    value <- x[[scale_name[1L]]]
    if (is.numeric(value)) value else numeric()
  }

  candidate_c <- unique(c(
    candidate_c,
    extract_scale(topology$scan),
    extract_scale(topology$extrema)
  ))
  candidate_c <- candidate_c[
    is.finite(candidate_c) & candidate_c >= topology$c_bounds[1L] &
      candidate_c <= topology$c_bounds[2L]
  ]

  candidate_rho <- vapply(candidate_c, fn, numeric(1))
  finite <- is.finite(candidate_rho)
  if (!any(finite)) {
    return(c(c = NA_real_, rho = NA_real_))
  }
  candidate_c <- candidate_c[finite]
  candidate_rho <- candidate_rho[finite]
  index <- which.max(candidate_rho)
  c(c = candidate_c[index], rho = candidate_rho[index])
}

.irtsimrel_topology_is_informative <- function(topology) {
  has_interior_extremum <- if (is.data.frame(topology$extrema) &&
                              "c" %in% names(topology$extrema)) {
    log_extrema <- log(topology$extrema$c)
    log_bounds <- log(topology$c_bounds)
    any(
      is.finite(log_extrema) &
        log_extrema > log_bounds[1L] + 1e-10 &
        log_extrema < log_bounds[2L] - 1e-10
    )
  } else {
    length(topology$extrema) > 0L
  }
  isTRUE(has_interior_extremum) ||
    !is.null(topology$topology_status) &&
      !identical(topology$topology_status, "resolved")
}

#' Check Feasibility of Target Reliability
#'
#' @description
#' Screens whether a given target reliability is attainable in the generated
#' finite Monte Carlo design for a particular test design (number of items,
#' model, latent distribution, item source) by computing the range of
#' empirical reliabilities across a range of scaling factors.
#'
#' This function is useful for determining whether a planned simulation
#' study is feasible before running the (potentially expensive) calibration
#' algorithms.
#'
#' @param n_items Integer. Number of items in the test form.
#' @param model Character. Measurement model: \code{"rasch"}, \code{"2pl"},
#'   or \code{"3pl"}. For 3PL, supply lower-asymptote generation controls in
#'   `item_params`; discrimination alone is varied under `D = 1`.
#' @param latent_shape Character. Shape argument passed to \code{sim_latentG()}.
#' @param item_source Character. Source argument passed to \code{sim_item_params()}.
#' @param c_bounds Numeric length-2 vector. Range of scaling factors to evaluate.
#'   Default: \code{c(0.1, 10)}.
#' @param M Integer. Monte Carlo sample size for theta. Default: 10000.
#' @param seed Optional integer for reproducibility.
#' @param latent_params List. Additional arguments passed to \code{sim_latentG()}.
#' @param item_params List. Additional arguments passed to \code{sim_item_params()}.
#' @param target_rho Optional numeric scalar in (0, 1). If supplied, the result
#'   includes status fields indicating whether the target is below, inside, or
#'   above the achievable range for each metric.
#' @param verbose Logical. If TRUE, print results.
#'
#' @return An object of class \code{"feasibility_check"} (a list) with:
#' \describe{
#'   \item{\code{rho_range_info}}{Numeric length-2 vector: empirical range of
#'     average-information reliability (\eqn{\tilde{\rho}}) computed by the
#'     configured finite scan; provisional when topology is unresolved.}
#'   \item{\code{rho_range_msem}}{Numeric length-2 vector: empirical range of
#'     MSEM-based reliability (\eqn{\bar{w}}) computed by the configured finite
#'     scan; provisional when topology is unresolved.}
#'   \item{\code{rho_bounds_info}}{Endpoint reliabilities for
#'     average-information reliability.}
#'   \item{\code{rho_bounds_msem}}{Endpoint reliabilities for MSEM-based
#'     reliability.}
#'   \item{\code{rho_info_max_c}}{Scaling factor at the maximum information
#'     reliability detected by the configured scan within \code{c_bounds};
#'     provisional when topology is unresolved.}
#'   \item{\code{rho_msem_max_c}}{Scaling factor at the maximum MSEM
#'     reliability detected by the configured scan within \code{c_bounds};
#'     provisional when topology is unresolved.}
#'   \item{\code{target_status_info}}{If \code{target_rho} is supplied, one of
#'     \code{"below_lower"}, \code{"boundary_lower"}, \code{"feasible"},
#'     \code{"boundary_upper"}, or \code{"above_upper"} for \eqn{\tilde{\rho}}.}
#'   \item{\code{target_status_msem}}{Analogous status for \eqn{\bar{w}}.}
#'   \item{\code{topology_info}, \code{topology_msem}}{Complete configured
#'     log-grid scans, polished extrema, detected roots, monotone branches,
#'     evaluation counts, and resolution status for each metric.}
#'   \item{\code{target_status_info_canonical},
#'     \code{target_status_msem_canonical}}{Canonical feasibility status:
#'     \code{"feasible"}, \code{"infeasible_below_range"},
#'     \code{"infeasible_above_range"}, or \code{"uncertain"}.}
#'   \item{\code{root_count_info}, \code{root_count_msem},
#'     \code{admissible_root_count_info},
#'     \code{admissible_root_count_msem}}{Detected roots and roots on
#'     increasing crossing/boundary branches.}
#'   \item{\code{best_achievable_info},
#'     \code{best_achievable_msem}}{Closest detected point to
#'     \code{target_rho}, including scale, reliability, residual, absolute
#'     error, and location type.}
#'   \item{\code{n_items}}{Number of items.}
#'   \item{\code{model}}{Model used.}
#'   \item{\code{latent_shape}}{Latent distribution shape.}
#'   \item{\code{c_bounds}}{Scaling factor bounds evaluated.}
#'   \item{\code{M}}{Monte Carlo sample size.}
#'   \item{\code{theta_var}}{Estimated latent variance.}
#' }
#'
#' @details
#' Either reliability metric can be non-monotone over a sufficiently broad
#' finite empirical integration interval. Both metrics therefore use the same
#' adaptive log-scale topology scan as the calibrators. Endpoint values remain
#' available in \code{rho_bounds_*}, while \code{rho_range_*} includes all
#' resolved interior extrema within the user-supplied bounds.
#'
#' These ranges and target classifications are conditional on the generated
#' finite theta sample. In particular, for the built-in
#' \code{latent_shape = "heavy_tail"}, the population MSEM functional can be
#' non-integrable even though a finite Monte Carlo sample yields a numeric
#' \code{rho_range_msem}. Treat that MSEM range as finite-grid sensitivity
#' evidence, not as a population feasibility guarantee.
#'
#' @examples
#' # Check feasibility for 25-item Rasch test
#' feas <- check_feasibility(n_items = 25, model = "rasch",
#'                           target_rho = 0.90, seed = 42,
#'                           M = 5000, verbose = FALSE)
#' print(feas)
#'
#' # Metric-specific target status
#' feas$target_status_info
#' feas$target_status_msem
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{rho_curve}}
#'
#' @export
check_feasibility <- function(n_items,
                              model = c("rasch", "2pl", "3pl"),
                              latent_shape = "normal",
                              item_source = "parametric",
                              c_bounds = c(0.1, 10),
                              M = 10000L,
                              seed = NULL,
                              latent_params = list(),
                              item_params = list(),
                              target_rho = NULL,
                              verbose = TRUE) {

  model <- match.arg(model)
  latent_shape <- .irtsimrel_match_latent_shape(latent_shape)

  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(!is.finite(c_bounds)) || any(c_bounds <= 0) ||
      c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  M <- .irtsimrel_validate_positive_integer_scalar(M, "M")
  if (M < 2L) {
    stop("`M` must be at least 2 to estimate latent variance.")
  }

  if (!is.null(target_rho) &&
      (!is.numeric(target_rho) || length(target_rho) != 1L ||
       is.na(target_rho) || !is.finite(target_rho) ||
       target_rho <= 0 || target_rho >= 1)) {
    stop("`target_rho` must be a scalar in (0, 1).")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)
  item_params <- .irtsimrel_normalize_item_params(item_params)

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  # Generate theta
  latent_args <- modifyList(
    list(n = M, shape = latent_shape),
    latent_params
  )
  theta_vec <- do.call(sim_latentG, latent_args)$theta
  theta_var <- var(theta_vec)
  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 1e-10) {
    stop("Feasibility latent variance must be finite and positive.")
  }

  # Generate items
  item_args <- modifyList(
    list(n_items = n_items, model = model, source = item_source, scale = 1),
    item_params
  )
  items <- do.call(sim_item_params, item_args)
  beta_vec <- items$data$beta
  guessing_vec <- if ("guessing" %in% names(items$data)) {
    items$data$guessing
  } else {
    rep.int(0, n_items)
  }

  if (model == "rasch") {
    lambda_base <- rep(1, n_items)
  } else {
    if ("lambda_unscaled" %in% names(items$data)) {
      lambda_base <- items$data$lambda_unscaled
    } else {
      lambda_base <- items$data$lambda
    }
  }

  # Both metrics share one memoized information evaluation. Each adaptive scan
  # nevertheless refines its own extrema and target crossings on log(c).
  rho_cache <- new.env(hash = TRUE, parent = emptyenv())
  rho_both_fn <- function(c_val) {
    key <- sprintf("%.17g", c_val)
    if (!exists(key, envir = rho_cache, inherits = FALSE)) {
      assign(
        key,
        compute_rho_both(
          c_val, theta_vec, beta_vec, lambda_base,
          theta_var = theta_var, guessing = guessing_vec
        ),
        envir = rho_cache
      )
    }
    get(key, envir = rho_cache, inherits = FALSE)
  }
  rho_tilde_fn <- function(c_val) rho_both_fn(c_val)$rho_tilde
  rho_bar_fn <- function(c_val) rho_both_fn(c_val)$rho_bar

  topology_info <- .irtsimrel_scan_topology(
    fn = rho_tilde_fn,
    c_bounds = c_bounds,
    target = target_rho
  )
  topology_msem <- .irtsimrel_scan_topology(
    fn = rho_bar_fn,
    c_bounds = c_bounds,
    target = target_rho
  )

  rho_bounds_info <- topology_info$rho_bounds
  rho_bounds_msem <- topology_msem$rho_bounds
  rho_range_info <- topology_info$rho_range
  rho_range_msem <- topology_msem$rho_range
  info_max <- .irtsimrel_topology_max_point(topology_info, rho_tilde_fn)
  msem_max <- .irtsimrel_topology_max_point(topology_msem, rho_bar_fn)

  result <- list(
    rho_range_info = rho_range_info,
    rho_range_msem = rho_range_msem,
    rho_bounds_info = rho_bounds_info,
    rho_bounds_msem = rho_bounds_msem,
    rho_info_max_c = unname(info_max[["c"]]),
    rho_msem_max_c = unname(msem_max[["c"]]),
    target_rho = target_rho,
    target_status_info = .classify_reliability_target(target_rho, rho_range_info),
    target_status_msem = .classify_reliability_target(target_rho, rho_range_msem),
    n_items = n_items,
    model = model,
    latent_shape = latent_shape,
    c_bounds = c_bounds,
    M = M,
    theta_var = theta_var,
    topology_info = topology_info,
    topology_msem = topology_msem,
    target_status_info_canonical = topology_info$target_status,
    target_status_msem_canonical = topology_msem$target_status,
    root_count_info = topology_info$root_count,
    root_count_msem = topology_msem$root_count,
    admissible_root_count_info = topology_info$admissible_root_count,
    admissible_root_count_msem = topology_msem$admissible_root_count,
    best_achievable_info = topology_info$best_achievable,
    best_achievable_msem = topology_msem$best_achievable
  )

  class(result) <- c("feasibility_check", "list")

  if (verbose) {
    print(result)
  }

  invisible(result)
}


#' @rdname check_feasibility
#' @param x An object of class \code{"feasibility_check"}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.feasibility_check <- function(x, digits = 4, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("  Feasibility Check: Achievable Reliability Range\n")
  cat("=======================================================\n\n")

  cat(sprintf("  Number of items  : %d\n", x$n_items))
  cat(sprintf("  Model            : %s\n", toupper(x$model)))
  cat(sprintf("  Latent shape     : %s\n", x$latent_shape))
  cat(sprintf("  Latent variance  : %.*f\n", digits, x$theta_var))
  cat(sprintf("  c range          : [%.2f, %.2f]\n", x$c_bounds[1], x$c_bounds[2]))
  cat(sprintf("  Monte Carlo M    : %d\n", x$M))
  cat("\n")

  cat("Achievable Reliability Ranges:\n")
  cat(sprintf("  rho_tilde (info) : [%.*f, %.*f]\n",
              digits, x$rho_range_info[1], digits, x$rho_range_info[2]))
  cat(sprintf("  rho_bar   (msem) : [%.*f, %.*f]\n",
              digits, x$rho_range_msem[1], digits, x$rho_range_msem[2]))
  cat("\n")

  if (!is.null(x$target_rho)) {
    cat(sprintf("Target rho*        : %.*f\n", digits, x$target_rho))
    cat(sprintf("  info status      : %s\n", x$target_status_info))
    cat(sprintf("  msem status      : %s\n", x$target_status_msem))
    cat("\n")
  }

  cat("Note: rho_tilde >= rho_bar on the same information grid (Jensen's inequality).\n")
  cat("  rho_tilde range screens EQC targets; root policy must still admit a root.\n")
  cat("  rho_bar range screens SAC targets; stable interior-branch preflight is also required.\n")
  if (identical(x$latent_shape, "heavy_tail")) {
    cat("  heavy_tail MSEM is finite-grid sensitivity only; population targeting is not certified.\n")
  }
  cat("\n")

  invisible(x)
}


# =============================================================================
# Reliability Curve
# =============================================================================

#' Compute Reliability as a Function of Scaling Factor
#'
#' @description
#' Computes and optionally plots the reliability curve \eqn{\rho(c)} across
#' a grid of scaling factor values. This visualization helps understand how
#' reliability varies with the discrimination scaling factor and aids in
#' selecting appropriate target reliability values.
#'
#' @param c_values Numeric vector. Grid of scaling factor values to evaluate.
#'   Default: \code{seq(0.1, 5, length.out = 50)}.
#' @param n_items Integer. Number of items in the test form.
#' @param model Character. Measurement model: \code{"rasch"}, \code{"2pl"},
#'   or \code{"3pl"}. For 3PL, supply lower-asymptote generation controls in
#'   `item_params`; discrimination alone is varied under `D = 1`.
#' @param latent_shape Character. Shape argument passed to \code{sim_latentG()}.
#' @param item_source Character. Source argument passed to \code{sim_item_params()}.
#' @param metric Character. Which reliability metric(s) to compute:
#'   \code{"both"}, \code{"info"}, or \code{"msem"}.
#' @param M Integer. Monte Carlo sample size. Default: 5000.
#' @param seed Optional integer for reproducibility.
#' @param latent_params List. Additional arguments passed to \code{sim_latentG()}.
#' @param item_params List. Additional arguments passed to \code{sim_item_params()}.
#' @param plot Logical. If TRUE (default), create a plot of the curve.
#'
#' @return A data frame of class \code{"rho_curve"} with columns:
#' \describe{
#'   \item{\code{c}}{Scaling factor values.}
#'   \item{\code{rho_tilde}}{Average-information reliability (if metric includes "info").}
#'   \item{\code{rho_bar}}{MSEM-based reliability (if metric includes "msem").}
#' }
#'
#' @details
#' The function generates a single set of theta and item parameters, then
#' evaluates the reliability at each value of \code{c_values}. When
#' \code{metric = "both"}, it uses \code{\link{compute_rho_both}} for
#' efficiency. For a curve with detected interior extrema or unresolved
#' regions, an independently scanned \code{topology_info} and/or
#' \code{topology_msem} attribute is attached. Topology is never inferred from
#' the display grid alone.
#'
#' Values are conditional on the generated finite theta sample. For the
#' built-in \code{latent_shape = "heavy_tail"}, a numeric MSEM curve is a
#' finite-grid sensitivity result; the corresponding population MSEM
#' functional can be non-integrable and is not certified by this function.
#'
#' @examples
#' # Basic usage: compute reliability curve for 25-item Rasch test
#' curve_data <- rho_curve(n_items = 25, model = "rasch", seed = 42,
#'                         M = 3000, plot = FALSE)
#' head(curve_data)
#'
#' @seealso
#' \code{\link{check_feasibility}}, \code{\link{compute_rho_both}}
#'
#' @export
rho_curve <- function(c_values = seq(0.1, 5, length.out = 50),
                      n_items,
                      model = c("rasch", "2pl", "3pl"),
                      latent_shape = "normal",
                      item_source = "parametric",
                      metric = c("both", "info", "msem"),
                      M = 5000L,
                      seed = NULL,
                      latent_params = list(),
                      item_params = list(),
                      plot = TRUE) {

  model  <- match.arg(model)
  metric <- match.arg(metric)
  latent_shape <- .irtsimrel_match_latent_shape(latent_shape)

  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")

  if (!is.numeric(c_values) || length(c_values) < 2L ||
      any(!is.finite(c_values)) || any(c_values <= 0)) {
    stop("`c_values` must be a finite numeric vector of positive values with length >= 2.")
  }

  M <- .irtsimrel_validate_positive_integer_scalar(M, "M")
  if (M < 2L) {
    stop("`M` must be at least 2 to estimate latent variance.")
  }

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)
  item_params <- .irtsimrel_normalize_item_params(item_params)

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  # Generate theta
  latent_args <- modifyList(
    list(n = M, shape = latent_shape),
    latent_params
  )
  theta_vec <- do.call(sim_latentG, latent_args)$theta
  theta_var <- var(theta_vec)

  # Generate items
  item_args <- modifyList(
    list(n_items = n_items, model = model, source = item_source, scale = 1),
    item_params
  )
  items <- do.call(sim_item_params, item_args)
  beta_vec <- items$data$beta
  guessing_vec <- if ("guessing" %in% names(items$data)) {
    items$data$guessing
  } else {
    rep.int(0, n_items)
  }

  if (model == "rasch") {
    lambda_base <- rep(1, n_items)
  } else {
    if ("lambda_unscaled" %in% names(items$data)) {
      lambda_base <- items$data$lambda_unscaled
    } else {
      lambda_base <- items$data$lambda
    }
  }

  # Compute reliability at each c value
  eta_base <- outer(theta_vec, beta_vec, "-")

  if (metric == "both") {
    results <- lapply(c_values, function(c_val) {
      compute_rho_both(
        c = c_val,
        theta_vec = theta_vec,
        beta_vec = beta_vec,
        lambda_base = lambda_base,
        theta_var = theta_var,
        guessing = guessing_vec
      )
    })
    df <- data.frame(
      c = c_values,
      rho_tilde = vapply(results, `[[`, numeric(1), "rho_tilde"),
      rho_bar   = vapply(results, `[[`, numeric(1), "rho_bar")
    )
  } else if (metric == "info") {
    rho_vals <- vapply(c_values, function(c_val) {
      .compute_rho_generic(
        c_val, theta_vec, beta_vec, lambda_base,
        theta_var = theta_var,
        metric_internal = "info",
        eta_base = eta_base,
        guessing = guessing_vec
      )
    }, numeric(1))
    df <- data.frame(c = c_values, rho_tilde = rho_vals)
  } else {
    rho_vals <- vapply(c_values, function(c_val) {
      .compute_rho_generic(
        c_val, theta_vec, beta_vec, lambda_base,
        theta_var = theta_var,
        metric_internal = "msem",
        eta_base = eta_base,
        guessing = guessing_vec
      )
    }, numeric(1))
    df <- data.frame(c = c_values, rho_bar = rho_vals)
  }

  class(df) <- c("rho_curve", "data.frame")
  attr(df, "metric") <- metric
  attr(df, "n_items") <- n_items
  attr(df, "model") <- model
  attr(df, "latent_shape") <- latent_shape
  attr(df, "theta_var") <- theta_var

  # Topology metadata, when present, comes from the same adaptive log-scale
  # scanner used by calibration and feasibility checks. The displayed points
  # are intentionally never used to infer extrema or branch structure.
  topology_cache <- new.env(hash = TRUE, parent = emptyenv())
  topology_both_fn <- function(c_val) {
    key <- sprintf("%.17g", c_val)
    if (!exists(key, envir = topology_cache, inherits = FALSE)) {
      assign(
        key,
        compute_rho_both(
          c_val, theta_vec, beta_vec, lambda_base,
          theta_var = theta_var, guessing = guessing_vec
        ),
        envir = topology_cache
      )
    }
    get(key, envir = topology_cache, inherits = FALSE)
  }
  curve_bounds <- range(c_values)
  scan_curve_topology <- curve_bounds[1L] < curve_bounds[2L]
  if (scan_curve_topology && metric %in% c("both", "info")) {
    topology_info <- .irtsimrel_scan_topology(
      fn = function(c_val) topology_both_fn(c_val)$rho_tilde,
      c_bounds = curve_bounds
    )
    if (.irtsimrel_topology_is_informative(topology_info)) {
      attr(df, "topology_info") <- topology_info
    }
  }
  if (scan_curve_topology && metric %in% c("both", "msem")) {
    topology_msem <- .irtsimrel_scan_topology(
      fn = function(c_val) topology_both_fn(c_val)$rho_bar,
      c_bounds = curve_bounds
    )
    if (.irtsimrel_topology_is_informative(topology_msem)) {
      attr(df, "topology_msem") <- topology_msem
    }
  }

  # Plot if requested
  if (plot) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      .plot_rho_curve_ggplot(df, metric, n_items, model, latent_shape)
    } else {
      .plot_rho_curve_base(df, metric, n_items, model, latent_shape)
    }
  }

  invisible(df)
}


#' @rdname rho_curve
#' @param x An object of class \code{"rho_curve"}.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.rho_curve <- function(x, ...) {
  metric <- attr(x, "metric")
  n_items <- attr(x, "n_items")
  model <- attr(x, "model")

  cat("Reliability Curve\n")
  cat("=================\n")
  cat(sprintf("  Items: %d | Model: %s | Metric: %s\n",
              n_items, toupper(model), metric))
  cat(sprintf("  c range: [%.2f, %.2f] (%d points)\n",
              min(x$c), max(x$c), nrow(x)))

  if ("rho_tilde" %in% names(x)) {
    cat(sprintf("  rho_tilde range: [%.4f, %.4f]\n",
                min(x$rho_tilde), max(x$rho_tilde)))
  }
  if ("rho_bar" %in% names(x)) {
    cat(sprintf("  rho_bar range  : [%.4f, %.4f]\n",
                min(x$rho_bar), max(x$rho_bar)))
  }
  cat("\n")

  # Print first few rows
  print(utils::head(as.data.frame(x), 6))
  if (nrow(x) > 6) {
    cat(sprintf("  ... (%d more rows)\n", nrow(x) - 6))
  }

  invisible(x)
}


#' Internal: Plot rho_curve with ggplot2
#' @noRd
.plot_rho_curve_ggplot <- function(df, metric, n_items, model, latent_shape) {
  # Build long-format data for ggplot
  if (metric == "both") {
    df_long <- data.frame(
      c = rep(df$c, 2),
      rho = c(df$rho_tilde, df$rho_bar),
      metric = rep(c("rho_tilde (info)", "rho_bar (msem)"), each = nrow(df))
    )
  } else if (metric == "info") {
    df_long <- data.frame(c = df$c, rho = df$rho_tilde,
                          metric = "rho_tilde (info)")
  } else {
    df_long <- data.frame(c = df$c, rho = df$rho_bar,
                          metric = "rho_bar (msem)")
  }

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = c, y = rho, color = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = sprintf("Reliability Curve: %d-item %s Test",
                       n_items, toupper(model)),
      subtitle = sprintf("Latent shape: %s", latent_shape),
      x = "Scaling Factor c",
      y = "Reliability",
      color = "Metric"
    ) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}


#' Internal: Plot rho_curve with base R
#' @noRd
.plot_rho_curve_base <- function(df, metric, n_items, model, latent_shape) {
  oldpar <- par(mar = c(5, 4, 4, 2) + 0.1)
  on.exit(par(oldpar))

  if (metric == "both") {
    plot(df$c, df$rho_tilde, type = "l", col = "steelblue", lwd = 2,
         xlab = "Scaling Factor c", ylab = "Reliability",
         main = sprintf("Reliability Curve: %d-item %s", n_items, toupper(model)),
         ylim = c(0, 1))
    lines(df$c, df$rho_bar, col = "coral", lwd = 2, lty = 2)
    legend("bottomright",
           legend = c("rho_tilde (info)", "rho_bar (msem)"),
           col = c("steelblue", "coral"), lty = c(1, 2), lwd = 2)
  } else if (metric == "info") {
    plot(df$c, df$rho_tilde, type = "l", col = "steelblue", lwd = 2,
         xlab = "Scaling Factor c", ylab = "Reliability (rho_tilde)",
         main = sprintf("Reliability Curve: %d-item %s", n_items, toupper(model)),
         ylim = c(0, 1))
  } else {
    plot(df$c, df$rho_bar, type = "l", col = "coral", lwd = 2,
         xlab = "Scaling Factor c", ylab = "Reliability (rho_bar)",
         main = sprintf("Reliability Curve: %d-item %s", n_items, toupper(model)),
         ylim = c(0, 1))
  }
}
