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
#'
#' @return Numeric scalar reliability value.
#'
#' @details
#' The computation proceeds as:
#' \enumerate{
#'   \item Form scaled discriminations \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}.
#'   \item Compute item response probabilities
#'     \eqn{p_{mi} = \text{logit}^{-1}\{\lambda_i(c)(\theta_m - \beta_i)\}}.
#'   \item Item information:
#'     \eqn{\mathcal{J}_{mi} = \lambda_i(c)^2 p_{mi}(1-p_{mi})}.
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
#' A small floor (\code{1e-10}) is applied to test information to avoid
#' numerical problems when taking reciprocals.
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
#' @export
compute_rho_bar <- function(c, theta_vec, beta_vec, lambda_base, theta_var = NULL) {
  .compute_rho_generic(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    theta_var = theta_var,
    metric_internal = "msem"
  )
}


#' @rdname compute_rho_bar
#' @export
compute_rho_tilde <- function(c, theta_vec, beta_vec, lambda_base, theta_var = NULL) {
  .compute_rho_generic(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base,
    theta_var = theta_var,
    metric_internal = "info"
  )
}


#' Internal: Resolve Theta Variance
#'
#' @noRd
.resolve_theta_var <- function(theta_vec, theta_var = NULL) {
  if (is.null(theta_var)) {
    theta_var <- var(theta_vec)
  }

  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 1e-10) {
    warning("Invalid theta_var. Falling back to 1.0.")
    theta_var <- 1.0
  }

  theta_var
}


#' Internal: Compute Test Information
#'
#' @noRd
.compute_test_information <- function(c,
                                      theta_vec,
                                      beta_vec,
                                      lambda_base,
                                      eta_base = NULL) {
  lambda_curr <- c * lambda_base

  if (is.null(eta_base)) {
    eta_base <- outer(theta_vec, beta_vec, "-")
  } else {
    eta_base <- as.matrix(eta_base)
    if (!identical(dim(eta_base), c(length(theta_vec), length(beta_vec)))) {
      stop("eta_base must have dimensions length(theta_vec) by length(beta_vec).")
    }
  }

  eta_mat <- sweep(eta_base, 2, lambda_curr, FUN = "*")
  p_mat <- plogis(eta_mat)
  info_mat <- sweep(p_mat * (1 - p_mat), 2, lambda_curr^2, FUN = "*")
  rowSums(info_mat)
}

.validate_reliability_inputs <- function(c, theta_vec, beta_vec, lambda_base) {
  c <- as.numeric(c)
  theta_vec <- as.numeric(theta_vec)
  beta_vec <- as.numeric(beta_vec)
  lambda_base <- as.numeric(lambda_base)

  if (length(c) != 1L || is.na(c) || !is.finite(c)) {
    stop("c must be a finite numeric scalar.")
  }
  if (length(theta_vec) < 2L || any(is.na(theta_vec)) ||
      any(!is.finite(theta_vec))) {
    stop("theta_vec must be a finite numeric vector with length >= 2.")
  }
  if (length(beta_vec) == 0L || any(is.na(beta_vec)) ||
      any(!is.finite(beta_vec))) {
    stop("beta_vec must be a non-empty finite numeric vector.")
  }
  if (length(lambda_base) == 0L || any(is.na(lambda_base)) ||
      any(!is.finite(lambda_base))) {
    stop("lambda_base must be a non-empty finite numeric vector.")
  }
  if (length(beta_vec) != length(lambda_base)) {
    stop("beta_vec and lambda_base must have the same length.")
  }

  list(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base
  )
}


#' Internal: Convert Test Information to Reliability
#'
#' @noRd
.compute_reliability_from_test_info <- function(test_info,
                                                theta_var,
                                                metric_internal = c("msem", "info")) {
  metric_internal <- match.arg(metric_internal)

  if (metric_internal == "info") {
    J_bar <- mean(test_info)
    return((theta_var * J_bar) / (theta_var * J_bar + 1))
  }

  test_info_safe <- pmax(test_info, 1e-10)
  msem <- mean(1 / test_info_safe)
  theta_var / (theta_var + msem)
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
                                 eta_base = NULL) {

  metric_internal <- match.arg(metric_internal)

  validated <- .validate_reliability_inputs(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base
  )
  c <- validated$c
  theta_vec <- validated$theta_vec
  beta_vec <- validated$beta_vec
  lambda_base <- validated$lambda_base

  if (c <= 0) {
    return(0)
  }

  theta_var <- .resolve_theta_var(theta_vec, theta_var)
  test_info <- .compute_test_information(c, theta_vec, beta_vec, lambda_base, eta_base)

  .compute_reliability_from_test_info(test_info, theta_var, metric_internal)
}


#' Compute Both Reliability Metrics in a Single Pass
#'
#' @description
#' Computes both the average-information reliability (\eqn{\tilde{\rho}}) and
#' the MSEM-based marginal reliability (\eqn{\bar{w}}) from a single set of
#' test information values, avoiding redundant matrix computation.
#'
#' This is a performance optimization over calling \code{compute_rho_tilde()}
#' and \code{compute_rho_bar()} separately, since both share the same
#' \eqn{M \times I} matrix computations.
#'
#' @param c Numeric scalar. Global discrimination scaling factor.
#' @param theta_vec Numeric vector of abilities \eqn{\theta_m}.
#' @param beta_vec Numeric vector of item difficulties \eqn{\beta_i}.
#' @param lambda_base Numeric vector of baseline discriminations
#'   \eqn{\lambda_{i,0}} (before scaling by \code{c}).
#' @param theta_var Optional numeric. Pre-calculated variance of theta.
#'   If NULL, computed from \code{theta_vec}.
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
#' @export
compute_rho_both <- function(c, theta_vec, beta_vec, lambda_base, theta_var = NULL) {

  validated <- .validate_reliability_inputs(
    c = c,
    theta_vec = theta_vec,
    beta_vec = beta_vec,
    lambda_base = lambda_base
  )
  c <- validated$c
  theta_vec <- validated$theta_vec
  beta_vec <- validated$beta_vec
  lambda_base <- validated$lambda_base

  if (c <= 0) {
    return(list(rho_tilde = 0, rho_bar = 0))
  }

  theta_var <- .resolve_theta_var(theta_vec, theta_var)
  test_info <- .compute_test_information(c, theta_vec, beta_vec, lambda_base)

  rho_tilde <- .compute_reliability_from_test_info(test_info, theta_var, "info")
  rho_bar <- .compute_reliability_from_test_info(test_info, theta_var, "msem")

  list(rho_tilde = rho_tilde, rho_bar = rho_bar)
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

#' Check Feasibility of Target Reliability
#'
#' @description
#' Screens whether a given target reliability is achievable for a particular
#' test design (number of items, model, latent distribution, item source)
#' by computing the range of achievable reliabilities across a range of
#' scaling factors.
#'
#' This function is useful for determining whether a planned simulation
#' study is feasible before running the (potentially expensive) calibration
#' algorithms.
#'
#' @param n_items Integer. Number of items in the test form.
#' @param model Character. Measurement model: \code{"rasch"} or \code{"2pl"}.
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
#'   \item{\code{rho_range_info}}{Numeric length-2 vector: achievable range of
#'     average-information reliability (\eqn{\tilde{\rho}}).}
#'   \item{\code{rho_range_msem}}{Numeric length-2 vector: achievable range of
#'     MSEM-based reliability (\eqn{\bar{w}}).}
#'   \item{\code{rho_bounds_info}}{Endpoint reliabilities for
#'     average-information reliability.}
#'   \item{\code{rho_bounds_msem}}{Endpoint reliabilities for MSEM-based
#'     reliability.}
#'   \item{\code{rho_info_max_c}}{Scaling factor at the maximum sampled info
#'     reliability within \code{c_bounds}.}
#'   \item{\code{rho_msem_max_c}}{Scaling factor at the maximum sampled MSEM
#'     reliability within \code{c_bounds}.}
#'   \item{\code{target_status_info}}{If \code{target_rho} is supplied, one of
#'     \code{"below_lower"}, \code{"boundary_lower"}, \code{"feasible"},
#'     \code{"boundary_upper"}, or \code{"above_upper"} for \eqn{\tilde{\rho}}.}
#'   \item{\code{target_status_msem}}{Analogous status for \eqn{\bar{w}}.}
#'   \item{\code{n_items}}{Number of items.}
#'   \item{\code{model}}{Model used.}
#'   \item{\code{latent_shape}}{Latent distribution shape.}
#'   \item{\code{c_bounds}}{Scaling factor bounds evaluated.}
#'   \item{\code{M}}{Monte Carlo sample size.}
#'   \item{\code{theta_var}}{Estimated latent variance.}
#' }
#'
#' @details
#' For the average-information metric (\eqn{\tilde{\rho}}), practical EQC
#' search intervals are usually monotone. With very wide finite empirical
#' quadrature bounds, however, logistic probabilities can saturate numerically.
#' The function therefore reports endpoint reliabilities separately and uses an
#' interior maximization pass for the feasible upper range.
#'
#' For the MSEM-based metric (\eqn{\bar{w}}), the reliability can be
#' non-monotone at extreme scaling factors. The function uses
#' \code{\link[stats]{optimize}} to find the maximum within the bounds, and
#' the range is \eqn{[\min(\bar{w}(c_{min}), \bar{w}(c_{max})), \bar{w}_{max}]}.
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
                              model = c("rasch", "2pl"),
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

  if (model == "rasch") {
    lambda_base <- rep(1, n_items)
  } else {
    if ("lambda_unscaled" %in% names(items$data)) {
      lambda_base <- items$data$lambda_unscaled
    } else {
      lambda_base <- items$data$lambda
    }
  }

  # Compute rho_tilde at bounds and estimate the interior maximum. Extremely
  # wide finite empirical quadrature intervals can saturate in the far tail, so
  # endpoint-only checks can understate the achievable upper reliability.
  rho_tilde_low  <- compute_rho_tilde(c_bounds[1], theta_vec, beta_vec,
                                       lambda_base, theta_var = theta_var)
  rho_tilde_high <- compute_rho_tilde(c_bounds[2], theta_vec, beta_vec,
                                       lambda_base, theta_var = theta_var)
  rho_tilde_fn <- function(c_val) {
    compute_rho_tilde(c_val, theta_vec, beta_vec, lambda_base,
                      theta_var = theta_var)
  }
  opt_info <- optimize(rho_tilde_fn, interval = c_bounds, maximum = TRUE)
  rho_tilde_max <- opt_info$objective
  rho_tilde_max_c <- opt_info$maximum
  if (rho_tilde_high >= rho_tilde_max - 1e-8) {
    rho_tilde_max <- rho_tilde_high
    rho_tilde_max_c <- c_bounds[2]
  }
  rho_bounds_info <- c(lower = rho_tilde_low, upper = rho_tilde_high)
  rho_range_info <- c(
    lower = min(rho_tilde_low, rho_tilde_high),
    upper = rho_tilde_max
  )

  # Compute rho_bar at bounds and find maximum (non-monotone!)
  rho_bar_low  <- compute_rho_bar(c_bounds[1], theta_vec, beta_vec,
                                   lambda_base, theta_var = theta_var)
  rho_bar_high <- compute_rho_bar(c_bounds[2], theta_vec, beta_vec,
                                   lambda_base, theta_var = theta_var)

  # Use optimize() to find maximum of rho_bar in c_bounds
  rho_bar_fn <- function(c_val) {
    compute_rho_bar(c_val, theta_vec, beta_vec, lambda_base, theta_var = theta_var)
  }
  opt_result <- optimize(rho_bar_fn, interval = c_bounds, maximum = TRUE)
  rho_bar_max <- opt_result$objective

  rho_bounds_msem <- c(lower = rho_bar_low, upper = rho_bar_high)
  rho_range_msem <- c(lower = min(rho_bar_low, rho_bar_high), upper = rho_bar_max)

  result <- list(
    rho_range_info = rho_range_info,
    rho_range_msem = rho_range_msem,
    rho_bounds_info = rho_bounds_info,
    rho_bounds_msem = rho_bounds_msem,
    rho_info_max_c = rho_tilde_max_c,
    rho_msem_max_c = opt_result$maximum,
    target_rho = target_rho,
    target_status_info = .classify_reliability_target(target_rho, rho_range_info),
    target_status_msem = .classify_reliability_target(target_rho, rho_range_msem),
    n_items = n_items,
    model = model,
    latent_shape = latent_shape,
    c_bounds = c_bounds,
    M = M,
    theta_var = theta_var
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
  cat("  Use rho_tilde range for EQC targets.\n")
  cat("  Use rho_bar range for SAC targets.\n")
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
#' @param model Character. Measurement model: \code{"rasch"} or \code{"2pl"}.
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
#' efficiency.
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
                      model = c("rasch", "2pl"),
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
      test_info <- .compute_test_information(
        c_val, theta_vec, beta_vec, lambda_base, eta_base = eta_base
      )
      list(
        rho_tilde = .compute_reliability_from_test_info(test_info, theta_var, "info"),
        rho_bar = .compute_reliability_from_test_info(test_info, theta_var, "msem")
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
        eta_base = eta_base
      )
    }, numeric(1))
    df <- data.frame(c = c_values, rho_tilde = rho_vals)
  } else {
    rho_vals <- vapply(c_values, function(c_val) {
      .compute_rho_generic(
        c_val, theta_vec, beta_vec, lambda_base,
        theta_var = theta_var,
        metric_internal = "msem",
        eta_base = eta_base
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
