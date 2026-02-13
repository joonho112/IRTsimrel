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
# Author: JoonHo Lee
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
                                 metric_internal = c("msem", "info")) {

  metric_internal <- match.arg(metric_internal)

  # Coerce to numeric
  theta_vec <- as.numeric(theta_vec)
  beta_vec  <- as.numeric(beta_vec)
  lambda_base <- as.numeric(lambda_base)

  # Input validation
  if (length(beta_vec) != length(lambda_base)) {
    stop("beta_vec and lambda_base must have the same length.")
  }
  if (length(theta_vec) < 2L) {
    stop("theta_vec must have length >= 2.")
  }
  if (c <= 0) {
    return(0)
  }

  # Use provided theta_var or compute from sample
  if (is.null(theta_var)) {
    theta_var <- var(theta_vec)
  }

  # Validate theta_var
  if (is.na(theta_var) || theta_var <= 1e-10) {
    warning("Invalid theta_var. Falling back to 1.0.")
    theta_var <- 1.0
  }

  I <- length(beta_vec)
  M <- length(theta_vec)

  # Scale discriminations
  lambda_curr <- c * lambda_base

  # Compute test information (vectorized)
  # eta = lambda * (theta - beta)
  eta_base <- outer(theta_vec, beta_vec, "-")              # M x I matrix
  eta_mat  <- sweep(eta_base, 2, lambda_curr, FUN = "*")   # M x I matrix
  p_mat    <- plogis(eta_mat)                              # Response probabilities

  # Item information: lambda^2 * p(1-p)
  info_mat  <- sweep(p_mat * (1 - p_mat), 2, lambda_curr^2, FUN = "*")
  test_info <- rowSums(info_mat)  # Test information at each theta

  if (metric_internal == "info") {
    # Average-information reliability (tilde)
    J_bar <- mean(test_info)
    rho   <- (theta_var * J_bar) / (theta_var * J_bar + 1)
  } else {
    # MSEM-based reliability (bar/w)
    test_info_safe <- pmax(test_info, 1e-10)  # Floor to avoid division by zero
    msem <- mean(1 / test_info_safe)
    rho  <- theta_var / (theta_var + msem)
  }

  rho
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
#' By Jensen's inequality, \eqn{\tilde{\rho} \geq \bar{w}} always holds.
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

  # Coerce to numeric
  theta_vec <- as.numeric(theta_vec)
  beta_vec  <- as.numeric(beta_vec)
  lambda_base <- as.numeric(lambda_base)

  # Input validation
  if (length(beta_vec) != length(lambda_base)) {
    stop("beta_vec and lambda_base must have the same length.")
  }
  if (length(theta_vec) < 2L) {
    stop("theta_vec must have length >= 2.")
  }
  if (c <= 0) {
    return(list(rho_tilde = 0, rho_bar = 0))
  }

  # Use provided theta_var or compute from sample
  if (is.null(theta_var)) {
    theta_var <- var(theta_vec)
  }

  # Validate theta_var
  if (is.na(theta_var) || theta_var <= 1e-10) {
    warning("Invalid theta_var. Falling back to 1.0.")
    theta_var <- 1.0
  }

  # Scale discriminations
  lambda_curr <- c * lambda_base

  # Compute test information (vectorized) -- single pass

  eta_base <- outer(theta_vec, beta_vec, "-")              # M x I matrix
  eta_mat  <- sweep(eta_base, 2, lambda_curr, FUN = "*")   # M x I matrix
  p_mat    <- plogis(eta_mat)                              # Response probabilities

  # Item information: lambda^2 * p(1-p)
  info_mat  <- sweep(p_mat * (1 - p_mat), 2, lambda_curr^2, FUN = "*")
  test_info <- rowSums(info_mat)  # Test information at each theta

  # Average-information reliability (tilde)
  J_bar <- mean(test_info)
  rho_tilde <- (theta_var * J_bar) / (theta_var * J_bar + 1)

  # MSEM-based reliability (bar/w)
  test_info_safe <- pmax(test_info, 1e-10)
  msem <- mean(1 / test_info_safe)
  rho_bar <- theta_var / (theta_var + msem)

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
#' @param verbose Logical. If TRUE, print results.
#'
#' @return An object of class \code{"feasibility_check"} (a list) with:
#' \describe{
#'   \item{\code{rho_range_info}}{Numeric length-2 vector: achievable range of
#'     average-information reliability (\eqn{\tilde{\rho}}).}
#'   \item{\code{rho_range_msem}}{Numeric length-2 vector: achievable range of
#'     MSEM-based reliability (\eqn{\bar{w}}).}
#'   \item{\code{n_items}}{Number of items.}
#'   \item{\code{model}}{Model used.}
#'   \item{\code{latent_shape}}{Latent distribution shape.}
#'   \item{\code{c_bounds}}{Scaling factor bounds evaluated.}
#'   \item{\code{M}}{Monte Carlo sample size.}
#'   \item{\code{theta_var}}{Estimated latent variance.}
#' }
#'
#' @details
#' For the average-information metric (\eqn{\tilde{\rho}}), the reliability
#' is monotone in \eqn{c}, so the range is simply
#' \eqn{[\tilde{\rho}(c_{min}), \tilde{\rho}(c_{max})]}.
#'
#' For the MSEM-based metric (\eqn{\bar{w}}), the reliability can be
#' non-monotone at extreme scaling factors. The function uses
#' \code{\link[stats]{optimize}} to find the maximum within the bounds, and
#' the range is \eqn{[\min(\bar{w}(c_{min}), \bar{w}(c_{max})), \bar{w}_{max}]}.
#'
#' @examples
#' # Check feasibility for 25-item Rasch test
#' feas <- check_feasibility(n_items = 25, model = "rasch", seed = 42, M = 5000)
#' print(feas)
#'
#' # Can we achieve rho = 0.90?
#' 0.90 >= feas$rho_range_info[1] && 0.90 <= feas$rho_range_info[2]
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
                              verbose = TRUE) {

  model <- match.arg(model)

  if (!is.numeric(n_items) || length(n_items) != 1L || n_items <= 0) {
    stop("`n_items` must be a positive integer.")
  }
  n_items <- as.integer(n_items)

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(c_bounds <= 0) || c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  if (!is.numeric(M) || length(M) != 1L || M <= 0) {
    stop("`M` must be a positive integer.")
  }
  M <- as.integer(M)

  # Save/restore RNG state
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv(), inherits = FALSE)
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

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

  # Compute rho_tilde at bounds (monotone, so range is straightforward)
  rho_tilde_low  <- compute_rho_tilde(c_bounds[1], theta_vec, beta_vec,
                                       lambda_base, theta_var = theta_var)
  rho_tilde_high <- compute_rho_tilde(c_bounds[2], theta_vec, beta_vec,
                                       lambda_base, theta_var = theta_var)
  rho_range_info <- c(rho_tilde_low, rho_tilde_high)

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

  rho_range_msem <- c(min(rho_bar_low, rho_bar_high), rho_bar_max)

  result <- list(
    rho_range_info = rho_range_info,
    rho_range_msem = rho_range_msem,
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

  cat("Note: rho_tilde >= rho_bar always (Jensen's inequality).\n")
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
#' # Basic usage: plot reliability curve for 25-item Rasch test
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

  if (!is.numeric(n_items) || length(n_items) != 1L || n_items <= 0) {
    stop("`n_items` must be a positive integer.")
  }
  n_items <- as.integer(n_items)

  if (!is.numeric(c_values) || length(c_values) < 2L || any(c_values <= 0)) {
    stop("`c_values` must be a numeric vector of positive values with length >= 2.")
  }

  if (!is.numeric(M) || length(M) != 1L || M <= 0) {
    stop("`M` must be a positive integer.")
  }
  M <- as.integer(M)

  # Save/restore RNG state
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv(), inherits = FALSE)
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }

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
  n_c <- length(c_values)

  if (metric == "both") {
    results <- lapply(c_values, function(c_val) {
      compute_rho_both(c_val, theta_vec, beta_vec, lambda_base, theta_var = theta_var)
    })
    df <- data.frame(
      c = c_values,
      rho_tilde = vapply(results, `[[`, numeric(1), "rho_tilde"),
      rho_bar   = vapply(results, `[[`, numeric(1), "rho_bar")
    )
  } else if (metric == "info") {
    rho_vals <- vapply(c_values, function(c_val) {
      compute_rho_tilde(c_val, theta_vec, beta_vec, lambda_base, theta_var = theta_var)
    }, numeric(1))
    df <- data.frame(c = c_values, rho_tilde = rho_vals)
  } else {
    rho_vals <- vapply(c_values, function(c_val) {
      compute_rho_bar(c_val, theta_vec, beta_vec, lambda_base, theta_var = theta_var)
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
