# =============================================================================
# reliability_utils.R
# =============================================================================
# Exported Low-Level Reliability Computation Functions
# These functions can be shared by EQC and SPC for consistency.
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
#' inside both `eqc_calibrate()` and `spc_calibrate()`:
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
#' # With pre-calculated theta variance (recommended for SPC)
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
