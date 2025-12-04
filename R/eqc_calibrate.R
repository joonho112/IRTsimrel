# =============================================================================
# eqc_calibrate.R
# =============================================================================
# Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)
#
# Contents:
#   - eqc_calibrate(): Main calibration function
#   - print.eqc_result(): Print method
#   - summary.eqc_result(): Summary method
#
# Dependencies:
#   - sim_latentG() from sim_latentG.R
#   - sim_item_params() from sim_item_params.R
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================


#' Empirical Quadrature Calibration (Algorithm 1: EQC/SQC)
#'
#' @description
#' `eqc_calibrate()` implements Algorithm 1 (Empirical / Stochastic
#' Quadrature Calibration, EQC/SQC) for reliability-targeted IRT simulation.
#'
#' Given a target marginal reliability \eqn{\rho^*}, a latent distribution
#' generator `sim_latentG()` (for \eqn{G}) and an item parameter generator
#' `sim_item_params()` (for \eqn{H}), the function searches for a global
#' discrimination scale \eqn{c^* > 0} such that the population reliability
#' \eqn{\rho(c)} of the Rasch/2PL model is approximately equal to \eqn{\rho^*}.
#'
#' The key idea is to:
#' \enumerate{
#'   \item Draw a large fixed "quadrature" sample
#'     \eqn{\{\theta_m\}_{m=1}^M \sim G} and item parameters
#'     \eqn{\{(\beta_i, \lambda_{i,0})\}_{i=1}^I \sim H} once.
#'   \item For any scale \eqn{c}, form \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}
#'     and compute the empirical approximation to population reliability
#'     \eqn{\hat\rho_M(c)} from the test information function.
#'   \item Solve the scalar equation \eqn{\hat\rho_M(c^*) = \rho^*} using
#'     deterministic root-finding (Brent's method via \code{uniroot()}).
#' }
#'
#' @param target_rho Numeric in (0, 1). Target marginal reliability \eqn{\rho^*}.
#'
#' @param n_items Integer. Number of items in the test form.
#'
#' @param model Character. Measurement model: \code{"rasch"} or \code{"2pl"}.
#'   For \code{"rasch"}, all baseline discriminations are set to 1 before scaling.
#'
#' @param latent_shape Character. Shape argument passed to \code{sim_latentG()}
#'   (e.g. \code{"normal"}, \code{"bimodal"}, \code{"heavy_tail"}, ...).
#'
#' @param item_source Character. Source argument passed to \code{sim_item_params()}
#'   (e.g. \code{"irw"}, \code{"parametric"}, \code{"hierarchical"}, \code{"custom"}).
#'
#' @param latent_params List. Additional arguments passed to \code{sim_latentG()}.
#'
#' @param item_params List. Additional arguments passed to \code{sim_item_params()}.
#'
#' @param reliability_metric Character. Reliability definition used inside EQC:
#'   \describe{
#'     \item{\code{"msem"}}{MSEM-based marginal reliability (default, theoretically exact).}
#'     \item{\code{"info"}}{Average-information reliability (faster, more stable).}
#'   }
#'   Synonyms: \code{"bar"} for \code{"msem"}, \code{"tilde"} for \code{"info"}.
#'
#' @param M Integer. Size of the empirical quadrature sample (default: 10000).
#'
#' @param c_bounds Numeric length-2 vector. Search bounds for \eqn{c}. Default: c(0.3, 3).
#'
#' @param tol Numeric. Tolerance for \code{uniroot()}. Default: 1e-4.
#'
#' @param seed Optional integer for reproducibility.
#'
#' @param verbose Logical. If TRUE, print progress messages.
#'
#' @return An object of class \code{"eqc_result"} (a list) with elements:
#' \describe{
#'   \item{\code{c_star}}{Calibrated discrimination scale \eqn{c^*}.}
#'   \item{\code{target_rho}}{Target reliability \eqn{\rho^*}.}
#'   \item{\code{achieved_rho}}{Empirical quadrature estimate \eqn{\hat\rho_M(c^*)}.}
#'   \item{\code{metric}}{Reliability metric used.}
#'   \item{\code{theta_quad}}{Length-M vector of quadrature abilities.}
#'   \item{\code{theta_var}}{Sample variance of theta_quad.}
#'   \item{\code{items_base}}{item_params object with scale = 1 (baseline).}
#'   \item{\code{items_calib}}{item_params object with discriminations scaled by c_star.}
#' }
#'
#' @details
#' ## Reliability Metrics
#'
#' The function supports two reliability definitions:
#'
#' - **MSEM-based** (\code{"msem"}/\code{"bar"}): Uses the harmonic mean of test information,
#'   \eqn{\bar{w}(c) = \sigma^2_\theta / (\sigma^2_\theta + E[1/\mathcal{J}(\theta;c)])}.
#'   This is theoretically exact but may have a lower ceiling for high reliability.
#'
#' - **Average-information** (\code{"info"}/\code{"tilde"}): Uses the arithmetic mean,
#'   \eqn{\tilde{\rho}(c) = \sigma^2_\theta \bar{\mathcal{J}}(c) / (\sigma^2_\theta \bar{\mathcal{J}}(c) + 1)}.
#'   By Jensen's inequality, \eqn{\tilde{\rho} \geq \bar{w}}, so this metric typically
#'   yields higher reliability values.
#'
#' ## WLE vs EAP Reliability Interpretation
#'
#' When validating with TAM, note that EAP reliability is systematically higher than
#' WLE reliability. This is not a bug but a mathematical property of TAM's definitions.
#' EAP reliability more directly corresponds to the MSEM-based population reliability
#' targeted by EQC. For conservative inference, treat WLE as a lower bound and EAP
#' as an upper bound for true measurement precision.
#'
#' @seealso
#' \code{\link{spc_calibrate}} for the stochastic approximation alternative,
#' \code{\link{compute_rho_bar}} and \code{\link{compute_rho_tilde}} for
#' reliability computation utilities,
#' \code{\link{compute_reliability_tam}} for TAM validation.
#'
#' @examples
#' \dontrun{
#' # Basic EQC calibration
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   latent_shape = "normal",
#'   item_source = "irw",
#'   seed = 42,
#'   verbose = TRUE
#' )
#' print(eqc_result)
#' }
#'
#' @export
eqc_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl"),
                          latent_shape = "normal",
                          item_source = "irw",
                          latent_params = list(),
                          item_params = list(),
                          reliability_metric = c("msem", "info", "bar", "tilde"),
                          M = 10000L,
                          c_bounds = c(0.3, 3),
                          tol = 1e-4,
                          seed = NULL,
                          verbose = FALSE) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================

  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a scalar in (0, 1).")
  }

  if (!is.numeric(n_items) || length(n_items) != 1L || n_items <= 0) {
    stop("`n_items` must be a positive integer.")
  }
  n_items <- as.integer(n_items)

  model <- match.arg(model)

  reliability_metric <- match.arg(reliability_metric)
  metric_internal <- if (reliability_metric %in% c("msem", "bar")) "msem" else "info"

  if (!is.numeric(M) || length(M) != 1L || M <= 0) {
    stop("`M` must be a positive integer.")
  }
  M <- as.integer(M)

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(c_bounds <= 0) || c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # ===========================================================================
  # Step 1: Generate Quadrature Samples from G and H
  # ===========================================================================

  if (verbose) cat("Step 1: Generating quadrature samples...\n")

  # --- Latent abilities ---
  latent_args <- modifyList(
    list(n = M, shape = latent_shape),
    latent_params
  )
  quad_latent <- do.call(sim_latentG, latent_args)
  theta_quad  <- quad_latent$theta
  theta_var   <- var(theta_quad)

  # --- Item parameters (baseline scale = 1) ---
  item_args <- modifyList(
    list(
      n_items = n_items,
      model   = model,
      source  = item_source,
      scale   = 1
    ),
    item_params
  )
  items_base <- do.call(sim_item_params, item_args)
  beta_vec   <- items_base$data$beta

  # CRITICAL: Get unscaled lambda for proper calibration
  if (model == "rasch") {
    lambda_base <- rep(1, n_items)
  } else {
    # Prefer lambda_unscaled if available, else use lambda
    if ("lambda_unscaled" %in% names(items_base$data)) {
      lambda_base <- items_base$data$lambda_unscaled
    } else {
      lambda_base <- items_base$data$lambda
    }
  }

  if (length(lambda_base) != n_items) {
    stop("Length of baseline discrimination vector does not match n_items.")
  }

  if (verbose) {
    cat(sprintf("  M (quad persons) = %d\n", M))
    cat(sprintf("  I (items)        = %d\n", n_items))
    cat(sprintf("  theta: mean = %.3f, sd = %.3f, var = %.3f\n",
                mean(theta_quad), sd(theta_quad), theta_var))
    cat(sprintf("  beta:  mean = %.3f, sd = %.3f\n",
                mean(beta_vec), sd(beta_vec)))
    cat(sprintf("  lambda_base: mean = %.3f, sd = %.3f\n",
                mean(lambda_base), sd(lambda_base)))
    cat(sprintf("  metric = %s\n", metric_internal))
  }

  # Precompute theta - beta grid (M x I)
  eta_base <- outer(theta_quad, beta_vec, "-")  # theta_p - beta_i

  # ===========================================================================
  # Step 2: Define Reliability Function rho(c)
  # ===========================================================================

  compute_rho <- function(c) {
    lambda_curr <- c * lambda_base

    # eta = lambda_i * (theta - beta_i)
    eta_mat <- sweep(eta_base, 2, lambda_curr, FUN = "*")
    p_mat   <- plogis(eta_mat)

    # Item information: lambda^2 * p(1-p)
    info_mat <- sweep(p_mat * (1 - p_mat), 2, lambda_curr^2, FUN = "*")
    test_info <- rowSums(info_mat)

    if (metric_internal == "info") {
      # Average test information reliability (tilde)
      J_bar <- mean(test_info)
      rho   <- (theta_var * J_bar) / (theta_var * J_bar + 1)
    } else {
      # MSEM-based reliability (bar/w)
      test_info_safe <- pmax(test_info, 1e-10)
      msem <- mean(1 / test_info_safe)
      rho  <- theta_var / (theta_var + msem)
    }

    rho
  }

  objective_fn <- function(c) compute_rho(c) - target_rho

  # ===========================================================================
  # Step 3: Root-Finding with Bracket Checks
  # ===========================================================================

  if (verbose) cat("Step 2: Running root-finding algorithm...\n")

  val_low  <- objective_fn(c_bounds[1])
  val_high <- objective_fn(c_bounds[2])

  rho_at_low  <- val_low + target_rho
  rho_at_high <- val_high + target_rho

  if (verbose) {
    cat(sprintf("  At c = %.3f: rho = %.4f, g = %.4f\n", c_bounds[1], rho_at_low, val_low))
    cat(sprintf("  At c = %.3f: rho = %.4f, g = %.4f\n", c_bounds[2], rho_at_high, val_high))
  }

  root_status <- "ok"

  if (val_low > 0 && val_high > 0) {
    # Target is too LOW: even at minimum c, reliability exceeds target
    warning(sprintf(
      paste0(
        "Target rho* = %.3f is below the minimum achievable reliability.\n",
        "  - At c = %.2f (lower bound): rho = %.4f\n",
        "Suggestion: Lower c_bounds[1] or increase target_rho.\n",
        "Returning c_star = %.2f (lower bound)."
      ),
      target_rho, c_bounds[1], rho_at_low, c_bounds[1]
    ))
    c_star <- c_bounds[1]
    root_status <- "below_lower"

  } else if (val_low < 0 && val_high < 0) {
    # Target is too HIGH: even at maximum c, reliability falls short
    warning(sprintf(
      paste0(
        "Target rho* = %.3f exceeds the maximum achievable reliability for this configuration.\n",
        "  - Items: %d\n",
        "  - At c = %.2f (upper bound): rho = %.4f\n",
        "  - Gap: %.4f\n\n",
        "Suggestions:\n",
        "  1. Use reliability_metric = 'info' (typically yields higher values)\n",
        "  2. Increase c_bounds[2] beyond %.1f\n",
        "  3. Increase n_items for higher achievable reliability\n",
        "  4. Accept achieved rho = %.4f as maximum for this design\n\n",
        "Returning c_star = %.2f (upper bound)."
      ),
      target_rho, n_items, c_bounds[2], rho_at_high,
      target_rho - rho_at_high, c_bounds[2], rho_at_high, c_bounds[2]
    ))
    c_star <- c_bounds[2]
    root_status <- "above_upper"

  } else {
    # Normal case: root exists within bounds
    root_res <- uniroot(objective_fn, interval = c_bounds, tol = tol)
    c_star   <- root_res$root
    root_status <- "uniroot_success"
  }

  achieved_rho <- compute_rho(c_star)

  # Informative message if target is near the achievable maximum
  if (root_status == "uniroot_success" && target_rho > rho_at_high * 0.95) {
    message(sprintf(
      "Note: Target rho* = %.3f is near the achievable maximum (%.3f) for this configuration.",
      target_rho, rho_at_high
    ))
  }

  if (verbose) {
    cat(sprintf("  c* = %.6f\n", c_star))
    cat(sprintf("  Target rho    = %.4f\n", target_rho))
    cat(sprintf("  Achieved rho  = %.4f\n", achieved_rho))
    cat(sprintf("  Root status   = %s\n", root_status))
  }

  # ===========================================================================
  # Step 4: Construct Calibrated Item Parameters
  # ===========================================================================

  items_calib <- items_base
  items_calib$scale <- items_base$scale * c_star
  items_calib$data$lambda <- lambda_base * c_star

  # ===========================================================================
  # Assemble Result
  # ===========================================================================

  res <- list(
    c_star       = c_star,
    target_rho   = target_rho,
    achieved_rho = achieved_rho,
    metric       = metric_internal,
    model        = model,
    n_items      = n_items,
    M            = M,
    theta_quad   = theta_quad,
    theta_var    = theta_var,
    beta_vec     = beta_vec,
    lambda_base  = lambda_base,
    lambda_scaled = lambda_base * c_star,
    items_base   = items_base,
    items_calib  = items_calib,
    call         = match.call(),
    misc         = list(
      c_bounds    = c_bounds,
      tol         = tol,
      root_status = root_status,
      val_low     = val_low,
      val_high    = val_high,
      rho_bounds  = c(rho_L = rho_at_low, rho_U = rho_at_high)
    )
  )

  class(res) <- c("eqc_result", "list")
  res
}


# =============================================================================
# S3 Methods for eqc_result
# =============================================================================

#' @export
print.eqc_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("  Empirical Quadrature Calibration (EQC) Results\n")
  cat("=======================================================\n\n")

  cat("Calibration Summary:\n")
  cat(sprintf("  Model                        : %s\n", toupper(x$model)))
  cat(sprintf("  Target reliability (rho*)    : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved reliability         : %.*f\n", digits, x$achieved_rho))
  cat(sprintf("  Absolute error               : %.2e\n", abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor (c*)          : %.*f\n", digits, x$c_star))
  cat("\n")

  cat("Design Parameters:\n")
  cat(sprintf("  Number of items (I)          : %d\n", x$n_items))
  cat(sprintf("  Quadrature points (M)        : %d\n", x$M))
  cat(sprintf("  Reliability metric           : %s\n",
              ifelse(x$metric == "info", "Average-information (tilde)", "MSEM-based (bar/w)")))
  cat(sprintf("  Latent variance              : %.4f\n", x$theta_var))
  cat("\n")

  cat("Convergence:\n")
  cat(sprintf("  Root status                  : %s\n", x$misc$root_status))
  cat(sprintf("  Search bracket               : [%.3f, %.3f]\n",
              x$misc$c_bounds[1], x$misc$c_bounds[2]))
  cat(sprintf("  Bracket reliabilities        : [%.4f, %.4f]\n",
              x$misc$rho_bounds[1], x$misc$rho_bounds[2]))
  cat("\n")

  cat("Parameter Summaries:\n")
  cat(sprintf("  theta:        mean = %.3f, sd = %.3f\n",
              mean(x$theta_quad), sd(x$theta_quad)))
  cat(sprintf("  beta:         mean = %.3f, sd = %.3f, range = [%.2f, %.2f]\n",
              mean(x$beta_vec), sd(x$beta_vec), min(x$beta_vec), max(x$beta_vec)))
  cat(sprintf("  lambda_base:  mean = %.3f, sd = %.3f\n",
              mean(x$lambda_base), sd(x$lambda_base)))
  cat(sprintf("  lambda_scaled: mean = %.3f, sd = %.3f\n",
              mean(x$lambda_scaled), sd(x$lambda_scaled)))
  cat("\n")

  invisible(x)
}


#' @export
summary.eqc_result <- function(object, ...) {
  print(object, ...)
}
