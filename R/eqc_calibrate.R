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
#'   (e.g. \code{"parametric"}, \code{"irw"}, \code{"hierarchical"}, \code{"custom"}).
#'   Defaults to \code{"parametric"} since the \pkg{irw} package is an optional
#'   dependency (listed in Suggests). Use \code{"irw"} for empirically-grounded
#'   difficulties when the \pkg{irw} package is installed.
#'
#' @param latent_params List. Additional arguments passed to \code{sim_latentG()}.
#'
#' @param item_params List. Additional arguments passed to \code{sim_item_params()}.
#'
#' @param reliability_metric Character. Reliability definition used inside EQC:
#'   \describe{
#'     \item{\code{"info"}}{Average-information reliability (default, recommended for EQC).
#'       Targets \eqn{\tilde{\rho}}, which is guaranteed to be monotone in \eqn{c}.}
#'     \item{\code{"msem"}}{MSEM-based marginal reliability (theoretically exact, but may
#'       have a non-monotone objective for EQC; see Details).}
#'   }
#'   Synonyms: \code{"tilde"} for \code{"info"}, \code{"bar"} for \code{"msem"}.
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
#' - **Average-information** (\code{"info"}/\code{"tilde"}, **default**): Uses the arithmetic mean,
#'   \eqn{\tilde{\rho}(c) = \sigma^2_\theta \bar{\mathcal{J}}(c) / (\sigma^2_\theta \bar{\mathcal{J}}(c) + 1)}.
#'   By Jensen's inequality, \eqn{\tilde{\rho} \geq \bar{w}}, so this metric typically
#'   yields higher reliability values. **This is the recommended default for EQC** because
#'   the objective function \eqn{\tilde{\rho}(c) - \rho^*} is guaranteed to be monotone
#'   in \eqn{c}, ensuring that \code{uniroot()} will find the unique root.
#'
#' - **MSEM-based** (\code{"msem"}/\code{"bar"}): Uses the harmonic mean of test information,
#'   \eqn{\bar{w}(c) = \sigma^2_\theta / (\sigma^2_\theta + E[1/\mathcal{J}(\theta;c)])}.
#'   This is theoretically exact but the objective \eqn{\bar{w}(c) - \rho^*} may be
#'   non-monotone (see Lee, 2025, Section 4.3), which can cause \code{uniroot()} to
#'   fail or find an incorrect root. Use \code{\link{sac_calibrate}} if you need
#'   to target \eqn{\bar{w}} directly.
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
#' \code{\link{sac_calibrate}} for the stochastic approximation alternative,
#' \code{\link{compute_rho_bar}} and \code{\link{compute_rho_tilde}} for
#' reliability computation utilities,
#' \code{\link{compute_reliability_tam}} for TAM validation.
#'
#' @examples
#' # Basic EQC calibration with parametric items (fast)
#' \donttest{
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   latent_shape = "normal",
#'   item_source = "parametric",
#'   M = 5000L,
#'   seed = 42
#' )
#' print(eqc_result)
#' }
#'
#' \dontrun{
#' # EQC with IRW difficulties (requires irw package)
#' eqc_result2 <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   item_source = "irw",
#'   seed = 42,
#'   verbose = TRUE
#' )
#' }
#'
#' @export
eqc_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl"),
                          latent_shape = "normal",
                          item_source = "parametric",
                          latent_params = list(),
                          item_params = list(),
                          reliability_metric = c("info", "tilde", "msem", "bar"),
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

  # Warn about non-monotone MSEM objective for EQC (Lee, 2025, Section 4.3)
  if (metric_internal == "msem") {
    warning(
      "MSEM-based reliability (w-bar) may have a non-monotone objective for EQC ",
      "(see Lee, 2025, Section 4.3). Consider using 'info' (rho-tilde) for EQC, ",
      "or use sac_calibrate() which handles w-bar targeting correctly."
    )
  }

  if (!is.numeric(M) || length(M) != 1L || M <= 0) {
    stop("`M` must be a positive integer.")
  }
  M <- as.integer(M)

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(c_bounds <= 0) || c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

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

  # ===========================================================================
  # Step 1: Generate Quadrature Samples from G and H
  # ===========================================================================

  # ===========================================================================
  # Convenience parsing: auto-wrap shape params in latent_params
  # ===========================================================================
  known_shape_params <- c("delta", "df", "k", "mu_mix", "sigma_mix", "weights",
                          "w0", "m", "m1", "m2", "w_inner",
                          "w_floor", "m_floor", "w_ceil", "m_ceil")
  top_level_shape <- intersect(names(latent_params), known_shape_params)
  if (length(top_level_shape) > 0 && is.null(latent_params$shape_params)) {
    shape_vals <- latent_params[top_level_shape]
    latent_params[top_level_shape] <- NULL
    latent_params$shape_params <- shape_vals
    message(
      "Auto-wrapping shape parameter(s) {",
      paste(top_level_shape, collapse = ", "),
      "} into latent_params$shape_params."
    )
  }

  if (verbose) message("Step 1: Generating quadrature samples...")

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
    message(sprintf("  M (quad persons) = %d", M))
    message(sprintf("  I (items)        = %d", n_items))
    message(sprintf("  theta: mean = %.3f, sd = %.3f, var = %.3f",
                mean(theta_quad), sd(theta_quad), theta_var))
    message(sprintf("  beta:  mean = %.3f, sd = %.3f",
                mean(beta_vec), sd(beta_vec)))
    message(sprintf("  lambda_base: mean = %.3f, sd = %.3f",
                mean(lambda_base), sd(lambda_base)))
    message(sprintf("  metric = %s", metric_internal))
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

  if (verbose) message("Step 2: Running root-finding algorithm...")

  val_low  <- objective_fn(c_bounds[1])
  val_high <- objective_fn(c_bounds[2])

  rho_at_low  <- val_low + target_rho
  rho_at_high <- val_high + target_rho

  if (verbose) {
    message(sprintf("  At c = %.3f: rho = %.4f, g = %.4f", c_bounds[1], rho_at_low, val_low))
    message(sprintf("  At c = %.3f: rho = %.4f, g = %.4f", c_bounds[2], rho_at_high, val_high))
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
    uniroot_result <- root_res
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
    message(sprintf("  c* = %.6f", c_star))
    message(sprintf("  Target rho    = %.4f", target_rho))
    message(sprintf("  Achieved rho  = %.4f", achieved_rho))
    message(sprintf("  Root status   = %s", root_status))
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
      rho_bounds  = c(rho_L = rho_at_low, rho_U = rho_at_high),
      uniroot_result = if (exists("uniroot_result", inherits = FALSE)) uniroot_result else NULL
    )
  )

  class(res) <- c("eqc_result", "list")
  res
}


# =============================================================================
# S3 Methods for eqc_result
# =============================================================================

#' @rdname eqc_calibrate
#' @param x An object of class \code{"eqc_result"}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments passed to or from other methods.
#' @return The input object, invisibly.
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


#' @rdname eqc_calibrate
#' @param object An object of class \code{"eqc_result"}.
#' @return An object of class \code{"summary.eqc_result"} containing key
#'   calibration results.
#' @export
summary.eqc_result <- function(object, ...) {
  out <- list(
    c_star       = object$c_star,
    target_rho   = object$target_rho,
    achieved_rho = object$achieved_rho,
    metric       = object$metric,
    model        = object$model,
    n_items      = object$n_items,
    M            = object$M,
    root_status  = object$misc$root_status,
    theta_var    = object$theta_var
  )
  class(out) <- "summary.eqc_result"
  out
}


#' Print Method for summary.eqc_result Objects
#'
#' @param x A \code{summary.eqc_result} object from \code{summary.eqc_result()}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.summary.eqc_result <- function(x, digits = 4, ...) {
  cat("Summary: Empirical Quadrature Calibration (EQC)\n")
  cat("================================================\n")
  cat(sprintf("  Model            : %s\n", toupper(x$model)))
  cat(sprintf("  Metric           : %s\n",
              ifelse(x$metric == "info", "Average-information (tilde)", "MSEM-based (bar/w)")))
  cat(sprintf("  Number of items  : %d\n", x$n_items))
  cat(sprintf("  Quadrature (M)   : %d\n", x$M))
  cat(sprintf("  Latent variance  : %.*f\n", digits, x$theta_var))
  cat("\nCalibration Results:\n")
  cat(sprintf("  Target rho*      : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved rho     : %.*f\n", digits, x$achieved_rho))
  cat(sprintf("  Absolute error   : %.2e\n", abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor c*: %.*f\n", digits, x$c_star))
  cat(sprintf("  Root status      : %s\n", x$root_status))
  invisible(x)
}


#' Extract Calibrated Item Parameters from EQC Results
#'
#' @description
#' Returns the calibrated item parameters as a data frame.
#'
#' @param object An object of class \code{"eqc_result"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{item_id}}{Item identifier (1 to I).}
#'   \item{\code{beta}}{Item difficulty.}
#'   \item{\code{lambda_base}}{Baseline (unscaled) discrimination.}
#'   \item{\code{lambda_scaled}}{Scaled discrimination (\code{lambda_base * c*}).}
#'   \item{\code{c_star}}{Calibrated scaling factor (same for all items).}
#' }
#'
#' @examples
#' \donttest{
#' eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 25,
#'                           model = "rasch", seed = 42, M = 5000)
#' coef(eqc_res)
#' }
#'
#' @export
coef.eqc_result <- function(object, ...) {
  data.frame(
    item_id = seq_len(object$n_items),
    beta = object$beta_vec,
    lambda_base = object$lambda_base,
    lambda_scaled = object$lambda_scaled,
    c_star = rep(object$c_star, object$n_items)
  )
}


#' Predict Reliability at New Scaling Factor Values (EQC)
#'
#' @description
#' If \code{newdata} is NULL, returns the achieved reliability from calibration.
#' If \code{newdata} is a numeric vector of scaling factor values, computes
#' the reliability at each value using the stored quadrature sample and item
#' parameters.
#'
#' @param object An object of class \code{"eqc_result"}.
#' @param newdata Optional numeric vector of scaling factor values.
#'   If NULL, returns \code{object$achieved_rho}.
#' @param ... Additional arguments (ignored).
#'
#' @return If \code{newdata} is NULL, a single numeric value (achieved reliability).
#'   If \code{newdata} is a numeric vector, a named numeric vector of reliability
#'   values at each scaling factor.
#'
#' @examples
#' \donttest{
#' eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 25,
#'                           model = "rasch", seed = 42, M = 5000)
#' predict(eqc_res)
#' predict(eqc_res, newdata = c(0.5, 1.0, 1.5, 2.0))
#' }
#'
#' @export
predict.eqc_result <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(object$achieved_rho)
  }

  c_values <- as.numeric(newdata)
  if (any(c_values <= 0)) {
    stop("All values in `newdata` must be positive.")
  }

  rho_fn <- if (object$metric == "info") compute_rho_tilde else compute_rho_bar

  result <- vapply(c_values, function(c_val) {
    rho_fn(c_val, object$theta_quad, object$beta_vec,
           object$lambda_base, theta_var = object$theta_var)
  }, numeric(1))

  names(result) <- paste0("c=", format(c_values, digits = 4))
  result
}
