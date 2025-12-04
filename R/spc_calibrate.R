# =============================================================================
# spc_calibrate.R
# =============================================================================
# Stochastic Approximation Calibration (Algorithm 2: SPC/SAC)
#
# Contents:
#   - spc_calibrate(): Main calibration function
#   - sac_calibrate(): Alias for spc_calibrate
#   - print.spc_result(): Print method
#   - summary.spc_result(): Summary method
#   - plot.spc_result(): Plot method
#
# Dependencies:
#   - sim_latentG() from sim_latentG.R
#   - sim_item_params() from sim_item_params.R
#   - .compute_rho_generic() from reliability_utils.R
#   - compute_apc_init() from reliability_utils.R
#
# Author: JoonHo Lee
# Date: December 2025
# =============================================================================


#' Stochastic Approximation Calibration (Algorithm 2: SPC/SAC)
#'
#' @description
#' `spc_calibrate()` implements Algorithm 2 (Stochastic Approximation
#' Calibration, SPC/SAC) for reliability-targeted IRT simulation using the
#' Robbins-Monro stochastic approximation algorithm.
#'
#' Given a target marginal reliability \eqn{\rho^*}, a latent distribution
#' generator `sim_latentG()` (for \eqn{G}) and an item parameter generator
#' `sim_item_params()` (for \eqn{H}), the function iteratively searches for a
#' global discrimination scale \eqn{c^* > 0} such that the population
#' reliability \eqn{\rho(c)} of the Rasch/2PL model is approximately equal to
#' \eqn{\rho^*}.
#'
#' SPC complements EQC (Algorithm 1) by:
#' \enumerate{
#'   \item Providing an independent validation of EQC calibration results.
#'   \item Enabling calibration to the exact marginal reliability \eqn{\bar{w}}
#'     (not just the average-information approximation \eqn{\tilde{\rho}}).
#'   \item Handling complex data-generating processes where analytic information
#'     functions may be unavailable.
#' }
#'
#' The algorithm uses the Robbins-Monro update rule:
#' \deqn{c_{n+1} = c_n - a_n \cdot (\hat{\rho}_n - \rho^*)}
#'
#' where \eqn{a_n = a / (n + A)^\gamma} is a decreasing step size sequence
#' satisfying \eqn{\sum a_n = \infty} and \eqn{\sum a_n^2 < \infty}.
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
#' @param reliability_metric Character. Reliability definition used inside SPC:
#'   \describe{
#'     \item{\code{"msem"}}{MSEM-based marginal reliability (theoretically exact, targets \eqn{\bar{w}}).}
#'     \item{\code{"info"}}{Average-information reliability (faster, targets \eqn{\tilde{\rho}}).}
#'   }
#'   Synonyms: \code{"bar"} for \code{"msem"}, \code{"tilde"} for \code{"info"}.
#'
#' @param c_init Numeric, \code{eqc_result} object, or NULL. Initial value for
#'   the scaling factor \eqn{c_0}.
#'   \itemize{
#'     \item If an \code{eqc_result} object is provided, its \code{c_star} is used (warm start).
#'     \item If a numeric value is provided, it is used directly.
#'     \item If NULL, initialized using Analytic Pre-Calibration (APC).
#'   }
#'   Providing a warm start from EQC greatly accelerates convergence.
#'
#' @param M_per_iter Integer. Number of Monte Carlo samples per iteration for
#'   estimating reliability. Default: 500. Larger values reduce variance but
#'   increase computation time.
#'
#' @param M_pre Integer. Number of Monte Carlo samples for pre-calculating the
#'   latent variance \eqn{\sigma^2_\theta}. Default: 10000. This variance is
#'   fixed throughout the iterations for stability. This is a CRITICAL parameter
#'   for numerical stability.
#'
#' @param n_iter Integer. Total number of Robbins-Monro iterations. Default: 300.
#'
#' @param burn_in Integer. Number of initial iterations to discard before
#'   Polyak-Ruppert averaging. Default: \code{floor(n_iter / 2)}.
#'
#' @param step_params List. Parameters controlling the step size sequence:
#'   \describe{
#'     \item{\code{a}}{Base step size (default: 1.0)}
#'     \item{\code{A}}{Stabilization constant (default: 50)}
#'     \item{\code{gamma}}{Decay exponent (default: 0.67, i.e., 2/3)}
#'   }
#'
#' @param c_bounds Numeric length-2 vector. Projection bounds for \eqn{c}.
#'   Iterates are clipped to this interval after each update. Default: c(0.01, 20).
#'
#' @param resample_items Logical. If TRUE (default), resample item parameters
#'   at each iteration. If FALSE, fix item parameters across all iterations
#'   (reduces variance but may introduce bias).
#'
#' @param seed Optional integer for reproducibility.
#'
#' @param verbose Logical or integer. If TRUE or >= 1, print progress messages.
#'   If >= 2, print detailed iteration-level output.
#'
#' @return An object of class \code{"spc_result"} (a list) with elements:
#' \describe{
#'   \item{\code{c_star}}{Calibrated discrimination scale (Polyak-Ruppert average).}
#'   \item{\code{c_final}}{Final iterate \eqn{c_{n_{iter}}}.}
#'   \item{\code{target_rho}}{Target reliability \eqn{\rho^*}.}
#'   \item{\code{achieved_rho}}{Estimated reliability at \eqn{c^*} (post-calibration).}
#'   \item{\code{theta_var}}{Pre-calculated latent variance used throughout.}
#'   \item{\code{trajectory}}{Numeric vector of all iterates.}
#'   \item{\code{rho_trajectory}}{Numeric vector of reliability estimates.}
#'   \item{\code{init_method}}{Character indicating initialization method.}
#'   \item{\code{metric}}{Reliability metric used.}
#'   \item{\code{convergence}}{List with convergence diagnostics.}
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}} for the faster deterministic Algorithm 1,
#' \code{\link{compute_rho_bar}} and \code{\link{compute_rho_tilde}} for
#' reliability computation utilities,
#' \code{\link{compare_eqc_spc}} for comparing EQC and SPC results.
#'
#' @references
#' Robbins, H., & Monro, S. (1951). A stochastic approximation method.
#'   \emph{The Annals of Mathematical Statistics, 22}(3), 400--407.
#'
#' Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
#'   approximation by averaging. \emph{SIAM Journal on Control and Optimization,
#'   30}(4), 838--855.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic SPC calibration
#' spc_result <- spc_calibrate(
#'   target_rho = 0.75,
#'   n_items = 20,
#'   model = "rasch",
#'   n_iter = 200,
#'   seed = 12345,
#'   verbose = TRUE
#' )
#' print(spc_result)
#' plot(spc_result)
#'
#' # Example 2: Warm start from EQC (RECOMMENDED)
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "2pl",
#'   M = 10000,
#'   seed = 42
#' )
#'
#' spc_result <- spc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "2pl",
#'   c_init = eqc_result,  # Direct EQC object passing!
#'   n_iter = 100,
#'   seed = 42
#' )
#'
#' # Compare EQC and SPC results
#' cat(sprintf("EQC c* = %.4f, SPC c* = %.4f\n",
#'             eqc_result$c_star, spc_result$c_star))
#' }
#'
#' @export
spc_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl"),
                          latent_shape = "normal",
                          item_source = "irw",
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

  if (!is.numeric(M_per_iter) || length(M_per_iter) != 1L || M_per_iter <= 0) {
    stop("`M_per_iter` must be a positive integer.")
  }
  M_per_iter <- as.integer(M_per_iter)

  if (!is.numeric(M_pre) || length(M_pre) != 1L || M_pre <= 0) {
    stop("`M_pre` must be a positive integer.")
  }
  M_pre <- as.integer(M_pre)

  if (!is.numeric(n_iter) || length(n_iter) != 1L || n_iter <= 0) {
    stop("`n_iter` must be a positive integer.")
  }
  n_iter <- as.integer(n_iter)

  # Set default burn-in
  if (is.null(burn_in)) {
    burn_in <- floor(n_iter / 2)
  }
  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
      burn_in < 0 || burn_in >= n_iter) {
    stop("`burn_in` must be a non-negative integer less than `n_iter`.")
  }
  burn_in <- as.integer(burn_in)

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(c_bounds <= 0) || c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  # Set default step parameters
  default_step <- list(a = 1.0, A = 50, gamma = 0.67)
  step_params <- modifyList(default_step, step_params)

  if (step_params$a <= 0) stop("step_params$a must be positive.")
  if (step_params$A < 0) stop("step_params$A must be non-negative.")
  if (step_params$gamma <= 0.5 || step_params$gamma > 1) {
    warning("step_params$gamma outside recommended range (0.5, 1]. Using ", step_params$gamma)
  }

  verbose_level <- if (is.logical(verbose)) {
    as.integer(verbose)
  } else {
    as.integer(verbose)
  }

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # ===========================================================================
  # Step 0: Pre-Calculate Latent Variance (CRITICAL for stability)
  # ===========================================================================

  if (verbose_level >= 1) {
    cat("\n")
    cat("================================================================\n")
    cat("  Stochastic Approximation Calibration (SPC/SAC)\n")
    cat("================================================================\n\n")
    cat("Configuration:\n")
    cat(sprintf("  Target reliability    : %.4f\n", target_rho))
    cat(sprintf("  Number of items       : %d\n", n_items))
    cat(sprintf("  Model                 : %s\n", toupper(model)))
    cat(sprintf("  Reliability metric    : %s\n", metric_internal))
    cat(sprintf("  M per iteration       : %d\n", M_per_iter))
    cat(sprintf("  M for variance pre-calc: %d\n", M_pre))
    cat(sprintf("  Total iterations      : %d\n", n_iter))
    cat(sprintf("  Burn-in               : %d\n", burn_in))
    cat(sprintf("  Resample items        : %s\n", ifelse(resample_items, "Yes", "No")))
    cat(sprintf("  Step params: a=%.2f, A=%.0f, gamma=%.2f\n",
                step_params$a, step_params$A, step_params$gamma))
    cat("\n")
    cat("Step 0: Pre-calculating latent variance...\n")
  }

  # Generate large sample for stable variance estimation
  latent_args_pre <- modifyList(
    list(n = M_pre, shape = latent_shape),
    latent_params
  )
  theta_pre <- do.call(sim_latentG, latent_args_pre)$theta
  theta_var <- var(theta_pre)

  if (verbose_level >= 1) {
    cat(sprintf("  Estimated theta_var = %.4f (from M_pre = %d samples)\n", theta_var, M_pre))
  }

  # ===========================================================================
  # Step 1: Initialize c_0 (warm start or default)
  # ===========================================================================

  if (verbose_level >= 1) {
    cat("\nStep 1: Initializing c_0...\n")
  }

  # Handle different c_init types
  if (is.null(c_init)) {
    # APC warm start
    c_init_val <- compute_apc_init(target_rho, n_items)
    init_method <- "apc_warm_start"
    if (verbose_level >= 1) {
      cat(sprintf("  c_0 = %.4f (APC warm start)\n", c_init_val))
    }
  } else if (inherits(c_init, "eqc_result")) {
    # EQC object warm start (NEW FEATURE)
    c_init_val <- c_init$c_star
    init_method <- "eqc_warm_start"
    if (verbose_level >= 1) {
      cat(sprintf("  c_0 = %.4f (EQC warm start from eqc_result object)\n", c_init_val))
    }
  } else if (is.numeric(c_init) && length(c_init) == 1L && c_init > 0) {
    # User-specified numeric value
    c_init_val <- c_init
    init_method <- "user_specified"
    if (verbose_level >= 1) {
      cat(sprintf("  c_0 = %.4f (user specified)\n", c_init_val))
    }
  } else {
    stop("c_init must be NULL, a positive scalar, or an eqc_result object.")
  }

  # Clip to bounds
  c_init_val <- max(c_bounds[1], min(c_bounds[2], c_init_val))

  # ===========================================================================
  # Step 2: Generate Fixed Item Parameters (if not resampling)
  # ===========================================================================

  fixed_items <- NULL
  if (!resample_items) {
    if (verbose_level >= 1) cat("Step 2: Generating fixed item parameters...\n")

    item_args <- modifyList(
      list(
        n_items = n_items,
        model   = model,
        source  = item_source,
        scale   = 1
      ),
      item_params
    )
    fixed_items <- do.call(sim_item_params, item_args)
  }

  # ===========================================================================
  # Step 3: Robbins-Monro Iteration
  # ===========================================================================

  if (verbose_level >= 1) {
    cat(sprintf("Step 3: Running %d Robbins-Monro iterations...\n", n_iter))
  }

  c_current <- c_init_val
  trajectory <- numeric(n_iter)
  rho_trajectory <- numeric(n_iter)

  # Progress tracking
  progress_checkpoints <- floor(seq(0.1, 1, by = 0.1) * n_iter)

  for (n in seq_len(n_iter)) {
    # --- Step 3a: Draw fresh samples ---

    # Latent abilities (using M_per_iter, NOT M_pre)
    latent_args <- modifyList(
      list(n = M_per_iter, shape = latent_shape),
      latent_params
    )
    theta_vec <- do.call(sim_latentG, latent_args)$theta

    # Item parameters
    if (resample_items) {
      item_args <- modifyList(
        list(
          n_items = n_items,
          model   = model,
          source  = item_source,
          scale   = 1
        ),
        item_params
      )
      items <- do.call(sim_item_params, item_args)
    } else {
      items <- fixed_items
    }

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

    # --- Step 3b: Compute reliability estimate at c_current ---
    # CRITICAL: Use pre-calculated theta_var (not var(theta_vec))
    rho_hat <- .compute_rho_generic(
      c = c_current,
      theta_vec = theta_vec,
      beta_vec = beta_vec,
      lambda_base = lambda_base,
      theta_var = theta_var,  # FIXED pre-calculated value!
      metric_internal = metric_internal
    )

    # --- Step 3c: Compute step size ---
    a_n <- step_params$a / (n + step_params$A)^step_params$gamma

    # --- Step 3d: Robbins-Monro update ---
    gradient <- rho_hat - target_rho
    c_new <- c_current - a_n * gradient

    # --- Step 3e: Project to feasible region ---
    c_new <- max(c_bounds[1], min(c_bounds[2], c_new))

    # Store trajectory
    trajectory[n] <- c_new
    rho_trajectory[n] <- rho_hat

    # Verbose output
    if (verbose_level >= 2) {
      cat(sprintf("  Iter %4d: c = %.6f, rho = %.4f, grad = %+.4f, a_n = %.4f\n",
                  n, c_new, rho_hat, gradient, a_n))
    } else if (verbose_level >= 1 && n %in% progress_checkpoints) {
      pct <- round(100 * n / n_iter)
      cat(sprintf("  [%3d%%] Iter %4d: c = %.4f, rho = %.4f\n", pct, n, c_new, rho_hat))
    }

    c_current <- c_new
  }

  # ===========================================================================
  # Step 4: Polyak-Ruppert Averaging
  # ===========================================================================

  if (verbose_level >= 1) cat("\nStep 4: Computing Polyak-Ruppert average...\n")

  # Average over post-burn-in iterates
  avg_indices <- (burn_in + 1):n_iter
  c_star <- mean(trajectory[avg_indices])

  # Clip final estimate to bounds
  c_star <- max(c_bounds[1], min(c_bounds[2], c_star))

  if (verbose_level >= 1) {
    cat(sprintf("  Polyak-Ruppert c* = %.6f (averaging %d iterates)\n",
                c_star, length(avg_indices)))
  }

  # ===========================================================================
  # Step 5: Post-Calibration Reliability Estimate
  # ===========================================================================

  if (verbose_level >= 1) cat("\nStep 5: Computing post-calibration reliability...\n")

  # Use larger sample for final estimate
  M_final <- max(M_per_iter * 5, 5000)
  latent_args_final <- modifyList(
    list(n = M_final, shape = latent_shape),
    latent_params
  )
  theta_final <- do.call(sim_latentG, latent_args_final)$theta

  # Generate fresh items for final estimate
  item_args_final <- modifyList(
    list(
      n_items = n_items,
      model   = model,
      source  = item_source,
      scale   = 1
    ),
    item_params
  )
  items_final <- do.call(sim_item_params, item_args_final)
  beta_final <- items_final$data$beta

  if (model == "rasch") {
    lambda_base_final <- rep(1, n_items)
  } else {
    if ("lambda_unscaled" %in% names(items_final$data)) {
      lambda_base_final <- items_final$data$lambda_unscaled
    } else {
      lambda_base_final <- items_final$data$lambda
    }
  }

  # Use pre-calculated theta_var for consistency
  achieved_rho <- .compute_rho_generic(
    c = c_star,
    theta_vec = theta_final,
    beta_vec = beta_final,
    lambda_base = lambda_base_final,
    theta_var = theta_var,
    metric_internal = metric_internal
  )

  # ===========================================================================
  # Step 6: Convergence Diagnostics
  # ===========================================================================

  if (verbose_level >= 1) cat("\nStep 6: Computing convergence diagnostics...\n")

  # Split post-burn-in trajectory for convergence check
  n_post_burn <- length(avg_indices)
  half_point <- floor(n_post_burn / 2)
  first_half <- trajectory[avg_indices[1:half_point]]
  second_half <- trajectory[avg_indices[(half_point + 1):n_post_burn]]

  convergence <- list(
    mean_first_half = mean(first_half),
    mean_second_half = mean(second_half),
    sd_post_burn = sd(trajectory[avg_indices]),
    range_post_burn = range(trajectory[avg_indices]),
    final_gradient = rho_trajectory[n_iter] - target_rho,
    hit_lower_bound = any(trajectory == c_bounds[1]),
    hit_upper_bound = any(trajectory == c_bounds[2]),
    converged = abs(mean(first_half) - mean(second_half)) < 0.05
  )

  # Convergence warning
  if (!convergence$converged) {
    warning(sprintf(
      paste0(
        "SPC may not have fully converged.\n",
        "  First half mean:  %.4f\n",
        "  Second half mean: %.4f\n",
        "  Difference:       %.4f\n",
        "Consider increasing n_iter or adjusting step_params."
      ),
      convergence$mean_first_half,
      convergence$mean_second_half,
      abs(convergence$mean_first_half - convergence$mean_second_half)
    ))
  }

  # Bound hit warnings
  if (convergence$hit_lower_bound) {
    warning("Some iterates hit the lower bound. Consider decreasing c_bounds[1].")
  }
  if (convergence$hit_upper_bound) {
    warning("Some iterates hit the upper bound. Consider increasing c_bounds[2].")
  }

  # ===========================================================================
  # Assemble Result
  # ===========================================================================

  if (verbose_level >= 1) {
    cat("\n")
    cat("================================================================\n")
    cat("  SPC Calibration Complete\n")
    cat("================================================================\n")
    cat(sprintf("  Target reliability   : %.4f\n", target_rho))
    cat(sprintf("  Achieved reliability : %.4f\n", achieved_rho))
    cat(sprintf("  Absolute error       : %.4f\n", abs(achieved_rho - target_rho)))
    cat(sprintf("  Calibrated c*        : %.6f\n", c_star))
    cat(sprintf("  Final iterate c_n    : %.6f\n", c_current))
    cat(sprintf("  Init method          : %s\n", init_method))
    cat(sprintf("  Converged            : %s\n", ifelse(convergence$converged, "Yes", "No")))
    cat("\n")
  }

  res <- list(
    c_star         = c_star,
    c_final        = c_current,
    target_rho     = target_rho,
    achieved_rho   = achieved_rho,
    theta_var      = theta_var,
    trajectory     = trajectory,
    rho_trajectory = rho_trajectory,
    metric         = metric_internal,
    model          = model,
    n_items        = n_items,
    n_iter         = n_iter,
    burn_in        = burn_in,
    M_per_iter     = M_per_iter,
    M_pre          = M_pre,
    step_params    = step_params,
    c_bounds       = c_bounds,
    c_init         = c_init_val,
    init_method    = init_method,
    convergence    = convergence,
    call           = match.call()
  )

  class(res) <- c("spc_result", "list")
  res
}


# =============================================================================
# SAC Alias (for nomenclature consistency)
# =============================================================================

#' @rdname spc_calibrate
#' @export
sac_calibrate <- spc_calibrate


# =============================================================================
# S3 Methods for spc_result
# =============================================================================

#' @export
print.spc_result <- function(x, digits = 4, ...) {
  cat("\n")
  cat("=======================================================\n")
  cat("  Stochastic Approximation Calibration (SPC) Results\n")
  cat("=======================================================\n\n")

  cat("Calibration Summary:\n")
  cat(sprintf("  Model                        : %s\n", toupper(x$model)))
  cat(sprintf("  Target reliability (rho*)    : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved reliability         : %.*f\n", digits, x$achieved_rho))
  cat(sprintf("  Absolute error               : %.2e\n", abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor (c*)          : %.*f\n", digits, x$c_star))
  cat("\n")

  cat("Algorithm Settings:\n")
  cat(sprintf("  Number of items (I)          : %d\n", x$n_items))
  cat(sprintf("  M per iteration              : %d\n", x$M_per_iter))
  cat(sprintf("  M for variance pre-calc      : %d\n", x$M_pre))
  cat(sprintf("  Total iterations             : %d\n", x$n_iter))
  cat(sprintf("  Burn-in                      : %d\n", x$burn_in))
  cat(sprintf("  Reliability metric           : %s\n",
              ifelse(x$metric == "info", "Average-information (tilde)", "MSEM-based (bar/w)")))
  cat(sprintf("  Step params: a=%.2f, A=%.0f, gamma=%.2f\n",
              x$step_params$a, x$step_params$A, x$step_params$gamma))
  cat("\n")

  cat("Convergence Diagnostics:\n")
  cat(sprintf("  Initialization method        : %s\n", x$init_method))
  cat(sprintf("  Initial c_0                  : %.4f\n", x$c_init))
  cat(sprintf("  Final iterate c_n            : %.4f\n", x$c_final))
  cat(sprintf("  Polyak-Ruppert c*            : %.4f\n", x$c_star))
  cat(sprintf("  Pre-calculated theta_var     : %.4f\n", x$theta_var))
  cat(sprintf("  Converged                    : %s\n",
              ifelse(x$convergence$converged, "Yes", "No")))
  cat(sprintf("  Post-burn-in SD              : %.4f\n", x$convergence$sd_post_burn))
  cat(sprintf("  Final gradient (rho - rho*)  : %+.4f\n", x$convergence$final_gradient))
  if (x$convergence$hit_lower_bound) {
    cat("  WARNING: Lower bound was hit during iteration.\n")
  }
  if (x$convergence$hit_upper_bound) {
    cat("  WARNING: Upper bound was hit during iteration.\n")
  }
  cat("\n")

  invisible(x)
}


#' @export
summary.spc_result <- function(object, ...) {
  print(object, ...)
}


#' Plot SPC Convergence Trajectory
#'
#' @description
#' Visualizes the Robbins-Monro iteration trajectory.
#'
#' @param x An \code{spc_result} object.
#' @param type Character. Plot type: \code{"trajectory"}, \code{"rho"},
#'   \code{"c"}, or \code{"both"}.
#' @param ... Additional arguments (unused).
#'
#' @return A ggplot object (if ggplot2 is available) or NULL (base R fallback).
#'
#' @export
plot.spc_result <- function(x, type = c("both", "trajectory", "c", "rho"), ...) {

  type <- match.arg(type)

  # Treat "trajectory" and "c" as synonyms
  if (type == "trajectory") type <- "c"

  n_iter <- x$n_iter
  burn_in <- x$burn_in

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    if (type == "both") {
      par(mfrow = c(2, 1))
    }

    if (type %in% c("both", "c")) {
      # Trajectory plot
      plot(1:n_iter, x$trajectory, type = "l", col = "steelblue",
           xlab = "Iteration", ylab = "c",
           main = "SPC Convergence Trajectory")
      abline(h = x$c_star, col = "red", lty = 2, lwd = 2)
      abline(v = burn_in, col = "gray", lty = 3)
      legend("topright", legend = c("Trajectory", "c* (average)", "Burn-in"),
             col = c("steelblue", "red", "gray"), lty = c(1, 2, 3), lwd = c(1, 2, 1))
    }

    if (type %in% c("both", "rho")) {
      # Reliability plot
      plot(1:n_iter, x$rho_trajectory, type = "l", col = "darkorange",
           xlab = "Iteration", ylab = expression(hat(rho)),
           main = "Reliability Estimates")
      abline(h = x$target_rho, col = "red", lty = 2, lwd = 2)
      abline(v = burn_in, col = "gray", lty = 3)
    }

    if (type == "both") {
      par(mfrow = c(1, 1))
    }
    return(invisible(NULL))
  }

  # ggplot2 implementation
  library(ggplot2)

  df <- data.frame(
    iteration = 1:n_iter,
    c = x$trajectory,
    rho = x$rho_trajectory,
    phase = factor(ifelse(1:n_iter <= burn_in, "Burn-in", "Averaging"),
                   levels = c("Burn-in", "Averaging"))
  )

  # Trajectory plot
  p1 <- ggplot(df, aes(x = iteration, y = c, color = phase)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = x$c_star, color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = burn_in, color = "gray50", linetype = "dotted") +
    scale_color_manual(values = c("Burn-in" = "steelblue", "Averaging" = "darkblue")) +
    labs(
      title = "SPC Convergence Trajectory",
      subtitle = sprintf("Target rho* = %.3f | c* = %.4f | Init: %s",
                         x$target_rho, x$c_star, x$init_method),
      x = "Iteration", y = "Scaling Factor c",
      color = "Phase"
    ) +
    annotate("text", x = n_iter * 0.95, y = x$c_star,
             label = sprintf("c* = %.3f", x$c_star),
             vjust = -0.5, color = "red", size = 3.5) +
    theme_minimal() +
    theme(legend.position = "bottom")

  # Reliability plot
  p2 <- ggplot(df, aes(x = iteration, y = rho, color = phase)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = x$target_rho, color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = burn_in, color = "gray50", linetype = "dotted") +
    scale_color_manual(values = c("Burn-in" = "darkorange", "Averaging" = "darkred")) +
    labs(
      title = "Reliability Estimates Across Iterations",
      subtitle = sprintf("Target rho* = %.3f | Achieved = %.3f",
                         x$target_rho, x$achieved_rho),
      x = "Iteration", y = expression(hat(rho)),
      color = "Phase"
    ) +
    annotate("text", x = n_iter * 0.95, y = x$target_rho,
             label = sprintf("rho* = %.3f", x$target_rho),
             vjust = -0.5, color = "red", size = 3.5) +
    theme_minimal() +
    theme(legend.position = "bottom")

  if (type == "c") {
    return(p1)
  } else if (type == "rho") {
    return(p2)
  } else {
    # Combine with patchwork if available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      library(patchwork)
      return(p1 / p2)
    } else {
      print(p1)
      print(p2)
      return(invisible(list(p1 = p1, p2 = p2)))
    }
  }
}
