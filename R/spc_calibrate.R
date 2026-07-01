# =============================================================================
# spc_calibrate.R
# =============================================================================
# Stochastic Approximation Calibration (Algorithm 2: SAC)
#
# Contents:
#   - sac_calibrate(): Main calibration function (primary)
#   - spc_calibrate(): Deprecated alias for sac_calibrate
#   - print.sac_result(): Print method
#   - summary.sac_result(): Summary method
#   - plot.sac_result(): Plot method
#
# Dependencies:
#   - sim_latentG() from sim_latentG.R
#   - sim_item_params() from sim_item_params.R
#   - .compute_rho_generic() from reliability_utils.R
#   - compute_apc_init() from reliability_utils.R
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# =============================================================================


#' Stochastic Approximation Calibration (Algorithm 2: SAC)
#'
#' @description
#' `sac_calibrate()` implements Algorithm 2 (Stochastic Approximation
#' Calibration, SAC) for reliability-targeted IRT simulation using the
#' Robbins-Monro stochastic approximation algorithm.
#'
#' Given a target marginal reliability \eqn{\rho^*}, a latent distribution
#' generator `sim_latentG()` (for \eqn{G}) and an item parameter generator
#' `sim_item_params()` (for \eqn{H}), the function iteratively searches for a
#' global discrimination scale \eqn{c^* > 0} such that the population
#' reliability \eqn{\rho(c)} of the Rasch/2PL model is approximately equal to
#' \eqn{\rho^*}.
#'
#' SAC complements EQC (Algorithm 1) by:
#' \enumerate{
#'   \item Providing an independent validation of EQC calibration results.
#'   \item Enabling calibration to the MSEM-based reciprocal-information
#'     reliability \eqn{\bar{w}} (not just the average-information metric
#'     \eqn{\tilde{\rho}}).
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
#'   (e.g. \code{"parametric"}, \code{"irw"}, \code{"hierarchical"}, \code{"custom"}).
#'   Defaults to \code{"parametric"} so core workflows do not require external
#'   item-pool data. Use \code{"irw"} for empirically-grounded difficulties
#'   when the Item Response Warehouse package is installed.
#'
#' @param latent_params List. Additional arguments passed to \code{sim_latentG()}.
#'
#' @param item_params List. Additional arguments passed to \code{sim_item_params()}.
#'
#' @param reliability_metric Character. Reliability definition used inside SAC:
#'   \describe{
#'     \item{\code{"msem"}}{MSEM-based reciprocal-information reliability (targets \eqn{\bar{w}}).}
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
#' @return An object of class \code{"sac_result"} (a list) with elements:
#' \describe{
#'   \item{\code{c_star}}{Calibrated discrimination scale (Polyak-Ruppert average).}
#'   \item{\code{c_final}}{Final iterate \eqn{c_{n_{iter}}}.}
#'   \item{\code{target_rho}}{Target reliability \eqn{\rho^*}.}
#'   \item{\code{achieved_rho}}{Estimated reliability at \eqn{c^*} (post-calibration).}
#'   \item{\code{theta_var}}{Pre-calculated latent variance used throughout.}
#'   \item{\code{trajectory}}{Numeric vector of all iterates.}
#'   \item{\code{rho_trajectory}}{Numeric vector of reliability estimates.}
#'   \item{\code{raw_trajectory}}{Unprojected proposal values before clipping to
#'     \code{c_bounds}.}
#'   \item{\code{step_size_trajectory}}{Robbins-Monro step size \eqn{a_n} at
#'     each iteration.}
#'   \item{\code{gradient_trajectory}}{Stochastic gradient estimate
#'     \eqn{\hat{\rho}_n - \rho^*} at each iteration.}
#'   \item{\code{projected}}{Logical vector indicating whether projection was
#'     applied at each iteration.}
#'   \item{\code{projection_side}}{Character vector indicating lower/upper/no
#'     projection at each iteration.}
#'   \item{\code{projection_count}}{Number of iterations where projection onto
#'     \code{c_bounds} was applied.}
#'   \item{\code{projection_rate}}{Proportion of iterations where projection was
#'     applied.}
#'   \item{\code{M_final}}{Monte Carlo sample size used for the post-calibration
#'     reliability estimate.}
#'   \item{\code{metric}}{Reliability metric used.}
#'   \item{\code{calibration_status}}{Canonical status label copied from
#'     \code{convergence$status}.}
#'   \item{\code{model}}{Measurement model.}
#'   \item{\code{n_items}}{Number of items.}
#'   \item{\code{n_iter}}{Total number of Robbins-Monro iterations.}
#'   \item{\code{burn_in}}{Number of initial iterations excluded from averaging.}
#'   \item{\code{M_per_iter}}{Monte Carlo sample size per iteration.}
#'   \item{\code{M_pre}}{Monte Carlo sample size used to estimate latent variance.}
#'   \item{\code{step_params}}{Step-size parameters used by the run.}
#'   \item{\code{c_bounds}}{Projection bounds for \eqn{c}.}
#'   \item{\code{c_init}}{Initial scaling factor after any projection.}
#'   \item{\code{init_method}}{Character indicating initialization method.}
#'   \item{\code{convergence}}{List with convergence diagnostics.}
#'   \item{\code{item_design}}{Whether stored item fields are a
#'     \code{"post_calibration_draw"} or the \code{"fixed_iteration_items"}.}
#'   \item{\code{beta_vec}}{Item difficulties from the stored item design.}
#'   \item{\code{lambda_base}}{Baseline (unscaled) item discriminations.}
#'   \item{\code{lambda_scaled}}{Scaled item discriminations (\code{lambda_base * c_star}).}
#'   \item{\code{items_base}}{Baseline item_params object (scale = 1).}
#'   \item{\code{items_calib}}{Calibrated item_params object (discriminations scaled by c_star).}
#'   \item{\code{theta_quad}}{Theta sample used for post-calibration reliability estimate.}
#'   \item{\code{call}}{Matched function call.}
#' }
#'
#' The stored item fields use the final post-calibration item draw at baseline
#' scale 1 when \code{resample_items = TRUE}. When \code{resample_items = FALSE},
#' they reuse the fixed item form used throughout the Robbins-Monro iterations.
#' In both cases, the calibrated design applies the Polyak-Ruppert average
#' \code{c_star}, not the last Robbins-Monro iterate \code{c_final}. Thus
#' \code{lambda_scaled} and \code{items_calib$data$lambda} equal
#' \code{lambda_base * c_star}.
#'
#' For same-estimand comparisons with an \code{eqc_result} warm start, set
#' \code{reliability_metric = "info"} in \code{sac_calibrate()}. SAC defaults
#' to direct MSEM targeting (\code{"msem"}), while EQC targets the
#' average-information metric (\code{"info"}).
#'
#' \code{spc_calibrate()} is a deprecated backward-compatible alias for
#' \code{sac_calibrate()}.
#'
#' @seealso
#' \code{\link{eqc_calibrate}} for the faster deterministic Algorithm 1,
#' \code{\link{compute_rho_bar}} and \code{\link{compute_rho_tilde}} for
#' reliability computation utilities,
#' \code{\link{compare_eqc_sac}} for comparing EQC and SAC results.
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
#' \donttest{
#' # Example 1: Basic SAC calibration
#' sac_result <- sac_calibrate(
#'   target_rho = 0.75,
#'   n_items = 20,
#'   model = "rasch",
#'   n_iter = 200,
#'   seed = 12345,
#'   verbose = TRUE
#' )
#' print(sac_result)
#' plot(sac_result)
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
#' sac_result <- sac_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "2pl",
#'   reliability_metric = "info",
#'   c_init = eqc_result,  # Direct EQC object passing!
#'   n_iter = 100,
#'   seed = 42
#' )
#'
#' # Compare EQC and SAC results
#' cat(sprintf("EQC c* = %.4f, SAC c* = %.4f\n",
#'             eqc_result$c_star, sac_result$c_star))
#' }
#'
#' @export
sac_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl"),
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
                          verbose = FALSE) {

  # ===========================================================================
  # Input Validation
  # ===========================================================================

  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      is.na(target_rho) || !is.finite(target_rho) ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a scalar in (0, 1).")
  }

  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")

  model <- match.arg(model)

  reliability_metric <- match.arg(reliability_metric)
  metric_internal <- if (reliability_metric %in% c("msem", "bar")) "msem" else "info"

  M_per_iter <- .irtsimrel_validate_positive_integer_scalar(
    M_per_iter,
    "M_per_iter"
  )

  M_pre <- .irtsimrel_validate_positive_integer_scalar(M_pre, "M_pre")
  if (M_pre < 2L) {
    stop("`M_pre` must be at least 2 to estimate latent variance.")
  }

  n_iter <- .irtsimrel_validate_positive_integer_scalar(n_iter, "n_iter")

  # Set default burn-in
  if (is.null(burn_in)) {
    burn_in <- floor(n_iter / 2)
  }
  if (!is.numeric(burn_in) || length(burn_in) != 1L ||
      is.na(burn_in) || !is.finite(burn_in) ||
      burn_in < 0 || burn_in >= n_iter || burn_in %% 1 != 0) {
    stop("`burn_in` must be a non-negative integer less than `n_iter`.")
  }
  burn_in <- as.integer(burn_in)

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(!is.finite(c_bounds)) || any(c_bounds <= 0) ||
      c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  if (!is.logical(resample_items) || length(resample_items) != 1L ||
      is.na(resample_items)) {
    stop("`resample_items` must be TRUE or FALSE.")
  }

  if (inherits(c_init, "eqc_result")) {
    .irtsimrel_validate_eqc_result_object(c_init, "c_init", "SAC warm start")
  }

  if (!(is.logical(verbose) || is.numeric(verbose)) ||
      length(verbose) != 1L || is.na(verbose) ||
      (is.numeric(verbose) && !is.finite(verbose))) {
    stop("`verbose` must be a logical or non-negative numeric scalar.")
  }
  verbose_level <- as.integer(verbose)
  if (verbose_level < 0) {
    stop("`verbose` must be a logical or non-negative numeric scalar.")
  }

  # Set default step parameters
  default_step <- list(a = 1.0, A = 50, gamma = 0.67)
  if (is.null(step_params)) {
    step_params <- list()
  }
  if (!is.list(step_params)) {
    stop("`step_params` must be a list.")
  }
  if (length(step_params) > 0L) {
    if (is.null(names(step_params)) || any(!nzchar(names(step_params)))) {
      stop("All entries in `step_params` must be named.")
    }
    unknown_step_params <- setdiff(names(step_params), names(default_step))
    if (length(unknown_step_params) > 0L) {
      stop(
        "Unknown `step_params` field(s): ",
        paste(unknown_step_params, collapse = ", "),
        ". Allowed fields are: a, A, gamma."
      )
    }
  }
  step_params <- modifyList(default_step, step_params)

  validate_step_scalar <- function(name) {
    value <- step_params[[name]]
    if (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value)) {
      stop("`step_params$", name, "` must be a finite numeric scalar.")
    }
    value
  }
  step_params$a <- validate_step_scalar("a")
  step_params$A <- validate_step_scalar("A")
  step_params$gamma <- validate_step_scalar("gamma")

  if (step_params$a <= 0) stop("`step_params$a` must be positive.")
  if (step_params$A < 0) stop("`step_params$A` must be non-negative.")
  if (step_params$gamma <= 0.5 || step_params$gamma > 1) {
    stop(
      "`step_params$gamma` must be in (0.5, 1] to satisfy the ",
      "Robbins-Monro step-size conditions."
    )
  }

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)
  item_params <- .irtsimrel_normalize_item_params(item_params)

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  # ===========================================================================
  # Step 0: Pre-Calculate Latent Variance (CRITICAL for stability)
  # ===========================================================================

  if (verbose_level >= 1) {
    message("")
    message("================================================================")
    message("  Stochastic Approximation Calibration (SAC)")
    message("================================================================")
    message("")
    message("Configuration:")
    message(sprintf("  Target reliability    : %.4f", target_rho))
    message(sprintf("  Number of items       : %d", n_items))
    message(sprintf("  Model                 : %s", toupper(model)))
    message(sprintf("  Reliability metric    : %s", metric_internal))
    message(sprintf("  M per iteration       : %d", M_per_iter))
    message(sprintf("  M for variance pre-calc: %d", M_pre))
    message(sprintf("  Total iterations      : %d", n_iter))
    message(sprintf("  Burn-in               : %d", burn_in))
    message(sprintf("  Resample items        : %s", ifelse(resample_items, "Yes", "No")))
    message(sprintf("  Step params: a=%.2f, A=%.0f, gamma=%.2f",
                step_params$a, step_params$A, step_params$gamma))
    message("")
    message("Step 0: Pre-calculating latent variance...")
  }

  # Generate large sample for stable variance estimation
  latent_args_pre <- modifyList(
    list(n = M_pre, shape = latent_shape),
    latent_params
  )
  theta_pre <- do.call(sim_latentG, latent_args_pre)$theta
  theta_var <- var(theta_pre)
  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 1e-10) {
    stop("Pre-calculated latent variance must be finite and positive.")
  }

  if (verbose_level >= 1) {
    message(sprintf("  Estimated theta_var = %.4f (from M_pre = %d samples)", theta_var, M_pre))
  }

  # ===========================================================================
  # Step 1: Initialize c_0 (warm start or default)
  # ===========================================================================

  if (verbose_level >= 1) {
    message("\nStep 1: Initializing c_0...")
  }

  # Handle different c_init types
  if (is.null(c_init)) {
    # APC warm start
    c_init_val <- compute_apc_init(target_rho, n_items)
    init_method <- "apc_warm_start"
    if (verbose_level >= 1) {
      message(sprintf("  c_0 = %.4f (APC warm start)", c_init_val))
    }
  } else if (inherits(c_init, "eqc_result")) {
    .irtsimrel_validate_eqc_result_object(c_init, "c_init", "SAC warm start")
    # EQC object warm start (NEW FEATURE)
    c_init_val <- c_init$c_star
    init_method <- "eqc_warm_start"
    if (verbose_level >= 1) {
      message(sprintf("  c_0 = %.4f (EQC warm start from eqc_result object)", c_init_val))
    }
  } else if (is.numeric(c_init) && length(c_init) == 1L &&
             is.finite(c_init) && c_init > 0) {
    # User-specified numeric value
    c_init_val <- c_init
    init_method <- "user_specified"
    if (verbose_level >= 1) {
      message(sprintf("  c_0 = %.4f (user specified)", c_init_val))
    }
  } else {
    stop("c_init must be NULL, a positive scalar, or an eqc_result object.")
  }

  # Clip to bounds
  c_init_raw <- c_init_val
  c_init_val <- max(c_bounds[1], min(c_bounds[2], c_init_val))
  c_init_projected <- !isTRUE(all.equal(c_init_raw, c_init_val, tolerance = 0))

  if (verbose_level >= 1 && c_init_projected) {
    message(sprintf(
      "  c_0 projected from %.4f to %.4f to satisfy c_bounds.",
      c_init_raw, c_init_val
    ))
  }

  # ===========================================================================
  # Step 2: Generate Fixed Item Parameters (if not resampling)
  # ===========================================================================

  fixed_items <- NULL
  if (!resample_items) {
    if (verbose_level >= 1) message("Step 2: Generating fixed item parameters...")

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
    message(sprintf("Step 3: Running %d Robbins-Monro iterations...", n_iter))
  }

  c_current <- c_init_val
  trajectory <- numeric(n_iter)
  rho_trajectory <- numeric(n_iter)
  raw_trajectory <- numeric(n_iter)
  step_size_trajectory <- numeric(n_iter)
  gradient_trajectory <- numeric(n_iter)
  projected <- logical(n_iter)
  projection_side <- rep("none", n_iter)

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
    c_raw <- c_current - a_n * gradient

    # --- Step 3e: Project to feasible region ---
    c_new <- max(c_bounds[1], min(c_bounds[2], c_raw))
    projected[n] <- !isTRUE(all.equal(c_raw, c_new, tolerance = 0))
    if (projected[n]) {
      projection_side[n] <- if (c_raw < c_bounds[1]) "lower" else "upper"
    }

    # Store trajectory
    trajectory[n] <- c_new
    rho_trajectory[n] <- rho_hat
    raw_trajectory[n] <- c_raw
    step_size_trajectory[n] <- a_n
    gradient_trajectory[n] <- gradient

    # Verbose output
    if (verbose_level >= 2) {
      message(sprintf("  Iter %4d: c = %.6f, rho = %.4f, grad = %+.4f, a_n = %.4f",
                  n, c_new, rho_hat, gradient, a_n))
    } else if (verbose_level >= 1 && n %in% progress_checkpoints) {
      pct <- round(100 * n / n_iter)
      message(sprintf("  [%3d%%] Iter %4d: c = %.4f, rho = %.4f", pct, n, c_new, rho_hat))
    }

    c_current <- c_new
  }

  # ===========================================================================
  # Step 4: Polyak-Ruppert Averaging
  # ===========================================================================

  if (verbose_level >= 1) message("\nStep 4: Computing Polyak-Ruppert average...")

  # Average over post-burn-in iterates
  avg_indices <- (burn_in + 1):n_iter
  c_star <- mean(trajectory[avg_indices])

  # Clip final estimate to bounds
  c_star <- max(c_bounds[1], min(c_bounds[2], c_star))

  if (verbose_level >= 1) {
    message(sprintf("  Polyak-Ruppert c* = %.6f (averaging %d iterates)",
                c_star, length(avg_indices)))
  }

  # ===========================================================================
  # Step 5: Post-Calibration Reliability Estimate
  # ===========================================================================

  if (verbose_level >= 1) message("\nStep 5: Computing post-calibration reliability...")

  # Use larger sample for final estimate
  M_final <- max(M_per_iter * 5, 5000)
  latent_args_final <- modifyList(
    list(n = M_final, shape = latent_shape),
    latent_params
  )
  theta_final <- do.call(sim_latentG, latent_args_final)$theta

  # Generate a representative final form only when the algorithm resampled
  # items. For fixed-item SAC, the result must describe the calibrated fixed
  # form used during the Robbins-Monro iterations.
  if (resample_items) {
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
    item_design <- "post_calibration_draw"
  } else {
    items_final <- fixed_items
    item_design <- "fixed_iteration_items"
  }
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
  # Step 5b: Store Calibrated Item Parameters
  # ===========================================================================

  items_calib_final <- .irtsimrel_apply_item_scale(
    items_final,
    lambda_base_final,
    c_star
  )

  # ===========================================================================
  # Step 6: Convergence Diagnostics
  # ===========================================================================

  if (verbose_level >= 1) message("\nStep 6: Computing convergence diagnostics...")

  n_post_burn <- length(avg_indices)
  post_burn_trajectory <- trajectory[avg_indices]
  min_post_burn <- 4L
  insufficient_post_burn <- n_post_burn < min_post_burn
  # Split-half drift is measured on the c scale, so use a relative tolerance.
  split_mean_rel_tolerance <- 0.05
  split_mean_scale <- max(abs(c_star), .Machine$double.eps)
  split_mean_tolerance <- split_mean_rel_tolerance * split_mean_scale

  if (insufficient_post_burn) {
    first_half <- post_burn_trajectory
    second_half <- post_burn_trajectory
    split_mean_diff <- NA_real_
    split_mean_rel_diff <- NA_real_
    converged <- FALSE
  } else {
    half_point <- floor(n_post_burn / 2)
    first_half <- trajectory[avg_indices[seq_len(half_point)]]
    second_half <- trajectory[avg_indices[(half_point + 1):n_post_burn]]
    split_mean_diff <- abs(mean(first_half) - mean(second_half))
    split_mean_rel_diff <- split_mean_diff / split_mean_scale
    converged <- split_mean_rel_diff < split_mean_rel_tolerance
  }

  final_iter_gradient <- rho_trajectory[n_iter] - target_rho
  post_calibration_gradient <- achieved_rho - target_rho
  achieved_gap_abs <- abs(post_calibration_gradient)
  achieved_gap_tolerance <- 0.05
  large_achieved_gap <- achieved_gap_abs > achieved_gap_tolerance
  trajectory_converged <- converged
  converged <- trajectory_converged && !large_achieved_gap
  projection_count <- sum(projected)
  projection_rate <- projection_count / n_iter
  hit_lower_bound <- any(trajectory == c_bounds[1]) ||
    (isTRUE(c_init_projected) && c_init_val == c_bounds[1])
  hit_upper_bound <- any(trajectory == c_bounds[2]) ||
    (isTRUE(c_init_projected) && c_init_val == c_bounds[2])

  convergence <- list(
    mean_first_half = mean(first_half),
    mean_second_half = mean(second_half),
    split_mean_diff = split_mean_diff,
    split_mean_rel_diff = split_mean_rel_diff,
    split_mean_scale = split_mean_scale,
    split_mean_tolerance = split_mean_tolerance,
    split_mean_rel_tolerance = split_mean_rel_tolerance,
    sd_post_burn = if (n_post_burn > 1L) stats::sd(post_burn_trajectory) else 0,
    range_post_burn = range(post_burn_trajectory),
    n_post_burn = n_post_burn,
    min_post_burn = min_post_burn,
    insufficient_post_burn = insufficient_post_burn,
    trajectory_converged = trajectory_converged,
    achieved_gap = post_calibration_gradient,
    achieved_gap_abs = achieved_gap_abs,
    achieved_gap_tolerance = achieved_gap_tolerance,
    large_achieved_gap = large_achieved_gap,
    final_gradient = final_iter_gradient,
    final_pre_update_gradient = final_iter_gradient,
    final_iter_gradient = final_iter_gradient,
    gradient_at_c_star = post_calibration_gradient,
    post_calibration_gradient = post_calibration_gradient,
    projection_count = projection_count,
    projection_rate = projection_rate,
    c_init_raw = c_init_raw,
    c_init_projected = c_init_projected,
    n_lower_hits = sum(projection_side == "lower"),
    n_upper_hits = sum(projection_side == "upper"),
    first_projection_iter = if (projection_count > 0L) which(projected)[1] else NA_integer_,
    post_burn_projection_rate = mean(projected[avg_indices]),
    hit_lower_bound = hit_lower_bound,
    hit_upper_bound = hit_upper_bound,
    converged = converged
  )

  status_flags <- character()
  if (!convergence$converged) {
    status_flags <- c(status_flags, "not_converged")
  }
  if (convergence$insufficient_post_burn) {
    status_flags <- c(status_flags, "insufficient_post_burn")
  }
  if (convergence$large_achieved_gap) {
    status_flags <- c(status_flags, "large_achieved_gap")
  }
  if (convergence$hit_lower_bound) {
    status_flags <- c(status_flags, "hit_lower_bound")
  }
  if (convergence$hit_upper_bound) {
    status_flags <- c(status_flags, "hit_upper_bound")
  }
  if (projection_count > 0L || c_init_projected) {
    status_flags <- c(status_flags, "projection_applied")
  }
  if (length(status_flags) == 0L) {
    status_flags <- "ok"
  }

  calibration_status <- if (!convergence$converged) {
    "not_converged"
  } else if (convergence$hit_lower_bound && convergence$hit_upper_bound) {
    "hit_both_bounds"
  } else if (convergence$hit_lower_bound) {
    "hit_lower_bound"
  } else if (convergence$hit_upper_bound) {
    "hit_upper_bound"
  } else {
    "ok"
  }
  convergence$status <- calibration_status
  convergence$status_flags <- status_flags

  # Convergence warning
  if (!convergence$converged) {
    warning(sprintf(
      paste0(
        "SAC may not have fully converged.\n",
        "  First half mean:  %.4f\n",
        "  Second half mean: %.4f\n",
        "  Difference:       %s (tolerance %.4f)\n",
        "  Achieved gap:     %+.4f (tolerance %.4f)\n",
        "Consider increasing n_iter or adjusting step_params."
      ),
      convergence$mean_first_half,
      convergence$mean_second_half,
      if (is.na(convergence$split_mean_diff)) {
        "NA (insufficient post-burn-in iterates)"
      } else {
        sprintf("%.4f", convergence$split_mean_diff)
      },
      convergence$split_mean_tolerance,
      convergence$achieved_gap,
      convergence$achieved_gap_tolerance
    ))
  }

  # Bound hit warnings
  if (convergence$hit_lower_bound) {
    warning(
      "SAC projection hit the lower bound. ",
      "Lower projections: ", convergence$n_lower_hits,
      "; post-burn projection rate: ",
      sprintf("%.1f%%", 100 * convergence$post_burn_projection_rate),
      "; gradient at c*: ", sprintf("%+.4f", convergence$gradient_at_c_star),
      ". Consider decreasing c_bounds[1], using a smaller step_params$a, ",
      "or checking whether the target reliability is too low for this design."
    )
  }
  if (convergence$hit_upper_bound) {
    warning(
      "SAC projection hit the upper bound. ",
      "Upper projections: ", convergence$n_upper_hits,
      "; post-burn projection rate: ",
      sprintf("%.1f%%", 100 * convergence$post_burn_projection_rate),
      "; gradient at c*: ", sprintf("%+.4f", convergence$gradient_at_c_star),
      ". Consider increasing c_bounds[2], using a smaller step_params$a, ",
      "or checking whether the target reliability is too high for this design."
    )
  }

  # ===========================================================================
  # Assemble Result
  # ===========================================================================

  if (verbose_level >= 1) {
    message("")
    message("================================================================")
    message("  SAC Calibration Complete")
    message("================================================================")
    message(sprintf("  Target reliability   : %.4f", target_rho))
    message(sprintf("  Achieved reliability : %.4f", achieved_rho))
    message(sprintf("  Absolute error       : %.4f", abs(achieved_rho - target_rho)))
    message(sprintf("  Calibrated c*        : %.6f", c_star))
    message(sprintf("  Final iterate c_n    : %.6f", c_current))
    message(sprintf("  Init method          : %s", init_method))
    message(sprintf("  Converged            : %s", ifelse(convergence$converged, "Yes", "No")))
    message("")
  }

  res <- list(
    c_star             = c_star,
    c_final            = c_current,
    target_rho         = target_rho,
    achieved_rho       = achieved_rho,
    theta_var          = theta_var,
    trajectory         = trajectory,
    rho_trajectory     = rho_trajectory,
    raw_trajectory     = raw_trajectory,
    step_size_trajectory = step_size_trajectory,
    gradient_trajectory = gradient_trajectory,
    projected          = projected,
    projection_side    = projection_side,
    projection_count   = projection_count,
    projection_rate    = projection_rate,
    M_final            = M_final,
    metric             = metric_internal,
    calibration_status = calibration_status,
    model              = model,
    n_items            = n_items,
    n_iter             = n_iter,
    burn_in            = burn_in,
    M_per_iter         = M_per_iter,
    M_pre              = M_pre,
    step_params        = step_params,
    c_bounds           = c_bounds,
    c_init             = c_init_val,
    init_method        = init_method,
    item_design        = item_design,
    convergence        = convergence,
    beta_vec           = beta_final,
    lambda_base        = lambda_base_final,
    lambda_scaled      = lambda_base_final * c_star,
    items_base         = items_final,
    items_calib        = items_calib_final,
    theta_quad         = theta_final,
    call               = match.call()
  )

  class(res) <- c("sac_result", "list")
  res
}


# =============================================================================
# SPC Deprecated Alias
# =============================================================================

#' Deprecated Alias for SAC Calibration
#'
#' @description
#' `spc_calibrate()` is a deprecated backward-compatible alias for
#' `sac_calibrate()`. It is retained for scripts written against the v0.1.x
#' naming convention and will issue a deprecation warning when called.
#'
#' @param ... Arguments passed to \code{sac_calibrate()}.
#'
#' @return The result of \code{\link{sac_calibrate}}.
#'
#' @seealso \code{\link{sac_calibrate}}
#'
#' @export
spc_calibrate <- function(...) {
  .Deprecated("sac_calibrate")
  sac_calibrate(...)
}


# =============================================================================
# S3 Methods for sac_result
# =============================================================================

#' @rdname sac_calibrate
#' @param x An object of class \code{"sac_result"}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments passed to or from other methods.
#' @return The input object, invisibly.
#' @export
print.sac_result <- function(x, digits = 4, ...) {
  .irtsimrel_validate_sac_result_object(x, "x", "SAC print")

  cat("\n")
  cat("=======================================================\n")
  cat("  Stochastic Approximation Calibration (SAC) Results\n")
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
  cat(sprintf("  Final iter gradient          : %+.4f\n", x$convergence$final_gradient))
  if (!is.null(x$convergence$gradient_at_c_star)) {
    cat(sprintf("  Gradient at c*               : %+.4f\n",
                x$convergence$gradient_at_c_star))
  }
  if (!is.null(x$convergence$projection_count)) {
    cat(sprintf("  Projection count             : %d (%.1f%%)\n",
                x$convergence$projection_count,
                100 * x$convergence$projection_rate))
  }
  if (!is.null(x$convergence$status_flags)) {
    cat(sprintf("  Status flags                 : %s\n",
                paste(x$convergence$status_flags, collapse = ", ")))
  }
  if (x$convergence$hit_lower_bound) {
    cat("  WARNING: Lower bound was hit during iteration.\n")
  }
  if (x$convergence$hit_upper_bound) {
    cat("  WARNING: Upper bound was hit during iteration.\n")
  }
  cat("\n")

  invisible(x)
}


#' @rdname sac_calibrate
#' @param object An object of class \code{"sac_result"}.
#' @return An object of class \code{"summary.sac_result"} containing key
#'   calibration results.
#' @export
summary.sac_result <- function(object, ...) {
  .irtsimrel_validate_sac_result_object(object, "object", "SAC summary")

  out <- list(
    c_star       = object$c_star,
    target_rho   = object$target_rho,
    achieved_rho = object$achieved_rho,
    metric       = object$metric,
    model        = object$model,
    n_items      = object$n_items,
    n_iter       = object$n_iter,
    M_per_iter   = object$M_per_iter,
    M_pre        = object$M_pre,
    theta_var    = object$theta_var,
    init_method  = object$init_method,
    convergence  = object$convergence,
    burn_in      = object$burn_in
  )
  class(out) <- "summary.sac_result"
  out
}


#' Print Method for summary.sac_result Objects
#'
#' @param x A \code{summary.sac_result} object from \code{summary.sac_result()}.
#' @param digits Integer. Number of decimal places for printing.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @export
print.summary.sac_result <- function(x, digits = 4, ...) {
  .irtsimrel_validate_summary_object(
    x,
    "summary.sac_result",
    c(
      "c_star", "target_rho", "achieved_rho", "metric", "model",
      "n_items", "n_iter", "M_per_iter", "M_pre", "theta_var",
      "init_method", "convergence", "burn_in"
    ),
    "x",
    "SAC summary print"
  )

  cat("Summary: Stochastic Approximation Calibration (SAC)\n")
  cat("====================================================\n")
  cat(sprintf("  Model            : %s\n", toupper(x$model)))
  cat(sprintf("  Metric           : %s\n",
              ifelse(x$metric == "info", "Average-information (tilde)", "MSEM-based (bar/w)")))
  cat(sprintf("  Number of items  : %d\n", x$n_items))
  cat(sprintf("  Iterations       : %d\n", x$n_iter))
  cat(sprintf("  Burn-in          : %d\n", x$burn_in))
  cat(sprintf("  M per iteration  : %d\n", x$M_per_iter))
  cat(sprintf("  M pre-calc       : %d\n", x$M_pre))
  cat(sprintf("  Latent variance  : %.*f\n", digits, x$theta_var))
  cat(sprintf("  Init method      : %s\n", x$init_method))
  cat("\nCalibration Results:\n")
  cat(sprintf("  Target rho*      : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved rho     : %.*f\n", digits, x$achieved_rho))
  cat(sprintf("  Absolute error   : %.2e\n", abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor c*: %.*f\n", digits, x$c_star))
  cat(sprintf("  Converged        : %s\n",
              ifelse(x$convergence$converged, "Yes", "No")))
  cat(sprintf("  Post-burn-in SD  : %.*f\n", digits, x$convergence$sd_post_burn))
  invisible(x)
}


#' Plot SAC Convergence Trajectory
#'
#' @description
#' Visualizes the Robbins-Monro iteration trajectory.
#'
#' @param x An \code{sac_result} object from \code{\link{sac_calibrate}}.
#' @param type Character. Plot type: \code{"trajectory"}, \code{"rho"},
#'   \code{"c"}, or \code{"both"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A \code{ggplot} object if \pkg{ggplot2} is available for
#'   \code{type = "c"} or \code{type = "rho"}. For \code{type = "both"},
#'   returns a patchwork object when \pkg{patchwork} is available; otherwise
#'   prints both ggplots and returns them invisibly as a list. Returns
#'   \code{NULL} invisibly when using base R graphics fallback.
#'
#' @export
plot.sac_result <- function(x, type = c("both", "trajectory", "c", "rho"), ...) {
  .irtsimrel_validate_sac_result_object(x, "x", "SAC plot")

  type <- match.arg(type)

  # Treat "trajectory" and "c" as synonyms
  if (type == "trajectory") type <- "c"

  n_iter <- x$n_iter
  burn_in <- x$burn_in

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    # Base R fallback
    if (type == "both") {
      oldpar <- par(mfrow = c(2, 1))
      on.exit(par(oldpar))
    }

    if (type %in% c("both", "c")) {
      # Trajectory plot
      plot(1:n_iter, x$trajectory, type = "l", col = "steelblue",
           xlab = "Iteration", ylab = "c",
           main = "SAC Convergence Trajectory")
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

    return(invisible(NULL))
  }

  # ggplot2 implementation
  df <- data.frame(
    iteration = 1:n_iter,
    c = x$trajectory,
    rho = x$rho_trajectory,
    phase = factor(ifelse(1:n_iter <= burn_in, "Burn-in", "Averaging"),
                   levels = c("Burn-in", "Averaging"))
  )

  # Trajectory plot
  p1 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = c, color = phase)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = x$c_star, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = burn_in, color = "gray50", linetype = "dotted") +
    ggplot2::scale_color_manual(values = c("Burn-in" = "steelblue", "Averaging" = "darkblue")) +
    ggplot2::labs(
      title = "SAC Convergence Trajectory",
      subtitle = sprintf("Target rho* = %.3f | c* = %.4f | Init: %s",
                         x$target_rho, x$c_star, x$init_method),
      x = "Iteration", y = "Scaling Factor c",
      color = "Phase"
    ) +
    ggplot2::annotate("text", x = n_iter * 0.95, y = x$c_star,
             label = sprintf("c* = %.3f", x$c_star),
             vjust = -0.5, color = "red", size = 3.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  # Reliability plot
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = rho, color = phase)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = x$target_rho, color = "red", linetype = "dashed", linewidth = 1) +
    ggplot2::geom_vline(xintercept = burn_in, color = "gray50", linetype = "dotted") +
    ggplot2::scale_color_manual(values = c("Burn-in" = "darkorange", "Averaging" = "darkred")) +
    ggplot2::labs(
      title = "Reliability Estimates Across Iterations",
      subtitle = sprintf("Target rho* = %.3f | Achieved = %.3f",
                         x$target_rho, x$achieved_rho),
      x = "Iteration", y = expression(hat(rho)),
      color = "Phase"
    ) +
    ggplot2::annotate("text", x = n_iter * 0.95, y = x$target_rho,
             label = sprintf("rho* = %.3f", x$target_rho),
             vjust = -0.5, color = "red", size = 3.5) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")

  if (type == "c") {
    return(p1)
  } else if (type == "rho") {
    return(p2)
  } else {
    # Combine with patchwork if available
    if (requireNamespace("patchwork", quietly = TRUE)) {
      return(patchwork::wrap_plots(p1, p2, ncol = 1))
    } else {
      print(p1)
      print(p2)
      return(invisible(list(p1 = p1, p2 = p2)))
    }
  }
}


#' Extract Calibrated Item Parameters from SAC Results
#'
#' @description
#' Returns the calibrated item parameters as a data frame.
#'
#' @param object An object of class \code{"sac_result"}.
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
#' sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
#'                           model = "rasch", n_iter = 200, seed = 42)
#' coef(sac_res)
#' }
#'
#' @export
coef.sac_result <- function(object, ...) {
  .irtsimrel_validate_sac_result_object(object, "object", "SAC coefficient extraction")

  data.frame(
    item_id = seq_len(object$n_items),
    beta = object$beta_vec,
    lambda_base = object$lambda_base,
    lambda_scaled = object$lambda_scaled,
    c_star = rep(object$c_star, object$n_items)
  )
}


#' Predict Reliability at New Scaling Factor Values (SAC)
#'
#' @description
#' If \code{newdata} is NULL, returns the achieved reliability from calibration.
#' If \code{newdata} is a numeric vector of scaling factor values, computes
#' the reliability at each value using the stored item parameters and either
#' the stored theta sample or a user-supplied theta sample.
#' By default, prediction uses the stored SAC theta sample with the
#' pre-calculated latent variance from calibration. If \code{theta_vec} is
#' supplied, reliability is recomputed on that ability sample and
#' \code{theta_var} is set to \code{stats::var(theta_vec)}.
#'
#' @param object An object of class \code{"sac_result"}.
#' @param newdata Optional numeric vector of scaling factor values.
#'   If NULL, returns \code{object$achieved_rho}.
#' @param theta_vec Optional numeric vector of abilities. If not provided,
#'   the stored \code{theta_quad} from the SAC result is used with
#'   \code{object$theta_var}. If provided, \code{theta_vec} defines both the
#'   ability sample and the latent-variance basis: reliability is recomputed
#'   with \code{theta_var = stats::var(theta_vec)}.
#' @param ... Additional arguments (ignored).
#'
#' @return If \code{newdata} is NULL, a single numeric value (achieved reliability).
#'   If \code{newdata} is a numeric vector, a named numeric vector of reliability
#'   values at each scaling factor.
#'
#' @examples
#' \donttest{
#' sac_res <- sac_calibrate(target_rho = 0.75, n_items = 20,
#'                           model = "rasch", n_iter = 200, seed = 42)
#' predict(sac_res)
#' predict(sac_res, newdata = c(0.5, 1.0, 1.5, 2.0))
#' }
#'
#' @export
predict.sac_result <- function(object, newdata = NULL, theta_vec = NULL, ...) {
  .irtsimrel_validate_sac_result_object(object, "object", "SAC prediction")

  if (is.null(newdata)) {
    return(object$achieved_rho)
  }

  c_values <- .irtsimrel_validate_c_values(newdata)

  # Use provided theta_vec or the stored theta_quad
  theta_var <- object$theta_var
  if (is.null(theta_vec)) {
    if (!is.null(object$theta_quad)) {
      theta_vec <- object$theta_quad
    } else {
      stop("No theta sample available. Provide `theta_vec` argument.")
    }
  } else {
    theta_vec <- .irtsimrel_validate_theta_vector(theta_vec)
    theta_var <- stats::var(theta_vec)
  }

  rho_fn <- if (object$metric == "info") compute_rho_tilde else compute_rho_bar

  result <- vapply(c_values, function(c_val) {
    rho_fn(c_val, theta_vec, object$beta_vec,
           object$lambda_base, theta_var = theta_var)
  }, numeric(1))

  names(result) <- paste0("c=", format(c_values, digits = 4))
  result
}
