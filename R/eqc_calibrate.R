# =============================================================================
# eqc_calibrate.R
# =============================================================================
# Empirical Quadrature Calibration (Algorithm 1: EQC)
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
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# =============================================================================


#' Empirical Quadrature Calibration (Algorithm 1: EQC)
#'
#' @description
#' `eqc_calibrate()` implements Algorithm 1 (Empirical Quadrature
#' Calibration, EQC) for reliability-targeted IRT simulation.
#'
#' Given a target marginal reliability \eqn{\rho^*}, a latent distribution
#' generator `sim_latentG()` (for \eqn{G}) and an item parameter generator
#' `sim_item_params()` (for \eqn{H}), the function searches for a global
#' discrimination scale \eqn{c^* > 0} such that a selected empirical
#' quadrature approximation \eqn{\hat\rho_M(c)} for the Rasch, 2PL, or 3PL
#' model matches \eqn{\rho^*}. Population accuracy is a Monte Carlo
#' approximation and should be assessed on independent draws.
#'
#' The key idea is to:
#' \enumerate{
#'   \item Draw a large fixed "quadrature" sample
#'     \eqn{\{\theta_m\}_{m=1}^M \sim G} and item parameters
#'     \eqn{\{(\beta_i, \lambda_{i,0})\}_{i=1}^I \sim H} once.
#'   \item For any scale \eqn{c}, form \eqn{\lambda_i(c) = c \cdot \lambda_{i,0}}
#'     and compute the empirical approximation to population reliability
#'     \eqn{\hat\rho_M(c)} from the test information function.
#'   \item Scan the explicit scale interval on a logarithmic grid, enumerate
#'     all detected roots of \eqn{\hat\rho_M(c^*) = \rho^*}, and apply the
#'     requested branch-selection policy.
#' }
#'
#' @param target_rho Numeric in (0, 1). Target marginal reliability \eqn{\rho^*}.
#'
#' @param n_items Integer. Number of items in the test form.
#'
#' @param model Character. Measurement model: \code{"rasch"}, \code{"2pl"},
#'   or \code{"3pl"}.
#'   For \code{"rasch"}, all baseline discriminations are set to 1 before scaling.
#'   For \code{"3pl"}, \code{item_params$guessing_params} controls the item
#'   lower asymptotes; calibration scales discrimination only and uses the
#'   fixed logistic convention `D = 1`.
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
#' @param reliability_metric Character. Reliability definition used inside EQC:
#'   \describe{
#'     \item{\code{"info"}}{Average-information reliability (default, recommended for EQC).
#'       Targets \eqn{\tilde{\rho}}, which is well behaved on practical EQC
#'       search intervals and is bracket-checked before root-finding.}
#'     \item{\code{"tilde"}}{Synonym for \code{"info"}.}
#'     \item{\code{"msem"}, \code{"bar"}}{Rejected by \code{eqc_calibrate()} because
#'       direct MSEM targeting can be non-monotone under EQC root-finding. Use
#'       \code{\link{sac_calibrate}} for direct MSEM-based calibration.}
#'   }
#'
#' @param M Integer. Size of the empirical quadrature sample (default: 10000).
#'
#' @param c_bounds Numeric length-2 vector. Search bounds for \eqn{c}. Default: c(0.3, 3).
#'
#' @param tol Numeric. Root-polishing tolerance. Default: 1e-4.
#'
#' @param seed Optional integer for reproducibility.
#'
#' @param verbose Logical. If TRUE, print progress messages.
#'
#' @param root_policy Character. Root-selection policy used when the fixed-form
#'   reliability curve has multiple roots. The default,
#'   \code{"lowest_increasing"}, selects the lowest-scale root with positive
#'   local slope. \code{"nearest_increasing"} uses the unscaled baseline
#'   \code{c = 1} as its reference. \code{"lowest_any"} and
#'   \code{"highest_any"} select the lowest- or highest-scale admissible root
#'   regardless of slope direction. Tangent and plateau roots remain
#'   inadmissible under the \code{"*_any"} policies unless
#'   \code{allow_tangent} or \code{allow_plateau}, respectively, is enabled.
#'
#' @param root_controls Named list of advanced topology-scanner controls. The
#'   supplied \code{c_bounds} are never expanded.
#'
#' @param allow_tangent Logical. If \code{TRUE}, a detected tangent root may be
#'   selected by an \code{"*_any"} policy. The default is \code{FALSE}.
#'
#' @param allow_plateau Logical. If \code{TRUE}, the geometric midpoint of a
#'   detected target plateau may be selected by an \code{"*_any"} policy. The
#'   complete plateau interval remains recorded in the root table. The default
#'   is \code{FALSE}.
#'
#' @return An object of class \code{"eqc_result"} (a list) with elements:
#' \describe{
#'   \item{\code{c_star}}{Calibrated discrimination scale \eqn{c^*}.}
#'   \item{\code{target_rho}}{Target reliability \eqn{\rho^*}.}
#'   \item{\code{achieved_rho}}{Empirical quadrature estimate \eqn{\hat\rho_M(c^*)}.}
#'   \item{\code{metric}}{Reliability metric used.}
#'   \item{\code{model}}{Measurement model.}
#'   \item{\code{n_items}}{Number of items.}
#'   \item{\code{M}}{Empirical quadrature sample size.}
#'   \item{\code{theta_quad}}{Length-M vector of quadrature abilities.}
#'   \item{\code{theta_var}}{Sample variance of theta_quad.}
#'   \item{\code{beta_vec}}{Item difficulties from the baseline item design.}
#'   \item{\code{lambda_base}}{Baseline (unscaled) item discriminations.}
#'   \item{\code{lambda_scaled}}{Scaled item discriminations
#'     (\code{lambda_base * c_star}).}
#'   \item{\code{guessing_vec}}{Item lower asymptotes. New Rasch/2PL objects
#'     store an all-zero vector; 3PL objects store the generated guessing
#'     parameters.}
#'   \item{\code{items_base}}{item_params object with scale = 1 (baseline).}
#'   \item{\code{items_calib}}{item_params object with discriminations scaled by c_star.}
#'   \item{\code{schema_version}, \code{item_scope}}{Additive result-schema and
#'     fixed-form estimand metadata.}
#'   \item{\code{estimand_signature}, \code{design_signature}}{Canonical
#'     estimand and realized-design provenance used by comparison utilities.}
#'   \item{\code{call}}{Matched function call.}
#'   \item{\code{misc}}{List of root-finding diagnostics, including bounds,
#'     endpoint reliabilities, global achievable range, root topology, branch
#'     selection, and a legacy-compatible \code{uniroot_result} for an interior
#'     crossing.}
#' }
#'
#' @details
#' ## Reliability Metrics
#'
#' The function calibrates the average-information reliability definition:
#'
#' - **Average-information** (\code{"info"}/\code{"tilde"}, **default**): Uses the arithmetic mean,
#'   \eqn{\tilde{\rho}(c) = \sigma^2_\theta \bar{\mathcal{J}}(c) / (\sigma^2_\theta \bar{\mathcal{J}}(c) + 1)}.
#'   By Jensen's inequality, \eqn{\tilde{\rho} \geq \bar{w}}, so this metric typically
#'   yields higher reliability values. **This is the recommended default for EQC** because
#'   the objective function \eqn{\tilde{\rho}(c) - \rho^*} is stable on practical
#'   search intervals. The implementation nevertheless scans the complete
#'   user-supplied interval on \code{log(c)}, polishes detected local extrema
#'   and successfully bracketed sign-changing crossings, records any
#'   near-target fallback separately, and applies an explicit root policy. It
#'   does not assume that the empirical curve is globally monotone or unimodal.
#'
#' **MSEM-based** (\code{"msem"}/\code{"bar"}) reliability uses the harmonic mean of test information,
#'   \eqn{\bar{w}(c) = \sigma^2_\theta / (\sigma^2_\theta + E[1/\mathcal{J}(\theta;c)])}.
#' This metric is rejected by \code{eqc_calibrate()} because the objective
#' \eqn{\bar{w}(c) - \rho^*} may be non-monotone (see the reliability-targeted
#' simulation manuscript, Section 4.3). Although the shared scanner can
#' diagnose such topology, MSEM remains outside EQC's calibration contract.
#' Use \code{\link{sac_calibrate}} if you need to target \eqn{\bar{w}} directly.
#'
#' ## WLE vs EAP Reliability Interpretation
#'
#' WLE and EAP reliability use different estimators and variance bases, so they
#' are complementary external diagnostics rather than ordered versions of one
#' estimand. EAP is more closely related to MSEM-based population reliability;
#' neither coefficient is identical to EQC's information-based target. The
#' public [compute_reliability_tam()] helper supports Rasch/2PL only. Phase 7's
#' direct TAM 3PL check was EAP-only because a validated TAM 3PL WLE path was
#' unavailable.
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
#'
#' # 3PL: guessing is generated once with the fixed item form and is not scaled.
#' eqc_3pl <- eqc_calibrate(
#'   target_rho = 0.70,
#'   n_items = 12L,
#'   model = "3pl",
#'   item_params = list(
#'     guessing_params = list(distribution = "fixed", value = 0.20)
#'   ),
#'   M = 1500L,
#'   c_bounds = c(0.1, 5),
#'   seed = 5102L
#' )
#' stopifnot(all(eqc_3pl$guessing_vec == 0.20))
#' }
#'
#' \dontrun{
#' # EQC with IRW difficulties (requires irw package)
#' if (requireNamespace("irw", quietly = TRUE)) {
#'   eqc_result2 <- eqc_calibrate(
#'     target_rho = 0.80,
#'     n_items = 25,
#'     model = "rasch",
#'     item_source = "irw",
#'     seed = 42,
#'     verbose = TRUE
#'   )
#' }
#' }
#'
#' @export
eqc_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl", "3pl"),
                          latent_shape = "normal",
                          item_source = "parametric",
                          latent_params = list(),
                          item_params = list(),
                          reliability_metric = "info",
                          M = 10000L,
                          c_bounds = c(0.3, 3),
                          tol = 1e-4,
                          seed = NULL,
                          verbose = FALSE,
                          root_policy = "lowest_increasing",
                          root_controls = list(),
                          allow_tangent = FALSE,
                          allow_plateau = FALSE) {

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
  latent_shape <- .irtsimrel_match_latent_shape(latent_shape)

  valid_metrics <- c("info", "tilde", "msem", "bar")
  if (!is.character(reliability_metric) || length(reliability_metric) != 1L ||
      is.na(reliability_metric) || !(reliability_metric %in% valid_metrics)) {
    stop(
      "`reliability_metric` must be one of \"info\", \"tilde\", ",
      "\"msem\", or \"bar\".",
      call. = FALSE
    )
  }
  if (reliability_metric %in% c("msem", "bar")) {
    stop(
      "`eqc_calibrate()` supports `reliability_metric = \"info\"` ",
      "(or \"tilde\") only. MSEM-based targets (\"msem\"/\"bar\") can be ",
      "non-monotone under EQC root-finding; use ",
      "`sac_calibrate(reliability_metric = \"msem\")` for direct w-bar ",
      "calibration, or `compute_rho_bar()`/`rho_curve(metric = \"msem\")` ",
      "for diagnostics.",
      call. = FALSE
    )
  }
  metric_internal <- "info"

  M <- .irtsimrel_validate_positive_integer_scalar(M, "M")
  if (M < 2L) {
    stop("`M` must be at least 2 to estimate latent variance.")
  }

  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      any(!is.finite(c_bounds)) || any(c_bounds <= 0) ||
      c_bounds[1] >= c_bounds[2]) {
    stop("`c_bounds` must be a numeric vector (c_min, c_max) with 0 < c_min < c_max.")
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite numeric scalar.")
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  valid_root_policies <- c(
    "lowest_increasing", "nearest_increasing", "lowest_any", "highest_any"
  )
  if (!is.character(root_policy) || length(root_policy) != 1L ||
      is.na(root_policy) || !(root_policy %in% valid_root_policies)) {
    stop(
      "`root_policy` must be one of ",
      paste(sprintf("'%s'", valid_root_policies), collapse = ", "),
      "."
    )
  }
  root_controls <- .irtsimrel_check_list_arg(root_controls, "root_controls")
  if (!is.logical(allow_tangent) || length(allow_tangent) != 1L ||
      is.na(allow_tangent)) {
    stop("`allow_tangent` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(allow_plateau) || length(allow_plateau) != 1L ||
      is.na(allow_plateau)) {
    stop("`allow_plateau` must be TRUE or FALSE.", call. = FALSE)
  }

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)
  item_params <- .irtsimrel_normalize_item_params(item_params)

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  # ===========================================================================
  # Step 1: Generate Quadrature Samples from G and H
  # ===========================================================================

  if (verbose) message("Step 1: Generating quadrature samples...")

  # --- Latent abilities ---
  latent_args <- modifyList(
    list(n = M, shape = latent_shape),
    latent_params
  )
  quad_latent <- do.call(sim_latentG, latent_args)
  theta_quad  <- quad_latent$theta
  theta_var   <- var(theta_quad)
  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 1e-10) {
    stop("Quadrature latent variance must be finite and positive.")
  }

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

  guessing_vec <- if (identical(model, "3pl")) {
    if (!"guessing" %in% names(items_base$data)) {
      stop("The generated 3PL item design is missing `guessing` parameters.")
    }
    as.numeric(items_base$data$guessing)
  } else {
    rep(0, n_items)
  }

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

  # ===========================================================================
  # Step 2: Define Reliability Function rho(c)
  # ===========================================================================

  compute_rho <- function(c) {
    .compute_rho_generic(
      c = c,
      theta_vec = theta_quad,
      beta_vec = beta_vec,
      lambda_base = lambda_base,
      theta_var = theta_var,
      metric_internal = metric_internal,
      guessing = guessing_vec
    )
  }

  # ===========================================================================
  # Step 3: Root-Finding with Bracket Checks
  # ===========================================================================

  if (verbose) message("Step 2: Running root-finding algorithm...")

  topology <- .irtsimrel_scan_topology(
    fn = compute_rho,
    c_bounds = c_bounds,
    target = target_rho,
    tol = tol,
    controls = root_controls
  )

  if (!identical(topology$topology_status, "resolved")) {
    .irtsimrel_abort_topology(
      paste0(
        "EQC root topology could not be resolved within `c_bounds` (status: ",
        topology$topology_status, "; stop reason: ", topology$stop_reason,
        "). Increase the scanner budget through `root_controls` or revise ",
        "the search interval."
      ),
      subclass = "irtsimrel_topology_unresolved",
      topology = topology
    )
  }

  rho_bounds <- topology$rho_bounds
  rho_at_low <- unname(rho_bounds[[1L]])
  rho_at_high <- unname(rho_bounds[[2L]])
  val_low <- rho_at_low - target_rho
  val_high <- rho_at_high - target_rho
  rho_range <- topology$rho_range
  target_status <- topology$target_status
  bracket_method <- "adaptive_log_grid"

  extrema_candidates <- topology$scan[
    topology$scan$finite,
    c("c", "rho"),
    drop = FALSE
  ]
  if (is.data.frame(topology$extrema) &&
      all(c("c", "rho") %in% names(topology$extrema))) {
    extrema_candidates <- rbind(
      extrema_candidates,
      topology$extrema[, c("c", "rho"), drop = FALSE]
    )
  }
  max_index <- which.max(extrema_candidates$rho)
  rho_max <- unname(rho_range[[2L]])
  c_at_rho_max <- extrema_candidates$c[[max_index]]

  if (verbose) {
    message(sprintf("  At c = %.3f: rho = %.4f, g = %.4f", c_bounds[1], rho_at_low, val_low))
    message(sprintf("  At c = %.3f: rho = %.4f, g = %.4f", c_bounds[2], rho_at_high, val_high))
  }

  root_status <- "ok"
  calibration_status <- "ok"
  selected_root_id <- NA_integer_
  selected_root_kind <- NA_character_
  selected_root_direction <- NA_character_
  uniroot_result <- NULL

  if (target_status == "infeasible_below_range") {
    c_star <- topology$best_achievable$c
    root_status <- "below_lower"
    calibration_status <- "infeasible"
    warning(sprintf(
      paste0(
        "Target rho* = %.3f is below the minimum achievable reliability.\n",
        "  - Best detected c in bounds: %.4f\n",
        "  - Best detected rho: %.4f\n",
        "Suggestion: Decrease the lower search bound or increase target_rho.\n",
        "Returning the best achievable c_star = %.4f."
      ),
      target_rho, c_star, topology$best_achievable$rho, c_star
    ))

  } else if (target_status == "infeasible_above_range") {
    c_star <- topology$best_achievable$c
    root_status <- "above_upper"
    calibration_status <- "infeasible"
    warning(sprintf(
      paste0(
        "Target rho* = %.3f exceeds the maximum achievable reliability for this configuration.\n",
        "  - Items: %d\n",
        "  - Best detected c in bounds: %.4f\n",
        "  - Best detected rho: %.4f\n",
        "  - At c = %.2f (upper bound): rho = %.4f\n",
        "  - Gap: %.4f\n\n",
        "Suggestions:\n",
        "  1. Revise the explicit c_bounds if scientifically justified\n",
        "  2. Increase n_items for higher achievable reliability\n",
        "  3. Accept achieved rho = %.4f as the maximum detected for this design\n\n",
        "Returning c_star = %.4f."
      ),
      target_rho, n_items, c_at_rho_max, rho_max, c_bounds[2], rho_at_high,
      target_rho - rho_max, rho_max, c_star
    ))

  } else if (target_status == "uncertain") {
    .irtsimrel_abort_topology(
      paste0(
        "EQC target feasibility is uncertain because the reliability ",
        "topology contains unresolved or non-finite regions."
      ),
      subclass = "irtsimrel_topology_unresolved",
      topology = topology
    )

  } else {
    selection <- .irtsimrel_select_root(
      topology = topology,
      root_policy = root_policy,
      reference = if (identical(root_policy, "nearest_increasing")) 1 else NULL,
      allow_tangent = allow_tangent,
      allow_plateau = allow_plateau
    )
    topology <- selection$topology
    c_star <- selection$c
    uniroot_result <- selection$uniroot_result
    selected_root_id <- selection$root_id
    selected_root_kind <- selection$root$kind[[1L]]
    selected_root_direction <- selection$root$direction[[1L]]

    at_lower <- abs(log(c_star / c_bounds[1])) <= max(tol, 1e-8)
    at_upper <- abs(log(c_star / c_bounds[2])) <= max(tol, 1e-8)
    if (identical(selected_root_kind, "boundary") && at_lower) {
      root_status <- "boundary_lower"
      calibration_status <- "boundary"
      uniroot_result <- NULL
    } else if (identical(selected_root_kind, "boundary") && at_upper) {
      root_status <- "boundary_upper"
      calibration_status <- "boundary"
      uniroot_result <- NULL
    } else if (identical(selected_root_kind, "tangent")) {
      root_status <- "tangent_root"
      calibration_status <- "tangent"
      uniroot_result <- NULL
    } else if (identical(selected_root_kind, "plateau")) {
      root_status <- "plateau_root"
      calibration_status <- "plateau"
      uniroot_result <- NULL
    } else {
      root_status <- "uniroot_success"
      calibration_status <- if (topology$root_count > 1L) {
        "ok_multiple_roots"
      } else {
        "ok"
      }
    }
  }

  # exp(log(boundary)) can differ from the supplied endpoint by one ulp.
  # Clamp only that floating-point representation error; never expand bounds.
  c_star <- min(max(c_star, c_bounds[[1L]]), c_bounds[[2L]])

  achieved_rho <- compute_rho(c_star)

  if (verbose && root_status == "uniroot_success" && target_rho > rho_max * 0.95) {
    message(sprintf(
      "Note: Target rho* = %.3f is near the achievable maximum (%.3f) for this configuration.",
      target_rho, rho_max
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

  items_calib <- .irtsimrel_apply_item_scale(items_base, lambda_base, c_star)

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
      rho_bounds_named = rho_bounds,
      rho_range = rho_range,
      rho_max = rho_max,
      c_at_rho_max = c_at_rho_max,
      bracket_method = bracket_method,
      target_status = target_status,
      topology_status = topology$topology_status,
      calibration_status = calibration_status,
      root_count = topology$root_count,
      admissible_root_count = topology$admissible_root_count,
      root_policy = root_policy,
      root_reference = if (identical(root_policy, "nearest_increasing")) 1 else NULL,
      allow_tangent = allow_tangent,
      allow_plateau = allow_plateau,
      selected_root_id = selected_root_id,
      selected_root_kind = selected_root_kind,
      selected_root_direction = selected_root_direction,
      topology_evaluations = topology$eval_count,
      topology_stop_reason = topology$stop_reason,
      roots = topology$roots,
      extrema = topology$extrema,
      topology = topology,
      uniroot_result = uniroot_result
    ),
    guessing_vec = guessing_vec,
    schema_version = 1L,
    item_scope = "fixed_form",
    estimand_signature = list(
      metric = metric_internal,
      theta_measure = "population",
      approximation_method = "fixed_mc",
      item_measure = "fixed_form"
    ),
    design_signature = list(
      schema_version = 1L,
      model = model,
      metric = metric_internal,
      n_items = n_items,
      item_scope = "fixed_form",
      latent_shape = latent_shape,
      latent_params = latent_params,
      beta = beta_vec,
      lambda_base = lambda_base,
      guessing = guessing_vec
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
  .irtsimrel_validate_eqc_result_object(x, "x", "EQC print")

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
  cat("  Reliability metric           : Average-information (tilde)\n")
  cat(sprintf("  Latent variance              : %.4f\n", x$theta_var))
  cat("\n")

  cat("Convergence:\n")
  cat(sprintf("  Root status                  : %s\n", x$misc$root_status))
  if (!is.null(x$misc$calibration_status)) {
    cat(sprintf("  Calibration status           : %s\n", x$misc$calibration_status))
  }
  if (!is.null(x$misc$root_count)) {
    cat(sprintf("  Roots detected               : %d\n", x$misc$root_count))
  }
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
  if (identical(x$model, "3pl")) {
    design <- .irtsimrel_normalize_result_item_design(x, "x", "EQC print")
    cat(sprintf(
      "  guessing:      mean = %.3f, sd = %.3f, range = [%.2f, %.2f]\n",
      mean(design$guessing), sd(design$guessing),
      min(design$guessing), max(design$guessing)
    ))
  }
  cat("\n")

  invisible(x)
}


#' @rdname eqc_calibrate
#' @param object An object of class \code{"eqc_result"}.
#' @return An object of class \code{"summary.eqc_result"} containing key
#'   calibration results.
#' @export
summary.eqc_result <- function(object, ...) {
  .irtsimrel_validate_eqc_result_object(object, "object", "EQC summary")

  out <- list(
    c_star       = object$c_star,
    target_rho   = object$target_rho,
    achieved_rho = object$achieved_rho,
    metric       = object$metric,
    model        = object$model,
    n_items      = object$n_items,
    M            = object$M,
    root_status  = object$misc$root_status,
    theta_var    = object$theta_var,
    calibration_status = object$misc$calibration_status,
    root_count = object$misc$root_count,
    root_policy = object$misc$root_policy
  )
  if (identical(object$model, "3pl")) {
    design <- .irtsimrel_normalize_result_item_design(
      object, "object", "EQC summary"
    )
    out$guessing_mean <- mean(design$guessing)
    out$guessing_range <- range(design$guessing)
  }
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
  .irtsimrel_validate_summary_object(
    x,
    "summary.eqc_result",
    c(
      "c_star", "target_rho", "achieved_rho", "metric", "model",
      "n_items", "M", "root_status", "theta_var"
    ),
    "x",
    "EQC summary print"
  )

  cat("Summary: Empirical Quadrature Calibration (EQC)\n")
  cat("================================================\n")
  cat(sprintf("  Model            : %s\n", toupper(x$model)))
  cat("  Metric           : Average-information (tilde)\n")
  cat(sprintf("  Number of items  : %d\n", x$n_items))
  cat(sprintf("  Quadrature (M)   : %d\n", x$M))
  cat(sprintf("  Latent variance  : %.*f\n", digits, x$theta_var))
  cat("\nCalibration Results:\n")
  cat(sprintf("  Target rho*      : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved rho     : %.*f\n", digits, x$achieved_rho))
  cat(sprintf("  Absolute error   : %.2e\n", abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor c*: %.*f\n", digits, x$c_star))
  cat(sprintf("  Root status      : %s\n", x$root_status))
  if (!is.null(x$calibration_status)) {
    cat(sprintf("  Calibration      : %s\n", x$calibration_status))
  }
  if (!is.null(x$root_count)) {
    cat(sprintf("  Roots detected   : %d\n", x$root_count))
  }
  if (identical(x$model, "3pl") && !is.null(x$guessing_mean)) {
    cat(sprintf(
      "  Guessing         : mean %.*f, range [%.*f, %.*f]\n",
      digits, x$guessing_mean,
      digits, x$guessing_range[1], digits, x$guessing_range[2]
    ))
  }
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
#'   \item{\code{guessing}}{For 3PL results only, the item lower asymptote.
#'     This column is omitted for Rasch and 2PL results.}
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
  .irtsimrel_validate_eqc_result_object(object, "object", "EQC coefficient extraction")

  out <- data.frame(
    item_id = seq_len(object$n_items),
    beta = object$beta_vec,
    lambda_base = object$lambda_base,
    lambda_scaled = object$lambda_scaled,
    c_star = rep(object$c_star, object$n_items)
  )
  if (identical(object$model, "3pl")) {
    design <- .irtsimrel_normalize_result_item_design(
      object, "object", "EQC coefficient extraction"
    )
    out$guessing <- design$guessing
  }
  out
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
  .irtsimrel_validate_eqc_result_object(object, "object", "EQC prediction")

  if (is.null(newdata)) {
    return(object$achieved_rho)
  }

  c_values <- .irtsimrel_validate_c_values(newdata)
  design <- .irtsimrel_normalize_result_item_design(
    object, "object", "EQC prediction"
  )

  rho_fn <- if (object$metric == "info") compute_rho_tilde else compute_rho_bar

  result <- vapply(c_values, function(c_val) {
    rho_fn(c_val, object$theta_quad, object$beta_vec,
           object$lambda_base, theta_var = object$theta_var,
           guessing = design$guessing)
  }, numeric(1))

  names(result) <- paste0("c=", format(c_values, digits = 4))
  result
}
