# =============================================================================
# validation_utils.R
# =============================================================================
# Validation and Comparison Utilities for EQC/SAC Calibration
#
# Contents:
#   - compute_reliability_tam(): Compute WLE/EAP reliability using TAM
#   - simulate_response_data(): Generate item response data from calibration results
#   - compare_eqc_sac(): Compare EQC and SAC calibration results (primary)
#   - compare_eqc_spc(): Deprecated alias for compare_eqc_sac
#
# Author: JoonHo Lee (jlee296@ua.edu)
# Date: December 2025
# =============================================================================


#' Compute WLE and EAP Reliability Using TAM
#'
#' @description
#' Fits a Rasch or 2PL model using \pkg{TAM} and computes WLE and EAP
#' reliability using the official \code{WLErel()} and \code{EAPrel()} functions.
#'
#' @param resp Matrix or data.frame of item responses (0/1).
#' @param model Character. \code{"rasch"} or \code{"2pl"}.
#' @param verbose Logical. If TRUE, print fitting messages.
#' @param ... Additional arguments passed to TAM fitting functions.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{rel_wle}}{WLE reliability.}
#'   \item{\code{rel_eap}}{EAP reliability.}
#'   \item{\code{mod}}{Fitted TAM model object.}
#'   \item{\code{wle}}{Output from \code{TAM::tam.wle()}.}
#' }
#'
#' @details
#' ## WLE vs EAP Reliability
#'
#' TAM defines these reliability coefficients differently:
#'
#' - **WLE reliability**: \eqn{1 - \bar{s}^2 / V_{WLE}}, based on design effect
#' - **EAP reliability**: \eqn{V_{EAP} / (V_{EAP} + \bar{\sigma}^2)}, based on posterior variance
#'
#' In many practical TAM fits, EAP reliability is greater than or equal to WLE
#' reliability, but the two coefficients use different estimators and variance
#' bases. EAP reliability more closely corresponds to MSEM-based population
#' reliability. For conservative inference, inspect WLE alongside EAP.
#'
#' @examples
#' \dontrun{
#' # Simulate response data from calibration results
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   seed = 42
#' )
#' sim_data <- simulate_response_data(result = eqc_result, n_persons = 500)
#'
#' # Compute TAM reliability if TAM is installed
#' if (requireNamespace("TAM", quietly = TRUE)) {
#'   tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
#'   cat(sprintf("WLE reliability: %.4f\n", tam_rel$rel_wle))
#'   cat(sprintf("EAP reliability: %.4f\n", tam_rel$rel_eap))
#' }
#' }
#'
#' @seealso
#' \code{\link{simulate_response_data}} for generating test data,
#' \code{\link{eqc_calibrate}} for calibration.
#'
#' @export
compute_reliability_tam <- function(resp,
                                    model = c("rasch", "2pl"),
                                    verbose = FALSE,
                                    ...) {
  if (!requireNamespace("TAM", quietly = TRUE)) {
    stop("Package 'TAM' is required. Please install it first.")
  }

  model <- match.arg(model)
  resp  <- as.matrix(resp)

  # Fit IRT model
  if (model == "rasch") {
    mod <- TAM::tam.mml(resp = resp, verbose = verbose, ...)
  } else {
    mod <- TAM::tam.mml.2pl(resp = resp, irtmodel = "2PL", verbose = verbose, ...)
  }

  # --- WLE reliability ---
  wle <- TAM::tam.wle(mod)
  if (!all(c("theta", "error") %in% colnames(wle))) {
    stop("`TAM::tam.wle()` output must contain 'theta' and 'error'.")
  }
  rel_wle <- TAM::WLErel(theta = wle$theta, error = wle$error)

  # --- EAP reliability ---
  person <- mod$person
  if (!all(c("EAP", "SD.EAP") %in% colnames(person))) {
    stop("TAM model object does not contain 'EAP' and 'SD.EAP' in $person.")
  }
  rel_eap <- TAM::EAPrel(theta = person$EAP, error = person$SD.EAP)

  list(
    rel_wle = rel_wle,
    rel_eap = rel_eap,
    mod     = mod,
    wle     = wle
  )
}


# Internal validation helpers --------------------------------------------------

.irtsimrel_missing_required_fields <- function(x, fields, arg_name, operation) {
  missing_fields <- setdiff(fields, names(x))
  if (length(missing_fields) > 0L) {
    stop(
      "`", arg_name, "` is missing required field(s) for ", operation, ": ",
      .irtsimrel_backtick_collapse(missing_fields), "."
    )
  }
  invisible(NULL)
}

.irtsimrel_validate_positive_integer_scalar <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1L ||
      is.na(x) || !is.finite(x) || x <= 0 || x %% 1 != 0) {
    stop("`", arg_name, "` must be a positive integer scalar.")
  }
  as.integer(x)
}

.irtsimrel_validate_numeric_scalar_field <- function(x,
                                                     field,
                                                     arg_name,
                                                     positive = FALSE,
                                                     bounds = NULL) {
  value <- x[[field]]
  if (!is.numeric(value) || length(value) != 1L ||
      is.na(value) || !is.finite(value)) {
    stop("`", arg_name, "$", field, "` must be a finite numeric scalar.")
  }
  if (positive && value <= 0) {
    stop("`", arg_name, "$", field, "` must be positive.")
  }
  if (!is.null(bounds) && (value <= bounds[1] || value >= bounds[2])) {
    stop(
      "`", arg_name, "$", field, "` must be in (",
      bounds[1], ", ", bounds[2], ")."
    )
  }
  value
}

.irtsimrel_validate_canonical_metric <- function(metric, arg_name) {
  if (!is.character(metric) || length(metric) != 1L ||
      is.na(metric) || !metric %in% c("info", "msem")) {
    stop(
      "`", arg_name,
      "` must be a canonical stored metric: 'info' or 'msem'."
    )
  }
  metric
}

.irtsimrel_primary_result_class <- function(result) {
  classes <- class(result)
  class_match <- intersect(classes, c("eqc_result", "sac_result", "spc_result"))
  if (length(class_match) == 0L) {
    return(classes[1])
  }
  class_match[1]
}

.irtsimrel_result_call <- function(result) {
  if (is.null(result$call)) {
    return(NULL)
  }
  paste(deparse(result$call), collapse = " ")
}

.irtsimrel_result_item_source <- function(result) {
  if (!is.null(result$items_base) && !is.null(result$items_base$source)) {
    return(result$items_base$source)
  }
  NA_character_
}

.irtsimrel_validate_result_design <- function(result, arg_name, operation) {
  result_class <- .irtsimrel_primary_result_class(result)
  .irtsimrel_missing_required_fields(
    result,
    c(
      "c_star", "target_rho", "achieved_rho", "metric", "model", "n_items",
      "theta_var", "beta_vec", "lambda_base", "lambda_scaled", "items_calib",
      "call"
    ),
    arg_name,
    operation
  )

  .irtsimrel_validate_numeric_scalar_field(result, "c_star", arg_name, positive = TRUE)
  .irtsimrel_validate_numeric_scalar_field(result, "target_rho", arg_name, bounds = c(0, 1))
  achieved_rho <- .irtsimrel_validate_numeric_scalar_field(result, "achieved_rho", arg_name)
  if (achieved_rho < 0 || achieved_rho > 1) {
    stop("`", arg_name, "$achieved_rho` must be in [0, 1].")
  }
  .irtsimrel_validate_numeric_scalar_field(result, "theta_var", arg_name, positive = TRUE)
  .irtsimrel_validate_canonical_metric(result$metric, paste0(arg_name, "$metric"))

  if (!is.character(result$model) || length(result$model) != 1L ||
      is.na(result$model) || !result$model %in% c("rasch", "2pl")) {
    stop("`", arg_name, "$model` must be either 'rasch' or '2pl'.")
  }
  n_items <- .irtsimrel_validate_positive_integer_scalar(
    result$n_items,
    paste0(arg_name, "$n_items")
  )

  if (result_class == "eqc_result") {
    if (is.null(result$misc) || is.null(result$misc$root_status)) {
      stop(
        "`", arg_name,
        "$misc$root_status` is required for response simulation provenance."
      )
    }
  } else if (result_class %in% c("sac_result", "spc_result")) {
    sac_status <- .irtsimrel_result_status(result, "sac")
    if (is.na(sac_status)) {
      stop(
        "`", arg_name,
        "$calibration_status` or `", arg_name,
        "$convergence$status` is required for response simulation provenance."
      )
    }
  }

  beta <- result$beta_vec
  lambda <- result$lambda_scaled

  if (!is.numeric(beta) || length(beta) == 0L || any(!is.finite(beta))) {
    stop("`", arg_name, "$beta_vec` must be a non-empty finite numeric vector.")
  }
  if (!is.numeric(lambda) || length(lambda) == 0L || any(!is.finite(lambda))) {
    stop("`", arg_name, "$lambda_scaled` must be a non-empty finite numeric vector.")
  }
  if (length(beta) != length(lambda)) {
    stop(
      "`", arg_name, "$beta_vec` and `", arg_name,
      "$lambda_scaled` must have the same length."
    )
  }
  if (any(lambda <= 0)) {
    stop("`", arg_name, "$lambda_scaled` must contain positive values.")
  }

  if (n_items != length(beta)) {
    stop(
      "`", arg_name, "$n_items` must match the length of `",
      arg_name, "$beta_vec` and `", arg_name, "$lambda_scaled`."
    )
  }

  lambda_base <- result$lambda_base
  if (!is.numeric(lambda_base) || length(lambda_base) != length(lambda) ||
      any(!is.finite(lambda_base)) || any(lambda_base <= 0)) {
    stop(
      "`", arg_name,
      "$lambda_base` must be a positive finite numeric vector matching ",
      "`", arg_name, "$lambda_scaled`."
    )
  }
  if (!isTRUE(all.equal(lambda_base * result$c_star, lambda, tolerance = 1e-10))) {
    stop(
      "`", arg_name,
      "$lambda_scaled` must match `", arg_name,
      "$lambda_base * ", arg_name, "$c_star` for response simulation."
    )
  }

  if (is.null(result$items_calib$data) ||
      !"lambda" %in% names(result$items_calib$data)) {
    stop(
      "`", arg_name,
      "$items_calib$data$lambda` is required for response simulation provenance."
    )
  }
  item_lambda <- result$items_calib$data$lambda
  if (!is.numeric(item_lambda) || length(item_lambda) != length(lambda) ||
      any(!is.finite(item_lambda))) {
    stop(
      "`", arg_name,
      "$items_calib$data$lambda` must be a finite numeric vector matching ",
      "`", arg_name, "$lambda_scaled`."
    )
  }
  if (!isTRUE(all.equal(item_lambda, lambda, tolerance = 1e-10))) {
    stop(
      "`", arg_name,
      "$items_calib$data$lambda` must match `",
      arg_name,
      "$lambda_scaled` for response simulation."
    )
  }

  list(beta = beta, lambda = lambda, n_items = length(beta))
}

.irtsimrel_validate_comparison_result <- function(result,
                                                  arg_name,
                                                  expected_class,
                                                  operation) {
  if (!inherits(result, expected_class)) {
    stop(
      "`", arg_name, "` must be an object of class '", expected_class, "'."
    )
  }

  .irtsimrel_missing_required_fields(
    result,
    c("c_star", "target_rho", "achieved_rho", "metric", "model", "n_items"),
    arg_name,
    operation
  )

  .irtsimrel_validate_numeric_scalar_field(result, "c_star", arg_name, positive = TRUE)
  .irtsimrel_validate_numeric_scalar_field(result, "target_rho", arg_name, bounds = c(0, 1))
  achieved_rho <- .irtsimrel_validate_numeric_scalar_field(result, "achieved_rho", arg_name)
  if (achieved_rho < 0 || achieved_rho > 1) {
    stop("`", arg_name, "$achieved_rho` must be in [0, 1].")
  }

  .irtsimrel_validate_canonical_metric(result$metric, paste0(arg_name, "$metric"))
  if (!is.character(result$model) || length(result$model) != 1L ||
      is.na(result$model) || !result$model %in% c("rasch", "2pl")) {
    stop("`", arg_name, "$model` must be either 'rasch' or '2pl'.")
  }
  .irtsimrel_validate_positive_integer_scalar(result$n_items, paste0(arg_name, "$n_items"))

  invisible(result)
}

.irtsimrel_result_status <- function(result, kind) {
  if (kind == "eqc") {
    if (!is.null(result$misc) &&
        "root_status" %in% names(result$misc) &&
        !is.null(result$misc[["root_status"]])) {
      return(result$misc[["root_status"]])
    }
    return(NA_character_)
  }

  if ("calibration_status" %in% names(result) &&
      !is.null(result[["calibration_status"]])) {
    return(result[["calibration_status"]])
  }
  if (!is.null(result$convergence) &&
      "status" %in% names(result$convergence) &&
      !is.null(result$convergence[["status"]])) {
    return(result$convergence[["status"]])
  }
  NA_character_
}

.irtsimrel_sac_status_flags <- function(result) {
  if (!is.null(result$convergence) &&
      "status_flags" %in% names(result$convergence) &&
      !is.null(result$convergence[["status_flags"]])) {
    return(result$convergence[["status_flags"]])
  }
  if ("calibration_status" %in% names(result) &&
      !is.null(result[["calibration_status"]])) {
    return(result[["calibration_status"]])
  }
  NA_character_
}

.irtsimrel_result_status_flags <- function(result, result_class, result_status) {
  flags <- result_status

  if (result_class %in% c("sac_result", "spc_result") &&
      !is.null(result$convergence) &&
      !is.null(result$convergence$status_flags)) {
    flags <- result$convergence$status_flags
  }

  flags <- as.character(flags)
  if (length(flags) == 0L || all(is.na(flags))) {
    return(as.character(result_status))
  }
  flags[!is.na(flags)]
}

.irtsimrel_validate_c_values <- function(newdata, arg_name = "newdata") {
  if (!is.numeric(newdata) || length(newdata) == 0L ||
      any(is.na(newdata)) || any(!is.finite(newdata)) ||
      any(newdata <= 0)) {
    stop("`", arg_name, "` must be a positive finite numeric vector.")
  }
  newdata
}

.irtsimrel_validate_theta_vector <- function(theta_vec, arg_name = "theta_vec") {
  if (!is.numeric(theta_vec) || length(theta_vec) < 2L ||
      any(is.na(theta_vec)) || any(!is.finite(theta_vec))) {
    stop("`", arg_name, "` must be a finite numeric vector with at least two values.")
  }
  if (stats::var(theta_vec) <= 1e-10) {
    stop("`", arg_name, "` must have positive variance.")
  }
  theta_vec
}

.irtsimrel_validate_eqc_result_object <- function(object,
                                                  arg_name = "object",
                                                  operation = "EQC S3 method") {
  if (!inherits(object, "eqc_result")) {
    stop("`", arg_name, "` must be an object of class 'eqc_result'.")
  }

  .irtsimrel_validate_result_design(object, arg_name, operation)
  .irtsimrel_missing_required_fields(
    object,
    c("M", "theta_quad", "misc"),
    arg_name,
    operation
  )

  .irtsimrel_validate_positive_integer_scalar(object$M, paste0(arg_name, "$M"))
  .irtsimrel_validate_theta_vector(object$theta_quad, paste0(arg_name, "$theta_quad"))

  if (!is.list(object$misc)) {
    stop("`", arg_name, "$misc` must be a list.")
  }
  misc_missing <- setdiff(
    c("root_status", "c_bounds", "rho_bounds"),
    names(object$misc)
  )
  if (length(misc_missing) > 0L) {
    stop(
      "`", arg_name,
      "$misc` is missing required field(s) for ", operation, ": ",
      .irtsimrel_backtick_collapse(misc_missing), "."
    )
  }
  if (!is.character(object$misc$root_status) ||
      length(object$misc$root_status) != 1L ||
      is.na(object$misc$root_status) ||
      !nzchar(object$misc$root_status)) {
    stop("`", arg_name, "$misc$root_status` must be a non-empty character scalar.")
  }
  if (!is.numeric(object$misc$c_bounds) || length(object$misc$c_bounds) != 2L ||
      any(!is.finite(object$misc$c_bounds))) {
    stop("`", arg_name, "$misc$c_bounds` must be a finite numeric vector of length 2.")
  }
  if (!is.numeric(object$misc$rho_bounds) || length(object$misc$rho_bounds) != 2L ||
      any(!is.finite(object$misc$rho_bounds))) {
    stop("`", arg_name, "$misc$rho_bounds` must be a finite numeric vector of length 2.")
  }

  invisible(object)
}

.irtsimrel_validate_sac_result_object <- function(object,
                                                  arg_name = "object",
                                                  operation = "SAC S3 method") {
  if (!inherits(object, "sac_result")) {
    stop("`", arg_name, "` must be an object of class 'sac_result'.")
  }

  .irtsimrel_validate_result_design(object, arg_name, operation)
  .irtsimrel_missing_required_fields(
    object,
    c(
      "c_final", "trajectory", "rho_trajectory", "theta_quad",
      "raw_trajectory", "step_size_trajectory", "gradient_trajectory",
      "projected", "projection_side", "projection_count", "projection_rate",
      "n_iter", "burn_in", "M_per_iter", "M_pre", "step_params",
      "c_bounds", "c_init", "init_method", "calibration_status",
      "convergence"
    ),
    arg_name,
    operation
  )

  .irtsimrel_validate_numeric_scalar_field(object, "c_final", arg_name, positive = TRUE)
  .irtsimrel_validate_theta_vector(object$theta_quad, paste0(arg_name, "$theta_quad"))
  n_iter <- .irtsimrel_validate_positive_integer_scalar(
    object$n_iter,
    paste0(arg_name, "$n_iter")
  )
  .irtsimrel_validate_positive_integer_scalar(object$M_per_iter, paste0(arg_name, "$M_per_iter"))
  .irtsimrel_validate_positive_integer_scalar(object$M_pre, paste0(arg_name, "$M_pre"))
  if (!is.numeric(object$burn_in) || length(object$burn_in) != 1L ||
      is.na(object$burn_in) || !is.finite(object$burn_in) ||
      object$burn_in < 0 || object$burn_in >= n_iter ||
      object$burn_in %% 1 != 0) {
    stop("`", arg_name, "$burn_in` must be a non-negative integer less than `", arg_name, "$n_iter`.")
  }
  if (!is.character(object$init_method) || length(object$init_method) != 1L ||
      is.na(object$init_method) || !nzchar(object$init_method)) {
    stop("`", arg_name, "$init_method` must be a non-empty character scalar.")
  }
  if (!is.numeric(object$c_init) || length(object$c_init) != 1L ||
      is.na(object$c_init) || !is.finite(object$c_init) || object$c_init <= 0) {
    stop("`", arg_name, "$c_init` must be a positive finite numeric scalar.")
  }
  if (!is.numeric(object$c_bounds) || length(object$c_bounds) != 2L ||
      any(!is.finite(object$c_bounds))) {
    stop("`", arg_name, "$c_bounds` must be a finite numeric vector of length 2.")
  }
  if (!is.numeric(object$trajectory) || length(object$trajectory) != n_iter ||
      any(!is.finite(object$trajectory))) {
    stop("`", arg_name, "$trajectory` must be a finite numeric vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.numeric(object$rho_trajectory) || length(object$rho_trajectory) != n_iter ||
      any(!is.finite(object$rho_trajectory))) {
    stop("`", arg_name, "$rho_trajectory` must be a finite numeric vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.numeric(object$raw_trajectory) || length(object$raw_trajectory) != n_iter ||
      any(!is.finite(object$raw_trajectory))) {
    stop("`", arg_name, "$raw_trajectory` must be a finite numeric vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.numeric(object$step_size_trajectory) ||
      length(object$step_size_trajectory) != n_iter ||
      any(!is.finite(object$step_size_trajectory)) ||
      any(object$step_size_trajectory <= 0)) {
    stop("`", arg_name, "$step_size_trajectory` must be a positive finite numeric vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.numeric(object$gradient_trajectory) ||
      length(object$gradient_trajectory) != n_iter ||
      any(!is.finite(object$gradient_trajectory))) {
    stop("`", arg_name, "$gradient_trajectory` must be a finite numeric vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.logical(object$projected) || length(object$projected) != n_iter ||
      any(is.na(object$projected))) {
    stop("`", arg_name, "$projected` must be a logical vector of length `", arg_name, "$n_iter`.")
  }
  if (!is.character(object$projection_side) ||
      length(object$projection_side) != n_iter ||
      any(is.na(object$projection_side)) ||
      any(!object$projection_side %in% c("none", "lower", "upper"))) {
    stop("`", arg_name, "$projection_side` must contain 'none', 'lower', or 'upper' for each iteration.")
  }
  if (!is.numeric(object$projection_count) ||
      length(object$projection_count) != 1L ||
      is.na(object$projection_count) || !is.finite(object$projection_count) ||
      object$projection_count < 0 || object$projection_count %% 1 != 0) {
    stop("`", arg_name, "$projection_count` must be a non-negative integer scalar.")
  }
  if (!is.numeric(object$projection_rate) || length(object$projection_rate) != 1L ||
      is.na(object$projection_rate) || !is.finite(object$projection_rate) ||
      object$projection_rate < 0 || object$projection_rate > 1) {
    stop("`", arg_name, "$projection_rate` must be a finite numeric scalar in [0, 1].")
  }
  if (!is.list(object$step_params) ||
      !all(c("a", "A", "gamma") %in% names(object$step_params))) {
    stop("`", arg_name, "$step_params` must be a list with fields `a`, `A`, and `gamma`.")
  }
  if (!is.character(object$calibration_status) ||
      length(object$calibration_status) != 1L ||
      is.na(object$calibration_status) ||
      !object$calibration_status %in% c(
        "ok", "not_converged", "hit_lower_bound",
        "hit_upper_bound", "hit_both_bounds"
      )) {
    stop("`", arg_name, "$calibration_status` must be a known SAC status.")
  }
  if (!is.list(object$convergence)) {
    stop("`", arg_name, "$convergence` must be a list.")
  }
  convergence_missing <- setdiff(
    c(
      "converged", "sd_post_burn", "hit_lower_bound", "hit_upper_bound",
      "achieved_gap_abs", "achieved_gap_tolerance", "large_achieved_gap",
      "status_flags"
    ),
    names(object$convergence)
  )
  if (length(convergence_missing) > 0L) {
    stop(
      "`", arg_name,
      "$convergence` is missing required field(s) for ", operation, ": ",
      .irtsimrel_backtick_collapse(convergence_missing), "."
    )
  }

  invisible(object)
}

.irtsimrel_validate_summary_object <- function(x,
                                               class_name,
                                               fields,
                                               arg_name = "x",
                                               operation = "summary print") {
  if (!inherits(x, class_name)) {
    stop("`", arg_name, "` must be an object of class '", class_name, "'.")
  }
  .irtsimrel_missing_required_fields(x, fields, arg_name, operation)
  invisible(x)
}


#' Simulate Item Response Data from Calibration Results
#'
#' @description
#' Generates item response data using the calibrated parameters from
#' \code{eqc_calibrate()} or \code{sac_calibrate()} (formerly
#' \code{spc_calibrate()}).
#'
#' @param result A calibration result object of class \code{"eqc_result"},
#'   \code{"sac_result"}, or \code{"spc_result"} (for backward compatibility),
#'   as returned by \code{eqc_calibrate()} or \code{sac_calibrate()}.
#' @param n_persons Integer. Number of persons to simulate.
#' @param latent_shape Character. Shape argument for \code{sim_latentG()}.
#' @param latent_params List. Additional arguments for \code{sim_latentG()}.
#' @param seed Optional integer for reproducibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{response_matrix}}{N x I matrix of binary responses}
#'   \item{\code{theta}}{True abilities (N x 1)}
#'   \item{\code{beta}}{Item difficulties (I x 1)}
#'   \item{\code{lambda}}{Item discriminations (I x 1)}
#'   \item{\code{provenance}}{List describing the calibration result and
#'     simulation settings used to generate the response data, including
#'     result class, target and achieved reliability, metric, model, test
#'     length, scalar \code{calibration_status}, character-vector
#'     \code{status_flags}, item source/design, calibration call, simulation
#'     seed, sample size, latent shape, and latent parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Using EQC calibration result
#' eqc_result <- eqc_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   seed = 42
#' )
#'
#' sim_data <- simulate_response_data(
#'   result = eqc_result,
#'   n_persons = 1000,
#'   latent_shape = "normal",
#'   seed = 123
#' )
#'
#' # Example 2: Using SAC calibration result
#' sac_result <- sac_calibrate(
#'   target_rho = 0.80,
#'   n_items = 25,
#'   model = "rasch",
#'   reliability_metric = "info",
#'   c_init = eqc_result,
#'   n_iter = 200,
#'   seed = 42
#' )
#'
#' sim_data2 <- simulate_response_data(
#'   result = sac_result,
#'   n_persons = 1000,
#'   latent_shape = "normal",
#'   seed = 123
#' )
#'
#' # Use with TAM for validation when TAM is installed
#' if (requireNamespace("TAM", quietly = TRUE)) {
#'   tam_rel <- compute_reliability_tam(sim_data$response_matrix, model = "rasch")
#' }
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}},
#' \code{\link{compute_reliability_tam}}
#'
#' @export
simulate_response_data <- function(result,
                                   n_persons,
                                   latent_shape = "normal",
                                   latent_params = list(),
                                   seed = NULL) {

  if (!inherits(result, c("eqc_result", "sac_result", "spc_result"))) {
    stop("'result' must be an 'eqc_result', 'sac_result', or legacy 'spc_result' object from eqc_calibrate() or sac_calibrate().")
  }
  n_persons <- .irtsimrel_validate_positive_integer_scalar(n_persons, "n_persons")
  design <- .irtsimrel_validate_result_design(
    result,
    "result",
    "response simulation"
  )

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  # Generate theta
  latent_args <- modifyList(
    list(n = n_persons, shape = latent_shape),
    latent_params
  )
  theta <- do.call(sim_latentG, latent_args)$theta

  beta   <- design$beta
  lambda <- design$lambda
  I <- design$n_items

  # Generate responses (vectorized)
  lin_pred <- outer(theta, beta, "-")
  lin_pred <- sweep(lin_pred, 2, lambda, "*")
  prob_mat <- plogis(lin_pred)

  response_matrix <- matrix(
    rbinom(n_persons * I, size = 1, prob = as.vector(prob_mat)),
    nrow = n_persons, ncol = I
  )

  colnames(response_matrix) <- paste0("item", 1:I)

  result_class <- .irtsimrel_primary_result_class(result)
  result_status <- if (result_class == "eqc_result") {
    .irtsimrel_result_status(result, "eqc")
  } else {
    .irtsimrel_result_status(result, "sac")
  }

  list(
    response_matrix = response_matrix,
    theta = theta,
    beta = beta,
    lambda = lambda,
    provenance = list(
      result_class = result_class,
      c_star = result$c_star,
      target_rho = result$target_rho,
      achieved_rho = result$achieved_rho,
      metric = result$metric,
      model = result$model,
      n_items = result$n_items,
      calibration_status = result_status,
      status_flags = .irtsimrel_result_status_flags(result, result_class, result_status),
      item_design = if (!is.null(result$item_design)) result$item_design else NA_character_,
      item_source = .irtsimrel_result_item_source(result),
      calibration_call = .irtsimrel_result_call(result),
      simulation_seed = seed,
      n_persons = n_persons,
      latent_shape = latent_shape,
      latent_params = latent_params
    )
  )
}


# =============================================================================
# Comparison Function for EQC and SAC Results
# =============================================================================

#' Compare EQC and SAC Calibration Results
#'
#' @description
#' Compares the calibration results from EQC and SAC algorithms.
#' The target reliability must match exactly. Differences in model, test
#' length, or stored reliability metric are reported as warnings because they
#' make \code{c*} agreement a weaker diagnostic.
#'
#' \code{compare_eqc_spc()} is a deprecated backward-compatible alias for
#' \code{compare_eqc_sac()}.
#'
#' @param eqc_result An object of class \code{"eqc_result"}.
#' @param sac_result An object of class \code{"sac_result"} (or \code{"spc_result"}
#'   for backward compatibility).
#' @param verbose Logical. If TRUE, print comparison summary.
#'
#' @return A list with comparison statistics (invisibly):
#' \describe{
#'   \item{\code{c_eqc}}{Calibrated c* from EQC.}
#'   \item{\code{c_sac}}{Calibrated c* from SAC.}
#'   \item{\code{diff_abs}}{Absolute difference between c* values.}
#'   \item{\code{diff_pct}}{Percent difference relative to EQC.}
#'   \item{\code{agreement}}{Logical. TRUE if difference is < 5%.}
#'   \item{\code{target_rho}}{Target reliability.}
#'   \item{\code{achieved_eqc}}{Achieved reliability from EQC.}
#'   \item{\code{achieved_sac}}{Achieved reliability from SAC.}
#'   \item{\code{achieved_diff_abs}}{Absolute achieved-reliability difference.}
#'   \item{\code{metric_eqc}, \code{metric_sac}}{Stored canonical reliability
#'     metrics.}
#'   \item{\code{model_eqc}, \code{model_sac}}{Measurement models.}
#'   \item{\code{n_items_eqc}, \code{n_items_sac}}{Test lengths.}
#'   \item{\code{eqc_status}}{EQC root/calibration status when available.}
#'   \item{\code{sac_status}}{SAC canonical calibration status when available.}
#'   \item{\code{sac_status_flags}}{SAC multi-condition status flags when
#'     available.}
#' }
#'
#' @examples
#' \donttest{
#' # Run both algorithms
#' eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
#' sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
#'                             reliability_metric = "info",
#'                             c_init = eqc_result, seed = 42)
#'
#' # Compare results
#' compare_eqc_sac(eqc_result, sac_result)
#' }
#'
#' @seealso
#' \code{\link{eqc_calibrate}}, \code{\link{sac_calibrate}}
#'
#' @export
compare_eqc_sac <- function(eqc_result, sac_result, verbose = TRUE) {

  if (!inherits(sac_result, c("sac_result", "spc_result"))) {
    stop("sac_result must be an object of class 'sac_result' (or 'spc_result' for backward compatibility).")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.")
  }

  .irtsimrel_validate_comparison_result(
    eqc_result,
    "eqc_result",
    "eqc_result",
    "EQC/SAC comparison"
  )
  .irtsimrel_validate_comparison_result(
    sac_result,
    "sac_result",
    class(sac_result)[1],
    "EQC/SAC comparison"
  )

  # Configuration validation
  if (!isTRUE(all.equal(eqc_result$target_rho, sac_result$target_rho))) {
    stop(sprintf(
      "target_rho differs between EQC (%.4f) and SAC (%.4f). Results are not comparable.",
      eqc_result$target_rho, sac_result$target_rho
    ))
  }
  if (!identical(eqc_result$model, sac_result$model)) {
    warning(sprintf(
      "Model differs between EQC ('%s') and SAC ('%s').",
      eqc_result$model, sac_result$model
    ))
  }
  if (!identical(eqc_result$n_items, sac_result$n_items)) {
    warning(sprintf(
      "n_items differs between EQC (%d) and SAC (%d).",
      eqc_result$n_items, sac_result$n_items
    ))
  }
  if (!identical(eqc_result$metric, sac_result$metric)) {
    warning(sprintf(
      "Reliability metric differs between EQC ('%s') and SAC ('%s').",
      eqc_result$metric, sac_result$metric
    ))
  }

  eqc_status <- .irtsimrel_result_status(eqc_result, "eqc")
  sac_status <- .irtsimrel_result_status(sac_result, "sac")
  sac_status_flags <- .irtsimrel_sac_status_flags(sac_result)
  if (!is.na(eqc_status) &&
      !eqc_status %in% c("ok", "uniroot_success", "boundary_lower", "boundary_upper")) {
    warning(sprintf(
      "EQC root status is '%s'; comparison may reflect a boundary or infeasible target.",
      eqc_status
    ))
  }
  if (!is.na(sac_status) && sac_status != "ok") {
    warning(sprintf(
      "SAC calibration status is '%s'; comparison may reflect non-convergence or boundary projection.",
      sac_status
    ))
  }

  # Extract key values
  c_eqc <- eqc_result$c_star
  c_sac <- sac_result$c_star

  diff_abs <- abs(c_eqc - c_sac)
  diff_pct <- 100 * diff_abs / c_eqc
  achieved_diff_abs <- abs(eqc_result$achieved_rho - sac_result$achieved_rho)

  agreement <- diff_pct < 5  # Within 5%

  comparison <- list(
    c_eqc = c_eqc,
    c_sac = c_sac,
    diff_abs = diff_abs,
    diff_pct = diff_pct,
    agreement = agreement,
    target_rho = sac_result$target_rho,
    achieved_eqc = eqc_result$achieved_rho,
    achieved_sac = sac_result$achieved_rho,
    achieved_diff_abs = achieved_diff_abs,
    metric_eqc = eqc_result$metric,
    metric_sac = sac_result$metric,
    model_eqc = eqc_result$model,
    model_sac = sac_result$model,
    n_items_eqc = eqc_result$n_items,
    n_items_sac = sac_result$n_items,
    eqc_status = eqc_status,
    sac_status = sac_status,
    sac_status_flags = sac_status_flags
  )

  if (verbose) {
    message("")
    message("=======================================================")
    message("  EQC vs SAC Comparison")
    message("=======================================================")
    message("")
    message(sprintf("  Target reliability  : %.4f", comparison$target_rho))
    message(sprintf("  EQC c*              : %.6f", c_eqc))
    message(sprintf("  SAC c*              : %.6f", c_sac))
    message(sprintf("  Absolute difference : %.6f", diff_abs))
    message(sprintf("  Percent difference  : %.2f%%", diff_pct))
    message(sprintf("  EQC achieved rho    : %.4f", comparison$achieved_eqc))
    message(sprintf("  SAC achieved rho    : %.4f", comparison$achieved_sac))
    message(sprintf("  Agreement (< 5%%)    : %s", ifelse(agreement, "YES", "NO")))
    message(sprintf("  EQC status          : %s", comparison$eqc_status))
    message(sprintf("  SAC status          : %s", comparison$sac_status))
    message(sprintf("  SAC status flags    : %s",
                    paste(comparison$sac_status_flags, collapse = ", ")))
    message("")
  }

  invisible(comparison)
}


#' Deprecated Alias for EQC/SAC Comparison
#'
#' @description
#' `compare_eqc_spc()` is a deprecated backward-compatible alias for
#' `compare_eqc_sac()`. It is retained for scripts written against the v0.1.x
#' naming convention and will issue a deprecation warning when called.
#'
#' @param ... Arguments passed to \code{compare_eqc_sac()}.
#'
#' @return The result of \code{\link{compare_eqc_sac}}.
#'
#' @seealso \code{\link{compare_eqc_sac}}
#'
#' @export
compare_eqc_spc <- function(...) {
  .Deprecated("compare_eqc_sac")
  compare_eqc_sac(...)
}
