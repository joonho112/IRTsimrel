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
#' reliability using `TAM::WLErel()` and `TAM::EAPrel()`. This helper does not
#' fit or score a 3PL model.
#'
#' @param resp Matrix or data.frame of item responses (0/1).
#' @param model Character. `"rasch"` or `"2pl"`. A 3PL value is deliberately
#'   unsupported because the package's validated TAM contract requires both
#'   WLE and EAP reliability, while the tested TAM 3PL path does not provide
#'   the required WLE workflow.
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
#' WLE and EAP use different estimators and variance bases, so neither is a
#' universal upper or lower bound for the other. EAP is the closer external
#' analogue to the MSEM-based population estimand, but it is not identical to
#' the information-based EQC estimand.
#'
#' ## 3PL limitation
#'
#' `compute_reliability_tam()` intentionally rejects `model = "3pl"`. The
#' package's Phase 7 external 3PL oracle used a direct TAM EAP-only fit; TAM 3PL
#' WLE was unavailable in that validated path. That oracle is evidence about
#' the probability/information implementation, not a public 3PL WLE API.
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

.irtsimrel_validate_item_design_vector <- function(value,
                                                    label,
                                                    n_items,
                                                    kind = c("beta", "lambda", "guessing")) {
  kind <- match.arg(kind)
  if (!is.numeric(value) || length(value) != n_items ||
      anyNA(value) || any(!is.finite(value))) {
    stop(
      "`", label, "` must be a finite numeric vector of length ",
      n_items, "."
    )
  }
  value <- as.numeric(value)

  if (kind == "lambda" && any(value <= 0)) {
    stop("`", label, "` must contain positive values.")
  }
  if (kind == "guessing" && any(value < 0 | value >= 1)) {
    stop("`", label, "` values must satisfy 0 <= guessing < 1.")
  }

  value
}

.irtsimrel_assert_item_design_match <- function(reference,
                                                 candidate,
                                                 reference_label,
                                                 candidate_label) {
  if (!isTRUE(all.equal(
    as.numeric(reference),
    as.numeric(candidate),
    tolerance = 1e-10
  ))) {
    stop(
      "`", candidate_label, "` must match `", reference_label,
      "` for response simulation."
    )
  }
  invisible(NULL)
}

# Normalize new and legacy stored item designs without modifying the input.
# Rasch/2PL objects written before guessing was stored receive an implicit
# zero vector. A 3PL result must carry guessing in at least one authoritative
# location, and all locations that are present must agree.
.irtsimrel_normalize_result_item_design <- function(result,
                                                     arg_name = "result",
                                                     operation = "response simulation") {
  result_class <- .irtsimrel_primary_result_class(result)
  if (!is.null(result$schema_version)) {
    if (!is.numeric(result$schema_version) ||
        length(result$schema_version) != 1L ||
        is.na(result$schema_version) ||
        !is.finite(result$schema_version) ||
        result$schema_version != floor(result$schema_version) ||
        result$schema_version < 1L) {
      stop("`", arg_name, "$schema_version` must be a positive integer.")
    }
    if (result$schema_version > 1L) {
      stop(
        "`", arg_name, "$schema_version` is newer than this IRTsimrel ",
        "version supports. Upgrade IRTsimrel before using this object."
      )
    }
  }
  design_fields <- c(
    "beta_vec", "lambda_base", "lambda_scaled", "items_calib"
  )

  if (identical(result_class, "spc_result") &&
      any(!design_fields %in% names(result))) {
    stop(
      "`", arg_name,
      "` is a legacy 'spc_result' without the stored item design required for ",
      operation, ". Rerun calibration with `sac_calibrate()` using the ",
      "current IRTsimrel version."
    )
  }

  .irtsimrel_missing_required_fields(
    result,
    design_fields,
    arg_name,
    operation
  )

  n_items <- .irtsimrel_validate_positive_integer_scalar(
    result$n_items,
    paste0(arg_name, "$n_items")
  )
  if (is.numeric(result$beta_vec) && is.numeric(result$lambda_scaled) &&
      length(result$beta_vec) != length(result$lambda_scaled)) {
    stop(
      "`", arg_name, "$beta_vec` and `", arg_name,
      "$lambda_scaled` must have the same length."
    )
  }
  if (is.numeric(result$beta_vec) && is.numeric(result$lambda_scaled) &&
      (n_items != length(result$beta_vec) ||
       n_items != length(result$lambda_scaled))) {
    stop(
      "`", arg_name, "$n_items` must match the length of `",
      arg_name, "$beta_vec` and `", arg_name, "$lambda_scaled`."
    )
  }
  beta <- .irtsimrel_validate_item_design_vector(
    result$beta_vec,
    paste0(arg_name, "$beta_vec"),
    n_items,
    "beta"
  )
  lambda <- .irtsimrel_validate_item_design_vector(
    result$lambda_scaled,
    paste0(arg_name, "$lambda_scaled"),
    n_items,
    "lambda"
  )
  lambda_base <- .irtsimrel_validate_item_design_vector(
    result$lambda_base,
    paste0(arg_name, "$lambda_base"),
    n_items,
    "lambda"
  )

  .irtsimrel_assert_item_design_match(
    lambda_base * result$c_star,
    lambda,
    paste0(arg_name, "$lambda_base * ", arg_name, "$c_star"),
    paste0(arg_name, "$lambda_scaled")
  )

  if (!is.list(result$items_calib) ||
      !is.data.frame(result$items_calib$data)) {
    stop(
      "`", arg_name,
      "$items_calib$data` must be a data frame containing calibrated item ",
      "parameters for response simulation."
    )
  }
  item_data <- result$items_calib$data
  if (!all(c("beta", "lambda") %in% names(item_data))) {
    stop(
      "`", arg_name,
      "$items_calib$data` must contain `beta` and `lambda` for response ",
      "simulation provenance."
    )
  }
  item_beta <- .irtsimrel_validate_item_design_vector(
    item_data$beta,
    paste0(arg_name, "$items_calib$data$beta"),
    n_items,
    "beta"
  )
  item_lambda <- .irtsimrel_validate_item_design_vector(
    item_data$lambda,
    paste0(arg_name, "$items_calib$data$lambda"),
    n_items,
    "lambda"
  )
  .irtsimrel_assert_item_design_match(
    beta,
    item_beta,
    paste0(arg_name, "$beta_vec"),
    paste0(arg_name, "$items_calib$data$beta")
  )
  .irtsimrel_assert_item_design_match(
    lambda,
    item_lambda,
    paste0(arg_name, "$lambda_scaled"),
    paste0(arg_name, "$items_calib$data$lambda")
  )

  beta_sources <- c("beta_vec", "items_calib$data$beta")
  lambda_sources <- c(
    "lambda_scaled", "lambda_base * c_star", "items_calib$data$lambda"
  )

  if (!is.null(result$items_base)) {
    if (!is.list(result$items_base) ||
        !is.data.frame(result$items_base$data)) {
      stop("`", arg_name, "$items_base$data` must be a data frame when present.")
    }
    base_data <- result$items_base$data
    if ("beta" %in% names(base_data)) {
      base_beta <- .irtsimrel_validate_item_design_vector(
        base_data$beta,
        paste0(arg_name, "$items_base$data$beta"),
        n_items,
        "beta"
      )
      .irtsimrel_assert_item_design_match(
        beta,
        base_beta,
        paste0(arg_name, "$beta_vec"),
        paste0(arg_name, "$items_base$data$beta")
      )
      beta_sources <- c(beta_sources, "items_base$data$beta")
    }
    base_lambda_name <- if ("lambda_unscaled" %in% names(base_data)) {
      "lambda_unscaled"
    } else if ("lambda" %in% names(base_data)) {
      "lambda"
    } else {
      NULL
    }
    if (!is.null(base_lambda_name)) {
      base_lambda <- .irtsimrel_validate_item_design_vector(
        base_data[[base_lambda_name]],
        paste0(arg_name, "$items_base$data$", base_lambda_name),
        n_items,
        "lambda"
      )
      .irtsimrel_assert_item_design_match(
        lambda_base,
        base_lambda,
        paste0(arg_name, "$lambda_base"),
        paste0(arg_name, "$items_base$data$", base_lambda_name)
      )
      lambda_sources <- c(
        lambda_sources,
        paste0("items_base$data$", base_lambda_name)
      )
    }
  }

  guessing_candidates <- list()
  if (!is.null(result$guessing_vec)) {
    guessing_candidates[["guessing_vec"]] <- result$guessing_vec
  }
  if ("guessing" %in% names(item_data)) {
    guessing_candidates[["items_calib$data$guessing"]] <- item_data$guessing
  }
  if (!is.null(result$items_base) &&
      is.data.frame(result$items_base$data) &&
      "guessing" %in% names(result$items_base$data)) {
    guessing_candidates[["items_base$data$guessing"]] <-
      result$items_base$data$guessing
  }

  if (length(guessing_candidates) == 0L) {
    if (identical(result$model, "3pl")) {
      stop(
        "`", arg_name,
        "` is a 3PL result but has no stored guessing parameters. Provide ",
        "`", arg_name, "$guessing_vec` or `", arg_name,
        "$items_calib$data$guessing`, or rerun calibration."
      )
    }
    guessing <- rep(0, n_items)
    guessing_sources <- "implicit_zero_for_legacy_rasch_2pl"
  } else {
    candidate_names <- names(guessing_candidates)
    guessing_candidates <- lapply(
      candidate_names,
      function(source_name) {
        .irtsimrel_validate_item_design_vector(
          guessing_candidates[[source_name]],
          paste0(arg_name, "$", source_name),
          n_items,
          "guessing"
        )
      }
    )
    names(guessing_candidates) <- candidate_names
    guessing <- guessing_candidates[[1L]]
    if (length(guessing_candidates) > 1L) {
      for (source_name in names(guessing_candidates)[-1L]) {
        .irtsimrel_assert_item_design_match(
          guessing,
          guessing_candidates[[source_name]],
          paste0(arg_name, "$", names(guessing_candidates)[1L]),
          paste0(arg_name, "$", source_name)
        )
      }
    }
    guessing_sources <- names(guessing_candidates)
  }

  if (!identical(result$model, "3pl") && any(guessing != 0)) {
    stop(
      "`", arg_name,
      "` stores nonzero guessing parameters but model is '", result$model,
      "'. Rasch and 2PL designs require zero guessing."
    )
  }

  # The response bank stores already-calibrated slopes and therefore uses
  # scale = 1. A calibrated Rasch form is represented internally as a 2PL bank
  # with a common effective discrimination.
  bank_model <- if (identical(result$model, "3pl")) "3pl" else "2pl"
  bank <- .irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda,
    guessing = guessing,
    model = bank_model
  )

  list(
    beta = beta,
    lambda = lambda,
    guessing = guessing,
    n_items = n_items,
    bank = bank,
    parameter_provenance = list(
      beta = beta_sources,
      lambda = lambda_sources,
      guessing = guessing_sources
    )
  )
}

.irtsimrel_validate_result_design <- function(result, arg_name, operation) {
  result_class <- .irtsimrel_primary_result_class(result)
  .irtsimrel_missing_required_fields(
    result,
    c(
      "c_star", "target_rho", "achieved_rho", "metric", "model", "n_items",
      "theta_var", "call"
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
      is.na(result$model) || !result$model %in% c("rasch", "2pl", "3pl")) {
    stop("`", arg_name, "$model` must be 'rasch', '2pl', or '3pl'.")
  }
  n_items <- .irtsimrel_validate_positive_integer_scalar(
    result$n_items,
    paste0(arg_name, "$n_items")
  )

  # Normalize the item design before status provenance so literal v0.1 SPC
  # objects receive the actionable rerun-calibration diagnostic even though
  # they also predate the current calibration-status schema.
  design <- .irtsimrel_normalize_result_item_design(result, arg_name, operation)

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

  design
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
      is.na(result$model) || !result$model %in% c("rasch", "2pl", "3pl")) {
    stop("`", arg_name, "$model` must be 'rasch', '2pl', or '3pl'.")
  }
  .irtsimrel_validate_positive_integer_scalar(result$n_items, paste0(arg_name, "$n_items"))

  invisible(result)
}

# Comparison-contract helpers -------------------------------------------------

# Comparability is deliberately based on a small, versioned semantic contract
# instead of an opaque object hash.  This keeps result objects readable and
# lets us explain precisely why two calibrations do (or do not) target the same
# design.  Named-list order is cosmetic, while functions and environments are
# intentionally treated as unstable signatures.
.irtsimrel_signature_value <- function(x) {
  if (is.null(x)) {
    return(list(value = NULL, stable = TRUE))
  }
  if (is.function(x) || is.environment(x) ||
      typeof(x) %in% c("externalptr", "weakref")) {
    return(list(value = NULL, stable = FALSE))
  }
  if (is.language(x)) {
    return(list(value = paste(deparse(x), collapse = " "), stable = TRUE))
  }
  if (is.list(x)) {
    values <- lapply(x, .irtsimrel_signature_value)
    if (any(!vapply(values, `[[`, logical(1), "stable"))) {
      return(list(value = NULL, stable = FALSE))
    }
    out <- lapply(values, `[[`, "value")
    if (!is.null(names(out)) && all(nzchar(names(out)))) {
      out <- out[order(names(out))]
    }
    return(list(value = out, stable = TRUE))
  }
  if (is.atomic(x)) {
    return(list(value = x, stable = TRUE))
  }
  list(value = NULL, stable = FALSE)
}

.irtsimrel_signature_equal <- function(x, y) {
  x_norm <- .irtsimrel_signature_value(x)
  y_norm <- .irtsimrel_signature_value(y)
  isTRUE(x_norm$stable) && isTRUE(y_norm$stable) &&
    identical(x_norm$value, y_norm$value)
}

.irtsimrel_first_contract_value <- function(...) {
  candidates <- list(...)
  for (candidate in candidates) {
    if (!is.null(candidate)) return(candidate)
  }
  NULL
}

.irtsimrel_call_argument <- function(result, argument) {
  result_call <- result$call
  if (is.null(result_call) || !is.call(result_call)) return(NULL)
  call_list <- as.list(result_call)
  if (!argument %in% names(call_list)) return(NULL)
  call_list[[argument]]
}

.irtsimrel_literal_call_value <- function(x) {
  if (is.null(x)) return(list(value = NULL, stable = FALSE))
  if (is.atomic(x)) return(list(value = x, stable = TRUE))
  if (is.call(x) && identical(x[[1L]], quote(list))) {
    values <- as.list(x)[-1L]
    normalized <- lapply(values, .irtsimrel_literal_call_value)
    if (any(!vapply(normalized, `[[`, logical(1), "stable"))) {
      return(list(value = NULL, stable = FALSE))
    }
    out <- lapply(normalized, `[[`, "value")
    names(out) <- names(values)
    return(list(value = out, stable = TRUE))
  }
  list(value = NULL, stable = FALSE)
}

.irtsimrel_contract_item_scope <- function(result, design_signature) {
  scope <- .irtsimrel_first_contract_value(
    result$item_scope,
    design_signature$item_scope,
    if (is.list(result$calibration_design)) result$calibration_design$item_scope else NULL
  )
  if (!is.null(scope)) return(scope)

  # Safe adapters for schema-complete pre-v0.3 SAC objects.  A literal v0.1
  # SPC object has neither stored item parameters nor this provenance and is
  # therefore left unresolved rather than guessed.
  if (inherits(result, "eqc_result")) return("fixed_form")
  if (identical(result$item_design, "fixed_iteration_items")) {
    return("fixed_form")
  }
  if (identical(result$item_design, "post_calibration_draw")) {
    return("item_superpopulation")
  }
  resample_arg <- .irtsimrel_literal_call_value(
    .irtsimrel_call_argument(result, "resample_items")
  )
  if (isTRUE(resample_arg$stable) && is.logical(resample_arg$value) &&
      length(resample_arg$value) == 1L && !is.na(resample_arg$value)) {
    return(if (resample_arg$value) "item_superpopulation" else "fixed_form")
  }
  NULL
}

.irtsimrel_contract_latent_spec <- function(result, design_signature) {
  canonicalize <- function(specification) {
    if (!is.list(specification) || is.null(specification$shape)) {
      return(specification)
    }
    shape <- tryCatch(
      .irtsimrel_match_latent_shape(specification$shape),
      error = function(e) NULL
    )
    if (is.null(shape)) return(NULL)
    specification$shape <- shape
    specification
  }

  embedded <- .irtsimrel_first_contract_value(
    design_signature$latent_specification,
    design_signature$latent_spec,
    if (is.list(result$calibration_design)) {
      result$calibration_design$latent_specification
    } else NULL
  )
  if (!is.null(embedded)) return(canonicalize(embedded))

  shape <- .irtsimrel_first_contract_value(
    design_signature$latent_shape,
    result$latent_shape,
    if (is.list(result$calibration_design)) result$calibration_design$latent_shape else NULL
  )
  params <- .irtsimrel_first_contract_value(
    design_signature$latent_params,
    result$latent_params,
    if (is.list(result$calibration_design)) result$calibration_design$latent_params else NULL
  )

  # match.call() stores only supplied arguments.  Omission is safe to map to
  # the public defaults; a symbol or other unevaluated expression is not.
  if (is.null(shape) && !is.null(result$call) && is.call(result$call)) {
    shape_call <- .irtsimrel_call_argument(result, "latent_shape")
    if (is.null(shape_call)) {
      shape <- "normal"
    } else {
      shape_literal <- .irtsimrel_literal_call_value(shape_call)
      if (isTRUE(shape_literal$stable)) shape <- shape_literal$value
    }
  }
  if (is.null(params) && !is.null(result$call) && is.call(result$call)) {
    params_call <- .irtsimrel_call_argument(result, "latent_params")
    if (is.null(params_call)) {
      params <- list()
    } else {
      params_literal <- .irtsimrel_literal_call_value(params_call)
      if (isTRUE(params_literal$stable)) params <- params_literal$value
    }
  }

  if (is.null(shape) || is.null(params)) return(NULL)
  canonicalize(list(shape = shape, params = params))
}

.irtsimrel_contract_scale_convention <- function(result, design_signature) {
  explicit <- .irtsimrel_first_contract_value(
    design_signature$scale_convention,
    result$scale_convention,
    if (is.list(result$calibration_design)) {
      result$calibration_design$scale_convention
    } else NULL
  )
  if (!is.null(explicit)) {
    # Accept the canonical string directly.  A structured representation with
    # the same semantics is normalized so signature format does not determine
    # comparability.
    if (identical(explicit, "global_discrimination_multiplier_D1")) {
      return(explicit)
    }
    if (is.list(explicit)) {
      d_value <- .irtsimrel_first_contract_value(explicit$D, explicit$d)
      target <- .irtsimrel_first_contract_value(
        explicit$target,
        explicit$scaling_target,
        explicit$parameter
      )
      if (is.numeric(d_value) && length(d_value) == 1L && d_value == 1 &&
          is.character(target) && length(target) == 1L &&
          target %in% c("lambda", "discrimination")) {
        return("global_discrimination_multiplier_D1")
      }
    }
    return(explicit)
  }

  # All package-produced EQC/SAC objects use D = 1 and multiply only the base
  # discrimination by c.  This is a safe compatibility inference from the
  # stored lambda_base/lambda_scaled relationship, not a guess about custom
  # calibration code.
  if (inherits(result, c("eqc_result", "sac_result", "spc_result")) &&
      is.numeric(result$lambda_base) && is.numeric(result$lambda_scaled) &&
      is.numeric(result$c_star) && length(result$c_star) == 1L &&
      isTRUE(all.equal(
        as.numeric(result$lambda_scaled),
        as.numeric(result$lambda_base) * result$c_star,
        tolerance = 1e-10
      ))) {
    return("global_discrimination_multiplier_D1")
  }
  NULL
}

.irtsimrel_contract_fixed_bank <- function(result, design_signature) {
  beta <- .irtsimrel_first_contract_value(
    design_signature$beta,
    design_signature$beta_vec,
    result$beta_vec
  )
  lambda <- .irtsimrel_first_contract_value(
    design_signature$lambda_base,
    result$lambda_base
  )
  guessing <- .irtsimrel_first_contract_value(
    design_signature$guessing,
    design_signature$guessing_vec,
    result$guessing_vec
  )
  if (is.null(guessing) && !identical(result$model, "3pl") &&
      !is.null(result$n_items)) {
    guessing <- rep(0, as.integer(result$n_items))
  }
  if (is.null(beta) || is.null(lambda) || is.null(guessing)) return(NULL)
  list(
    beta = as.numeric(beta),
    lambda_base = as.numeric(lambda),
    guessing = as.numeric(guessing)
  )
}

.irtsimrel_contract_generator <- function(result, design_signature) {
  nested <- .irtsimrel_first_contract_value(
    design_signature$item_generator,
    design_signature$generator_spec,
    design_signature$item_generator_spec,
    if (is.list(result$calibration_design)) result$calibration_design$item_generator else NULL,
    if (is.list(result$calibration_design)) result$calibration_design$generator_spec else NULL,
    result$item_generator_spec
  )
  if (!is.null(nested)) return(nested)

  item_source <- .irtsimrel_first_contract_value(
    design_signature$item_source,
    result$item_source,
    if (is.list(result$calibration_design)) result$calibration_design$item_source else NULL
  )
  item_params <- .irtsimrel_first_contract_value(
    design_signature$item_params,
    result$item_params,
    if (is.list(result$calibration_design)) result$calibration_design$item_params else NULL
  )
  if (is.null(item_source) || is.null(item_params)) return(NULL)
  list(item_source = item_source, item_params = item_params)
}

.irtsimrel_comparison_contract <- function(result) {
  design_signature <- if (is.list(result$design_signature)) {
    result$design_signature
  } else {
    list()
  }
  estimand_signature <- if (is.list(result$estimand_signature)) {
    result$estimand_signature
  } else {
    list()
  }
  item_scope <- .irtsimrel_contract_item_scope(result, design_signature)
  theta_measure <- .irtsimrel_first_contract_value(
    estimand_signature$theta_measure,
    design_signature$theta_measure,
    if (is.list(result$calibration_design)) result$calibration_design$theta_measure else NULL
  )
  if (is.null(theta_measure) && !is.null(result$theta_quad) &&
      !is.null(result$theta_var)) {
    theta_measure <- "population"
  }

  signature_state <- .irtsimrel_signature_value(result$design_signature)
  generator <- if (identical(item_scope, "item_superpopulation")) {
    .irtsimrel_contract_generator(result, design_signature)
  } else {
    NULL
  }
  generator_state <- .irtsimrel_signature_value(generator)

  list(
    model = result$model,
    metric = result$metric,
    n_items = as.integer(result$n_items),
    scale_convention = .irtsimrel_contract_scale_convention(result, design_signature),
    theta_measure = theta_measure,
    latent_specification = .irtsimrel_contract_latent_spec(result, design_signature),
    item_scope = item_scope,
    fixed_bank = if (identical(item_scope, "fixed_form")) {
      .irtsimrel_contract_fixed_bank(result, design_signature)
    } else NULL,
    generator = generator,
    generator_stable = isTRUE(generator_state$stable),
    signature_present = is.list(result$design_signature),
    signature_stable = isTRUE(signature_state$stable),
    signature_schema = .irtsimrel_first_contract_value(
      design_signature$schema_version,
      result$schema_version
    ),
    signature_inconsistent = any(c(
      !is.null(estimand_signature$metric) &&
        !identical(estimand_signature$metric, result$metric),
      !is.null(design_signature$model) &&
        !identical(design_signature$model, result$model),
      !is.null(design_signature$metric) &&
        !identical(design_signature$metric, result$metric),
      !is.null(design_signature$n_items) &&
        !identical(as.integer(design_signature$n_items), as.integer(result$n_items)),
      !is.null(design_signature$item_scope) && !is.null(result$item_scope) &&
        !identical(design_signature$item_scope, result$item_scope),
      !is.null(design_signature$beta) && !is.null(result$beta_vec) &&
        !identical(as.numeric(design_signature$beta), as.numeric(result$beta_vec)),
      !is.null(design_signature$lambda_base) && !is.null(result$lambda_base) &&
        !identical(
          as.numeric(design_signature$lambda_base),
          as.numeric(result$lambda_base)
        ),
      !is.null(design_signature$guessing) && !is.null(result$guessing_vec) &&
        !identical(
          as.numeric(design_signature$guessing),
          as.numeric(result$guessing_vec)
        ),
      !is.null(design_signature$scale_convention) &&
        !is.null(result$scale_convention) &&
        !.irtsimrel_signature_equal(
          design_signature$scale_convention,
          result$scale_convention
        )
    ))
  )
}

.irtsimrel_comparability_reasons <- function(eqc, sac) {
  reasons <- character()
  add_reason <- function(reason) {
    if (!reason %in% reasons) reasons <<- c(reasons, reason)
  }

  compare_required <- function(field, missing_code, mismatch_code) {
    if (is.null(eqc[[field]]) || is.null(sac[[field]])) {
      add_reason(missing_code)
    } else if (!.irtsimrel_signature_equal(eqc[[field]], sac[[field]])) {
      add_reason(mismatch_code)
    }
  }

  if (!identical(eqc$model, sac$model)) add_reason("model_mismatch")
  if (!identical(eqc$metric, sac$metric)) add_reason("metric_mismatch")
  if (!identical(eqc$n_items, sac$n_items)) add_reason("n_items_mismatch")
  compare_required(
    "scale_convention", "scale_convention_missing",
    "scale_convention_mismatch"
  )
  compare_required(
    "theta_measure", "theta_measure_missing", "theta_measure_mismatch"
  )
  compare_required(
    "latent_specification", "latent_specification_missing",
    "latent_specification_mismatch"
  )
  compare_required("item_scope", "item_scope_missing", "item_scope_mismatch")

  if (isTRUE(eqc$signature_inconsistent)) {
    add_reason("eqc_design_signature_inconsistent")
  }
  if (isTRUE(sac$signature_inconsistent)) {
    add_reason("sac_design_signature_inconsistent")
  }
  if (!isTRUE(eqc$signature_stable)) add_reason("eqc_design_signature_unstable")
  if (!isTRUE(sac$signature_stable)) add_reason("sac_design_signature_unstable")
  if (!is.null(eqc$signature_schema) && !is.null(sac$signature_schema) &&
      !identical(as.integer(eqc$signature_schema), as.integer(sac$signature_schema))) {
    add_reason("design_signature_schema_mismatch")
  }

  if (!is.null(eqc$item_scope) && identical(eqc$item_scope, sac$item_scope)) {
    if (identical(eqc$item_scope, "fixed_form")) {
      if (is.null(eqc$fixed_bank) || is.null(sac$fixed_bank)) {
        add_reason("fixed_bank_signature_missing")
      } else {
        if (!identical(eqc$fixed_bank$beta, sac$fixed_bank$beta)) {
          add_reason("fixed_bank_beta_mismatch")
        }
        if (!identical(eqc$fixed_bank$lambda_base, sac$fixed_bank$lambda_base)) {
          add_reason("fixed_bank_lambda_mismatch")
        }
        if (!identical(eqc$fixed_bank$guessing, sac$fixed_bank$guessing)) {
          add_reason("fixed_bank_guessing_mismatch")
        }
      }
    } else if (identical(eqc$item_scope, "item_superpopulation")) {
      if (is.null(eqc$generator) || is.null(sac$generator)) {
        add_reason("generator_signature_missing")
      } else if (!isTRUE(eqc$generator_stable) || !isTRUE(sac$generator_stable)) {
        add_reason("generator_signature_unstable")
      } else if (!.irtsimrel_signature_equal(eqc$generator, sac$generator)) {
        add_reason("generator_signature_mismatch")
      }
    } else {
      add_reason("item_scope_unknown")
    }
  }

  # An absent explicit signature is acceptable only when every required
  # semantic component above was safely reconstructed from stored provenance.
  if (!isTRUE(eqc$signature_present) &&
      (is.null(eqc$latent_specification) || is.null(eqc$item_scope) ||
       (identical(eqc$item_scope, "fixed_form") && is.null(eqc$fixed_bank)) ||
       (identical(eqc$item_scope, "item_superpopulation") && is.null(eqc$generator)))) {
    add_reason("eqc_design_signature_missing")
  }
  if (!isTRUE(sac$signature_present) &&
      (is.null(sac$latent_specification) || is.null(sac$item_scope) ||
       (identical(sac$item_scope, "fixed_form") && is.null(sac$fixed_bank)) ||
       (identical(sac$item_scope, "item_superpopulation") && is.null(sac$generator)))) {
    add_reason("sac_design_signature_missing")
  }

  reasons
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
        "hit_upper_bound", "hit_both_bounds", "branch_lost",
        "topology_unresolved", "topology_uncertain"
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

  # v0.3 schema fields are additive.  Pre-schema Rasch/2PL objects retain the
  # lazy zero-guessing adapter in .irtsimrel_normalize_result_item_design(),
  # while objects that opt into the versioned schema must carry the complete
  # comparison contract.
  if (!is.null(object$schema_version)) {
    .irtsimrel_missing_required_fields(
      object,
      c(
        "guessing_vec", "item_scope", "estimand_signature",
        "design_signature"
      ),
      arg_name,
      operation
    )
    if (!is.character(object$item_scope) || length(object$item_scope) != 1L ||
        is.na(object$item_scope) ||
        !object$item_scope %in% c("fixed_form", "item_superpopulation")) {
      stop(
        "`", arg_name,
        "$item_scope` must be 'fixed_form' or 'item_superpopulation'."
      )
    }
    if (!is.list(object$estimand_signature)) {
      stop("`", arg_name, "$estimand_signature` must be a list.")
    }
    if (!is.list(object$design_signature)) {
      stop("`", arg_name, "$design_signature` must be a list.")
    }
    .irtsimrel_missing_required_fields(
      object$estimand_signature,
      c("metric", "theta_measure", "item_measure"),
      paste0(arg_name, "$estimand_signature"),
      operation
    )
    .irtsimrel_missing_required_fields(
      object$design_signature,
      c("schema_version", "model", "metric", "n_items", "item_scope"),
      paste0(arg_name, "$design_signature"),
      operation
    )
    if (!identical(object$estimand_signature$metric, object$metric) ||
        !identical(object$design_signature$metric, object$metric) ||
        !identical(object$design_signature$model, object$model) ||
        !identical(
          as.integer(object$design_signature$n_items),
          as.integer(object$n_items)
        ) ||
        !identical(object$design_signature$item_scope, object$item_scope)) {
      stop(
        "`", arg_name,
        "` has inconsistent top-level and signature design fields."
      )
    }
  }

  if (!is.null(object$rho_scale_trajectory) &&
      (!is.numeric(object$rho_scale_trajectory) ||
       length(object$rho_scale_trajectory) != n_iter ||
       any(!is.finite(object$rho_scale_trajectory)) ||
       any(object$rho_scale_trajectory <= 0))) {
    stop(
      "`", arg_name,
      "$rho_scale_trajectory` must be a positive finite numeric vector of ",
      "length `", arg_name, "$n_iter`."
    )
  }
  if (!is.null(object$rho_update_trajectory) &&
      (!is.numeric(object$rho_update_trajectory) ||
       length(object$rho_update_trajectory) != n_iter ||
       any(!is.finite(object$rho_update_trajectory)) ||
       any(object$rho_update_trajectory < 0 | object$rho_update_trajectory > 1))) {
    stop(
      "`", arg_name,
      "$rho_update_trajectory` must be a finite numeric vector in [0, 1] ",
      "with length `", arg_name, "$n_iter`."
    )
  }
  if (!is.null(object$evaluation_trajectory)) {
    evaluation_n <- if (is.data.frame(object$evaluation_trajectory) ||
                        is.matrix(object$evaluation_trajectory)) {
      nrow(object$evaluation_trajectory)
    } else {
      length(object$evaluation_trajectory)
    }
    if (!identical(as.integer(evaluation_n), as.integer(n_iter))) {
      stop(
        "`", arg_name,
        "$evaluation_trajectory` must contain one entry per SAC iteration."
      )
    }
  }
  if (!is.null(object$iteration_trace) &&
      (!is.data.frame(object$iteration_trace) ||
       nrow(object$iteration_trace) != n_iter)) {
    stop(
      "`", arg_name,
      "$iteration_trace` must be a data frame with `", arg_name,
      "$n_iter` rows."
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
#' The function uses the normalized item design stored in the calibration
#' result. For a fixed-form result this is the calibrated fixed form. For a
#' SAC result with \code{item_scope = "item_superpopulation"}, it is one stored
#' representative evaluation form; this function does not redraw a new item
#' form or simulate the aggregate item-superpopulation estimand.
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
#'   \item{\code{guessing}}{Item lower asymptotes (I x 1). This additive
#'     component is a zero vector for Rasch and 2PL designs.}
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
#'   resample_items = FALSE,
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
  latent_shape <- .irtsimrel_match_latent_shape(latent_shape)
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
  guessing <- design$guessing
  I <- design$n_items

  # Generate responses through the common IRT kernel. Item-column chunks keep
  # the same column-major RNG order as one vectorized rbinom() call while
  # bounding the temporary probability matrix for large simulations.
  max_probability_cells <- 1000000L
  items_per_chunk <- max(
    1L,
    min(I, as.integer(floor(max_probability_cells / n_persons)))
  )
  response_matrix <- matrix(NA_integer_, nrow = n_persons, ncol = I)
  item_starts <- seq.int(1L, I, by = items_per_chunk)

  for (item_start in item_starts) {
    item_end <- min(I, item_start + items_per_chunk - 1L)
    item_rows <- item_start:item_end
    response_bank <- .irtsimrel_item_bank(
      beta = beta[item_rows],
      lambda_base = lambda[item_rows],
      guessing = guessing[item_rows],
      model = if (identical(result$model, "3pl")) "3pl" else "2pl"
    )
    prob_mat <- .irtsimrel_probability(
      theta = theta,
      bank = response_bank,
      scale = 1,
      chunk_size = min(n_persons, max_probability_cells)
    )
    response_matrix[, item_rows] <- matrix(
      rbinom(length(prob_mat), size = 1, prob = as.vector(prob_mat)),
      nrow = n_persons,
      ncol = length(item_rows)
    )
  }

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
      latent_params = latent_params,
      schema_version = if (!is.null(result$schema_version)) {
        result$schema_version
      } else {
        0L
      },
      item_scope = if (!is.null(result$item_scope)) {
        result$item_scope
      } else {
        NA_character_
      },
      design_signature = if (!is.null(result$design_signature)) {
        result$design_signature
      } else {
        NULL
      },
      item_parameter_provenance = design$parameter_provenance
    ),
    guessing = guessing
  )
}


# =============================================================================
# Comparison Function for EQC and SAC Results
# =============================================================================

#' Compare EQC and SAC Calibration Results
#'
#' @description
#' Compares the calibration results from EQC and SAC algorithms. The target
#' reliability must match exactly. Numeric estimates are always shown side by
#' side, but agreement is evaluated only when the stored estimand and design
#' contracts match: model, reliability metric, theta population, latent
#' specification, item scope, test length, scale convention, and either the
#' complete fixed item bank or normalized superpopulation generator.
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
#'   \item{\code{agreement}}{Logical. \code{TRUE} if the objects are comparable,
#'     both calibrations are valid for agreement, and the difference is below
#'     5 percent; \code{FALSE} if evaluated and above threshold; \code{NA} when
#'     agreement is withheld.}
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
#'   \item{\code{comparable}}{Logical indicating whether both objects target
#'     the same estimand and design.}
#'   \item{\code{comparability_reasons}}{Stable reason codes explaining a
#'     non-comparable result.}
#'   \item{\code{agreement_status}, \code{agreement_reasons}}{Whether agreement
#'     was evaluated, and why it was withheld when it was not.}
#'   \item{\code{eqc_contract}, \code{sac_contract}}{Normalized estimand and
#'     design contracts used for the comparability decision.}
#' }
#'
#' @examples
#' \donttest{
#' # Run both algorithms
#' eqc_result <- eqc_calibrate(target_rho = 0.80, n_items = 25, seed = 42)
#' sac_result <- sac_calibrate(target_rho = 0.80, n_items = 25,
#'                             reliability_metric = "info",
#'                             c_init = eqc_result,
#'                             resample_items = FALSE,
#'                             seed = 42)
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

  # Target equality remains a hard precondition: unlike the other design
  # dimensions, comparing calibrations for different requested targets has no
  # useful side-by-side interpretation.
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

  eqc_contract <- .irtsimrel_comparison_contract(eqc_result)
  sac_contract <- .irtsimrel_comparison_contract(sac_result)
  comparability_reasons <- .irtsimrel_comparability_reasons(
    eqc_contract,
    sac_contract
  )
  comparable <- length(comparability_reasons) == 0L

  if (!comparable) {
    warning(
      "EQC and SAC designs are not comparable (",
      paste(comparability_reasons, collapse = ", "),
      "). Numeric estimates are returned, but agreement is NA."
    )
  }

  # Extract key values.  These legacy scalar fields remain available even
  # when design comparison fails.
  c_eqc <- eqc_result$c_star
  c_sac <- sac_result$c_star

  diff_abs <- abs(c_eqc - c_sac)
  diff_pct <- 100 * diff_abs / c_eqc
  achieved_diff_abs <- abs(eqc_result$achieved_rho - sac_result$achieved_rho)

  eqc_success <- !is.na(eqc_status) &&
    eqc_status %in% c("ok", "uniroot_success")
  sac_success <- !is.na(sac_status) && identical(sac_status, "ok")
  agreement_reasons <- character()
  if (!comparable) {
    agreement_status <- "not_comparable"
    agreement_reasons <- comparability_reasons
    agreement <- NA
  } else if (!eqc_success || !sac_success) {
    agreement_status <- "calibration_unsuccessful"
    if (!eqc_success) {
      agreement_reasons <- c(
        agreement_reasons,
        "eqc_calibration_unsuccessful"
      )
    }
    if (!sac_success) {
      agreement_reasons <- c(
        agreement_reasons,
        "sac_calibration_unsuccessful"
      )
    }
    agreement <- NA
  } else {
    agreement_status <- "evaluated"
    agreement <- diff_pct < 5
  }

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
    sac_status_flags = sac_status_flags,
    comparable = comparable,
    comparability_reasons = comparability_reasons,
    agreement_status = agreement_status,
    agreement_reasons = agreement_reasons,
    eqc_contract = eqc_contract,
    sac_contract = sac_contract
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
    message(sprintf("  Comparable           : %s", ifelse(comparable, "YES", "NO")))
    if (!comparable) {
      message(sprintf(
        "  Comparability reason : %s",
        paste(comparability_reasons, collapse = ", ")
      ))
    }
    message(sprintf(
      "  Agreement (< 5%%)    : %s",
      if (is.na(agreement)) "NOT EVALUATED" else if (agreement) "YES" else "NO"
    ))
    message(sprintf("  Agreement status     : %s", agreement_status))
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
