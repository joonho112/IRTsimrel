# =============================================================================
# stochastic_methods.R
# =============================================================================
# Canonical v0.3 SAC result methods. The deprecated alias and plot method
# remain in spc_calibrate.R; no duplicate method definitions are retained.
# =============================================================================


.irtsimrel_sac_method_scope <- function(object) {
  if (!is.null(object$item_scope)) {
    return(object$item_scope)
  }
  if (is.list(object$calibration_design) &&
      !is.null(object$calibration_design$item_scope)) {
    return(object$calibration_design$item_scope)
  }
  if (identical(object$item_design, "fixed_iteration_items")) {
    return("fixed_form")
  }
  if (identical(object$item_design, "post_calibration_draw")) {
    return("item_superpopulation")
  }
  "legacy_unspecified"
}


.irtsimrel_sac_method_distribution <- function(object) {
  distribution <- object$achieved_distribution
  if (!is.null(distribution) && !is.list(distribution)) {
    stop("`object$achieved_distribution` must be a list when present.")
  }

  achieved_se <- object$achieved_se
  if (is.null(achieved_se) && is.list(distribution)) {
    achieved_se <- distribution$se
  }
  if (!is.null(achieved_se) &&
      (!is.numeric(achieved_se) || length(achieved_se) != 1L ||
       is.nan(achieved_se) ||
       (!is.na(achieved_se) && (!is.finite(achieved_se) || achieved_se < 0)))) {
    stop("Stored achieved-reliability SE must be a non-negative finite scalar.")
  }

  n_forms <- if (is.list(distribution)) distribution$n_forms else NULL
  if (is.null(n_forms) && is.list(object$evaluation_design)) {
    n_forms <- object$evaluation_design$n_forms
  }
  if (!is.null(n_forms) &&
      (!is.numeric(n_forms) || length(n_forms) != 1L || is.na(n_forms) ||
       !is.finite(n_forms) || n_forms < 1 || n_forms != floor(n_forms))) {
    stop("Stored evaluation form count must be a positive integer scalar.")
  }

  list(
    distribution = distribution,
    achieved_se = if (is.null(achieved_se) || is.na(achieved_se)) {
      NA_real_
    } else {
      as.numeric(achieved_se)
    },
    n_forms = if (is.null(n_forms)) NA_integer_ else as.integer(n_forms)
  )
}


.irtsimrel_sac_method_branch <- function(object) {
  branch <- if (is.list(object$branch)) object$branch else list()
  preflight <- if (is.list(object$preflight)) object$preflight else list()

  selected_root <- if (!is.null(branch$selected_root)) {
    branch$selected_root
  } else {
    preflight$selected_root
  }
  if (is.null(selected_root)) selected_root <- NA_real_

  selected_branch <- if (!is.null(branch$selected_branch)) {
    branch$selected_branch
  } else {
    preflight$selected_branch
  }
  root_policy <- if (!is.null(branch$root_policy)) {
    branch$root_policy
  } else {
    preflight$root_policy
  }
  if (is.null(root_policy)) root_policy <- NA_character_

  topology <- if (is.list(preflight$topology)) {
    preflight$topology
  } else if (is.list(preflight$full_topology)) {
    preflight$full_topology
  } else {
    list()
  }
  preflight_status <- topology$topology_status
  if (is.null(preflight_status)) preflight_status <- NA_character_

  list(
    selected_root = selected_root,
    selected_branch = selected_branch,
    root_policy = root_policy,
    preflight_status = preflight_status
  )
}


.irtsimrel_sac_guessing_summary <- function(guessing) {
  c(
    min = min(guessing),
    mean = mean(guessing),
    max = max(guessing)
  )
}


.irtsimrel_sac_branch_label <- function(selected_branch) {
  if (is.null(selected_branch)) return(NA_character_)
  if (is.data.frame(selected_branch) && nrow(selected_branch) == 1L) {
    branch_id <- if ("branch_id" %in% names(selected_branch)) {
      selected_branch$branch_id[[1L]]
    } else {
      NA
    }
    direction <- if ("direction" %in% names(selected_branch)) {
      selected_branch$direction[[1L]]
    } else {
      NA_character_
    }
    bounds <- if (all(c("lower", "upper") %in% names(selected_branch))) {
      sprintf(
        "[%.4g, %.4g]",
        selected_branch$lower[[1L]],
        selected_branch$upper[[1L]]
      )
    } else {
      NA_character_
    }
    pieces <- character()
    if (!is.na(branch_id)) pieces <- c(pieces, paste0("id=", branch_id))
    if (!is.na(direction)) pieces <- c(pieces, direction)
    if (!is.na(bounds)) pieces <- c(pieces, bounds)
    if (length(pieces) > 0L) return(paste(pieces, collapse = ", "))
  }
  if (is.atomic(selected_branch) && length(selected_branch) == 1L &&
      !is.na(selected_branch)) {
    return(as.character(selected_branch))
  }
  "available"
}


#' @rdname sac_calibrate
#' @export
print.sac_result <- function(x, digits = 4, ...) {
  .irtsimrel_validate_sac_result_object(x, "x", "SAC print")
  design <- .irtsimrel_normalize_result_item_design(x, "x", "SAC print")
  distribution <- .irtsimrel_sac_method_distribution(x)
  branch <- .irtsimrel_sac_method_branch(x)
  scope <- .irtsimrel_sac_method_scope(x)

  cat("\n")
  cat("=======================================================\n")
  cat("  Stochastic Approximation Calibration (SAC) Results\n")
  cat("=======================================================\n\n")

  cat("Calibration Summary:\n")
  cat(sprintf("  Model                        : %s\n", toupper(x$model)))
  cat(sprintf("  Item scope                   : %s\n", scope))
  cat(sprintf("  Target reliability (rho*)    : %.*f\n", digits, x$target_rho))
  cat(sprintf("  Achieved reliability         : %.*f\n", digits, x$achieved_rho))
  if (is.finite(distribution$achieved_se)) {
    cat(sprintf("  Achieved reliability SE      : %.*f\n",
                digits, distribution$achieved_se))
  }
  if (!is.na(distribution$n_forms)) {
    cat(sprintf("  Evaluation forms             : %d\n", distribution$n_forms))
  }
  cat(sprintf("  Absolute error               : %.2e\n",
              abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor (c*)          : %.*f\n", digits, x$c_star))
  cat(sprintf("  Calibration status           : %s\n", x$calibration_status))
  if (is.finite(branch$selected_root)) {
    cat(sprintf("  Selected preflight root      : %.*f\n",
                digits, branch$selected_root))
  }
  branch_label <- .irtsimrel_sac_branch_label(branch$selected_branch)
  if (!is.na(branch_label)) {
    cat(sprintf("  Selected branch              : %s\n", branch_label))
  }
  if (!is.na(branch$root_policy)) {
    cat(sprintf("  Root policy                  : %s\n", branch$root_policy))
  }
  if (identical(x$model, "3pl")) {
    guessing_summary <- .irtsimrel_sac_guessing_summary(design$guessing)
    cat(sprintf(
      "  Guessing (min/mean/max)      : %.*f / %.*f / %.*f\n",
      digits, guessing_summary[["min"]],
      digits, guessing_summary[["mean"]],
      digits, guessing_summary[["max"]]
    ))
  }
  cat("\n")

  cat("Algorithm Settings:\n")
  cat(sprintf("  Number of items (I)          : %d\n", x$n_items))
  cat(sprintf("  M per iteration              : %d\n", x$M_per_iter))
  cat(sprintf("  M for variance pre-calc      : %d\n", x$M_pre))
  cat(sprintf("  Total iterations             : %d\n", x$n_iter))
  cat(sprintf("  Burn-in                      : %d\n", x$burn_in))
  cat(sprintf(
    "  Reliability metric           : %s\n",
    ifelse(
      x$metric == "info",
      "Average-information (tilde)",
      "MSEM-based (bar/w)"
    )
  ))
  cat(sprintf(
    "  Step params: a=%.2f, A=%.0f, gamma=%.2f\n",
    x$step_params$a, x$step_params$A, x$step_params$gamma
  ))
  cat("\n")

  cat("Convergence Diagnostics:\n")
  cat(sprintf("  Initialization method        : %s\n", x$init_method))
  cat(sprintf("  Initial c_0                  : %.4f\n", x$c_init))
  cat(sprintf("  Final iterate c_n            : %.4f\n", x$c_final))
  cat(sprintf("  Polyak-Ruppert c*            : %.4f\n", x$c_star))
  cat(sprintf("  Pre-calculated theta_var     : %.4f\n", x$theta_var))
  cat(sprintf(
    "  Converged                    : %s\n",
    ifelse(x$convergence$converged, "Yes", "No")
  ))
  cat(sprintf("  Post-burn-in SD              : %.4f\n",
              x$convergence$sd_post_burn))
  cat(sprintf("  Final iter gradient          : %+.4f\n",
              x$convergence$final_gradient))
  if (!is.null(x$convergence$gradient_at_c_star)) {
    cat(sprintf("  Gradient at c*               : %+.4f\n",
                x$convergence$gradient_at_c_star))
  }
  if (!is.null(x$convergence$projection_count)) {
    cat(sprintf(
      "  Projection count             : %d (%.1f%%)\n",
      x$convergence$projection_count,
      100 * x$convergence$projection_rate
    ))
  }
  if (!is.null(x$convergence$status_flags)) {
    cat(sprintf(
      "  Status flags                 : %s\n",
      paste(x$convergence$status_flags, collapse = ", ")
    ))
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
#' @export
summary.sac_result <- function(object, ...) {
  .irtsimrel_validate_sac_result_object(object, "object", "SAC summary")
  design <- .irtsimrel_normalize_result_item_design(
    object,
    "object",
    "SAC summary"
  )
  distribution <- .irtsimrel_sac_method_distribution(object)
  branch <- .irtsimrel_sac_method_branch(object)

  # The historical field order is a public compatibility contract. New
  # diagnostics are appended only after this exact prefix.
  out <- list(
    c_star = object$c_star,
    target_rho = object$target_rho,
    achieved_rho = object$achieved_rho,
    metric = object$metric,
    model = object$model,
    n_items = object$n_items,
    n_iter = object$n_iter,
    M_per_iter = object$M_per_iter,
    M_pre = object$M_pre,
    theta_var = object$theta_var,
    init_method = object$init_method,
    convergence = object$convergence,
    burn_in = object$burn_in,
    item_scope = .irtsimrel_sac_method_scope(object),
    achieved_se = distribution$achieved_se,
    achieved_distribution_n = distribution$n_forms,
    achieved_distribution = distribution$distribution,
    calibration_status = object$calibration_status,
    root_policy = branch$root_policy,
    selected_root = branch$selected_root,
    selected_branch = branch$selected_branch,
    preflight_status = branch$preflight_status
  )
  if (identical(object$model, "3pl")) {
    out$guessing_summary <- .irtsimrel_sac_guessing_summary(design$guessing)
  }
  class(out) <- "summary.sac_result"
  out
}


#' @rdname sac_calibrate
#' @export
print.summary.sac_result <- function(x, digits = 4, ...) {
  legacy_fields <- c(
    "c_star", "target_rho", "achieved_rho", "metric", "model",
    "n_items", "n_iter", "M_per_iter", "M_pre", "theta_var",
    "init_method", "convergence", "burn_in"
  )
  .irtsimrel_validate_summary_object(
    x,
    "summary.sac_result",
    legacy_fields,
    "x",
    "SAC summary print"
  )

  cat("Summary: Stochastic Approximation Calibration (SAC)\n")
  cat("====================================================\n")
  cat(sprintf("  Model            : %s\n", toupper(x$model)))
  cat(sprintf(
    "  Metric           : %s\n",
    ifelse(
      x$metric == "info",
      "Average-information (tilde)",
      "MSEM-based (bar/w)"
    )
  ))
  if (!is.null(x$item_scope)) {
    cat(sprintf("  Item scope       : %s\n", x$item_scope))
  }
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
  if (!is.null(x$achieved_se) && is.finite(x$achieved_se)) {
    cat(sprintf("  Achieved SE      : %.*f\n", digits, x$achieved_se))
  }
  if (!is.null(x$achieved_distribution_n) &&
      !is.na(x$achieved_distribution_n)) {
    cat(sprintf("  Evaluation forms : %d\n", x$achieved_distribution_n))
  }
  cat(sprintf("  Absolute error   : %.2e\n",
              abs(x$achieved_rho - x$target_rho)))
  cat(sprintf("  Scaling factor c*: %.*f\n", digits, x$c_star))
  if (!is.null(x$calibration_status)) {
    cat(sprintf("  Status           : %s\n", x$calibration_status))
  }
  if (!is.null(x$selected_root) && is.finite(x$selected_root)) {
    cat(sprintf("  Selected root    : %.*f\n", digits, x$selected_root))
  }
  branch_label <- .irtsimrel_sac_branch_label(x$selected_branch)
  if (!is.na(branch_label)) {
    cat(sprintf("  Selected branch  : %s\n", branch_label))
  }
  if (!is.null(x$root_policy) && !is.na(x$root_policy)) {
    cat(sprintf("  Root policy      : %s\n", x$root_policy))
  }
  if (!is.null(x$guessing_summary)) {
    cat(sprintf(
      "  Guessing min/mean/max: %.*f / %.*f / %.*f\n",
      digits, x$guessing_summary[["min"]],
      digits, x$guessing_summary[["mean"]],
      digits, x$guessing_summary[["max"]]
    ))
  }
  cat(sprintf(
    "  Converged        : %s\n",
    ifelse(x$convergence$converged, "Yes", "No")
  ))
  cat(sprintf("  Post-burn-in SD  : %.*f\n",
              digits, x$convergence$sd_post_burn))
  invisible(x)
}


#' @rdname sac_calibrate
#' @export
coef.sac_result <- function(object, ...) {
  .irtsimrel_validate_sac_result_object(
    object,
    "object",
    "SAC coefficient extraction"
  )
  design <- .irtsimrel_normalize_result_item_design(
    object,
    "object",
    "SAC coefficient extraction"
  )

  out <- data.frame(
    item_id = seq_len(object$n_items),
    beta = design$beta,
    lambda_base = object$lambda_base,
    lambda_scaled = design$lambda,
    c_star = rep(object$c_star, object$n_items)
  )
  if (identical(object$model, "3pl")) {
    out$guessing <- design$guessing
  }
  out
}


.irtsimrel_sac_predict_fixed <- function(object,
                                         design,
                                         c_values,
                                         theta_vec) {
  if (is.null(theta_vec)) {
    theta_vec <- object$theta_quad
    theta_var <- object$theta_var
  } else {
    theta_vec <- .irtsimrel_validate_theta_vector(theta_vec)
    theta_var <- stats::var(theta_vec)
  }

  vapply(
    c_values,
    function(c_value) {
      .compute_rho_generic(
        c = c_value,
        theta_vec = theta_vec,
        beta_vec = design$beta,
        lambda_base = object$lambda_base,
        theta_var = theta_var,
        metric_internal = object$metric,
        guessing = design$guessing
      )
    },
    numeric(1)
  )
}


.irtsimrel_sac_predict_superpopulation <- function(object,
                                                   c_values,
                                                   theta_vec) {
  evaluation_design <- object$evaluation_design
  if (!is.list(evaluation_design)) {
    stop(
      "`object$evaluation_design` is required for item-superpopulation prediction."
    )
  }
  item_banks <- evaluation_design$item_banks
  if (!is.list(item_banks) || length(item_banks) == 0L ||
      any(!vapply(item_banks, inherits, logical(1), "irtsimrel_item_bank"))) {
    stop(
      "`object$evaluation_design$item_banks` must be a non-empty list of normalized item banks."
    )
  }
  if (any(vapply(item_banks, function(bank) {
    !identical(as.integer(bank$n_items), as.integer(object$n_items))
  }, logical(1)))) {
    stop(
      "Every stored evaluation item bank must match `object$n_items`."
    )
  }
  if (!is.null(evaluation_design$n_forms) &&
      !identical(
        as.integer(evaluation_design$n_forms),
        as.integer(length(item_banks))
      )) {
    stop(
      "`object$evaluation_design$n_forms` must match the number of stored item banks."
    )
  }

  if (is.null(theta_vec)) {
    theta_blocks <- evaluation_design$theta_blocks
    if (!is.list(theta_blocks) || length(theta_blocks) != length(item_banks)) {
      stop(
        "`object$evaluation_design$theta_blocks` must correspond one-to-one with stored item banks."
      )
    }
    theta_blocks <- lapply(
      seq_along(theta_blocks),
      function(index) {
        .irtsimrel_validate_theta_vector(
          theta_blocks[[index]],
          paste0("object$evaluation_design$theta_blocks[[", index, "]]" )
        )
      }
    )
    theta_vars <- rep(object$theta_var, length(theta_blocks))
  } else {
    theta_vec <- .irtsimrel_validate_theta_vector(theta_vec)
    theta_blocks <- rep(list(theta_vec), length(item_banks))
    theta_vars <- rep(stats::var(theta_vec), length(item_banks))
  }

  out <- vapply(
    c_values,
    function(c_value) {
      mean(vapply(
        seq_along(item_banks),
        function(index) {
          bank <- item_banks[[index]]
          .compute_rho_generic(
            c = c_value,
            theta_vec = theta_blocks[[index]],
            beta_vec = bank$beta,
            lambda_base = bank$lambda_base,
            theta_var = theta_vars[[index]],
            metric_internal = object$metric,
            guessing = bank$guessing
          )
        },
        numeric(1)
      ))
    },
    numeric(1)
  )
  attr(out, "prediction_scope") <- "item_superpopulation"
  out
}


#' @rdname sac_calibrate
#' @export
predict.sac_result <- function(object, newdata = NULL, theta_vec = NULL, ...) {
  .irtsimrel_validate_sac_result_object(object, "object", "SAC prediction")
  design <- .irtsimrel_normalize_result_item_design(
    object,
    "object",
    "SAC prediction"
  )

  if (is.null(newdata)) {
    return(object$achieved_rho)
  }
  c_values <- .irtsimrel_validate_c_values(newdata)

  result <- if (identical(object$item_scope, "item_superpopulation")) {
    .irtsimrel_sac_predict_superpopulation(object, c_values, theta_vec)
  } else {
    .irtsimrel_sac_predict_fixed(object, design, c_values, theta_vec)
  }
  names(result) <- paste0("c=", format(c_values, digits = 4))
  result
}
