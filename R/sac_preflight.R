# =============================================================================
# sac_preflight.R
# =============================================================================
# Internal deterministic guard for stochastic-approximation calibration (SAC).
# The guard establishes one scientifically admissible increasing reliability
# branch before a Robbins--Monro trajectory is allowed to start.  Fixed-form
# and item-superpopulation objectives share the same theta nodes, item banks,
# and scale cache throughout the full and split scans.
# =============================================================================


.irtsimrel_sac_preflight_validate <- function(theta,
                                               theta_var,
                                               item_banks,
                                               target_rho,
                                               metric_internal,
                                               c_bounds,
                                               root_policy,
                                               root_controls,
                                               split_check,
                                               split_log_tolerance) {
  if (!is.numeric(theta) || length(theta) < 2L || anyNA(theta) ||
      any(!is.finite(theta))) {
    stop("`theta` must be a finite numeric vector with length >= 2.",
         call. = FALSE)
  }
  if (!is.numeric(theta_var) || length(theta_var) != 1L ||
      is.na(theta_var) || !is.finite(theta_var) || theta_var <= 0) {
    stop("`theta_var` must be a finite positive numeric scalar.",
         call. = FALSE)
  }
  if (!is.list(item_banks) || length(item_banks) == 0L) {
    stop("`item_banks` must be a non-empty list of normalized item banks.",
         call. = FALSE)
  }

  required_bank_fields <- c(
    "beta", "lambda_base", "guessing", "model", "n_items"
  )
  valid_bank <- vapply(
    item_banks,
    function(bank) {
      inherits(bank, "irtsimrel_item_bank") && is.list(bank) &&
        all(required_bank_fields %in% names(bank)) &&
        is.numeric(bank$beta) && is.numeric(bank$lambda_base) &&
        is.numeric(bank$guessing) &&
        length(bank$beta) > 0L &&
        length(bank$beta) == length(bank$lambda_base) &&
        length(bank$beta) == length(bank$guessing) &&
        identical(as.integer(bank$n_items), as.integer(length(bank$beta))) &&
        !anyNA(bank$beta) && !anyNA(bank$lambda_base) &&
        !anyNA(bank$guessing) && all(is.finite(bank$beta)) &&
        all(is.finite(bank$lambda_base)) && all(is.finite(bank$guessing)) &&
        all(bank$lambda_base > 0) &&
        all(bank$guessing >= 0 & bank$guessing < 1)
    },
    logical(1)
  )
  if (!all(valid_bank)) {
    invalid <- paste(which(!valid_bank), collapse = ", ")
    stop(
      "Every element of `item_banks` must be created by ",
      "`.irtsimrel_item_bank()` and remain internally valid; invalid ",
      "position(s): ", invalid, ".",
      call. = FALSE
    )
  }

  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      is.na(target_rho) || !is.finite(target_rho) ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a finite scalar strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      anyNA(c_bounds) || any(!is.finite(c_bounds)) ||
      any(c_bounds <= 0) || c_bounds[[1L]] >= c_bounds[[2L]]) {
    stop("`c_bounds` must contain finite values with 0 < lower < upper.",
         call. = FALSE)
  }
  valid_policies <- c(
    "lowest_increasing", "nearest_increasing", "lowest_any", "highest_any"
  )
  if (!is.character(root_policy) || length(root_policy) != 1L ||
      is.na(root_policy) || !root_policy %in% valid_policies) {
    stop(
      "`root_policy` must be one of ",
      paste(sprintf("'%s'", valid_policies), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.list(root_controls)) {
    stop("`root_controls` must be a named list.", call. = FALSE)
  }
  if (!is.logical(split_check) || length(split_check) != 1L ||
      is.na(split_check)) {
    stop("`split_check` must be TRUE or FALSE.", call. = FALSE)
  }
  if (split_check && length(item_banks) < 4L) {
    stop("`split_check = TRUE` requires at least four item banks.",
         call. = FALSE)
  }
  if (!is.numeric(split_log_tolerance) ||
      length(split_log_tolerance) != 1L ||
      is.na(split_log_tolerance) || !is.finite(split_log_tolerance) ||
      split_log_tolerance < 0) {
    stop("`split_log_tolerance` must be a finite nonnegative scalar.",
         call. = FALSE)
  }

  invisible(TRUE)
}


.irtsimrel_sac_abort_preflight <- function(message,
                                            subclass,
                                            topology) {
  .irtsimrel_abort_topology(
    paste0("SAC preflight: ", message),
    subclass = subclass,
    topology = topology
  )
}


.irtsimrel_sac_select_preflight_root <- function(topology,
                                                  root_policy,
                                                  c_bounds,
                                                  failure_subclass) {
  if (!identical(topology$topology_status, "resolved")) {
    .irtsimrel_sac_abort_preflight(
      paste0(
        "reliability topology is not resolved (status: ",
        topology$topology_status, "; stop reason: ",
        topology$stop_reason, ")."
      ),
      subclass = "irtsimrel_topology_uncertain",
      topology = topology
    )
  }
  if (!identical(topology$target_status, "feasible")) {
    .irtsimrel_sac_abort_preflight(
      paste0(
        "target reliability has no resolved feasible root in `c_bounds` ",
        "(target status: ", topology$target_status, ")."
      ),
      subclass = failure_subclass,
      topology = topology
    )
  }

  selection <- tryCatch(
    .irtsimrel_select_root(
      topology = topology,
      root_policy = root_policy,
      reference = if (identical(root_policy, "nearest_increasing")) 1 else
        NULL,
      allow_tangent = FALSE,
      allow_plateau = FALSE
    ),
    irtsimrel_selection_unavailable = function(error) {
      .irtsimrel_sac_abort_preflight(
        paste0(
          "no root is available under `root_policy = \"", root_policy,
          "\"` that can initialize an increasing SAC update."
        ),
        subclass = failure_subclass,
        topology = topology
      )
    }
  )

  root <- selection$root
  interior_tolerance <- max(topology$tol, 1e-8)
  interior <- root$log_c[[1L]] > log(c_bounds[[1L]]) + interior_tolerance &&
    root$log_c[[1L]] < log(c_bounds[[2L]]) - interior_tolerance
  admissible <- identical(root$kind[[1L]], "crossing") &&
    identical(root$direction[[1L]], "increasing") && interior
  if (!admissible) {
    .irtsimrel_sac_abort_preflight(
      paste0(
        "selected root is not an interior increasing crossing (kind: ",
        root$kind[[1L]], "; direction: ", root$direction[[1L]],
        "; c: ", format(root$c[[1L]], digits = 8), ")."
      ),
      subclass = failure_subclass,
      topology = selection$topology
    )
  }

  branch_id <- root$branch_id[[1L]]
  branch <- if (is.finite(branch_id)) {
    selection$topology$branches[
      selection$topology$branches$branch_id == branch_id,
      ,
      drop = FALSE
    ]
  } else {
    selection$topology$branches[FALSE, , drop = FALSE]
  }
  if (nrow(branch) != 1L ||
      !identical(branch$direction[[1L]], "increasing")) {
    .irtsimrel_sac_abort_preflight(
      "selected crossing is not attached to one resolved increasing branch.",
      subclass = failure_subclass,
      topology = selection$topology
    )
  }

  effective_bounds <- c(
    lower = max(c_bounds[[1L]], branch$lower[[1L]]),
    upper = min(c_bounds[[2L]], branch$upper[[1L]])
  )
  root_inside_branch <- selection$c > effective_bounds[["lower"]] &&
    selection$c < effective_bounds[["upper"]]
  if (!all(is.finite(effective_bounds)) ||
      effective_bounds[["lower"]] >= effective_bounds[["upper"]] ||
      !root_inside_branch) {
    .irtsimrel_sac_abort_preflight(
      "the selected branch does not yield a nonempty interior projection interval.",
      subclass = failure_subclass,
      topology = selection$topology
    )
  }

  list(
    topology = selection$topology,
    selection = selection,
    selected_root = unname(selection$c),
    selected_root_record = root,
    selected_branch = branch,
    effective_branch_bounds = effective_bounds
  )
}


#' Internal: Establish an Increasing SAC Calibration Branch
#'
#' @param theta Fixed ability nodes used at every evaluated scale.
#' @param theta_var Fixed positive theta variance.
#' @param item_banks A non-empty list of normalized item-bank objects. A
#'   one-element list represents a fixed form; multiple banks represent an
#'   item-superpopulation pilot sample.
#' @param target_rho Target reliability.
#' @param metric_internal Either `"msem"` or `"info"`.
#' @param c_bounds Positive scale-search bounds.
#' @param root_policy Root-selection policy understood by
#'   `.irtsimrel_select_root()`.
#' @param root_controls Controls passed to `.irtsimrel_scan_topology()`.
#' @param split_check Whether to require odd/even form-split stability.
#' @param split_log_tolerance Maximum absolute log ratio between the two split
#'   roots.
#'
#' @return An object of class `irtsimrel_sac_preflight`.
#' @noRd
.irtsimrel_sac_preflight <- function(theta,
                                      theta_var,
                                      item_banks,
                                      target_rho,
                                      metric_internal = c("msem", "info"),
                                      c_bounds,
                                      root_policy = "lowest_increasing",
                                      root_controls = list(),
                                      split_check = length(item_banks) >= 4L,
                                      split_log_tolerance = 0.5) {
  metric_internal <- match.arg(metric_internal)
  .irtsimrel_sac_preflight_validate(
    theta = theta,
    theta_var = theta_var,
    item_banks = item_banks,
    target_rho = target_rho,
    metric_internal = metric_internal,
    c_bounds = c_bounds,
    root_policy = root_policy,
    root_controls = root_controls,
    split_check = split_check,
    split_log_tolerance = split_log_tolerance
  )

  theta <- as.numeric(theta)
  theta_var <- as.numeric(theta_var)
  c_bounds <- as.numeric(c_bounds)
  n_forms <- length(item_banks)
  item_scope <- if (n_forms == 1L) {
    "fixed_form"
  } else {
    "item_superpopulation"
  }

  reliability_cache <- new.env(hash = TRUE, parent = emptyenv())
  request_counts <- new.env(hash = TRUE, parent = emptyenv())
  scale_values <- new.env(hash = TRUE, parent = emptyenv())
  cache_eval_count <- 0L
  cache_hit_count <- 0L

  evaluate_forms <- function(c) {
    key <- sprintf("%.17g", as.numeric(c))
    requests <- if (exists(key, envir = request_counts, inherits = FALSE)) {
      get(key, envir = request_counts, inherits = FALSE)
    } else {
      0L
    }
    assign(key, requests + 1L, envir = request_counts)

    if (exists(key, envir = reliability_cache, inherits = FALSE)) {
      cache_hit_count <<- cache_hit_count + 1L
      return(get(key, envir = reliability_cache, inherits = FALSE))
    }

    values <- vapply(
      item_banks,
      function(bank) {
        value <- tryCatch(
          .compute_rho_generic(
            c = c,
            theta_vec = theta,
            beta_vec = bank$beta,
            lambda_base = bank$lambda_base,
            theta_var = theta_var,
            metric_internal = metric_internal,
            guessing = bank$guessing
          ),
          error = function(error) NA_real_
        )
        if (!is.numeric(value) || length(value) != 1L ||
            is.na(value) || !is.finite(value)) {
          NA_real_
        } else {
          as.numeric(value)
        }
      },
      numeric(1)
    )
    cache_eval_count <<- cache_eval_count + 1L
    assign(key, values, envir = reliability_cache)
    assign(key, as.numeric(c), envir = scale_values)
    values
  }

  objective_for <- function(indices) {
    force(indices)
    function(c) {
      values <- evaluate_forms(c)[indices]
      if (length(values) == 0L || any(!is.finite(values))) {
        return(NA_real_)
      }
      mean(values)
    }
  }

  full_topology <- .irtsimrel_scan_topology(
    fn = objective_for(seq_len(n_forms)),
    c_bounds = c_bounds,
    target = target_rho,
    tol = 1e-6,
    controls = root_controls
  )
  full <- .irtsimrel_sac_select_preflight_root(
    topology = full_topology,
    root_policy = root_policy,
    c_bounds = c_bounds,
    failure_subclass = "irtsimrel_branch_unavailable"
  )

  split <- list(
    checked = FALSE,
    stable = NA,
    reason = "not_requested",
    odd_indices = integer(),
    even_indices = integer(),
    roots = c(odd = NA_real_, even = NA_real_),
    absolute_log_root_ratio = NA_real_,
    log_tolerance = split_log_tolerance,
    odd = NULL,
    even = NULL
  )

  if (split_check) {
    odd_indices <- seq.int(1L, n_forms, by = 2L)
    even_indices <- seq.int(2L, n_forms, by = 2L)
    odd_topology <- .irtsimrel_scan_topology(
      fn = objective_for(odd_indices),
      c_bounds = c_bounds,
      target = target_rho,
      tol = 1e-6,
      controls = root_controls
    )
    even_topology <- .irtsimrel_scan_topology(
      fn = objective_for(even_indices),
      c_bounds = c_bounds,
      target = target_rho,
      tol = 1e-6,
      controls = root_controls
    )

    split_context <- list(
      full = full$topology,
      odd = odd_topology,
      even = even_topology,
      split_stable = FALSE,
      split_reason = "split_root_unavailable"
    )
    class(split_context) <- c("irtsimrel_sac_split_topology", "list")

    odd <- tryCatch(
      .irtsimrel_sac_select_preflight_root(
        topology = odd_topology,
        root_policy = root_policy,
        c_bounds = c_bounds,
        failure_subclass = "irtsimrel_topology_uncertain"
      ),
      irtsimrel_topology_uncertain = function(error) {
        .irtsimrel_sac_abort_preflight(
          "odd/even split has no stable admissible odd-form root.",
          subclass = "irtsimrel_topology_uncertain",
          topology = split_context
        )
      }
    )
    even <- tryCatch(
      .irtsimrel_sac_select_preflight_root(
        topology = even_topology,
        root_policy = root_policy,
        c_bounds = c_bounds,
        failure_subclass = "irtsimrel_topology_uncertain"
      ),
      irtsimrel_topology_uncertain = function(error) {
        .irtsimrel_sac_abort_preflight(
          "odd/even split has no stable admissible even-form root.",
          subclass = "irtsimrel_topology_uncertain",
          topology = split_context
        )
      }
    )

    full_bounds <- full$effective_branch_bounds
    roots <- c(odd = odd$selected_root, even = even$selected_root)
    roots_inside_full_branch <- all(
      roots > full_bounds[["lower"]] & roots < full_bounds[["upper"]]
    )
    absolute_log_root_ratio <- abs(log(roots[["odd"]] / roots[["even"]]))
    log_stable <- is.finite(absolute_log_root_ratio) &&
      absolute_log_root_ratio <= split_log_tolerance

    if (!roots_inside_full_branch || !log_stable) {
      split_context$odd <- odd$topology
      split_context$even <- even$topology
      split_context$split_roots <- roots
      split_context$full_branch_bounds <- full_bounds
      split_context$absolute_log_root_ratio <- absolute_log_root_ratio
      split_context$split_log_tolerance <- split_log_tolerance
      split_context$split_reason <- if (!roots_inside_full_branch) {
        "split_root_outside_full_branch"
      } else {
        "split_log_root_mismatch"
      }
      .irtsimrel_sac_abort_preflight(
        paste0(
          "odd/even split roots are unstable (reason: ",
          split_context$split_reason, "; |log(c_odd/c_even)|: ",
          format(absolute_log_root_ratio, digits = 6), "; tolerance: ",
          format(split_log_tolerance, digits = 6), ")."
        ),
        subclass = "irtsimrel_topology_uncertain",
        topology = split_context
      )
    }

    split <- list(
      checked = TRUE,
      stable = TRUE,
      reason = "stable",
      odd_indices = odd_indices,
      even_indices = even_indices,
      roots = roots,
      absolute_log_root_ratio = absolute_log_root_ratio,
      log_tolerance = split_log_tolerance,
      odd = odd,
      even = even
    )
  }

  cache_keys <- ls(reliability_cache, all.names = TRUE)
  per_scale <- if (length(cache_keys) == 0L) {
    data.frame(
      c = numeric(), request_count = integer(), cache_eval_count = integer(),
      form_reliability_eval_count = integer(), stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      c = vapply(
        cache_keys,
        function(key) get(key, envir = scale_values, inherits = FALSE),
        numeric(1)
      ),
      request_count = vapply(
        cache_keys,
        function(key) get(key, envir = request_counts, inherits = FALSE),
        integer(1)
      ),
      cache_eval_count = rep.int(1L, length(cache_keys)),
      form_reliability_eval_count = rep.int(n_forms, length(cache_keys)),
      stringsAsFactors = FALSE
    )
  }
  per_scale <- per_scale[order(per_scale$c), , drop = FALSE]
  rownames(per_scale) <- NULL

  split_topologies <- list(
    odd = if (isTRUE(split$checked)) split$odd$topology else NULL,
    even = if (isTRUE(split$checked)) split$even$topology else NULL
  )
  split_selections <- list(
    odd = if (isTRUE(split$checked)) split$odd$selection else NULL,
    even = if (isTRUE(split$checked)) split$even$selection else NULL
  )
  eval_diagnostics <- list(
    cache_eval_count = cache_eval_count,
    cache_hit_count = cache_hit_count,
    form_reliability_eval_count = cache_eval_count * n_forms,
    n_forms = n_forms,
    per_scale = per_scale,
    scanner_eval_count = list(
      full = full$topology$eval_count,
      odd = if (isTRUE(split$checked)) split$odd$topology$eval_count else 0L,
      even = if (isTRUE(split$checked)) split$even$topology$eval_count else 0L
    )
  )

  result <- list(
    selected_root = full$selected_root,
    c_star = full$selected_root,
    selected_root_record = full$selected_root_record,
    selected_branch = full$selected_branch,
    selected_branch_bounds = full$effective_branch_bounds,
    effective_branch_bounds = full$effective_branch_bounds,
    topology = full$topology,
    full_topology = full$topology,
    full_selection = full$selection,
    split_topologies = split_topologies,
    split_selections = split_selections,
    split = split,
    split_stable = split$stable,
    split_reason = split$reason,
    eval_diagnostics = eval_diagnostics,
    cache_eval_count = cache_eval_count,
    item_scope = item_scope,
    n_forms = n_forms,
    metric = metric_internal,
    target_rho = target_rho,
    c_bounds = c_bounds,
    root_policy = root_policy
  )
  class(result) <- c("irtsimrel_sac_preflight", "list")
  result
}
