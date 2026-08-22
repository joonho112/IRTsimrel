# =============================================================================
# stochastic_calibration.R
# =============================================================================
# v0.3 SAC implementation.  The single public sac_calibrate() definition lives
# here.  The deprecated SPC forwarding alias and SAC plot method remain in
# spc_calibrate.R; the other SAC S3 methods live in stochastic_methods.R.
# =============================================================================


.irtsimrel_sac_abort <- function(message, subclass, diagnostics = NULL) {
  condition <- structure(
    list(message = message, call = NULL, diagnostics = diagnostics),
    class = c(subclass, "irtsimrel_sac_error", "error", "condition")
  )
  stop(condition)
}


.irtsimrel_sac_named_list <- function(x, name) {
  if (!is.list(x)) {
    stop("`", name, "` must be a named list.", call. = FALSE)
  }
  if (length(x) > 0L &&
      (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop("Every entry in `", name, "` must be named.", call. = FALSE)
  }
  x
}


.irtsimrel_sac_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum) {
    stop("`", name, "` must be an integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(x)
}


.irtsimrel_sac_scalar <- function(x,
                                  name,
                                  lower = -Inf,
                                  upper = Inf,
                                  lower_open = FALSE,
                                  upper_open = FALSE) {
  valid <- is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
  if (valid) {
    valid <- if (lower_open) x > lower else x >= lower
  }
  if (valid) {
    valid <- if (upper_open) x < upper else x <= upper
  }
  if (!valid) {
    stop("`", name, "` is outside its permitted finite range.", call. = FALSE)
  }
  as.numeric(x)
}


.irtsimrel_sac_controls <- function(preflight_controls,
                                    evaluation_controls,
                                    M_pre,
                                    M_per_iter,
                                    resample_items) {
  preflight_controls <- .irtsimrel_sac_named_list(
    preflight_controls,
    "preflight_controls"
  )
  evaluation_controls <- .irtsimrel_sac_named_list(
    evaluation_controls,
    "evaluation_controls"
  )

  preflight_defaults <- list(
    n_forms = 8L,
    M = min(M_pre, max(500L, M_per_iter)),
    split_check = isTRUE(resample_items),
    split_log_tolerance = log(1.5),
    root_controls = list(),
    branch_margin = 0.01,
    max_consecutive_hits = 5L,
    max_post_burn_hit_rate = 0.10
  )
  evaluation_defaults <- list(
    n_forms = 20L,
    M = max(1000L, M_per_iter),
    conf_level = 0.95,
    probs = c(0.025, 0.5, 0.975)
  )

  unknown_preflight <- setdiff(
    names(preflight_controls),
    names(preflight_defaults)
  )
  if (length(unknown_preflight) > 0L) {
    stop(
      "Unknown `preflight_controls` field(s): ",
      paste(unknown_preflight, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  unknown_evaluation <- setdiff(
    names(evaluation_controls),
    names(evaluation_defaults)
  )
  if (length(unknown_evaluation) > 0L) {
    stop(
      "Unknown `evaluation_controls` field(s): ",
      paste(unknown_evaluation, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  preflight <- modifyList(preflight_defaults, preflight_controls)
  evaluation <- modifyList(evaluation_defaults, evaluation_controls)

  preflight$n_forms <- .irtsimrel_sac_integer(
    preflight$n_forms,
    "preflight_controls$n_forms",
    if (resample_items) 2L else 1L
  )
  preflight$M <- .irtsimrel_sac_integer(
    preflight$M,
    "preflight_controls$M",
    2L
  )
  if (!is.logical(preflight$split_check) ||
      length(preflight$split_check) != 1L || is.na(preflight$split_check)) {
    stop("`preflight_controls$split_check` must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (preflight$split_check && preflight$n_forms < 4L) {
    stop(
      "`preflight_controls$n_forms` must be at least 4 when split checking.",
      call. = FALSE
    )
  }
  preflight$split_log_tolerance <- .irtsimrel_sac_scalar(
    preflight$split_log_tolerance,
    "preflight_controls$split_log_tolerance",
    lower = 0
  )
  preflight$root_controls <- .irtsimrel_sac_named_list(
    preflight$root_controls,
    "preflight_controls$root_controls"
  )
  preflight$branch_margin <- .irtsimrel_sac_scalar(
    preflight$branch_margin,
    "preflight_controls$branch_margin",
    lower = 0,
    upper = 0.49
  )
  preflight$max_consecutive_hits <- .irtsimrel_sac_integer(
    preflight$max_consecutive_hits,
    "preflight_controls$max_consecutive_hits",
    1L
  )
  preflight$max_post_burn_hit_rate <- .irtsimrel_sac_scalar(
    preflight$max_post_burn_hit_rate,
    "preflight_controls$max_post_burn_hit_rate",
    lower = 0,
    upper = 1
  )

  evaluation$n_forms <- .irtsimrel_sac_integer(
    evaluation$n_forms,
    "evaluation_controls$n_forms",
    1L
  )
  evaluation$M <- .irtsimrel_sac_integer(
    evaluation$M,
    "evaluation_controls$M",
    2L
  )
  evaluation$conf_level <- .irtsimrel_sac_scalar(
    evaluation$conf_level,
    "evaluation_controls$conf_level",
    lower = 0,
    upper = 1,
    lower_open = TRUE,
    upper_open = TRUE
  )
  if (!is.numeric(evaluation$probs) || length(evaluation$probs) == 0L ||
      anyNA(evaluation$probs) || any(!is.finite(evaluation$probs)) ||
      any(evaluation$probs < 0 | evaluation$probs > 1)) {
    stop("`evaluation_controls$probs` must contain probabilities in [0, 1].",
         call. = FALSE)
  }
  evaluation$probs <- as.numeric(evaluation$probs)

  list(preflight = preflight, evaluation = evaluation)
}


.irtsimrel_sac_capture_rng <- function() {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  list(
    had_seed = had_seed,
    seed = if (had_seed) {
      get(".Random.seed", envir = globalenv(), inherits = FALSE)
    } else {
      NULL
    }
  )
}


.irtsimrel_sac_restore_rng <- function(state) {
  if (isTRUE(state$had_seed)) {
    assign(".Random.seed", state$seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  invisible(NULL)
}


.irtsimrel_sac_stream_seed <- function(master_seed, state, stream_id) {
  modulus <- 2147483646
  if (!is.null(master_seed)) {
    base <- abs(as.double(master_seed)) %% modulus
  } else if (isTRUE(state$had_seed) && length(state$seed) > 0L) {
    base <- 0
    for (value in head(as.double(state$seed), 64L)) {
      base <- (base * 1664525 + (value %% modulus) + 1013904223) %% modulus
    }
  } else {
    base <- 104729
  }
  as.integer((base + as.double(stream_id) * 1000003) %% modulus + 1)
}


.irtsimrel_sac_isolated <- function(stream_seed, fn) {
  state <- .irtsimrel_sac_capture_rng()
  on.exit(.irtsimrel_sac_restore_rng(state), add = TRUE)
  set.seed(stream_seed)
  fn()
}


.irtsimrel_sac_draw_theta <- function(n, latent_shape, latent_params) {
  args <- modifyList(
    list(n = n, shape = latent_shape),
    latent_params
  )
  do.call(sim_latentG, args)$theta
}


.irtsimrel_sac_form_from_items <- function(items, model, n_items) {
  if (!inherits(items, "item_params") || !is.data.frame(items$data) ||
      nrow(items$data) != n_items || !"beta" %in% names(items$data)) {
    stop("Stored item parameters do not define the requested item form.",
         call. = FALSE)
  }
  beta <- as.numeric(items$data$beta)
  lambda_base <- if (identical(model, "rasch")) {
    rep(1, n_items)
  } else if ("lambda_unscaled" %in% names(items$data)) {
    as.numeric(items$data$lambda_unscaled)
  } else {
    as.numeric(items$data$lambda)
  }
  guessing <- if (identical(model, "3pl")) {
    if (!"guessing" %in% names(items$data)) {
      stop("A 3PL item form must contain `guessing`.", call. = FALSE)
    }
    as.numeric(items$data$guessing)
  } else {
    rep(0, n_items)
  }
  bank <- .irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda_base,
    guessing = guessing,
    model = model
  )
  list(
    items = items,
    bank = bank,
    beta = beta,
    lambda_base = lambda_base,
    guessing = guessing
  )
}


.irtsimrel_sac_draw_form <- function(n_items,
                                     model,
                                     item_source,
                                     item_params) {
  args <- modifyList(
    list(
      n_items = n_items,
      model = model,
      source = item_source,
      scale = 1
    ),
    item_params
  )
  items <- do.call(sim_item_params, args)
  .irtsimrel_sac_form_from_items(items, model, n_items)
}


.irtsimrel_sac_form_rho <- function(scale,
                                    theta,
                                    form,
                                    theta_var,
                                    metric_internal) {
  .compute_rho_generic(
    c = scale,
    theta_vec = theta,
    beta_vec = form$beta,
    lambda_base = form$lambda_base,
    theta_var = theta_var,
    metric_internal = metric_internal,
    guessing = form$guessing
  )
}


.irtsimrel_sac_effective_bounds <- function(preflight, margin_fraction) {
  branch_bounds <- as.numeric(preflight$effective_branch_bounds)
  log_bounds <- log(branch_bounds)
  log_root <- log(preflight$selected_root)
  span <- diff(log_bounds)
  requested_margin <- margin_fraction * span
  safe_margin <- min(
    requested_margin,
    0.45 * (log_root - log_bounds[[1L]]),
    0.45 * (log_bounds[[2L]] - log_root)
  )
  if (!is.finite(safe_margin) || safe_margin < 0) safe_margin <- 0
  effective <- exp(c(
    lower = log_bounds[[1L]] + safe_margin,
    upper = log_bounds[[2L]] - safe_margin
  ))
  if (effective[[1L]] >= preflight$selected_root ||
      effective[[2L]] <= preflight$selected_root) {
    .irtsimrel_sac_abort(
      "The selected SAC root is not interior to the guarded branch.",
      "irtsimrel_sac_branch_unavailable",
      preflight
    )
  }
  effective
}


.irtsimrel_sac_max_run <- function(x) {
  if (!any(x)) return(0L)
  runs <- rle(x)
  as.integer(max(runs$lengths[runs$values]))
}


.irtsimrel_sac_has_function <- function(x) {
  if (is.function(x)) return(TRUE)
  if (!is.list(x)) return(FALSE)
  any(vapply(x, .irtsimrel_sac_has_function, logical(1)))
}


.irtsimrel_sac_calibrate_impl <- function(target_rho,
                                          n_items,
                                          model,
                                          latent_shape,
                                          item_source,
                                          latent_params,
                                          item_params,
                                          reliability_metric,
                                          c_init,
                                          M_per_iter,
                                          M_pre,
                                          n_iter,
                                          burn_in,
                                          step_params,
                                          c_bounds,
                                          resample_items,
                                          seed,
                                          verbose,
                                          root_policy,
                                          preflight_controls,
                                          evaluation_controls,
                                          .call,
                                          .item_source_missing,
                                          .item_params_missing) {
  # Public input contract ----------------------------------------------------
  if (!is.numeric(target_rho) || length(target_rho) != 1L ||
      is.na(target_rho) || !is.finite(target_rho) ||
      target_rho <= 0 || target_rho >= 1) {
    stop("`target_rho` must be a scalar in (0, 1).", call. = FALSE)
  }
  n_items <- .irtsimrel_validate_positive_integer_scalar(n_items, "n_items")
  model <- match.arg(model, c("rasch", "2pl", "3pl"))
  latent_shape <- .irtsimrel_match_latent_shape(latent_shape)
  reliability_metric <- match.arg(
    reliability_metric,
    c("msem", "info", "bar", "tilde")
  )
  metric_internal <- if (reliability_metric %in% c("msem", "bar")) {
    "msem"
  } else {
    "info"
  }
  if (identical(metric_internal, "msem") &&
      identical(latent_shape, "heavy_tail")) {
    .irtsimrel_sac_abort(
      paste0(
        "Population MSEM reliability is not finite for the built-in ",
        "Student-t heavy-tail distribution. Use `reliability_metric = ",
        "\"info\"` or an explicitly truncated empirical theta design."
      ),
      "irtsimrel_sac_nonintegrable"
    )
  }

  M_per_iter <- .irtsimrel_validate_positive_integer_scalar(
    M_per_iter,
    "M_per_iter"
  )
  M_pre <- .irtsimrel_validate_positive_integer_scalar(M_pre, "M_pre")
  if (M_pre < 2L) stop("`M_pre` must be at least 2.", call. = FALSE)
  n_iter <- .irtsimrel_validate_positive_integer_scalar(n_iter, "n_iter")
  if (is.null(burn_in)) burn_in <- floor(n_iter / 2)
  if (!is.numeric(burn_in) || length(burn_in) != 1L || is.na(burn_in) ||
      !is.finite(burn_in) || burn_in < 0 || burn_in >= n_iter ||
      burn_in != floor(burn_in)) {
    stop("`burn_in` must be a non-negative integer less than `n_iter`.",
         call. = FALSE)
  }
  burn_in <- as.integer(burn_in)
  if (!is.numeric(c_bounds) || length(c_bounds) != 2L || anyNA(c_bounds) ||
      any(!is.finite(c_bounds)) || any(c_bounds <= 0) ||
      c_bounds[[1L]] >= c_bounds[[2L]]) {
    stop("`c_bounds` must satisfy 0 < lower < upper.", call. = FALSE)
  }
  requested_c_bounds <- as.numeric(c_bounds)
  if (!is.logical(resample_items) || length(resample_items) != 1L ||
      is.na(resample_items)) {
    stop("`resample_items` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!(is.logical(verbose) || is.numeric(verbose)) || length(verbose) != 1L ||
      is.na(verbose) || (is.numeric(verbose) && !is.finite(verbose))) {
    stop("`verbose` must be a logical or non-negative numeric scalar.",
         call. = FALSE)
  }
  verbose_level <- as.integer(verbose)
  if (verbose_level < 0L) {
    stop("`verbose` must be a logical or non-negative numeric scalar.",
         call. = FALSE)
  }
  if (!is.character(root_policy) || length(root_policy) != 1L ||
      is.na(root_policy) || !root_policy %in% c(
        "lowest_increasing", "nearest_increasing", "lowest_any", "highest_any"
      )) {
    stop("`root_policy` is not a supported topology-root policy.",
         call. = FALSE)
  }

  default_step <- list(a = 1, A = 50, gamma = 0.67)
  if (is.null(step_params)) step_params <- list()
  step_params <- .irtsimrel_sac_named_list(step_params, "step_params")
  unknown_step <- setdiff(names(step_params), names(default_step))
  if (length(unknown_step) > 0L) {
    stop("Unknown `step_params` field(s): ",
         paste(unknown_step, collapse = ", "), ".", call. = FALSE)
  }
  step_params <- modifyList(default_step, step_params)
  for (field in names(default_step)) {
    step_params[[field]] <- .irtsimrel_sac_scalar(
      step_params[[field]],
      paste0("step_params$", field)
    )
  }
  if (step_params$a <= 0) stop("`step_params$a` must be positive.", call. = FALSE)
  if (step_params$A < 0) {
    stop("`step_params$A` must be non-negative.", call. = FALSE)
  }
  if (step_params$gamma <= 0.5 || step_params$gamma > 1) {
    stop(
      "`step_params$gamma` must be in (0.5, 1] to satisfy the ",
      "Robbins-Monro step-size conditions.",
      call. = FALSE
    )
  }

  latent_params <- .irtsimrel_normalize_latent_params(latent_params)
  item_params <- .irtsimrel_normalize_item_params(item_params)
  controls <- .irtsimrel_sac_controls(
    preflight_controls,
    evaluation_controls,
    M_pre,
    M_per_iter,
    resample_items
  )
  pre_controls <- controls$preflight
  eval_controls <- controls$evaluation
  item_scope <- if (resample_items) "item_superpopulation" else "fixed_form"

  if (inherits(c_init, "eqc_result")) {
    .irtsimrel_validate_eqc_result_object(c_init, "c_init", "SAC warm start")
    if (!isTRUE(all.equal(c_init$target_rho, target_rho, tolerance = 1e-12))) {
      .irtsimrel_sac_abort(
        "The EQC warm start has a different target reliability.",
        "irtsimrel_sac_warm_start_mismatch",
        c_init
      )
    }
    if (!identical(c_init$model, model) ||
        as.integer(c_init$n_items) != as.integer(n_items)) {
      .irtsimrel_sac_abort(
        "The EQC warm start model or item count does not match SAC.",
        "irtsimrel_sac_warm_start_mismatch",
        c_init
      )
    }
    if (!identical(c_init$metric, metric_internal)) {
      .irtsimrel_sac_abort(
        "The EQC warm start targets a different reliability metric.",
        "irtsimrel_sac_warm_start_mismatch",
        c_init
      )
    }
    if (resample_items) {
      .irtsimrel_sac_abort(
        paste0(
          "An EQC warm start is conditional on its fixed item form; use ",
          "`resample_items = FALSE` for the same item-scope contract."
        ),
        "irtsimrel_sac_warm_start_mismatch",
        c_init
      )
    }
  } else if (!is.null(c_init) &&
             (!is.numeric(c_init) || length(c_init) != 1L || is.na(c_init) ||
              !is.finite(c_init) || c_init <= 0)) {
    stop("`c_init` must be NULL, a positive scalar, or an `eqc_result`.",
         call. = FALSE)
  }

  restore_seed <- .irtsimrel_set_seed(seed)
  if (!is.null(restore_seed)) on.exit(restore_seed(), add = TRUE)

  if (verbose_level >= 1L) {
    message("SAC: estimating latent variance and establishing a branch guard.")
  }

  theta_pre <- .irtsimrel_sac_draw_theta(M_pre, latent_shape, latent_params)
  theta_var <- stats::var(theta_pre)
  if (!is.numeric(theta_var) || length(theta_var) != 1L || is.na(theta_var) ||
      !is.finite(theta_var) || theta_var <= 0) {
    stop("Pre-calculated latent variance must be finite and positive.",
         call. = FALSE)
  }

  warm_start <- list(
    supplied = inherits(c_init, "eqc_result"),
    metric_match = if (inherits(c_init, "eqc_result")) {
      identical(c_init$metric, metric_internal)
    } else {
      NA
    },
    item_scope_match = if (inherits(c_init, "eqc_result")) {
      identical(item_scope, "fixed_form")
    } else {
      NA
    },
    form_reused = FALSE,
    comparable = FALSE
  )

  fixed_form <- NULL
  reuse_eqc_form <- inherits(c_init, "eqc_result") && !resample_items &&
    isTRUE(.item_source_missing) && isTRUE(.item_params_missing)
  if (!resample_items) {
    if (reuse_eqc_form) {
      fixed_form <- .irtsimrel_sac_form_from_items(
        c_init$items_base,
        model,
        n_items
      )
      warm_start$form_reused <- TRUE
      warm_start$comparable <- isTRUE(warm_start$metric_match)
    } else {
      fixed_form <- .irtsimrel_sac_draw_form(
        n_items,
        model,
        item_source,
        item_params
      )
    }
  }
  if (inherits(c_init, "eqc_result") && !reuse_eqc_form) {
    warm_design <- .irtsimrel_normalize_result_item_design(
      c_init,
      "c_init",
      "SAC warm start"
    )
    same_bank <- isTRUE(all.equal(
      fixed_form$beta,
      warm_design$beta,
      tolerance = 1e-10
    )) && isTRUE(all.equal(
      fixed_form$lambda_base,
      c_init$lambda_base,
      tolerance = 1e-10
    )) && isTRUE(all.equal(
      fixed_form$guessing,
      warm_design$guessing,
      tolerance = 1e-10
    ))
    if (!same_bank) {
      .irtsimrel_sac_abort(
        "The EQC warm start and SAC fixed item forms are not the same design.",
        "irtsimrel_sac_warm_start_mismatch",
        list(eqc = warm_design, sac = fixed_form$bank)
      )
    }
    warm_start$comparable <- TRUE
  }

  calibration_rng_state <- .irtsimrel_sac_capture_rng()
  preflight_seed <- .irtsimrel_sac_stream_seed(seed, calibration_rng_state, 1L)
  evaluation_seed <- .irtsimrel_sac_stream_seed(seed, calibration_rng_state, 2L)

  preflight_payload <- .irtsimrel_sac_isolated(preflight_seed, function() {
    theta <- .irtsimrel_sac_draw_theta(
      pre_controls$M,
      latent_shape,
      latent_params
    )
    forms <- if (resample_items) {
      lapply(seq_len(pre_controls$n_forms), function(index) {
        .irtsimrel_sac_draw_form(n_items, model, item_source, item_params)
      })
    } else {
      list(fixed_form)
    }
    list(theta = theta, forms = forms)
  })
  preflight <- .irtsimrel_sac_preflight(
    theta = preflight_payload$theta,
    theta_var = theta_var,
    item_banks = lapply(preflight_payload$forms, `[[`, "bank"),
    target_rho = target_rho,
    metric_internal = metric_internal,
    c_bounds = requested_c_bounds,
    root_policy = root_policy,
    root_controls = pre_controls$root_controls,
    split_check = isTRUE(pre_controls$split_check) && resample_items,
    split_log_tolerance = pre_controls$split_log_tolerance
  )
  effective_c_bounds <- .irtsimrel_sac_effective_bounds(
    preflight,
    pre_controls$branch_margin
  )

  # Initializer --------------------------------------------------------------
  if (is.null(c_init)) {
    if (identical(model, "3pl")) {
      c_init_raw <- preflight$selected_root
      c_init_val <- c_init_raw
      init_method <- "preflight_root"
    } else {
      c_init_raw <- compute_apc_init(target_rho, n_items)
      if (c_init_raw > effective_c_bounds[[1L]] &&
          c_init_raw < effective_c_bounds[[2L]]) {
        c_init_val <- c_init_raw
        init_method <- "apc_warm_start"
      } else {
        c_init_val <- preflight$selected_root
        init_method <- "preflight_root_fallback"
      }
    }
  } else if (inherits(c_init, "eqc_result")) {
    direction <- c_init$misc$selected_root_direction
    if (!is.null(direction) && !identical(direction, "increasing")) {
      .irtsimrel_sac_abort(
        "The EQC warm start is not on an increasing root branch.",
        "irtsimrel_sac_warm_start_mismatch",
        c_init$misc$topology
      )
    }
    c_init_raw <- as.numeric(c_init$c_star)
    if (c_init_raw <= effective_c_bounds[[1L]] ||
        c_init_raw >= effective_c_bounds[[2L]]) {
      .irtsimrel_sac_abort(
        "The EQC warm start is outside the selected SAC branch.",
        "irtsimrel_sac_warm_start_mismatch",
        list(warm_start = c_init_raw, branch = effective_c_bounds)
      )
    }
    c_init_val <- c_init_raw
    init_method <- "eqc_warm_start"
  } else {
    c_init_raw <- as.numeric(c_init)
    c_init_val <- min(
      max(c_init_raw, effective_c_bounds[[1L]]),
      effective_c_bounds[[2L]]
    )
    init_method <- "user_specified"
  }
  c_init_projected <- !isTRUE(all.equal(c_init_raw, c_init_val, tolerance = 0))

  # Robbins--Monro core ------------------------------------------------------
  trajectory <- numeric(n_iter)
  rho_trajectory <- numeric(n_iter)
  rho_update_trajectory <- numeric(n_iter)
  evaluation_trajectory <- numeric(n_iter)
  raw_trajectory <- numeric(n_iter)
  step_size_trajectory <- numeric(n_iter)
  gradient_trajectory <- numeric(n_iter)
  projected <- logical(n_iter)
  projection_side <- rep("none", n_iter)
  projection_source <- rep("none", n_iter)
  c_current <- c_init_val

  progress_checkpoints <- unique(floor(seq(0.1, 1, by = 0.1) * n_iter))
  for (iteration in seq_len(n_iter)) {
    theta_iteration <- .irtsimrel_sac_draw_theta(
      M_per_iter,
      latent_shape,
      latent_params
    )
    form <- if (resample_items) {
      .irtsimrel_sac_draw_form(n_items, model, item_source, item_params)
    } else {
      fixed_form
    }

    rho_update <- .irtsimrel_sac_form_rho(
      c_current,
      theta_iteration,
      form,
      theta_var,
      metric_internal
    )
    step_size <- step_params$a /
      (iteration + step_params$A)^step_params$gamma
    gradient <- rho_update - target_rho
    c_raw <- c_current - step_size * gradient
    c_updated <- min(
      max(c_raw, effective_c_bounds[[1L]]),
      effective_c_bounds[[2L]]
    )
    was_projected <- !isTRUE(all.equal(c_raw, c_updated, tolerance = 0))
    if (was_projected) {
      side <- if (c_raw < effective_c_bounds[[1L]]) "lower" else "upper"
      projection_side[[iteration]] <- side
      requested_same <- isTRUE(all.equal(
        effective_c_bounds[[if (side == "lower") 1L else 2L]],
        requested_c_bounds[[if (side == "lower") 1L else 2L]],
        tolerance = 1e-12
      ))
      projection_source[[iteration]] <- paste0(
        if (requested_same) "requested_" else "branch_",
        side
      )
    }
    rho_aligned <- .irtsimrel_sac_form_rho(
      c_updated,
      theta_iteration,
      form,
      theta_var,
      metric_internal
    )

    evaluation_trajectory[[iteration]] <- c_current
    rho_update_trajectory[[iteration]] <- rho_update
    trajectory[[iteration]] <- c_updated
    rho_trajectory[[iteration]] <- rho_aligned
    raw_trajectory[[iteration]] <- c_raw
    step_size_trajectory[[iteration]] <- step_size
    gradient_trajectory[[iteration]] <- gradient
    projected[[iteration]] <- was_projected
    c_current <- c_updated

    if (verbose_level >= 2L ||
        (verbose_level >= 1L && iteration %in% progress_checkpoints)) {
      message(sprintf(
        "SAC %d/%d: c=%.5f, rho(update)=%.4f, rho(aligned)=%.4f",
        iteration,
        n_iter,
        c_updated,
        rho_update,
        rho_aligned
      ))
    }
  }

  avg_indices <- seq.int(burn_in + 1L, n_iter)
  c_star <- mean(trajectory[avg_indices])
  c_star <- min(
    max(c_star, effective_c_bounds[[1L]]),
    effective_c_bounds[[2L]]
  )

  iteration_trace <- data.frame(
    iteration = seq_len(n_iter),
    c_evaluated = evaluation_trajectory,
    rho_update = rho_update_trajectory,
    gradient = gradient_trajectory,
    step_size = step_size_trajectory,
    c_raw = raw_trajectory,
    c_updated = trajectory,
    rho_updated = rho_trajectory,
    projected = projected,
    projection_side = projection_side,
    projection_source = projection_source,
    stringsAsFactors = FALSE
  )

  # Independent post-calibration evaluation --------------------------------
  evaluation <- .irtsimrel_sac_isolated(evaluation_seed, function() {
    forms <- if (resample_items) {
      lapply(seq_len(eval_controls$n_forms), function(index) {
        .irtsimrel_sac_draw_form(n_items, model, item_source, item_params)
      })
    } else {
      rep(list(fixed_form), eval_controls$n_forms)
    }
    theta_blocks <- lapply(seq_len(eval_controls$n_forms), function(index) {
      .irtsimrel_sac_draw_theta(eval_controls$M, latent_shape, latent_params)
    })
    form_reliability <- vapply(
      seq_len(eval_controls$n_forms),
      function(index) {
        .irtsimrel_sac_form_rho(
          c_star,
          theta_blocks[[index]],
          forms[[index]],
          theta_var,
          metric_internal
        )
      },
      numeric(1)
    )
    list(
      forms = forms,
      theta_blocks = theta_blocks,
      form_reliability = form_reliability
    )
  })

  if (resample_items) {
    achieved_rho <- mean(evaluation$form_reliability)
    theta_final <- evaluation$theta_blocks[[1L]]
    items_final <- evaluation$forms[[1L]]
    aggregation <- "mean_of_form_reliabilities"
  } else {
    theta_final <- unlist(evaluation$theta_blocks, use.names = FALSE)
    items_final <- fixed_form
    achieved_rho <- .irtsimrel_sac_form_rho(
      c_star,
      theta_final,
      fixed_form,
      theta_var,
      metric_internal
    )
    aggregation <- "fixed_form_concatenated_theta"
  }
  representative_achieved_rho <- evaluation$form_reliability[[1L]]
  achieved_sd <- stats::sd(evaluation$form_reliability)
  achieved_se <- achieved_sd / sqrt(eval_controls$n_forms)
  achieved_quantiles <- stats::quantile(
    evaluation$form_reliability,
    probs = eval_controls$probs,
    names = TRUE,
    type = 7
  )
  critical <- stats::qnorm((1 + eval_controls$conf_level) / 2)
  achieved_interval <- c(
    lower = max(0, achieved_rho - critical * achieved_se),
    upper = min(1, achieved_rho + critical * achieved_se)
  )

  items_calib_final <- .irtsimrel_apply_item_scale(
    items_final$items,
    items_final$lambda_base,
    c_star
  )

  # Convergence and branch diagnostics -------------------------------------
  post_burn <- trajectory[avg_indices]
  n_post_burn <- length(post_burn)
  insufficient_post_burn <- n_post_burn < 4L
  split_mean_scale <- max(abs(c_star), .Machine$double.eps)
  split_mean_rel_tolerance <- 0.05
  split_mean_tolerance <- split_mean_rel_tolerance * split_mean_scale
  if (insufficient_post_burn) {
    first_half <- post_burn
    second_half <- post_burn
    split_mean_diff <- NA_real_
    split_mean_rel_diff <- NA_real_
    trajectory_converged <- FALSE
  } else {
    half <- floor(n_post_burn / 2)
    first_half <- post_burn[seq_len(half)]
    second_half <- post_burn[(half + 1L):n_post_burn]
    split_mean_diff <- abs(mean(first_half) - mean(second_half))
    split_mean_rel_diff <- split_mean_diff / split_mean_scale
    trajectory_converged <- split_mean_rel_diff < split_mean_rel_tolerance
  }

  branch_hits <- projected & startsWith(projection_source, "branch_")
  max_consecutive_branch_hits <- .irtsimrel_sac_max_run(branch_hits)
  post_burn_branch_hit_rate <- mean(branch_hits[avg_indices])
  c_slope_delta <- 1e-3
  slope_lower <- max(
    effective_c_bounds[[1L]],
    c_star * exp(-c_slope_delta)
  )
  slope_upper <- min(
    effective_c_bounds[[2L]],
    c_star * exp(c_slope_delta)
  )
  mean_preflight_rho <- function(scale) {
    mean(vapply(
      preflight_payload$forms,
      function(form) {
        .irtsimrel_sac_form_rho(
          scale,
          preflight_payload$theta,
          form,
          theta_var,
          metric_internal
        )
      },
      numeric(1)
    ))
  }
  final_branch_slope <- if (slope_upper > slope_lower) {
    (mean_preflight_rho(slope_upper) - mean_preflight_rho(slope_lower)) /
      (log(slope_upper) - log(slope_lower))
  } else {
    NA_real_
  }
  final_inside_branch <- c_star > effective_c_bounds[[1L]] &&
    c_star < effective_c_bounds[[2L]]
  branch_lost <- max_consecutive_branch_hits >=
    pre_controls$max_consecutive_hits ||
    post_burn_branch_hit_rate > pre_controls$max_post_burn_hit_rate ||
    !isTRUE(final_inside_branch) || !is.finite(final_branch_slope) ||
    final_branch_slope <= 0

  achieved_gap <- achieved_rho - target_rho
  achieved_gap_abs <- abs(achieved_gap)
  achieved_gap_tolerance <- 0.05
  large_achieved_gap <- achieved_gap_abs > achieved_gap_tolerance
  converged <- trajectory_converged && !large_achieved_gap && !branch_lost
  projection_count <- sum(projected)
  projection_rate <- projection_count / n_iter
  hit_lower_bound <- any(trajectory == effective_c_bounds[[1L]]) ||
    (c_init_projected && c_init_val == effective_c_bounds[[1L]])
  hit_upper_bound <- any(trajectory == effective_c_bounds[[2L]]) ||
    (c_init_projected && c_init_val == effective_c_bounds[[2L]])

  calibration_status <- if (branch_lost) {
    "branch_lost"
  } else if (!converged) {
    "not_converged"
  } else if (hit_lower_bound && hit_upper_bound) {
    "hit_both_bounds"
  } else if (hit_lower_bound) {
    "hit_lower_bound"
  } else if (hit_upper_bound) {
    "hit_upper_bound"
  } else {
    "ok"
  }
  status_flags <- character()
  if (!converged) status_flags <- c(status_flags, "not_converged")
  if (insufficient_post_burn) {
    status_flags <- c(status_flags, "insufficient_post_burn")
  }
  if (large_achieved_gap) status_flags <- c(status_flags, "large_achieved_gap")
  if (branch_lost) status_flags <- c(status_flags, "branch_lost")
  if (hit_lower_bound) status_flags <- c(status_flags, "hit_lower_bound")
  if (hit_upper_bound) status_flags <- c(status_flags, "hit_upper_bound")
  if (projection_count > 0L || c_init_projected) {
    status_flags <- c(status_flags, "projection_applied")
  }
  if (length(status_flags) == 0L) status_flags <- "ok"

  convergence <- list(
    mean_first_half = mean(first_half),
    mean_second_half = mean(second_half),
    split_mean_diff = split_mean_diff,
    split_mean_rel_diff = split_mean_rel_diff,
    split_mean_scale = split_mean_scale,
    split_mean_tolerance = split_mean_tolerance,
    split_mean_rel_tolerance = split_mean_rel_tolerance,
    sd_post_burn = if (n_post_burn > 1L) stats::sd(post_burn) else 0,
    range_post_burn = range(post_burn),
    n_post_burn = n_post_burn,
    min_post_burn = 4L,
    insufficient_post_burn = insufficient_post_burn,
    trajectory_converged = trajectory_converged,
    achieved_gap = achieved_gap,
    achieved_gap_abs = achieved_gap_abs,
    achieved_gap_tolerance = achieved_gap_tolerance,
    large_achieved_gap = large_achieved_gap,
    final_gradient = utils::tail(rho_update_trajectory, 1L) - target_rho,
    final_pre_update_gradient = utils::tail(rho_update_trajectory, 1L) -
      target_rho,
    final_iter_gradient = utils::tail(rho_trajectory, 1L) - target_rho,
    gradient_at_c_star = achieved_gap,
    post_calibration_gradient = achieved_gap,
    projection_count = projection_count,
    projection_rate = projection_rate,
    c_init_raw = c_init_raw,
    c_init_projected = c_init_projected,
    n_lower_hits = sum(projection_side == "lower"),
    n_upper_hits = sum(projection_side == "upper"),
    first_projection_iter = if (projection_count > 0L) {
      which(projected)[[1L]]
    } else {
      NA_integer_
    },
    post_burn_projection_rate = mean(projected[avg_indices]),
    post_burn_branch_hit_rate = post_burn_branch_hit_rate,
    max_consecutive_branch_hits = max_consecutive_branch_hits,
    final_branch_slope = final_branch_slope,
    final_inside_branch = final_inside_branch,
    branch_lost = branch_lost,
    hit_lower_bound = hit_lower_bound,
    hit_upper_bound = hit_upper_bound,
    converged = converged,
    status = calibration_status,
    status_flags = unique(status_flags)
  )

  if (!converged) {
    warning(
      "SAC may not have fully converged (status: ", calibration_status,
      "; achieved gap: ", sprintf("%+.4f", achieved_gap), ").",
      call. = FALSE
    )
  }
  if (hit_lower_bound) {
    warning("SAC projection hit the lower effective branch bound.",
            call. = FALSE)
  }
  if (hit_upper_bound) {
    warning("SAC projection hit the upper effective branch bound.",
            call. = FALSE)
  }

  representative_form <- items_final
  achieved_distribution <- list(
    aggregation = aggregation,
    form_reliability = evaluation$form_reliability,
    form_reliabilities = evaluation$form_reliability,
    mean_form_reliability = mean(evaluation$form_reliability),
    sd = achieved_sd,
    se = achieved_se,
    quantiles = achieved_quantiles,
    confidence_level = eval_controls$conf_level,
    confidence_interval = achieved_interval,
    n_forms = eval_controls$n_forms,
    M_per_form = eval_controls$M,
    probs = eval_controls$probs
  )
  branch <- list(
    root_policy = root_policy,
    selected_root = preflight$selected_root,
    selected_root_record = preflight$selected_root_record,
    selected_branch = preflight$selected_branch,
    raw_bounds = preflight$effective_branch_bounds,
    effective_bounds = effective_c_bounds,
    requested_bounds = requested_c_bounds,
    margin_fraction = pre_controls$branch_margin,
    split = preflight$split,
    final_slope = final_branch_slope,
    lost = branch_lost
  )
  calibration_design <- list(
    item_scope = item_scope,
    resample_items = resample_items,
    latent_shape = latent_shape,
    latent_params = latent_params,
    item_source = item_source,
    item_params = item_params,
    theta_measure = "population",
    scale_convention = "global_discrimination_multiplier_D1",
    preflight_n_forms = length(preflight_payload$forms),
    preflight_M = length(preflight_payload$theta)
  )
  evaluation_design <- list(
    item_scope = item_scope,
    aggregation = aggregation,
    independent_from_calibration = TRUE,
    n_forms = eval_controls$n_forms,
    M_per_form = eval_controls$M,
    form_ids = seq_len(eval_controls$n_forms),
    item_banks = lapply(evaluation$forms, `[[`, "bank"),
    theta_blocks = evaluation$theta_blocks,
    representative_form_id = 1L
  )
  estimand_signature <- list(
    metric = metric_internal,
    theta_measure = "population",
    approximation_method = "stochastic_mc",
    item_measure = item_scope
  )
  design_signature <- list(
    schema_version = 1L,
    model = model,
    metric = metric_internal,
    n_items = n_items,
    item_scope = item_scope,
    latent_shape = latent_shape,
    latent_params = latent_params,
    scale_convention = "global_discrimination_multiplier_D1"
  )
  if (identical(item_scope, "fixed_form")) {
    design_signature$beta <- representative_form$beta
    design_signature$lambda_base <- representative_form$lambda_base
    design_signature$guessing <- representative_form$guessing
    design_signature$verifiable <- TRUE
  } else {
    design_signature$item_source <- item_source
    design_signature$item_params <- item_params
    design_signature$verifiable <- !.irtsimrel_sac_has_function(item_params)
  }

  rng_provenance <- list(
    master_seed = seed,
    preflight_stream_seed = preflight_seed,
    evaluation_stream_seed = evaluation_seed,
    preflight_seed = preflight_seed,
    evaluation_seed = evaluation_seed,
    calibration_stream = "outer_stream_after_theta_pre_and_fixed_form",
    stream_overlap = FALSE,
    isolated_streams = TRUE,
    caller_state_restored_when_seeded = !is.null(seed),
    caller_state_restored = !is.null(seed),
    rng_kind = RNGkind(),
    r_version = R.version.string
  )

  if (verbose_level >= 1L) {
    message(sprintf(
      "SAC complete: c*=%.6f, achieved rho=%.4f (SE %.4f), status=%s",
      c_star,
      achieved_rho,
      achieved_se,
      calibration_status
    ))
  }

  result <- list(
    c_star = c_star,
    c_final = c_current,
    target_rho = target_rho,
    achieved_rho = achieved_rho,
    theta_var = theta_var,
    trajectory = trajectory,
    rho_trajectory = rho_trajectory,
    raw_trajectory = raw_trajectory,
    step_size_trajectory = step_size_trajectory,
    gradient_trajectory = gradient_trajectory,
    projected = projected,
    projection_side = projection_side,
    projection_count = projection_count,
    projection_rate = projection_rate,
    M_final = sum(lengths(evaluation$theta_blocks)),
    metric = metric_internal,
    calibration_status = calibration_status,
    model = model,
    n_items = n_items,
    n_iter = n_iter,
    burn_in = burn_in,
    M_per_iter = M_per_iter,
    M_pre = M_pre,
    step_params = step_params,
    c_bounds = effective_c_bounds,
    c_init = c_init_val,
    init_method = init_method,
    item_design = if (resample_items) {
      "post_calibration_draw"
    } else {
      "fixed_iteration_items"
    },
    convergence = convergence,
    beta_vec = representative_form$beta,
    lambda_base = representative_form$lambda_base,
    lambda_scaled = representative_form$lambda_base * c_star,
    items_base = representative_form$items,
    items_calib = items_calib_final,
    theta_quad = theta_final,
    call = .call,
    guessing_vec = representative_form$guessing,
    schema_version = 1L,
    item_scope = item_scope,
    estimand_signature = estimand_signature,
    design_signature = design_signature,
    calibration_design = calibration_design,
    evaluation_design = evaluation_design,
    achieved_distribution = achieved_distribution,
    achieved_se = achieved_se,
    representative_achieved_rho = representative_achieved_rho,
    rho_update_trajectory = rho_update_trajectory,
    evaluation_trajectory = evaluation_trajectory,
    rho_scale_trajectory = evaluation_trajectory,
    iteration_trace = iteration_trace,
    preflight = preflight,
    branch = branch,
    requested_c_bounds = requested_c_bounds,
    projection_source = projection_source,
    rng_provenance = rng_provenance,
    warm_start = warm_start,
    preflight_controls = pre_controls,
    evaluation_controls = eval_controls
  )
  class(result) <- c("sac_result", "list")
  result
}


#' Stochastic Approximation Calibration (Algorithm 2: SAC)
#'
#' `sac_calibrate()` uses branch-guarded Robbins--Monro stochastic
#' approximation to choose a positive global discrimination multiplier that
#' targets a marginal reliability estimand for a unidimensional Rasch, 2PL, or
#' 3PL item generator.
#'
#' The first 18 arguments preserve the v0.2 positional API. The v0.3 branch
#' and evaluation controls are appended after `verbose`.
#'
#' @param target_rho Finite numeric scalar in `(0, 1)`. Target marginal
#'   reliability.
#' @param n_items Positive integer. Number of items per generated form.
#' @param model Character. One of `"rasch"`, `"2pl"`, or `"3pl"`.
#'   Rasch fixes baseline discriminations at one. In the 3PL model,
#'   `item_params$guessing_params` (or `custom_params$guessing` for a custom
#'   source) supplies the lower asymptotes. SAC scales discrimination only and
#'   uses the fixed logistic convention `D = 1`.
#' @param latent_shape Character shape passed to [sim_latentG()].
#' @param item_source Character source passed to [sim_item_params()]. The
#'   default, `"parametric"`, has no external data dependency; `"irw"`
#'   requires the optional Item Response Warehouse package.
#' @param latent_params Named list of additional arguments passed to
#'   [sim_latentG()].
#' @param item_params Named list of additional arguments passed to
#'   [sim_item_params()].
#' @param reliability_metric Character reliability estimand. `"msem"` (the
#'   default) and its synonym `"bar"` target reciprocal-information MSEM
#'   reliability, \eqn{\bar w}. `"info"` and its synonym `"tilde"` target
#'   average-information reliability, \eqn{\tilde\rho}. The stored `metric`
#'   is canonicalized to `"msem"` or `"info"`. Population MSEM is rejected
#'   for the built-in untruncated `"heavy_tail"` distribution because its
#'   reciprocal-information expectation is not finite.
#' @param c_init `NULL`, a positive numeric scalar, or an `eqc_result`. With
#'   `NULL`, 3PL uses the selected preflight root; Rasch/2PL uses APC when APC
#'   lies inside the selected branch and otherwise uses the preflight root. A
#'   numeric initializer is projected into the effective branch. An EQC warm
#'   start is accepted only when target, model, item count, canonical metric,
#'   fixed-form scope, realized item form, and increasing branch agree.
#' @param M_per_iter Positive integer. Ability draws per Robbins--Monro
#'   iteration. Default `500L`.
#' @param M_pre Integer at least two. Draws used once to estimate the population
#'   latent variance held fixed during calibration. Default `10000L`.
#' @param n_iter Positive integer. Robbins--Monro iterations. Default `300L`.
#' @param burn_in Non-negative integer less than `n_iter`. Iterations through
#'   `burn_in` are excluded from Polyak--Ruppert averaging. The default is
#'   `floor(n_iter / 2)`.
#' @param step_params Named list overriding `a = 1`, `A = 50`, or
#'   `gamma = 0.67` in \eqn{a_n=a/(n+A)^\gamma}. `a` must be positive, `A`
#'   non-negative, and `gamma` in `(0.5, 1]`.
#' @param c_bounds Positive increasing length-two numeric vector. The preflight
#'   scanner never expands these requested bounds. Iterates are projected to
#'   the selected increasing branch after applying its safety margin.
#' @param resample_items Logical. `FALSE` calibrates one realized fixed form and
#'   sets `item_scope = "fixed_form"`. `TRUE` (the default) draws a form at each
#'   iteration, targets the item-superpopulation mean, and sets
#'   `item_scope = "item_superpopulation"`.
#' @param seed Optional integer seed. When supplied, the caller's RNG state is
#'   restored on exit; isolated deterministic substreams are used for
#'   preflight and independent post-calibration evaluation.
#' @param verbose Logical or non-negative numeric verbosity level. Level one
#'   reports checkpoints; level two reports every iteration.
#' @param root_policy Character root-selection policy: `"lowest_increasing"`,
#'   `"nearest_increasing"`, `"lowest_any"`, or `"highest_any"`. The default
#'   selects the lowest admissible increasing root. Tangent and plateau roots
#'   are not made admissible by SAC's public controls.
#' @param preflight_controls Named list of branch-preflight overrides:
#'   `n_forms` (default `8L` for item-superpopulation calibration), `M`
#'   (default `min(M_pre, max(500L, M_per_iter))`), `split_check` (default
#'   `resample_items`), `split_log_tolerance` (default `log(1.5)`),
#'   `root_controls` (a named list passed to the topology scanner),
#'   `branch_margin` (default `0.01`), `max_consecutive_hits` (default `5L`),
#'   and `max_post_burn_hit_rate` (default `0.10`). Unknown names are rejected.
#' @param evaluation_controls Named list controlling the independent final
#'   evaluation: `n_forms` (default `20L`), `M` per form (default
#'   `max(1000L, M_per_iter)`), `conf_level` (default `0.95`), and `probs`
#'   (default `c(0.025, 0.5, 0.975)`). Unknown names are rejected. For a fixed
#'   form, theta blocks are concatenated before computing `achieved_rho`; for
#'   item-superpopulation calibration, form reliabilities are averaged.
#' @param x A `sac_result` or `summary.sac_result`, as appropriate.
#' @param digits Number of decimal places used by print methods.
#' @param object A `sac_result`.
#' @param newdata Optional numeric vector of positive discrimination scales for
#'   `predict.sac_result()`. `NULL` returns the stored `achieved_rho`.
#' @param theta_vec Optional finite ability vector for prediction. When
#'   supplied, its sample variance defines the prediction variance basis.
#' @param ... Additional arguments passed to or from methods. Current SAC
#'   methods do not otherwise use them.
#'
#' @return `sac_calibrate()` returns a `sac_result` list. Its historical fields
#'   remain present and v0.3 metadata are appended:
#' \describe{
#'   \item{Calibration}{`c_star` is the Polyak--Ruppert scale; `c_final` is the
#'     last iterate; `target_rho`, `achieved_rho`, `achieved_se`,
#'     `achieved_distribution`, `metric`, and `calibration_status` describe the
#'     independently evaluated result. `M_final` is the total number of final
#'     evaluation theta draws.}
#'   \item{Trajectories}{`trajectory` stores updated scales and
#'     `rho_trajectory` stores reliabilities recomputed at those same scales.
#'     `evaluation_trajectory`, `rho_update_trajectory`,
#'     `raw_trajectory`, `step_size_trajectory`, `gradient_trajectory`,
#'     `projected`, `projection_side`, `projection_source`, and
#'     `iteration_trace` preserve the full update audit trail.}
#'   \item{Algorithm settings}{`theta_var`, `model`, `n_items`, `n_iter`,
#'     `burn_in`, `M_per_iter`, `M_pre`, `step_params`, `c_bounds`,
#'     `requested_c_bounds`, `c_init`, `init_method`, `preflight_controls`, and
#'     `evaluation_controls` record the effective run settings.}
#'   \item{Item design}{`beta_vec`, `lambda_base`, `lambda_scaled`,
#'     `guessing_vec`, `items_base`, `items_calib`, `item_design`, and
#'     `item_scope` describe the representative stored form. Rasch/2PL results
#'     store zero guessing; 3PL results preserve the generated lower
#'     asymptotes. For item-superpopulation runs the representative form is not
#'     itself the estimand.}
#'   \item{Contracts and diagnostics}{`schema_version`, `estimand_signature`,
#'     `design_signature`, `calibration_design`, `evaluation_design`,
#'     `convergence`, `preflight`, `branch`, `rng_provenance`, `warm_start`,
#'     `representative_achieved_rho`, `theta_quad`, and `call` preserve
#'     estimand, topology, evaluation, RNG, and provenance diagnostics.}
#' }
#'
#' `print()` returns its input invisibly; `summary()` returns a
#' `summary.sac_result`; `coef()` returns calibrated item parameters (including
#' `guessing` for 3PL); and `predict()` returns stored or recomputed reliability.
#'
#' @details
#' ## Estimand and branch contract
#'
#' SAC first scans the requested log-scale interval and must resolve an
#' admissible increasing root branch. A classed preflight error is raised when
#' the target is infeasible, a stable branch cannot be resolved, or only an
#' inadmissible boundary/decreasing branch is available. The effective branch
#' bounds are then used as projection guards. Repeated branch-bound contact or
#' a non-positive final branch slope produces `calibration_status =
#' "branch_lost"` rather than silently treating the run as converged.
#'
#' `reliability_metric` and `item_scope` are both part of the estimand. Thus an
#' EQC result (`metric = "info"`, fixed form) is directly comparable to SAC
#' only when SAC also uses `reliability_metric = "info"` and
#' `resample_items = FALSE` on that same form. Passing the compatible EQC object
#' as `c_init` without overriding item-generation arguments reuses its form.
#'
#' ## 3PL and external validation
#'
#' The 3PL probability and Fisher-information kernels use `D = 1`; the global
#' scale multiplies discrimination and never guessing. The package's public
#' [compute_reliability_tam()] helper intentionally supports Rasch and 2PL only
#' because the validated TAM workflow requires both WLE and EAP outputs. Phase
#' 7's 3PL external oracle is an EAP-only validation path; TAM 3PL WLE is not a
#' supported public-package workflow. Do not interpret WLE and EAP as ordered
#' versions of the same estimand.
#'
#' `spc_calibrate()` remains a deprecated forwarding alias.
#'
#' @examples
#' \donttest{
#' # Small fixed-form 3PL calibration with explicit lower asymptotes.
#' form <- list(
#'   custom_params = list(
#'     beta = c(-1.1, -0.55, -0.1, 0.35, 0.8, 1.25),
#'     lambda = c(0.75, 0.9, 1.05, 1.2, 1.35, 1.5),
#'     guessing = c(0.08, 0.12, 0.16, 0.20, 0.24, 0.28)
#'   ),
#'   center_difficulties = FALSE
#' )
#' fit <- sac_calibrate(
#'   target_rho = 0.40,
#'   n_items = 6L,
#'   model = "3pl",
#'   item_source = "custom",
#'   item_params = form,
#'   reliability_metric = "info",
#'   M_per_iter = 128L,
#'   M_pre = 512L,
#'   n_iter = 8L,
#'   burn_in = 4L,
#'   c_bounds = c(0.1, 3),
#'   resample_items = FALSE,
#'   seed = 7601L,
#'   preflight_controls = list(
#'     M = 256L,
#'     split_check = FALSE,
#'     root_controls = list(min_grid = 65L, max_evals = 257L)
#'   ),
#'   evaluation_controls = list(n_forms = 4L, M = 512L)
#' )
#' print(fit)
#' coef(fit)
#' predict(fit, newdata = c(0.5, 1))
#' }
#'
#' @seealso [eqc_calibrate()], [compute_rho_bar()], [compute_rho_tilde()],
#'   [compare_eqc_sac()], [sim_item_params()]
#'
#' @references
#' Robbins, H., & Monro, S. (1951). A stochastic approximation method.
#' *The Annals of Mathematical Statistics, 22*(3), 400--407.
#'
#' Polyak, B. T., & Juditsky, A. B. (1992). Acceleration of stochastic
#' approximation by averaging. *SIAM Journal on Control and Optimization,
#' 30*(4), 838--855.
#'
#' @export
# The public wrapper intentionally keeps the historical 18-argument prefix.
# New v0.3 controls are appended so positional code remains valid.
sac_calibrate <- function(target_rho,
                          n_items,
                          model = c("rasch", "2pl", "3pl"),
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
                          verbose = FALSE,
                          root_policy = "lowest_increasing",
                          preflight_controls = list(),
                          evaluation_controls = list()) {
  item_source_missing <- missing(item_source)
  item_params_missing <- missing(item_params)
  .irtsimrel_sac_calibrate_impl(
    target_rho = target_rho,
    n_items = n_items,
    model = model,
    latent_shape = latent_shape,
    item_source = item_source,
    latent_params = latent_params,
    item_params = item_params,
    reliability_metric = reliability_metric,
    c_init = c_init,
    M_per_iter = M_per_iter,
    M_pre = M_pre,
    n_iter = n_iter,
    burn_in = burn_in,
    step_params = step_params,
    c_bounds = c_bounds,
    resample_items = resample_items,
    seed = seed,
    verbose = verbose,
    root_policy = root_policy,
    preflight_controls = preflight_controls,
    evaluation_controls = evaluation_controls,
    .call = match.call(),
    .item_source_missing = item_source_missing,
    .item_params_missing = item_params_missing
  )
}
