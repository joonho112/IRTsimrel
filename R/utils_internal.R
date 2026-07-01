# Internal seed helper ---------------------------------------------------------

.irtsimrel_set_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) || length(seed) != 1L ||
      is.na(seed) || !is.finite(seed) || seed %% 1 != 0) {
    stop("`seed` must be a single integer value.")
  }
  if (seed < -.Machine$integer.max || seed > .Machine$integer.max) {
    stop(
      "`seed` must be within the signed 32-bit integer range [",
      -.Machine$integer.max, ", ", .Machine$integer.max, "]."
    )
  }

  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }

  set.seed(as.integer(seed))

  function() {
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    invisible(NULL)
  }
}


# Internal argument normalization ---------------------------------------------

.irtsimrel_backtick_collapse <- function(x) {
  paste0("`", x, "`", collapse = ", ")
}

.irtsimrel_check_list_arg <- function(x, arg_name) {
  if (!is.list(x)) {
    stop("`", arg_name, "` must be a list.")
  }
  if (length(x) > 0L && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop("`", arg_name, "` entries must all be named.")
  }
  x
}

.irtsimrel_shape_param_names <- function() {
  c(
    "delta", "df", "k", "mu_mix", "sigma_mix", "weights",
    "w0", "m", "m1", "m2", "w_inner",
    "w_floor", "m_floor", "w_ceil", "m_ceil"
  )
}

.irtsimrel_normalize_latent_params <- function(latent_params,
                                               arg_name = "latent_params",
                                               message_auto_wrap = TRUE) {
  latent_params <- .irtsimrel_check_list_arg(latent_params, arg_name)

  reserved <- c("n", "shape", "seed")
  found_reserved <- intersect(names(latent_params), reserved)
  if (length(found_reserved) > 0L) {
    stop(
      "`", arg_name, "` cannot contain reserved argument(s): ",
      .irtsimrel_backtick_collapse(found_reserved),
      ". Use the corresponding top-level argument instead."
    )
  }

  top_level_shape <- intersect(names(latent_params), .irtsimrel_shape_param_names())
  if (length(top_level_shape) > 0L) {
    if (!is.null(latent_params$shape_params)) {
      stop(
        "`", arg_name, "` contains shape-specific argument(s) ",
        .irtsimrel_backtick_collapse(top_level_shape),
        " and `shape_params`. Put shape-specific arguments inside `",
        arg_name, "$shape_params` only once."
      )
    }

    shape_vals <- latent_params[top_level_shape]
    latent_params[top_level_shape] <- NULL
    latent_params$shape_params <- shape_vals

    if (message_auto_wrap) {
      message(
        "Auto-wrapping shape parameter(s) {",
        paste(top_level_shape, collapse = ", "),
        "} into ", arg_name, "$shape_params."
      )
    }
  }

  latent_params
}

.irtsimrel_normalize_item_params <- function(item_params,
                                             arg_name = "item_params") {
  item_params <- .irtsimrel_check_list_arg(item_params, arg_name)

  if (!is.null(item_params$n_forms)) {
    if (!is.numeric(item_params$n_forms) || length(item_params$n_forms) != 1L ||
        item_params$n_forms != 1L) {
      stop(
        "`", arg_name, "$n_forms` > 1 is only supported in direct calls to ",
        "`sim_item_params()`. Calibration and diagnostic functions are ",
        "single-form only."
      )
    }
    item_params$n_forms <- NULL
  }

  reserved <- c("n_items", "model", "source", "scale", "seed")
  found_reserved <- intersect(names(item_params), reserved)
  if (length(found_reserved) > 0L) {
    stop(
      "`", arg_name, "` cannot contain reserved argument(s): ",
      .irtsimrel_backtick_collapse(found_reserved),
      ". Use top-level arguments for test length, model, item source, ",
      "scale, and seed. Multiple forms (`n_forms > 1`) are only supported ",
      "by direct calls to `sim_item_params()`."
    )
  }

  item_params
}


# Internal item parameter helpers ---------------------------------------------

.irtsimrel_cor_or_na <- function(x, y, method = "pearson") {
  if (length(x) < 2L || length(y) < 2L ||
      stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }

  suppressWarnings(stats::cor(x, y, method = method))
}

.irtsimrel_recompute_item_achieved <- function(data, model) {
  compute_one <- function(df) {
    lambda_unscaled <- if ("lambda_unscaled" %in% names(df)) {
      df$lambda_unscaled
    } else {
      df$lambda
    }

    list(
      beta_mean = mean(df$beta),
      beta_sd = stats::sd(df$beta),
      beta_range = range(df$beta),
      lambda_mean = mean(df$lambda),
      lambda_sd = stats::sd(df$lambda),
      lambda_range = range(df$lambda),
      cor_pearson = if (model == "2pl") {
        .irtsimrel_cor_or_na(df$beta, log(lambda_unscaled))
      } else {
        NA_real_
      },
      cor_spearman = if (model == "2pl") {
        .irtsimrel_cor_or_na(df$beta, log(lambda_unscaled), method = "spearman")
      } else {
        NA_real_
      }
    )
  }

  lambda_unscaled <- if ("lambda_unscaled" %in% names(data)) {
    data$lambda_unscaled
  } else {
    data$lambda
  }

  list(
    by_form = lapply(split(data, data$form_id), compute_one),
    overall = list(
      n_total = nrow(data),
      beta_mean = mean(data$beta),
      beta_sd = stats::sd(data$beta),
      lambda_mean = mean(data$lambda),
      lambda_sd = stats::sd(data$lambda),
      cor_pearson_pooled = if (model == "2pl") {
        .irtsimrel_cor_or_na(data$beta, log(lambda_unscaled))
      } else {
        NA_real_
      },
      cor_spearman_pooled = if (model == "2pl") {
        .irtsimrel_cor_or_na(data$beta, log(lambda_unscaled), method = "spearman")
      } else {
        NA_real_
      }
    )
  )
}

.irtsimrel_apply_item_scale <- function(items, lambda_base, scale) {
  if (!inherits(items, "item_params") || is.null(items$data)) {
    stop("`items` must be an item_params object.")
  }
  if (!is.numeric(lambda_base) || length(lambda_base) != nrow(items$data) ||
      any(!is.finite(lambda_base)) || any(lambda_base <= 0)) {
    stop("`lambda_base` must be a positive numeric vector matching item rows.")
  }
  if (!is.numeric(scale) || length(scale) != 1L ||
      !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a positive finite scalar.")
  }

  out <- items
  out$scale <- scale
  out$data$lambda_unscaled <- lambda_base
  out$data$lambda <- lambda_base * scale
  out$achieved <- .irtsimrel_recompute_item_achieved(out$data, out$model)
  out
}
