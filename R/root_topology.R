# =============================================================================
# root_topology.R
# =============================================================================
# Internal adaptive topology scanner for reliability curves.  Search is carried
# out on log(c), never outside the caller-supplied bounds.  The scanner is
# shared by EQC, feasibility screening, rho-curve diagnostics, and (later) SAC
# preflight so that every public entry point uses the same definition of a
# detected branch and root.
# =============================================================================


.irtsimrel_empty_extrema <- function() {
  data.frame(
    extremum_id = integer(),
    c = numeric(),
    log_c = numeric(),
    rho = numeric(),
    kind = character(),
    component_id = integer(),
    stringsAsFactors = FALSE
  )
}


.irtsimrel_empty_branches <- function() {
  data.frame(
    branch_id = integer(),
    lower = numeric(),
    upper = numeric(),
    log_lower = numeric(),
    log_upper = numeric(),
    rho_lower = numeric(),
    rho_upper = numeric(),
    direction = character(),
    component_id = integer(),
    stringsAsFactors = FALSE
  )
}


.irtsimrel_empty_roots <- function() {
  data.frame(
    root_id = integer(),
    c = numeric(),
    log_c = numeric(),
    rho = numeric(),
    residual = numeric(),
    kind = character(),
    d_rho_d_log_c = numeric(),
    direction = character(),
    branch_id = integer(),
    branch_lower = numeric(),
    branch_upper = numeric(),
    near_extremum = logical(),
    selected = logical(),
    bracket_lower = numeric(),
    bracket_upper = numeric(),
    iterations = integer(),
    estimated_precision = numeric(),
    stringsAsFactors = FALSE
  )
}


.irtsimrel_validate_topology_controls <- function(controls, c_bounds, tol) {
  if (!is.list(controls)) {
    stop("`controls` must be a named list.", call. = FALSE)
  }
  if (length(controls) > 0L &&
      (is.null(names(controls)) || any(names(controls) == ""))) {
    stop("Every element of `controls` must be named.", call. = FALSE)
  }

  defaults <- list(
    points_per_decade = 32,
    min_grid = 65L,
    max_refine = 8L,
    max_evals = 2049L,
    curve_tol = max(10 * tol, 1e-5),
    jump_tol = 0.05,
    root_tol = tol,
    target_tol = max(tol, 1e-8),
    plateau_tol = max(tol, 1e-8),
    slope_tol = max(sqrt(.Machine$double.eps), tol / 10)
  )
  unknown <- setdiff(names(controls), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      "Unknown topology control(s): ", paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  out <- modifyList(defaults, controls)

  positive_scalar <- function(x, name, allow_zero = FALSE) {
    valid <- is.numeric(x) && length(x) == 1L && !is.na(x) &&
      is.finite(x) && if (allow_zero) x >= 0 else x > 0
    if (!valid) {
      qualifier <- if (allow_zero) "nonnegative" else "positive"
      stop("`controls$", name, "` must be a finite ", qualifier,
           " scalar.", call. = FALSE)
    }
    as.numeric(x)
  }
  integer_scalar <- function(x, name, minimum) {
    x <- positive_scalar(x, name, allow_zero = minimum == 0L)
    if (x != floor(x) || x < minimum) {
      stop("`controls$", name, "` must be an integer >= ", minimum,
           ".", call. = FALSE)
    }
    as.integer(x)
  }

  out$points_per_decade <- positive_scalar(
    out$points_per_decade, "points_per_decade"
  )
  out$min_grid <- integer_scalar(out$min_grid, "min_grid", 3L)
  out$max_refine <- integer_scalar(out$max_refine, "max_refine", 0L)
  out$max_evals <- integer_scalar(out$max_evals, "max_evals", 3L)
  out$curve_tol <- positive_scalar(out$curve_tol, "curve_tol")
  out$jump_tol <- positive_scalar(out$jump_tol, "jump_tol")
  out$root_tol <- positive_scalar(out$root_tol, "root_tol")
  out$target_tol <- positive_scalar(out$target_tol, "target_tol")
  out$plateau_tol <- positive_scalar(out$plateau_tol, "plateau_tol")
  out$slope_tol <- positive_scalar(
    out$slope_tol, "slope_tol", allow_zero = TRUE
  )

  decades <- log10(c_bounds[2L] / c_bounds[1L])
  out$requested_initial_grid <- max(
    out$min_grid,
    as.integer(ceiling(out$points_per_decade * decades) + 1L)
  )
  out
}


.irtsimrel_component_ids <- function(finite) {
  ids <- rep.int(NA_integer_, length(finite))
  current <- 0L
  previous_finite <- FALSE
  for (i in seq_along(finite)) {
    if (isTRUE(finite[[i]])) {
      if (!previous_finite) current <- current + 1L
      ids[[i]] <- current
      previous_finite <- TRUE
    } else {
      previous_finite <- FALSE
    }
  }
  ids
}


.irtsimrel_direction <- function(slope, slope_tol) {
  if (!is.finite(slope) || abs(slope) <= slope_tol) {
    "flat"
  } else if (slope > 0) {
    "increasing"
  } else {
    "decreasing"
  }
}


.irtsimrel_dedupe_locations <- function(x, tolerance) {
  if (length(x) <= 1L) return(seq_along(x))
  ordering <- order(x)
  keep <- ordering[[1L]]
  last_x <- x[[ordering[[1L]]]]
  for (index in ordering[-1L]) {
    if (abs(x[[index]] - last_x) > tolerance) {
      keep <- c(keep, index)
      last_x <- x[[index]]
    }
  }
  keep
}


.irtsimrel_unique_grid <- function(x, tolerance = 1e-12) {
  x <- sort(x)
  if (length(x) <= 1L) return(x)
  keep <- c(TRUE, diff(x) > tolerance)
  x[keep]
}


.irtsimrel_abort_topology <- function(message, subclass, topology = NULL) {
  condition <- structure(
    list(message = message, call = NULL, topology = topology),
    class = c(subclass, "irtsimrel_topology_error", "error", "condition")
  )
  stop(condition)
}


#' Internal: Scan Reliability-Curve Topology on log(c)
#'
#' @param fn Scalar reliability function of a positive scale.
#' @param c_bounds Strictly positive lower and upper bounds.
#' @param target Optional target reliability.
#' @param anchors Optional positive scales that must occur in the initial grid.
#' @param tol Root tolerance.
#' @param controls Named scanner-control list.
#'
#' @return An object of class `irtsimrel_root_topology`.
#' @noRd
.irtsimrel_scan_topology <- function(fn,
                                     c_bounds,
                                     target = NULL,
                                     anchors = NULL,
                                     tol = 1e-6,
                                     controls = list()) {
  if (!is.function(fn)) {
    stop("`fn` must be a function.", call. = FALSE)
  }
  if (!is.numeric(c_bounds) || length(c_bounds) != 2L ||
      anyNA(c_bounds) || any(!is.finite(c_bounds)) ||
      any(c_bounds <= 0) || c_bounds[[1L]] >= c_bounds[[2L]]) {
    stop(
      "`c_bounds` must contain finite values with 0 < lower < upper.",
      call. = FALSE
    )
  }
  c_bounds <- as.numeric(c_bounds)
  if (!is.null(target) &&
      (!is.numeric(target) || length(target) != 1L || is.na(target) ||
       !is.finite(target))) {
    stop("`target` must be NULL or a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a finite positive scalar.", call. = FALSE)
  }
  if (!is.null(anchors)) {
    if (!is.numeric(anchors) || anyNA(anchors) ||
        any(!is.finite(anchors)) || any(anchors <= 0)) {
      stop("`anchors` must be NULL or finite positive scales.", call. = FALSE)
    }
    anchors <- as.numeric(anchors)
  }

  controls <- .irtsimrel_validate_topology_controls(controls, c_bounds, tol)
  x_bounds <- log(c_bounds)
  cache <- new.env(hash = TRUE, parent = emptyenv())
  eval_count <- 0L
  budget_hit <- FALSE

  evaluate_one <- function(x) {
    key <- sprintf("%.17g", x)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    if (eval_count >= controls$max_evals) {
      budget_hit <<- TRUE
      return(NA_real_)
    }
    eval_count <<- eval_count + 1L
    value <- tryCatch(fn(exp(x)), error = function(e) NA_real_)
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value)) {
      value <- NA_real_
    } else {
      value <- as.numeric(value)
    }
    assign(key, value, envir = cache)
    value
  }
  evaluate_many <- function(x) {
    vapply(x, evaluate_one, numeric(1))
  }

  required_x <- x_bounds
  if (c_bounds[[1L]] <= 1 && c_bounds[[2L]] >= 1) {
    required_x <- c(required_x, 0)
  }
  if (length(anchors) > 0L) {
    inside <- anchors >= c_bounds[[1L]] & anchors <= c_bounds[[2L]]
    required_x <- c(required_x, log(anchors[inside]))
  }
  required_x <- sort(unique(required_x))

  requested_n <- controls$requested_initial_grid
  available_regular <- max(3L, controls$max_evals - length(required_x) + 2L)
  initial_n <- min(requested_n, available_regular)
  x_grid <- .irtsimrel_unique_grid(c(
    seq(x_bounds[[1L]], x_bounds[[2L]], length.out = initial_n),
    required_x
  ))
  if (length(x_grid) > controls$max_evals) {
    # Keep all required anchors, filling the remaining budget uniformly.
    regular <- seq(x_bounds[[1L]], x_bounds[[2L]],
                   length.out = controls$max_evals)
    distance_to_required <- vapply(
      regular,
      function(x) min(abs(x - required_x)),
      numeric(1)
    )
    remaining <- max(0L, controls$max_evals - length(required_x))
    fill <- if (remaining > 0L) {
      regular[order(distance_to_required, decreasing = TRUE)[
        seq_len(remaining)
      ]]
    } else {
      numeric()
    }
    x_grid <- .irtsimrel_unique_grid(c(required_x, fill))
  }
  y_grid <- evaluate_many(x_grid)
  if (initial_n < requested_n) budget_hit <- TRUE

  stop_reason <- if (controls$max_refine == 0L) {
    "refinement_disabled"
  } else {
    "max_refine_reached"
  }

  # Refinement concentrates on target brackets, slope reversals, and large
  # local changes. The deterministic initial grid already supplies at least
  # 32 points per log-decade; local refinement retains most of the evaluation
  # budget for polishing the actually observed topology.
  if (controls$max_refine > 0L && !budget_hit) {
    for (iteration in seq_len(controls$max_refine)) {
      n <- length(x_grid)
      if (n < 2L) break
      midpoints <- (x_grid[-n] + x_grid[-1L]) / 2

      score <- rep.int(0, length(midpoints))
      pair_finite <- is.finite(y_grid[-n]) & is.finite(y_grid[-1L])
      jump <- abs(diff(y_grid))
      score[pair_finite & jump > controls$jump_tol] <- 20

      if (!is.null(target)) {
        residual <- y_grid - target
        pair_slope <- diff(y_grid) / diff(x_grid)
        flat_near_pair <- pair_finite &
          abs(residual[-n]) <= controls$plateau_tol &
          abs(residual[-1L]) <= controls$plateau_tol &
          abs(pair_slope) <= controls$slope_tol
        crossing <- pair_finite &
          residual[-n] * residual[-1L] <= 0 & !flat_near_pair
        score[crossing] <- pmax(score[crossing], 100)
      }

      slopes <- diff(y_grid) / diff(x_grid)
      if (length(slopes) >= 2L) {
        left_slopes <- slopes[-length(slopes)]
        right_slopes <- slopes[-1L]
        turns <- which(
          is.finite(left_slopes) & is.finite(right_slopes) &
            ((left_slopes > controls$slope_tol &
                right_slopes < -controls$slope_tol) |
               (left_slopes < -controls$slope_tol &
                  right_slopes > controls$slope_tol))
        )
        for (turn in turns) {
          score[turn] <- max(score[turn], 80)
          score[turn + 1L] <- max(score[turn + 1L], 80)
        }

        # A small, capped curvature probe in the first two rounds can expose a
        # narrow feature before it becomes an observed slope reversal.  The
        # cap prevents a smooth curved objective from doubling the full grid.
        if (iteration <= 2L && n >= 3L) {
          interior <- 2:(n - 1L)
          finite_triples <- y_grid[interior - 1L]
          finite_triples <- is.finite(finite_triples) &
            is.finite(y_grid[interior]) & is.finite(y_grid[interior + 1L])
          fraction <- (x_grid[interior] - x_grid[interior - 1L]) /
            (x_grid[interior + 1L] - x_grid[interior - 1L])
          linear_prediction <- y_grid[interior - 1L] + fraction *
            (y_grid[interior + 1L] - y_grid[interior - 1L])
          curvature_error <- abs(y_grid[interior] - linear_prediction)
          curved <- which(finite_triples &
                            curvature_error > controls$curve_tol)
          if (length(curved) > 0L) {
            curved <- curved[order(
              curvature_error[curved], decreasing = TRUE
            )]
            curved <- head(curved, 4L)
            for (position in interior[curved]) {
              score[position - 1L] <- max(score[position - 1L], 40)
              score[position] <- max(score[position], 40)
            }
          }
        }
      }

      chosen <- which(score > 0)
      if (length(chosen) == 0L) {
        stop_reason <- "local_resolution_stable"
        break
      }
      candidate_x <- midpoints[chosen]
      candidate_score <- score[chosen]
      remaining <- controls$max_evals - eval_count
      if (remaining <= 0L) {
        budget_hit <- TRUE
        stop_reason <- "evaluation_budget_exhausted"
        break
      }
      if (length(candidate_x) > remaining) {
        ordering <- order(candidate_score, decreasing = TRUE, candidate_x)
        candidate_x <- candidate_x[ordering[seq_len(remaining)]]
        budget_hit <- TRUE
        stop_reason <- "evaluation_budget_exhausted"
      }
      candidate_y <- evaluate_many(candidate_x)
      x_grid <- c(x_grid, candidate_x)
      y_grid <- c(y_grid, candidate_y)
      ordering <- order(x_grid)
      x_grid <- x_grid[ordering]
      y_grid <- y_grid[ordering]
      distinct <- c(TRUE, diff(x_grid) > 1e-12)
      x_grid <- x_grid[distinct]
      y_grid <- y_grid[distinct]
      if (budget_hit) break
    }
  }

  scan <- data.frame(
    c = exp(x_grid),
    log_c = x_grid,
    rho = y_grid,
    finite = is.finite(y_grid),
    stringsAsFactors = FALSE
  )
  scan$component_id <- .irtsimrel_component_ids(scan$finite)

  # Detect and polish every sampled slope reversal independently.  Using all
  # local extrema rather than one global optimize() is essential for the
  # empirically observed four-root 2PL and 3PL curves.
  extrema_rows <- list()
  if (nrow(scan) >= 3L) {
    slopes <- diff(scan$rho) / diff(scan$log_c)
    for (i in 2:(nrow(scan) - 1L)) {
      left_slope <- slopes[[i - 1L]]
      right_slope <- slopes[[i]]
      if (!all(scan$finite[c(i - 1L, i, i + 1L)]) ||
          !is.finite(left_slope) || !is.finite(right_slope)) {
        next
      }
      kind <- if (left_slope > controls$slope_tol &&
                  right_slope < -controls$slope_tol) {
        "maximum"
      } else if (left_slope < -controls$slope_tol &&
                 right_slope > controls$slope_tol) {
        "minimum"
      } else {
        NULL
      }
      if (is.null(kind)) next

      interval <- scan$log_c[c(i - 1L, i + 1L)]
      objective <- function(x) {
        value <- evaluate_one(x)
        if (!is.finite(value)) return(.Machine$double.xmax / 100)
        if (identical(kind, "maximum")) -value else value
      }
      polished <- if (eval_count < controls$max_evals) {
        tryCatch(
          optimize(objective, interval = interval, tol = controls$root_tol),
          error = function(e) NULL
        )
      } else {
        NULL
      }
      if (is.null(polished)) {
        x_value <- scan$log_c[[i]]
        rho_value <- scan$rho[[i]]
      } else {
        x_value <- polished$minimum
        rho_value <- evaluate_one(x_value)
      }
      if (!is.finite(rho_value)) next
      extrema_rows[[length(extrema_rows) + 1L]] <- data.frame(
        extremum_id = NA_integer_,
        c = exp(x_value),
        log_c = x_value,
        rho = rho_value,
        kind = kind,
        component_id = scan$component_id[[i]],
        stringsAsFactors = FALSE
      )
    }
  }
  extrema <- if (length(extrema_rows) == 0L) {
    .irtsimrel_empty_extrema()
  } else {
    do.call(rbind, extrema_rows)
  }
  if (nrow(extrema) > 0L) {
    extrema <- extrema[order(extrema$log_c), , drop = FALSE]
    keep <- .irtsimrel_dedupe_locations(
      extrema$log_c,
      max(controls$root_tol * 10, 1e-7)
    )
    extrema <- extrema[keep, , drop = FALSE]
    extrema$extremum_id <- seq_len(nrow(extrema))
    rownames(extrema) <- NULL
  }

  # Branches are the finite monotone segments separated by polished extrema.
  branch_rows <- list()
  component_values <- sort(unique(scan$component_id[scan$finite]))
  for (component in component_values) {
    component_scan <- scan[scan$component_id == component & scan$finite,
                           , drop = FALSE]
    if (nrow(component_scan) == 0L) next
    component_extrema <- extrema[extrema$component_id == component,
                                 , drop = FALSE]
    breaks_x <- sort(unique(c(
      min(component_scan$log_c),
      component_extrema$log_c,
      max(component_scan$log_c)
    )))
    if (length(breaks_x) == 1L) {
      breaks_x <- rep(breaks_x, 2L)
    }
    for (i in seq_len(length(breaks_x) - 1L)) {
      lower_x <- breaks_x[[i]]
      upper_x <- breaks_x[[i + 1L]]
      lower_rho <- evaluate_one(lower_x)
      upper_rho <- evaluate_one(upper_x)
      slope <- if (upper_x > lower_x) {
        (upper_rho - lower_rho) / (upper_x - lower_x)
      } else {
        0
      }
      branch_rows[[length(branch_rows) + 1L]] <- data.frame(
        branch_id = NA_integer_,
        lower = exp(lower_x),
        upper = exp(upper_x),
        log_lower = lower_x,
        log_upper = upper_x,
        rho_lower = lower_rho,
        rho_upper = upper_rho,
        direction = .irtsimrel_direction(slope, controls$slope_tol),
        component_id = component,
        stringsAsFactors = FALSE
      )
    }
  }
  branches <- if (length(branch_rows) == 0L) {
    .irtsimrel_empty_branches()
  } else {
    do.call(rbind, branch_rows)
  }
  if (nrow(branches) > 0L) {
    branches$branch_id <- seq_len(nrow(branches))
    rownames(branches) <- NULL
  }

  derivative_at <- function(x) {
    span <- diff(x_bounds)
    h <- min(max(sqrt(controls$root_tol), 1e-5), span / 100)
    left <- max(x_bounds[[1L]], x - h)
    right <- min(x_bounds[[2L]], x + h)
    if (right <= left) return(NA_real_)
    values <- c(evaluate_one(left), evaluate_one(right))
    if (all(is.finite(values))) {
      return((values[[2L]] - values[[1L]]) / (right - left))
    }
    index <- which.min(abs(scan$log_c - x))
    candidates <- integer()
    if (index > 1L) candidates <- c(candidates, index - 1L)
    if (index < nrow(scan)) candidates <- c(candidates, index)
    for (candidate in candidates) {
      if (all(scan$finite[c(candidate, candidate + 1L)])) {
        return(
          diff(scan$rho[c(candidate, candidate + 1L)]) /
            diff(scan$log_c[c(candidate, candidate + 1L)])
        )
      }
    }
    NA_real_
  }

  root_rows <- list()
  add_root <- function(x,
                       kind,
                       bracket = c(NA_real_, NA_real_),
                       iterations = NA_integer_,
                       precision = NA_real_,
                       forced_direction = NULL,
                       plateau_bounds = NULL) {
    rho_value <- evaluate_one(x)
    if (!is.finite(rho_value)) return(invisible(NULL))
    slope <- derivative_at(x)
    direction <- if (is.null(forced_direction)) {
      .irtsimrel_direction(slope, controls$slope_tol)
    } else {
      forced_direction
    }
    branch_id <- NA_integer_
    branch_lower <- NA_real_
    branch_upper <- NA_real_
    if (!is.null(plateau_bounds)) {
      branch_lower <- exp(plateau_bounds[[1L]])
      branch_upper <- exp(plateau_bounds[[2L]])
    } else if (nrow(branches) > 0L) {
      contained <- which(
        x >= branches$log_lower - controls$root_tol * 10 &
          x <= branches$log_upper + controls$root_tol * 10
      )
      if (length(contained) > 0L) {
        if (direction == "increasing") {
          matching <- contained[branches$direction[contained] == "increasing"]
        } else if (direction == "decreasing") {
          matching <- contained[branches$direction[contained] == "decreasing"]
        } else {
          matching <- contained
        }
        chosen <- if (length(matching) > 0L) matching[[1L]] else contained[[1L]]
        branch_id <- branches$branch_id[[chosen]]
        branch_lower <- branches$lower[[chosen]]
        branch_upper <- branches$upper[[chosen]]
      }
    }
    near_extremum <- nrow(extrema) > 0L &&
      any(abs(extrema$log_c - x) <= max(controls$root_tol * 10, 1e-6))
    root_rows[[length(root_rows) + 1L]] <<- data.frame(
      root_id = NA_integer_,
      c = exp(x),
      log_c = x,
      rho = rho_value,
      residual = rho_value - target,
      kind = kind,
      d_rho_d_log_c = slope,
      direction = direction,
      branch_id = branch_id,
      branch_lower = branch_lower,
      branch_upper = branch_upper,
      near_extremum = near_extremum,
      selected = FALSE,
      bracket_lower = exp(bracket[[1L]]),
      bracket_upper = exp(bracket[[2L]]),
      iterations = as.integer(iterations),
      estimated_precision = precision,
      stringsAsFactors = FALSE
    )
    invisible(NULL)
  }

  if (!is.null(target)) {
    residual <- scan$rho - target
    near_target <- scan$finite & abs(residual) <= controls$plateau_tol

    # Collapse a genuinely flat, non-trivial near-target run to one plateau.
    # Dense refinement around an ordinary crossing also creates several
    # near-target points, so those groups are bracketed and polished instead.
    near_indices <- which(near_target)
    handled_indices <- integer()
    if (length(near_indices) > 0L) {
      groups <- split(
        near_indices,
        cumsum(c(1L, diff(near_indices) != 1L))
      )
      for (group in groups) {
        handled_indices <- c(handled_indices, group)
        group_x <- scan$log_c[range(group)]
        group_span <- diff(group_x)
        group_slope <- if (group_span > 0) {
          diff(scan$rho[range(group)]) / group_span
        } else {
          Inf
        }
        plateau <- length(group) >= 2L &&
          group_span >= max(100 * controls$root_tol, 1e-5) &&
          is.finite(group_slope) && abs(group_slope) <= controls$slope_tol
        if (plateau) {
          add_root(
            mean(group_x),
            kind = "plateau",
            forced_direction = "flat",
            plateau_bounds = group_x
          )
          next
        }

        before <- min(group) - 1L
        after <- max(group) + 1L
        has_before <- before >= 1L && scan$finite[[before]]
        has_after <- after <= nrow(scan) && scan$finite[[after]]
        if (eval_count < controls$max_evals && has_before && has_after &&
            residual[[before]] * residual[[after]] < 0) {
          interval <- scan$log_c[c(before, after)]
          root_result <- tryCatch(
            uniroot(
              function(x) evaluate_one(x) - target,
              interval = interval,
              tol = controls$root_tol
            ),
            error = function(e) NULL
          )
          if (!is.null(root_result)) {
            add_root(
              root_result$root,
              kind = "crossing",
              bracket = interval,
              iterations = root_result$iter,
              precision = root_result$estim.prec
            )
            next
          }
        }

        closest <- group[[which.min(abs(residual[group]))]]
        boundary <- min(group) == 1L || max(group) == nrow(scan)
        tangent <- has_before && has_after &&
          sign(residual[[before]]) == sign(residual[[after]])
        add_root(
          scan$log_c[[closest]],
          kind = if (boundary) "boundary" else if (tangent) "tangent" else
            "crossing",
          forced_direction = if (tangent && !boundary) "flat" else NULL
        )
      }
    }

    if (nrow(scan) >= 2L) {
      for (i in seq_len(nrow(scan) - 1L)) {
        if (eval_count >= controls$max_evals) break
        if (!all(scan$finite[c(i, i + 1L)]) ||
            near_target[[i]] || near_target[[i + 1L]] ||
            residual[[i]] * residual[[i + 1L]] >= 0) {
          next
        }
        interval <- scan$log_c[c(i, i + 1L)]
        root_result <- tryCatch(
          uniroot(
            function(x) evaluate_one(x) - target,
            interval = interval,
            tol = controls$root_tol
          ),
          error = function(e) NULL
        )
        if (is.null(root_result)) next
        add_root(
          root_result$root,
          kind = if (abs(root_result$root - x_bounds[[1L]]) <=
                     controls$root_tol ||
                     abs(root_result$root - x_bounds[[2L]]) <=
                     controls$root_tol) "boundary" else "crossing",
          bracket = interval,
          iterations = root_result$iter,
          precision = root_result$estim.prec
        )
      }
    }

    # A polished extremum can touch the target between sampled points without
    # changing sign.  Record it explicitly as a tangent candidate.
    if (nrow(extrema) > 0L) {
      tangent_extrema <- which(
        abs(extrema$rho - target) <= controls$target_tol
      )
      for (i in tangent_extrema) {
        add_root(
          extrema$log_c[[i]],
          kind = "tangent",
          forced_direction = "flat"
        )
      }
    }
  }

  roots <- if (length(root_rows) == 0L) {
    .irtsimrel_empty_roots()
  } else {
    do.call(rbind, root_rows)
  }
  if (nrow(roots) > 0L) {
    # Prefer crossings over grid hits, and tangencies over duplicate near-hits.
    kind_priority <- match(
      roots$kind,
      c("boundary", "crossing", "tangent", "plateau")
    )
    ordering <- order(roots$log_c, kind_priority)
    roots <- roots[ordering, , drop = FALSE]
    keep <- .irtsimrel_dedupe_locations(
      roots$log_c,
      max(controls$root_tol * 10, 1e-7)
    )
    roots <- roots[keep, , drop = FALSE]
    roots <- roots[order(roots$log_c), , drop = FALSE]
    roots$root_id <- seq_len(nrow(roots))
    rownames(roots) <- NULL
  }

  finite_values <- scan$rho[scan$finite]
  if (nrow(extrema) > 0L) finite_values <- c(finite_values, extrema$rho)
  rho_range <- if (length(finite_values) > 0L) {
    c(lower = min(finite_values), upper = max(finite_values))
  } else {
    c(lower = NA_real_, upper = NA_real_)
  }
  rho_bounds <- c(lower = scan$rho[[1L]], upper = scan$rho[[nrow(scan)]])

  topology_status <- if (any(!scan$finite)) {
    "nonfinite"
  } else if (budget_hit) {
    "resolution_limited"
  } else {
    "resolved"
  }
  if (budget_hit) stop_reason <- "evaluation_budget_exhausted"

  target_status <- NULL
  best_achievable <- NULL
  if (!is.null(target)) {
    if (!identical(topology_status, "resolved") ||
        any(!is.finite(rho_range))) {
      target_status <- "uncertain"
    } else if (target < rho_range[["lower"]] - controls$target_tol) {
      target_status <- "infeasible_below_range"
    } else if (target > rho_range[["upper"]] + controls$target_tol) {
      target_status <- "infeasible_above_range"
    } else {
      target_status <- "feasible"
    }

    candidates <- data.frame(
      c = scan$c[scan$finite],
      rho = scan$rho[scan$finite],
      location = rep.int("scan", sum(scan$finite)),
      stringsAsFactors = FALSE
    )
    if (nrow(extrema) > 0L) {
      candidates <- rbind(
        candidates,
        data.frame(
          c = extrema$c,
          rho = extrema$rho,
          location = "extremum",
          stringsAsFactors = FALSE
        )
      )
    }
    if (nrow(roots) > 0L) {
      candidates <- rbind(
        candidates,
        data.frame(
          c = roots$c,
          rho = roots$rho,
          location = roots$kind,
          stringsAsFactors = FALSE
        )
      )
    }
    if (nrow(candidates) == 0L) {
      best_achievable <- list(
        c = NA_real_,
        rho = NA_real_,
        residual = NA_real_,
        absolute_error = NA_real_,
        location = "unavailable"
      )
    } else {
      candidates$residual <- candidates$rho - target
      candidates$absolute_error <- abs(candidates$residual)
      ordering <- order(candidates$absolute_error, candidates$c)
      best <- candidates[ordering[[1L]], , drop = FALSE]
      best_achievable <- list(
        c = best$c[[1L]],
        rho = best$rho[[1L]],
        residual = best$residual[[1L]],
        absolute_error = best$absolute_error[[1L]],
        location = best$location[[1L]]
      )
    }
  }

  admissible <- nrow(roots) > 0L &
    roots$direction == "increasing" &
    roots$kind %in% c("crossing", "boundary")
  result <- list(
    c_bounds = c_bounds,
    target = target,
    tol = tol,
    controls = controls,
    scan = scan,
    extrema = extrema,
    branches = branches,
    roots = roots,
    rho_bounds = rho_bounds,
    rho_range = rho_range,
    best_achievable = best_achievable,
    topology_status = topology_status,
    target_status = target_status,
    root_count = nrow(roots),
    admissible_root_count = sum(admissible),
    eval_count = eval_count,
    stop_reason = stop_reason
  )
  class(result) <- c("irtsimrel_root_topology", "list")
  result
}


#' Internal: Select One Root Under an Explicit Branch Policy
#'
#' @noRd
.irtsimrel_select_root <- function(topology,
                                   root_policy = "lowest_increasing",
                                   reference = NULL,
                                   allow_tangent = FALSE,
                                   allow_plateau = FALSE) {
  if (!inherits(topology, "irtsimrel_root_topology")) {
    stop("`topology` must be an irtsimrel_root_topology object.",
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
  if (!is.logical(allow_tangent) || length(allow_tangent) != 1L ||
      is.na(allow_tangent)) {
    stop("`allow_tangent` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(allow_plateau) || length(allow_plateau) != 1L ||
      is.na(allow_plateau)) {
    stop("`allow_plateau` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!identical(topology$topology_status, "resolved")) {
    .irtsimrel_abort_topology(
      paste0(
        "Cannot select a root from topology status '",
        topology$topology_status, "'."
      ),
      subclass = "irtsimrel_topology_unresolved",
      topology = topology
    )
  }

  roots <- topology$roots
  selectable_kind <- roots$kind %in% c("crossing", "boundary")
  if (allow_tangent) selectable_kind <- selectable_kind | roots$kind == "tangent"
  if (allow_plateau) selectable_kind <- selectable_kind | roots$kind == "plateau"
  increasing <- roots$direction == "increasing"

  candidates <- switch(
    root_policy,
    lowest_increasing = which(selectable_kind & increasing),
    nearest_increasing = which(selectable_kind & increasing),
    lowest_any = which(selectable_kind),
    highest_any = which(selectable_kind)
  )
  if (length(candidates) == 0L) {
    .irtsimrel_abort_topology(
      paste0(
        "No admissible root is available under root_policy = '",
        root_policy,
        "'. Detected roots must satisfy the requested direction and kind."
      ),
      subclass = "irtsimrel_selection_unavailable",
      topology = topology
    )
  }

  if (identical(root_policy, "nearest_increasing")) {
    if (!is.numeric(reference) || length(reference) != 1L ||
        is.na(reference) || !is.finite(reference) || reference <= 0) {
      stop(
        "`reference` must be a finite positive scalar for nearest_increasing.",
        call. = FALSE
      )
    }
    chosen <- candidates[[which.min(abs(log(roots$c[candidates] / reference)))]]
  } else if (identical(root_policy, "highest_any")) {
    chosen <- candidates[[which.max(roots$c[candidates])]]
  } else {
    chosen <- candidates[[which.min(roots$c[candidates])]]
  }

  topology$roots$selected <- FALSE
  topology$roots$selected[[chosen]] <- TRUE
  root <- topology$roots[chosen, , drop = FALSE]
  rationale <- paste0(
    "Selected root ", root$root_id[[1L]], " under ", root_policy,
    " (", root$direction[[1L]], ", ", root$kind[[1L]], ")."
  )
  topology$selection <- list(
    root_policy = root_policy,
    reference = reference,
    allow_tangent = allow_tangent,
    allow_plateau = allow_plateau,
    root_id = root$root_id[[1L]],
    rationale = rationale
  )
  uniroot_result <- if (identical(root$kind[[1L]], "crossing")) {
    list(
      root = root$c[[1L]],
      f.root = root$residual[[1L]],
      iter = root$iterations[[1L]],
      init.it = NA_integer_,
      estim.prec = root$estimated_precision[[1L]]
    )
  } else {
    NULL
  }

  list(
    topology = topology,
    c = root$c[[1L]],
    root_id = root$root_id[[1L]],
    root = root,
    uniroot_result = uniroot_result,
    rationale = rationale
  )
}
