# =============================================================================
# model_engine.R
# =============================================================================
# Internal, numerically stable unidimensional IRT model engine.
#
# The package-facing parameterization is
#
#   P_i(theta; scale) = g_i + (1 - g_i) * logistic{
#     scale * lambda_i * (theta - beta_i)
#   },
#
# with D = 1.  The normalized bank representation also carries the equivalent
# slope/intercept form so that a future multidimensional engine can consume the
# same adapter without changing the public beta/lambda/guessing vocabulary.
# =============================================================================


#' Internal: Normalize a Unidimensional IRT Item Bank
#'
#' @param beta Numeric item difficulties.
#' @param lambda_base Numeric positive baseline discriminations.
#' @param guessing Optional lower-asymptote parameter(s). A scalar is recycled.
#' @param model Optional model label: `"rasch"`, `"2pl"`, or `"3pl"`.
#'
#' @return A normalized internal item-bank list.
#' @noRd
.irtsimrel_item_bank <- function(beta,
                                 lambda_base,
                                 guessing = NULL,
                                 model = NULL) {
  if (!is.numeric(beta) || length(beta) == 0L ||
      anyNA(beta) || any(!is.finite(beta))) {
    stop("beta must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (!is.numeric(lambda_base) || length(lambda_base) == 0L ||
      anyNA(lambda_base) || any(!is.finite(lambda_base)) ||
      any(lambda_base <= 0)) {
    stop(
      "lambda_base must be a non-empty finite numeric vector with all values > 0.",
      call. = FALSE
    )
  }
  if (length(beta) != length(lambda_base)) {
    stop("beta and lambda_base must have the same length.", call. = FALSE)
  }

  beta <- as.numeric(beta)
  lambda_base <- as.numeric(lambda_base)
  n_items <- length(beta)

  if (!is.null(guessing)) {
    if (!is.numeric(guessing) || length(guessing) == 0L ||
        anyNA(guessing) || any(!is.finite(guessing))) {
      stop("guessing must contain finite numeric values.", call. = FALSE)
    }
    if (length(guessing) == 1L) {
      guessing <- rep(as.numeric(guessing), n_items)
    } else if (length(guessing) != n_items) {
      stop(
        "guessing must have length 1 or the same length as beta.",
        call. = FALSE
      )
    } else {
      guessing <- as.numeric(guessing)
    }
    if (any(guessing < 0 | guessing >= 1)) {
      stop("guessing values must satisfy 0 <= guessing < 1.", call. = FALSE)
    }
  }

  if (is.null(model)) {
    model <- if (!is.null(guessing) && any(guessing > 0)) "3pl" else "2pl"
  } else {
    if (!is.character(model) || length(model) != 1L || is.na(model)) {
      stop("model must be one of 'rasch', '2pl', or '3pl'.", call. = FALSE)
    }
    model <- tolower(model)
    if (!model %in% c("rasch", "2pl", "3pl")) {
      stop("model must be one of 'rasch', '2pl', or '3pl'.", call. = FALSE)
    }
  }

  if (identical(model, "3pl")) {
    if (is.null(guessing)) {
      stop("guessing is required when model = '3pl'.", call. = FALSE)
    }
  } else {
    if (!is.null(guessing) && any(guessing != 0)) {
      stop(
        "guessing must be NULL or zero for Rasch and 2PL models.",
        call. = FALSE
      )
    }
    guessing <- rep(0, n_items)
  }

  if (identical(model, "rasch") && any(lambda_base != 1)) {
    stop("lambda_base must equal 1 for every Rasch item.", call. = FALSE)
  }

  if (is.null(guessing)) {
    guessing <- rep(0, n_items)
  }

  d0 <- -lambda_base * beta
  if (any(!is.finite(d0))) {
    stop(
      "beta and lambda_base imply non-finite slope/intercept parameters.",
      call. = FALSE
    )
  }

  structure(
    list(
      model = model,
      beta = beta,
      lambda_base = lambda_base,
      guessing = guessing,
      A0 = matrix(lambda_base, ncol = 1L),
      d0 = d0,
      lower = guessing,
      upper = rep(1, n_items),
      D = 1,
      n_items = n_items
    ),
    class = "irtsimrel_item_bank"
  )
}


#' Internal: Validate Model-Engine Evaluation Inputs
#'
#' @noRd
.irtsimrel_engine_inputs <- function(theta, bank, scale, chunk_size) {
  if (!is.numeric(theta) || length(theta) == 0L ||
      anyNA(theta) || any(!is.finite(theta))) {
    stop("theta must be a non-empty finite numeric vector.", call. = FALSE)
  }
  if (!inherits(bank, "irtsimrel_item_bank") ||
      !is.list(bank) || is.null(bank$n_items)) {
    stop("bank must be created by .irtsimrel_item_bank().", call. = FALSE)
  }
  if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) ||
      !is.finite(scale) || scale < 0) {
    stop("scale must be a finite nonnegative numeric scalar.", call. = FALSE)
  }

  if (is.null(chunk_size)) {
    chunk_size <- length(theta)
  } else {
    if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
        is.na(chunk_size) || !is.finite(chunk_size) || chunk_size < 1 ||
        chunk_size != floor(chunk_size)) {
      stop("chunk_size must be NULL or a positive integer.", call. = FALSE)
    }
    chunk_size <- as.integer(min(chunk_size, length(theta)))
  }

  list(
    theta = as.numeric(theta),
    bank = bank,
    scale = as.numeric(scale),
    chunk_size = chunk_size
  )
}


#' Internal: Apply a Model-Engine Kernel in Row Chunks
#'
#' @noRd
.irtsimrel_chunked_matrix <- function(theta, n_items, chunk_size, kernel) {
  out <- matrix(NA_real_, nrow = length(theta), ncol = n_items)
  starts <- seq.int(1L, length(theta), by = chunk_size)

  for (start in starts) {
    end <- min(length(theta), start + chunk_size - 1L)
    rows <- start:end
    out[rows, ] <- kernel(theta[rows])
  }

  out
}


#' Internal: Stable Softplus
#'
#' @noRd
.irtsimrel_softplus <- function(x) {
  positive <- x > 0
  out <- numeric(length(x))
  out[positive] <- x[positive] + log1p(exp(-x[positive]))
  out[!positive] <- log1p(exp(x[!positive]))
  dim(out) <- dim(x)
  out
}


#' Internal: Linear Predictors for a Normalized Bank
#'
#' @noRd
.irtsimrel_eta_matrix <- function(theta, bank, scale) {
  if (scale == 0) {
    return(matrix(0, nrow = length(theta), ncol = bank$n_items))
  }

  delta <- outer(theta, bank$beta, FUN = "-")
  effective_lambda <- scale * bank$lambda_base
  eta <- sweep(
    delta,
    MARGIN = 2L,
    STATS = effective_lambda,
    FUN = "*"
  )

  # Recover products whose effective slope overflowed or underflowed before it
  # could be combined with a compensating theta difference, as well as
  # differences that overflowed before multiplication by a tiny slope.
  # Ordinary cells retain the direct multiply path and its historical rounding
  # behavior.
  unstable <- !is.finite(effective_lambda) | effective_lambda == 0
  recover <- !is.finite(eta)
  if (any(unstable)) recover[, unstable] <- TRUE
  if (any(recover)) {
    log_abs_delta <- log(abs(delta))
    nonfinite_delta <- !is.finite(delta)
    if (any(nonfinite_delta)) {
      for (item in which(colSums(nonfinite_delta) > 0L)) {
        rows <- which(nonfinite_delta[, item])
        abs_theta <- abs(theta[rows])
        abs_beta <- abs(bank$beta[[item]])
        high <- pmax(abs_theta, abs_beta)
        low <- pmin(abs_theta, abs_beta)
        same_sign <- sign(theta[rows]) == sign(bank$beta[[item]])
        adjustment <- ifelse(
          same_sign,
          log1p(-low / high),
          log1p(low / high)
        )
        log_abs_delta[rows, item] <- log(high) + adjustment
      }
    }
    log_abs_eta <- sweep(
      log_abs_delta,
      MARGIN = 2L,
      STATS = log(scale) + log(bank$lambda_base),
      FUN = "+"
    )
    recovered_eta <- sign(delta) * exp(log_abs_eta)
    eta[recover] <- recovered_eta[recover]
  }

  eta
}


#' Internal: Response Probabilities for a Normalized Bank
#'
#' @param theta Numeric ability vector.
#' @param bank Normalized bank from `.irtsimrel_item_bank()`.
#' @param scale Nonnegative global discrimination scale.
#' @param chunk_size Optional number of theta rows evaluated per chunk.
#'
#' @return A length(theta)-by-n_items probability matrix.
#' @noRd
.irtsimrel_probability <- function(theta,
                                   bank,
                                   scale = 1,
                                   chunk_size = NULL) {
  inputs <- .irtsimrel_engine_inputs(theta, bank, scale, chunk_size)
  theta <- inputs$theta
  bank <- inputs$bank
  scale <- inputs$scale

  .irtsimrel_chunked_matrix(
    theta = theta,
    n_items = bank$n_items,
    chunk_size = inputs$chunk_size,
    kernel = function(theta_chunk) {
      q <- plogis(.irtsimrel_eta_matrix(theta_chunk, bank, scale))

      # This branch deliberately preserves exact 2PL behavior when g = 0.
      if (all(bank$guessing == 0)) {
        return(q)
      }

      probability <- q
      has_guessing <- bank$guessing > 0
      probability[, has_guessing] <- sweep(
        q[, has_guessing, drop = FALSE],
        MARGIN = 2L,
        STATS = 1 - bank$guessing[has_guessing],
        FUN = "*"
      )
      probability[, has_guessing] <- sweep(
        probability[, has_guessing, drop = FALSE],
        MARGIN = 2L,
        STATS = bank$guessing[has_guessing],
        FUN = "+"
      )
      probability
    }
  )
}


#' Internal: Stable Log Item Fisher Information
#'
#' @inheritParams .irtsimrel_probability
#'
#' @return A length(theta)-by-n_items log-information matrix.
#' @noRd
.irtsimrel_log_item_information <- function(theta,
                                            bank,
                                            scale = 1,
                                            chunk_size = NULL) {
  inputs <- .irtsimrel_engine_inputs(theta, bank, scale, chunk_size)
  theta <- inputs$theta
  bank <- inputs$bank
  scale <- inputs$scale

  if (scale == 0) {
    return(matrix(-Inf, nrow = length(theta), ncol = bank$n_items))
  }

  log_effective_lambda_sq <- 2 * (log(scale) + log(bank$lambda_base))

  .irtsimrel_chunked_matrix(
    theta = theta,
    n_items = bank$n_items,
    chunk_size = inputs$chunk_size,
    kernel = function(theta_chunk) {
      eta <- .irtsimrel_eta_matrix(theta_chunk, bank, scale)
      log_q <- -.irtsimrel_softplus(-eta)
      log_one_minus_q <- -.irtsimrel_softplus(eta)

      # Exact 2PL/Rasch identity: log(a^2 q (1-q)).
      log_information <- sweep(
        log_q + log_one_minus_q,
        MARGIN = 2L,
        STATS = log_effective_lambda_sq,
        FUN = "+"
      )

      has_guessing <- bank$guessing > 0
      if (!any(has_guessing)) {
        return(log_information)
      }

      g <- bank$guessing[has_guessing]
      log_q_g <- log_q[, has_guessing, drop = FALSE]
      log_one_minus_q_g <- log_one_minus_q[, has_guessing, drop = FALSE]

      # For the 3PL,
      #   I = a^2 (1-g) q^2 (1-q) / {g + (1-g)q}.
      # log(P) is evaluated with logaddexp to avoid loss in the lower tail.
      log_g <- log(g)
      log_one_minus_g <- log1p(-g)
      log_scaled_q <- sweep(
        log_q_g,
        MARGIN = 2L,
        STATS = log_one_minus_g,
        FUN = "+"
      )
      log_probability <- sweep(
        log_scaled_q,
        MARGIN = 2L,
        STATS = log_g,
        FUN = function(x, y) pmax(x, y) + log1p(exp(-abs(x - y)))
      )

      log_information[, has_guessing] <-
        sweep(
          2 * log_q_g + log_one_minus_q_g - log_probability,
          MARGIN = 2L,
          STATS = log_effective_lambda_sq[has_guessing] + log_one_minus_g,
          FUN = "+"
        )

      log_information
    }
  )
}


#' Internal: Row Log-Sum-Exp
#'
#' @noRd
.irtsimrel_row_logsumexp <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  row_max <- apply(x, 1L, max)
  out <- rep(-Inf, nrow(x))
  finite_rows <- is.finite(row_max)

  if (any(finite_rows)) {
    centered <- sweep(
      x[finite_rows, , drop = FALSE],
      MARGIN = 1L,
      STATS = row_max[finite_rows],
      FUN = "-"
    )
    out[finite_rows] <- row_max[finite_rows] + log(rowSums(exp(centered)))
  }

  out
}


#' Internal: Default Row Budget for Reliability Reduction
#'
#' Reliability reducers only return one test-information value per theta node,
#' so retaining a full M-by-I item-information matrix is unnecessary.  Bound
#' each transient item matrix to roughly 250,000 cells unless the caller asks
#' for an explicit chunk size.  Matrix-returning kernels keep their historical
#' behavior because their full output is part of the contract.
#'
#' @noRd
.irtsimrel_reducer_chunk_size <- function(n_nodes,
                                          n_items,
                                          max_cells = 250000L) {
  if (!is.numeric(n_nodes) || length(n_nodes) != 1L || is.na(n_nodes) ||
      !is.finite(n_nodes) || n_nodes < 1 || n_nodes != floor(n_nodes) ||
      !is.numeric(n_items) || length(n_items) != 1L || is.na(n_items) ||
      !is.finite(n_items) || n_items < 1 || n_items != floor(n_items) ||
      !is.numeric(max_cells) || length(max_cells) != 1L || is.na(max_cells) ||
      !is.finite(max_cells) || max_cells < 1 || max_cells != floor(max_cells)) {
    stop(
      "n_nodes, n_items, and max_cells must be positive integer scalars.",
      call. = FALSE
    )
  }
  rows <- max(1, floor(max_cells / n_items))
  as.integer(min(n_nodes, rows))
}


#' Internal: Stable Log Test Fisher Information
#'
#' @inheritParams .irtsimrel_probability
#'
#' @return A numeric vector of log test information, one value per theta.
#' @noRd
.irtsimrel_log_test_information <- function(theta,
                                            bank,
                                            scale = 1,
                                            chunk_size = NULL) {
  chunk_size_requested <- chunk_size
  inputs <- .irtsimrel_engine_inputs(theta, bank, scale, chunk_size)
  theta <- inputs$theta
  bank <- inputs$bank
  scale <- inputs$scale
  chunk_size <- if (is.null(chunk_size_requested)) {
    .irtsimrel_reducer_chunk_size(length(theta), bank$n_items)
  } else {
    inputs$chunk_size
  }

  out <- rep(-Inf, length(theta))
  starts <- seq.int(1L, length(theta), by = chunk_size)

  for (start in starts) {
    end <- min(length(theta), start + chunk_size - 1L)
    rows <- start:end
    log_item_information <- .irtsimrel_log_item_information(
      theta = theta[rows],
      bank = bank,
      scale = scale,
      chunk_size = length(rows)
    )
    out[rows] <- .irtsimrel_row_logsumexp(log_item_information)
  }

  out
}


#' Internal: Test Fisher Information on the Natural Scale
#'
#' This compatibility helper intentionally performs no flooring. Values whose
#' exact natural-scale representation is below double precision may underflow
#' to zero; callers that need tail-stable reductions should consume the log
#' information returned by `.irtsimrel_log_test_information()`.
#'
#' @inheritParams .irtsimrel_probability
#'
#' @return A numeric vector of test information, one value per theta.
#' @noRd
.irtsimrel_test_information <- function(theta,
                                        bank,
                                        scale = 1,
                                        chunk_size = NULL) {
  exp(
    .irtsimrel_log_test_information(
      theta = theta,
      bank = bank,
      scale = scale,
      chunk_size = chunk_size
    )
  )
}
