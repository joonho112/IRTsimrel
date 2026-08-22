# =============================================================================
# test-reliability-v03.R
# =============================================================================
# Version 0.3 reliability-contract tests: log-domain reduction, 3PL guessing,
# weighted empirical nodes, diagnostics, and strict validation.
# =============================================================================

test_that("v0.3 reliability wrappers preserve their positional API", {
  expected <- c("c", "theta_vec", "beta_vec", "lambda_base", "theta_var")

  expect_identical(head(names(formals(compute_rho_bar)), 5L), expected)
  expect_identical(head(names(formals(compute_rho_tilde)), 5L), expected)
  expect_identical(head(names(formals(compute_rho_both)), 5L), expected)

  theta <- c(-1, 0, 1)
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.2)
  expect_type(compute_rho_bar(1, theta, beta, lambda), "double")
  expect_type(compute_rho_tilde(1, theta, beta, lambda), "double")
  expect_identical(
    names(compute_rho_both(1, theta, beta, lambda)),
    c("rho_tilde", "rho_bar")
  )
})

test_that("log-domain 2PL reducers match direct moderate-range arithmetic", {
  theta <- c(-1.5, -0.2, 0.7, 1.8)
  beta <- c(-0.8, 0.1, 1.2)
  lambda <- c(0.7, 1.1, 1.4)
  scale <- 1.3
  theta_var <- stats::var(theta)

  eta <- outer(theta, beta, "-")
  slope <- scale * lambda
  probability <- stats::plogis(sweep(eta, 2L, slope, `*`))
  item_info <- sweep(probability * (1 - probability), 2L, slope^2, `*`)
  test_info <- rowSums(item_info)
  expected_info <- theta_var * mean(test_info) /
    (theta_var * mean(test_info) + 1)
  expected_msem <- theta_var /
    (theta_var + mean(1 / test_info))

  expect_equal(
    compute_rho_tilde(scale, theta, beta, lambda, theta_var),
    expected_info,
    tolerance = 1e-13
  )
  expect_equal(
    compute_rho_bar(scale, theta, beta, lambda, theta_var),
    expected_msem,
    tolerance = 1e-13
  )
})

test_that("3PL reliability agrees with the fixed-D Fisher information formula", {
  theta <- c(-2, -0.5, 0.5, 2)
  beta <- c(-0.7, 0.9)
  lambda <- c(0.8, 1.3)
  guessing <- c(0.15, 0.25)
  scale <- 1.4
  theta_var <- 1.2

  eta <- sweep(outer(theta, beta, "-"), 2L, scale * lambda, `*`)
  q <- stats::plogis(eta)
  p <- sweep(q, 2L, 1 - guessing, `*`)
  p <- sweep(p, 2L, guessing, `+`)
  dp <- sweep(q * (1 - q), 2L, (scale * lambda) * (1 - guessing), `*`)
  test_info <- rowSums(dp^2 / (p * (1 - p)))

  expected_info <- stats::plogis(log(theta_var) + log(mean(test_info)))
  expected_msem <- stats::plogis(log(theta_var) - log(mean(1 / test_info)))

  expect_equal(
    compute_rho_tilde(
      scale, theta, beta, lambda,
      theta_var = theta_var,
      guessing = guessing
    ),
    expected_info,
    tolerance = 1e-13
  )
  expect_equal(
    compute_rho_bar(
      scale, theta, beta, lambda,
      theta_var = theta_var,
      guessing = guessing
    ),
    expected_msem,
    tolerance = 1e-13
  )
})

test_that("zero guessing is the exact 2PL reliability special case", {
  theta <- seq(-3, 3, length.out = 51)
  beta <- c(-1, 0, 1)
  lambda <- c(0.75, 1, 1.25)

  expect_identical(
    compute_rho_tilde(1.2, theta, beta, lambda),
    compute_rho_tilde(1.2, theta, beta, lambda, guessing = 0)
  )
  expect_identical(
    compute_rho_bar(1.2, theta, beta, lambda),
    compute_rho_bar(1.2, theta, beta, lambda, guessing = rep(0, 3))
  )
})

test_that("weighted nodes use normalized weights for both estimands", {
  theta <- c(-1, 0, 1)
  beta <- c(-0.4, 0.6)
  lambda <- c(0.9, 1.2)
  weights <- c(2, 3, 5)
  weights_norm <- weights / sum(weights)
  theta_var <- 1

  q <- stats::plogis(outer(theta, beta, "-") *
    matrix(lambda, nrow = length(theta), ncol = length(lambda), byrow = TRUE))
  test_info <- rowSums(sweep(q * (1 - q), 2L, lambda^2, `*`))
  expected_info <- stats::plogis(log(theta_var) + log(sum(weights_norm * test_info)))
  expected_msem <- stats::plogis(log(theta_var) - log(sum(weights_norm / test_info)))

  both <- compute_rho_both(
    1, theta, beta, lambda,
    theta_var = theta_var,
    weights = weights
  )
  expect_equal(both$rho_tilde, expected_info, tolerance = 1e-13)
  expect_equal(both$rho_bar, expected_msem, tolerance = 1e-13)
})


test_that("finite extreme weight ratios remain positive in log reductions", {
  theta <- c(-1000, 0)
  weights <- c(1e-200, 1e200)
  details <- compute_rho_bar(
    1,
    theta,
    beta_vec = 0,
    lambda_base = 1,
    theta_var = 1,
    weights = weights,
    return_diagnostics = TRUE
  )

  log_raw_weights <- log(weights)
  pivot_weights <- max(log_raw_weights)
  log_total_weights <- pivot_weights +
    log(sum(exp(log_raw_weights - pivot_weights)))
  log_weights <- log_raw_weights - log_total_weights
  log_information <- c(-1000, -log(4))
  log_terms <- log_weights - log_information
  pivot_terms <- max(log_terms)
  expected_log_msem <- pivot_terms + log(sum(exp(log_terms - pivot_terms)))

  expect_equal(details$log_msem, expected_log_msem, tolerance = 1e-12)
  expect_equal(details$rho, stats::plogis(-expected_log_msem), tolerance = 1e-12)
  expect_gt(details$rho, 0)
  expect_identical(details$status, "ok")

  # A tiny relative mass can still make an order-one variance contribution
  # when its node is correspondingly far from the weighted mean.
  weighted_variance <- .resolve_theta_var(
    c(1e200, 0),
    weights = weights
  )
  expect_equal(as.numeric(weighted_variance), 1, tolerance = 1e-12)
  expect_identical(attr(weighted_variance, "source"),
                   "weighted_population_variance")
})

test_that("zero-weight nodes do not affect log-domain reductions", {
  theta <- c(-1000, -1, 0, 1, 1000)
  weights <- c(0, 1, 2, 1, 0)
  beta <- c(-0.25, 0.25)
  lambda <- c(0.9, 1.1)

  with_zeros <- compute_rho_both(
    1, theta, beta, lambda,
    theta_var = 1,
    weights = weights
  )
  without_zeros <- compute_rho_both(
    1, theta[2:4], beta, lambda,
    theta_var = 1,
    weights = weights[2:4]
  )

  expect_equal(with_zeros, without_zeros, tolerance = 1e-15)
})

test_that("extreme finite nodes are not replaced by the legacy information floor", {
  theta <- c(-40, 0, 40)
  beta <- 0
  lambda <- 1

  details <- compute_rho_bar(
    1, theta, beta, lambda,
    theta_var = 1,
    return_diagnostics = TRUE
  )

  expect_gt(details$rho, 0)
  expect_lt(details$rho, 1e-15)
  expect_gt(details$log_msem, 35)
  expect_identical(details$below_legacy_floor_count, 2L)
  expect_equal(details$below_legacy_floor_rate, 2 / 3)
})

test_that("diagnostic mode reports estimand and numerical provenance", {
  theta <- c(-2, -1, 0, 1, 2)
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.2)

  details <- compute_rho_tilde(
    1, theta, beta, lambda,
    weights = c(1, 2, 4, 2, 1),
    return_diagnostics = TRUE
  )
  expected_names <- c(
    "rho", "rho_logit", "theta_var", "theta_var_source",
    "log_mean_information", "log_msem", "log_info_range",
    "below_legacy_floor_count", "below_legacy_floor_rate", "status",
    "theta_measure", "approximation", "metric"
  )

  expect_true(all(expected_names %in% names(details)))
  expect_identical(details$theta_var_source, "weighted_population_variance")
  expect_identical(details$theta_measure, "empirical")
  expect_identical(details$approximation, "weighted_nodes")
  expect_identical(details$metric, "info")
  expect_identical(details$status, "ok")
  expect_length(details$log_info_range, 2L)
  expect_equal(details$rho, stats::plogis(details$rho_logit))

  both <- compute_rho_both(
    1, theta, beta, lambda,
    return_diagnostics = TRUE
  )
  expect_equal(unname(both$rho), c(both$rho_tilde, both$rho_bar))
  expect_length(both$rho_logit, 2L)
})

test_that("zero scale is valid but negative scale and invalid inputs error", {
  theta <- c(-1, 0, 1)
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.2)

  expect_identical(compute_rho_bar(0, theta, beta, lambda), 0)
  expect_identical(compute_rho_tilde(0, theta, beta, lambda), 0)
  expect_identical(
    compute_rho_both(0, theta, beta, lambda),
    list(rho_tilde = 0, rho_bar = 0)
  )
  expect_identical(
    compute_rho_bar(0, theta, beta, lambda, return_diagnostics = TRUE)$status,
    "zero_scale"
  )

  expect_error(compute_rho_bar(-1, theta, beta, lambda), "non-negative")
  expect_error(compute_rho_tilde(1, theta, beta, c(1, 0)), "positive")
  expect_error(
    compute_rho_bar(1, theta, beta, lambda, theta_var = NA_real_),
    "theta_var"
  )
  expect_error(
    compute_rho_bar(1, theta, beta, lambda, theta_var = 0),
    "theta_var"
  )
  expect_error(
    compute_rho_bar(1, theta, beta, lambda, guessing = c(0.1, 1)),
    "guessing"
  )
  expect_error(
    compute_rho_bar(1, theta, beta, lambda, weights = c(1, -1, 1)),
    "weights"
  )
  expect_error(
    compute_rho_bar(1, theta, beta, lambda, weights = c(0, 0, 0)),
    "positive finite sum"
  )
})
