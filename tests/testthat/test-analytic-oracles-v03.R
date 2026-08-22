# =============================================================================
# test-analytic-oracles-v03.R
# =============================================================================
# Phase 7 independent mathematical oracles for the D = 1 dichotomous 3PL.
#
# Independence boundary:
#   * the probability oracle below does not call stats::plogis() or any
#     IRTsimrel probability helper;
#   * the finite-difference derivative uses a five-point stencil applied only
#     to that independent probability oracle;
#   * the stable information oracle is assembled from the Bernoulli identity
#       I(theta) = {dP(theta)/dtheta}^2 / {P(theta) [1 - P(theta)]},
#     rather than copying the package's algebraically reduced 3PL expression.
#
# The tests are deterministic and intentionally use internal production
# kernels only as the system under test.
# =============================================================================

.oracle_softplus <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

.oracle_logistic <- function(eta) {
  out <- numeric(length(eta))
  nonnegative <- eta >= 0
  out[nonnegative] <- 1 / (1 + exp(-eta[nonnegative]))
  exp_eta <- exp(eta[!nonnegative])
  out[!nonnegative] <- exp_eta / (1 + exp_eta)
  out
}

.oracle_logspace_add_pair <- function(x, y) {
  pivot <- pmax(x, y)
  pivot + log(exp(x - pivot) + exp(y - pivot))
}

.oracle_probability <- function(theta, beta, effective_a, guessing) {
  eta <- effective_a * (theta - beta)
  q <- .oracle_logistic(eta)
  guessing + (1 - guessing) * q
}

.oracle_probability_from_eta <- function(eta, guessing) {
  q <- .oracle_logistic(eta)
  guessing + (1 - guessing) * q
}

.oracle_log_information_bernoulli <- function(theta,
                                              beta,
                                              effective_a,
                                              guessing) {
  eta <- effective_a * (theta - beta)
  log_q <- -.oracle_softplus(-eta)
  log_one_minus_q <- -.oracle_softplus(eta)
  log_one_minus_g <- log1p(-guessing)

  # Construct P, 1-P, and P' separately from the Bernoulli definition.
  log_probability <- .oracle_logspace_add_pair(
    log(guessing),
    log_one_minus_g + log_q
  )
  log_one_minus_probability <- log_one_minus_g + log_one_minus_q
  log_abs_derivative <- log(effective_a) + log_one_minus_g +
    log_q + log_one_minus_q

  2 * log_abs_derivative - log_probability - log_one_minus_probability
}

.oracle_five_point_derivative <- function(theta,
                                          beta,
                                          effective_a,
                                          guessing,
                                          h) {
  probability_at <- function(x) {
    .oracle_probability(x, beta, effective_a, guessing)
  }
  (
    -probability_at(theta + 2 * h) +
      8 * probability_at(theta + h) -
      8 * probability_at(theta - h) +
      probability_at(theta - 2 * h)
  ) / (12 * h)
}

.oracle_item_matrix <- function(theta, beta, effective_a, guessing) {
  eta <- outer(theta, beta, "-")
  eta <- sweep(eta, 2L, effective_a, `*`)
  q <- matrix(.oracle_logistic(as.vector(eta)), nrow = nrow(eta))
  probability <- sweep(q, 2L, 1 - guessing, `*`)
  probability <- sweep(probability, 2L, guessing, `+`)
  derivative <- sweep(q * (1 - q), 2L, effective_a * (1 - guessing), `*`)
  derivative^2 / (probability * (1 - probability))
}

test_that("g = 0 is an exact 2PL reduction on a randomized grid", {
  set.seed(71001)
  theta <- sort(c(seq(-20, 20, length.out = 81), runif(20, -6, 6)))
  beta <- runif(24, -3, 3)
  lambda <- exp(runif(24, log(0.25), log(3)))
  scale <- 1.7

  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, model = "2pl"
  )
  bank_3pl_zero <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, guessing = rep(0, length(beta)), model = "3pl"
  )

  probability_2pl <- IRTsimrel:::.irtsimrel_probability(
    theta, bank_2pl, scale = scale
  )
  probability_3pl <- IRTsimrel:::.irtsimrel_probability(
    theta, bank_3pl_zero, scale = scale
  )
  log_information_2pl <- IRTsimrel:::.irtsimrel_log_item_information(
    theta, bank_2pl, scale = scale
  )
  log_information_3pl <- IRTsimrel:::.irtsimrel_log_item_information(
    theta, bank_3pl_zero, scale = scale
  )

  expect_identical(probability_3pl, probability_2pl)
  expect_identical(log_information_3pl, log_information_2pl)
  expect_identical(
    IRTsimrel:::.irtsimrel_log_test_information(
      theta, bank_3pl_zero, scale = scale
    ),
    IRTsimrel:::.irtsimrel_log_test_information(theta, bank_2pl, scale = scale)
  )
  expect_identical(
    compute_rho_both(
      scale, theta, beta, lambda,
      theta_var = 1.25,
      guessing = rep(0, length(beta))
    ),
    compute_rho_both(scale, theta, beta, lambda, theta_var = 1.25)
  )

  eta <- sweep(outer(theta, beta, "-"), 2L, scale * lambda, `*`)
  probability_oracle <- matrix(
    .oracle_logistic(as.vector(eta)),
    nrow = length(theta)
  )
  log_information_oracle <-
    sweep(
      -.oracle_softplus(-eta) - .oracle_softplus(eta),
      2L,
      2 * log(scale * lambda),
      `+`
    )

  expect_equal(probability_2pl, probability_oracle, tolerance = 1e-15)
  expect_equal(log_information_2pl, log_information_oracle, tolerance = 2e-13)
})

test_that("theta = beta satisfies the closed-form 3PL midpoint identities", {
  beta <- c(-2, -0.4, 0, 1.2, 3, 4.5)
  lambda <- c(0.3, 0.7, 1, 1.6, 2.5, 4)
  guessing <- c(0, 1e-12, 0.05, 0.2, 0.45, 0.999999)
  scale <- 1.35
  effective_a <- scale * lambda
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, guessing = guessing, model = "3pl"
  )

  probability_diagonal <- diag(IRTsimrel:::.irtsimrel_probability(
    beta, bank, scale = scale
  ))
  log_information_diagonal <- diag(
    IRTsimrel:::.irtsimrel_log_item_information(beta, bank, scale = scale)
  )

  expected_probability <- (1 + guessing) / 2
  expected_information <- effective_a^2 * (1 - guessing) /
    (4 * (1 + guessing))

  expect_equal(probability_diagonal, expected_probability, tolerance = 1e-15)
  expect_equal(
    exp(log_information_diagonal),
    expected_information,
    tolerance = 2e-14
  )
  expect_equal(
    log_information_diagonal,
    .oracle_log_information_bernoulli(
      beta, beta, effective_a, guessing
    ),
    tolerance = 2e-14
  )
})

test_that("extreme eta remains finite in log information without a floor", {
  theta <- c(-1000, -750, -100, 0, 100, 750, 1000)
  guessing <- 0.2
  bank_3pl <- IRTsimrel:::.irtsimrel_item_bank(
    beta = 0,
    lambda_base = 1,
    guessing = guessing,
    model = "3pl"
  )
  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(
    beta = 0,
    lambda_base = 1,
    model = "2pl"
  )

  probability <- drop(IRTsimrel:::.irtsimrel_probability(theta, bank_3pl))
  log_information_3pl <- drop(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank_3pl)
  )
  log_information_2pl <- drop(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank_2pl)
  )
  oracle_3pl <- .oracle_log_information_bernoulli(
    theta, 0, 1, guessing
  )
  oracle_2pl <- .oracle_log_information_bernoulli(theta, 0, 1, 0)

  expect_equal(probability[c(1, 7)], c(guessing, 1), tolerance = 0)
  expect_true(all(probability >= guessing & probability <= 1))
  expect_true(all(is.finite(log_information_3pl)))
  expect_true(all(is.finite(log_information_2pl)))
  expect_equal(log_information_3pl, oracle_3pl, tolerance = 3e-13)
  expect_equal(log_information_2pl, oracle_2pl, tolerance = 3e-13)
  expect_lt(log_information_3pl[1], -1900)
  expect_lt(log_information_3pl[7], -900)
  expect_identical(
    IRTsimrel:::.irtsimrel_test_information(theta[1], bank_3pl),
    0
  )
})

test_that("3PL Fisher information never exceeds its matching 2PL information", {
  theta <- c(-1000, -100, -20, -8, -3, -1, 0, 1, 3, 8, 20, 100, 1000)
  beta <- c(-2.5, -0.5, 0, 1.25, 3)
  lambda <- c(0.25, 0.8, 1, 1.75, 3)
  guessing <- c(1e-12, 0.05, 0.2, 0.45, 0.99)
  scale <- 1.4
  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, model = "2pl"
  )
  bank_3pl <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, guessing = guessing, model = "3pl"
  )

  log_information_2pl <- IRTsimrel:::.irtsimrel_log_item_information(
    theta, bank_2pl, scale = scale
  )
  log_information_3pl <- IRTsimrel:::.irtsimrel_log_item_information(
    theta, bank_3pl, scale = scale
  )

  eta <- sweep(outer(theta, beta, "-"), 2L, scale * lambda, `*`)
  log_q <- -.oracle_softplus(-eta)
  log_probability <- matrix(NA_real_, nrow(eta), ncol(eta))
  for (j in seq_along(guessing)) {
    log_probability[, j] <- .oracle_logspace_add_pair(
      log(guessing[j]),
      log1p(-guessing[j]) + log_q[, j]
    )
  }
  log_ratio_oracle <- sweep(log_q, 2L, log1p(-guessing), `+`) -
    log_probability

  expect_true(all(log_information_3pl <= log_information_2pl + 5e-13))
  expect_true(all(log_ratio_oracle <= 5e-16))
  expect_equal(
    log_information_3pl - log_information_2pl,
    log_ratio_oracle,
    tolerance = 5e-13
  )
})

test_that("five-point finite differences recover the 3PL derivative and information", {
  set.seed(71002)
  n_cases <- 400L
  effective_a <- exp(runif(n_cases, log(0.2), log(3.5)))
  beta <- runif(n_cases, -3, 3)
  eta <- runif(n_cases, -8, 8)
  theta <- beta + eta / effective_a
  guessing <- runif(n_cases, 0, 0.49)
  guessing[seq.int(1L, n_cases, by = 37L)] <- 0
  h <- 2e-4 / pmax(effective_a, 1)

  derivative_fd <- .oracle_five_point_derivative(
    theta, beta, effective_a, guessing, h
  )
  q <- .oracle_logistic(eta)
  derivative_closed <- (1 - guessing) * effective_a * q * (1 - q)
  probability <- .oracle_probability_from_eta(eta, guessing)
  information_fd <- derivative_fd^2 / (probability * (1 - probability))

  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta,
    effective_a,
    guessing = guessing,
    model = "3pl"
  )
  information_production <- exp(diag(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank)
  ))
  probability_production <- diag(
    IRTsimrel:::.irtsimrel_probability(theta, bank)
  )

  derivative_relative_error <- abs(derivative_fd - derivative_closed) /
    pmax(abs(derivative_closed), 1e-14)
  information_relative_error <- abs(information_fd - information_production) /
    pmax(information_production, 1e-14)

  expect_lt(max(derivative_relative_error), 2e-8)
  expect_lt(max(information_relative_error), 4e-8)
  expect_equal(probability_production, probability, tolerance = 3e-16)
})

test_that("stable Bernoulli oracle agrees over randomized edge grids", {
  set.seed(71003)
  theta <- c(-60, -30, -12, -5, -2, -0.5, 0, 0.5, 2, 5, 12, 30, 60)
  n_items <- 31L
  beta <- runif(n_items, -4, 4)
  lambda <- exp(runif(n_items, log(0.15), log(4)))
  guessing <- runif(n_items, 1e-10, 0.499)
  guessing[c(1, 11, 21, 31)] <- 0
  scales <- c(0.03, 0.2, 1, 3, 8)
  max_error <- 0

  for (scale in scales) {
    bank <- IRTsimrel:::.irtsimrel_item_bank(
      beta, lambda, guessing = guessing, model = "3pl"
    )
    observed <- IRTsimrel:::.irtsimrel_log_item_information(
      theta, bank, scale = scale
    )
    expected <- matrix(NA_real_, nrow = length(theta), ncol = n_items)
    for (j in seq_len(n_items)) {
      expected[, j] <- .oracle_log_information_bernoulli(
        theta,
        beta[j],
        scale * lambda[j],
        guessing[j]
      )
    }
    max_error <- max(max_error, abs(observed - expected))
    expect_true(all(is.finite(observed)))
  }

  expect_lt(max_error, 2e-12)
})

test_that("Jensen ordering holds on randomized finite empirical supports", {
  set.seed(71004)
  n_designs <- 120L
  minimum_gap <- Inf
  maximum_formula_error <- 0
  minimum_jensen_margin <- Inf

  for (design_id in seq_len(n_designs)) {
    n_nodes <- sample(5:25, 1)
    n_items <- sample(1:12, 1)
    # Keep this direct natural-scale Jensen oracle away from probability
    # saturation. Extreme tails are covered separately in the log-domain test.
    theta <- sort(runif(n_nodes, -3, 3))
    beta <- runif(n_items, -2, 2)
    lambda <- exp(runif(n_items, log(0.3), log(2)))
    guessing <- runif(n_items, 0, 0.49)
    scale <- exp(runif(1, log(0.15), log(1.5)))
    theta_var <- exp(runif(1, log(0.05), log(5)))
    weights <- rgamma(n_nodes, shape = 1.5)
    if (design_id %% 4L == 0L) weights[1] <- 0
    weights <- weights / sum(weights)

    item_information <- .oracle_item_matrix(
      theta,
      beta,
      scale * lambda,
      guessing
    )
    test_information <- rowSums(item_information)
    mean_information <- sum(weights * test_information)
    mean_inverse_information <- sum(weights / test_information)
    rho_info_oracle <- theta_var * mean_information /
      (1 + theta_var * mean_information)
    rho_msem_oracle <- theta_var /
      (theta_var + mean_inverse_information)

    observed <- compute_rho_both(
      scale,
      theta,
      beta,
      lambda,
      theta_var = theta_var,
      guessing = guessing,
      weights = weights
    )
    minimum_gap <- min(minimum_gap, observed$rho_tilde - observed$rho_bar)
    maximum_formula_error <- max(
      maximum_formula_error,
      abs(observed$rho_tilde - rho_info_oracle),
      abs(observed$rho_bar - rho_msem_oracle)
    )
    minimum_jensen_margin <- min(
      minimum_jensen_margin,
      mean_inverse_information - 1 / mean_information
    )
  }

  expect_gte(minimum_jensen_margin, -2e-13)
  expect_gte(minimum_gap, -2e-14)
  expect_lt(maximum_formula_error, 2e-12)
})
