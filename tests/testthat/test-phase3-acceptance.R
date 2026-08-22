# Independent Phase 3 acceptance tests
#
# These tests deliberately exercise the model engine and reliability reducers
# from hand calculations rather than reusing package implementation helpers.

test_that("Rasch and 2PL kernels match hand probability and information", {
  theta <- c(-1, 0.5)
  beta <- c(-0.25, 1.25)

  rasch_bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = beta,
    lambda_base = c(1, 1),
    model = "rasch"
  )
  rasch_eta <- outer(theta, beta, "-")
  rasch_probability <- plogis(rasch_eta)
  rasch_item_info <- rasch_probability * (1 - rasch_probability)

  expect_equal(
    IRTsimrel:::.irtsimrel_probability(theta, rasch_bank),
    rasch_probability,
    tolerance = 1e-15
  )
  expect_equal(
    exp(IRTsimrel:::.irtsimrel_log_item_information(theta, rasch_bank)),
    rasch_item_info,
    tolerance = 1e-15
  )
  expect_equal(
    IRTsimrel:::.irtsimrel_test_information(theta, rasch_bank),
    rowSums(rasch_item_info),
    tolerance = 1e-15
  )

  lambda <- c(0.8, 1.7)
  scale <- 1.3
  twopl_bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda,
    model = "2pl"
  )
  effective_lambda <- scale * lambda
  twopl_eta <- sweep(outer(theta, beta, "-"), 2, effective_lambda, "*")
  twopl_probability <- plogis(twopl_eta)
  twopl_item_info <- sweep(
    twopl_probability * (1 - twopl_probability),
    2,
    effective_lambda^2,
    "*"
  )

  expect_equal(
    IRTsimrel:::.irtsimrel_probability(theta, twopl_bank, scale = scale),
    twopl_probability,
    tolerance = 1e-15
  )
  expect_equal(
    exp(IRTsimrel:::.irtsimrel_log_item_information(
      theta,
      twopl_bank,
      scale = scale
    )),
    twopl_item_info,
    tolerance = 1e-14
  )
})

test_that("3PL information matches a finite-difference Bernoulli oracle", {
  theta <- c(-1.2, 0.2, 1.1)
  beta <- 0.2
  lambda <- 1.4
  guessing <- 0.25
  scale <- 0.9
  h <- 1e-6
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda,
    guessing = guessing,
    model = "3pl"
  )

  probability <- drop(IRTsimrel:::.irtsimrel_probability(
    theta,
    bank,
    scale = scale
  ))
  probability_plus <- drop(IRTsimrel:::.irtsimrel_probability(
    theta + h,
    bank,
    scale = scale
  ))
  probability_minus <- drop(IRTsimrel:::.irtsimrel_probability(
    theta - h,
    bank,
    scale = scale
  ))
  derivative_fd <- (probability_plus - probability_minus) / (2 * h)
  information_fd <- derivative_fd^2 / (probability * (1 - probability))
  information_kernel <- drop(exp(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank, scale = scale)
  ))

  expect_equal(information_kernel, information_fd, tolerance = 2e-9)

  effective_lambda <- scale * lambda
  probability_at_beta <- (1 + guessing) / 2
  closed_form_at_beta <- effective_lambda^2 * (1 - guessing)^2 /
    (16 * probability_at_beta * (1 - probability_at_beta))
  expect_equal(information_kernel[2], closed_form_at_beta, tolerance = 1e-14)
})

test_that("extreme eta retains bounded probabilities and meaningful log diagnostics", {
  theta <- c(-1000, 1000)
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = 0,
    lambda_base = 1,
    guessing = 0.2,
    model = "3pl"
  )

  probability <- drop(IRTsimrel:::.irtsimrel_probability(theta, bank))
  log_item_info <- drop(IRTsimrel:::.irtsimrel_log_item_information(theta, bank))
  log_test_info <- IRTsimrel:::.irtsimrel_log_test_information(theta, bank)

  expect_equal(probability, c(0.2, 1), tolerance = 0)
  expect_true(all(is.finite(log_item_info)))
  expect_true(all(is.finite(log_test_info)))
  expect_true(all(log_item_info < -900))

  diagnostics <- compute_rho_both(
    1,
    theta,
    beta_vec = 0,
    lambda_base = 1,
    theta_var = 1,
    guessing = 0.2,
    return_diagnostics = TRUE
  )
  expect_true(all(is.finite(diagnostics$rho_logit)))
  expect_true(all(is.finite(diagnostics$log_info_range)))
  expect_identical(diagnostics$below_legacy_floor_count, 2L)
  expect_identical(diagnostics$status, "ok")
})

test_that("zero guessing is an exact 2PL reduction", {
  theta <- seq(-4, 4, length.out = 29)
  beta <- c(-1.1, 0.3, 1.7)
  lambda <- c(0.7, 1.2, 1.8)
  scale <- 1.15
  twopl <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, model = "2pl")
  threepl_zero <- IRTsimrel:::.irtsimrel_item_bank(
    beta,
    lambda,
    guessing = rep(0, length(beta)),
    model = "3pl"
  )

  expect_identical(
    IRTsimrel:::.irtsimrel_probability(theta, twopl, scale),
    IRTsimrel:::.irtsimrel_probability(theta, threepl_zero, scale)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_item_information(theta, twopl, scale),
    IRTsimrel:::.irtsimrel_log_item_information(theta, threepl_zero, scale)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_test_information(theta, twopl, scale),
    IRTsimrel:::.irtsimrel_log_test_information(theta, threepl_zero, scale)
  )

  expect_identical(
    compute_rho_both(scale, theta, beta, lambda),
    compute_rho_both(scale, theta, beta, lambda, guessing = 0)
  )
})

test_that("positive guessing cannot increase pointwise Fisher information", {
  theta <- seq(-8, 8, length.out = 101)
  beta <- c(-1, 0.75)
  lambda <- c(0.9, 1.6)
  twopl <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, model = "2pl")
  threepl <- IRTsimrel:::.irtsimrel_item_bank(
    beta,
    lambda,
    guessing = c(0.1, 0.3),
    model = "3pl"
  )

  log_info_2pl <- IRTsimrel:::.irtsimrel_log_item_information(theta, twopl)
  log_info_3pl <- IRTsimrel:::.irtsimrel_log_item_information(theta, threepl)
  expect_true(all(log_info_3pl <= log_info_2pl + 1e-14))
  expect_true(any(log_info_3pl < log_info_2pl - 1e-6))
})

test_that("chunked and unchunked model kernels are equivalent", {
  theta <- seq(-3.5, 3.5, length.out = 31)
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = c(-1.2, 0.1, 1.5),
    lambda_base = c(0.65, 1.1, 1.9),
    guessing = c(0, 0.15, 0.35),
    model = "3pl"
  )

  for (chunk_size in c(1L, 7L, 30L, 100L)) {
    expect_identical(
      IRTsimrel:::.irtsimrel_probability(theta, bank, chunk_size = NULL),
      IRTsimrel:::.irtsimrel_probability(theta, bank, chunk_size = chunk_size)
    )
    expect_identical(
      IRTsimrel:::.irtsimrel_log_item_information(
        theta,
        bank,
        chunk_size = NULL
      ),
      IRTsimrel:::.irtsimrel_log_item_information(
        theta,
        bank,
        chunk_size = chunk_size
      )
    )
    expect_identical(
      IRTsimrel:::.irtsimrel_log_test_information(
        theta,
        bank,
        chunk_size = NULL
      ),
      IRTsimrel:::.irtsimrel_log_test_information(
        theta,
        bank,
        chunk_size = chunk_size
      )
    )
  }
})

test_that("weighted reliability matches a direct hand calculation", {
  theta <- c(-1, 0, 2)
  beta <- c(-0.5, 0.75)
  lambda <- c(0.9, 1.2)
  guessing <- c(0.1, 0.25)
  weights <- c(1, 2, 3)
  weights_normalized <- weights / sum(weights)
  theta_var <- 1.5
  scale <- 0.8

  eta <- sweep(outer(theta, beta, "-"), 2, scale * lambda, "*")
  q <- plogis(eta)
  probability <- sweep(q, 2, 1 - guessing, "*")
  probability <- sweep(probability, 2, guessing, "+")
  derivative <- sweep(q * (1 - q), 2, (1 - guessing) * scale * lambda, "*")
  item_info <- derivative^2 / (probability * (1 - probability))
  test_info <- rowSums(item_info)
  mean_info <- sum(weights_normalized * test_info)
  mean_inverse_info <- sum(weights_normalized / test_info)
  rho_tilde_hand <- theta_var * mean_info / (1 + theta_var * mean_info)
  rho_bar_hand <- theta_var / (theta_var + mean_inverse_info)

  both <- compute_rho_both(
    scale,
    theta,
    beta,
    lambda,
    theta_var = theta_var,
    guessing = guessing,
    weights = weights
  )
  expect_equal(both$rho_tilde, rho_tilde_hand, tolerance = 1e-14)
  expect_equal(both$rho_bar, rho_bar_hand, tolerance = 1e-14)
  expect_gte(both$rho_tilde, both$rho_bar)
})

test_that("strict invalid inputs error while zero scale remains defined", {
  theta <- c(-1, 0, 1)
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.2)

  expect_error(compute_rho_tilde(-0.1, theta, beta, lambda), "non-negative")
  expect_error(compute_rho_bar(-0.1, theta, beta, lambda), "non-negative")
  expect_error(compute_rho_both(-0.1, theta, beta, lambda), "non-negative")

  bank <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, model = "2pl")
  expect_error(
    IRTsimrel:::.irtsimrel_probability(theta, bank, scale = -0.1),
    "nonnegative"
  )

  for (bad_variance in list(0, -1, Inf, NA_real_)) {
    expect_error(
      compute_rho_tilde(
        1,
        theta,
        beta,
        lambda,
        theta_var = bad_variance
      ),
      "theta_var"
    )
  }

  expect_identical(compute_rho_tilde(0, theta, beta, lambda), 0)
  expect_identical(
    compute_rho_both(0, theta, beta, lambda),
    list(rho_tilde = 0, rho_bar = 0)
  )
})

test_that("public reliability defaults preserve legacy shapes and names", {
  expected_prefix <- c(
    "c",
    "theta_vec",
    "beta_vec",
    "lambda_base",
    "theta_var"
  )
  expect_identical(names(formals(compute_rho_tilde))[1:5], expected_prefix)
  expect_identical(names(formals(compute_rho_bar))[1:5], expected_prefix)
  expect_identical(names(formals(compute_rho_both))[1:5], expected_prefix)

  theta <- c(-1, -0.2, 0.4, 1.1)
  beta <- c(-0.4, 0.8)
  lambda <- c(0.9, 1.3)
  rho_tilde <- compute_rho_tilde(1, theta, beta, lambda)
  rho_bar <- compute_rho_bar(1, theta, beta, lambda)
  both <- compute_rho_both(1, theta, beta, lambda)

  expect_type(rho_tilde, "double")
  expect_length(rho_tilde, 1L)
  expect_null(names(rho_tilde))
  expect_type(rho_bar, "double")
  expect_length(rho_bar, 1L)
  expect_null(names(rho_bar))
  expect_identical(names(both), c("rho_tilde", "rho_bar"))
  expect_true(all(vapply(both, is.double, logical(1))))
  expect_gte(both$rho_tilde, both$rho_bar)
})
