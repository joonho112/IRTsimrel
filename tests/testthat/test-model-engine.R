# Tests for the internal common unidimensional model engine.

test_that("item-bank adapter normalizes public parameters", {
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = c(-1, 0.5),
    lambda_base = c(0.8, 1.2),
    guessing = c(0.1, 0.25),
    model = "3PL"
  )

  expect_s3_class(bank, "irtsimrel_item_bank")
  expect_identical(bank$model, "3pl")
  expect_equal(bank$beta, c(-1, 0.5))
  expect_equal(bank$lambda_base, c(0.8, 1.2))
  expect_equal(bank$guessing, c(0.1, 0.25))
  expect_equal(bank$A0, matrix(c(0.8, 1.2), ncol = 1))
  expect_equal(bank$d0, c(0.8, -0.6))
  expect_equal(bank$lower, bank$guessing)
  expect_equal(bank$upper, c(1, 1))
  expect_identical(bank$D, 1)
  expect_identical(bank$n_items, 2L)
})

test_that("item-bank adapter enforces model and parameter contracts", {
  adapter <- IRTsimrel:::.irtsimrel_item_bank

  expect_error(adapter(numeric(), 1), "non-empty")
  expect_error(adapter(c(0, NA), c(1, 1)), "finite")
  expect_error(adapter(c(0, 1), 1), "same length")
  expect_error(adapter(0, 0), "all values > 0", fixed = TRUE)
  expect_error(adapter(0, -1), "all values > 0", fixed = TRUE)
  expect_error(adapter(0, Inf), "finite")
  expect_error(adapter(0, 1, model = "4pl"), "rasch")
  expect_error(adapter(0, 1, model = "3pl"), "guessing is required")
  expect_error(adapter(0, 1, guessing = 1, model = "3pl"), "0 <= guessing < 1")
  expect_error(adapter(0, 1, guessing = -0.1, model = "3pl"), "0 <= guessing < 1")
  expect_error(
    adapter(c(0, 1), c(1, 1), guessing = c(0.1, 0.2, 0.3), model = "3pl"),
    "length 1 or"
  )
  expect_error(adapter(0, 1, guessing = 0.2, model = "2pl"), "NULL or zero")
  expect_error(adapter(0, 1.1, model = "rasch"), "must equal 1")

  expect_equal(
    adapter(c(-1, 1), c(1, 1), guessing = 0.2, model = "3pl")$guessing,
    c(0.2, 0.2)
  )
})

test_that("probabilities implement the Rasch, 2PL, and 3PL formulas", {
  theta <- c(-1, 0, 1)
  beta <- c(-0.5, 0.5)
  lambda <- c(0.8, 1.4)
  scale <- 1.7

  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, model = "2pl")
  observed_2pl <- IRTsimrel:::.irtsimrel_probability(theta, bank_2pl, scale)
  expected_q <- plogis(scale * outer(theta, beta, "-") *
    matrix(lambda, nrow = length(theta), ncol = length(lambda), byrow = TRUE))
  expect_equal(observed_2pl, expected_q, tolerance = 1e-15)

  g <- c(0.15, 0.25)
  bank_3pl <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, g, "3pl")
  observed_3pl <- IRTsimrel:::.irtsimrel_probability(theta, bank_3pl, scale)
  expected_3pl <- sweep(expected_q, 2, 1 - g, "*")
  expected_3pl <- sweep(expected_3pl, 2, g, "+")
  expect_equal(observed_3pl, expected_3pl, tolerance = 1e-15)

  bank_rasch <- IRTsimrel:::.irtsimrel_item_bank(beta, c(1, 1), model = "rasch")
  expect_equal(
    IRTsimrel:::.irtsimrel_probability(theta, bank_rasch),
    plogis(outer(theta, beta, "-")),
    tolerance = 1e-15
  )
})

test_that("g = 0 follows the exact 2PL path", {
  theta <- c(-1000, -2, 0, 2, 1000)
  beta <- c(-0.3, 0.7)
  lambda <- c(0.9, 1.3)
  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, model = "2pl")
  bank_3pl_zero <- IRTsimrel:::.irtsimrel_item_bank(
    beta, lambda, guessing = c(0, 0), model = "3pl"
  )

  expect_identical(
    IRTsimrel:::.irtsimrel_probability(theta, bank_3pl_zero, scale = 2),
    IRTsimrel:::.irtsimrel_probability(theta, bank_2pl, scale = 2)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank_3pl_zero, scale = 2),
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank_2pl, scale = 2)
  )
})

test_that("extreme predictors remain stable without an information floor", {
  theta <- c(-1000, 1000)
  bank_2pl <- IRTsimrel:::.irtsimrel_item_bank(0, 1, model = "2pl")
  p_2pl <- IRTsimrel:::.irtsimrel_probability(theta, bank_2pl)
  log_i_2pl <- IRTsimrel:::.irtsimrel_log_item_information(theta, bank_2pl)

  expect_equal(as.numeric(p_2pl), c(0, 1))
  expect_true(all(is.finite(log_i_2pl)))
  expect_equal(as.numeric(log_i_2pl), c(-1000, -1000), tolerance = 1e-12)
  expect_true(all(log_i_2pl < log(1e-10)))

  bank_3pl <- IRTsimrel:::.irtsimrel_item_bank(0, 1, 0.2, "3pl")
  p_3pl <- IRTsimrel:::.irtsimrel_probability(theta, bank_3pl)
  log_i_3pl <- IRTsimrel:::.irtsimrel_log_item_information(theta, bank_3pl)

  expect_equal(as.numeric(p_3pl), c(0.2, 1), tolerance = 1e-15)
  expect_true(all(is.finite(log_i_3pl)))
  expect_true(log_i_3pl[1] < -1900)
  expect_equal(log_i_3pl[2], -1000 + log(0.8), tolerance = 1e-12)
})


test_that("finite slope factors avoid intermediate overflow and underflow", {
  overflow_bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = 0,
    lambda_base = 1e308,
    model = "2pl"
  )
  probability <- IRTsimrel:::.irtsimrel_probability(
    theta = 2,
    bank = overflow_bank,
    scale = 1e-308
  )
  log_information <- IRTsimrel:::.irtsimrel_log_item_information(
    theta = 2,
    bank = overflow_bank,
    scale = 1e-308
  )
  expected_eta <- 2
  expect_equal(as.numeric(probability), plogis(expected_eta), tolerance = 1e-15)
  expect_equal(
    as.numeric(log_information),
    2 * (log(1e-308) + log(1e308)) -
      IRTsimrel:::.irtsimrel_softplus(-expected_eta) -
      IRTsimrel:::.irtsimrel_softplus(expected_eta),
    tolerance = 1e-12
  )

  underflow_bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = 0,
    lambda_base = 1e-200,
    model = "2pl"
  )
  underflow_log_information <- IRTsimrel:::.irtsimrel_log_item_information(
    theta = 0,
    bank = underflow_bank,
    scale = 1e-200
  )
  expect_equal(
    as.numeric(underflow_log_information),
    2 * (log(1e-200) + log(1e-200)) - log(4),
    tolerance = 1e-12
  )
  expect_true(is.finite(underflow_log_information))

  compensated_bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = -1e308,
    lambda_base = 1e-308,
    model = "2pl"
  )
  compensated_probability <- IRTsimrel:::.irtsimrel_probability(
    theta = 1e308,
    bank = compensated_bank,
    scale = 1
  )
  compensated_log_information <- IRTsimrel:::.irtsimrel_log_item_information(
    theta = 1e308,
    bank = compensated_bank,
    scale = 1
  )
  expect_equal(
    as.numeric(compensated_probability),
    plogis(2),
    tolerance = 2e-14
  )
  expect_equal(
    as.numeric(compensated_log_information),
    2 * log(1e-308) -
      IRTsimrel:::.irtsimrel_softplus(-2) -
      IRTsimrel:::.irtsimrel_softplus(2),
    tolerance = 1e-12
  )
})


test_that("linear predictors are invariant to representable translations", {
  base_bank <- IRTsimrel:::.irtsimrel_item_bank(0, 1, model = "2pl")
  shifted_bank <- IRTsimrel:::.irtsimrel_item_bank(1e16, 1, model = "2pl")

  expect_identical(
    IRTsimrel:::.irtsimrel_probability(2, base_bank),
    IRTsimrel:::.irtsimrel_probability(1e16 + 2, shifted_bank)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_item_information(2, base_bank),
    IRTsimrel:::.irtsimrel_log_item_information(1e16 + 2, shifted_bank)
  )
})

test_that("log information agrees with the analytic Fisher information", {
  theta <- seq(-2, 2, length.out = 11)
  beta <- c(-0.7, 0.4, 1.1)
  lambda <- c(0.8, 1.2, 1.6)
  guessing <- c(0.1, 0.2, 0.3)
  scale <- 1.35
  bank <- IRTsimrel:::.irtsimrel_item_bank(beta, lambda, guessing, "3pl")

  q <- plogis(scale * sweep(outer(theta, beta, "-"), 2, lambda, "*"))
  probability <- sweep(sweep(q, 2, 1 - guessing, "*"), 2, guessing, "+")
  derivative <- sweep(
    q * (1 - q),
    2,
    scale * lambda * (1 - guessing),
    "*"
  )
  expected_item_info <- derivative^2 / (probability * (1 - probability))

  observed_log_item <- IRTsimrel:::.irtsimrel_log_item_information(
    theta, bank, scale
  )
  observed_log_test <- IRTsimrel:::.irtsimrel_log_test_information(
    theta, bank, scale
  )

  expect_equal(exp(observed_log_item), expected_item_info, tolerance = 1e-13)
  expect_equal(exp(observed_log_test), rowSums(expected_item_info), tolerance = 1e-13)
  expect_equal(
    IRTsimrel:::.irtsimrel_test_information(theta, bank, scale),
    rowSums(expected_item_info),
    tolerance = 1e-13
  )
})

test_that("chunked and unchunked model-engine results are equivalent", {
  set.seed(912)
  theta <- rnorm(37)
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    beta = rnorm(7),
    lambda_base = exp(rnorm(7, sd = 0.2)),
    guessing = runif(7, 0.05, 0.3),
    model = "3pl"
  )

  expect_identical(
    IRTsimrel:::.irtsimrel_probability(theta, bank, chunk_size = 4),
    IRTsimrel:::.irtsimrel_probability(theta, bank)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank, chunk_size = 4),
    IRTsimrel:::.irtsimrel_log_item_information(theta, bank)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_log_test_information(theta, bank, chunk_size = 4),
    IRTsimrel:::.irtsimrel_log_test_information(theta, bank)
  )
  expect_identical(
    IRTsimrel:::.irtsimrel_test_information(theta, bank, chunk_size = 4),
    IRTsimrel:::.irtsimrel_test_information(theta, bank)
  )
})

test_that("scale zero and invalid evaluation inputs have explicit behavior", {
  bank <- IRTsimrel:::.irtsimrel_item_bank(c(-1, 1), c(0.8, 1.2), 0.2, "3pl")
  theta <- c(-2, 0, 2)

  expected_at_zero <- matrix(0.6, nrow = 3, ncol = 2)
  expect_equal(
    IRTsimrel:::.irtsimrel_probability(theta, bank, scale = 0),
    expected_at_zero
  )
  expect_equal(
    IRTsimrel:::.irtsimrel_log_test_information(theta, bank, scale = 0),
    rep(-Inf, 3)
  )
  expect_equal(
    IRTsimrel:::.irtsimrel_test_information(theta, bank, scale = 0),
    rep(0, 3)
  )

  probability <- IRTsimrel:::.irtsimrel_probability
  expect_error(probability(numeric(), bank), "non-empty")
  expect_error(probability(c(0, NA), bank), "finite")
  expect_error(probability(theta, list()), "item_bank")
  expect_error(probability(theta, bank, scale = -1), "nonnegative")
  expect_error(probability(theta, bank, scale = Inf), "finite")
  expect_error(probability(theta, bank, chunk_size = 0), "positive integer")
  expect_error(probability(theta, bank, chunk_size = 1.5), "positive integer")
})
