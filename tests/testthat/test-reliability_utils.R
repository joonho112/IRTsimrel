# =============================================================================
# test-reliability_utils.R
# =============================================================================
# Comprehensive tests for compute_rho_bar(), compute_rho_tilde(),
# .compute_rho_generic() (via the public wrappers), and compute_apc_init()
# =============================================================================

# ---------------------------------------------------------------------------
# Jensen's inequality: rho_tilde >= rho_bar ALWAYS holds
# ---------------------------------------------------------------------------

describe("Jensen's inequality: rho_tilde >= rho_bar", {

  test_that("inequality holds at c=0.5, 1.0, 2.0 for Rasch items", {
    set.seed(42)
    theta <- rnorm(2000)
    beta  <- rnorm(20)
    lambda <- rep(1, 20)

    for (c_val in c(0.5, 1.0, 2.0)) {
      rho_bar   <- compute_rho_bar(c_val, theta, beta, lambda)
      rho_tilde <- compute_rho_tilde(c_val, theta, beta, lambda)
      expect_true(
        rho_tilde >= rho_bar - 1e-10,
        info = paste("c =", c_val,
                     "| tilde:", round(rho_tilde, 6),
                     "| bar:", round(rho_bar, 6))
      )
    }
  })

  test_that("inequality holds with 2PL (varying lambda)", {
    set.seed(7)
    theta <- rnorm(2000)
    beta  <- rnorm(20)
    lambda_2pl <- exp(rnorm(20, 0, 0.3))

    for (c_val in c(0.5, 1.0, 3.0)) {
      rho_bar   <- compute_rho_bar(c_val, theta, beta, lambda_2pl)
      rho_tilde <- compute_rho_tilde(c_val, theta, beta, lambda_2pl)
      expect_true(rho_tilde >= rho_bar - 1e-10)
    }
  })
})


# ---------------------------------------------------------------------------
# Monotonicity of rho_tilde in c
# ---------------------------------------------------------------------------

test_that("rho_tilde is strictly monotonically increasing in c", {
  set.seed(42)
  theta <- rnorm(2000)
  beta  <- rnorm(20)
  lambda <- rep(1, 20)

  c_vals <- c(0.3, 0.5, 1.0, 2.0, 5.0)
  rho_vals <- sapply(c_vals, function(cv) {
    compute_rho_tilde(cv, theta, beta, lambda)
  })

  for (i in 2:length(rho_vals)) {
    expect_true(
      rho_vals[i] > rho_vals[i - 1],
      info = paste("c =", c_vals[i], "rho =", round(rho_vals[i], 6),
                   "should exceed c =", c_vals[i - 1],
                   "rho =", round(rho_vals[i - 1], 6))
    )
  }
})


# ---------------------------------------------------------------------------
# Small c (near 0): reliability near 0
# ---------------------------------------------------------------------------

test_that("very small c (0.001) gives reliability near 0", {
  set.seed(42)
  theta <- rnorm(2000)
  beta  <- rnorm(20)
  lambda <- rep(1, 20)

  rho_bar   <- compute_rho_bar(0.001, theta, beta, lambda)
  rho_tilde <- compute_rho_tilde(0.001, theta, beta, lambda)

  expect_true(rho_bar < 0.01,
              info = paste("rho_bar at c=0.001:", round(rho_bar, 8)))
  expect_true(rho_tilde < 0.01,
              info = paste("rho_tilde at c=0.001:", round(rho_tilde, 8)))
})


# ---------------------------------------------------------------------------
# Large c: reliability near 1 but < 1
# ---------------------------------------------------------------------------

test_that("large c gives rho_tilde near 1 but strictly < 1", {
  set.seed(42)
  theta <- rnorm(2000)
  beta  <- rnorm(20)
  lambda <- rep(1, 20)

  rho_tilde <- compute_rho_tilde(50, theta, beta, lambda)
  expect_true(rho_tilde > 0.99)
  expect_true(rho_tilde < 1)
})

test_that("rho_bar may decrease at very large c (known non-monotone behavior)", {
  # MSEM-based rho_bar is NON-MONOTONE for large c due to the harmonic mean.

  # At extreme c, p(1-p) -> 0 for extreme theta, making test_info -> 0,
  # causing 1/test_info -> Inf. This is the documented reason why EQC
  # restricts to rho_tilde (the "info" metric) by default.
  set.seed(42)
  theta <- rnorm(2000)
  beta  <- rnorm(20)
  lambda <- rep(1, 20)

  rho_bar_moderate <- compute_rho_bar(1.5, theta, beta, lambda)
  rho_bar_extreme  <- compute_rho_bar(50, theta, beta, lambda)

  # At moderate c, rho_bar should be reasonable
  expect_true(rho_bar_moderate > 0.5)
  # At extreme c, rho_bar may actually decrease (non-monotone!)
  expect_true(rho_bar_extreme < 1)
  expect_true(rho_bar_extreme >= 0)
})


# ---------------------------------------------------------------------------
# c <= 0 returns 0
# ---------------------------------------------------------------------------

test_that("c=0 returns reliability 0", {
  set.seed(42)
  theta <- rnorm(100)
  beta  <- rnorm(10)
  lambda <- rep(1, 10)

  expect_equal(compute_rho_bar(0, theta, beta, lambda), 0)
  expect_equal(compute_rho_tilde(0, theta, beta, lambda), 0)
})

test_that("c=-1 returns reliability 0", {
  set.seed(42)
  theta <- rnorm(100)
  beta  <- rnorm(10)
  lambda <- rep(1, 10)

  expect_equal(compute_rho_bar(-1, theta, beta, lambda), 0)
  expect_equal(compute_rho_tilde(-1, theta, beta, lambda), 0)
})


# ---------------------------------------------------------------------------
# Pre-calculated theta_var matches computed
# ---------------------------------------------------------------------------

test_that("pre-calculated theta_var gives identical result to auto-computed", {
  set.seed(42)
  theta <- rnorm(2000)
  beta  <- rnorm(20)
  lambda <- rep(1, 20)

  tv <- var(theta)
  rho1 <- compute_rho_tilde(1, theta, beta, lambda)
  rho2 <- compute_rho_tilde(1, theta, beta, lambda, theta_var = tv)
  expect_equal(rho1, rho2, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# theta_var = NA or near-zero triggers warning and fallback to 1.0
# ---------------------------------------------------------------------------

test_that("theta_var=NA triggers warning and fallback to 1.0", {
  set.seed(42)
  theta <- rnorm(100)
  beta  <- rnorm(10)
  lambda <- rep(1, 10)

  expect_warning(
    compute_rho_bar(1, theta, beta, lambda, theta_var = NA),
    "Invalid theta_var"
  )
})

test_that("theta_var near-zero triggers warning and fallback", {
  set.seed(42)
  theta <- rnorm(100)
  beta  <- rnorm(10)
  lambda <- rep(1, 10)

  expect_warning(
    compute_rho_bar(1, theta, beta, lambda, theta_var = 1e-12),
    "Invalid theta_var"
  )
})


# ---------------------------------------------------------------------------
# compute_apc_init()
# ---------------------------------------------------------------------------

describe("compute_apc_init()", {

  test_that("result is positive and in [0.1, 10]", {
    c_init <- compute_apc_init(0.80, 25)
    expect_true(c_init >= 0.1)
    expect_true(c_init <= 10)
    expect_true(c_init > 0)
  })

  test_that("higher target requires larger c_init", {
    c_low  <- compute_apc_init(0.50, 25)
    c_high <- compute_apc_init(0.90, 25)
    expect_true(c_high > c_low)
  })

  test_that("more items require smaller c_init for same target", {
    c_few  <- compute_apc_init(0.80, 10)
    c_many <- compute_apc_init(0.80, 50)
    expect_true(c_few > c_many)
  })

  test_that("APC formula is correctly implemented", {
    target <- 0.80
    I <- 25
    sigma_beta <- 1.0
    sigma_sq <- 1 + sigma_beta^2
    kappa <- 0.25 / sqrt(1 + sigma_sq * pi^2 / 3)
    expected <- sqrt(target / (I * kappa * (1 - target)))
    expected <- max(0.1, min(10, expected))

    result <- compute_apc_init(target, I, sigma_beta)
    expect_equal(result, expected, tolerance = 1e-10)
  })

  test_that("APC init is bounded: very high target hits ceiling of 10", {
    c_init <- compute_apc_init(0.999, 3, sigma_beta = 1.0)
    expect_equal(c_init, 10)
  })

  test_that("APC init with many items and low target hits floor of 0.1", {
    c_init <- compute_apc_init(0.01, 200, sigma_beta = 0.1)
    expect_equal(c_init, 0.1)
  })
})


# ---------------------------------------------------------------------------
# Input validation for reliability functions
# ---------------------------------------------------------------------------

describe("reliability utility input validation", {

  test_that("mismatched beta/lambda lengths errors", {
    set.seed(42)
    theta <- rnorm(100)
    expect_error(
      compute_rho_bar(1, theta, rnorm(10), rep(1, 5)),
      "beta_vec and lambda_base must have the same length"
    )
  })

  test_that("theta_vec length 1 errors", {
    expect_error(
      compute_rho_tilde(1, c(0.5), rnorm(10), rep(1, 10)),
      "theta_vec must have length >= 2"
    )
  })

  test_that("empty theta_vec errors", {
    expect_error(
      compute_rho_bar(1, numeric(0), rnorm(10), rep(1, 10)),
      "theta_vec must have length >= 2"
    )
  })
})
