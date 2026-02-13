# =============================================================================
# test-sim_latentG.R
# =============================================================================
# Comprehensive tests for sim_latentG() and related helper functions
# =============================================================================

# ---------------------------------------------------------------------------
# Pre-standardization verification for all 12 shapes
# Each built-in shape is mathematically constructed to have E[z]=0, Var[z]=1
# ---------------------------------------------------------------------------

describe("Pre-standardization: all 12 shapes have E[z]=0, Var[z]=1", {

  # We use n=10000 and seed=42 for all shapes; tolerance accounts for
  # Monte Carlo sampling variability.
  n_large <- 10000L
  seed    <- 42L
  mean_tol <- 0.05
  sd_tol   <- 0.10

  test_that("normal shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "normal", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("bimodal shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "bimodal", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("trimodal shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "trimodal", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("multimodal shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "multimodal", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("skew_pos shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "skew_pos", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("skew_neg shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "skew_neg", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("heavy_tail shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "heavy_tail", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("light_tail shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "light_tail", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("uniform shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "uniform", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("floor shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "floor", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("ceiling shape: mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(n = n_large, shape = "ceiling", seed = seed)
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })

  test_that("custom shape (standardize_custom=TRUE): mean(z) ~ 0 and sd(z) ~ 1", {
    g <- sim_latentG(
      n = n_large, shape = "custom", seed = seed,
      mixture_spec = list(
        weights = c(0.3, 0.5, 0.2),
        means = c(-1.5, 0, 2),
        sds = c(0.5, 0.7, 0.5)
      ),
      standardize_custom = TRUE
    )
    expect_true(abs(mean(g$z)) < mean_tol)
    expect_true(abs(sd(g$z) - 1) < sd_tol)
  })
})


# ---------------------------------------------------------------------------
# sigma parameter: controls the SD of theta
# ---------------------------------------------------------------------------

test_that("sigma=2 produces sd(theta) approximately 2", {
  g <- sim_latentG(n = 10000L, shape = "normal", sigma = 2, seed = 42)
  expect_true(abs(sd(g$theta) - 2) < 0.15)
})

test_that("sigma=0.5 produces sd(theta) approximately 0.5", {
  g <- sim_latentG(n = 10000L, shape = "normal", sigma = 0.5, seed = 42)
  expect_true(abs(sd(g$theta) - 0.5) < 0.08)
})


# ---------------------------------------------------------------------------
# mu parameter: shifts the mean of theta
# ---------------------------------------------------------------------------

test_that("mu=3 produces mean(theta) approximately 3", {
  g <- sim_latentG(n = 10000L, shape = "normal", mu = 3, seed = 42)
  expect_true(abs(mean(g$theta) - 3) < 0.1)
})

test_that("mu=-5 produces mean(theta) approximately -5", {
  g <- sim_latentG(n = 10000L, shape = "normal", mu = -5, seed = 42)
  expect_true(abs(mean(g$theta) - (-5)) < 0.1)
})


# ---------------------------------------------------------------------------
# Covariate effects via xcov and beta
# ---------------------------------------------------------------------------

test_that("covariate with beta=1 shifts group means by approx 1", {
  n <- 6000L
  set.seed(99)
  group <- c(rep(0L, n / 2), rep(1L, n / 2))
  g <- sim_latentG(
    n = n, shape = "normal",
    xcov = data.frame(group = group),
    beta = 1,
    seed = 42
  )
  grp0 <- mean(g$theta[group == 0])
  grp1 <- mean(g$theta[group == 1])
  expect_true(abs((grp1 - grp0) - 1) < 0.2)
})

test_that("eta_cov equals xcov %*% beta exactly", {
  n <- 100L
  set.seed(10)
  x <- data.frame(v1 = rnorm(n))
  g <- sim_latentG(n = n, shape = "normal", seed = 42,
                   xcov = x, beta = 2)
  expect_equal(g$eta_cov, 2 * x$v1, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Seed reproducibility: same seed yields identical output
# ---------------------------------------------------------------------------

test_that("same seed produces identical theta and z vectors", {
  g1 <- sim_latentG(n = 500L, shape = "bimodal", seed = 123)
  g2 <- sim_latentG(n = 500L, shape = "bimodal", seed = 123)
  expect_identical(g1$theta, g2$theta)
  expect_identical(g1$z, g2$z)
})

test_that("seed save/restore preserves external RNG state", {
  set.seed(999)
  state_before <- .Random.seed
  g <- sim_latentG(n = 100L, shape = "normal", seed = 42)
  state_after <- .Random.seed
  # sim_latentG with explicit seed should restore the external state
  expect_identical(state_before, state_after)
})


# ---------------------------------------------------------------------------
# Custom mixture: standardize_custom = FALSE returns raw values
# ---------------------------------------------------------------------------

test_that("custom with standardize_custom=FALSE returns non-unit-variance z", {
  spec <- list(
    weights = c(0.5, 0.5),
    means = c(-2, 2),
    sds = c(0.5, 0.5)
  )
  g <- sim_latentG(
    n = 10000L, shape = "custom",
    mixture_spec = spec,
    standardize_custom = FALSE,
    seed = 42
  )
  # The raw mixture has theoretical variance = 0.5*(0.25+4) + 0.5*(0.25+4) - 0 = 4.25
  # So sd should be ~ 2.06, NOT 1
  expect_true(sd(g$z) > 1.5)
})


# ---------------------------------------------------------------------------
# Edge cases: bimodal delta extremes
# ---------------------------------------------------------------------------

test_that("bimodal delta=0.05 (nearly unimodal) still has unit variance z", {
  g <- sim_latentG(
    n = 10000L, shape = "bimodal",
    shape_params = list(delta = 0.05),
    seed = 42
  )
  expect_true(abs(mean(g$z)) < 0.05)
  expect_true(abs(sd(g$z) - 1) < 0.1)
})

test_that("bimodal delta=0.95 (strongly bimodal) still has unit variance z", {
  g <- sim_latentG(
    n = 10000L, shape = "bimodal",
    shape_params = list(delta = 0.95),
    seed = 42
  )
  expect_true(abs(mean(g$z)) < 0.05)
  expect_true(abs(sd(g$z) - 1) < 0.1)
})


# ---------------------------------------------------------------------------
# Edge case: heavy_tail with df near minimum (>2)
# ---------------------------------------------------------------------------

test_that("heavy_tail df=2.5 produces finite mean ~ 0", {
  g <- sim_latentG(
    n = 10000L, shape = "heavy_tail",
    shape_params = list(df = 2.5),
    seed = 42
  )
  # Very heavy tails: wider tolerance needed
  expect_true(abs(mean(g$z)) < 0.15)
  # Variance may deviate more in finite samples for df near 2
  expect_true(abs(sd(g$z) - 1) < 0.35)
})


# ---------------------------------------------------------------------------
# Edge case: skew_pos with k=1 (extreme skew)
# ---------------------------------------------------------------------------

test_that("skew_pos k=1 (extreme skew) still has mean ~ 0 and sd ~ 1", {
  g <- sim_latentG(
    n = 10000L, shape = "skew_pos",
    shape_params = list(k = 1),
    seed = 42
  )
  expect_true(abs(mean(g$z)) < 0.1)
  expect_true(abs(sd(g$z) - 1) < 0.15)
})


# ---------------------------------------------------------------------------
# Input validation: error conditions
# ---------------------------------------------------------------------------

describe("Input validation errors in sim_latentG()", {

  test_that("n=0 errors", {
    expect_error(sim_latentG(n = 0), "`n` must be a positive integer")
  })

  test_that("n=-5 errors", {
    expect_error(sim_latentG(n = -5), "`n` must be a positive integer")
  })

  test_that("n=2.5 (non-integer) errors", {
    expect_error(sim_latentG(n = 2.5), "`n` must be a positive integer")
  })

  test_that("n='abc' (character) errors", {
    expect_error(sim_latentG(n = "abc"), "`n` must be a positive integer")
  })

  test_that("sigma=-1 errors", {
    expect_error(sim_latentG(n = 10, sigma = -1), "`sigma` must be a positive scalar")
  })

  test_that("sigma=0 errors", {
    expect_error(sim_latentG(n = 10, sigma = 0), "`sigma` must be a positive scalar")
  })

  test_that("mu=Inf errors", {
    expect_error(sim_latentG(n = 10, mu = Inf), "`mu` must be a finite numeric scalar")
  })

  test_that("mu=NaN errors", {
    expect_error(sim_latentG(n = 10, mu = NaN), "`mu` must be a finite numeric scalar")
  })

  test_that("bimodal delta=0 errors", {
    expect_error(
      sim_latentG(n = 10, shape = "bimodal", shape_params = list(delta = 0)),
      "delta"
    )
  })

  test_that("bimodal delta=1 errors", {
    expect_error(
      sim_latentG(n = 10, shape = "bimodal", shape_params = list(delta = 1)),
      "delta"
    )
  })

  test_that("heavy_tail df=2 errors (must be > 2)", {
    expect_error(
      sim_latentG(n = 10, shape = "heavy_tail", shape_params = list(df = 2)),
      "df"
    )
  })

  test_that("custom without mixture_spec errors", {
    expect_error(
      sim_latentG(n = 10, shape = "custom"),
      "mixture_spec"
    )
  })

  test_that("custom with weights not summing to 1 errors", {
    expect_error(
      sim_latentG(
        n = 10, shape = "custom",
        mixture_spec = list(
          weights = c(0.3, 0.3), means = c(0, 1), sds = c(1, 1)
        )
      ),
      "sum to 1"
    )
  })

  test_that("xcov without beta errors", {
    expect_error(
      sim_latentG(n = 10, shape = "normal", xcov = data.frame(x = 1:10)),
      "beta"
    )
  })

  test_that("trimodal with (1-w0)*m^2 >= 1 errors", {
    expect_error(
      sim_latentG(n = 10, shape = "trimodal",
                  shape_params = list(w0 = 0.1, m = 3)),
      "must have"
    )
  })
})


# ---------------------------------------------------------------------------
# return_z parameter
# ---------------------------------------------------------------------------

test_that("return_z=FALSE omits z; return_z=TRUE includes z", {
  g_no_z <- sim_latentG(n = 100L, shape = "normal", return_z = FALSE, seed = 1)
  g_with_z <- sim_latentG(n = 100L, shape = "normal", return_z = TRUE, seed = 1)

  expect_null(g_no_z$z)
  expect_false(is.null(g_with_z$z))
  expect_length(g_with_z$z, 100L)
})


# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("sim_latentG returns latent_G class with all expected fields", {
  g <- sim_latentG(n = 200L, shape = "bimodal", seed = 1)

  expect_s3_class(g, "latent_G")
  expect_true(all(
    c("theta", "z", "mu", "sigma", "eta_cov", "shape",
      "shape_params", "n", "sample_moments") %in% names(g)
  ))
  expect_length(g$theta, 200L)
  expect_length(g$z, 200L)
  expect_equal(g$n, 200L)
  expect_equal(g$shape, "bimodal")
  expect_equal(g$mu, 0)
  expect_equal(g$sigma, 1)
})

test_that("sample_moments contains mean, sd, skewness, and kurtosis", {
  g <- sim_latentG(n = 100L, shape = "normal", seed = 1)
  sm <- g$sample_moments
  expect_true(all(c("mean", "sd", "skewness", "kurtosis") %in% names(sm)))
  expect_true(is.numeric(sm$mean))
  expect_true(is.numeric(sm$sd))
  expect_true(is.numeric(sm$skewness))
  expect_true(is.numeric(sm$kurtosis))
})

test_that("custom shape result includes mixture_spec in output", {
  g <- sim_latentG(
    n = 100L, shape = "custom", seed = 1,
    mixture_spec = list(
      weights = c(0.5, 0.5), means = c(-1, 1), sds = c(0.5, 0.5)
    )
  )
  expect_true("mixture_spec" %in% names(g))
  expect_equal(g$mixture_spec$weights, c(0.5, 0.5))
})
