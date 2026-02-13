# =============================================================================
# test-sac_calibrate.R
# =============================================================================
# Comprehensive tests for sac_calibrate() -- Algorithm 2: SAC
# All tests use small n_iter (50-80) and M_per_iter (200-300) for speed.
# =============================================================================

# ---------------------------------------------------------------------------
# Basic convergence: c_star in reasonable range, trajectory correct length
# ---------------------------------------------------------------------------

test_that("SAC converges and c_star is in reasonable range", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    n_iter = 80L,
    M_per_iter = 300L,
    M_pre = 3000L,
    seed = 42
  )))

  # c_star should be in a reasonable range
  expect_true(sac$c_star > 0.1)
  expect_true(sac$c_star < 10)
  # Trajectory should have correct length
  expect_length(sac$trajectory, 80)
  expect_length(sac$rho_trajectory, 80)
})


# ---------------------------------------------------------------------------
# EQC warm start: init_method should be "eqc_warm_start"
# ---------------------------------------------------------------------------

test_that("SAC with EQC warm start sets init_method to 'eqc_warm_start'", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 20, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))

  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    c_init = eqc,
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  expect_equal(sac$init_method, "eqc_warm_start")
  expect_equal(sac$c_init, eqc$c_star, tolerance = 1e-6)
})


# ---------------------------------------------------------------------------
# APC init (c_init=NULL): init_method should be "apc_warm_start"
# ---------------------------------------------------------------------------

test_that("SAC with c_init=NULL uses APC init method", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    c_init = NULL,
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  expect_equal(sac$init_method, "apc_warm_start")
})


# ---------------------------------------------------------------------------
# Numeric c_init: init_method should be "user_specified"
# ---------------------------------------------------------------------------

test_that("SAC with numeric c_init=1.0 uses user_specified init", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    c_init = 1.0,
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  expect_equal(sac$init_method, "user_specified")
})


# ---------------------------------------------------------------------------
# Trajectory: length equals n_iter
# ---------------------------------------------------------------------------

test_that("SAC trajectory and rho_trajectory have length n_iter", {
  n_iter <- 60L
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    n_iter = n_iter,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  expect_length(sac$trajectory, n_iter)
  expect_length(sac$rho_trajectory, n_iter)
})


# ---------------------------------------------------------------------------
# Polyak-Ruppert: c_star = mean(trajectory[(burn_in+1):n_iter])
# ---------------------------------------------------------------------------

test_that("c_star equals Polyak-Ruppert average of post-burn-in iterates", {
  n_iter <- 80L
  burn_in <- 40L
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    n_iter = n_iter,
    burn_in = burn_in,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  expected_mean <- mean(sac$trajectory[(burn_in + 1):n_iter])
  expected_mean <- max(sac$c_bounds[1], min(sac$c_bounds[2], expected_mean))
  expect_equal(sac$c_star, expected_mean, tolerance = 1e-10)
})


# ---------------------------------------------------------------------------
# Default metric is "msem"
# ---------------------------------------------------------------------------

test_that("SAC default metric is 'msem'", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  expect_equal(sac$metric, "msem")
})


# ---------------------------------------------------------------------------
# SAC with "info" metric
# ---------------------------------------------------------------------------

test_that("SAC can use metric='info' successfully", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "info",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  expect_equal(sac$metric, "info")
  expect_true(sac$c_star > 0)
})


# ---------------------------------------------------------------------------
# Deprecated spc_calibrate alias
# ---------------------------------------------------------------------------

test_that("spc_calibrate triggers deprecation warning", {
  # .Deprecated produces a warning containing "deprecated"
  # We catch it with expect_warning; other warnings from SAC convergence
  # are allowed via the class argument
  expect_warning(
    suppressMessages(spc_calibrate(
      target_rho = 0.75,
      n_items = 20,
      model = "rasch",
      item_source = "parametric",
      n_iter = 50L,
      M_per_iter = 200L,
      M_pre = 2000L,
      seed = 42
    )),
    "deprecated"
  )
})


# ---------------------------------------------------------------------------
# Return structure: class and fields
# ---------------------------------------------------------------------------

test_that("SAC return object has correct class 'sac_result' and all fields", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  expect_s3_class(sac, "sac_result")

  expected_fields <- c(
    "c_star", "c_final", "target_rho", "achieved_rho",
    "theta_var", "trajectory", "rho_trajectory",
    "metric", "model", "n_items", "n_iter", "burn_in",
    "M_per_iter", "M_pre", "step_params", "c_bounds",
    "c_init", "init_method", "convergence",
    "beta_vec", "lambda_base", "lambda_scaled",
    "items_base", "items_calib", "theta_quad", "call"
  )
  for (f in expected_fields) {
    expect_true(f %in% names(sac), info = paste("Missing field:", f))
  }

  # convergence sub-structure
  expect_true("converged" %in% names(sac$convergence))
  expect_true("mean_first_half" %in% names(sac$convergence))
  expect_true("mean_second_half" %in% names(sac$convergence))
  expect_true("sd_post_burn" %in% names(sac$convergence))
  expect_true("hit_lower_bound" %in% names(sac$convergence))
  expect_true("hit_upper_bound" %in% names(sac$convergence))

  # Numerical checks
  expect_true(sac$c_star > 0)
  expect_true(sac$theta_var > 0)
  expect_length(sac$beta_vec, 15)
  expect_length(sac$lambda_base, 15)
  expect_length(sac$lambda_scaled, 15)
  expect_equal(sac$n_items, 15L)
})


# ---------------------------------------------------------------------------
# Parameter auto-wrapping (Issue #27)
# ---------------------------------------------------------------------------

test_that("SAC auto-wraps shape params in latent_params", {
  expect_message(
    suppressWarnings(sac_calibrate(
      target_rho = 0.75, n_items = 15, model = "rasch",
      item_source = "parametric",
      latent_shape = "bimodal",
      latent_params = list(delta = 0.9),
      n_iter = 50L, M_per_iter = 200L, M_pre = 2000L, seed = 42
    )),
    "Auto-wrapping"
  )
})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

describe("SAC input validation", {

  test_that("target_rho=0 errors", {
    expect_error(
      sac_calibrate(target_rho = 0, n_items = 20, item_source = "parametric"),
      "target_rho"
    )
  })

  test_that("target_rho=1 errors", {
    expect_error(
      sac_calibrate(target_rho = 1, n_items = 20, item_source = "parametric"),
      "target_rho"
    )
  })

  test_that("n_items=0 errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 0, item_source = "parametric"),
      "n_items"
    )
  })

  test_that("n_iter=0 errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                     item_source = "parametric", n_iter = 0),
      "n_iter"
    )
  })

  test_that("M_per_iter=0 errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                     item_source = "parametric", M_per_iter = 0),
      "M_per_iter"
    )
  })

  test_that("invalid c_init type errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                     item_source = "parametric", c_init = "bad"),
      "c_init"
    )
  })

  test_that("burn_in >= n_iter errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                     item_source = "parametric",
                     n_iter = 50, burn_in = 50),
      "burn_in"
    )
  })
})
