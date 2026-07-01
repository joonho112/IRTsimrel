# =============================================================================
# test-eqc_calibrate.R
# =============================================================================
# Comprehensive tests for eqc_calibrate() -- Algorithm 1: EQC
# All tests use item_source="parametric" and small M for speed.
# =============================================================================

# ---------------------------------------------------------------------------
# Accuracy: Rasch model
# ---------------------------------------------------------------------------

test_that("EQC Rasch: achieved_rho within 0.02 of target 0.80", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80,
    n_items = 25,
    model = "rasch",
    item_source = "parametric",
    M = 5000L,
    seed = 42
  ))
  expect_true(
    abs(eqc$achieved_rho - 0.80) < 0.02,
    info = paste("Achieved rho:", round(eqc$achieved_rho, 5))
  )
  expect_true(eqc$c_star > 0)
})


# ---------------------------------------------------------------------------
# Accuracy: 2PL model
# ---------------------------------------------------------------------------

test_that("EQC 2PL: achieved_rho within 0.03 of target 0.80", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80,
    n_items = 25,
    model = "2pl",
    item_source = "parametric",
    M = 5000L,
    seed = 42
  ))
  expect_true(
    abs(eqc$achieved_rho - 0.80) < 0.03,
    info = paste("Achieved rho:", round(eqc$achieved_rho, 5))
  )
})


# ---------------------------------------------------------------------------
# Default metric "info" does NOT produce warning
# ---------------------------------------------------------------------------

test_that("EQC with default metric 'info' produces no metric-policy warning", {
  # expect_no_warning may not exist in all testthat versions,
  # so we use expect_silent for the warning-specific part.
  # Note: messages are still allowed, so suppress those.
  expect_no_warning(
    suppressMessages(eqc_calibrate(
      target_rho = 0.75,
      n_items = 20,
      model = "rasch",
      item_source = "parametric",
      reliability_metric = "info",
      M = 3000L,
      seed = 42
    ))
  )
})


# ---------------------------------------------------------------------------
# Metrics "msem"/"bar" are rejected by EQC policy
# ---------------------------------------------------------------------------

test_that("EQC with metric 'msem' errors before calibration", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.75,
      n_items = 20,
      model = "rasch",
      item_source = "parametric",
      reliability_metric = "msem",
      M = 3000L,
      seed = 42
    ),
    "sac_calibrate.*msem"
  )
})

test_that("EQC rejects the documented non-monotone MSEM failure-probe path", {
  failure_probe_beta <- seq(-2, 2, length.out = 20)

  expect_error(
    eqc_calibrate(
      target_rho = 0.75,
      n_items = 20,
      model = "rasch",
      item_source = "custom",
      item_params = list(custom_params = list(beta = failure_probe_beta)),
      reliability_metric = "msem",
      M = 5000L,
      c_bounds = c(0.1, 50),
      seed = 42
    ),
    "non-monotone.*sac_calibrate"
  )
})


# ---------------------------------------------------------------------------
# Synonyms: "tilde" maps to "info"; "bar" is rejected with MSEM policy
# ---------------------------------------------------------------------------

test_that("metric synonym 'tilde' maps to metric='info' internally", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.75,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "tilde",
    M = 3000L,
    seed = 42
  ))
  expect_equal(eqc$metric, "info")
})

test_that("metric synonym 'bar' errors before item generation", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.75,
      n_items = 20,
      model = "rasch",
      item_source = "custom",
      item_params = list(),
      reliability_metric = "bar",
      M = 3000L,
      seed = 42
    ),
    "sac_calibrate.*msem"
  )
})


# ---------------------------------------------------------------------------
# Determinism: same seed -> identical c_star
# ---------------------------------------------------------------------------

test_that("EQC with same seed produces identical c_star and achieved_rho", {
  eqc1 <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))
  eqc2 <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))
  expect_equal(eqc1$c_star, eqc2$c_star)
  expect_equal(eqc1$achieved_rho, eqc2$achieved_rho)
})


# ---------------------------------------------------------------------------
# Return structure: all expected fields present
# ---------------------------------------------------------------------------

test_that("EQC return object has correct class and all expected fields", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 20, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))

  expect_s3_class(eqc, "eqc_result")

  expected_fields <- c(
    "c_star", "target_rho", "achieved_rho", "metric",
    "model", "n_items", "M", "theta_quad", "theta_var",
    "beta_vec", "lambda_base", "lambda_scaled",
    "items_base", "items_calib", "call", "misc"
  )
  for (f in expected_fields) {
    expect_true(f %in% names(eqc), info = paste("Missing field:", f))
  }

  # Dimensional checks
  expect_length(eqc$theta_quad, 3000)
  expect_length(eqc$beta_vec, 20)
  expect_length(eqc$lambda_base, 20)
  expect_length(eqc$lambda_scaled, 20)

  # Numerical sanity
  expect_true(eqc$c_star > 0)
  expect_equal(eqc$target_rho, 0.80)
  expect_equal(eqc$n_items, 20L)
  expect_true(eqc$theta_var > 0)

  # misc sub-structure
  expect_true("root_status" %in% names(eqc$misc))
  expect_true("c_bounds" %in% names(eqc$misc))
  expect_true("rho_bounds" %in% names(eqc$misc))
})


# ---------------------------------------------------------------------------
# Boundary: high target with few items -> upper bound warning
# ---------------------------------------------------------------------------

test_that("EQC: target_rho=0.99 with n_items=5 triggers upper bound warning", {
  eqc <- NULL
  expect_warning(
    eqc <- suppressMessages(eqc_calibrate(
      target_rho = 0.99,
      n_items = 5,
      model = "rasch",
      item_source = "parametric",
      M = 3000L,
      seed = 42
    )),
    "exceeds the maximum achievable"
  )
  expect_equal(eqc$misc$root_status, "above_upper")
  expect_equal(eqc$misc$target_status, "above_upper")
  expect_equal(eqc$c_star, eqc$misc$c_bounds[2])
  expect_null(eqc$misc$uniroot_result)
})


# ---------------------------------------------------------------------------
# Boundary: very low target -> lower bound warning
# ---------------------------------------------------------------------------

test_that("EQC: very low target_rho triggers below-lower warning", {
  eqc <- NULL
  expect_warning(
    eqc <- suppressMessages(eqc_calibrate(
      target_rho = 0.01,
      n_items = 50,
      model = "rasch",
      item_source = "parametric",
      M = 3000L,
      seed = 42,
      c_bounds = c(0.3, 3)
    )),
    "below the minimum achievable"
  )
  expect_equal(eqc$misc$root_status, "below_lower")
  expect_equal(eqc$misc$target_status, "below_lower")
  expect_equal(eqc$c_star, eqc$misc$c_bounds[1])
  expect_null(eqc$misc$uniroot_result)
})

test_that("EQC exact endpoint targets return boundary status without warning", {
  base <- suppressMessages(eqc_calibrate(
    target_rho = 0.80,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    M = 2000L,
    c_bounds = c(0.3, 3),
    seed = 42,
    verbose = FALSE
  ))

  lower_target <- unname(base$misc$rho_bounds_named["lower"])
  upper_target <- unname(base$misc$rho_bounds_named["upper"])

  eqc_lower <- expect_no_warning(suppressMessages(eqc_calibrate(
    target_rho = lower_target,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    M = 2000L,
    c_bounds = c(0.3, 3),
    seed = 42,
    verbose = FALSE
  )))
  expect_equal(eqc_lower$misc$root_status, "boundary_lower")
  expect_equal(eqc_lower$misc$target_status, "boundary_lower")
  expect_equal(eqc_lower$c_star, 0.3)
  expect_equal(eqc_lower$achieved_rho, lower_target, tolerance = 1e-12)
  expect_null(eqc_lower$misc$uniroot_result)

  eqc_upper <- expect_no_warning(suppressMessages(eqc_calibrate(
    target_rho = upper_target,
    n_items = 20,
    model = "rasch",
    item_source = "parametric",
    M = 2000L,
    c_bounds = c(0.3, 3),
    seed = 42,
    verbose = FALSE
  )))
  expect_equal(eqc_upper$misc$root_status, "boundary_upper")
  expect_equal(eqc_upper$misc$target_status, "boundary_upper")
  expect_equal(eqc_upper$c_star, 3)
  expect_equal(eqc_upper$achieved_rho, upper_target, tolerance = 1e-12)
  expect_null(eqc_upper$misc$uniroot_result)
})


# ---------------------------------------------------------------------------
# Root status: normal case should be "uniroot_success"
# ---------------------------------------------------------------------------

test_that("EQC normal case has root_status = 'uniroot_success'", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))
  expect_equal(eqc$misc$root_status, "uniroot_success")
})


# ---------------------------------------------------------------------------
# Full uniroot diagnostics stored in misc$uniroot_result (Issue #37)
# ---------------------------------------------------------------------------

test_that("EQC stores full uniroot diagnostics in misc$uniroot_result", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80, n_items = 25, model = "rasch",
    item_source = "parametric", M = 3000L, seed = 42
  ))
  expect_true(!is.null(eqc$misc$uniroot_result))
  ur <- eqc$misc$uniroot_result
  expect_true("root" %in% names(ur))
  expect_true("f.root" %in% names(ur))
  expect_true("iter" %in% names(ur))
  expect_true("init.it" %in% names(ur))
  expect_true("estim.prec" %in% names(ur))
  expect_equal(ur$root, eqc$c_star)
})

test_that("EQC boundary case has NULL uniroot_result", {
  eqc <- suppressWarnings(suppressMessages(eqc_calibrate(
    target_rho = 0.99,
    n_items = 5,
    model = "rasch",
    item_source = "parametric",
    M = 3000L,
    seed = 42
  )))
  expect_null(eqc$misc$uniroot_result)
})

test_that("EQC finds an interior root when a very wide upper bound saturates", {
  eqc <- suppressMessages(eqc_calibrate(
    target_rho = 0.80,
    n_items = 25,
    model = "rasch",
    item_source = "parametric",
    M = 1000L,
    c_bounds = c(0.3, 1e6),
    seed = 42,
    verbose = FALSE
  ))

  expect_equal(eqc$misc$root_status, "uniroot_success")
  expect_equal(eqc$misc$target_status, "feasible")
  expect_equal(eqc$misc$bracket_method, "interior_max")
  expect_lt(unname(eqc$misc$rho_bounds_named["upper"]), eqc$target_rho)
  expect_gt(eqc$misc$rho_max, eqc$target_rho)
  expect_lt(eqc$c_star, 10)
  expect_equal(eqc$achieved_rho, eqc$target_rho, tolerance = 1e-3)
})


# ---------------------------------------------------------------------------
# Parameter auto-wrapping (Issue #27)
# ---------------------------------------------------------------------------

test_that("EQC auto-wraps shape params in latent_params", {
  expect_message(
    suppressWarnings(eqc_calibrate(
      target_rho = 0.80, n_items = 20, model = "rasch",
      item_source = "parametric",
      latent_shape = "bimodal",
      latent_params = list(delta = 0.9),
      M = 2000L, seed = 42
    )),
    "Auto-wrapping"
  )
})


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

describe("EQC input validation", {

  test_that("target_rho=0 errors", {
    expect_error(
      eqc_calibrate(target_rho = 0, n_items = 25, item_source = "parametric"),
      "target_rho"
    )
  })

  test_that("target_rho=1 errors", {
    expect_error(
      eqc_calibrate(target_rho = 1, n_items = 25, item_source = "parametric"),
      "target_rho"
    )
  })

  test_that("target_rho=1.5 errors", {
    expect_error(
      eqc_calibrate(target_rho = 1.5, n_items = 25, item_source = "parametric"),
      "target_rho"
    )
  })

  test_that("n_items=0 errors", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 0, item_source = "parametric"),
      "n_items"
    )
  })

  test_that("n_items fractional errors", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 25.5,
                    item_source = "parametric"),
      "n_items"
    )
  })

  test_that("M=0 errors", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 25,
                     item_source = "parametric", M = 0),
      "M"
    )
  })

  test_that("M must be sufficient to estimate variance", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 25,
                    item_source = "parametric", M = 1),
      "M.*at least 2"
    )
  })

  test_that("invalid c_bounds errors", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 25,
                     item_source = "parametric", c_bounds = c(3, 1)),
      "c_bounds"
    )
  })

  test_that("invalid reliability_metric errors", {
    expect_error(
      eqc_calibrate(target_rho = 0.80, n_items = 25,
                    item_source = "parametric",
                    reliability_metric = "foobar"),
      "reliability_metric"
    )
  })
})
