# =============================================================================
# test-edge_cases.R
# =============================================================================
# Focused edge-case coverage for scalar validation, optional dependency guards,
# nonfinite calibrated designs, and low-boundary simulation sizes.
# =============================================================================

test_that("direct sim_item_params rejects fractional and nonfinite design scalars", {
  expect_error(
    sim_item_params(n_items = 5.5, model = "rasch", source = "parametric"),
    "n_items"
  )

  expect_error(
    sim_item_params(
      n_items = 5L,
      model = "rasch",
      source = "parametric",
      n_forms = 1.5
    ),
    "n_forms"
  )

  expect_error(
    sim_item_params(
      n_items = 5L,
      model = "rasch",
      source = "parametric",
      n_forms = NA_real_
    ),
    "n_forms"
  )

  expect_error(
    sim_item_params(
      n_items = 5L,
      model = "rasch",
      source = "parametric",
      scale = Inf
    ),
    "scale"
  )
})

test_that("rho_curve rejects invalid grid and Monte Carlo boundary inputs", {
  expect_error(
    rho_curve(
      c_values = c(0.5, 1),
      n_items = 5.5,
      item_source = "parametric",
      M = 50L,
      plot = FALSE
    ),
    "n_items"
  )

  expect_error(
    rho_curve(
      c_values = c(0.5, Inf),
      n_items = 5L,
      item_source = "parametric",
      M = 50L,
      plot = FALSE
    ),
    "c_values"
  )

  expect_error(
    rho_curve(
      c_values = c(0.5, 1),
      n_items = 5L,
      item_source = "parametric",
      M = 1L,
      plot = FALSE
    ),
    "at least 2"
  )

  expect_error(
    rho_curve(
      c_values = c(0.5, 1),
      n_items = 5L,
      item_source = "parametric",
      M = 50.5,
      plot = FALSE
    ),
    "M"
  )
})

test_that("rho_curve accepts the documented finite lower Monte Carlo boundary", {
  curve <- rho_curve(
    c_values = c(0.5, 1),
    n_items = 4L,
    item_source = "parametric",
    M = 2L,
    seed = 42,
    plot = FALSE
  )

  expect_s3_class(curve, "rho_curve")
  expect_equal(nrow(curve), 2L)
})

test_that("IRW optional dependency absence fails with an actionable error", {
  if (requireNamespace("irw", quietly = TRUE)) {
    skip("irw is installed; absence path cannot be tested in this library.")
  }

  expect_error(
    sim_item_params(n_items = 5L, model = "rasch", source = "irw"),
    "Package 'irw' is required"
  )
})

test_that("latent_params and item_params must not route unnamed values positionally", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 5L,
      item_source = "parametric",
      latent_params = list(2),
      M = 50L,
      seed = 1,
      verbose = FALSE
    ),
    "latent_params.*named"
  )

  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 5L,
      item_source = "parametric",
      item_params = list("independent"),
      M = 50L,
      seed = 1,
      verbose = FALSE
    ),
    "item_params.*named"
  )
})

test_that("low-level reliability utilities reject nonfinite core inputs", {
  theta <- c(-1, 0, 1)
  beta <- c(-0.5, 0.5)
  lambda <- c(1, 1.2)

  expect_error(
    compute_rho_tilde(Inf, theta, beta, lambda),
    "c"
  )

  expect_error(
    compute_rho_bar(1, c(NA_real_, 1), beta, lambda),
    "theta_vec"
  )

  expect_error(
    compute_rho_both(1, theta, beta, c(1, Inf)),
    "lambda_base"
  )
})

test_that("simulate_response_data rejects nonfinite calibrated designs", {
  result <- suppressMessages(eqc_calibrate(
    target_rho = 0.50,
    n_items = 4L,
    item_source = "parametric",
    M = 100L,
    seed = 11,
    verbose = FALSE
  ))

  nonfinite_lambda <- result
  nonfinite_lambda$lambda_scaled[1] <- Inf
  expect_error(
    simulate_response_data(nonfinite_lambda, n_persons = 10L, seed = 1),
    "lambda_scaled.*finite"
  )

  nonfinite_item <- result
  nonfinite_item$items_calib$data$lambda[1] <- NaN
  expect_error(
    simulate_response_data(nonfinite_item, n_persons = 10L, seed = 1),
    "items_calib\\$data\\$lambda.*finite"
  )
})

test_that("check_feasibility classifies exact upper boundary targets", {
  base <- check_feasibility(
    n_items = 8L,
    item_source = "parametric",
    seed = 7,
    M = 300L,
    verbose = FALSE
  )
  upper <- check_feasibility(
    n_items = 8L,
    item_source = "parametric",
    target_rho = base$rho_range_info["upper"],
    seed = 7,
    M = 300L,
    verbose = FALSE
  )

  expect_equal(upper$target_status_info, "boundary_upper")
})
