# =============================================================================
# test-argument_normalization.R
# =============================================================================
# latent_params/item_params ownership and normalization
# =============================================================================

test_that("sim_item_params defaults to parametric source", {
  items <- sim_item_params(n_items = 8L, model = "rasch", seed = 1)

  expect_s3_class(items, "item_params")
  expect_equal(items$source, "parametric")
  expect_equal(nrow(items$data), 8L)
})

test_that("direct sim_item_params still supports multiple forms", {
  items <- sim_item_params(
    n_items = 5L,
    model = "rasch",
    source = "parametric",
    n_forms = 3L,
    seed = 1
  )

  expect_equal(items$n_forms, 3L)
  expect_equal(nrow(items$data), 15L)
  expect_equal(length(unique(items$data$form_id)), 3L)
})

test_that("latent shape parameters are auto-wrapped consistently", {
  expect_message(
    suppressWarnings(eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      latent_shape = "bimodal",
      latent_params = list(delta = 0.9),
      M = 500L,
      seed = 1,
      verbose = FALSE
    )),
    "Auto-wrapping shape parameter"
  )

  expect_message(
    suppressWarnings(check_feasibility(
      n_items = 8L,
      latent_shape = "bimodal",
      latent_params = list(delta = 0.9),
      M = 300L,
      seed = 1,
      verbose = FALSE
    )),
    "Auto-wrapping shape parameter"
  )
})

test_that("latent_params cannot override top-level owner arguments", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      latent_params = list(n = 999L),
      M = 100L,
      seed = 1
    ),
    "latent_params.*reserved.*`n`"
  )

  expect_error(
    sac_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      latent_params = list(shape = "uniform"),
      M_per_iter = 50L,
      M_pre = 100L,
      n_iter = 5L,
      burn_in = 2L,
      seed = 1
    ),
    "latent_params.*reserved.*`shape`"
  )

  result <- suppressWarnings(suppressMessages(eqc_calibrate(
    target_rho = 0.70,
    n_items = 4L,
    item_source = "parametric",
    M = 100L,
    seed = 1
  )))

  expect_error(
    simulate_response_data(
      result,
      n_persons = 10L,
      latent_params = list(seed = 2L),
      seed = 1
    ),
    "latent_params.*reserved.*`seed`"
  )
})

test_that("shape params cannot be specified twice", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      latent_shape = "bimodal",
      latent_params = list(delta = 0.9, shape_params = list(delta = 0.8)),
      M = 100L,
      seed = 1
    ),
    "shape-specific argument"
  )
})

test_that("item_params cannot override top-level owner arguments", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = list(n_items = 99L),
      M = 100L,
      seed = 1
    ),
    "item_params.*reserved.*`n_items`"
  )

  expect_error(
    sac_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = list(model = "2pl"),
      M_per_iter = 50L,
      M_pre = 100L,
      n_iter = 5L,
      burn_in = 2L,
      seed = 1
    ),
    "item_params.*reserved.*`model`"
  )

  expect_error(
    check_feasibility(
      n_items = 8L,
      item_source = "parametric",
      item_params = list(source = "irw"),
      M = 100L,
      seed = 1,
      verbose = FALSE
    ),
    "item_params.*reserved.*`source`"
  )

  expect_error(
    rho_curve(
      c_values = c(0.5, 1, 1.5),
      n_items = 8L,
      item_source = "parametric",
      item_params = list(scale = 2),
      M = 100L,
      seed = 1,
      plot = FALSE
    ),
    "item_params.*reserved.*`scale`"
  )
})

test_that("n_forms is direct-only and rejected in calibrator item_params", {
  eqc <- suppressWarnings(suppressMessages(eqc_calibrate(
    target_rho = 0.70,
    n_items = 8L,
    item_source = "parametric",
    item_params = list(n_forms = 1L),
    M = 100L,
    seed = 1,
    verbose = FALSE
  )))

  expect_s3_class(eqc, "eqc_result")
  expect_equal(eqc$items_base$n_forms, 1L)

  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = list(n_forms = 2L),
      M = 100L,
      seed = 1
    ),
    "n_forms.*direct calls to `sim_item_params\\(\\)`"
  )

  expect_error(
    sac_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = list(n_forms = 2L),
      M_per_iter = 50L,
      M_pre = 100L,
      n_iter = 5L,
      burn_in = 2L,
      seed = 1
    ),
    "n_forms.*direct calls to `sim_item_params\\(\\)`"
  )

  expect_error(
    check_feasibility(
      n_items = 8L,
      item_source = "parametric",
      item_params = list(n_forms = 2L),
      M = 100L,
      seed = 1,
      verbose = FALSE
    ),
    "n_forms.*direct calls to `sim_item_params\\(\\)`"
  )

  expect_error(
    rho_curve(
      c_values = c(0.5, 1, 1.5),
      n_items = 8L,
      item_source = "parametric",
      item_params = list(n_forms = 2L),
      M = 100L,
      seed = 1,
      plot = FALSE
    ),
    "n_forms.*direct calls to `sim_item_params\\(\\)`"
  )
})

test_that("item_params cannot carry a nested seed", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = list(seed = 2L),
      M = 100L,
      seed = 1
    ),
    "item_params.*reserved.*`seed`"
  )
})

test_that("generator-specific item parameters still pass through normalization", {
  eqc <- suppressWarnings(suppressMessages(eqc_calibrate(
    target_rho = 0.70,
    n_items = 8L,
    model = "2pl",
    item_source = "parametric",
    item_params = list(
      method = "independent",
      difficulty_params = list(mu = 0.5, sigma = 0.2),
      discrimination_params = list(mu_log = 0, sigma_log = 0.1),
      center_difficulties = FALSE
    ),
    M = 200L,
    seed = 1,
    verbose = FALSE
  )))

  expect_s3_class(eqc, "eqc_result")
  expect_equal(eqc$items_base$method, "independent")
  expect_false(eqc$items_base$centered)
})

test_that("item_params and latent_params must be lists", {
  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      latent_params = "delta",
      M = 100L,
      seed = 1
    ),
    "`latent_params` must be a list"
  )

  expect_error(
    eqc_calibrate(
      target_rho = 0.70,
      n_items = 8L,
      item_source = "parametric",
      item_params = "difficulty",
      M = 100L,
      seed = 1
    ),
    "`item_params` must be a list"
  )
})
