test_that("EQC preserves the v0.2 formal prefix and appends topology controls", {
  expected_prefix <- c(
    "target_rho", "n_items", "model", "latent_shape", "item_source",
    "latent_params", "item_params", "reliability_metric", "M", "c_bounds",
    "tol", "seed", "verbose"
  )
  observed <- names(formals(eqc_calibrate))

  expect_identical(observed[seq_along(expected_prefix)], expected_prefix)
  expect_identical(
    tail(observed, 4L),
    c("root_policy", "root_controls", "allow_tangent", "allow_plateau")
  )
})

test_that("EQC near-maximum diagnostic reports the detected global maximum", {
  messages <- capture_messages(
    eqc_calibrate(
      target_rho = 0.98,
      n_items = 1,
      model = "rasch",
      latent_shape = "normal",
      item_source = "custom",
      item_params = list(
        custom_params = list(beta = 0),
        center_difficulties = FALSE
      ),
      M = 100,
      c_bounds = c(0.01, 1e6),
      seed = 1,
      verbose = TRUE
    )
  )

  expect_true(any(grepl(
    "near the achievable maximum \\(1\\.000\\)", messages
  )))
  expect_false(any(grepl(
    "achievable maximum \\(0\\.000\\)", messages
  )))
})

test_that("EQC 3PL with zero guessing is the exact 2PL special case", {
  common <- list(
    target_rho = 0.70,
    n_items = 12L,
    item_source = "parametric",
    M = 1000L,
    c_bounds = c(0.1, 5),
    seed = 5101
  )
  twopl <- suppressMessages(do.call(
    eqc_calibrate,
    c(common, list(model = "2pl"))
  ))
  threepl_zero <- suppressMessages(do.call(
    eqc_calibrate,
    c(common, list(
      model = "3pl",
      item_params = list(
        guessing_params = list(distribution = "fixed", value = 0)
      )
    ))
  ))

  expect_equal(threepl_zero$beta_vec, twopl$beta_vec, tolerance = 0)
  expect_equal(threepl_zero$lambda_base, twopl$lambda_base, tolerance = 0)
  expect_equal(threepl_zero$guessing_vec, rep(0, 12L), tolerance = 0)
  expect_equal(threepl_zero$c_star, twopl$c_star, tolerance = 1e-10)
  expect_equal(threepl_zero$achieved_rho, twopl$achieved_rho, tolerance = 1e-12)
})

test_that("EQC calibrates a 3PL fixed form and records additive provenance", {
  out <- suppressMessages(eqc_calibrate(
    target_rho = 0.70,
    n_items = 12L,
    model = "3pl",
    item_source = "parametric",
    item_params = list(
      guessing_params = list(distribution = "fixed", value = 0.20)
    ),
    M = 1500L,
    c_bounds = c(0.1, 5),
    seed = 5102
  ))

  expect_s3_class(out, "eqc_result")
  expect_lt(abs(out$achieved_rho - out$target_rho), 5e-4)
  expect_identical(out$schema_version, 1L)
  expect_identical(out$item_scope, "fixed_form")
  expect_equal(out$guessing_vec, rep(0.20, 12L), tolerance = 0)
  expect_equal(out$items_base$data$guessing, out$guessing_vec, tolerance = 0)
  expect_equal(out$items_calib$data$guessing, out$guessing_vec, tolerance = 0)
  expect_equal(out$items_calib$data$beta, out$items_base$data$beta, tolerance = 0)
  expect_equal(out$estimand_signature$metric, "info")
  expect_equal(out$estimand_signature$approximation_method, "fixed_mc")
  expect_equal(out$estimand_signature$item_measure, "fixed_form")
  expect_equal(out$design_signature$guessing, out$guessing_vec, tolerance = 0)
  expect_identical(out$misc$bracket_method, "adaptive_log_grid")
  expect_identical(out$misc$root_policy, "lowest_increasing")
  expect_true(out$misc$root_count >= 1L)
})

test_that("EQC recovers a known 3PL scale on the selected branch", {
  calibration_seed <- 5191L
  pilot <- suppressMessages(eqc_calibrate(
    target_rho = 0.65,
    n_items = 15L,
    model = "3pl",
    M = 1200L,
    c_bounds = c(0.2, 3),
    seed = calibration_seed,
    verbose = FALSE
  ))
  known_scale <- 0.8
  known_target <- compute_rho_tilde(
    known_scale,
    pilot$theta_quad,
    pilot$beta_vec,
    pilot$lambda_base,
    theta_var = pilot$theta_var,
    guessing = pilot$guessing_vec
  )

  recovered <- suppressMessages(eqc_calibrate(
    target_rho = known_target,
    n_items = 15L,
    model = "3pl",
    M = 1200L,
    c_bounds = c(0.2, 3),
    seed = calibration_seed,
    verbose = FALSE
  ))

  expect_equal(recovered$c_star, known_scale, tolerance = 2e-5)
  expect_equal(recovered$achieved_rho, known_target, tolerance = 1e-7)
  expect_identical(recovered$misc$selected_root_direction, "increasing")
})

test_that("EQC coef and predict expose 3PL guessing without changing 2PL columns", {
  threepl <- suppressMessages(eqc_calibrate(
    target_rho = 0.65,
    n_items = 8L,
    model = "3pl",
    M = 700L,
    c_bounds = c(0.1, 5),
    seed = 5103
  ))
  twopl <- suppressMessages(eqc_calibrate(
    target_rho = 0.65,
    n_items = 8L,
    model = "2pl",
    M = 700L,
    c_bounds = c(0.1, 5),
    seed = 5103
  ))

  expect_identical(
    names(coef(twopl)),
    c("item_id", "beta", "lambda_base", "lambda_scaled", "c_star")
  )
  expect_identical(
    names(coef(threepl)),
    c(
      "item_id", "beta", "lambda_base", "lambda_scaled", "c_star",
      "guessing"
    )
  )
  expect_equal(coef(threepl)$guessing, threepl$guessing_vec, tolerance = 0)

  scales <- c(0.5, threepl$c_star, 2)
  expected <- vapply(scales, function(scale) {
    compute_rho_tilde(
      scale,
      threepl$theta_quad,
      threepl$beta_vec,
      threepl$lambda_base,
      theta_var = threepl$theta_var,
      guessing = threepl$guessing_vec
    )
  }, numeric(1))
  expect_equal(unname(predict(threepl, newdata = scales)), expected, tolerance = 1e-12)
})

test_that("EQC reports canonical infeasibility and returns the best scanned scale", {
  out <- NULL
  expect_warning(
    out <- suppressMessages(eqc_calibrate(
      target_rho = 0.995,
      n_items = 3L,
      model = "3pl",
      M = 700L,
      c_bounds = c(0.1, 3),
      seed = 5104
    )),
    "maximum achievable|outside|infeasible"
  )

  expect_identical(out$misc$target_status, "infeasible_above_range")
  expect_identical(out$misc$calibration_status, "infeasible")
  expect_true(out$c_star >= 0.1 && out$c_star <= 3)
  expect_equal(out$achieved_rho, out$misc$rho_max, tolerance = 1e-8)
  expect_null(out$misc$uniroot_result)
})

test_that("EQC validates and records an explicit root policy", {
  out <- suppressMessages(eqc_calibrate(
    target_rho = 0.65,
    n_items = 10L,
    model = "3pl",
    M = 700L,
    c_bounds = c(0.1, 5),
    root_policy = "lowest_increasing",
    seed = 5105
  ))
  expect_identical(out$misc$root_policy, "lowest_increasing")
  expect_true(out$misc$selected_root_id >= 1L)

  nearest <- suppressMessages(eqc_calibrate(
    target_rho = 0.65,
    n_items = 10L,
    model = "3pl",
    M = 300L,
    c_bounds = c(0.1, 5),
    root_policy = "nearest_increasing",
    seed = 5105
  ))
  expect_identical(nearest$misc$root_policy, "nearest_increasing")
  expect_identical(nearest$misc$root_reference, 1)

  expect_error(
    eqc_calibrate(
      target_rho = 0.65,
      n_items = 10L,
      model = "3pl",
      M = 100L,
      root_policy = "not_a_policy",
      seed = 5105
    ),
    "root_policy|policy"
  )
  expect_error(
    eqc_calibrate(
      target_rho = 0.65,
      n_items = 10L,
      model = "3pl",
      M = 100L,
      allow_tangent = NA,
      seed = 5105
    ),
    "allow_tangent"
  )
})
