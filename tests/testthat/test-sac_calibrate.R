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

test_that("SAC metric aliases are stored canonically", {
  sac_bar <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.70,
    n_items = 6,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "bar",
    n_iter = 8L,
    burn_in = 4L,
    M_per_iter = 60L,
    M_pre = 120L,
    seed = 3
  )))

  sac_tilde <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.70,
    n_items = 6,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "tilde",
    n_iter = 8L,
    burn_in = 4L,
    M_per_iter = 60L,
    M_pre = 120L,
    seed = 3
  )))

  expect_equal(sac_bar$metric, "msem")
  expect_equal(sac_tilde$metric, "info")
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
    "raw_trajectory", "step_size_trajectory", "gradient_trajectory",
    "projected", "projection_side", "projection_count", "projection_rate",
    "M_final", "metric", "model", "n_items", "n_iter", "burn_in",
    "M_per_iter", "M_pre", "step_params", "c_bounds",
    "c_init", "init_method", "item_design", "calibration_status", "convergence",
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
  expect_true("split_mean_diff" %in% names(sac$convergence))
  expect_true("split_mean_rel_diff" %in% names(sac$convergence))
  expect_true("split_mean_scale" %in% names(sac$convergence))
  expect_true("split_mean_tolerance" %in% names(sac$convergence))
  expect_true("split_mean_rel_tolerance" %in% names(sac$convergence))
  expect_true("sd_post_burn" %in% names(sac$convergence))
  expect_true("n_post_burn" %in% names(sac$convergence))
  expect_true("insufficient_post_burn" %in% names(sac$convergence))
  expect_true("trajectory_converged" %in% names(sac$convergence))
  expect_true("achieved_gap" %in% names(sac$convergence))
  expect_true("achieved_gap_abs" %in% names(sac$convergence))
  expect_true("achieved_gap_tolerance" %in% names(sac$convergence))
  expect_true("large_achieved_gap" %in% names(sac$convergence))
  expect_true("projection_count" %in% names(sac$convergence))
  expect_true("projection_rate" %in% names(sac$convergence))
  expect_true("n_lower_hits" %in% names(sac$convergence))
  expect_true("n_upper_hits" %in% names(sac$convergence))
  expect_true("status_flags" %in% names(sac$convergence))
  expect_true("hit_lower_bound" %in% names(sac$convergence))
  expect_true("hit_upper_bound" %in% names(sac$convergence))
  expect_true("final_pre_update_gradient" %in% names(sac$convergence))
  expect_true("final_iter_gradient" %in% names(sac$convergence))
  expect_true("gradient_at_c_star" %in% names(sac$convergence))
  expect_true("post_calibration_gradient" %in% names(sac$convergence))
  expect_true("status" %in% names(sac$convergence))
  expect_equal(sac$calibration_status, sac$convergence$status)
  expect_true(sac$calibration_status %in% c(
    "ok", "not_converged", "hit_lower_bound", "hit_upper_bound", "hit_both_bounds"
  ))

  # Numerical checks
  expect_true(sac$c_star > 0)
  expect_true(sac$theta_var > 0)
  expect_length(sac$beta_vec, 15)
  expect_length(sac$lambda_base, 15)
  expect_length(sac$lambda_scaled, 15)
  expect_length(sac$raw_trajectory, sac$n_iter)
  expect_length(sac$step_size_trajectory, sac$n_iter)
  expect_length(sac$gradient_trajectory, sac$n_iter)
  expect_length(sac$projected, sac$n_iter)
  expect_length(sac$projection_side, sac$n_iter)
  expect_equal(sac$n_items, 15L)
  expect_equal(sac$c_final, utils::tail(sac$trajectory, 1))
  expect_true(all(sac$trajectory >= sac$c_bounds[1]))
  expect_true(all(sac$trajectory <= sac$c_bounds[2]))
  expect_equal(length(sac$theta_quad), sac$M_final)
  expect_equal(sac$M_final, max(sac$M_per_iter * 5, 5000))
  expect_equal(sac$projection_count, sum(sac$projected))
  expect_equal(sac$projection_rate, mean(sac$projected))
  expect_equal(sac$convergence$projection_count, sac$projection_count)
  expect_equal(sac$convergence$projection_rate, sac$projection_rate)
  expect_true(all(sac$step_size_trajectory > 0))
  expect_true(all(diff(sac$step_size_trajectory) < 0))
  expect_equal(sac$gradient_trajectory,
               sac$rho_trajectory - sac$target_rho)
  expect_true(sac$item_design %in% c(
    "post_calibration_draw", "fixed_iteration_items"
  ))
  expect_equal(sac$convergence$final_gradient,
               utils::tail(sac$rho_trajectory, 1) - sac$target_rho)
  expect_equal(sac$convergence$final_pre_update_gradient,
               sac$convergence$final_gradient)
  expect_equal(sac$convergence$final_iter_gradient,
               sac$convergence$final_gradient)
  expect_equal(sac$convergence$gradient_at_c_star,
               sac$achieved_rho - sac$target_rho)
  expect_equal(sac$convergence$post_calibration_gradient,
               sac$convergence$gradient_at_c_star)
  expect_equal(sac$convergence$achieved_gap,
               sac$achieved_rho - sac$target_rho)
  expect_equal(sac$convergence$achieved_gap_abs,
               abs(sac$achieved_rho - sac$target_rho))
  expect_equal(sac$convergence$large_achieved_gap,
               sac$convergence$achieved_gap_abs >
                 sac$convergence$achieved_gap_tolerance)
  expect_equal(sac$convergence$split_mean_scale,
               max(abs(sac$c_star), .Machine$double.eps))
  expect_equal(sac$convergence$split_mean_tolerance,
               sac$convergence$split_mean_rel_tolerance *
                 sac$convergence$split_mean_scale)
  expect_equal(sac$convergence$split_mean_rel_tolerance, 0.05)
  expect_equal(sac$convergence$split_mean_rel_diff,
               sac$convergence$split_mean_diff /
                 sac$convergence$split_mean_scale)
  expect_equal(sac$convergence$trajectory_converged,
               sac$convergence$split_mean_rel_diff <
                 sac$convergence$split_mean_rel_tolerance)
  expect_equal(sac$convergence$n_post_burn, sac$n_iter - sac$burn_in)
  expect_false(sac$convergence$insufficient_post_burn)
})

test_that("SAC step size defaults, overrides, and validation are explicit", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 10,
    model = "rasch",
    item_source = "parametric",
    n_iter = 20L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))
  expect_equal(sac$step_params$a, 1.0)
  expect_equal(sac$step_params$A, 50)
  expect_equal(sac$step_params$gamma, 0.67)

  sac_override <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 10,
    model = "rasch",
    item_source = "parametric",
    step_params = list(a = 0.5, A = 20, gamma = 0.8),
    n_iter = 20L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))
  expect_equal(sac_override$step_params$a, 0.5)
  expect_equal(sac_override$step_params$A, 20)
  expect_equal(sac_override$step_params$gamma, 0.8)

  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = "bad"),
    "step_params"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(a = c(1, 2))),
    "step_params\\$a"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(a = NA_real_)),
    "step_params\\$a"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(a = Inf)),
    "step_params\\$a"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(a = -1)),
    "step_params\\$a"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(A = -1)),
    "step_params\\$A"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(gamma = 0.5)),
    "Robbins-Monro"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(gamma = 1.1)),
    "Robbins-Monro"
  )
  expect_error(
    sac_calibrate(0.75, 10, item_source = "parametric",
                  step_params = list(foo = 1)),
    "Unknown `step_params`"
  )
})

test_that("SAC projection diagnostics record initial clipping and bound hits", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 10,
    model = "rasch",
    item_source = "parametric",
    c_init = 100,
    c_bounds = c(0.2, 0.3),
    step_params = list(a = 2),
    n_iter = 20L,
    burn_in = 10L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))

  expect_equal(sac$c_init, 0.3)
  expect_true(sac$convergence$c_init_projected)
  expect_equal(sac$convergence$c_init_raw, 100)
  expect_true(sac$convergence$hit_upper_bound)
  expect_true(sac$projection_count > 0)
  expect_equal(sac$projection_count, sac$convergence$projection_count)
  expect_equal(sac$projection_rate, sac$convergence$projection_rate)
  expect_true("projection_applied" %in% sac$convergence$status_flags)
  expect_true(all(sac$raw_trajectory != sac$trajectory | !sac$projected))
  expect_true(all(sac$projection_side %in% c("none", "lower", "upper")))
})

test_that("SAC emits non-convergence warning when achieved gap is large", {
  sac <- NULL
  expect_warning(
    sac <- suppressMessages(sac_calibrate(
      target_rho = 0.90,
      n_items = 8,
      model = "rasch",
      item_source = "parametric",
      c_init = 0.01,
      step_params = list(a = 1e-8),
      n_iter = 20L,
      burn_in = 10L,
      M_per_iter = 80L,
      M_pre = 800L,
      seed = 42
    )),
    "SAC may not have fully converged"
  )

  expect_true(sac$convergence$large_achieved_gap)
  expect_true("large_achieved_gap" %in% sac$convergence$status_flags)
  expect_equal(sac$projection_count, 0)
})

test_that("SAC emits upper projection warning when upper bound is active", {
  sac <- NULL
  expect_warning(
    sac <- suppressMessages(sac_calibrate(
      target_rho = 0.20,
      n_items = 10,
      model = "rasch",
      item_source = "parametric",
      c_init = 100,
      c_bounds = c(0.2, 0.3),
      step_params = list(a = 2),
      n_iter = 20L,
      burn_in = 10L,
      M_per_iter = 100L,
      M_pre = 1000L,
      seed = 42
    )),
    "SAC projection hit the upper bound"
  )

  expect_false(sac$convergence$hit_lower_bound)
  expect_true(sac$convergence$hit_upper_bound)
  expect_equal(sac$projection_count, 20L)
  expect_true("projection_applied" %in% sac$convergence$status_flags)
})

test_that("SAC emits lower projection warning when lower bound is active", {
  sac <- NULL
  expect_warning(
    sac <- suppressMessages(sac_calibrate(
      target_rho = 0.05,
      n_items = 10,
      model = "rasch",
      item_source = "parametric",
      c_init = 0.001,
      c_bounds = c(0.2, 0.3),
      step_params = list(a = 2),
      n_iter = 20L,
      burn_in = 10L,
      M_per_iter = 100L,
      M_pre = 1000L,
      seed = 42
    )),
    "SAC projection hit the lower bound"
  )

  expect_true(sac$convergence$hit_lower_bound)
  expect_false(sac$convergence$hit_upper_bound)
  expect_equal(sac$projection_count, 20L)
  expect_true("projection_applied" %in% sac$convergence$status_flags)
})

test_that("SAC flags insufficient post-burn-in diagnostics", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 8,
    model = "rasch",
    item_source = "parametric",
    n_iter = 3L,
    burn_in = 2L,
    M_per_iter = 80L,
    M_pre = 800L,
    seed = 42
  )))

  expect_equal(sac$convergence$n_post_burn, 1L)
  expect_true(sac$convergence$insufficient_post_burn)
  expect_true("insufficient_post_burn" %in% sac$convergence$status_flags)
  expect_equal(sac$calibration_status, "not_converged")
  expect_false(sac$convergence$converged)
  expect_equal(sac$convergence$sd_post_burn, 0)
  expect_true(is.na(sac$convergence$split_mean_diff))
  expect_true(is.na(sac$convergence$split_mean_rel_diff))
  expect_true(is.finite(sac$convergence$split_mean_scale))
  expect_true(is.finite(sac$convergence$split_mean_tolerance))
  expect_equal(sac$convergence$split_mean_rel_tolerance, 0.05)
})

test_that("SAC does not report ok when achieved reliability misses target", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.90,
    n_items = 8,
    model = "rasch",
    item_source = "parametric",
    c_init = 0.01,
    step_params = list(a = 1e-8),
    n_iter = 20L,
    burn_in = 10L,
    M_per_iter = 80L,
    M_pre = 800L,
    seed = 42
  )))

  expect_true(sac$convergence$large_achieved_gap)
  expect_true("large_achieved_gap" %in% sac$convergence$status_flags)
  expect_equal(sac$calibration_status, "not_converged")
  expect_false(sac$convergence$converged)
  expect_true(abs(sac$achieved_rho - sac$target_rho) >
                sac$convergence$achieved_gap_tolerance)
})

test_that("SAC calibrated item object uses c_star on the post-calibration draw", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 12,
    model = "2pl",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 1500L,
    seed = 42
  )))

  expect_equal(sac$items_base$scale, 1)
  expect_equal(sac$items_calib$scale, sac$c_star)
  expect_equal(sac$lambda_scaled, sac$lambda_base * sac$c_star,
               tolerance = 1e-10)
  expect_equal(sac$items_calib$data$lambda, sac$lambda_scaled,
               tolerance = 1e-10)
  expect_equal(sac$items_calib$data$lambda_unscaled, sac$lambda_base,
               tolerance = 1e-10)
  expect_equal(sac$items_calib$achieved$overall$lambda_mean,
               mean(sac$lambda_scaled), tolerance = 1e-10)
  expect_equal(sac$items_calib$achieved$overall$lambda_sd,
               stats::sd(sac$lambda_scaled), tolerance = 1e-10)
})

test_that("SAC with resample_items=FALSE stores the fixed iteration item form", {
  beta_call_count <- 0L
  beta_fun <- function(n) {
    beta_call_count <<- beta_call_count + 1L
    seq_len(n) + 100 * beta_call_count
  }

  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 6,
    model = "rasch",
    item_source = "custom",
    item_params = list(
      custom_params = list(beta = beta_fun),
      center_difficulties = FALSE
    ),
    resample_items = FALSE,
    n_iter = 30L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))

  expect_equal(beta_call_count, 1L)
  expect_equal(sac$item_design, "fixed_iteration_items")
  expect_equal(sac$beta_vec, seq_len(6) + 100)
  expect_equal(sac$items_base$data$beta, sac$beta_vec)
  expect_equal(sac$items_calib$data$beta, sac$beta_vec)
  expect_equal(sac$items_calib$data$lambda, sac$lambda_scaled)
})

test_that("SAC predict at c_star reproduces achieved_rho from stored final draw", {
  sac <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.75,
    n_items = 12,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 1500L,
    seed = 42
  )))

  pred_at_cstar <- predict(sac, newdata = sac$c_star)
  expect_equal(unname(pred_at_cstar), sac$achieved_rho, tolerance = 1e-12)
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

  test_that("M_pre must be sufficient to estimate variance", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                    item_source = "parametric", M_pre = 1),
      "M_pre.*at least 2"
    )
  })

  test_that("integer-like arguments reject fractional values", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20.5,
                    item_source = "parametric"),
      "n_items"
    )
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                    item_source = "parametric", n_iter = 10.5),
      "n_iter"
    )
  })

  test_that("invalid c_init type errors", {
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                     item_source = "parametric", c_init = "bad"),
      "c_init"
    )
  })

  test_that("malformed EQC warm starts error before extraction", {
    bad_eqc <- structure(list(c_star = 1), class = "eqc_result")
    expect_error(
      sac_calibrate(target_rho = 0.75, n_items = 20,
                    item_source = "parametric", c_init = bad_eqc),
      "c_init.*missing required field"
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
