# =============================================================================
# test-enhancements.R
# =============================================================================
# Tests for Phase 6 Enhancement Functions
#
# Contents:
#   - compute_rho_both()
#   - check_feasibility()
#   - rho_curve()
#   - coef.eqc_result / predict.eqc_result
#   - coef.sac_result / predict.sac_result
# =============================================================================


# =============================================================================
# Helper: shared fixtures
# =============================================================================
set.seed(42)
theta_fix <- rnorm(2000)
beta_fix  <- rnorm(20)
lambda_fix <- rep(1, 20)
theta_var_fix <- var(theta_fix)


# =============================================================================
# Tests for compute_rho_both()
# =============================================================================

test_that("compute_rho_both returns list with correct names", {
  result <- compute_rho_both(1, theta_fix, beta_fix, lambda_fix)
  expect_type(result, "list")
  expect_named(result, c("rho_tilde", "rho_bar"))
})

test_that("compute_rho_both matches individual functions", {
  c_val <- 1.5
  both <- compute_rho_both(c_val, theta_fix, beta_fix, lambda_fix,
                            theta_var = theta_var_fix)
  tilde <- compute_rho_tilde(c_val, theta_fix, beta_fix, lambda_fix,
                              theta_var = theta_var_fix)
  bar   <- compute_rho_bar(c_val, theta_fix, beta_fix, lambda_fix,
                            theta_var = theta_var_fix)

  expect_equal(both$rho_tilde, tilde, tolerance = 1e-12)
  expect_equal(both$rho_bar, bar, tolerance = 1e-12)
})

test_that("compute_rho_both satisfies Jensen's inequality", {
  for (c_val in c(0.5, 1.0, 2.0, 3.0)) {
    both <- compute_rho_both(c_val, theta_fix, beta_fix, lambda_fix)
    expect_gte(both$rho_tilde, both$rho_bar)
  }
})

test_that("compute_rho_both returns values in (0, 1)", {
  both <- compute_rho_both(1, theta_fix, beta_fix, lambda_fix)
  expect_gt(both$rho_tilde, 0)
  expect_lt(both$rho_tilde, 1)
  expect_gt(both$rho_bar, 0)
  expect_lt(both$rho_bar, 1)
})

test_that("compute_rho_both returns zeros for c <= 0", {
  both <- compute_rho_both(0, theta_fix, beta_fix, lambda_fix)
  expect_equal(both$rho_tilde, 0)
  expect_equal(both$rho_bar, 0)

  both_neg <- compute_rho_both(-1, theta_fix, beta_fix, lambda_fix)
  expect_equal(both_neg$rho_tilde, 0)
  expect_equal(both_neg$rho_bar, 0)
})

test_that("compute_rho_both validates inputs", {
  expect_error(compute_rho_both(1, theta_fix, beta_fix[1:5], lambda_fix),
               "beta_vec and lambda_base must have the same length")
  expect_error(compute_rho_both(1, 5, beta_fix, lambda_fix),
               "theta_vec must have length >= 2")
})

test_that("compute_rho_both works with theta_var argument", {
  both1 <- compute_rho_both(1, theta_fix, beta_fix, lambda_fix)
  both2 <- compute_rho_both(1, theta_fix, beta_fix, lambda_fix,
                             theta_var = theta_var_fix)
  # Should be very close (same theta sample used for both)
  expect_equal(both1$rho_tilde, both2$rho_tilde, tolerance = 1e-10)
  expect_equal(both1$rho_bar, both2$rho_bar, tolerance = 1e-10)
})


# =============================================================================
# Tests for check_feasibility()
# =============================================================================

test_that("check_feasibility returns correct class and structure", {
  feas <- check_feasibility(n_items = 20, model = "rasch", seed = 42,
                             M = 2000, verbose = FALSE)
  expect_s3_class(feas, "feasibility_check")
  expect_true(is.list(feas))
  expect_named(feas, c("rho_range_info", "rho_range_msem", "n_items", "model",
                        "latent_shape", "c_bounds", "M", "theta_var"),
               ignore.order = TRUE)
})

test_that("check_feasibility rho ranges are valid", {
  feas <- check_feasibility(n_items = 25, model = "rasch", seed = 42,
                             M = 3000, verbose = FALSE)

  # All rho values in (0, 1)
  expect_true(all(feas$rho_range_info > 0 & feas$rho_range_info < 1))
  expect_true(all(feas$rho_range_msem > 0 & feas$rho_range_msem < 1))

  # rho_tilde range should be ordered (monotone)
  expect_lte(feas$rho_range_info[1], feas$rho_range_info[2])

  # msem range should be ordered
  expect_lte(feas$rho_range_msem[1], feas$rho_range_msem[2])
})

test_that("check_feasibility: Jensen's inequality holds", {
  feas <- check_feasibility(n_items = 25, model = "rasch", seed = 42,
                             M = 3000, verbose = FALSE)

  # Upper bound of info metric should be >= upper bound of msem
  expect_gte(feas$rho_range_info[2], feas$rho_range_msem[2])
})

test_that("check_feasibility stores correct metadata", {
  feas <- check_feasibility(n_items = 30, model = "2pl",
                             latent_shape = "bimodal", seed = 42,
                             M = 2000, verbose = FALSE)
  expect_equal(feas$n_items, 30)
  expect_equal(feas$model, "2pl")
  expect_equal(feas$latent_shape, "bimodal")
  expect_equal(feas$M, 2000)
  expect_gt(feas$theta_var, 0)
})

test_that("check_feasibility validates inputs", {
  expect_error(check_feasibility(n_items = -5, verbose = FALSE),
               "`n_items` must be a positive integer")
  expect_error(check_feasibility(n_items = 20, c_bounds = c(5, 1), verbose = FALSE),
               "`c_bounds`")
  expect_error(check_feasibility(n_items = 20, M = -10, verbose = FALSE),
               "`M` must be a positive integer")
})

test_that("check_feasibility print method works", {
  feas <- check_feasibility(n_items = 20, model = "rasch", seed = 42,
                             M = 2000, verbose = FALSE)
  expect_output(print(feas), "Feasibility Check")
  expect_output(print(feas), "rho_tilde")
  expect_output(print(feas), "rho_bar")
})

test_that("check_feasibility with wider c_bounds gives wider or equal range", {
  feas_narrow <- check_feasibility(n_items = 20, c_bounds = c(0.5, 3),
                                    seed = 42, M = 2000, verbose = FALSE)
  feas_wide   <- check_feasibility(n_items = 20, c_bounds = c(0.1, 10),
                                    seed = 42, M = 2000, verbose = FALSE)
  # Wider bounds should give wider or equal info range
  expect_lte(feas_wide$rho_range_info[1], feas_narrow$rho_range_info[1] + 1e-6)
  expect_gte(feas_wide$rho_range_info[2], feas_narrow$rho_range_info[2] - 1e-6)
})


# =============================================================================
# Tests for rho_curve()
# =============================================================================

test_that("rho_curve returns correct class and structure", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", seed = 42,
                  M = 2000, plot = FALSE)
  expect_s3_class(rc, "rho_curve")
  expect_s3_class(rc, "data.frame")
  expect_true("c" %in% names(rc))
})

test_that("rho_curve with metric='both' has both columns", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", metric = "both",
                  seed = 42, M = 2000, plot = FALSE)
  expect_true("rho_tilde" %in% names(rc))
  expect_true("rho_bar" %in% names(rc))
  expect_equal(nrow(rc), 10)
})

test_that("rho_curve with metric='info' only has rho_tilde", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", metric = "info",
                  seed = 42, M = 2000, plot = FALSE)
  expect_true("rho_tilde" %in% names(rc))
  expect_false("rho_bar" %in% names(rc))
})

test_that("rho_curve with metric='msem' only has rho_bar", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", metric = "msem",
                  seed = 42, M = 2000, plot = FALSE)
  expect_true("rho_bar" %in% names(rc))
  expect_false("rho_tilde" %in% names(rc))
})

test_that("rho_curve: rho_tilde is monotone increasing", {
  rc <- rho_curve(c_values = seq(0.3, 5, length.out = 30),
                  n_items = 25, model = "rasch", metric = "info",
                  seed = 42, M = 3000, plot = FALSE)
  diffs <- diff(rc$rho_tilde)
  expect_true(all(diffs >= -1e-8),
              info = "rho_tilde should be monotone increasing in c")
})

test_that("rho_curve: all values in (0, 1)", {
  rc <- rho_curve(c_values = seq(0.3, 5, length.out = 20),
                  n_items = 25, model = "rasch", metric = "both",
                  seed = 42, M = 2000, plot = FALSE)
  expect_true(all(rc$rho_tilde > 0 & rc$rho_tilde < 1))
  expect_true(all(rc$rho_bar > 0 & rc$rho_bar < 1))
})

test_that("rho_curve: Jensen's inequality holds pointwise", {
  rc <- rho_curve(c_values = seq(0.3, 5, length.out = 20),
                  n_items = 25, model = "rasch", metric = "both",
                  seed = 42, M = 2000, plot = FALSE)
  expect_true(all(rc$rho_tilde >= rc$rho_bar - 1e-10))
})

test_that("rho_curve validates inputs", {
  expect_error(rho_curve(c_values = -1, n_items = 20, plot = FALSE),
               "`c_values`")
  expect_error(rho_curve(c_values = seq(0.5, 3, length.out = 10),
                         n_items = -5, plot = FALSE),
               "`n_items`")
})

test_that("rho_curve print method works", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", seed = 42,
                  M = 2000, plot = FALSE)
  expect_output(print(rc), "Reliability Curve")
})

test_that("rho_curve attributes are set correctly", {
  rc <- rho_curve(c_values = seq(0.5, 3, length.out = 10),
                  n_items = 20, model = "rasch", metric = "both",
                  seed = 42, M = 2000, plot = FALSE)
  expect_equal(attr(rc, "metric"), "both")
  expect_equal(attr(rc, "n_items"), 20)
  expect_equal(attr(rc, "model"), "rasch")
  expect_equal(attr(rc, "latent_shape"), "normal")
  expect_gt(attr(rc, "theta_var"), 0)
})


# =============================================================================
# Tests for coef.eqc_result
# =============================================================================

test_that("coef.eqc_result returns correct structure", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  cf <- coef(eqc_res)
  expect_s3_class(cf, "data.frame")
  expect_named(cf, c("item_id", "beta", "lambda_base", "lambda_scaled", "c_star"))
  expect_equal(nrow(cf), 20)
})

test_that("coef.eqc_result values match stored values", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  cf <- coef(eqc_res)
  expect_equal(cf$beta, eqc_res$beta_vec)
  expect_equal(cf$lambda_base, eqc_res$lambda_base)
  expect_equal(cf$lambda_scaled, eqc_res$lambda_scaled)
  expect_equal(unique(cf$c_star), eqc_res$c_star)
  expect_equal(cf$item_id, 1:20)
})

test_that("coef.eqc_result: lambda_scaled = lambda_base * c_star", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  cf <- coef(eqc_res)
  expect_equal(cf$lambda_scaled, cf$lambda_base * cf$c_star, tolerance = 1e-10)
})


# =============================================================================
# Tests for predict.eqc_result
# =============================================================================

test_that("predict.eqc_result returns achieved_rho when newdata is NULL", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  pred <- predict(eqc_res)
  expect_equal(pred, eqc_res$achieved_rho)
})

test_that("predict.eqc_result returns valid reliabilities for newdata", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  c_vals <- c(0.5, 1.0, 1.5, 2.0)
  pred <- predict(eqc_res, newdata = c_vals)

  expect_length(pred, 4)
  expect_true(all(pred > 0 & pred < 1))
  expect_true(!is.null(names(pred)))
})

test_that("predict.eqc_result is monotone for info metric", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  c_vals <- seq(0.3, 3, length.out = 10)
  pred <- predict(eqc_res, newdata = c_vals)
  diffs <- diff(pred)
  expect_true(all(diffs >= -1e-8), info = "Predictions should be monotone for info metric")
})

test_that("predict.eqc_result errors for negative c values", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  expect_error(predict(eqc_res, newdata = c(-1, 1)),
               "All values in `newdata` must be positive")
})

test_that("predict.eqc_result at c_star approximately equals achieved_rho", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  pred_at_cstar <- predict(eqc_res, newdata = eqc_res$c_star)
  expect_equal(unname(pred_at_cstar), eqc_res$achieved_rho, tolerance = 1e-6)
})


# =============================================================================
# Tests for coef.sac_result
# =============================================================================

test_that("coef.sac_result returns correct structure", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  cf <- coef(sac_res)
  expect_s3_class(cf, "data.frame")
  expect_named(cf, c("item_id", "beta", "lambda_base", "lambda_scaled", "c_star"))
  expect_equal(nrow(cf), 15)
})

test_that("coef.sac_result values match stored values", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  cf <- coef(sac_res)
  expect_equal(cf$beta, sac_res$beta_vec)
  expect_equal(cf$lambda_base, sac_res$lambda_base)
  expect_equal(cf$lambda_scaled, sac_res$lambda_scaled)
  expect_equal(unique(cf$c_star), sac_res$c_star)
})

test_that("coef.sac_result: lambda_scaled = lambda_base * c_star", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  cf <- coef(sac_res)
  expect_equal(cf$lambda_scaled, cf$lambda_base * cf$c_star, tolerance = 1e-10)
})


# =============================================================================
# Tests for predict.sac_result
# =============================================================================

test_that("predict.sac_result returns achieved_rho when newdata is NULL", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  pred <- predict(sac_res)
  expect_equal(pred, sac_res$achieved_rho)
})

test_that("predict.sac_result returns valid reliabilities for newdata", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  c_vals <- c(0.5, 1.0, 1.5, 2.0)
  pred <- predict(sac_res, newdata = c_vals)

  expect_length(pred, 4)
  expect_true(all(pred > 0 & pred < 1))
  expect_true(!is.null(names(pred)))
})

test_that("predict.sac_result errors for negative c values", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  expect_error(predict(sac_res, newdata = c(-1, 1)),
               "All values in `newdata` must be positive")
})

test_that("predict.sac_result accepts custom theta_vec", {
  skip_on_cran()
  sac_res <- sac_calibrate(target_rho = 0.75, n_items = 15,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  set.seed(99)
  custom_theta <- rnorm(5000)
  pred <- predict(sac_res, newdata = c(1.0, 2.0), theta_vec = custom_theta)
  expect_length(pred, 2)
  expect_true(all(pred > 0 & pred < 1))
})


# =============================================================================
# Integration test: coef/predict consistency between EQC and SAC
# =============================================================================

test_that("coef format is consistent between EQC and SAC", {
  skip_on_cran()
  eqc_res <- eqc_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", seed = 42, M = 3000)
  sac_res <- sac_calibrate(target_rho = 0.80, n_items = 20,
                             model = "rasch", n_iter = 100,
                             seed = 42, M_per_iter = 300, M_pre = 2000)
  cf_eqc <- coef(eqc_res)
  cf_sac <- coef(sac_res)
  expect_equal(names(cf_eqc), names(cf_sac))
  expect_equal(nrow(cf_eqc), nrow(cf_sac))
})
