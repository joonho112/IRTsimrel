# =============================================================================
# test-s3_methods.R
# =============================================================================
# Comprehensive tests for all S3 methods across the four main classes:
#   latent_G, item_params, eqc_result, sac_result
# Also tests compare_shapes() utility.
# =============================================================================

# ---------------------------------------------------------------------------
# Setup: create objects for all S3 method tests
# ---------------------------------------------------------------------------

latent_obj <- sim_latentG(n = 200L, shape = "bimodal", seed = 1)

item_obj <- sim_item_params(
  n_items = 15, model = "2pl", source = "parametric",
  method = "copula", seed = 1
)

eqc_obj <- suppressMessages(eqc_calibrate(
  target_rho = 0.80, n_items = 15, model = "rasch",
  item_source = "parametric", M = 2000L, seed = 42
))

sac_obj <- suppressWarnings(suppressMessages(sac_calibrate(
  target_rho = 0.75, n_items = 15, model = "rasch",
  item_source = "parametric",
  n_iter = 50L, M_per_iter = 200L, M_pre = 2000L, seed = 42
)))


# ===========================================================================
# latent_G S3 methods
# ===========================================================================

describe("latent_G S3 methods", {

  test_that("print.latent_G does not error and returns input invisibly", {
    out <- expect_invisible(print(latent_obj))
    expect_identical(out, latent_obj)
  })

  test_that("print.latent_G produces informative output", {
    output <- capture.output(print(latent_obj))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Latent", output, ignore.case = TRUE)))
  })

  test_that("summary.latent_G returns object of class 'summary.latent_G'", {
    s <- summary(latent_obj)
    expect_s3_class(s, "summary.latent_G")
    expect_true("shape" %in% names(s))
    expect_true("n" %in% names(s))
    expect_true("sample_mean" %in% names(s))
    expect_true("sample_sd" %in% names(s))
    expect_true("quantiles" %in% names(s))
  })

  test_that("print.summary.latent_G does not error", {
    s <- summary(latent_obj)
    output <- capture.output(print(s))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Summary", output, ignore.case = TRUE)))
  })

  test_that("as.numeric.latent_G returns numeric vector of correct length", {
    # Note: as.numeric() primitive doesn't dispatch for list subclasses in R;
    # the S3 method must be called directly or via UseMethod dispatch.
    num <- as.numeric.latent_G(latent_obj)
    expect_true(is.numeric(num))
    expect_length(num, 200L)
    expect_identical(num, latent_obj$theta)
  })

  test_that("plot.latent_G with ggplot2 returns ggplot object", {
    skip_if_not_installed("ggplot2")
    p <- plot(latent_obj)
    expect_s3_class(p, "gg")
  })

  test_that("plot.latent_G base R fallback does not error", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    result <- tryCatch(plot(latent_obj), error = function(e) e)
    expect_false(inherits(result, "error"))
  })
})


# ===========================================================================
# item_params S3 methods
# ===========================================================================

describe("item_params S3 methods", {

  test_that("print.item_params does not error and returns invisibly", {
    out <- expect_invisible(print(item_obj))
    expect_identical(out, item_obj)
  })

  test_that("print.item_params produces informative output", {
    output <- capture.output(print(item_obj))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Item Parameters", output, ignore.case = TRUE)))
  })

  test_that("summary.item_params returns object of class 'summary.item_params'", {
    s <- summary(item_obj)
    expect_s3_class(s, "summary.item_params")
    expect_true("model" %in% names(s))
    expect_true("source" %in% names(s))
    expect_true("n_items" %in% names(s))
    expect_true("n_forms" %in% names(s))
    expect_true("scale" %in% names(s))
    expect_true("centered" %in% names(s))
    expect_true("beta_summary" %in% names(s))
    expect_true("lambda_summary" %in% names(s))
    expect_true("achieved_cors" %in% names(s))
    # Check beta_summary subfields
    expect_true("mean" %in% names(s$beta_summary))
    expect_true("sd" %in% names(s$beta_summary))
    expect_true("quantiles" %in% names(s$beta_summary))
  })

  test_that("print.summary.item_params does not error", {
    s <- summary(item_obj)
    output <- capture.output(print(s))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Summary", output, ignore.case = TRUE)))
  })

  test_that("summary.item_params for rasch model has NULL lambda_summary", {
    item_rasch <- sim_item_params(
      n_items = 10, model = "rasch", source = "parametric", seed = 1
    )
    s <- summary(item_rasch)
    expect_null(s$lambda_summary)
    expect_null(s$achieved_cors)
  })

  test_that("as.data.frame.item_params returns data.frame with correct columns", {
    df <- as.data.frame(item_obj)
    expect_s3_class(df, "data.frame")
    expected_cols <- c("form_id", "item_id", "beta", "lambda", "lambda_unscaled")
    expect_true(all(expected_cols %in% names(df)))
    expect_equal(nrow(df), 15)
  })

  test_that("as.data.frame.item_params with multiple forms returns all rows", {
    items_multi <- sim_item_params(
      n_items = 10, model = "rasch", source = "parametric",
      n_forms = 3, seed = 1
    )
    df <- as.data.frame(items_multi)
    expect_equal(nrow(df), 30)
    expect_equal(length(unique(df$form_id)), 3)
  })

  test_that("plot.item_params scatter type works with ggplot2", {
    skip_if_not_installed("ggplot2")
    p <- plot(item_obj, type = "scatter")
    expect_s3_class(p, "gg")
  })

  test_that("plot.item_params density type works with ggplot2", {
    skip_if_not_installed("ggplot2")
    p <- plot(item_obj, type = "density")
    # Returns ggplot or patchwork object
    expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
  })

  test_that("plot.item_params base R fallback does not error", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    result <- tryCatch(plot(item_obj), error = function(e) e)
    expect_false(inherits(result, "error"))
  })
})


# ===========================================================================
# eqc_result S3 methods
# ===========================================================================

describe("eqc_result S3 methods", {

  test_that("print.eqc_result does not error and returns invisibly", {
    out <- expect_invisible(print(eqc_obj))
    expect_identical(out, eqc_obj)
  })

  test_that("print.eqc_result produces informative output", {
    output <- capture.output(print(eqc_obj))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Empirical Quadrature Calibration|EQC", output,
                          ignore.case = TRUE)))
  })

  test_that("summary.eqc_result returns object of class 'summary.eqc_result'", {
    s <- summary(eqc_obj)
    expect_s3_class(s, "summary.eqc_result")
    expect_true("c_star" %in% names(s))
    expect_true("target_rho" %in% names(s))
    expect_true("achieved_rho" %in% names(s))
    expect_true("metric" %in% names(s))
    expect_true("model" %in% names(s))
    expect_true("n_items" %in% names(s))
    expect_true("M" %in% names(s))
    expect_true("root_status" %in% names(s))
    expect_true("theta_var" %in% names(s))
  })

  test_that("print.summary.eqc_result does not error and returns invisibly", {
    s <- summary(eqc_obj)
    out <- expect_invisible(print(s))
    expect_s3_class(out, "summary.eqc_result")
    output <- capture.output(print(s))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Summary.*EQC", output, ignore.case = TRUE)))
  })
})


# ===========================================================================
# sac_result S3 methods
# ===========================================================================

describe("sac_result S3 methods", {

  test_that("print.sac_result does not error and returns invisibly", {
    out <- expect_invisible(print(sac_obj))
    expect_identical(out, sac_obj)
  })

  test_that("print.sac_result produces informative output", {
    output <- capture.output(print(sac_obj))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Stochastic Approximation Calibration|SAC", output,
                          ignore.case = TRUE)))
  })

  test_that("summary.sac_result returns object of class 'summary.sac_result'", {
    s <- summary(sac_obj)
    expect_s3_class(s, "summary.sac_result")
    expect_true("c_star" %in% names(s))
    expect_true("target_rho" %in% names(s))
    expect_true("achieved_rho" %in% names(s))
    expect_true("metric" %in% names(s))
    expect_true("model" %in% names(s))
    expect_true("n_items" %in% names(s))
    expect_true("n_iter" %in% names(s))
    expect_true("M_per_iter" %in% names(s))
    expect_true("M_pre" %in% names(s))
    expect_true("init_method" %in% names(s))
    expect_true("convergence" %in% names(s))
    expect_true("burn_in" %in% names(s))
  })

  test_that("print.summary.sac_result does not error and returns invisibly", {
    s <- summary(sac_obj)
    out <- expect_invisible(print(s))
    expect_s3_class(out, "summary.sac_result")
    output <- capture.output(print(s))
    expect_true(length(output) > 0)
    expect_true(any(grepl("Summary.*SAC", output, ignore.case = TRUE)))
  })

  test_that("plot.sac_result type='c' with ggplot2 returns ggplot", {
    skip_if_not_installed("ggplot2")
    p <- plot(sac_obj, type = "c")
    expect_s3_class(p, "gg")
  })

  test_that("plot.sac_result type='rho' with ggplot2 returns ggplot", {
    skip_if_not_installed("ggplot2")
    p <- plot(sac_obj, type = "rho")
    expect_s3_class(p, "gg")
  })

  test_that("plot.sac_result type='both' with ggplot2 works", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("patchwork")
    p <- plot(sac_obj, type = "both")
    expect_true(inherits(p, "gg") || inherits(p, "patchwork"))
  })

  test_that("plot.sac_result base R fallback does not error", {
    pdf(NULL)
    on.exit(dev.off(), add = TRUE)
    result <- tryCatch(plot(sac_obj, type = "c"), error = function(e) e)
    expect_false(inherits(result, "error"))
  })
})


# ===========================================================================
# compare_shapes() utility
# ===========================================================================

describe("compare_shapes()", {

  test_that("compare_shapes returns ggplot with ggplot2 installed", {
    skip_if_not_installed("ggplot2")
    p <- compare_shapes(
      n = 500,
      shapes = c("normal", "bimodal", "skew_pos"),
      seed = 42
    )
    expect_s3_class(p, "gg")
  })

  test_that("compare_shapes with 2 shapes works", {
    skip_if_not_installed("ggplot2")
    p <- compare_shapes(n = 100, shapes = c("normal", "uniform"), seed = 1)
    expect_s3_class(p, "gg")
  })

  test_that("compare_shapes preserves external RNG state when seed is set", {
    skip_if_not_installed("ggplot2")
    set.seed(111)
    state_before <- .Random.seed
    p <- compare_shapes(n = 100, shapes = c("normal"), seed = 42)
    state_after <- .Random.seed
    expect_identical(state_before, state_after)
  })
})


# ===========================================================================
# Cross-class integration: objects used in the right contexts
# ===========================================================================

describe("Cross-class integration", {

  test_that("eqc_result can be passed to simulate_response_data", {
    sim <- simulate_response_data(
      result = eqc_obj, n_persons = 50, latent_shape = "normal", seed = 1
    )
    expect_equal(nrow(sim$response_matrix), 50)
    expect_equal(ncol(sim$response_matrix), 15)
  })

  test_that("sac_result can be passed to simulate_response_data", {
    sim <- simulate_response_data(
      result = sac_obj, n_persons = 50, latent_shape = "normal", seed = 1
    )
    expect_equal(nrow(sim$response_matrix), 50)
    expect_equal(ncol(sim$response_matrix), 15)
  })

  test_that("eqc_result can be passed as c_init to sac_calibrate (warm start)", {
    sac2 <- suppressWarnings(suppressMessages(sac_calibrate(
      target_rho = 0.80, n_items = 15, model = "rasch",
      item_source = "parametric",
      c_init = eqc_obj,
      n_iter = 50L, M_per_iter = 200L, M_pre = 2000L, seed = 99
    )))
    expect_s3_class(sac2, "sac_result")
    expect_equal(sac2$init_method, "eqc_warm_start")
  })
})
