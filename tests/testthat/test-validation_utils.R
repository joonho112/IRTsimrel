# =============================================================================
# test-validation_utils.R
# =============================================================================
# Comprehensive tests for simulate_response_data(), compare_eqc_sac(),
# compare_eqc_spc(), and compute_reliability_tam()
# =============================================================================

# ---------------------------------------------------------------------------
# Helper: create a quick EQC result for testing
# Each test_that block that needs it constructs its own data, but we create
# a shared EQC result at the top level to keep tests fast.
# ---------------------------------------------------------------------------

eqc_for_val <- suppressMessages(eqc_calibrate(
  target_rho = 0.80,
  n_items = 15,
  model = "rasch",
  item_source = "parametric",
  M = 2000L,
  seed = 42
))


# ---------------------------------------------------------------------------
# simulate_response_data: dimensions and content
# ---------------------------------------------------------------------------

describe("simulate_response_data() dimensions and content", {

  test_that("response matrix has correct dimensions: n_persons x n_items", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 200,
      latent_shape = "normal",
      seed = 123
    )
    expect_equal(nrow(sim$response_matrix), 200)
    expect_equal(ncol(sim$response_matrix), 15)
  })

  test_that("response matrix is binary: only 0s and 1s", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 100,
      latent_shape = "normal",
      seed = 123
    )
    expect_true(all(sim$response_matrix %in% c(0L, 1L)))
  })

  test_that("column names are paste0('item', 1:n_items)", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 50,
      latent_shape = "normal",
      seed = 123
    )
    expected_names <- paste0("item", 1:15)
    expect_equal(colnames(sim$response_matrix), expected_names)
  })

  test_that("theta, beta, and lambda have correct lengths", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 300,
      latent_shape = "normal",
      seed = 123
    )
    expect_length(sim$theta, 300)
    expect_length(sim$beta, 15)
    expect_length(sim$lambda, 15)
  })

  test_that("beta and lambda come from the result object", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 50,
      latent_shape = "normal",
      seed = 123
    )
    # beta should match the result's beta_vec
    expect_equal(sim$beta, eqc_for_val$beta_vec)
    # lambda should match the result's lambda_scaled
    expect_equal(sim$lambda, eqc_for_val$lambda_scaled)
  })
})


# ---------------------------------------------------------------------------
# Seed reproducibility
# ---------------------------------------------------------------------------

test_that("simulate_response_data with same seed gives identical results", {
  sim1 <- simulate_response_data(
    result = eqc_for_val,
    n_persons = 100,
    latent_shape = "normal",
    seed = 789
  )
  sim2 <- simulate_response_data(
    result = eqc_for_val,
    n_persons = 100,
    latent_shape = "normal",
    seed = 789
  )
  expect_identical(sim1$response_matrix, sim2$response_matrix)
  expect_identical(sim1$theta, sim2$theta)
})


# ---------------------------------------------------------------------------
# simulate_response_data works with sac_result class
# ---------------------------------------------------------------------------

test_that("simulate_response_data works with sac_result class", {
  sac_for_val <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  sim <- simulate_response_data(
    result = sac_for_val,
    n_persons = 100,
    latent_shape = "normal",
    seed = 456
  )

  expect_equal(nrow(sim$response_matrix), 100)
  expect_equal(ncol(sim$response_matrix), 15)
  expect_true(all(sim$response_matrix %in% c(0L, 1L)))
  expect_length(sim$theta, 100)
})


# ---------------------------------------------------------------------------
# simulate_response_data: invalid result class errors
# ---------------------------------------------------------------------------

test_that("simulate_response_data errors on invalid result class", {
  fake_result <- list(c_star = 1.0, beta_vec = rnorm(10))
  class(fake_result) <- "not_a_valid_result"

  expect_error(
    simulate_response_data(result = fake_result, n_persons = 100),
    "eqc_result.*sac_result"
  )
})


# ---------------------------------------------------------------------------
# compare_eqc_sac: correct fields and computation
# ---------------------------------------------------------------------------

describe("compare_eqc_sac()", {

  # Create both results for comparison
  sac_for_comp <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    c_init = eqc_for_val,
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  test_that("returns list with all expected fields", {
    comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = FALSE))
    expected_fields <- c("c_eqc", "c_sac", "diff_abs", "diff_pct",
                          "agreement", "target_rho")
    for (f in expected_fields) {
      expect_true(f %in% names(comp), info = paste("Missing:", f))
    }
  })

  test_that("c_eqc and c_sac match the input objects", {
    comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = FALSE))
    expect_equal(comp$c_eqc, eqc_for_val$c_star)
    expect_equal(comp$c_sac, sac_for_comp$c_star)
  })

  test_that("diff_abs and diff_pct are computed correctly", {
    comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = FALSE))
    expect_equal(comp$diff_abs, abs(comp$c_eqc - comp$c_sac))
    expect_equal(comp$diff_pct, 100 * comp$diff_abs / comp$c_eqc)
  })

  test_that("agreement is TRUE when diff_pct < 5", {
    comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = FALSE))
    expect_equal(comp$agreement, comp$diff_pct < 5)
  })

  test_that("errors with wrong class for eqc_result", {
    fake <- list(c_star = 1)
    class(fake) <- "not_eqc"
    expect_error(
      compare_eqc_sac(fake, sac_for_comp),
      "eqc_result"
    )
  })

  test_that("errors with wrong class for sac_result", {
    fake <- list(c_star = 1)
    class(fake) <- "not_sac"
    expect_error(
      compare_eqc_sac(eqc_for_val, fake),
      "sac_result"
    )
  })

  test_that("errors when target_rho differs between EQC and SAC", {
    # Create a fake SAC result with different target_rho
    sac_diff <- sac_for_comp
    sac_diff$target_rho <- 0.90
    expect_error(
      compare_eqc_sac(eqc_for_val, sac_diff, verbose = FALSE),
      "target_rho differs"
    )
  })

  test_that("warns when model differs between EQC and SAC", {
    sac_diff_model <- sac_for_comp
    sac_diff_model$model <- "2pl"
    # Muffle non-target warnings (e.g., metric mismatch), only let model warning through
    expect_warning(
      withCallingHandlers(
        compare_eqc_sac(eqc_for_val, sac_diff_model, verbose = FALSE),
        warning = function(w) {
          if (!grepl("Model differs", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      "Model differs"
    )
  })

  test_that("warns when n_items differs between EQC and SAC", {
    sac_diff_items <- sac_for_comp
    sac_diff_items$n_items <- 30L
    # Muffle non-target warnings (e.g., metric mismatch), only let n_items warning through
    expect_warning(
      withCallingHandlers(
        compare_eqc_sac(eqc_for_val, sac_diff_items, verbose = FALSE),
        warning = function(w) {
          if (!grepl("n_items differs", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      ),
      "n_items differs"
    )
  })

  test_that("warns when metric differs between EQC and SAC", {
    # SAC default metric is msem, EQC default is info, so they already differ
    expect_warning(
      compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = FALSE),
      "metric differs"
    )
  })
})


# ---------------------------------------------------------------------------
# compare_eqc_sac: verbose output
# ---------------------------------------------------------------------------

test_that("compare_eqc_sac with verbose=TRUE prints output via message()", {
  sac_v <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  msgs <- capture.output(
    comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, sac_v, verbose = TRUE)),
    type = "message"
  )
  expect_true(length(msgs) > 0)
  expect_true(any(grepl("EQC vs SAC", msgs)))
})


# ---------------------------------------------------------------------------
# Deprecated compare_eqc_spc alias
# ---------------------------------------------------------------------------

test_that("compare_eqc_spc triggers deprecation warning", {
  sac_dep <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  # Expect deprecation warning; suppress other warnings (e.g., metric mismatch)
  expect_warning(
    withCallingHandlers(
      compare_eqc_spc(eqc_for_val, sac_dep, verbose = FALSE),
      warning = function(w) {
        if (!grepl("deprecated", conditionMessage(w), ignore.case = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "deprecated"
  )
})


# ---------------------------------------------------------------------------
# compute_reliability_tam (skip if TAM not available)
# ---------------------------------------------------------------------------

test_that("compute_reliability_tam works and returns WLE and EAP reliabilities", {
  skip_if_not_installed("TAM")

  sim <- simulate_response_data(
    result = eqc_for_val,
    n_persons = 300,
    latent_shape = "normal",
    seed = 99
  )
  tam_rel <- compute_reliability_tam(sim$response_matrix, model = "rasch")

  expect_true("rel_wle" %in% names(tam_rel))
  expect_true("rel_eap" %in% names(tam_rel))
  expect_true("mod" %in% names(tam_rel))
  expect_true("wle" %in% names(tam_rel))
  expect_true(tam_rel$rel_wle > 0 && tam_rel$rel_wle < 1)
  expect_true(tam_rel$rel_eap > 0 && tam_rel$rel_eap < 1)
  # EAP >= WLE (mathematical property)
  expect_true(tam_rel$rel_eap >= tam_rel$rel_wle - 0.01)
})
