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

  test_that("provenance records calibration and simulation settings", {
    sim <- simulate_response_data(
      result = eqc_for_val,
      n_persons = 50,
      latent_shape = "normal",
      seed = 123
    )

    expect_true("provenance" %in% names(sim))
    expect_equal(sim$provenance$result_class, "eqc_result")
    expect_equal(sim$provenance$c_star, eqc_for_val$c_star)
    expect_equal(sim$provenance$target_rho, eqc_for_val$target_rho)
    expect_equal(sim$provenance$achieved_rho, eqc_for_val$achieved_rho)
    expect_equal(sim$provenance$metric, eqc_for_val$metric)
    expect_equal(sim$provenance$model, eqc_for_val$model)
    expect_equal(sim$provenance$n_items, eqc_for_val$n_items)
    expect_equal(sim$provenance$calibration_status,
                 eqc_for_val$misc$root_status)
    expect_type(sim$provenance$status_flags, "character")
    expect_equal(sim$provenance$status_flags,
                 as.character(eqc_for_val$misc$root_status))
    expect_equal(sim$provenance$simulation_seed, 123)
    expect_equal(sim$provenance$n_persons, 50L)
    expect_equal(sim$provenance$latent_shape, "normal")
    expect_true(is.character(sim$provenance$calibration_call))
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
  expect_equal(sim$beta, sac_for_val$beta_vec)
  expect_equal(sim$lambda, sac_for_val$lambda_scaled)
  expect_equal(sim$lambda, sac_for_val$items_calib$data$lambda)
  expect_equal(sim$provenance$result_class, "sac_result")
  expect_equal(sim$provenance$item_design, sac_for_val$item_design)
  expect_equal(sim$provenance$calibration_status,
               sac_for_val$calibration_status)
  expect_type(sim$provenance$status_flags, "character")
  expect_equal(sim$provenance$status_flags,
               sac_for_val$convergence$status_flags)
})

test_that("simulate_response_data preserves multi-flag SAC provenance", {
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
  sac_for_val$calibration_status <- "not_converged"
  sac_for_val$convergence$status <- "not_converged"
  sac_for_val$convergence$status_flags <- c("not_converged", "large_achieved_gap")

  sim <- simulate_response_data(
    result = sac_for_val,
    n_persons = 20,
    latent_shape = "normal",
    seed = 456
  )

  expect_length(sim$provenance$calibration_status, 1L)
  expect_equal(sim$provenance$calibration_status, "not_converged")
  expect_type(sim$provenance$status_flags, "character")
  expect_equal(sim$provenance$status_flags,
               c("not_converged", "large_achieved_gap"))
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

test_that("simulate_response_data validates n_persons", {
  expect_error(
    simulate_response_data(result = eqc_for_val, n_persons = 0),
    "n_persons"
  )
  expect_error(
    simulate_response_data(result = eqc_for_val, n_persons = 10.5),
    "n_persons"
  )
  expect_error(
    simulate_response_data(result = eqc_for_val, n_persons = NA_real_),
    "n_persons"
  )
})

test_that("simulate_response_data validates result design fields", {
  skeletal <- structure(
    list(beta_vec = rep(0, 4), lambda_scaled = rep(1, 4)),
    class = "eqc_result"
  )
  expect_error(
    simulate_response_data(result = skeletal, n_persons = 10),
    "missing required field"
  )

  missing_n_items <- eqc_for_val
  missing_n_items$n_items <- NULL
  expect_error(
    simulate_response_data(result = missing_n_items, n_persons = 10),
    "missing required field.*n_items"
  )

  missing_status <- eqc_for_val
  missing_status$misc$root_status <- NULL
  expect_error(
    simulate_response_data(result = missing_status, n_persons = 10),
    "misc\\$root_status"
  )

  missing_lambda <- eqc_for_val
  missing_lambda$lambda_scaled <- NULL
  expect_error(
    simulate_response_data(result = missing_lambda, n_persons = 10),
    "missing required field.*lambda_scaled"
  )

  mismatched_length <- eqc_for_val
  mismatched_length$lambda_scaled <- mismatched_length$lambda_scaled[-1]
  expect_error(
    simulate_response_data(result = mismatched_length, n_persons = 10),
    "same length"
  )

  nonpositive_lambda <- eqc_for_val
  nonpositive_lambda$lambda_scaled[1] <- 0
  expect_error(
    simulate_response_data(result = nonpositive_lambda, n_persons = 10),
    "positive values"
  )

  inconsistent_scale <- eqc_for_val
  inconsistent_scale$lambda_scaled <-
    inconsistent_scale$lambda_scaled * 1.2
  expect_error(
    simulate_response_data(result = inconsistent_scale, n_persons = 10),
    "lambda_base.*c_star"
  )

  mismatched_n_items <- eqc_for_val
  mismatched_n_items$n_items <- mismatched_n_items$n_items + 1L
  expect_error(
    simulate_response_data(result = mismatched_n_items, n_persons = 10),
    "n_items.*match"
  )

  inconsistent_items <- eqc_for_val
  inconsistent_items$items_calib$data$lambda[1] <-
    inconsistent_items$items_calib$data$lambda[1] * 2
  expect_error(
    simulate_response_data(result = inconsistent_items, n_persons = 10),
    "items_calib\\$data\\$lambda.*lambda_scaled"
  )

  noncanonical_metric <- eqc_for_val
  noncanonical_metric$metric <- "tilde"
  expect_error(
    simulate_response_data(result = noncanonical_metric, n_persons = 10),
    "canonical stored metric"
  )
})

test_that("simulate_response_data accepts schema-complete legacy spc_result", {
  spc_legacy <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 8,
    model = "rasch",
    item_source = "parametric",
    n_iter = 20L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))
  class(spc_legacy) <- c("spc_result", "list")

  sim <- simulate_response_data(spc_legacy, n_persons = 20, seed = 99)
  expect_equal(ncol(sim$response_matrix), 8)
  expect_equal(sim$provenance$result_class, "spc_result")
  expect_equal(sim$provenance$calibration_status,
               spc_legacy$calibration_status)
  expect_type(sim$provenance$status_flags, "character")
  expect_equal(sim$provenance$status_flags,
               spc_legacy$convergence$status_flags)
})

test_that("simulate_response_data falls back for empty legacy status_flags", {
  spc_legacy <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 8,
    model = "rasch",
    item_source = "parametric",
    n_iter = 20L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))
  class(spc_legacy) <- c("spc_result", "list")
  spc_legacy$calibration_status <- "not_converged"
  spc_legacy$convergence$status <- "not_converged"
  spc_legacy$convergence$status_flags <- character()

  sim <- simulate_response_data(spc_legacy, n_persons = 20, seed = 99)
  expect_equal(sim$provenance$calibration_status, "not_converged")
  expect_equal(sim$provenance$status_flags, "not_converged")
})

test_that("simulate_response_data rejects status-incomplete legacy spc_result", {
  spc_legacy <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 8,
    model = "rasch",
    item_source = "parametric",
    n_iter = 20L,
    M_per_iter = 100L,
    M_pre = 1000L,
    seed = 42
  )))
  class(spc_legacy) <- c("spc_result", "list")
  spc_legacy$calibration_status <- NULL
  spc_legacy$convergence$status <- NULL

  expect_error(
    simulate_response_data(spc_legacy, n_persons = 20, seed = 99),
    "calibration_status.*convergence\\$status"
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
                         "agreement", "target_rho",
                         "achieved_eqc", "achieved_sac",
                         "achieved_diff_abs", "metric_eqc", "metric_sac",
                         "model_eqc", "model_sac",
                         "n_items_eqc", "n_items_sac",
                         "eqc_status", "sac_status", "sac_status_flags")
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
    expect_equal(comp$achieved_diff_abs,
                 abs(comp$achieved_eqc - comp$achieved_sac))
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

  test_that("validates comparison result schemas before computing", {
    bad_eqc <- eqc_for_val
    bad_eqc$c_star <- NA_real_
    expect_error(
      compare_eqc_sac(bad_eqc, sac_for_comp, verbose = FALSE),
      "eqc_result\\$c_star"
    )

    bad_eqc <- eqc_for_val
    bad_eqc$achieved_rho <- 1.2
    expect_error(
      compare_eqc_sac(bad_eqc, sac_for_comp, verbose = FALSE),
      "achieved_rho"
    )

    bad_sac <- sac_for_comp
    bad_sac$metric <- NULL
    expect_error(
      compare_eqc_sac(eqc_for_val, bad_sac, verbose = FALSE),
      "missing required field.*metric"
    )

    bad_sac <- sac_for_comp
    bad_sac$model <- "3pl"
    expect_error(
      compare_eqc_sac(eqc_for_val, bad_sac, verbose = FALSE),
      "sac_result\\$model"
    )

    bad_sac <- sac_for_comp
    bad_sac$metric <- "bar"
    expect_error(
      compare_eqc_sac(eqc_for_val, bad_sac, verbose = FALSE),
      "canonical stored metric"
    )

    expect_error(
      compare_eqc_sac(eqc_for_val, sac_for_comp, verbose = NA),
      "verbose"
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

  test_that("warns when calibration statuses indicate boundary or instability", {
    eqc_boundary <- eqc_for_val
    eqc_boundary$misc$root_status <- "above_upper"
    sac_unstable <- sac_for_comp
    sac_unstable$calibration_status <- "hit_upper_bound"
    sac_unstable$convergence$status <- "hit_upper_bound"
    sac_unstable$convergence$status_flags <- "hit_upper_bound"

    warnings <- character()
    comp <- NULL
    withCallingHandlers(
      suppressMessages(
        comp <- compare_eqc_sac(eqc_boundary, sac_unstable, verbose = FALSE)
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    expect_true(any(grepl("EQC root status", warnings)))
    expect_true(any(grepl("SAC calibration status", warnings)))
    expect_equal(comp$eqc_status, "above_upper")
    expect_equal(comp$sac_status, "hit_upper_bound")
    expect_true("hit_upper_bound" %in% comp$sac_status_flags)
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

test_that("compare_eqc_spc forwards to compare_eqc_sac after deprecation", {
  sac_dep <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "info",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))

  primary <- compare_eqc_sac(eqc_for_val, sac_dep, verbose = FALSE)
  alias <- suppressWarnings(compare_eqc_spc(eqc_for_val, sac_dep,
                                            verbose = FALSE))

  expect_identical(alias, primary)
})

test_that("compare_eqc_sac accepts schema-complete legacy spc_result", {
  spc_for_comp <- suppressWarnings(suppressMessages(sac_calibrate(
    target_rho = 0.80,
    n_items = 15,
    model = "rasch",
    item_source = "parametric",
    reliability_metric = "info",
    n_iter = 50L,
    M_per_iter = 200L,
    M_pre = 2000L,
    seed = 42
  )))
  class(spc_for_comp) <- c("spc_result", "list")

  comp <- suppressWarnings(compare_eqc_sac(eqc_for_val, spc_for_comp,
                                           verbose = FALSE))

  expect_equal(comp$sac_status, spc_for_comp$calibration_status)
  expect_equal(comp$sac_status_flags,
               spc_for_comp$convergence$status_flags)
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
  tam_rel <- NULL
  capture.output(
    tam_rel <- compute_reliability_tam(sim$response_matrix, model = "rasch")
  )

  expect_true("rel_wle" %in% names(tam_rel))
  expect_true("rel_eap" %in% names(tam_rel))
  expect_true("mod" %in% names(tam_rel))
  expect_true("wle" %in% names(tam_rel))
  expect_true(tam_rel$rel_wle > 0 && tam_rel$rel_wle < 1)
  expect_true(tam_rel$rel_eap > 0 && tam_rel$rel_eap < 1)
  # Seeded fixture check only; TAM does not guarantee a universal EAP/WLE ordering.
  expect_true(tam_rel$rel_eap >= tam_rel$rel_wle - 0.01)
})
