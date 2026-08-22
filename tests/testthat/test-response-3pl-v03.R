# =============================================================================
# test-response-3pl-v03.R
# =============================================================================
# 3PL response generation and lazy result-design normalization.
# =============================================================================

make_response_result_v03 <- function(model = "3pl",
                                     guessing = c(0.15, 0.25),
                                     include_guessing = TRUE) {
  beta <- c(-0.4, 0.7)
  lambda_base <- c(0.8, 1.2)
  c_star <- 1.5
  lambda_scaled <- lambda_base * c_star

  base_data <- data.frame(
    form_id = 1L,
    item_id = seq_along(beta),
    beta = beta,
    lambda = lambda_base,
    lambda_unscaled = lambda_base
  )
  calib_data <- base_data
  calib_data$lambda <- lambda_scaled
  if (include_guessing) {
    base_data$guessing <- guessing
    calib_data$guessing <- guessing
  }

  result <- list(
    c_star = c_star,
    target_rho = 0.8,
    achieved_rho = 0.8,
    metric = "info",
    model = model,
    n_items = length(beta),
    theta_var = 1,
    beta_vec = beta,
    lambda_base = lambda_base,
    lambda_scaled = lambda_scaled,
    items_base = structure(
      list(data = base_data, source = "custom"),
      class = c("item_params", "list")
    ),
    items_calib = structure(
      list(data = calib_data, source = "custom"),
      class = c("item_params", "list")
    ),
    schema_version = 1L,
    item_scope = "fixed_form",
    design_signature = list(model = model, n_items = length(beta)),
    call = quote(eqc_calibrate()),
    misc = list(root_status = "uniroot_success")
  )
  if (include_guessing) {
    result$guessing_vec <- guessing
  }
  class(result) <- c("eqc_result", "list")
  result
}

test_that("3PL result normalization drives the common probability kernel", {
  result <- make_response_result_v03()
  design <- IRTsimrel:::.irtsimrel_normalize_result_item_design(result)
  theta <- c(-1, 0, 1)

  observed <- IRTsimrel:::.irtsimrel_probability(
    theta,
    design$bank,
    scale = 1
  )
  eta <- outer(theta, result$beta_vec, "-")
  eta <- sweep(eta, 2L, result$lambda_scaled, "*")
  q <- plogis(eta)
  expected <- sweep(q, 2L, 1 - result$guessing_vec, "*")
  expected <- sweep(expected, 2L, result$guessing_vec, "+")

  expect_equal(observed, expected, tolerance = 1e-14)
  expect_equal(design$guessing, result$guessing_vec)
})

test_that("3PL simulated response rates follow person-specific probabilities", {
  result <- make_response_result_v03()
  sim <- simulate_response_data(result, n_persons = 20000L, seed = 782)
  bank <- IRTsimrel:::.irtsimrel_item_bank(
    result$beta_vec,
    result$lambda_scaled,
    result$guessing_vec,
    model = "3pl"
  )
  expected_rates <- colMeans(IRTsimrel:::.irtsimrel_probability(sim$theta, bank))

  expect_equal(
    unname(colMeans(sim$response_matrix)),
    expected_rates,
    tolerance = 0.015
  )
  expect_equal(sim$guessing, result$guessing_vec)
  expect_equal(sim$provenance$schema_version, 1L)
  expect_equal(sim$provenance$item_scope, "fixed_form")
  expect_equal(sim$provenance$design_signature, result$design_signature)
})

test_that("zero-guessing 3PL response simulation is seed-identical to 2PL", {
  threepl <- make_response_result_v03(model = "3pl", guessing = c(0, 0))
  twopl <- threepl
  twopl$model <- "2pl"
  twopl$design_signature$model <- "2pl"
  twopl$guessing_vec <- NULL
  twopl$items_base$data$guessing <- NULL
  twopl$items_calib$data$guessing <- NULL

  sim_3pl <- simulate_response_data(threepl, n_persons = 500L, seed = 91)
  sim_2pl <- simulate_response_data(twopl, n_persons = 500L, seed = 91)

  expect_identical(sim_3pl$theta, sim_2pl$theta)
  expect_identical(sim_3pl$response_matrix, sim_2pl$response_matrix)
  expect_equal(sim_2pl$guessing, c(0, 0))
})

test_that("stored guessing and beta provenance must be valid and consistent", {
  malformed_guessing <- make_response_result_v03()
  malformed_guessing$guessing_vec <- 0.2
  expect_error(
    simulate_response_data(malformed_guessing, n_persons = 10L, seed = 1),
    "guessing_vec.*length 2"
  )

  mismatched_guessing <- make_response_result_v03()
  mismatched_guessing$items_calib$data$guessing[1] <- 0.3
  expect_error(
    simulate_response_data(mismatched_guessing, n_persons = 10L, seed = 1),
    "items_calib\\$data\\$guessing.*guessing_vec"
  )

  mismatched_beta <- make_response_result_v03()
  mismatched_beta$items_calib$data$beta[1] <- 99
  expect_error(
    simulate_response_data(mismatched_beta, n_persons = 10L, seed = 1),
    "items_calib\\$data\\$beta.*beta_vec"
  )

  missing_guessing <- make_response_result_v03(include_guessing = FALSE)
  expect_error(
    simulate_response_data(missing_guessing, n_persons = 10L, seed = 1),
    "3PL result.*no stored guessing"
  )
})

test_that("legacy Rasch and 2PL designs lazily receive zero guessing", {
  legacy <- make_response_result_v03(model = "2pl", include_guessing = FALSE)
  before <- serialize(legacy, NULL)

  design <- IRTsimrel:::.irtsimrel_normalize_result_item_design(legacy)
  sim <- simulate_response_data(legacy, n_persons = 40L, seed = 12)

  expect_equal(design$guessing, c(0, 0))
  expect_equal(sim$guessing, c(0, 0))
  expect_equal(
    sim$provenance$item_parameter_provenance$guessing,
    "implicit_zero_for_legacy_rasch_2pl"
  )
  expect_identical(serialize(legacy, NULL), before)
})

test_that("literal v0.1 SPC objects give an actionable response error", {
  legacy_spc <- structure(
    list(
      c_star = 1,
      target_rho = 0.8,
      achieved_rho = 0.8,
      metric = "info",
      model = "rasch",
      n_items = 5L,
      theta_var = 1,
      call = quote(spc_calibrate()),
      convergence = list(status = "ok")
    ),
    class = c("spc_result", "list")
  )

  expect_error(
    simulate_response_data(legacy_spc, n_persons = 10L, seed = 1),
    "legacy 'spc_result'.*Rerun calibration.*sac_calibrate"
  )
})

test_that("future or malformed result schemas are not guessed", {
  malformed <- make_response_result_v03()
  malformed$schema_version <- 0.5
  expect_error(
    simulate_response_data(malformed, n_persons = 10L),
    "schema_version.*positive integer"
  )

  future <- make_response_result_v03()
  future$schema_version <- 2L
  expect_error(
    simulate_response_data(future, n_persons = 10L),
    "newer.*Upgrade"
  )
})
