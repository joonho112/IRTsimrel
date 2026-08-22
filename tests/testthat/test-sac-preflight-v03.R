# =============================================================================
# SAC branch-aware deterministic preflight (v0.3)
# =============================================================================

.sac_preflight_theta <- seq(-2.5, 2.5, length.out = 101L)

.sac_preflight_bank <- function(guessing = NULL,
                                model = NULL,
                                beta_shift = 0,
                                lambda_multiplier = 1) {
  beta <- c(-1, -0.3, 0.3, 1) + beta_shift
  lambda <- c(0.8, 1, 1.2, 1.4) * lambda_multiplier
  .irtsimrel_item_bank(
    beta = beta,
    lambda_base = lambda,
    guessing = guessing,
    model = model
  )
}

.sac_preflight_rho <- function(c, bank, metric = "info") {
  .compute_rho_generic(
    c = c,
    theta_vec = .sac_preflight_theta,
    beta_vec = bank$beta,
    lambda_base = bank$lambda_base,
    theta_var = stats::var(.sac_preflight_theta),
    metric_internal = metric,
    guessing = bank$guessing
  )
}

.run_mocked_sac_preflight <- function(objective,
                                      target,
                                      banks = NULL,
                                      split_check = FALSE,
                                      split_log_tolerance = 0.5,
                                      root_controls = list()) {
  if (is.null(banks)) {
    banks <- list(.sac_preflight_bank(model = "2pl"))
  }
  local_mocked_bindings(
    .compute_rho_generic = objective,
    .package = "IRTsimrel",
    .env = parent.frame()
  )
  .irtsimrel_sac_preflight(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    item_banks = banks,
    target_rho = target,
    metric_internal = "info",
    c_bounds = c(0.1, 10),
    split_check = split_check,
    split_log_tolerance = split_log_tolerance,
    root_controls = root_controls
  )
}


test_that("fixed-form 2PL and 3PL preflights recover monotone roots", {
  bank_2pl <- .sac_preflight_bank(model = "2pl")
  bank_3pl <- .sac_preflight_bank(
    guessing = c(0.15, 0.20, 0.20, 0.25),
    model = "3pl"
  )

  target_2pl <- .sac_preflight_rho(0.8, bank_2pl)
  target_3pl <- .sac_preflight_rho(0.7, bank_3pl)
  result_2pl <- .irtsimrel_sac_preflight(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    item_banks = list(bank_2pl),
    target_rho = target_2pl,
    metric_internal = "info",
    c_bounds = c(0.1, 1.5),
    split_check = FALSE
  )
  result_3pl <- .irtsimrel_sac_preflight(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    item_banks = list(bank_3pl),
    target_rho = target_3pl,
    metric_internal = "info",
    c_bounds = c(0.1, 1.5),
    split_check = FALSE
  )
  target_3pl_msem <- .sac_preflight_rho(0.65, bank_3pl, metric = "msem")
  result_3pl_msem <- .irtsimrel_sac_preflight(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    item_banks = list(bank_3pl),
    target_rho = target_3pl_msem,
    metric_internal = "msem",
    c_bounds = c(0.1, 1.3),
    split_check = FALSE
  )

  expect_s3_class(result_2pl, "irtsimrel_sac_preflight")
  expect_s3_class(result_3pl, "irtsimrel_sac_preflight")
  expect_equal(result_2pl$selected_root, 0.8, tolerance = 2e-5)
  expect_equal(result_3pl$selected_root, 0.7, tolerance = 2e-5)
  expect_equal(result_3pl_msem$selected_root, 0.65, tolerance = 2e-5)
  expect_identical(result_3pl_msem$metric, "msem")
  expect_identical(result_2pl$item_scope, "fixed_form")
  expect_identical(result_3pl$selected_root_record$kind, "crossing")
  expect_identical(
    result_3pl$selected_root_record$direction,
    "increasing"
  )
  expect_true(
    result_3pl$selected_root > result_3pl$effective_branch_bounds[["lower"]]
  )
  expect_true(
    result_3pl$selected_root < result_3pl$effective_branch_bounds[["upper"]]
  )
})


test_that("3PL guessing zero is the exact 2PL preflight special case", {
  bank_2pl <- .sac_preflight_bank(model = "2pl")
  bank_3pl_zero <- .sac_preflight_bank(
    guessing = rep(0, 4L),
    model = "3pl"
  )
  target <- .sac_preflight_rho(0.8, bank_2pl)

  common <- list(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    target_rho = target,
    metric_internal = "info",
    c_bounds = c(0.1, 1.5),
    split_check = FALSE
  )
  result_2pl <- do.call(
    .irtsimrel_sac_preflight,
    c(common, list(item_banks = list(bank_2pl)))
  )
  result_3pl <- do.call(
    .irtsimrel_sac_preflight,
    c(common, list(item_banks = list(bank_3pl_zero)))
  )

  expect_identical(
    result_2pl$full_topology$scan$rho,
    result_3pl$full_topology$scan$rho
  )
  expect_identical(result_2pl$selected_root, result_3pl$selected_root)
  expect_identical(
    result_2pl$effective_branch_bounds,
    result_3pl$effective_branch_bounds
  )
})


test_that("superpopulation split preflight is stable, cached, and reproducible", {
  bank_a <- .sac_preflight_bank(model = "2pl")
  bank_b <- .sac_preflight_bank(
    model = "2pl", beta_shift = 0.08, lambda_multiplier = 1.04
  )
  # Odd and even halves contain the same two forms.
  banks <- list(bank_a, bank_a, bank_b, bank_b)
  target <- mean(vapply(
    banks,
    function(bank) .sac_preflight_rho(0.75, bank),
    numeric(1)
  ))
  args <- list(
    theta = .sac_preflight_theta,
    theta_var = stats::var(.sac_preflight_theta),
    item_banks = banks,
    target_rho = target,
    metric_internal = "info",
    c_bounds = c(0.1, 1.5)
  )

  first <- do.call(.irtsimrel_sac_preflight, args)
  second <- do.call(.irtsimrel_sac_preflight, args)

  expect_identical(first$item_scope, "item_superpopulation")
  expect_true(first$split$checked)
  expect_true(first$split$stable)
  expect_identical(first$split$reason, "stable")
  expect_equal(first$selected_root, 0.75, tolerance = 2e-5)
  expect_equal(first$split$roots, c(odd = 0.75, even = 0.75),
               tolerance = 2e-5)
  expect_identical(first$selected_root, second$selected_root)
  expect_identical(first$split$roots, second$split$roots)
  expect_s3_class(first$split_topologies$odd, "irtsimrel_root_topology")
  expect_s3_class(first$split_topologies$even, "irtsimrel_root_topology")
  expect_gt(first$eval_diagnostics$cache_hit_count, 0L)
  expect_true(all(
    first$eval_diagnostics$per_scale$cache_eval_count == 1L
  ))
  expect_true(all(
    first$eval_diagnostics$per_scale$form_reliability_eval_count == 4L
  ))
  expect_lt(
    first$cache_eval_count,
    sum(unlist(first$eval_diagnostics$scanner_eval_count))
  )
})


test_that("resolved infeasibility is a branch-unavailable condition", {
  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) 0.2 + 0.5 * plogis(log(c)),
      target = 0.9
    ),
    "no resolved feasible root",
    class = "irtsimrel_branch_unavailable"
  )
})


test_that("decreasing, tangent, and boundary-only roots are rejected", {
  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) plogis(-log(c)),
      target = 0.5
    ),
    "no root is available",
    class = "irtsimrel_branch_unavailable"
  )

  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) 0.7 - 0.05 * log(c)^2,
      target = 0.7
    ),
    "no root is available",
    class = "irtsimrel_branch_unavailable"
  )

  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) plogis(log(c)),
      target = plogis(log(0.1))
    ),
    "not an interior increasing crossing",
    class = "irtsimrel_branch_unavailable"
  )
})


test_that("split-root mismatch is topology uncertainty", {
  low_root_bank <- .sac_preflight_bank(model = "2pl", beta_shift = -2)
  high_root_bank <- .sac_preflight_bank(model = "2pl", beta_shift = 2)
  banks <- list(
    low_root_bank, high_root_bank, low_root_bank, high_root_bank
  )
  condition <- tryCatch(
    .run_mocked_sac_preflight(
      objective = function(c, beta_vec, ...) {
        if (beta_vec[[1L]] < -1.5) {
          plogis(log(c / 0.5))
        } else {
          plogis(log(c / 2))
        }
      },
      target = 0.5,
      banks = banks,
      split_check = TRUE,
      split_log_tolerance = 0.5
    ),
    irtsimrel_topology_uncertain = identity
  )

  expect_s3_class(condition, "irtsimrel_topology_uncertain")
  expect_identical(
    condition$topology$split_reason,
    "split_log_root_mismatch"
  )
  expect_equal(condition$topology$split_roots, c(odd = 0.5, even = 2),
               tolerance = 2e-5)
  expect_gt(
    condition$topology$absolute_log_root_ratio,
    condition$topology$split_log_tolerance
  )
})


test_that("nonfinite curves and exhausted scan budgets are uncertain", {
  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) {
        if (c > 0.8 && c < 1.2) return(NA_real_)
        plogis(log(c))
      },
      target = 0.5
    ),
    "not resolved",
    class = "irtsimrel_topology_uncertain"
  )

  expect_error(
    .run_mocked_sac_preflight(
      objective = function(c, ...) plogis(log(c)),
      target = 0.5,
      root_controls = list(min_grid = 65L, max_evals = 10L)
    ),
    "not resolved",
    class = "irtsimrel_topology_uncertain"
  )
})


test_that("preflight validates normalized-bank and split contracts", {
  bank <- .sac_preflight_bank(model = "2pl")
  expect_error(
    .irtsimrel_sac_preflight(
      theta = .sac_preflight_theta,
      theta_var = stats::var(.sac_preflight_theta),
      item_banks = list(unclass(bank)),
      target_rho = 0.5,
      c_bounds = c(0.1, 10),
      split_check = FALSE
    ),
    "created by"
  )
  expect_error(
    .irtsimrel_sac_preflight(
      theta = .sac_preflight_theta,
      theta_var = stats::var(.sac_preflight_theta),
      item_banks = list(bank),
      target_rho = 0.5,
      c_bounds = c(0.1, 10),
      split_check = TRUE
    ),
    "at least four"
  )
})
